// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/plugins/DefaultReportPlugin.cpp (HydroCouple/
// openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin, licensed
// under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The legacy-format .rpt writer, part 1: the plugin's lifecycle (update, finalize,
// write_summary and the destructor's "[Report interrupted ...]" footer), the file-local
// helpers, write_preamble (title, errors and warnings, the input summaries, the analysis
// options) and write_timing. write_results (the continuity blocks and the summary tables)
// is core/report_results.rs.
//
// Differences of form, never in result:
//   * The FILE* is a CFile (core/cprintf.rs); fprintf is cfprintf!, whose formatter is
//     glibc's printf (C locale).
//   * The destructor is Drop; it runs where upstream deletes the plugin (unload_all) or
//     where a temporary plugin goes out of scope (SWMMEngine::write_open_failure_report).
//   * write_timing prints ctime(3) for the wall-clock start and end. The time zone database
//     is not read: the time is printed as UTC (the gate machine's zone). These three lines
//     are the ones the .rpt gate masks.
//   * The 2D blocks of write_results (core/report_results.rs; 12 November 2026, phase 2a wave
//     5): upstream reaches the router through ctx.twod_io.surface_router; here, with the
//     `twod` feature, the engine passes it to write_summary.

use crate::cfprintf;
use crate::core::context::SimulationContext;
use crate::core::cprintf::CFile;
use crate::core::datetime;
use crate::core::plugins::{DefaultReportPlugin, PluginState};
use crate::core::unit_conversion as ucf;
use crate::data::gage::RainSource;
use crate::data::link::{LinkType, XsectShape};
use crate::data::table::TableType;
use crate::hydraulics::transect_types::{TransectData, N_TRANSECT_TBL};

/// OPENSWMM_VERSION_FULL (the build's generated version.h)
pub const OPENSWMM_VERSION_FULL: &[u8] = b"6.0.0-alpha.4";

pub(crate) const FLOW_UNIT_WORDS: [&[u8]; 6] = [b"CFS", b"GPM", b"MGD", b"CMS", b"LPS", b"MLD"];
const INFIL_MODEL_WORDS: [&[u8]; 5] = [b"HORTON", b"MODIFIED_HORTON", b"GREEN_AMPT", b"MODIFIED_GREEN_AMPT", b"CURVE_NUMBER"];
const SURCHARGE_WORDS: [&[u8]; 4] = [b"EXTRAN", b"SLOT", b"DYNAMIC_SLOT", b"TPA"];
const NODE_TYPE_WORDS: [&[u8]; 4] = [b"JUNCTION", b"OUTFALL", b"DIVIDER", b"STORAGE"];
const LINK_TYPE_WORDS: [&[u8]; 5] = [b"CONDUIT", b"PUMP", b"ORIFICE", b"WEIR", b"OUTLET"];
const RAIN_TYPE_WORDS: [&[u8]; 3] = [b"INTENSITY", b"VOLUME", b"CUMULATIVE"];
const XSECT_SHAPE_WORDS: [&[u8]; 26] = [
    b"CIRCULAR",
    b"FILLED_CIRCULAR",
    b"RECT_CLOSED",
    b"RECT_OPEN",
    b"TRAPEZOIDAL",
    b"TRIANGULAR",
    b"PARABOLIC",
    b"POWER",
    b"MODBASKETHANDLE",
    b"EGG",
    b"HORSESHOE",
    b"GOTHIC",
    b"CATENARY",
    b"SEMIELLIPTICAL",
    b"BASKETHANDLE",
    b"SEMICIRCULAR",
    b"RECT_TRIANGULAR",
    b"RECT_ROUND",
    b"HORIZ_ELLIPSE",
    b"VERT_ELLIPSE",
    b"ARCH",
    b"IRREGULAR",
    b"CUSTOM",
    b"FORCE_MAIN",
    b"STREET",
    b"DUMMY",
];

// ---------------------------------------------------------------------------
// Date/time helpers
// ---------------------------------------------------------------------------

/// dateToStr: MM/DD/YYYY, or N/A for a nonpositive date.
pub(crate) fn date_to_str(date: f64) -> Vec<u8> {
    if date <= 0.0 {
        return b"N/A".to_vec();
    }
    let (y, m, d) = datetime::decode_date(date);
    crate::cformat!(b"%02d/%02d/%04d", m, d, y)
}

/// timeToStr: HH:MM:SS.
pub(crate) fn time_to_str(date: f64) -> Vec<u8> {
    let (h, mn, s) = datetime::decode_time(date);
    crate::cformat!(b"%02d:%02d:%02d", h, mn, s)
}

/// secsToHMS: hours may exceed 23.
fn secs_to_hms(secs: i32) -> Vec<u8> {
    let h = secs / 3600;
    let m = (secs % 3600) / 60;
    let s = secs % 60;
    crate::cformat!(b"%02d:%02d:%02d", h, m, s)
}

/// elapsedToParts: days, hours and minutes since the report start, rounded to the second.
pub(crate) fn elapsed_to_parts(date: f64, report_start: f64) -> (i32, i32, i32) {
    let (mut days, mut hrs, mut mins) = (0, 0, 0);
    if date <= 0.0 {
        return (days, hrs, mins);
    }
    let x = date - report_start;
    if x <= 0.0 {
        return (days, hrs, mins);
    }
    days = crate::core::cmath::cvt_i32(x);
    let mut secs = crate::core::cmath::cvt_i32(((x - x.floor()) * 86400.0 + 0.5).floor());
    if secs >= 86400 {
        secs = 86399;
    }
    let total_mins = secs / 60;
    hrs = total_mins / 60;
    mins = total_mins % 60;
    if hrs > 23 {
        hrs = 0;
    }
    (days, hrs, mins)
}

pub(crate) fn nt_str(nt: i32) -> &'static [u8] {
    if (0..=3).contains(&nt) {
        NODE_TYPE_WORDS[nt as usize]
    } else {
        b"JUNCTION"
    }
}
pub(crate) fn lt_str(lt: i32) -> &'static [u8] {
    if (0..=4).contains(&lt) {
        LINK_TYPE_WORDS[lt as usize]
    } else {
        b"CONDUIT"
    }
}

/// WRITE(f, s): fprintf(f, "\n  %s", s)
pub(crate) fn wr(f: &mut CFile, s: &[u8]) {
    cfprintf!(f, b"\n  %s", s);
}

/// flowFmt: "%9.3f" for MGD/CMS, "%9.2f" otherwise.
pub(crate) fn flow_fmt(fu: i32) -> &'static [u8] {
    if fu == 2 || fu == 3 {
        b"%9.3f"
    } else {
        b"%9.2f"
    }
}

/// landAreaToFt2: 1 / UCF(LANDAREA) by the flow units' system.
pub(crate) fn land_area_to_ft2(flow_units: i32) -> f64 {
    let us = if flow_units >= 3 { 1 } else { 0 };
    1.0 / ucf::UCF[ucf::LANDAREA as usize][us]
}

/// unindented: the message after its leading blanks.
pub(crate) fn unindented(m: &[u8]) -> &[u8] {
    match m.iter().position(|&c| c != b' ') {
        Some(k) => &m[k..],
        None => &m[m.len()..],
    }
}

/// write_errors: the input errors from `from` on, capped at 100 as legacy input.c caps them.
pub(crate) fn write_errors(f: &mut CFile, errors: &[Vec<u8>], from: usize) {
    const K_MAX_ERRS: usize = 100;
    let mut i = from;
    while i < errors.len() && i < K_MAX_ERRS {
        cfprintf!(f, b"\n  %s", unindented(&errors[i]));
        i += 1;
    }
    if errors.len() > K_MAX_ERRS && from <= K_MAX_ERRS {
        cfprintf!(f, b"\n  \n  Maximum error count exceeded.");
    }
}

impl Drop for DefaultReportPlugin {
    /// ~DefaultReportPlugin: a report still open was not completed.
    fn drop(&mut self) {
        if let Some(mut f) = self.file_.take() {
            f.write_bytes("\n\n  [Report interrupted \u{2014} simulation did not complete normally]\n".as_bytes());
            f.flush();
        }
    }
}

impl DefaultReportPlugin {
    /// update(snapshot): UPDATING; nothing is written per period.
    pub fn update(&mut self, _snapshot: &crate::core::output_plugin::SimulationSnapshot) -> i32 {
        self.set_state(PluginState::Updating);
        0
    }

    /// finalize(ctx): flush the errors and warnings that accumulated during the run.
    pub fn finalize(&mut self, ctx: &SimulationContext) -> i32 {
        if let Some(f) = self.file_.as_mut() {
            write_errors(f, &ctx.errors, self.errors_written_);
            self.errors_written_ = ctx.errors.len();
            let mut i = self.warnings_written_;
            while i < ctx.warnings.len() {
                cfprintf!(f, b"\n  %s", unindented(&ctx.warnings[i]));
                i += 1;
            }
            self.warnings_written_ = ctx.warnings.len();
            f.flush();
        }
        self.set_state(PluginState::Finalized);
        0
    }

    /// write_summary(ctx): the result sections, the snow balance, the timing, then fclose.
    /// With `twod`, `router` is what upstream's ctx.twod_io.surface_router points at.
    pub fn write_summary(&mut self, ctx: &SimulationContext, #[cfg(feature = "twod")] router: Option<&crate::twod::surface_router_2d::SurfaceRouter2D>) -> i32 {
        if self.rpt_path().is_empty() {
            return 0;
        }
        if ctx.options.rpt_disabled {
            return 0;
        }
        let mut f = match self.file_.take() {
            Some(f) => f,
            None => {
                // prepare() was not called or its fopen failed: the whole report at once
                let mut f = match CFile::create(self.rpt_path(), false) {
                    Some(f) => f,
                    None => return -1,
                };
                self.write_preamble(&mut f, ctx);
                write_errors(&mut f, &ctx.errors, 0);
                for warn in &ctx.warnings {
                    cfprintf!(f, b"\n  %s", unindented(warn));
                }
                f
            }
        };

        #[cfg(not(feature = "twod"))]
        self.write_results(&mut f, ctx);
        #[cfg(feature = "twod")]
        self.write_results(&mut f, ctx, router);

        if !ctx.snowpacks.report_values.is_empty() {
            cfprintf!(f, b"\n\n  Snow Process Water and Energy Balance (SI)\n  -----------------------------------------\n");
            for row in &ctx.snowpacks.report_values {
                cfprintf!(f, b"\n  %-44s % .9g", row.0, row.1);
            }
            cfprintf!(f, b"\n  Negative sublimation is condensation; energy includes the active ground layer.\n");
        }

        self.write_timing(&mut f, ctx);
        drop(f); // fclose
        0
    }

    /// write_preamble(f, ctx): title, errors and warnings, input summaries, options.
    pub(crate) fn write_preamble(&mut self, f: &mut CFile, ctx: &SimulationContext) {
        let opt = &ctx.options;
        let mut fu = opt.flow_units as i32;
        if !(0..=5).contains(&fu) {
            fu = 0;
        }
        let du_opt = ucf::DisplayUnits::from(opt);

        cfprintf!(f, b"\n  OPENSWMM ENGINE - VERSION %s", OPENSWMM_VERSION_FULL);
        cfprintf!(f, b"\n  ------------------------------------------------------------");
        cfprintf!(f, b"\n");

        // Warnings raised while reading the input carry their own leading "\n  " and
        // precede the title.
        for warn in &ctx.warnings {
            if !warn.is_empty() && warn[0] == b'\n' {
                cfprintf!(f, b"\n  %s", warn);
            }
        }

        // [TITLE]: at most 3 lines, comment lines skipped
        {
            let mut n_title = 0;
            for line in &ctx.title_notes {
                let first = line.iter().position(|&c| c != b' ' && c != b'\t');
                match first {
                    None => continue,
                    Some(k) if line[k] == b';' => continue,
                    _ => {}
                }
                let was = n_title;
                n_title += 1;
                if was == 3 {
                    break;
                }
                cfprintf!(f, b"\n  %s ", line);
            }
        }

        write_errors(f, &ctx.errors, 0);
        self.errors_written_ = ctx.errors.len();

        for warn in &ctx.warnings {
            if warn.is_empty() || warn[0] != b'\n' {
                cfprintf!(f, b"\n  %s", unindented(warn));
            }
        }
        self.warnings_written_ = ctx.warnings.len();

        if opt.rpt_input {
            wr(f, b"");
            wr(f, b"*************");
            wr(f, b"Element Count");
            wr(f, b"*************");
            cfprintf!(f, b"\n  Number of rain gages ...... %d", ctx.n_gages());
            cfprintf!(f, b"\n  Number of subcatchments ... %d", ctx.n_subcatches());
            cfprintf!(f, b"\n  Number of nodes ........... %d", ctx.n_nodes());
            cfprintf!(f, b"\n  Number of links ........... %d", ctx.n_links());
            cfprintf!(f, b"\n  Number of pollutants ...... %d", ctx.n_pollutants());
            cfprintf!(f, b"\n  Number of land uses ....... %d", ctx.landuse_names.size());

            wr(f, b"");
            wr(f, b"");

            if ctx.n_pollutants() > 0 {
                const K_QUAL_UNITS: [&[u8]; 3] = [b"MG/L", b"UG/L", b"#/L"];
                wr(f, b"*****************");
                wr(f, b"Pollutant Summary");
                wr(f, b"*****************");
                cfprintf!(f, b"\n                               Ppt.      GW         Kdecay");
                cfprintf!(f, b"\n  Name                 Units   Concen.   Concen.    1/days    CoPollutant");
                cfprintf!(f, b"\n  -----------------------------------------------------------------------");
                let p = &ctx.pollutants;
                for i in 0..ctx.n_pollutants() {
                    let ui = i as usize;
                    let u = p.units[ui] as i32;
                    cfprintf!(
                        f,
                        b"\n  %-20s %5s%10.2f%10.2f%10.2f",
                        ctx.pollutant_names.name_of(i),
                        if (0..=2).contains(&u) { K_QUAL_UNITS[u as usize] } else { b"MG/L" as &[u8] },
                        p.c_rain[ui],
                        p.c_gw[ui],
                        p.k_decay[ui] * 86400.0
                    );
                    if p.co_pollut[ui] >= 0 {
                        cfprintf!(f, b"    %-s  (%.2f)", ctx.pollutant_names.name_of(p.co_pollut[ui]), p.co_frac[ui]);
                    }
                }
                wr(f, b"");
                wr(f, b"");
            }
            if ctx.landuse_names.size() > 0 {
                wr(f, b"***************");
                wr(f, b"Landuse Summary");
                wr(f, b"***************");
                cfprintf!(f, b"\n                         Sweeping   Maximum      Last");
                cfprintf!(f, b"\n  Name                   Interval   Removal     Swept");
                cfprintf!(f, b"\n  ---------------------------------------------------");
                for i in 0..ctx.landuses.count() {
                    let ui = i as usize;
                    cfprintf!(f, b"\n  %-20s %10.2f%10.2f%10.2f", ctx.landuse_names.name_of(i), ctx.landuses.sweep_interval[ui], ctx.landuses.sweep_removal[ui], ctx.landuses.last_swept[ui]);
                }
                wr(f, b"");
                wr(f, b"");
            }

            if ctx.n_gages() > 0 {
                wr(f, b"****************");
                wr(f, b"Raingage Summary");
                wr(f, b"****************");
                cfprintf!(f, b"\n                                                      Data       Recording");
                cfprintf!(f, b"\n  Name                 Data Source                    Type       Interval ");
                cfprintf!(f, b"\n  ------------------------------------------------------------------------");
                for i in 0..ctx.n_gages() {
                    let ui = i as usize;
                    let name = ctx.gage_names.name_of(i);
                    let rt = ctx.gages.rain_type[ui];
                    let interval_min = ctx.gages.interval_sec[ui] / 60;
                    let rt_str: &[u8] = if (0..=2).contains(&rt) { RAIN_TYPE_WORDS[rt as usize] } else { b"INTENSITY" };
                    if ctx.gages.source[ui] == RainSource::Timeseries {
                        cfprintf!(f, b"\n  %-20s %-30s %-10s %3d min.", name, ctx.gages.ts_name[ui], rt_str, interval_min);
                    } else {
                        cfprintf!(f, b"\n  %-20s %-30s %-10s %3d min.", name, ctx.gages.file_path[ui].original, rt_str, interval_min);
                    }
                }
                wr(f, b"");
                wr(f, b"");
            }

            if ctx.n_subcatches() > 0 {
                wr(f, b"********************");
                wr(f, b"Subcatchment Summary");
                wr(f, b"********************");
                cfprintf!(f, b"\n  Name                       Area     Width   %%Imperv    %%Slope Rain Gage            Outlet              ");
                cfprintf!(f, b"\n  -----------------------------------------------------------------------------------------------------------");
                for i in 0..ctx.n_subcatches() {
                    let ui = i as usize;
                    let name = ctx.subcatch_names.name_of(i);
                    let gage_idx = ctx.subcatches.gage[ui];
                    let gage_name: &[u8] = if gage_idx >= 0 { ctx.gage_names.name_of(gage_idx) } else { b"" };
                    let out_node = ctx.subcatches.outlet_node[ui];
                    let mut outlet_name: &[u8] = if out_node >= 0 { ctx.node_names.name_of(out_node) } else { b"" };
                    let out_sub = ctx.subcatches.outlet_subcatch[ui];
                    if out_node < 0 && out_sub >= 0 {
                        outlet_name = ctx.subcatch_names.name_of(out_sub);
                    }
                    cfprintf!(
                        f,
                        b"\n  %-20s %10.2f%10.2f%10.2f%10.4f %-20s %-20s",
                        name,
                        (ctx.subcatches.area[ui] / ucf::ucf(ucf::LANDAREA, opt)) * ucf::ucf(ucf::LANDAREA, opt),
                        ctx.subcatches.width[ui] * du_opt.length,
                        ctx.subcatches.frac_imperv[ui] * 100.0,
                        ctx.subcatches.slope[ui] * 100.0,
                        gage_name,
                        outlet_name
                    );
                }
                wr(f, b"");
                wr(f, b"");
            }

            // LID Control Summary (legacy lid_writeSummary), [LID_USAGE] order reversed
            {
                let lu = &ctx.lid_usage;
                let mut any_lid = false;
                let mut k = 0;
                while k < lu.count() && !any_lid {
                    any_lid = lu.area[k as usize] * lu.number[k as usize] as f64 > 0.0;
                    k += 1;
                }
                if any_lid {
                    let ucf_len = du_opt.length;
                    let ucf_len2 = ucf_len * ucf_len;
                    let la_ft2 = land_area_to_ft2(fu);
                    cfprintf!(f, b"\n  *******************");
                    cfprintf!(f, b"\n  LID Control Summary");
                    cfprintf!(f, b"\n  *******************");
                    cfprintf!(f, b"\n                                   No. of        Unit        Unit      %% Area    %% Imperv      %% Perv");
                    cfprintf!(f, b"\n  Subcatchment     LID Control      Units        Area       Width     Covered     Treated     Treated");
                    cfprintf!(f, b"\n  ---------------------------------------------------------------------------------------------------");
                    for j in 0..ctx.n_subcatches() {
                        let sub_area = ctx.subcatches.area[j as usize] * la_ft2;
                        let mut k = lu.count() - 1;
                        while k >= 0 {
                            let uk = k as usize;
                            if lu.subcatch_index[uk] != j {
                                k -= 1;
                                continue;
                            }
                            let area = lu.area[uk] / ucf_len2;
                            let width = lu.width[uk] / ucf_len;
                            let lid = lu.lid_index[uk];
                            cfprintf!(f, b"\n  %-16s %-16s", ctx.subcatch_names.name_of(j), if lid >= 0 { ctx.lid_names.name_of(lid) } else { b"" as &[u8] });
                            cfprintf!(
                                f,
                                b"%6d  %10.2f  %10.2f  %10.2f  %10.2f  %10.2f",
                                lu.number[uk],
                                area * ucf_len2,
                                width * ucf_len,
                                area * lu.number[uk] as f64 / sub_area * 100.0,
                                lu.from_imperv[uk] / 100.0 * 100.0,
                                lu.from_perv[uk] / 100.0 * 100.0
                            );
                            k -= 1;
                        }
                    }
                    wr(f, b"");
                    wr(f, b"");
                }
            }

            if ctx.n_nodes() > 0 {
                wr(f, b"************");
                wr(f, b"Node Summary");
                wr(f, b"************");
                cfprintf!(f, b"\n                                           Invert      Max.    Ponded    External");
                cfprintf!(f, b"\n  Name                 Type                 Elev.     Depth      Area    Inflow  ");
                cfprintf!(f, b"\n  -------------------------------------------------------------------------------");
                let nn = ctx.n_nodes();
                let mut has_inflow = vec![0u8; nn as usize];
                let mut mark = |node_idx: &[i32]| {
                    for &n in node_idx {
                        if n >= 0 && n < nn {
                            has_inflow[n as usize] = 1;
                        }
                    }
                };
                mark(&ctx.ext_inflows.node_idx);
                mark(&ctx.dwf_inflows.node_idx);
                mark(&ctx.rdii_assigns.node_idx);
                for i in 0..nn {
                    let ui = i as usize;
                    let nt = ctx.nodes.type_[ui] as i32;
                    cfprintf!(
                        f,
                        b"\n  %-20s %-16s%10.2f%10.2f%10.1f",
                        ctx.node_names.name_of(i),
                        nt_str(nt),
                        ctx.nodes.invert_elev[ui] * du_opt.length,
                        ctx.nodes.full_depth[ui] * du_opt.length,
                        ctx.nodes.ponded_area[ui] * du_opt.length * du_opt.length
                    );
                    if has_inflow[ui] != 0 {
                        cfprintf!(f, b"    Yes");
                    }
                }
                wr(f, b"");
                wr(f, b"");
            }

            if ctx.n_links() > 0 {
                wr(f, b"************");
                wr(f, b"Link Summary");
                wr(f, b"************");
                cfprintf!(f, b"\n  Name             From Node        To Node          Type            Length    %%Slope Roughness");
                cfprintf!(f, b"\n  ---------------------------------------------------------------------------------------------");
                for i in 0..ctx.n_links() {
                    let ui = i as usize;
                    let lt = ctx.links.type_[ui] as i32;
                    let (mut n1, mut n2) = (ctx.links.node1[ui], ctx.links.node2[ui]);
                    let dir = ctx.links.direction[ui];
                    if dir < 0 {
                        std::mem::swap(&mut n1, &mut n2);
                    }
                    let n1_name: &[u8] = if n1 >= 0 { ctx.node_names.name_of(n1) } else { b"" };
                    let n2_name: &[u8] = if n2 >= 0 { ctx.node_names.name_of(n2) } else { b"" };
                    cfprintf!(f, b"\n  %-16s %-16s %-16s ", ctx.link_names.name_of(i), n1_name, n2_name);

                    if lt == LinkType::Conduit as i32 {
                        let cr = ctx.link_subtypes.conduit_row(i);
                        let cd = &ctx.link_subtypes.conduits;
                        let mut n_rep = if cr >= 0 { cd.roughness[cr as usize] } else { 0.01 };
                        if cr >= 0 && ctx.links.xsect_shape[ui] == XsectShape::Irregular {
                            let ti = ctx.links.xsect_curve[ui];
                            if ti >= 0 && (ti as usize) < ctx.transects.n_channel.len() && ctx.transects.n_channel[ti as usize] > 0.0 {
                                n_rep = ctx.transects.n_channel[ti as usize];
                            }
                        }
                        cfprintf!(
                            f,
                            b"%-12s%10.1f%10.4f%10.4f",
                            b"CONDUIT",
                            (if cr >= 0 { cd.length[cr as usize] } else { 0.0 }) * du_opt.length,
                            (if cr >= 0 { cd.slope[cr as usize] } else { 0.0 }) * 100.0 * dir as f64,
                            n_rep
                        );
                    } else if lt == LinkType::Pump as i32 {
                        const K_PUMP_WORDS: [&[u8]; 5] = [b"TYPE1", b"TYPE2", b"TYPE3", b"TYPE4", b"TYPE5"];
                        let pr = ctx.link_subtypes.pump_row(i);
                        let ct = if pr >= 0 { ctx.link_subtypes.pumps.curve_type[pr as usize] } else { -1 };
                        cfprintf!(f, b"%-5s PUMP  ", if (1..=5).contains(&ct) { K_PUMP_WORDS[(ct - 1) as usize] } else { b"IDEAL" as &[u8] });
                    } else {
                        cfprintf!(f, b"%-12s", lt_str(lt));
                    }
                }
                wr(f, b"");
                wr(f, b"");
            }

            if ctx.n_links() > 0 {
                wr(f, b"*********************");
                wr(f, b"Cross Section Summary");
                wr(f, b"*********************");
                cfprintf!(f, b"\n                                        Full     Full     Hyd.     Max.   No. of     Full");
                cfprintf!(f, b"\n  Conduit          Shape               Depth     Area     Rad.    Width  Barrels     Flow");
                cfprintf!(f, b"\n  ---------------------------------------------------------------------------------------");
                let qcf_pre = ucf::QCF[fu as usize];
                for i in 0..ctx.n_links() {
                    let ui = i as usize;
                    if ctx.links.type_[ui] != LinkType::Conduit {
                        continue;
                    }
                    let shape = ctx.links.xsect_shape[ui] as i32;
                    let mut shape_str: &[u8] = if (0..=25).contains(&shape) { XSECT_SHAPE_WORDS[shape as usize] } else { b"CIRCULAR" };
                    let xs = ctx.links.xsect_shape[ui];
                    if (xs == XsectShape::Custom || xs == XsectShape::Irregular || xs == XsectShape::StreetXsect) && !ctx.links.pump_curve_name[ui].is_empty() {
                        shape_str = &ctx.links.pump_curve_name[ui];
                    }
                    let cr = ctx.link_subtypes.conduit_row(i);
                    let cd = &ctx.link_subtypes.conduits;
                    cfprintf!(
                        f,
                        b"\n  %-16s %-16s %8.2f %8.2f %8.2f %8.2f      %3d %8.2f",
                        ctx.link_names.name_of(i),
                        shape_str,
                        ctx.links.xsect_y_full[ui] * du_opt.length,
                        ctx.links.xsect_a_full[ui] * du_opt.length * du_opt.length,
                        ctx.links.xsect_r_full[ui] * du_opt.length,
                        ctx.links.xsect_w_max[ui] * du_opt.length,
                        if cr >= 0 { cd.barrels[cr as usize] } else { 1 },
                        (if cr >= 0 { cd.q_full[cr as usize] } else { 0.0 }) * qcf_pre
                    );
                }
                wr(f, b"");
                wr(f, b"");
            }

            // Shape / Transect / Street Summary: the normalized tables, five values a row
            {
                fn tables(f: &mut CFile, kind: &[u8], name: &[u8], td: &TransectData) {
                    cfprintf!(f, b"\n\n  %s %s", kind, name);
                    let tbl: [&[f64]; 3] = [&td.area_tbl, &td.hrad_tbl, &td.width_tbl];
                    let lbl: [&[u8]; 3] = [b"Area:  ", b"Hrad:  ", b"Width: "];
                    for t in 0..3 {
                        cfprintf!(f, b"\n  %s", lbl[t]);
                        for m in 1..N_TRANSECT_TBL as usize {
                            if m % 5 == 1 {
                                cfprintf!(f, b"\n          ");
                            }
                            cfprintf!(f, b"%10.4f ", tbl[t][m]);
                        }
                    }
                }

                let any_shape = ctx.tables.tables.iter().any(|t| t.type_ == TableType::CurveShape);
                if any_shape {
                    wr(f, b"*************");
                    wr(f, b"Shape Summary");
                    wr(f, b"*************");
                    for t in &ctx.tables.tables {
                        if t.type_ != TableType::CurveShape {
                            continue;
                        }
                        let mut td = TransectData::default();
                        if !t.x.is_empty() {
                            crate::hydraulics::transect::build_custom_tables(&mut td, 1.0, &t.x, &t.y, t.x.len() as i32);
                        }
                        tables(f, b"Shape", &t.id, &td);
                    }
                    wr(f, b"");
                    wr(f, b"");
                }

                let nt = ctx.transects.count();
                if nt > 0 {
                    wr(f, b"****************");
                    wr(f, b"Transect Summary");
                    wr(f, b"****************");
                    let mut i = 0;
                    while i < nt && (i as usize) < ctx.transect_tables.len() {
                        let td = &ctx.transect_tables[i as usize];
                        tables(f, b"Transect", &td.name, td);
                        i += 1;
                    }
                    wr(f, b"");
                    wr(f, b"");
                }

                if ctx.streets.count() > 0 {
                    wr(f, b"**************");
                    wr(f, b"Street Summary");
                    wr(f, b"**************");
                    let us = if fu >= 3 { 1 } else { 0 };
                    let ucf_len = ucf::UCF[ucf::LENGTH as usize][us];
                    for s in 0..ctx.streets.count() {
                        let su = s as usize;
                        let mut sp = crate::hydraulics::street::StreetParams::default();
                        sp.width = ctx.streets.t_crown[su] / ucf_len;
                        sp.curb_height = ctx.streets.h_curb[su] / ucf_len;
                        sp.slope = ctx.streets.sx[su] / 100.0;
                        sp.roughness = ctx.streets.n_road[su];
                        sp.gutter_depression = ctx.streets.gutter_depres[su] / ucf_len;
                        sp.gutter_width = ctx.streets.gutter_width[su] / ucf_len;
                        sp.sides = ctx.streets.sides[su];
                        sp.back_width = ctx.streets.back_width[su] / ucf_len;
                        sp.back_slope = ctx.streets.back_slope[su] / 100.0;
                        sp.back_roughness = ctx.streets.back_n[su];
                        let mut td = TransectData::default();
                        crate::hydraulics::street::build_transect(&sp, &mut td);
                        tables(f, b"Street", &ctx.streets.names[su], &td);
                    }
                    wr(f, b"");
                    wr(f, b"");
                }
            }
        } // end rpt_input

        // Analysis Options (legacy report_writeOptions)
        wr(f, b"");
        wr(f, b"****************");
        wr(f, b"Analysis Options");
        wr(f, b"****************");

        cfprintf!(f, b"\n  Flow Units ............... %s", FLOW_UNIT_WORDS[fu as usize]);

        let yn = |b: bool| -> &'static [u8] {
            if b {
                b"YES"
            } else {
                b"NO"
            }
        };
        cfprintf!(f, b"\n  Process Models:");
        cfprintf!(f, b"\n    Rainfall/Runoff ........ %s", yn(ctx.n_gages() > 0 && !opt.ignore_rainfall));
        cfprintf!(f, b"\n    RDII ................... %s", yn(ctx.unit_hyds.count() > 0 && !opt.ignore_rdii));
        cfprintf!(f, b"\n    Snowmelt ............... %s", yn(ctx.snowpack_names.size() > 0 && !opt.ignore_snow_melt));
        cfprintf!(f, b"\n    Groundwater ............ %s", yn(ctx.aquifer_names.size() > 0 && !opt.ignore_groundwater));
        cfprintf!(f, b"\n    Flow Routing ........... %s", yn(ctx.n_links() > 0 && !opt.ignore_routing));
        if ctx.n_links() > 0 {
            cfprintf!(f, b"\n    Ponding Allowed ........ %s", yn(opt.allow_ponding));
        }
        cfprintf!(f, b"\n    Water Quality .......... %s", yn(ctx.n_pollutants() > 0 && !opt.ignore_quality));

        // E2: the Domain x Species transport matrix, only when a species class legacy does
        // not have is in use
        if ctx.options.water_age || ctx.options.heat_transport || (ctx.reactions.configured && ctx.reactions.compiled) {
            let block = crate::transport::transport_policy::report_block(ctx);
            cfprintf!(f, b"\n%s", block);
        }

        if ctx.n_subcatches() > 0 {
            let mut im = opt.infiltration as i32;
            if !(0..=4).contains(&im) {
                im = 0;
            }
            cfprintf!(f, b"\n  Infiltration Method ...... %s", INFIL_MODEL_WORDS[im as usize]);
        }

        if ctx.n_links() > 0 {
            let rm = opt.routing_model as i32;
            let rm_name: &[u8] = if rm == 3 {
                b"FV"
            } else if rm == 2 {
                b"DYNWAVE"
            } else if rm == 1 {
                b"KINWAVE"
            } else {
                b"STEADY"
            };
            cfprintf!(f, b"\n  Flow Routing Method ...... %s", rm_name);
            if rm == 3 && opt.fv.pressure_closure == 1 {
                cfprintf!(f, b"\n  Pressure Closure ......... TPA");
            }
            if (rm == 2 || rm == 3) && opt.unsteady_friction != 0 {
                cfprintf!(f, b"\n  Unsteady Friction ........ VITKOVSKY (k3 = %g)", opt.uf_k3);
            }
        }

        if opt.routing_model as i32 == 2 {
            let sm = opt.surcharge_method;
            let sm_name: &[u8] = if (0..=3).contains(&sm) { SURCHARGE_WORDS[sm as usize] } else { b"EXTRAN" };
            cfprintf!(f, b"\n  Surcharge Method ......... %s", sm_name);
            if opt.node_continuity == crate::core::options::NodeContinuity::SemiImplicit {
                cfprintf!(f, b"\n  Node Continuity .......... SEMI_IMPLICIT");
            }
            if opt.anderson_accel {
                cfprintf!(f, b"\n  Anderson Acceleration .... YES");
            }
        }

        cfprintf!(f, b"\n  Starting Date ............ %s %s", date_to_str(opt.start_date), time_to_str(opt.start_date));
        cfprintf!(f, b"\n  Ending Date .............. %s %s", date_to_str(opt.end_date), time_to_str(opt.end_date));
        cfprintf!(f, b"\n  Antecedent Dry Days ...... %.1f", opt.dry_days);

        cfprintf!(f, b"\n  Report Time Step ......... %s", secs_to_hms(crate::core::cmath::cvt_i32(opt.report_step) % 86400));

        if ctx.n_subcatches() > 0 {
            cfprintf!(f, b"\n  Wet Time Step ............ %s", secs_to_hms(crate::core::cmath::cvt_i32(opt.wet_step) % 86400));
            cfprintf!(f, b"\n  Dry Time Step ............ %s", secs_to_hms(crate::core::cmath::cvt_i32(opt.dry_step) % 86400));
        }

        if ctx.n_links() > 0 {
            cfprintf!(f, b"\n  Routing Time Step ........ %.2f sec", opt.routing_step);
            if opt.routing_model as i32 == 2 {
                cfprintf!(f, b"\n  Variable Time Step ....... %s", yn(opt.variable_step > 0.0));
                cfprintf!(f, b"\n  Maximum Trials ........... %d", if opt.max_trials > 0 { opt.max_trials } else { 8 });
                cfprintf!(f, b"\n  Number of Threads ........ %d", opt.num_threads);
                cfprintf!(
                    f,
                    b"\n  Head Tolerance ........... %f %s",
                    if opt.head_tol == 0.0 { 0.005 * du_opt.length } else { opt.head_tol },
                    if du_opt.unit_system == 1 { b"m" as &[u8] } else { b"ft" }
                );
            }
        }

        wr(f, b"");
        wr(f, b"");

        let scaled_rain = ctx.gages.scale_factor.iter().any(|&factor| factor != 1.0);
        if scaled_rain {
            wr(f, b"************************");
            wr(f, b"Rainfall Scaling Factors");
            wr(f, b"************************");
            wr(f, b"Configured factors multiply rain gage input.");
            cfprintf!(f, b"\n  Rain Gage             Scale Factor");
            cfprintf!(f, b"\n  ----------------------------------");
            for i in 0..ctx.n_gages() {
                cfprintf!(f, b"\n  %-20s %.17g", ctx.gage_names.name_of(i), ctx.gages.scale_factor[i as usize]);
            }
            wr(f, b"");
            wr(f, b"");
        }
    }

    /// write_timing(f, ctx): the analysis's wall-clock window (legacy report_writeSysTime).
    pub(crate) fn write_timing(&mut self, f: &mut CFile, ctx: &SimulationContext) {
        let wall_end = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs() as i64).unwrap_or(0);
        let wall_start = if ctx.wall_start != 0 { ctx.wall_start } else { wall_end };
        let begin_str = ctime_utc(wall_start);
        let end_str = ctime_utc(wall_end);
        cfprintf!(f, b"\n\n  Analysis begun on:  %s", begin_str);
        cfprintf!(f, b"\n  Analysis ended on:  %s", end_str);
        cfprintf!(f, b"\n  Total elapsed time: ");
        let elapsed_secs = (wall_end - wall_start) as f64;
        if elapsed_secs < 1.0 {
            cfprintf!(f, b"< 1 sec");
        } else {
            let es = elapsed_secs as i64;
            let days = es / 86400;
            let rem = es % 86400;
            if days > 0 {
                cfprintf!(f, b"%ld.", days);
            }
            cfprintf!(f, b"%02ld:%02ld:%02ld", rem / 3600, (rem % 3600) / 60, rem % 60);
        }
        cfprintf!(f, b"\n");
    }
}

/// ctime(3) without its newline, for a time in seconds since the epoch, as UTC.
fn ctime_utc(t: i64) -> Vec<u8> {
    const WDAY: [&[u8]; 7] = [b"Thu", b"Fri", b"Sat", b"Sun", b"Mon", b"Tue", b"Wed"];
    const MON: [&[u8]; 12] = [b"Jan", b"Feb", b"Mar", b"Apr", b"May", b"Jun", b"Jul", b"Aug", b"Sep", b"Oct", b"Nov", b"Dec"];
    let days = t.div_euclid(86400);
    let secs = t.rem_euclid(86400);
    // civil_from_days (Howard Hinnant)
    let z = days + 719468;
    let era = z.div_euclid(146097);
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let mut y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    if m <= 2 {
        y += 1;
    }
    crate::cformat!(b"%s %s %2d %02d:%02d:%02d %d", WDAY[days.rem_euclid(7) as usize], MON[(m - 1) as usize], d, secs / 3600, (secs % 3600) / 60, secs % 60, y)
}
