// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/SWMMEngine.cpp (postOutputSnapshot, the
// file-local convertSnapshotToDisplay, fillSurfaceSnapshot, accumulateAvgResults,
// applyAvgResults, end, writeDueHotstartSaves, writeHotstartSave, legacyRunoffSection and
// the file-local infilToLegacy, report, close, write_open_failure_report),
// src/engine/core/HotStartManager.cpp (save_legacy_routing), src/engine/core/
// InterfaceFile.cpp (InterfaceManager::writeOutfallResults), src/engine/hydraulics/Inlet.cpp
// (InletSolver::gatherStats), src/engine/quality/NegativeSources.hpp
// (summarizeNegativeSourceClamps) and src/engine/core/openswmm_snow_impl.cpp
// (swmm_snow_get_ueb_ledger), src/engine/hydrology/surface/NativeSnow.cpp
// (nativeSnowResults) and UebSnow.cpp (bulkTemperature) (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The output side of the run and its end (1D build: OPENSWMM_HAS_2D off): the report-period
// snapshot and its writers, time-step averaging, SAVE HOTSTART, end(), report(), close().
//
// Differences of form, never in result:
//   * The IO thread (io_thread_.post, .stop) moves the plugin writes off the stepping
//     thread upstream; the writes happen in post order either way (one consumer, a FIFO
//     queue drained before end() finalizes), so here post_output_snapshot calls the
//     plugins' update directly.
//   * fillSurfaceSnapshot is empty without OPENSWMM_HAS_2D, as upstream's is; with the `twod`
//     feature it and the other OPENSWMM_HAS_2D blocks of these functions are in
//     core/engine_output_twod.rs (12 November 2026, phase 2a wave 5).
//   * The perf:: timers and the OPENSWMM_PERF stderr line only time the run; left out.
//   * The FV solver's run statistics (end()) and the FV router are not translated: a deck
//     with ROUTING FV stops in stepRouting before end() is reached, and end() returns
//     Err(stub) for an engine with an FV backend.
//   * NOT YET TRANSLATED, an Err(stub) where it is met: a NATIVE SAVE HOTSTART row
//     (HotStartManager::save and the engine checkpoint).

use crate::core::context::{EngineState, HotstartFormat, HotstartSaveEntry, SimulationContext};
use crate::core::datetime;
use crate::core::engine::{SwmmEngine, SWMM_ERR_PLUGIN, SWMM_ERR_WRONG_STATE, SWMM_OK};
use crate::core::engine_init_hydraulics::Terminated;
use crate::core::engine_step::stub;
use crate::core::options::Bytes;
use crate::core::output_plugin::SimulationSnapshot;
use crate::core::plugins::{DefaultReportPlugin, PluginFactory};
use crate::core::unit_conversion as ucf;
use crate::data::link::LinkType;
use crate::data::node::NodeType;
use crate::hydraulics::{link, timestep_controller as tc};

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut m = Vec::new();
    for p in parts {
        m.extend_from_slice(p);
    }
    m
}

/// plugin_io_error: the first plugin's last error after the fallback, or the fallback.
fn plugin_io_error(plugins: &PluginFactory, fallback: &[u8]) -> Bytes {
    let mut details: Vec<Bytes> = Vec::new();
    details.push(plugins.output_plugins().iter().map(|p| crate::core::cstdlib::c_str(p.last_error_message()).to_vec()).find(|m| !m.is_empty()).unwrap_or_default());
    details.push(plugins.report_plugins().iter().map(|p| crate::core::cstdlib::c_str(p.last_error_message()).to_vec()).find(|m| !m.is_empty()).unwrap_or_default());
    details.push(plugins.state_io_plugins().iter().map(|p| crate::core::cstdlib::c_str(p.last_error_message()).to_vec()).find(|m| !m.is_empty()).unwrap_or_default());
    for m in details {
        if !m.is_empty() {
            return cat(&[fallback, b": ", &m]);
        }
    }
    fallback.to_vec()
}

/// convertSnapshotToDisplay: the 1D snapshot from internal to display units, in place.
fn convert_snapshot_to_display(s: &mut SimulationSnapshot, du: &ucf::DisplayUnits) {
    fn scale(v: &mut [f64], f: f64) {
        if f == 1.0 {
            return;
        }
        for x in v.iter_mut() {
            *x *= f;
        }
    }
    scale(&mut s.subcatch.rainfall, du.rainfall);
    scale(&mut s.subcatch.snow_depth, du.raindepth);
    scale(&mut s.subcatch.evap, du.evaprate);
    scale(&mut s.subcatch.infil, du.rainfall);
    scale(&mut s.subcatch.runoff, du.flow);
    scale(&mut s.subcatch.gw_flow, du.flow);
    scale(&mut s.subcatch.gw_elev, du.length);
    scale(&mut s.nodes.depth, du.length);
    scale(&mut s.nodes.head, du.length);
    scale(&mut s.nodes.volume, du.volume);
    scale(&mut s.nodes.lateral_inflow, du.flow);
    scale(&mut s.nodes.total_inflow, du.flow);
    scale(&mut s.nodes.overflow, du.flow);
    scale(&mut s.links.flow, du.flow);
    scale(&mut s.links.depth, du.length);
    scale(&mut s.links.velocity, du.length);
    scale(&mut s.links.volume, du.volume);
    if s.has_subcatchments {
        s.sys_temperature = du.temperature(s.sys_temperature);
    }
    s.sys_rainfall *= du.rainfall;
    s.sys_snow_depth *= du.raindepth;
    s.sys_infil *= du.rainfall;
    s.sys_runoff *= du.flow;
    s.sys_dw_inflow *= du.flow;
    s.sys_gw_inflow *= du.flow;
    s.sys_ii_inflow *= du.flow;
    s.sys_ext_inflow *= du.flow;
    s.sys_flooding *= du.flow;
    s.sys_outflow *= du.flow;
    s.sys_storage *= du.volume;
    s.sys_evap *= du.evaprate;
    s.sys_pet *= du.evaprate;
}

/// infilToLegacy: the six native infiltration-state slots in the legacy runoff layout.
fn infil_to_legacy(model: i32, v6: &[f64; 6], lg: &mut [f64; 6]) {
    use crate::hydrology::infiltration::InfilModel;
    *lg = [0.0; 6];
    if model == InfilModel::Horton as i32 || model == InfilModel::ModHorton as i32 {
        lg[0] = v6[0];
        lg[1] = v6[1];
        lg[2] = v6[2];
    } else if model == InfilModel::GreenAmpt as i32 || model == InfilModel::ModGreenAmpt as i32 {
        lg[0] = v6[0];
        lg[1] = v6[1];
        lg[2] = v6[2];
        lg[3] = v6[4];
        lg[4] = v6[3];
    } else if model == InfilModel::CurveNum as i32 {
        lg[0] = v6[0];
        lg[1] = v6[2];
        lg[2] = v6[3];
        lg[3] = v6[5];
        lg[4] = v6[1];
        lg[5] = v6[4];
    }
}

/// HotStartManager::save_legacy_routing: a legacy "SWMM5-HOTSTART4" file (hotstart.c
/// saveHotstart + saveRunoff + saveRouting), with the optional reaction-species trailer.
/// The 2D-only refusals (layered groundwater, coupled thermal) are false in a 1D build.
fn save_legacy_routing(path: &[u8], ctx: &SimulationContext, runoff: Option<&Vec<f64>>) -> i32 {
    use crate::core::hot_start_manager as hsm;
    if !ctx.reactions.lid_controls.is_empty() {
        hsm::set_last_io_error(b"LID reaction inventories require a full native engine checkpoint".to_vec());
        return 1;
    }
    if hsm::mesh_lid_restart_unsupported(ctx) {
        hsm::set_last_io_error(b"Mesh LID hotstart requires the R5 joint owner restart extension".to_vec());
        return 1;
    }
    // std::ofstream(path, binary | trunc)
    let p = crate::core::file_io::path_from_bytes(crate::core::cstdlib::c_str(path));
    let file = std::fs::OpenOptions::new().write(true).create(true).truncate(true).open(p);
    let file = match file {
        Ok(f) => f,
        Err(_) => {
            hsm::set_last_io_error(cat(&[b"Cannot open hotstart save file '", path, b"'"]));
            return 1;
        }
    };
    let mut buf: Vec<u8> = Vec::new();
    buf.extend_from_slice(b"SWMM5-HOTSTART4");
    let n_sub = ctx.n_subcatches();
    let n_land = if runoff.is_some() { ctx.n_landuses() } else { 0 };
    let n_nodes = ctx.n_nodes();
    let n_links = ctx.n_links();
    let n_pollut = ctx.n_pollutants();
    let flow_units = ctx.options.flow_units as i32;
    for v in [n_sub, n_land, n_nodes, n_links, n_pollut, flow_units] {
        buf.extend_from_slice(&v.to_ne_bytes());
    }
    if let Some(r) = runoff {
        if r.len() != hsm::legacy_runoff_size(ctx, n_land) {
            // the header is already in the stream's buffer, so it reaches the file
            use std::io::Write;
            let mut file = file;
            let _ = file.write_all(&buf);
            hsm::set_last_io_error(cat(&[b"Hotstart runoff section has the wrong size for '", path, b"'"]));
            return 1;
        }
        for &x in r.iter() {
            buf.extend_from_slice(&x.to_ne_bytes());
        }
    }
    let nodes = &ctx.nodes;
    let links = &ctx.links;
    for i in 0..n_nodes {
        let ui = i as usize;
        buf.extend_from_slice(&(nodes.depth[ui] as f32).to_ne_bytes());
        buf.extend_from_slice(&(nodes.lat_flow[ui] as f32).to_ne_bytes());
        if nodes.type_[ui] == NodeType::Storage {
            let h: f32 = if ui < nodes.hrt.len() { nodes.hrt[ui] as f32 } else { 0.0 };
            buf.extend_from_slice(&h.to_ne_bytes());
        }
        for j in 0..n_pollut {
            let qi = ui * n_pollut as usize + j as usize;
            let q: f32 = if qi < nodes.conc.len() { nodes.conc[qi] as f32 } else { 0.0 };
            buf.extend_from_slice(&q.to_ne_bytes());
        }
    }
    for i in 0..n_links {
        let ui = i as usize;
        buf.extend_from_slice(&(links.flow[ui] as f32).to_ne_bytes());
        buf.extend_from_slice(&(links.depth[ui] as f32).to_ne_bytes());
        buf.extend_from_slice(&(links.setting[ui] as f32).to_ne_bytes());
        for j in 0..n_pollut {
            let qi = ui * n_pollut as usize + j as usize;
            let q: f32 = if qi < links.conc.len() { links.conc[qi] as f32 } else { 0.0 };
            buf.extend_from_slice(&q.to_ne_bytes());
        }
    }
    let rx = &ctx.reactions;
    let nm = if rx.configured { rx.n_species() } else { 0 };
    let unm = nm.max(0) as usize;
    if nm > 0 && rx.msx_node_conc.len() == n_nodes as usize * unm && rx.msx_link_conc.len() == n_links as usize * unm {
        buf.extend_from_slice(b"OPENSWMM-HSMSX1");
        buf.extend_from_slice(&nm.to_ne_bytes());
        for m in 0..unm {
            let s = &rx.species_name[m];
            let len = (s.len() + 1) as u32;
            buf.extend_from_slice(&len.to_ne_bytes());
            buf.extend_from_slice(s);
            buf.push(0);
        }
        buf.extend_from_slice(&n_nodes.to_ne_bytes());
        buf.extend_from_slice(&n_links.to_ne_bytes());
        for &x in &rx.msx_node_conc {
            buf.extend_from_slice(&x.to_ne_bytes());
        }
        for &x in &rx.msx_link_conc {
            buf.extend_from_slice(&x.to_ne_bytes());
        }
        let s = &rx.surface;
        let surf = s.active() && s.n_species == nm && s.n_subcatch == n_sub;
        buf.extend_from_slice(&(if surf { 1i32 } else { 0 }).to_ne_bytes());
        if surf {
            buf.extend_from_slice(&s.n_subcatch.to_ne_bytes());
            buf.extend_from_slice(&s.n_landuses.to_ne_bytes());
            for v in [&s.washoff_conc, &s.ponded_qual, &s.buildup] {
                for &x in v.iter() {
                    buf.extend_from_slice(&x.to_ne_bytes());
                }
            }
        }
    }
    use std::io::Write;
    let mut file = file;
    if file.write_all(&buf).and_then(|_| file.flush()).is_err() {
        hsm::set_last_io_error(cat(&[b"Hotstart body write failed for '", path, b"'"]));
        return 1;
    }
    hsm::set_last_io_error(Vec::new());
    0
}

/// InletSolver::gatherStats: publish the inlet statistics for the report.
fn inlet_gather_stats(inlet: &crate::hydraulics::inlet::InletSolver, ctx: &mut SimulationContext) {
    let soa = inlet.soa();
    let n = std::cmp::min(soa.count, ctx.inlet_usages.count());
    if n <= 0 {
        return;
    }
    let cnt = ctx.inlet_usages.count();
    ctx.inlet_usages.resize_stats(cnt);
    ctx.inlet_diag.resize(n);
    for i in 0..n as usize {
        ctx.inlet_usages.stat_capture_vol[i] = soa.stat_capture_vol[i];
        ctx.inlet_usages.stat_bypass_vol[i] = soa.stat_bypass_vol[i];
        ctx.inlet_usages.stat_backflow_vol[i] = soa.stat_backflow_vol[i];
        ctx.inlet_usages.stat_peak_flow[i] = soa.stat_peak_flow[i];
        let dg = &mut ctx.inlet_diag;
        dg.host_node[i] = soa.host_node[i];
        dg.up_link[i] = soa.up_link[i];
        dg.is_sag[i] = soa.is_sag[i];
        dg.num_inlets[i] = soa.num_inlets[i];
        dg.flow_periods[i] = soa.flow_periods[i];
        dg.capture_periods[i] = soa.capture_periods[i];
        dg.backflow_periods[i] = soa.backflow_periods[i];
        dg.peak_flow[i] = soa.peak_flow[i];
        dg.peak_flow_capture[i] = soa.peak_flow_capture[i];
        dg.avg_flow_capture[i] = soa.avg_flow_capture[i];
        dg.bypass_freq[i] = soa.bypass_freq[i];
    }
}

/// quality::summarizeNegativeSourceClamps: the end-of-run extraction-clamp warning.
fn summarize_negative_source_clamps(ctx: &mut SimulationContext) {
    let st = &ctx.negsrc;
    if st.clamp_events == 0 && st.age_clamp_events == 0 {
        return;
    }
    let mut m = b"D-NS1 summary: extraction exceeded available mass ".to_vec();
    m.extend_from_slice(st.clamp_events.to_string().as_bytes());
    m.extend_from_slice(b" time(s) (unmet extraction ");
    m.extend_from_slice(&crate::core::cprintf::to_string_f64(st.shortfall_mass));
    m.extend_from_slice(b" internal mass units)");
    if st.age_clamp_events > 0 {
        m.extend_from_slice(b", water-age extraction clamped ");
        m.extend_from_slice(st.age_clamp_events.to_string().as_bytes());
        m.extend_from_slice(b" time(s)");
    }
    m.extend_from_slice(b", first at element index ");
    m.extend_from_slice(st.first_node.to_string().as_bytes());
    m.extend_from_slice(
        b". Extraction clamps ROUTINELY while the system fills, when elements are near-empty; a nonzero count is not by itself a modelling error. Node-seam clamps are ledgered (the ledger carries the mass actually removed); ARD cell-source clamps are counted only \xE2\x80\x94 see bookNegativeCellSourceClamp.",
    );
    ctx.warnings.push(m);
}

/// ueb::bulkTemperature (hydrology/surface/UebSnow.cpp): the snow/active-soil bulk
/// temperature of a point state.
fn ueb_bulk_temperature(s: &crate::hydrology::ueb_types::State, p: &crate::hydrology::ueb_types::Parameters) -> f64 {
    const WATER_DENSITY: f64 = 1000.0;
    const LATENT_HEAT: f64 = 333.5;
    const ICE_HEAT: f64 = 2.09;
    let capacity = WATER_DENSITY * s.water * ICE_HEAT + p.soil_density * p.active_soil_depth * p.soil_heat_capacity;
    if s.energy <= 0.0 {
        s.energy / capacity
    } else {
        crate::core::cmath::cmax(0.0, s.energy - s.water * WATER_DENSITY * LATENT_HEAT) / capacity
    }
}

/// surface::nativeSnowResults (hydrology/surface/NativeSnow.cpp): one subcatchment's ten
/// snow-process columns, NaN where it has no UEB column.
fn native_snow_results(c: &SimulationContext, s: &crate::hydrology::snow::SnowSoA, sc: i32) -> [f64; 10] {
    let mut out = [f64::NAN; 10];
    let rt = &c.snowpacks.runtime;
    if sc < 0 || sc >= s.n_subcatch || !s.external[sc as usize] || (sc + 1) * 3 > rt.columns.len() as i32 {
        return out;
    }
    out = [0.0; 10];
    let mut advanced = false;
    for k in 0..3 {
        let i = (sc * 3 + k) as usize;
        let a = s.fArea[i];
        if a <= 0.0 {
            continue;
        }
        let l = &rt.columns[i].accepted;
        let p = &l.point;
        let f = &rt.last[i];
        out[9] = crate::core::cmath::cmax(out[9], rt.sample_end[i]);
        advanced = advanced || rt.sample_seconds[i] > 0.0;
        out[0] += a * l.water;
        out[1] += a * (if l.cover > 0.0 { l.cover } else { 1.0 }) * p.energy;
        let bulk = ueb_bulk_temperature(p, &c.snowpacks.config.parameters);
        out[2] += a * (if p.history.is_empty() { crate::core::cmath::cmin(0.0, bulk) } else { p.history.last().unwrap().surface });
        out[3] += a * bulk;
        out[4] += a * f.point.albedo;
        out[5] += a * l.cover;
        out[6] += a * s.imelt[i] * 0.3048;
        out[7] += a * s.external_net[i] * 0.3048;
        out[8] += a * (if rt.sample_seconds[i] > 0.0 { f.sublimation / rt.sample_seconds[i] } else { 0.0 });
    }
    if !advanced {
        out[4] = f64::NAN;
    }
    out
}

/// SWMM_SnowLedger
#[derive(Debug, Clone, Copy, Default)]
struct SnowLedger {
    elapsed_s: f64,
    initial_m3: f64,
    precipitation_m3: f64,
    liquid_m3: f64,
    throughfall_m3: f64,
    sublimation_m3: f64,
    removed_m3: f64,
    transferred_m3: f64,
    stored_m3: f64,
    residual_m3: f64,
    initial_energy_kj: f64,
    net_energy_kj: f64,
    geometry_energy_kj: f64,
    plowed_energy_kj: f64,
    stored_energy_kj: f64,
    energy_residual_kj: f64,
}

/// swmm_snow_get_ueb_ledger: the run's UEB water and energy balance (None: SWMM_ERR_LIFECYCLE,
/// no UEB runtime).
fn ueb_ledger(ctx: &SimulationContext, soa: &crate::hydrology::snow::SnowSoA) -> Option<SnowLedger> {
    if !ctx.snowpacks.runtime.active {
        return None;
    }
    let mut out = SnowLedger::default();
    let c = ctx;
    let ordinary = true;
    let runtime = &c.snowpacks.runtime;
    out.elapsed_s = crate::core::cmath::cmax(out.elapsed_s, runtime.completed_seconds);
    out.initial_m3 += runtime.initial_water;
    out.precipitation_m3 += runtime.precipitation;
    out.liquid_m3 += runtime.liquid;
    out.throughfall_m3 += runtime.throughfall;
    out.sublimation_m3 += runtime.sublimation;
    out.removed_m3 += runtime.removed;
    out.transferred_m3 += runtime.transferred_water;
    out.initial_energy_kj += runtime.initial_energy;
    out.net_energy_kj += runtime.net_energy;
    out.geometry_energy_kj += runtime.geometry_energy;
    out.plowed_energy_kj += runtime.plowed_energy;
    for sc in 0..soa.n_subcatch as usize {
        if soa.external[sc] {
            let lid = if sc < c.subcatches.total_lid_area_ft2.len() { c.subcatches.total_lid_area_ft2[sc] } else { 0.0 };
            let area = crate::core::cmath::cmax(0.0, c.subcatches.area[sc] / ucf::ucf(ucf::LANDAREA, &c.options) - lid) * 0.3048 * 0.3048;
            let mut stored = 0.0;
            let mut energy = 0.0;
            for k in 0..3 {
                let l = &runtime.columns[sc * 3 + k].accepted;
                let a = area * soa.fArea[sc * 3 + k];
                if !ordinary || soa.enabled[sc] != 0 {
                    out.elapsed_s = crate::core::cmath::cmax(out.elapsed_s, runtime.sample_end[sc * 3 + k]);
                }
                stored += a * l.water;
                energy += a * (if l.cover > 0.0 { l.cover } else { 1.0 }) * l.point.energy;
            }
            if ordinary && soa.enabled[sc] == 0 {
                out.initial_m3 -= stored;
                out.initial_energy_kj -= energy;
            } else {
                out.stored_m3 += stored;
                out.stored_energy_kj += energy;
            }
        }
    }
    out.energy_residual_kj = out.initial_energy_kj + out.net_energy_kj + out.geometry_energy_kj - out.plowed_energy_kj - out.stored_energy_kj;
    out.residual_m3 = out.initial_m3 + out.precipitation_m3 + out.transferred_m3 - out.liquid_m3 - out.throughfall_m3 - out.sublimation_m3 - out.removed_m3 - out.stored_m3;
    Some(out)
}

impl crate::core::interface_file::InterfaceManager {
    /// InterfaceManager::writeOutfallResults: one row per outlet node to the SAVE OUTFLOWS
    /// file at a report instant.
    pub fn write_outfall_results(&mut self, ctx: &SimulationContext, current_time: f64) {
        let f = match self.outfile_.as_mut() {
            Some(f) => f,
            None => return,
        };
        let (yr, mon, day) = datetime::decode_date(current_time);
        let (hr, mn, sec) = datetime::decode_time(current_time);
        let mut out: Vec<u8> = Vec::new();
        for i in 0..ctx.n_nodes() {
            let ui = i as usize;
            if !Self::is_outlet_node(ctx, i) {
                continue;
            }
            let name = &ctx.node_names.names()[ui];
            out.extend(crate::cformat!(b"\n%-16s", name));
            out.extend(crate::cformat!(b" %04d %02d  %02d  %02d  %02d  %02d ", yr, mon, day, hr, mn, sec));
            let mut fu = ctx.options.flow_units as i32;
            if !(0..=5).contains(&fu) {
                fu = 0;
            }
            out.extend(crate::cformat!(b" %-10f", ctx.nodes.inflow[ui] * ucf::QCF[fu as usize]));
            let np = ctx.n_pollutants();
            for p in 0..np {
                let qi = ui * np as usize + p as usize;
                let conc = if qi < ctx.nodes.conc.len() { ctx.nodes.conc[qi] } else { 0.0 };
                out.extend(crate::cformat!(b" %-10f", conc));
            }
        }
        use std::io::Write;
        let _ = f.write_all(&out);
    }
}

impl SwmmEngine {
    /// SWMMEngine::postOutputSnapshot: build and write the report period's snapshot when
    /// output is due. noexcept upstream; Err only for a call into untranslated code.
    pub(crate) fn post_output_snapshot(&mut self, _dt_step: f64) -> Result<(), Terminated> {
        if !tc::output_due(&self.ctx_) {
            return Ok(());
        }
        let report_date = datetime::add_seconds(self.ctx_.options.start_date, (self.ctx_.next_report_ms + 1.0) / 1000.0);

        // REPORT_START gate: no bytes before ReportStart, but the grid advances
        if report_date < self.ctx_.options.report_start {
            tc::reset_output_timer(&mut self.ctx_);
            return Ok(());
        }

        // legacy stats_updateMaxNodeDepth: the largest REPORTED depth of every node
        {
            let c = &mut self.ctx_;
            let span = c.elapsed_ms - c.old_elapsed_ms;
            let f = if span > 0.0 { (c.next_report_ms - c.old_elapsed_ms) / span } else { 1.0 };
            let f1 = 1.0 - f;
            let ucf_len = ucf::ucf(ucf::LENGTH, &c.options);
            for i in 0..c.n_nodes() as usize {
                let z = ((f1 * c.nodes.old_depth[i] + f * c.nodes.depth[i]) * ucf_len) as f32;
                if z as f64 > c.nodes.stat_max_rpt_depth[i] {
                    c.nodes.stat_max_rpt_depth[i] = z as f64;
                }
            }
        }

        // Routing interface file: one outfall row per reporting step
        if self.do_routing_ {
            self.iface_.write_outfall_results(&self.ctx_, report_date);
        }

        if self.save_results_ != 0 && !self.plugins_.is_empty() {
            let mut snap = SimulationSnapshot { sim_time: report_date, ..Default::default() };
            snap.node_count = self.ctx_.n_nodes();
            snap.link_count = self.ctx_.n_links();
            snap.subcatch_count = self.ctx_.n_subcatches();
            snap.gage_count = self.ctx_.n_gages();
            snap.pollut_count = self.ctx_.n_reported_species();
            snap.flow_units_code = self.ctx_.options.flow_units as i32;

            let route_span = self.ctx_.elapsed_ms - self.ctx_.old_elapsed_ms;
            let f_rt = if route_span > 0.0 { (self.ctx_.next_report_ms - self.ctx_.old_elapsed_ms) / route_span } else { 1.0 };
            let f1_rt = 1.0 - f_rt;

            // Node results (legacy node_getResults)
            {
                let n_n = self.ctx_.n_nodes() as usize;
                snap.nodes.depth.resize(n_n, 0.0);
                snap.nodes.head.resize(n_n, 0.0);
                snap.nodes.volume.resize(n_n, 0.0);
                snap.nodes.lateral_inflow.resize(n_n, 0.0);
                snap.nodes.total_inflow.resize(n_n, 0.0);
                for ui in 0..n_n {
                    let i = ui as i32;
                    let c = &self.ctx_;
                    let y = f1_rt * c.nodes.old_depth[ui] + f_rt * c.nodes.depth[ui];
                    snap.nodes.depth[ui] = y;
                    snap.nodes.head[ui] = if c.options.report_signed_heads != 0 { c.nodes.invert_elev[ui] + f1_rt * c.nodes.old_depth[ui] + f_rt * (c.nodes.head[ui] - c.nodes.invert_elev[ui]) } else { c.nodes.invert_elev[ui] + y };
                    snap.nodes.volume[ui] = f1_rt * self.reported_node_volume_at(i, c.nodes.old_depth[ui], c.nodes.old_volume[ui]) + f_rt * self.reported_node_volume_at(i, c.nodes.depth[ui], c.nodes.volume[ui]);
                    snap.nodes.lateral_inflow[ui] = f1_rt * c.nodes.old_lat_flow[ui] + f_rt * c.nodes.lat_flow[ui];
                    snap.nodes.total_inflow[ui] = f1_rt * c.nodes.old_inflow[ui] + f_rt * c.nodes.inflow[ui];
                }
            }
            snap.nodes.overflow = self.ctx_.nodes.overflow.clone();

            // Link results (legacy link_getResults)
            {
                let n_l = self.ctx_.n_links() as usize;
                snap.links.flow.resize(n_l, 0.0);
                snap.links.depth.resize(n_l, 0.0);
                snap.links.volume.resize(n_l, 0.0);
                snap.links.velocity.resize(n_l, 0.0);
                snap.links.capacity.resize(n_l, 0.0);
                self.ensure_xsp_cache();
                let c = &self.ctx_;
                for uj in 0..n_l {
                    let j = uj as i32;
                    let y = f1_rt * c.links.old_depth[uj] + f_rt * c.links.depth[uj];
                    let mut q = f1_rt * c.links.old_flow[uj] + f_rt * c.links.flow[uj];
                    let v = f1_rt * c.links.old_volume[uj] + f_rt * c.links.volume[uj];
                    let mut veloc = 0.0;
                    let cap;
                    let lt = c.links.type_[uj];
                    if lt == LinkType::Conduit {
                        let xs = link::build_xsect_params(&c.links, uj, Some(&c.transect_tables));
                        let cr = c.link_subtypes.conduit_row(j);
                        let nb = if cr >= 0 { c.link_subtypes.conduits.barrels[cr as usize] } else { 1 };
                        veloc = link::get_velocity(&xs, q, y, nb);
                        cap = link::get_capacity(&xs, y);
                    } else {
                        cap = c.links.setting[uj];
                    }
                    if lt == LinkType::Pump && c.links.old_flow[uj] * c.links.flow[uj] == 0.0 {
                        q = if f_rt >= f1_rt { c.links.flow[uj] } else { c.links.old_flow[uj] };
                    }
                    let dir = c.links.direction[uj] as f64;
                    snap.links.depth[uj] = y;
                    snap.links.flow[uj] = q * dir;
                    snap.links.velocity[uj] = veloc * dir;
                    snap.links.volume[uj] = v;
                    snap.links.capacity[uj] = cap;
                }
            }

            // Subcatchment state
            snap.subcatch.rainfall = self.ctx_.subcatches.rainfall.clone();
            snap.subcatch.evap = self.ctx_.subcatches.evap_loss.clone();
            snap.subcatch.infil = self.ctx_.subcatches.infil_loss.clone();
            snap.subcatch.runoff = self.ctx_.subcatches.runoff.clone();
            snap.subcatch.gw_flow = self.ctx_.subcatches.gw_flow.clone();

            // Legacy-parity reported losses on LID subcatchments (subcatch.c:740-744)
            {
                let c = &self.ctx_;
                let land2ft2 = 1.0 / ucf::ucf(ucf::LANDAREA, &c.options);
                let dt_last = (self.new_runoff_ms_ - self.old_runoff_ms_) / 1000.0;
                let rsoa = self.runoff_.soa();
                let mut lid_evap_cfs: Vec<f64> = Vec::new();
                let mut lid_infil_cfs: Vec<f64> = Vec::new();
                for t in 0..self.lid_.num_groups() {
                    let g = self.lid_.group(t);
                    for uu in 0..g.count as usize {
                        let sc = g.subcatch_idx[uu];
                        if sc < 0 || sc >= c.n_subcatches() || dt_last <= 0.0 {
                            continue;
                        }
                        if lid_evap_cfs.is_empty() {
                            let un = c.n_subcatches() as usize;
                            lid_evap_cfs.clear();
                            lid_evap_cfs.resize(un, 0.0);
                            lid_infil_cfs.clear();
                            lid_infil_cfs.resize(un, 0.0);
                        }
                        let usc = sc as usize;
                        lid_evap_cfs[usc] += g.evap_loss[uu] * g.area[uu] / dt_last;
                        lid_infil_cfs[usc] += g.infil_loss[uu] * g.area[uu] / dt_last;
                    }
                }
                for ui in 0..c.n_subcatches() as usize {
                    if c.subcatches.total_lid_area_ft2[ui] <= 0.0 {
                        continue;
                    }
                    let full_ft2 = c.subcatches.area[ui] * land2ft2;
                    if full_ft2 <= 0.0 {
                        continue;
                    }
                    let nonlid_ft2 = rsoa.area[ui];
                    let ev = if lid_evap_cfs.is_empty() { 0.0 } else { lid_evap_cfs[ui] };
                    let inn = if lid_infil_cfs.is_empty() { 0.0 } else { lid_infil_cfs[ui] };
                    snap.subcatch.evap[ui] = (c.subcatches.evap_loss[ui] * nonlid_ft2 + ev) / full_ft2;
                    snap.subcatch.infil[ui] = (c.subcatches.infil_loss[ui] * nonlid_ft2 + inn) / full_ft2;
                }
            }

            // Snow-process columns of a UEB-snow deck, on the snow's own accepted clock
            if crate::hydrology::native_snow::has_ueb_snow(&self.ctx_) {
                let n_s = self.ctx_.n_subcatches() as usize;
                snap.subcatch.snow_process.clear();
                snap.subcatch.snow_process.resize(n_s * crate::core::output_plugin::SNOW_PROCESS_COUNT, 0.0);
                let snow = self.snow_.state();
                for sc in 0..n_s {
                    let row = native_snow_results(&self.ctx_, snow, sc as i32);
                    snap.subcatch.snow_process[sc * row.len()..sc * row.len() + row.len()].copy_from_slice(&row);
                }
            }

            // Reported runoff, time-interpolated on the runoff clock (subcatch.c:865)
            {
                const MIN_RUNOFF: f64 = 2.31481e-8;
                let c = &self.ctx_;
                let span = self.new_runoff_ms_ - self.old_runoff_ms_;
                let f = if span > 0.0 { (c.next_report_ms - self.old_runoff_ms_) / span } else { 1.0 };
                let f1 = 1.0 - f;
                let land2ft2 = 1.0 / ucf::ucf(ucf::LANDAREA, &c.options);
                let gw = self.groundwater_.state();
                for ui in 0..c.n_subcatches() as usize {
                    let mut ro = f1 * c.subcatches.old_runoff[ui] + f * c.subcatches.runoff[ui];
                    if c.subcatches.total_lid_area_ft2[ui] > 0.0 {
                        ro += f1 * c.subcatches.old_lid_drain_flow[ui] + f * c.subcatches.lid_drain_flow[ui];
                    }
                    if ro < MIN_RUNOFF * (c.subcatches.area[ui] * land2ft2) {
                        ro = 0.0;
                    }
                    snap.subcatch.runoff[ui] = ro;
                    let a_ft2 = c.subcatches.area[ui] / ucf::ucf(ucf::LANDAREA, &c.options);
                    let r_new = if ui < gw.gw_flow.len() { gw.gw_flow[ui] } else { 0.0 };
                    let r_old = if ui < self.old_gw_rate_.len() { self.old_gw_rate_[ui] } else { 0.0 };
                    snap.subcatch.gw_flow[ui] = (f1 * r_old + f * r_new) * a_ft2;
                }
            }

            // GW elevation and soil moisture (legacy subcatch_getResults)
            {
                let c = &self.ctx_;
                let n_s = c.n_subcatches() as usize;
                snap.subcatch.gw_elev.resize(n_s, 0.0);
                snap.subcatch.soil_moist.resize(n_s, 0.0);
                let gw = self.groundwater_.state();
                for us in 0..n_s {
                    let aq = c.subcatches.gw_aquifer[us];
                    if aq >= 0 {
                        snap.subcatch.gw_elev[us] = gw.bottom_elev[us] + gw.lower_depth[us];
                        snap.subcatch.soil_moist[us] = gw.theta[us];
                    }
                }
            }

            // Per-subcatch snow depth on the runoff-clock weight (subcatch.c:888-890)
            {
                let c = &self.ctx_;
                let n_s = c.n_subcatches() as usize;
                snap.subcatch.snow_depth.resize(n_s, 0.0);
                let span = self.new_runoff_ms_ - self.old_runoff_ms_;
                let f = if span > 0.0 { (c.next_report_ms - self.old_runoff_ms_) / span } else { 1.0 };
                let f1 = 1.0 - f;
                for us in 0..n_s {
                    snap.subcatch.snow_depth[us] = f1 * c.subcatches.old_snow_depth[us] + f * c.subcatches.snow_depth[us];
                }
            }

            // System-level results
            let has_subcatchments = self.ctx_.n_subcatches() > 0;
            snap.has_subcatchments = has_subcatchments;
            snap.sys_temperature = if has_subcatchments { self.ctx_.climate_state.temperature } else { 0.0 };
            snap.sys_rainfall = 0.0;
            snap.sys_snow_depth = 0.0;
            snap.sys_evap = 0.0;
            snap.sys_infil = 0.0;
            snap.sys_runoff = 0.0;
            snap.sys_pet = if has_subcatchments { self.ctx_.climate_state.evap_rate } else { 0.0 };
            {
                let mb = &self.ctx_.mass_balance;
                snap.sys_dw_inflow = mb.step_dw_inflow;
                snap.sys_gw_inflow = mb.step_gw_inflow;
                snap.sys_ii_inflow = mb.step_rdii_inflow;
                snap.sys_ext_inflow = mb.step_ext_inflow;
                snap.sys_flooding = mb.step_flooding;
                snap.sys_outflow = mb.step_outflow;
            }

            // Total storage: the sum of the REPORTED volumes, in REAL4 display units
            {
                let vol_ucf = ucf::UCF[ucf::VOLUME as usize][ucf::get_unit_system(self.ctx_.options.flow_units as i32) as usize];
                let mut tot_store = 0.0f32;
                if self.ctx_.options.rpt_averages {
                    for j in 0..self.ctx_.n_nodes() {
                        tot_store += (self.reported_node_volume(j) * vol_ucf) as f32;
                    }
                    for uj in 0..self.ctx_.links.volume.len() {
                        tot_store += (self.ctx_.links.volume[uj] * vol_ucf) as f32;
                    }
                } else {
                    for uj in 0..snap.nodes.volume.len() {
                        tot_store += (snap.nodes.volume[uj] * vol_ucf) as f32;
                    }
                    for uj in 0..snap.links.volume.len() {
                        tot_store += (snap.links.volume[uj] * vol_ucf) as f32;
                    }
                }
                snap.sys_storage = tot_store as f64 / vol_ucf;
            }

            // OPENSWMM_HAS_2D: close a pending (COUPLING_SYNC-batched) 2D window and bring the
            // render fields up to date before the snapshot (the checkpoint profiles are off
            // without engine checkpointing). postOutputSnapshot is noexcept: an exception from
            // the batch terminates (a C++ throw becomes a termination here).
            #[cfg(feature = "twod")]
            if self.surface_router_.is_active() {
                let r = self.surface_router_.flush_pending_batch(&mut self.ctx_);
                if let Err(t) = r {
                    return Err(match t.exception.strip_prefix(crate::core::engine_step::THROWN_PREFIX) {
                        Some(ty) => Terminated { exception: ty, what: t.what },
                        None => t,
                    });
                }
                self.surface_router_.refresh_render_fields_if_stale();
            }

            // fillSurfaceSnapshot(snap): the 2D surface block (core/engine_output_twod.rs)
            #[cfg(feature = "twod")]
            self.fill_surface_snapshot(&mut snap);
            #[cfg(all(feature = "twod", feature = "gate"))]
            crate::core::engine_output_twod::gate_snapshot_hook(&snap);

            let du = ucf::DisplayUnits::from(&self.ctx_.options);
            convert_snapshot_to_display(&mut snap, &du);

            if self.ctx_.options.rpt_averages {
                self.apply_avg_results(&mut snap);
            }

            // legacy node_getResults NODE_HEAD: f32(depth) + f32(invert), in float
            if !self.ctx_.options.rpt_averages && self.ctx_.options.report_signed_heads == 0 {
                let nn = std::cmp::min(snap.nodes.head.len(), snap.nodes.depth.len());
                for i in 0..nn {
                    let depth_f32 = snap.nodes.depth[i] as f32;
                    let invert_f32 = (self.ctx_.nodes.invert_elev[i] * du.length) as f32;
                    snap.nodes.head[i] = (depth_f32 + invert_f32) as f64;
                }
            }

            // Subcatchment RAINFALL: the gage's report-time rainfall, native units
            {
                let c = &self.ctx_;
                let n_g = c.n_gages();
                let mut gage_report = vec![0.0f64; n_g as usize];
                for g in 0..n_g {
                    gage_report[g as usize] = crate::hydrology::gage_report::get_report_rainfall(c, g, report_date);
                }
                let n_s = snap.subcatch.rainfall.len();
                for s in 0..n_s {
                    let g = c.subcatches.gage[s];
                    snap.subcatch.rainfall[s] = if g >= 0 && g < n_g { gage_report[g as usize] } else { 0.0 };
                }
            }

            // Subcatchment-derived system results, from the REPORTED float32 values
            if has_subcatchments {
                let c = &self.ctx_;
                let mut total_area = 0.0f32;
                let (mut a_rain, mut a_snow, mut a_evap) = (0.0f32, 0.0f32, 0.0f32);
                let (mut a_infil, mut a_runoff) = (0.0f32, 0.0f32);
                let gw = self.groundwater_.state();
                let ucf_land = ucf::ucf(ucf::LANDAREA, &c.options);
                for ui in 0..c.n_subcatches() as usize {
                    let area = (c.subcatches.area[ui] / ucf_land) * ucf_land;
                    total_area += area as f32;
                    a_rain += ((snap.subcatch.rainfall[ui] as f32) as f64 * area) as f32;
                    a_snow += ((snap.subcatch.snow_depth[ui] as f32) as f64 * area) as f32;
                    a_evap += ((snap.subcatch.evap[ui] as f32) as f64 * area) as f32;
                    if c.subcatches.gw_aquifer[ui] >= 0 && ui < gw.upper_evap.len() {
                        a_evap += ((gw.upper_evap[ui] + gw.lower_evap[ui]) * du.evaprate * area) as f32;
                    }
                    a_infil += ((snap.subcatch.infil[ui] as f32) as f64 * area) as f32;
                    a_runoff += snap.subcatch.runoff[ui] as f32;
                }
                if total_area > 0.0 {
                    a_evap /= total_area;
                    a_rain /= total_area;
                    a_snow /= total_area;
                    a_infil /= total_area;
                }
                snap.sys_rainfall = a_rain as f64;
                snap.sys_snow_depth = a_snow as f64;
                snap.sys_evap = a_evap as f64;
                snap.sys_infil = a_infil as f64;
                snap.sys_runoff = a_runoff as f64;
            }

            // Pollutant concentrations
            if !self.ctx_.options.ignore_quality && self.ctx_.n_reported_species() > 0 {
                let c = &self.ctx_;
                let np_s = c.n_pollutants() as usize;
                let nr_s = c.n_reported_species() as usize;
                let n_n = c.n_nodes() as usize;
                let n_l = c.n_links() as usize;
                let n_s = c.n_subcatches() as usize;
                let age_col = c.options.water_age && nr_s > np_s;
                let temp_col = c.options.heat_transport;
                let temp_i = np_s + if c.options.water_age { 1 } else { 0 };
                let rx = &c.reactions;
                let unm = if rx.n_species() > 0 { rx.n_species() as usize } else { 0 };
                let msx_i0 = temp_i + if c.options.heat_transport { 1 } else { 0 };
                const K_SEC_PER_HOUR: f64 = 3600.0;
                const K_DRY_REPORT_DEPTH: f64 = 1.0e-9;

                snap.node_quality.clear();
                snap.node_quality.resize(n_n * nr_s, 0.0);
                for n in 0..n_n {
                    for p in 0..np_s {
                        let src = n * np_s + p;
                        if src < c.nodes.conc.len() && src < c.nodes.conc_old.len() {
                            snap.node_quality[n * nr_s + p] = f1_rt * c.nodes.conc_old[src] + f_rt * c.nodes.conc[src];
                        }
                    }
                    if age_col && n < c.water_age_state.node_age.len() {
                        let wet = n < snap.nodes.depth.len() && snap.nodes.depth[n] > K_DRY_REPORT_DEPTH;
                        snap.node_quality[n * nr_s + np_s] = if wet { c.water_age_state.node_age[n] / K_SEC_PER_HOUR } else { 0.0 };
                    }
                    if temp_col && n < c.heat_state.node_temp.len() {
                        snap.node_quality[n * nr_s + temp_i] = c.heat_state.node_temp[n];
                    }
                    for r in 0..rx.reported_species.len() {
                        let m = rx.reported_species[r] as usize;
                        let src = n * unm + m;
                        if src < rx.msx_node_conc.len() {
                            snap.node_quality[n * nr_s + msx_i0 + r] = rx.msx_node_conc[src];
                        }
                    }
                }

                snap.link_quality.clear();
                snap.link_quality.resize(n_l * nr_s, 0.0);
                for l in 0..n_l {
                    for p in 0..np_s {
                        let src = l * np_s + p;
                        if src < c.links.conc.len() && src < c.links.conc_old.len() {
                            snap.link_quality[l * nr_s + p] = f1_rt * c.links.conc_old[src] + f_rt * c.links.conc[src];
                        }
                    }
                    if age_col && l < c.water_age_state.link_age.len() {
                        let holds = l < c.links.volume.len() && c.links.volume[l] > crate::quality::quality_routing::ZERO_VOLUME;
                        let conveys = l < c.links.flow.len() && c.links.flow[l].abs() > crate::core::constants::TINY;
                        let wet = holds || conveys;
                        snap.link_quality[l * nr_s + np_s] = if wet { c.water_age_state.link_age[l] / K_SEC_PER_HOUR } else { 0.0 };
                    }
                    if temp_col && l < c.heat_state.link_temp.len() {
                        snap.link_quality[l * nr_s + temp_i] = c.heat_state.link_temp[l];
                    }
                    for r in 0..rx.reported_species.len() {
                        let m = rx.reported_species[r] as usize;
                        let src = l * unm + m;
                        if src < rx.msx_link_conc.len() {
                            snap.link_quality[l * nr_s + msx_i0 + r] = rx.msx_link_conc[src];
                        }
                    }
                }

                snap.subcatch_quality.clear();
                snap.subcatch_quality.resize(n_s * nr_s, 0.0);
                let ro_span_q = self.new_runoff_ms_ - self.old_runoff_ms_;
                let f_ro_q = if ro_span_q > 0.0 { (c.next_report_ms - self.old_runoff_ms_) / ro_span_q } else { 1.0 };
                let f1_ro_q = 1.0 - f_ro_q;
                for s in 0..n_s {
                    if temp_col && s < c.heat_state.subcatch_runoff_temp.len() {
                        snap.subcatch_quality[s * nr_s + temp_i] = c.heat_state.subcatch_runoff_temp[s];
                    }
                    let has_runoff = s < snap.subcatch.runoff.len() && snap.subcatch.runoff[s] != 0.0;
                    if !has_runoff {
                        continue;
                    }
                    for p in 0..np_s {
                        let src = s * np_s + p;
                        if src < c.subcatches.conc.len() && src < c.subcatches.conc_old.len() {
                            snap.subcatch_quality[s * nr_s + p] = f1_ro_q * c.subcatches.conc_old[src] + f_ro_q * c.subcatches.conc[src];
                        }
                    }
                    for r in 0..rx.reported_species.len() {
                        let m = rx.reported_species[r] as usize;
                        let src = s * unm + m;
                        if src < rx.surface.washoff_conc.len() && src < rx.surface.washoff_conc_old.len() {
                            snap.subcatch_quality[s * nr_s + msx_i0 + r] = f1_ro_q * rx.surface.washoff_conc_old[src] + f_ro_q * rx.surface.washoff_conc[src];
                        }
                    }
                    if age_col && s < c.water_age_state.subcatch_runoff_age.len() {
                        snap.subcatch_quality[s * nr_s + np_s] = c.water_age_state.subcatch_runoff_age[s] / K_SEC_PER_HOUR;
                    }
                }
            }

            // io_thread_.post(snap): the plugins' update, in order
            let _ = self.plugins_.update_all(&snap);
        }
        tc::reset_output_timer(&mut self.ctx_);
        Ok(())
    }

    /// SWMMEngine::accumulateAvgResults: one routing step into the float32 display-unit
    /// accumulators (legacy output_updateAvgResults).
    pub(crate) fn accumulate_avg_results(&mut self) -> Result<(), Terminated> {
        let du = ucf::DisplayUnits::from(&self.ctx_.options);
        for j in 0..self.ctx_.n_nodes() {
            let uj = j as usize;
            let c = &self.ctx_;
            let dd = (c.nodes.depth[uj] * du.length) as f32;
            let rv = self.reported_node_volume(j);
            let a = &mut self.avg_;
            a.node_depth[uj] += dd;
            a.node_head[uj] += dd + (c.nodes.invert_elev[uj] * du.length) as f32;
            a.node_volume[uj] += (rv * du.volume) as f32;
            a.node_lat_inflow[uj] += (c.nodes.lat_flow[uj] * du.flow) as f32;
            a.node_total_inflow[uj] += (c.nodes.inflow[uj] * du.flow) as f32;
            a.node_overflow[uj] += (c.nodes.overflow[uj] * du.flow) as f32;
        }
        self.ensure_xsp_cache();
        let c = &self.ctx_;
        let a = &mut self.avg_;
        for j in 0..c.n_links() {
            let uj = j as usize;
            let q = c.links.flow[uj];
            let d = c.links.depth[uj];
            let dir = c.links.direction[uj] as f64;
            a.link_flow[uj] += (q * du.flow * dir) as f32;
            a.link_depth[uj] += (d * du.length) as f32;
            a.link_volume[uj] += (c.links.volume[uj] * du.volume) as f32;
            let lt = c.links.type_[uj];
            if lt == LinkType::Conduit {
                let xs = link::build_xsect_params(&c.links, uj, Some(&c.transect_tables));
                let cr = c.link_subtypes.conduit_row(j);
                let nb = if cr >= 0 { c.link_subtypes.conduits.barrels[cr as usize] } else { 1 };
                a.link_velocity[uj] += (link::get_velocity(&xs, q, d, nb) * du.length * dir) as f32;
                a.link_capacity[uj] += link::get_capacity(&xs, d) as f32;
            } else {
                a.link_capacity[uj] = (c.links.setting[uj] as f32) * ((a.n_steps + 1) as f32);
            }
        }
        a.n_steps += 1;
        Ok(())
    }

    /// SWMMEngine::applyAvgResults: the averaged node/link values into a snapshot
    /// (legacy output_saveAvgResults), then a fresh accumulator.
    fn apply_avg_results(&mut self, snap: &mut SimulationSnapshot) {
        if self.avg_.n_steps <= 0 {
            return;
        }
        let n = self.avg_.n_steps as f32;
        let a = &self.avg_;
        for uj in 0..self.ctx_.n_nodes() as usize {
            snap.nodes.depth[uj] = (a.node_depth[uj] / n) as f64;
            snap.nodes.head[uj] = (a.node_head[uj] / n) as f64;
            snap.nodes.volume[uj] = (a.node_volume[uj] / n) as f64;
            snap.nodes.lateral_inflow[uj] = (a.node_lat_inflow[uj] / n) as f64;
            snap.nodes.total_inflow[uj] = (a.node_total_inflow[uj] / n) as f64;
            snap.nodes.overflow[uj] = (a.node_overflow[uj] / n) as f64;
        }
        for uj in 0..self.ctx_.n_links() as usize {
            snap.links.flow[uj] = (a.link_flow[uj] / n) as f64;
            snap.links.depth[uj] = (a.link_depth[uj] / n) as f64;
            snap.links.velocity[uj] = (a.link_velocity[uj] / n) as f64;
            snap.links.volume[uj] = (a.link_volume[uj] / n) as f64;
            snap.links.capacity[uj] = (a.link_capacity[uj] / n) as f64;
        }
        self.avg_.reset();
    }

    /// SWMMEngine::end().
    pub fn end(&mut self) -> Result<i32, Terminated> {
        if self.ctx_.state != EngineState::Running && self.ctx_.state != EngineState::Ended {
            if self.ctx_.state == EngineState::ErrorState {
                return Ok(if self.ctx_.error_code != SWMM_OK { self.ctx_.error_code } else { SWMM_ERR_WRONG_STATE });
            }
            self.set_error(SWMM_ERR_WRONG_STATE, b"swmm_engine_end: engine must be running or ended");
            return Ok(SWMM_ERR_WRONG_STATE);
        }

        // Native thermal hotstarts: the last accepted full engine checkpoint
        if self.ctx_.state == EngineState::Running && self.ctx_.options.heat_transport {
            let saves = self.ctx_.files.hotstart_saves.clone();
            for entry in &saves {
                if entry.datetime == 0.0 {
                    self.write_hotstart_save(entry)?;
                }
            }
        }

        // OPENSWMM_HAS_2D: the source histories close at the last accepted checkpoint
        #[cfg(feature = "twod")]
        self.detach_engine_source_history();

        // Final GW storage for the mass balance report
        {
            let gw = self.groundwater_.state();
            let c = &mut self.ctx_;
            c.mass_balance.gw_final_storage = 0.0;
            for ui in 0..c.n_subcatches() as usize {
                if gw.total_depth[ui] <= 0.0 {
                    continue;
                }
                let upper_d = gw.total_depth[ui] - gw.lower_depth[ui];
                let vol = gw.theta[ui] * upper_d + gw.porosity[ui] * gw.lower_depth[ui];
                let area = c.subcatches.area[ui] / ucf::ucf(ucf::LANDAREA, &c.options);
                c.mass_balance.gw_final_storage += vol * area;
            }
        }

        summarize_negative_source_clamps(&mut self.ctx_);

        if self.ctx_.routing_stats.step_intervals[0] <= 0.0 {
            self.ctx_.routing_stats.build_histogram();
        }

        // The FV solver's cumulative counters, published for the report when it ran.
        if let Some(fv) = self.router_.fv_solver_.as_ref() {
            if fv.run_stats().nsteps > 0 {
                let s = fv.run_stats();
                let rs = &mut self.ctx_.routing_stats;
                rs.fv_nsteps = s.nsteps;
                rs.fv_nflux = s.nflux;
                rs.fv_avg_h = s.avg_h;
                rs.fv_last_h = s.last_h;
                rs.fv_min_h = s.min_h;
                rs.fv_active_min = s.active_frac_min;
                rs.fv_active_mean = s.active_frac_mean;
                rs.fv_active_max = s.active_frac_max;
                rs.fv_n_tiers = s.n_tiers;
                let mut kk = 0;
                while kk < s.n_tiers && kk < 8 {
                    rs.fv_tier_cells[kk as usize] = s.tier_cells[kk as usize];
                    kk += 1;
                }
                rs.fv_dt_argmin_pressurized = s.dt_argmin_pressurized;
                rs.fv_dt_argmin_band = s.dt_argmin_band;
                rs.fv_dt_argmin_free = s.dt_argmin_free;
                rs.fv_dt_argmin_node = s.dt_argmin_node;
                rs.fv_macro_cycles = s.n_macro_cycles;
                rs.fv_macro_rejected = s.n_macro_rejected;
            }
        }

        // OPENSWMM_HAS_2D: finalize the 2D module before the max-stats finalize (its partial
        // batch flush books the last exchanges and the solver statistics the report reads)
        #[cfg(feature = "twod")]
        self.end_2d_finalize()?;

        self.ctx_.finalize_max_stats();

        // SAVE HOTSTART (end-of-run rows); a dated row the run never reached is reported
        let saves = self.ctx_.files.hotstart_saves.clone();
        for (k, entry) in saves.iter().enumerate() {
            if entry.datetime == 0.0 {
                if !self.ctx_.options.heat_transport {
                    self.write_hotstart_save(entry)?;
                }
            } else if k >= self.hotstart_save_done_.len() || self.hotstart_save_done_[k] == 0 {
                let (y, mo, d) = datetime::decode_date(entry.datetime);
                let (h, mi, sec) = datetime::decode_time(entry.datetime);
                let mut when = crate::cformat!(b"%02d/%02d/%04d %02d:%02d:%02d", mo, d, y, h, mi, sec);
                when.truncate(31);
                let p: &[u8] = if !entry.path.absolute.is_empty() { &entry.path.absolute } else { &entry.path.original };
                self.ctx_.warnings.push(cat(&[b"SAVE HOTSTART: '", p, b"' was not written; the run ended before its save time ", &when, b"."]));
            }
        }

        // io_thread_.stop(): nothing queued here

        let plugin_result = if self.plugins_.is_empty() { 0 } else { self.plugins_.finalize_all(&self.ctx_) };
        self.ctx_.state = EngineState::Ended;
        if plugin_result != 0 {
            let m = plugin_io_error(&self.plugins_, b"Cannot finish writing results");
            self.set_error(SWMM_ERR_PLUGIN, &m);
            return Ok(SWMM_ERR_PLUGIN);
        }
        Ok(SWMM_OK)
    }

    /// SWMMEngine::writeDueHotstartSaves: each dated row once, at the first step reaching it.
    pub(crate) fn write_due_hotstart_saves(&mut self) -> Result<(), Terminated> {
        let n = self.ctx_.files.hotstart_saves.len();
        for k in 0..n {
            let entry = self.ctx_.files.hotstart_saves[k].clone();
            if entry.datetime == 0.0 || k >= self.hotstart_save_done_.len() || self.hotstart_save_done_[k] != 0 {
                continue;
            }
            let due_ms = (entry.datetime - self.ctx_.options.start_date) * 86400000.0;
            if self.ctx_.elapsed_ms + 1.0 < due_ms {
                continue;
            }
            self.hotstart_save_done_[k] = 1;
            self.write_hotstart_save(&entry)?;
        }
        Ok(())
    }

    /// SWMMEngine::writeHotstartSave: one SAVE HOTSTART row (LEGACY .hsf here; NATIVE is a
    /// stub).
    fn write_hotstart_save(&mut self, entry: &HotstartSaveEntry) -> Result<(), Terminated> {
        let sp: Bytes = if !entry.path.absolute.is_empty() { entry.path.absolute.clone() } else { entry.path.original.clone() };
        if sp.is_empty() {
            return Ok(());
        }
        let lid_restart = self.ctx_.node_subtypes.storages.lid.iter().any(|c| c.control >= 0);
        #[allow(unused_mut)]
        let mut full_continuation = self.ctx_.options.heat_transport;
        #[cfg(feature = "twod")]
        {
            full_continuation |= self.hotstart_independent_clock();
        }
        let native = entry.format == HotstartFormat::Native || (entry.format == HotstartFormat::Auto && (lid_restart || full_continuation));
        let rc;
        if native {
            return Err(stub("SWMMEngine::writeHotstartSave (a NATIVE SAVE HOTSTART file: HotStartManager::save)"));
        } else {
            if lid_restart {
                self.ctx_.warnings.push(cat(&[b"SAVE HOTSTART: '", &sp, b"' is LEGACY; a legacy SWMM5 file carries no LID state."]));
            }
            let runoff = self.legacy_runoff_section();
            rc = save_legacy_routing(&sp, &self.ctx_, Some(&runoff));
        }
        if rc != 0 {
            let m = cat(&[b"SAVE HOTSTART: ", &crate::core::hot_start_manager::last_io_error()]);
            self.ctx_.warnings.push(m);
        }
        Ok(())
    }

    /// SWMMEngine::legacyRunoffSection: the legacy hot-start runoff block (saveRunoff).
    fn legacy_runoff_section(&self) -> Vec<f64> {
        let c = &self.ctx_;
        let nsc = c.n_subcatches();
        let np = c.n_pollutants();
        let nlu = c.n_landuses();
        let mut x: Vec<f64> = Vec::with_capacity(crate::core::hot_start_manager::legacy_runoff_size(c, nlu));
        let sub = &c.subcatches;
        let soa = self.runoff_.soa();
        let gw = self.groundwater_.state();
        let sn = self.snow_.state();
        for i in 0..nsc {
            let ui = i as usize;
            x.push(soa.depth_imperv0[ui]);
            x.push(soa.depth_imperv1[ui]);
            x.push(soa.depth_perv[ui]);
            x.push(sub.runoff[ui]);
            let mut model = 0;
            let mut v6 = [0.0f64; 6];
            let mut lg = [0.0f64; 6];
            self.runoff_.infil_get_state(i, &mut model, &mut v6);
            infil_to_legacy(model, &v6, &mut lg);
            x.extend_from_slice(&lg);
            if ui < sub.gw_aquifer.len() && sub.gw_aquifer[ui] >= 0 {
                x.push(gw.theta[ui]);
                x.push(gw.bottom_elev[ui] + gw.lower_depth[ui]);
                x.push(gw.gw_flow[ui]);
                x.push(gw.max_infil_vol[ui]);
            }
            if ui < sub.snowpack.len() && sub.snowpack[ui] >= 0 {
                for j in 0..crate::hydrology::snow::N_SUBAREAS as usize {
                    let pk = ui * crate::hydrology::snow::N_SUBAREAS as usize + j;
                    let has = pk < sn.wsnow.len();
                    x.push(if has { sn.wsnow[pk] } else { 0.0 });
                    x.push(if has { sn.fw[pk] } else { 0.0 });
                    x.push(if has { sn.coldc[pk] } else { 0.0 });
                    x.push(if has { sn.ati[pk] } else { 0.0 });
                    x.push(if has { sn.awe[pk] } else { 0.0 });
                }
            }
            if np > 0 {
                for p in 0..np as usize {
                    x.push(sub.conc[ui * np as usize + p]);
                }
                for p in 0..np as usize {
                    x.push(sub.ponded_qual[ui * np as usize + p]);
                }
                for lu in 0..nlu {
                    for p in 0..np {
                        let bu = self.surface_quality_.bu_idx(i, lu, p);
                        x.push(if bu < self.surface_quality_.buildup.len() { self.surface_quality_.buildup[bu] * crate::core::engine_initialize::buildup_norm(c, &self.landuse_solver_, ui, lu, p) } else { 0.0 });
                    }
                    let sw = ui * nlu as usize + lu as usize;
                    x.push(if sw < sub.sweep_last_date.len() {
                        sub.sweep_last_date[sw]
                    } else if sw < sub.sweep_last_swept.len() {
                        c.options.start_date - sub.sweep_last_swept[sw]
                    } else {
                        c.options.start_date
                    });
                }
            }
        }
        x
    }

    /// SWMMEngine::report().
    pub fn report(&mut self) -> Result<i32, Terminated> {
        if self.ctx_.state != EngineState::Ended {
            if self.ctx_.state == EngineState::ErrorState {
                return Ok(if self.ctx_.error_code != SWMM_OK { self.ctx_.error_code } else { SWMM_ERR_WRONG_STATE });
            }
            self.set_error(SWMM_ERR_WRONG_STATE, b"swmm_engine_report: must call end() first");
            return Ok(SWMM_ERR_WRONG_STATE);
        }

        // LID water balance results from the solver groups into lid_usage
        {
            let n_usage = self.ctx_.lid_usage.count();
            let n_groups = self.lid_.num_groups();
            for t in 0..n_groups {
                let g = self.lid_.group(t);
                for us in 0..g.count as usize {
                    let ctrl_idx = if us < g.control_idx.len() { g.control_idx[us] } else { -1 };
                    let sc_idx = if us < g.subcatch_idx.len() { g.subcatch_idx[us] } else { -1 };
                    if ctrl_idx < 0 || sc_idx < 0 {
                        continue;
                    }
                    let lu = &mut self.ctx_.lid_usage;
                    for uj in 0..n_usage as usize {
                        if lu.lid_index[uj] == ctrl_idx && lu.subcatch_index[uj] == sc_idx {
                            let gv = |v: &Vec<f64>| if us < v.len() { v[us] } else { 0.0 };
                            lu.wb_inflow[uj] = gv(&g.wb_inflow);
                            lu.wb_evap[uj] = gv(&g.wb_evap);
                            lu.wb_infil[uj] = gv(&g.wb_infil);
                            lu.wb_surf_flow[uj] = gv(&g.wb_surf_flow);
                            lu.wb_drain_flow[uj] = gv(&g.wb_drain_flow);
                            lu.wb_init_vol[uj] = gv(&g.wb_init_vol);
                            lu.wb_final_vol[uj] = gv(&g.wb_final_vol);
                            break;
                        }
                    }
                }
            }
        }

        // Gap #68: the inlet statistics for the report
        if self.ctx_.inlet_usages.count() > 0 {
            inlet_gather_stats(&self.inlet_, &mut self.ctx_);
        }

        self.ctx_.snowpacks.report_values.clear();
        if let Some(b) = ueb_ledger(&self.ctx_, self.snow_.state()) {
            let rows: [(&[u8], f64); 15] = [
                (b"Subcatchment UEB initial water (m3)", b.initial_m3),
                (b"Subcatchment UEB precipitation (m3)", b.precipitation_m3),
                (b"Subcatchment UEB liquid outflow (m3)", b.liquid_m3),
                (b"Subcatchment UEB throughfall (m3)", b.throughfall_m3),
                (b"Subcatchment UEB sublimation (m3)", b.sublimation_m3),
                (b"Subcatchment UEB removal (m3)", b.removed_m3),
                (b"Subcatchment UEB transfer in (m3)", b.transferred_m3),
                (b"Subcatchment UEB stored water (m3)", b.stored_m3),
                (b"Subcatchment UEB water residual (m3)", b.residual_m3),
                (b"Subcatchment UEB initial energy (kJ)", b.initial_energy_kj),
                (b"Subcatchment UEB net energy input (kJ)", b.net_energy_kj),
                (b"Subcatchment UEB depletion remapping (kJ)", b.geometry_energy_kj),
                (b"Subcatchment UEB plowed energy out (kJ)", b.plowed_energy_kj),
                (b"Subcatchment UEB stored energy (kJ)", b.stored_energy_kj),
                (b"Subcatchment UEB energy residual (kJ)", b.energy_residual_kj),
            ];
            self.ctx_.snowpacks.report_values = rows.iter().map(|(k, v)| (k.to_vec(), *v)).collect();
        }

        // OPENSWMM_HAS_2D: the mesh snow rows (phase 2b; refused)
        #[cfg(feature = "twod")]
        self.report_mesh_snow_rows()?;

        #[cfg(not(feature = "twod"))]
        let plugin_result = if self.plugins_.is_empty() { 0 } else { self.plugins_.write_summary_all(&self.ctx_) };
        // the report's 2D blocks read the router (upstream: ctx.twod_io.surface_router)
        #[cfg(feature = "twod")]
        let plugin_result = if self.plugins_.is_empty() { 0 } else { self.plugins_.write_summary_all(&self.ctx_, Some(&self.surface_router_)) };
        self.ctx_.state = EngineState::Reported;
        if plugin_result != 0 {
            let m = plugin_io_error(&self.plugins_, b"Cannot write report summary");
            self.set_error(SWMM_ERR_PLUGIN, &m);
            return Ok(SWMM_ERR_PLUGIN);
        }
        Ok(SWMM_OK)
    }

    /// SWMMEngine::close().
    pub fn close(&mut self) -> i32 {
        self.checkpoint_enabled_ = false;
        self.checkpoint_ready_ = false;
        self.checkpoint_output_ = false;
        self.binary_output_plugin_ = None;
        self.checkpoint_model_.clear();
        self.pending_hotstart_checkpoint_.clear();
        self.hotstart_checkpoint_error_.clear();
        #[cfg(feature = "twod")]
        self.close_2d();
        // perf:: dumps: timing only. io_thread_.stop(): nothing queued.

        if self.ctx_.state == EngineState::ErrorState {
            if !self.ctx_.error_message.is_empty() && !self.ctx_.errors.iter().any(|e| *e == self.ctx_.error_message) {
                let m = self.ctx_.error_message.clone();
                self.ctx_.errors.push(m);
            }
            if self.plugins_prepare_attempted_ {
                self.plugins_.finalize_all(&self.ctx_);
            } else {
                self.write_open_failure_report();
            }
        }
        self.plugins_prepare_attempted_ = false;

        self.iface_.close_files();
        self.rdii_iface_file_.close();
        self.close_runoff_iface();
        self.plugins_.unload_all();

        self.ctx_.state = EngineState::Closed;
        SWMM_OK
    }

    /// SWMMEngine::write_open_failure_report: the .rpt of a run that failed before start()
    /// prepared its plugins: a temporary report plugin's preamble (title, errors, warnings),
    /// its finalize, and its destructor's interrupted footer.
    pub(crate) fn write_open_failure_report(&mut self) {
        if self.rpt_path_.is_empty() || self.ctx_.options.rpt_disabled {
            return;
        }
        let mut rp = DefaultReportPlugin::new(self.rpt_path_.clone());
        rp.initialize(&[]);
        rp.prepare(&self.ctx_);
        rp.finalize(&self.ctx_);
        drop(rp);
    }
}

/// ~SWMMEngine: a running or ended engine is closed.
impl Drop for SwmmEngine {
    fn drop(&mut self) {
        if self.ctx_.state == EngineState::Running || self.ctx_.state == EngineState::Ended {
            self.close();
        }
    }
}

/// Gate only (cargo feature `gate`): what examples/outgate.rs needs beyond the loaders. Not
/// part of the translation.
#[cfg(feature = "gate")]
impl SwmmEngine {
    pub fn gate_old_gw_rate_mut(&mut self) -> &mut Vec<f64> {
        &mut self.old_gw_rate_
    }
    pub fn gate_plugins_mut(&mut self) -> &mut PluginFactory {
        &mut self.plugins_
    }
    /// fflush(nullptr) for the plugins' open files.
    pub fn gate_flush_plugin_files(&mut self) {
        for p in self.plugins_.output_plugins_mut().iter_mut() {
            if let Some(f) = p.out_file_.as_mut() {
                f.flush();
            }
        }
        for p in self.plugins_.report_plugins_mut().iter_mut() {
            if let Some(f) = p.file_.as_mut() {
                f.flush();
            }
        }
    }
    pub fn gate_post_output_snapshot(&mut self) -> Result<(), Terminated> {
        self.post_output_snapshot(0.0)
    }
}
