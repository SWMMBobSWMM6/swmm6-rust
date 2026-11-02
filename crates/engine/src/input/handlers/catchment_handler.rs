// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/handlers/CatchmentHandler.cpp and
// CatchmentHandler.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Section handlers for [SUBCATCHMENTS], [SUBAREAS], [INFILTRATION] and [RAINGAGES].
//
// Differences in form only: `ensure_gage_capacity`'s generic `grow` lambda is written out
// per column. Upstream's `Table{}` for a new rain_series row is Table::new(b"", TIMESERIES),
// the same value GageData's own sizing uses. Assigning a std::string to a FilePathPair is
// FilePathPair::set (original = token, absolute cleared).
//
// `split_series_file_token` is MultiColumnSeriesFile.cpp's, from
// crate::input::multi_column_series_file.

use crate::core::cmath::{cmin, cvt_i32};
use crate::core::context::SimulationContext;
use crate::core::error_codes::{format_error, ERR_DATETIME, ERR_INFILTRATION};
use crate::core::options::{Bytes, InfiltrationModel};
use crate::data::gage::{GageData, RainFileFormat, RainSource};
use crate::data::table::{Table, TableType};
use crate::input::parse_utils::{add_unique, parse_date, parse_time_seconds, to_double};
use crate::input::section_parser::parse_section;
use crate::input::tokenizer;

use crate::input::multi_column_series_file::split_series_file_token;

/// Extend base subcatchment columns through a registered row index.
fn ensure_subcatch_capacity(ctx: &mut SimulationContext, idx: i32) {
    ctx.subcatches.grow_to(idx + 1);
}

/// Grow rain-gage authoring/runtime columns through a registered index, keeping existing
/// rows.
fn ensure_gage_capacity(ctx: &mut SimulationContext, idx: i32) {
    let n = (idx + 1) as usize;
    let g = &mut ctx.gages;
    if g.source.len() < n {
        g.source.resize(n, RainSource::Timeseries);
    }
    if g.ts_index.len() < n {
        g.ts_index.resize(n, -1);
    }
    if g.ts_name.len() < n {
        g.ts_name.resize(n, Bytes::new());
    }
    if g.file_path.len() < n {
        g.file_path.resize(n, Default::default());
    }
    if g.col_name.len() < n {
        g.col_name.resize(n, Bytes::new());
    }
    if g.station_id.len() < n {
        g.station_id.resize(n, Bytes::new());
    }
    if g.rain_units.len() < n {
        g.rain_units.resize(n, 0);
    }
    if g.file_first_date.len() < n {
        g.file_first_date.resize(n, 0.0);
    }
    if g.file_last_date.len() < n {
        g.file_last_date.resize(n, 0.0);
    }
    if g.file_periods_precip.len() < n {
        g.file_periods_precip.resize(n, 0);
    }
    if g.file_start_date.len() < n {
        g.file_start_date.resize(n, 0.0);
    }
    if g.rain_series.len() < n {
        g.rain_series.resize_with(n, || Table::new(b"", TableType::Timeseries));
    }
    if g.file_format.len() < n {
        g.file_format.resize(n, RainFileFormat::Unknown);
    }
    if g.interval_sec.len() < n {
        g.interval_sec.resize(n, 3600);
    }
    if g.snow_factor.len() < n {
        g.snow_factor.resize(n, 1.0);
    }
    if g.scale_factor.len() < n {
        g.scale_factor.resize(n, 1.0);
    }
    if g.rain_type.len() < n {
        g.rain_type.resize(n, 0);
    }
    if g.rainfall.len() < n {
        g.rainfall.resize(n, 0.0);
    }
    if g.next_rainfall.len() < n {
        g.next_rainfall.resize(n, 0.0);
    }
    if g.api_rainfall.len() < n {
        g.api_rainfall.resize(n, -1.0); // -1 = no API override
    }
    if g.next_rain_date.len() < n {
        g.next_rain_date.resize(n, 0.0);
    }
    if g.is_raining.len() < n {
        g.is_raining.resize(n, false);
    }
    if g.st_cur.len() < n {
        g.st_cur.resize(n, -2);
    }
    if g.st_next.len() < n {
        g.st_next.resize(n, -1);
    }
    if g.st_rain.len() < n {
        g.st_rain.resize(n, 0.0);
    }
    if g.st_init.len() < n {
        g.st_init.resize(n, 0u8);
    }
    if g.st_used.len() < n {
        g.st_used.resize(n, 0u8);
    }

    // Past-rain tracking (used by updateAllGages for controls past-rain)
    if g.past_rain_accum.len() < n {
        g.past_rain_accum.resize(n, 0.0);
    }
    if g.past_rain_interval.len() < n {
        g.past_rain_interval.resize(n, 0.0);
    }
    if g.cumul_rain_accum.len() < n {
        g.cumul_rain_accum.resize(n, 0.0);
    }
    if g.co_gage_index.len() < n {
        g.co_gage_index.resize(n, -1);
    }
    if g.comments.len() < n {
        g.comments.resize(n, Bytes::new());
    }
    // past_rain is a flat 2-D array: [gage * MAXPASTRAIN + hour]
    {
        let nr = n * GageData::MAXPASTRAIN as usize;
        if g.past_rain.len() < nr {
            g.past_rain.resize(nr, 0.0);
        }
    }
}

// ============================================================================
// handle_subcatchments()
// ============================================================================

/// [SUBCATCHMENTS]: Name Gage Outlet Area %Imperv Width %Slope [CurbLen] [SnowPack]
/// [RainScale] [SnowScale].
pub fn handle_subcatchments(ctx: &mut SimulationContext, lines: &[Bytes]) {
    // Pre-reserve from the section's row count (an upper bound). Capacity only.
    let cnt = ctx.subcatches.count();
    ctx.subcatches.reserve_to(cnt + lines.len() as i32);
    let ns = ctx.subcatch_names.size() as usize;
    ctx.subcatch_names.reserve(ns + lines.len());
    for pl in parse_section(lines) {
        let tok = tokenizer::tokenize(&pl.data);
        if tok.len() < 7 {
            continue;
        }
        // Snowpack (tok 8) accepts "*" as a placeholder meaning "no snow pack".
        // RainScale (tok 9) and SnowScale (tok 10) default to 1.0.

        let name = &tok[0];
        let idx = add_unique(&mut ctx.subcatch_names, name, &mut ctx.errors);
        if idx < 0 {
            continue; // duplicate ID (ERR 207, legacy input.c parity)
        }
        let u = idx as usize;

        ensure_subcatch_capacity(ctx, idx);

        // Gage: may be -1 if not parsed yet; the name lets PostParseResolver re-resolve.
        ctx.subcatches.gage_name[u] = tok[1].clone();
        ctx.subcatches.gage[u] = ctx.gage_names.find(&tok[1]);

        // Outlet: a node or another subcatchment; the name is kept for deferred resolution.
        ctx.subcatches.outlet_name[u] = tok[2].clone();
        let node_idx = ctx.node_names.find(&tok[2]);
        if node_idx >= 0 {
            ctx.subcatches.outlet_node[u] = node_idx;
        } else {
            let sub_idx = ctx.subcatch_names.find(&tok[2]);
            if sub_idx >= 0 {
                ctx.subcatches.outlet_subcatch[u] = sub_idx;
            }
        }

        ctx.subcatches.area[u] = to_double(&tok[3], 0.0);
        // Legacy subcatch_readParams: fracImperv = MIN(x[4], 100) / 100.
        ctx.subcatches.frac_imperv[u] = cmin(to_double(&tok[4], 0.0), 100.0) / 100.0; // % -> fraction
        ctx.subcatches.width[u] = to_double(&tok[5], 0.0);
        ctx.subcatches.slope[u] = to_double(&tok[6], 0.0) / 100.0; // % -> fraction

        // Optional CurbLen (column 7)
        if tok.len() > 7 {
            ctx.subcatches.curb_length[u] = to_double(&tok[7], 0.0);
        }

        // Optional SnowPack (column 8); "*" is a positional placeholder.
        if tok.len() > 8 && !tok[8].is_empty() && tok[8] != b"*" {
            ctx.subcatches.snowpack_name[u] = tok[8].clone();
            ctx.subcatches.snowpack[u] = ctx.snowpack_names.find(&tok[8]);
        }

        // Optional rainfall / snowfall scale factors (columns 9, 10); non-positive ignored.
        if tok.len() > 9 {
            let rsf = to_double(&tok[9], 0.0);
            if rsf > 0.0 {
                ctx.subcatches.rain_scale_factor[u] = rsf;
            }
        }
        if tok.len() > 10 {
            let ssf = to_double(&tok[10], 0.0);
            if ssf > 0.0 {
                ctx.subcatches.snow_scale_factor[u] = ssf;
            }
        }

        if !pl.comment.is_empty() {
            ctx.subcatches.comments[u] = pl.comment.clone();
        }
    }
}

// ============================================================================
// handle_subareas()
// ============================================================================

/// [SUBAREAS]: Subcatch N-Imperv N-Perv S-Imperv S-Perv PctZero [RouteTo [PctRouted]].
pub fn handle_subareas(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 6 {
            continue;
        }

        let idx = ctx.subcatch_names.find(&tok[0]);
        if idx < 0 {
            continue;
        }
        let u = idx as usize;

        ctx.subcatches.n_imperv[u] = to_double(&tok[1], 0.0);
        ctx.subcatches.n_perv[u] = to_double(&tok[2], 0.0);
        ctx.subcatches.ds_imperv[u] = to_double(&tok[3], 0.0);
        ctx.subcatches.ds_perv[u] = to_double(&tok[4], 0.0);

        // PctZero: the authored percent and its 0-1 fraction are both kept.
        ctx.subcatches.pct_zero[u] = to_double(&tok[5], 0.0);
        ctx.subcatches.frac_imperv_no_store[u] = to_double(&tok[5], 0.0) / 100.0;

        // RouteTo: OUTLET (0), IMPERV (1), or PERV (2), keyword as a case-insensitive
        // PREFIX of the token (legacy findmatch).
        if tok.len() > 6 {
            let rt = tokenizer::to_upper(&tok[6]);
            let mut mode = 0; // OUTLET
            if rt.starts_with(b"OUTLET") {
                mode = 0;
            } else if rt.starts_with(b"IMPERV") {
                mode = 1;
            } else if rt.starts_with(b"PERV") {
                mode = 2;
            }
            ctx.subcatches.subarea_routing[u] = mode;

            // PctRouted defaults to 100%, then tok[7] overrides (every RouteTo).
            if tok.len() > 7 {
                ctx.subcatches.pct_routed[u] = to_double(&tok[7], 0.0) / 100.0;
            } else {
                ctx.subcatches.pct_routed[u] = 1.0;
            }
        }
    }
}

// ============================================================================
// handle_infiltration()
// ============================================================================

/// [INFILTRATION]: Subcatch Param1 Param2 Param3 [Param4] [Param5] [Method].
pub fn handle_infiltration(ctx: &mut SimulationContext, lines: &[Bytes]) {
    // Determine infiltration model from OPTIONS enum
    let model = ctx.options.infiltration as i32;

    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 4 {
            continue;
        }

        let idx = ctx.subcatch_names.find(&tok[0]);
        if idx < 0 {
            continue;
        }
        let u = idx as usize;

        // A per-subcatchment method override as the LAST token (legacy infil.c:150-155).
        let mut ntoks = tok.len();
        let mut row_model = model;
        {
            let last = tokenizer::to_upper(&tok[ntoks - 1]);
            if last == b"HORTON" {
                row_model = 0;
                ntoks -= 1;
            } else if last == b"MOD_HORTON" || last == b"MODIFIED_HORTON" {
                row_model = 1;
                ntoks -= 1;
            } else if last == b"GREEN_AMPT" {
                row_model = 2;
                ntoks -= 1;
            } else if last == b"MOD_GREEN_AMPT" || last == b"MODIFIED_GREEN_AMPT" {
                row_model = 3;
                ntoks -= 1;
            } else if last == b"CURVE_NUMBER" {
                row_model = 4;
                ntoks -= 1;
            }
        }

        ctx.subcatches.infil_model[u] = row_model;
        if ntoks > 1 {
            ctx.subcatches.infil_p1[u] = to_double(&tok[1], 0.0);
        }
        if ntoks > 2 {
            ctx.subcatches.infil_p2[u] = to_double(&tok[2], 0.0);
        }
        if ntoks > 3 {
            ctx.subcatches.infil_p3[u] = to_double(&tok[3], 0.0);
        }
        if ntoks > 4 {
            ctx.subcatches.infil_p4[u] = to_double(&tok[4], 0.0);
        }
        if ntoks > 5 {
            ctx.subcatches.infil_p5[u] = to_double(&tok[5], 0.0);
        }

        // Parameter validation on the RAW display-unit values (legacy ERROR 235).
        let p1 = ctx.subcatches.infil_p1[u];
        let p2 = ctx.subcatches.infil_p2[u];
        let p3 = ctx.subcatches.infil_p3[u];
        let p4 = ctx.subcatches.infil_p4[u];
        let p5 = ctx.subcatches.infil_p5[u];
        let mut ok = true;
        // switch (static_cast<InfiltrationModel>(row_model))
        if row_model == InfiltrationModel::Horton as i32 || row_model == InfiltrationModel::ModHorton as i32 {
            // horton_setParams (infil.c:328): no negatives, f0 >= fmin
            ok = p1 >= 0.0 && p2 >= 0.0 && p3 >= 0.0 && p4 >= 0.0 && p5 >= 0.0 && p1 >= p2;
        } else if row_model == InfiltrationModel::GreenAmpt as i32 || row_model == InfiltrationModel::ModGreenAmpt as i32 {
            // grnampt_setParams (infil.c:605): suction >= 0, Ksat > 0, IMD in [0, 1]
            ok = p1 >= 0.0 && p2 > 0.0 && p3 >= 0.0 && p3 <= 1.0;
        } else if row_model == InfiltrationModel::CurveNumber as i32 {
            // curvenum_setParams (infil.c): drying time (3rd param) > 0
            ok = p3 > 0.0;
        }
        if !ok {
            ctx.errors.push(format_error(ERR_INFILTRATION, b""));
        }
    }
}

// ============================================================================
// handle_raingages()
// ============================================================================

/// [RAINGAGES]: Name Format Interval SCF TIMESERIES name [Scale] | FILE path [Station
/// Units [StartDate] [Scale]].
pub fn handle_raingages(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for pl in parse_section(lines) {
        let tok = tokenizer::tokenize(&pl.data);
        if tok.len() < 5 {
            continue;
        }

        let name = &tok[0];
        let idx = add_unique(&mut ctx.gage_names, name, &mut ctx.errors);
        if idx < 0 {
            continue; // duplicate ID (ERR 207, legacy input.c parity)
        }
        let u = idx as usize;

        ensure_gage_capacity(ctx, idx);

        // Format: VOLUME or INTENSITY -> rain_type
        let fmt = tokenizer::to_upper(&tok[1]);
        if fmt == b"VOLUME" {
            ctx.gages.rain_type[u] = 1;
        } else if fmt == b"CUMULATIVE" {
            ctx.gages.rain_type[u] = 2;
        } else {
            ctx.gages.rain_type[u] = 0; // INTENSITY
        }

        // Recording interval: DECIMAL HOURS or "H:MM[:SS]" (legacy datetime_strToTime).
        let itok = &tok[2];
        let interval_secs = if itok.contains(&b':') {
            parse_time_seconds(itok) // H:MM[:SS] clock duration
        } else {
            to_double(itok, 0.0) * 3600.0 // decimal hours
        };
        ctx.gages.interval_sec[u] = cvt_i32(interval_secs + 0.5);

        // Snow correction factor
        ctx.gages.snow_factor[u] = to_double(&tok[3], 1.0);

        // Source: TIMESERIES <name>  or  FILE "<path>"  or  FILE "<path:col>"
        let src = tokenizer::to_upper(&tok[4]);

        if src == b"TIMESERIES" && tok.len() > 5 {
            ctx.gages.source[u] = RainSource::Timeseries;
            ctx.gages.ts_name[u] = tok[5].clone(); // name kept for deferred resolution
            ctx.gages.ts_index[u] = ctx.find_timeseries(&tok[5]);

            // Optional trailing rainfall scaling factor (legacy gage.c tok[6])
            if tok.len() > 6 {
                let sf = to_double(&tok[6], 1.0);
                if sf > 0.0 {
                    ctx.gages.scale_factor[u] = sf;
                }
            }
        } else if src == b"FILE" && tok.len() > 5 {
            ctx.gages.source[u] = RainSource::FileRain;

            // "path:COLUMN" by the shared split rule (MultiColumnSeriesFile).
            let mut file_only = Bytes::new();
            let mut file_col = Bytes::new();
            if split_series_file_token(&tok[5], &mut file_only, &mut file_col) {
                ctx.gages.file_path[u].set(&file_only);
                ctx.gages.col_name[u] = file_col;
                ctx.gages.file_format[u] = RainFileFormat::UserCsv;
            } else {
                ctx.gages.file_path[u].set(&file_only);
                ctx.gages.file_format[u] = RainFileFormat::StanPrcp;
            }

            // Compact FILE "path" [ScaleFactor] is a CSV/TSF gage on its first data column.
            if ctx.gages.file_format[u] == RainFileFormat::StanPrcp && tok.len() == 7 && to_double(&tok[6], -1.0) > 0.0 {
                ctx.gages.file_format[u] = RainFileFormat::UserCsv;
            }

            // Standard SWMM FILE grammar (legacy gage.c gage_readParams):
            //   Name Format Interval SCF FILE Fname Station Units [StartDate] [SCF]
            if ctx.gages.file_format[u] == RainFileFormat::StanPrcp {
                // '*' is the writer's placeholder for "no station".
                if tok.len() > 6 && tok[6] != b"*" {
                    ctx.gages.station_id[u] = tok[6].clone();
                }
                if tok.len() > 7 {
                    let un = tokenizer::to_upper(&tok[7]);
                    ctx.gages.rain_units[u] = if un == b"MM" { 1 } else { 0 }; // default IN
                }
                // Optional start date (legacy gage.c:207-212): a leading `*` means none;
                // anything else must be a date or it is ERROR 213.
                ctx.gages.file_start_date[u] = 0.0;
                if tok.len() > 8 && !tok[8].is_empty() && tok[8][0] != b'*' {
                    let start = parse_date(&tok[8]);
                    if start > 0.0 {
                        ctx.gages.file_start_date[u] = start;
                    } else {
                        ctx.errors.push(format_error(ERR_DATETIME, &tok[8]));
                    }
                }
                // The trailing scale factor, only as a positive number.
                if tok.len() > 9 {
                    let sf = to_double(&tok[9], 1.0);
                    if sf > 0.0 {
                        ctx.gages.scale_factor[u] = sf;
                    }
                }
            } else if tok.len() > 6 {
                let sf = to_double(&tok[6], 1.0);
                if sf > 0.0 {
                    ctx.gages.scale_factor[u] = sf;
                }
            }
        }
        if !pl.comment.is_empty() {
            ctx.gages.comments[u] = pl.comment.clone();
        }
    }
}
