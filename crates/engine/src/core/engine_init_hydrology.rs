// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/SWMMEngine.cpp (initHydrology and
// setNextEvapDate) and SWMMEngine.hpp (the hydrology solver members)
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// SWMMEngine::initHydrology: the LID, runoff (with infiltration), gage, snow, climate,
// radiation and groundwater setup, and the project checks it makes (dates, climate
// file, [GWF] expressions, aquifer and LID parameters). The solver members it fills are
// fields of SwmmEngine (core/engine.rs).
//
// Differences of form:
//   * RunoffSolver::init can throw only for reviewed-area overrides, which initHydrology
//     never passes; its Result is therefore always Ok here.
//   * mathexpr::parse can throw only where gwf_validate's drift guard (the same parse)
//     already failed, and `||` stops before parse then. Were it ever to throw here,
//     upstream's noexcept would terminate; this translation treats it as a failed
//     compile instead (the [GWF] error), and the case is unreachable.
//   * The accessors at the end are for the gate dumper (examples/inithydro.rs); upstream
//     reads the private members directly (its harness compiles with private as public).

use crate::core::constants;
use crate::core::engine::{SwmmEngine, SWMM_ERR_PARSE};
use crate::core::error_codes::{format_error, format_error_detail, ERR_AQUIFER_PARAMS, ERR_LID_PARAMS, ERR_MATH_EXPR, ERR_NO_CLIMATE_FILE, ERR_REPORT_DATE, ERR_START_DATE};
use crate::core::options::Bytes;
use crate::core::unit_conversion as ucf;
use crate::data::name_index::ieq;
use crate::hydrology::climate_file::ClimateFileReader;
use crate::hydrology::climate_types::EvapMethod;
use crate::hydrology::snow::SnowSolver;
use crate::hydrology::gage::GageState;
use crate::hydrology::groundwater::{gwf_validate, GWSolver, GW_VAR_NAMES, GWV_MAX};
use crate::hydrology::lid::LIDSolver;
use crate::hydrology::rdii::RDIISolver;
use crate::hydrology::runoff::RunoffSolver;
use crate::math::math_expr::{self as mathexpr, Expression};

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut m = Vec::new();
    for p in parts {
        m.extend_from_slice(p);
    }
    m
}

impl SwmmEngine {
    /// Initialize hydrology solvers: runoff, snow, groundwater, LID. Populates the
    /// runoff store from subcatchment properties, sizes the gage states, and
    /// initializes the snow, groundwater and LID solvers.
    pub(crate) fn init_hydrology(&mut self) {
        // 1. LID solver — must run before the runoff solver: it fills
        //    subcatches.total_lid_area_ft2, which RunoffSolver::init subtracts (Gap #23,
        //    issue #131).
        self.lid_.init(&mut self.ctx_);
        // A4: size the per-layer species block against the units the manager just
        // built, and seed whatever they already hold.
        crate::transport::water_age_lid::init_lid_layer_age(&mut self.ctx_, &self.lid_);
        crate::transport::heat_lid::init_lid_layer_temperature(&mut self.ctx_, &self.lid_);

        // 2. Runoff solver: populate the runoff store from subcatchment properties
        let _ = self.runoff_.init(&mut self.ctx_, &[]);

        // 2b. Store base values for pattern-adjusted subcatchment parameters
        {
            let ns = self.ctx_.n_subcatches();
            let uns = ns as usize;
            let ctx = &mut self.ctx_;
            ctx.base_n_perv.resize(uns, 0.0);
            ctx.base_ds_perv.resize(uns, 0.0);
            for i in 0..ns {
                let ui = i as usize;
                ctx.base_n_perv[ui] = ctx.subcatches.n_perv[ui];
                ctx.base_ds_perv[ui] = ctx.subcatches.ds_perv[ui];
            }
            // Check if any pattern assignments exist
            ctx.has_subcatch_adj_patterns = false;
            for i in 0..ctx.subcatch_n_perv_pattern.len() {
                if ctx.subcatch_n_perv_pattern[i] >= 0 {
                    ctx.has_subcatch_adj_patterns = true;
                    break;
                }
            }
            if !ctx.has_subcatch_adj_patterns {
                for i in 0..ctx.subcatch_d_store_pattern.len() {
                    if ctx.subcatch_d_store_pattern[i] >= 0 {
                        ctx.has_subcatch_adj_patterns = true;
                        break;
                    }
                }
            }
            if !ctx.has_subcatch_adj_patterns {
                for i in 0..ctx.subcatch_infil_pattern.len() {
                    if ctx.subcatch_infil_pattern[i] >= 0 {
                        ctx.has_subcatch_adj_patterns = true;
                        break;
                    }
                }
            }
        }

        // 3. Gage states
        self.gage_states_.resize(self.ctx_.n_gages() as usize, GageState::default());

        // Both ordinary runoff and completed source owners use one snow setup.
        {
            let error = crate::hydrology::native_snow::initialize_native_snow(&mut self.ctx_, &mut self.snow_);
            if !error.is_empty() {
                self.ctx_.errors.push(error.clone());
                self.set_error(SWMM_ERR_PARSE, &error);
            }
        }

        // Date sanity (legacy project_validate, project.c:161-169), gated on an explicit
        // end date.
        if self.ctx_.options.end_date > 0.0 && self.ctx_.options.end_date <= self.ctx_.options.start_date {
            self.ctx_.errors.push(format_error(ERR_START_DATE, b""));
            let m = self.ctx_.errors.last().unwrap().clone();
            self.set_error(SWMM_ERR_PARSE, &m);
        } else if self.ctx_.options.report_start > 0.0 && self.ctx_.options.end_date <= self.ctx_.options.report_start {
            self.ctx_.errors.push(format_error(ERR_REPORT_DATE, b""));
            let m = self.ctx_.errors.last().unwrap().clone();
            self.set_error(SWMM_ERR_PARSE, &m);
        }

        // Climate-file requirement (legacy climate_validate, climate.c:850-858): FILE
        // wind, FILE evaporation, or TEMPERATURE (Hargreaves) evaporation need the
        // [TEMPERATURE] FILE.
        if (self.ctx_.options.wind_type == 1 || self.ctx_.options.evap_type == 4 || self.ctx_.options.evap_type == 3) && self.ctx_.options.temp_file.is_empty() {
            self.ctx_.errors.push(format_error(ERR_NO_CLIMATE_FILE, b""));
            let m = self.ctx_.errors.last().unwrap().clone();
            self.set_error(SWMM_ERR_PARSE, &m);
        }

        // 5. Climate: transfer evaporation data from options to climate state
        {
            let evap_type = self.ctx_.options.evap_type;
            {
                let ctx = &mut self.ctx_;
                let evaprate = ucf::ucf(ucf::EVAPRATE, &ctx.options);
                if evap_type == 0 {
                    // CONSTANT: EVAPRATE units follow the project's unit system.
                    ctx.climate_state.evap_method = EvapMethod::Constant;
                    ctx.climate_state.evaprate_ucf = evaprate;
                    // legacy keeps the constant in monthlyEvap[0]
                    ctx.climate_state.monthly_evap[0] = ctx.options.evap_values[0];
                    ctx.climate_state.evap_rate = ctx.options.evap_values[0] / evaprate;
                } else if evap_type == 1 {
                    // MONTHLY
                    ctx.climate_state.evap_method = EvapMethod::Monthly;
                    ctx.climate_state.evaprate_ucf = evaprate;
                    for i in 0..12 {
                        ctx.climate_state.monthly_evap[i] = ctx.options.evap_values[i];
                    }
                } else if evap_type == 2 {
                    ctx.climate_state.evap_method = EvapMethod::Timeseries;
                } else if evap_type == 3 {
                    ctx.climate_state.evap_method = EvapMethod::Temperature;
                } else if evap_type == 4 {
                    ctx.climate_state.evap_method = EvapMethod::Pan;
                }

                // Transfer latitude for Hargreaves ET calculation
                ctx.climate_state.latitude = ctx.options.snow_lat;

                // Site elevation (legacy Temp.elev = x[3] / UCF(LENGTH)) and the rain/snow
                // dividing temperature, degF internally.
                let si = ucf::get_unit_system(ctx.options.flow_units as i32) == 1;
                ctx.climate_state.elev = ctx.options.snow_elev / ucf::ucf(ucf::LENGTH, &ctx.options);
                ctx.climate_state.snow_divt = if si { 9.0 / 5.0 * ctx.options.snow_divt + 32.0 } else { ctx.options.snow_divt };

                // Longitude/solar-time correction: minutes -> hours (Temp.dtlong = x[5]/60).
                ctx.climate_state.dtlong = ctx.options.snow_dtlong / 60.0;

                // Transfer monthly adjustment arrays
                for i in 0..12 {
                    // legacy climate_validate: an SI deck's degC adjustment is a degF
                    // DIFFERENCE, x 9/5.
                    ctx.climate_state.adjust_temp[i] = ctx.adjust_temp[i] * (if si { 9.0 / 5.0 } else { 1.0 });
                    ctx.climate_state.adjust_evap[i] = ctx.adjust_evap[i] / ucf::ucf(ucf::EVAPRATE, &ctx.options);
                    ctx.climate_state.adjust_rain[i] = ctx.adjust_rain[i];
                    ctx.climate_state.adjust_hydcon[i] = ctx.adjust_hydcon[i];
                }

                // Resolve timeseries names to table indices
                if ctx.options.temp_source == 1 && !ctx.options.temp_ts_name.is_empty() {
                    ctx.climate_state.temp_ts_index = ctx.find_timeseries(&ctx.options.temp_ts_name);
                }

                // With NO temperature source, legacy's air temperature stays at its
                // project_setDefaults value, 70 degF (project.c:941).
                if ctx.options.temp_source == 0 && ctx.climate_state.temp_ts_index < 0 {
                    ctx.climate_state.temperature_src = 70.0;
                    ctx.climate_state.temperature = 70.0;
                }
                if evap_type == 2 && !ctx.options.evap_ts_name.is_empty() {
                    ctx.climate_state.evap_ts_index = ctx.find_timeseries(&ctx.options.evap_ts_name);
                }
            }
            // legacy climate_initState: NextEvapDate / NextEvapRate (step-wise series).
            {
                let start_date = self.ctx_.options.start_date;
                let cs = &mut self.ctx_.climate_state;
                cs.next_evap_date = start_date;
                cs.next_evap_rate = 0.0;
                cs.evap_ts_pos = -1;
                if cs.evap_method == EvapMethod::Timeseries && cs.evap_ts_index >= 0 {
                    let tbl = &self.ctx_.tables.tables[cs.evap_ts_index as usize];
                    if !tbl.x.is_empty() {
                        cs.evap_ts_pos = 0;
                        cs.next_evap_date = tbl.x[0];
                        cs.next_evap_rate = tbl.y[0];
                    }
                    if cs.next_evap_date < start_date {
                        self.set_next_evap_date(start_date);
                    }
                    let evaprate = ucf::ucf(ucf::EVAPRATE, &self.ctx_.options);
                    let cs = &mut self.ctx_.climate_state;
                    cs.evap_rate = cs.next_evap_rate / evaprate;
                    let d = cs.next_evap_date;
                    self.set_next_evap_date(d);
                }
            }
            let ctx = &mut self.ctx_;
            if ctx.options.humidity_type == 2 && !ctx.options.humidity_ts_name.is_empty() {
                ctx.climate_state.humidity_ts_index = ctx.find_timeseries(&ctx.options.humidity_ts_name);
            }

            // Resolve recovery pattern name to pattern index (case-insensitive)
            if !ctx.options.evap_recovery_pat.is_empty() {
                let np = ctx.patterns.count();
                for i in 0..np {
                    if ieq(&ctx.patterns.names[i as usize], &ctx.options.evap_recovery_pat) {
                        ctx.climate_state.recovery_pat_index = i;
                        break;
                    }
                }
            }

            // Open the climate file if temperature or evaporation uses FILE source, by
            // its RESOLVED path (falling back to the token when the resolver never ran).
            if (ctx.options.temp_source == 2 || evap_type == 4) && !ctx.options.temp_file.is_empty() {
                let us = ucf::get_unit_system(ctx.options.flow_units as i32);
                let climate_path: Bytes = if !ctx.options.temp_file.absolute.is_empty() { ctx.options.temp_file.absolute.clone() } else { ctx.options.temp_file.original.clone() };
                if !self.climate_file_.open(&climate_path, ctx.options.temp_file_start, us, ctx.options.temp_units) {
                    ctx.warnings.push(cat(&[
                        b"[TEMPERATURE] climate file '",
                        &climate_path,
                        "' could not be opened or its format was not recognised \u{2014} no climate data will be applied".as_bytes(),
                    ]));
                }
                // A FILE start date shifts which file day each simulation day reads.
                self.climate_file_.set_day_offset(if ctx.options.temp_file_start > 0.0 { ctx.options.temp_file_start.floor() - ctx.options.start_date.floor() } else { 0.0 });
                // legacy climate_openFile seeds both held temperatures with Temp.ta.
                ctx.climate_state.file_tmin = 70.0;
                ctx.climate_state.file_tmax = 70.0;
            }
        }

        {
            let error = crate::hydrology::ueb_forcing::initialize_radiation(&mut self.ctx_);
            if !error.is_empty() {
                self.ctx_.errors.push(error.clone());
                self.set_error(SWMM_ERR_PARSE, &error);
            }
        }

        // 6. Groundwater solver: init and populate from parsed aquifer/GW data
        self.groundwater_.init(self.ctx_.n_subcatches());
        {
            let ctx = &self.ctx_;
            let gw = self.groundwater_.state_mut();
            let unit_sys = ucf::get_unit_system(ctx.options.flow_units as i32) as usize;
            let u = |q: i32| ucf::UCF[q as usize][unit_sys];
            let sc = &ctx.subcatches;
            let aq = &ctx.aquifers;
            for i in 0..ctx.n_subcatches() {
                let ui = i as usize;
                let aq_idx = sc.gw_aquifer[ui];
                if aq_idx < 0 {
                    continue;
                }
                let uaq = aq_idx as usize;

                // Copy aquifer properties (legacy gwater.c line 170-182 conversions)
                gw.porosity[ui] = aq.porosity[uaq];
                gw.field_cap[ui] = aq.field_capacity[uaq];
                gw.wilt_point[ui] = aq.wilting_point[uaq];
                gw.k_sat[ui] = aq.conductivity[uaq] / u(ucf::RAINFALL);
                gw.k_slope[ui] = aq.conduct_slope[uaq];
                gw.tension_slope[ui] = aq.tension_slope[uaq] / u(ucf::LENGTH);
                gw.upper_evap_frac[ui] = aq.upper_evap[uaq];
                // Resolve upper evaporation pattern name to index
                gw.upper_evap_pat[ui] = -1;
                let pat_name = &aq.upper_evap_pat[uaq];
                if !pat_name.is_empty() {
                    gw.upper_evap_pat[ui] = ctx.patterns.find(pat_name);
                }
                gw.lower_evap_depth[ui] = aq.lower_evap[uaq] / u(ucf::LENGTH);
                gw.lower_loss_coeff[ui] = aq.lower_loss[uaq] / u(ucf::RAINFALL);
                // The [GROUNDWATER] row's optional Ebot / Wgw / Umc override the aquifer's
                // (legacy gwater_validate); Egwt overrides the node invert as the
                // lateral-flow threshold. Elevations convert to ft first.
                let ucf_len = u(ucf::LENGTH);
                let bot_elev = if sc.gw_bot_elev[ui] != constants::MISSING { sc.gw_bot_elev[ui] / ucf_len } else { aq.bottom_elev[uaq] / ucf_len };
                let wt_elev = if sc.gw_wt_elev[ui] != constants::MISSING { sc.gw_wt_elev[ui] / ucf_len } else { aq.water_table_elev[uaq] / ucf_len };
                let surf_elev = sc.gw_surf_elev[ui] / ucf_len;
                gw.total_depth[ui] = surf_elev - bot_elev;
                gw.bottom_elev[ui] = bot_elev;

                // Copy GW lateral flow coefficients
                gw.a1[ui] = sc.gw_a1[ui];
                gw.b1[ui] = sc.gw_b1[ui];
                gw.a2[ui] = sc.gw_a2[ui];
                gw.b2[ui] = sc.gw_b2[ui];
                gw.a3[ui] = sc.gw_a3[ui];
                {
                    let node_elev = if sc.gw_hstar[ui] != constants::MISSING {
                        sc.gw_hstar[ui] / ucf_len
                    } else {
                        let mut gwn = sc.gw_node[ui];
                        if gwn < 0 {
                            gwn = sc.outlet_node[ui];
                        }
                        if gwn >= 0 && gwn < ctx.n_nodes() {
                            ctx.nodes.invert_elev[gwn as usize]
                        } else {
                            0.0
                        }
                    };
                    gw.h_star[ui] = node_elev - bot_elev;
                }

                // Initial conditions (legacy gwater_initState)
                let mut theta0 = if sc.gw_upper_moist[ui] != constants::MISSING { sc.gw_upper_moist[ui] } else { aq.upper_moist[uaq] };
                if theta0 >= gw.porosity[ui] {
                    theta0 = gw.porosity[ui] - constants::XTOL;
                }
                gw.theta[ui] = theta0;
                let mut lower0 = wt_elev - bot_elev;
                if lower0 >= gw.total_depth[ui] {
                    lower0 = gw.total_depth[ui] - constants::XTOL;
                }
                gw.lower_depth[ui] = lower0;
            }
        }

        // Parse the [GWF] custom expressions, stored in ext_options during input parsing
        // as "GWF:<subcatch_name>:LATERAL" / "GWF:<subcatch_name>:DEEP".
        for i in 0..self.ctx_.n_subcatches() {
            let ui = i as usize;
            if self.ctx_.subcatches.gw_aquifer[ui] < 0 {
                continue;
            }
            let sc_name = self.ctx_.subcatch_names.names()[ui].clone();
            for (type_, lateral) in [(b"LATERAL" as &[u8], true), (b"DEEP" as &[u8], false)] {
                // An expression that fails to compile is a model error (legacy
                // ERR_MATH_EXPR); gwf_validate supplies the diagnostic.
                let key = cat(&[b"GWF:", &sc_name, b":", type_]);
                let text = match self.ctx_.options.ext_options.get(&key) {
                    None => continue,
                    Some(t) if t.is_empty() => continue,
                    Some(t) => t.clone(),
                };
                let mut msg: Bytes = Bytes::new();
                let mut col: i32 = -1;
                let mut expr = Expression::default();
                let failed = gwf_validate(&text, &mut msg, &mut col) != 0 || mathexpr::parse(&text, &mut expr).map(|rc| rc != 0).unwrap_or(true) || !expr.valid;
                if failed {
                    let detail = cat(&[b"in [GWF] ", type_, b" for Subcatchment ", &sc_name, b": ", &msg]);
                    self.ctx_.errors.push(format_error_detail(ERR_MATH_EXPR, b"", &detail));
                    let m = self.ctx_.errors.last().unwrap().clone();
                    self.set_error(SWMM_ERR_PARSE, &m);
                    continue;
                }
                mathexpr::bind_variables(&mut expr, &GW_VAR_NAMES, GWV_MAX);
                let gw = self.groundwater_.state_mut();
                if lateral {
                    gw.lateral_expr[ui] = expr;
                } else {
                    gw.deep_expr[ui] = expr;
                }
            }
        }

        // Gap #81: Aquifer parameter validation
        {
            let naq = self.ctx_.aquifers.count();
            for a in 0..naq {
                let ua = a as usize;
                let aq = &self.ctx_.aquifers;
                let por = aq.porosity[ua];
                let fc = aq.field_capacity[ua];
                let wp = aq.wilting_point[ua];
                let ksat = aq.conductivity[ua];
                let uevap = aq.upper_evap[ua];
                let bad = por <= 0.0 || por > 1.0 || fc < 0.0 || fc >= por || wp < 0.0 || wp >= fc || ksat < 0.0 || uevap < 0.0 || uevap > 1.0;
                if bad {
                    let e = format_error(ERR_AQUIFER_PARAMS, &aq.names[ua]);
                    self.ctx_.errors.push(e);
                    let m = self.ctx_.errors.last().unwrap().clone();
                    self.set_error(SWMM_ERR_PARSE, &m);
                }
            }
        }

        // 7. LID solver: initialized first (step 1).

        // Gap #82: LID parameter validation
        {
            let nlid = self.ctx_.lid_controls.count();
            for j in 0..nlid {
                let uj = j as usize;
                let lc = &self.ctx_.lid_controls;
                let lname = lc.names[uj].clone();
                let mut bad: Vec<bool> = Vec::new();
                // Soil layer: [0]=thickness, [1]=porosity, [2]=fieldCap, [3]=wiltPt
                if uj < lc.soil.len() {
                    let s = &lc.soil[uj];
                    if s[0] > 0.0 && (s[1] <= 0.0 || s[1] > 1.0 || s[2] < 0.0 || s[2] >= s[1] || s[3] < 0.0 || s[3] >= s[2]) {
                        bad.push(true);
                    }
                }
                // Storage layer: [0]=thickness, [1]=void RATIO (issue #102)
                if uj < lc.storage.len() {
                    let s = &lc.storage[uj];
                    if s[0] > 0.0 && s[1] <= 0.0 {
                        bad.push(true);
                    }
                }
                // Pavement layer: [0]=thickness, [1]=void RATIO — same convention.
                if uj < lc.pavement.len() {
                    let p = &lc.pavement[uj];
                    if p[0] > 0.0 && p[1] <= 0.0 {
                        bad.push(true);
                    }
                }
                // Drain layer: [4]=hOpen, [5]=hClose — hOpen must be >= hClose
                if uj < lc.drain.len() {
                    let d = &lc.drain[uj];
                    if d[0] > 0.0 && d[4] < d[5] {
                        bad.push(true);
                    }
                }
                // Each failed layer check pushes its own ERR_LID_PARAMS, in layer order.
                for _ in bad {
                    self.ctx_.errors.push(format_error(ERR_LID_PARAMS, &lname));
                    let m = self.ctx_.errors.last().unwrap().clone();
                    self.set_error(SWMM_ERR_PARSE, &m);
                }
            }
        }
    }

    /// legacy climate.c setNextEvapDate: the date the evaporation source next changes.
    pub(crate) fn set_next_evap_date(&mut self, the_date: f64) {
        let cs = &mut self.ctx_.climate_state;
        if cs.next_evap_date > the_date {
            return;
        }
        match cs.evap_method {
            EvapMethod::Constant => {
                cs.next_evap_date = the_date + 365.0;
            }
            EvapMethod::Monthly => {
                let (mut yr, mut mon, _day) = crate::core::datetime::decode_date(the_date);
                if mon == 12 {
                    mon = 1;
                    yr += 1;
                } else {
                    mon += 1;
                }
                cs.next_evap_date = crate::core::datetime::encode_date(yr, mon, 1);
            }
            EvapMethod::Timeseries => {
                cs.next_evap_date = the_date + 365.0;
                if cs.evap_ts_index >= 0 {
                    let tbl = &self.ctx_.tables.tables[cs.evap_ts_index as usize];
                    let n = tbl.x.len() as i32;
                    let end_date = crate::core::datetime::add_seconds(self.ctx_.options.start_date, self.ctx_.options.total_duration_ms() / 1000.0);
                    let mut pos = cs.evap_ts_pos;
                    while pos + 1 < n {
                        pos += 1;
                        let d = tbl.x[pos as usize];
                        let e = tbl.y[pos as usize];
                        if d > end_date {
                            break;
                        }
                        if d >= the_date {
                            cs.next_evap_date = d;
                            cs.next_evap_rate = e;
                            break;
                        }
                    }
                    cs.evap_ts_pos = pos;
                }
            }
            EvapMethod::Pan => {
                cs.next_evap_date = the_date.floor() + 1.0;
            }
            // TEMPERATURE takes upstream's default branch
            _ => {
                cs.next_evap_date = the_date + 365.0;
            }
        }
    }

    // ---- accessors for the gate dumper ------------------------------------

    pub fn runoff_solver(&self) -> &RunoffSolver {
        &self.runoff_
    }
    pub fn snow_solver(&self) -> &SnowSolver {
        &self.snow_
    }
    pub fn groundwater_solver(&self) -> &GWSolver {
        &self.groundwater_
    }
    pub fn lid_solver(&self) -> &LIDSolver {
        &self.lid_
    }
    pub fn rdii_solver(&self) -> &RDIISolver {
        &self.rdii_
    }
    pub fn climate_file_reader(&mut self) -> &mut ClimateFileReader {
        &mut self.climate_file_
    }
    pub fn gage_states(&self) -> &Vec<GageState> {
        &self.gage_states_
    }
    pub fn last_melt_doy(&self) -> i32 {
        self.last_melt_doy_
    }
    /// gw_frac_perv_, gw_perv_evap_, gw_lid_perv_area_, gw_lid_infil_vol_,
    /// gw_lid_perv_evap_vol_, gw_2d_recharge_vol_, gw_infil_with_2d_, gw_infil_rate_,
    /// snow_rain_, snow_snow_, gw_q_interp_, old_gw_rate_, in that order.
    pub fn hydrology_vectors(&self) -> [&Vec<f64>; 12] {
        [
            &self.gw_frac_perv_,
            &self.gw_perv_evap_,
            &self.gw_lid_perv_area_,
            &self.gw_lid_infil_vol_,
            &self.gw_lid_perv_evap_vol_,
            &self.gw_2d_recharge_vol_,
            &self.gw_infil_with_2d_,
            &self.gw_infil_rate_,
            &self.snow_rain_,
            &self.snow_snow_,
            &self.gw_q_interp_,
            &self.old_gw_rate_,
        ]
    }
    pub fn gw_q_node(&self) -> &Vec<i32> {
        &self.gw_q_node_
    }
    /// SWMMEngine::initHydrology, for harnesses (upstream's is private).
    pub fn init_hydrology_for_gate(&mut self) {
        self.init_hydrology();
    }
    /// `rdii_.init(ctx_)`, the call upstream's initHydraulics makes (step 10a), for
    /// harnesses.
    pub fn init_rdii_for_gate(&mut self) {
        self.rdii_.init(&mut self.ctx_);
    }
}
