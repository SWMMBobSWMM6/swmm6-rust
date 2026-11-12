// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/SWMMEngine.cpp (HydroCouple/openswmm.engine,
// branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the
// Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The runoff (hydrology) step: SWMMEngine::stepRunoff and what it calls on the engine:
// computeRunoffTimestep (setNextEvapDate is in engine_init_hydrology.rs), assembleRunon,
// assembleGWCoupling, stepGroundwater, accumulateRunoffMassBalance, stepSurfaceQuality,
// bookWashoffLoads, saveRunoffIfaceStep and subcatchSnowDepth, and the LID blocks A6a/A6b
// and the volume-domain runoff fold. The solvers it drives: hydrology/runoff.rs (the
// nonlinear reservoir), infiltration.rs and infil_bank.rs, gage.rs, climate.rs, snow.rs
// (degree-day snow), groundwater.rs, rdii.rs, lid_kernel.rs (the LID process kernel),
// math/ode_solver.rs, quality/surface_species_kernel.rs (buildup and washoff) and
// transport/water_age_watershed.rs.
//
// The step's signatures are the agreed ones (PLAN.md, "The step's signatures"):
//   SwmmEngine::step_runoff(&mut self, dt_routing: f64) -> Result<(), Terminated>
//   SwmmEngine::step_groundwater(&mut self, dt_runoff: f64) -> Result<(), Terminated>
//   SwmmEngine::step_surface_quality(&mut self, dt_runoff: f64) -> Result<(), Terminated>
//
// Differences of form, never in result:
//   * stepRunoff throws std::runtime_error in three places (an hourly meteorology average
//     error, a UEB advance error, the end of a USE RUNOFF file); each is an Err(Terminated)
//     whose exception starts with engine_step::THROWN_PREFIX, which step() turns into
//     upstream's catch.
//   * The 1D build's spatial-source map is absent (ordinarySourceOwners() is null), so
//     ownsSpatialSource is false everywhere; the calls are kept where upstream makes them.
//   * The local `total_sec_clamp` at the top of stepRunoff is computed and never used
//     upstream; it is left out.
//   * Code a deck reaches whose solver is not translated yet returns Err(stub(...)) at the
//     point upstream would first act (see `stub_if`). NOT YET TRANSLATED: the UEB
//     EQUILIBRIUM surface-temperature kernel (refused before ueb::advanceSubcatchments)
//     and the runoff interface file (USE RUNOFF readResults, SAVE RUNOFF saveResults).

use crate::core::cmath::{cclamp, cmax, cmin};
use crate::core::context::FileMode;
use crate::core::datetime;
use crate::core::engine::SwmmEngine;
use crate::core::engine_init_hydraulics::Terminated;
use crate::core::engine_step::{stub, THROWN_PREFIX};
use crate::core::unit_conversion as ucf;
use crate::data::forcing::ForcingMode;
use crate::data::lid_layer_species::LidSpecies;
use crate::data::pollutant::MassUnits;
use crate::data::table::table_lookup_cursor;
use crate::data::water_age::WaterAgeSource;
use crate::hydrology::climate;
use crate::hydrology::climate_file::DailyClimateRecord;
use crate::hydrology::climate_types::EvapMethod;
use crate::hydrology::gage;
use crate::hydrology::snow;
use crate::quality::surface_species_kernel as sskernel;

fn stub_if(cond: bool, what: &'static str) -> Result<(), Terminated> {
    if cond {
        Err(stub(what))
    } else {
        Ok(())
    }
}

/// A std::runtime_error thrown toward step()'s catch.
fn thrown(what: &[u8]) -> Terminated {
    let w = String::from_utf8_lossy(what).into_owned();
    let _ = THROWN_PREFIX;
    Terminated { exception: "throw std::runtime_error", what: Box::leak(w.into_boxed_str()) }
}

/// transport::heat::saturationVapourPressure (HeatFluxModules/SurfaceExchange.cpp), kPa.
fn saturation_vapour_pressure(t_c: f64) -> f64 {
    0.61275 * (17.27 * t_c / (237.3 + t_c)).exp() // LIBM
}

impl SwmmEngine {
    /// SWMMEngine::stepRunoff: runoff sub-stepping on the runoff clock until it passes the
    /// end of this routing step, then the interpolation to the node inflow arrays.
    pub(crate) fn step_runoff(&mut self, dt_routing: f64) -> Result<(), Terminated> {
        // The trigger is legacy's, on the ms clocks (swmm5.c execRouting).
        let next_routing_ms = cmin(self.ctx_.elapsed_ms + 1000.0 * dt_routing, self.ctx_.options.total_duration_ms());
        while self.new_runoff_ms_ < next_routing_ms {
            self.runoff_substep()?;
        }
        self.interpolate_runoff_inflows();
        Ok(())
    }

    /// The body of stepRunoff's `while (new_runoff_ms_ < next_routing_ms)` loop.
    fn runoff_substep(&mut self) -> Result<(), Terminated> {
        // Save old runoff/runon/conc + GW/snow/LID-drain state for interpolation.
        self.ctx_.subcatches.save_state();
        if self.ctx_.reactions.surface.active() {
            self.ctx_.reactions.surface.washoff_conc_old = self.ctx_.reactions.surface.washoff_conc.clone();
        }
        for i in 0..self.ctx_.n_subcatches() {
            if self.owns_spatial_source(i) {
                continue;
            }
            let ui = i as usize;
            self.ctx_.subcatches.old_gw_flow[ui] = self.ctx_.subcatches.gw_flow[ui];
            if ui < self.old_gw_rate_.len() {
                self.old_gw_rate_[ui] = if ui < self.groundwater_.state().gw_flow.len() { self.groundwater_.state().gw_flow[ui] } else { 0.0 };
            }
            self.ctx_.subcatches.old_snow_depth[ui] = self.ctx_.subcatches.snow_depth[ui];
            self.ctx_.subcatches.old_lid_drain_flow[ui] = self.ctx_.subcatches.lid_drain_flow[ui];
            self.ctx_.subcatches.lid_drain_flow[ui] = 0.0;
        }
        // legacy lid_setOldGroupState
        for t in 0..self.lid_.num_groups() {
            let count = self.lid_.group(t).count;
            for u in 0..count {
                let uu = u as usize;
                if self.owns_spatial_source(self.lid_.group(t).subcatch_idx[uu]) {
                    continue;
                }
                let g = self.lid_.group_mut(t);
                g.old_drain_flow[uu] = g.drain_flow[uu];
            }
        }

        // Advance runoff clock (legacy runoff.c:177, 229).
        self.prev_runoff_step_sec_ = (self.new_runoff_ms_ - self.old_runoff_ms_) / 1000.0;
        self.old_runoff_time_ = self.new_runoff_time_;
        self.old_runoff_ms_ = self.new_runoff_ms_;
        // Run-on for THIS step (legacy runoff.c:247-253).
        self.assemble_runon()?;
        // legacy runoff.c:180 currentDate = StartDateTime + (elapsedMsec + 1)/1000 s
        let abs_time = datetime::add_seconds(self.ctx_.options.start_date, (self.old_runoff_ms_ + 1.0) / 1000.0);

        // Gap #26: clear LID drain quality node accumulators each runoff step.
        self.ctx_.nodes.lid_drain_inflow.fill(0.0);
        if self.ctx_.n_pollutants() > 0 || self.ctx_.reactions.n_species() > 0 || self.ctx_.options.water_age || self.ctx_.options.heat_transport {
            self.ctx_.nodes.lid_drain_qual_load_old = self.ctx_.nodes.lid_drain_qual_load.clone();
            msx_roll_lid_drain(&mut self.ctx_);
            self.ctx_.nodes.lid_drain_qual_vol_old = self.ctx_.nodes.lid_drain_qual_vol.clone();
            self.ctx_.nodes.lid_drain_qual_load.fill(0.0);
            self.ctx_.nodes.lid_drain_qual_vol.fill(0.0);
            self.ctx_.water_age_state.node_lid_drain_age_vol_in.fill(0.0);
            self.ctx_.heat_state.node_lid_drain_temp_vol_in.fill(0.0);
        }

        let mut is_raining = false;

        // A1. Update rain gages and detect rainfall
        gage::update_all_gages(&mut self.ctx_, abs_time);
        for g in 0..self.ctx_.n_gages() {
            if self.ctx_.gages.rainfall[g as usize] > 0.0 {
                is_raining = true;
            }
        }

        // A2. Update climate state
        let doy = datetime::day_of_year(abs_time);
        let mon = datetime::month_of_year(abs_time) - 1;
        let umon = mon as usize;
        let unit_sys = ucf::get_unit_system(self.ctx_.options.flow_units as i32);

        // A2a. Temperature source
        if self.ctx_.climate_state.temp_ts_index >= 0 {
            let idx = self.ctx_.climate_state.temp_ts_index as usize;
            let v = table_lookup_cursor(&mut self.ctx_.tables.tables[idx], abs_time);
            self.ctx_.climate_state.temperature_src = v;
            if unit_sys == 1 {
                self.ctx_.climate_state.temperature_src = (9. / 5.) * self.ctx_.climate_state.temperature_src + 32.0;
            }
            self.ctx_.climate_state.temperature_src += self.ctx_.climate_state.adjust_temp[umon];
        } else if self.ctx_.options.temp_source == 2 && self.climate_file_.is_open() {
            let mut rec = DailyClimateRecord::default();
            if self.climate_file_.get_record(abs_time, &mut rec) {
                let cs = &mut self.ctx_.climate_state;
                if !rec.tmin.is_nan() {
                    cs.file_tmin = rec.tmin;
                }
                if !rec.tmax.is_nan() {
                    cs.file_tmax = rec.tmax;
                }
                {
                    let tmin = cs.file_tmin + cs.adjust_temp[umon];
                    let tmax = cs.file_tmax + cs.adjust_temp[umon];
                    cs.temp_range = tmax - tmin;

                    if doy != cs.last_temp_doy {
                        cs.prev_tmax = if cs.has_minmax { cs.tmax_daily } else { tmax };
                        cs.tmin_daily = tmin;
                        cs.tmax_daily = tmax;
                        cs.has_minmax = true;
                        cs.last_temp_doy = doy;
                        climate::update_temp_times(cs, doy);
                    } else {
                        cs.tmin_daily = tmin;
                        cs.tmax_daily = tmax;
                    }

                    let hour = (abs_time - abs_time.floor()) * 24.0;
                    cs.temperature_src = climate::get_subdaily_temp(cs, hour);
                }
            }
        }

        // A2a'. Climate temperature forcing
        self.ctx_.climate_state.temperature = self.ctx_.forcing.effective_temperature(self.ctx_.climate_state.temperature_src);

        climate::update_daily_climate(&mut self.ctx_.climate_state, doy, mon);

        // A2b'. Seasonal snowmelt coefficients (legacy climate.c:1176-1180, once a day)
        if doy != self.last_melt_doy_ {
            self.last_melt_doy_ = doy;
            self.snow_.set_melt_coeffs(doy);
        }

        // A2b. Evaporation from timeseries or climate file (legacy setEvap)
        if self.ctx_.climate_state.evap_method == EvapMethod::Timeseries && self.ctx_.climate_state.evap_ts_index >= 0 {
            let cs = &mut self.ctx_.climate_state;
            if abs_time >= cs.next_evap_date {
                cs.evap_rate = cs.next_evap_rate / ucf::UCF[ucf::EVAPRATE as usize][unit_sys as usize];
            }
            cs.evap_rate += cs.adjust_evap[umon];
        } else if self.ctx_.climate_state.evap_method == EvapMethod::Pan && self.climate_file_.is_open() {
            let mut rec = DailyClimateRecord::default();
            if self.climate_file_.get_record(abs_time, &mut rec) && !rec.evap.is_nan() {
                let cs = &mut self.ctx_.climate_state;
                cs.evap_rate = rec.evap / ucf::UCF[ucf::EVAPRATE as usize][unit_sys as usize];
                cs.evap_rate *= self.ctx_.options.pan_coeff[umon];
                cs.evap_rate += cs.adjust_evap[umon];
            }
        }
        self.set_next_evap_date(abs_time);

        // A2b'. System-wide evaporation forcing
        self.ctx_.climate_state.evap_rate = self.ctx_.forcing.effective_climate_evap(self.ctx_.climate_state.evap_rate);

        // A2c. Wind speed lookup
        if self.ctx_.options.wind_type == 0 {
            self.ctx_.climate_state.wind_speed_src = self.ctx_.options.wind_speed[umon] / ucf::ucf(ucf::WINDSPEED, &self.ctx_.options);
        } else if self.ctx_.options.wind_type == 1 && self.climate_file_.is_open() {
            let mut rec = DailyClimateRecord::default();
            if self.climate_file_.get_record(abs_time, &mut rec) && !rec.wind.is_nan() {
                self.ctx_.climate_state.wind_speed_src = rec.wind;
            }
        }
        self.ctx_.climate_state.wind_speed = self.ctx_.forcing.effective_wind(self.ctx_.climate_state.wind_speed_src);

        // A2c''. H2: relative humidity
        {
            let mut h = self.ctx_.options.humidity[umon];
            if self.ctx_.climate_state.humidity_ts_index >= 0 {
                let idx = self.ctx_.climate_state.humidity_ts_index as usize;
                h = table_lookup_cursor(&mut self.ctx_.tables.tables[idx], abs_time);
            }
            if self.ctx_.options.humidity_var == 1 {
                let td_c = if unit_sys == 1 { h } else { (h - 32.0) * 5.0 / 9.0 };
                let ta_c = (self.ctx_.climate_state.temperature - 32.0) * 5.0 / 9.0;
                h = 100.0 * saturation_vapour_pressure(td_c) / saturation_vapour_pressure(ta_c);
            }
            self.ctx_.climate_state.humidity = cclamp(h, 0.0, 100.0);
        }

        // The hourly meteorology source precedes legacy temperature, wind and humidity.
        if !self.ctx_.snowpacks.meteorology.is_empty() {
            let mut met = [0.0f64; 5];
            let error = self.ctx_.snowpacks.meteorology.average(abs_time, 1.0, &mut met);
            if !error.is_empty() {
                return Err(thrown(&error));
            }
            self.ctx_.climate_state.temperature = self.ctx_.forcing.effective_temperature(met[0] * 1.8 + 32.0);
            self.ctx_.climate_state.wind_speed = self.ctx_.forcing.effective_wind(met[2] / 0.44704);
            self.ctx_.climate_state.humidity = met[1];
            let ta = self.ctx_.climate_state.temperature;
            self.ctx_.climate_state.ea = 8.1175e6 * (-7701.544 / (ta + 405.0265)).exp(); // LIBM
        }

        // A2d. Monthly adjustment factors
        self.ctx_.climate_state.infil_factor = self.ctx_.adjust_hydcon[umon];

        // A2e. Recovery pattern lookup
        if self.ctx_.climate_state.recovery_pat_index >= 0 {
            let ui = self.ctx_.climate_state.recovery_pat_index as usize;
            if ui < self.ctx_.patterns.factors.len() {
                let facs = &self.ctx_.patterns.factors[ui];
                self.ctx_.climate_state.recovery_factor = if umon < facs.len() { facs[umon] } else { 1.0 };
            }
        }

        // The 2D surface integrates exactly the rates this runoff step applies, from its start,
        // so the mesh and the subcatchments receive the same gage depth.
        #[cfg(feature = "twod")]
        if self.surface_router_.is_active() {
            let start = self.old_runoff_ms_ / 1000.0;
            self.surface_router_.book_gage_rates(&self.ctx_, start);
        }

        // A2z. [FILES] USE RUNOFF
        if let Some(f) = self.runoff_iface_file_.as_ref() {
            if f.is_open() && self.ctx_.files.runoff_mode == FileMode::Use {
                return Err(stub("RunoffInterfaceFile::readResults (USE RUNOFF)"));
            }
        }

        // Compute variable runoff timestep
        let mut dt_runoff = self.compute_runoff_timestep(abs_time, is_raining, self.has_runoff_, self.has_snow_);
        if dt_runoff <= 0.0 {
            dt_runoff = 1.0;
        }
        dt_runoff = crate::hydrology::ueb_subcatchments::next_snow_boundary(&self.ctx_, self.old_runoff_time_, dt_runoff);

        // Update runoff clock (legacy runoff.c:230-237)
        self.new_runoff_time_ = self.old_runoff_time_ + dt_runoff;
        self.new_runoff_ms_ += 1000.0 * dt_runoff;
        let total_ms = self.ctx_.options.total_duration_ms();
        let total_sec = total_ms / 1000.0;
        if self.new_runoff_time_ > total_sec {
            dt_runoff = (total_ms - self.old_runoff_ms_) / 1000.0;
            self.new_runoff_time_ = total_sec;
        }
        if self.new_runoff_ms_ > total_ms {
            self.new_runoff_ms_ = total_ms;
        }

        // Legacy runoff.c:241
        gage::update_past_rain(&mut self.ctx_, dt_runoff as i32);

        // RDII convolution at the wet-weather step
        if !self.ctx_.options.ignore_rdii && !self.ctx_.options.ignore_rainfall {
            if self.rdii_iface_file_.is_open() && !self.rdii_iface_file_.is_writing() {
                // [FILES] USE RDII: flows come from the interface file in stepRouting
            } else {
                let t = self.new_runoff_time_;
                self.rdii_.advance(&mut self.ctx_, t);
                // [FILES] SAVE RDII: export the freshly computed flows (legacy saveRdiiFlows)
                if self.rdii_iface_file_.is_open() && self.rdii_iface_file_.is_writing() {
                    self.rdii_iface_file_.save_flows(abs_time, self.rdii_.node_flows());
                }
            }
        }

        // A3. Snowmelt
        if !self.ctx_.options.ignore_snow_melt {
            self.snow_split_and_melt(dt_runoff)?;
        }

        // A4. Runoff
        let owners = self.ordinary_source_owners().cloned();
        let (evap, infil_f, rec_f) = (self.ctx_.climate_state.evap_rate, self.ctx_.climate_state.infil_factor, self.ctx_.climate_state.recovery_factor);
        if let Err(what) = self.runoff_.execute(&mut self.ctx_, dt_runoff, evap, infil_f, rec_f, mon, owners.as_ref()) {
            return Err(Terminated { exception: "throw std::invalid_argument", what });
        }

        // A3: subcatchment surface age; H5a: its temperature
        crate::transport::water_age_watershed::route_subcatchment_age(&mut self.ctx_, &self.runoff_.soa_, dt_runoff);
        crate::transport::heat_watershed::route_subcatchment_temperature(&mut self.ctx_, &self.runoff_.soa_, dt_runoff).map_err(Terminated::from)?;

        // HasRunoff / HasSnow for the NEXT step's timestep (legacy runoff.c:266-271)
        self.has_runoff_ = false;
        self.has_snow_ = false;
        for i in 0..self.ctx_.n_subcatches() {
            if self.owns_spatial_source(i) {
                continue;
            }
            let usi = i as usize;
            if self.ctx_.subcatches.area[usi] == 0.0 {
                continue;
            }
            let rsoa = self.runoff_.soa();
            if usi < rsoa.subarea_runoff_rate.len() && rsoa.subarea_runoff_rate[usi] > 0.0 {
                self.has_runoff_ = true;
            }
            if self.ctx_.subcatches.snow_depth[usi] > 0.0 {
                self.has_snow_ = true;
            }
        }

        // A4b'. SAVE RUNOFF
        self.save_runoff_iface_step(dt_runoff)?;

        // A4b'. Street sweeping decisions for this step (legacy runoff.c:283)
        self.plan_sweeping(abs_time);

        if !self.ctx_.options.ignore_quality {
            self.step_surface_quality(dt_runoff)?;
            // BW-MSX: the same buildup / washoff arithmetic over the reactions component's
            // species (separate store, same cadence).
            crate::transport::msx_surface::step(&mut self.ctx_, dt_runoff, self.old_runoff_time_)?;
            // S7: the buildup / washoff relations on the 2D cells' land uses (surface quality on
            // the mesh is a phase-2b refusal at initialize, so the store is never active here).
            #[cfg(feature = "twod")]
            if self.surface_router_.surface_quality().borrow().active() {
                return Err(stub("2D surface quality on the mesh (SurfaceQuality2D::step; phase 2b)"));
            }
        }

        // A6. LID performance
        self.lid_inflows(dt_runoff)?;
        self.lid_native_infil(dt_runoff);
        {
            let owners = self.ordinary_source_owners().cloned();
            let evap = self.ctx_.climate_state.evap_rate;
            if let Err(w) = self.lid_.execute(&self.ctx_, dt_runoff, 0.0, evap, owners.as_ref()) {
                return Err(if w.starts_with("Ordinary") { Terminated { exception: "throw std::invalid_argument", what: w } } else { stub(w) });
            }
        }
        // A4: LID layer ages; H5b: temperature column
        if !self.ctx_.options.ignore_quality {
            crate::quality::msx_lid::route_units(&mut self.ctx_, &self.lid_, dt_runoff)?;
        }
        crate::transport::water_age_lid::route_lid_layer_age(&mut self.ctx_, &self.lid_, dt_runoff);
        crate::transport::heat_lid::route_lid_layer_temperature(&mut self.ctx_, &self.lid_, dt_runoff).map_err(Terminated::from)?;

        // A6b. Route LID outputs back (no unit was evaluated above, or the step stopped)
        self.lid_outputs(dt_runoff)?;

        // A6b''. Volume-domain runoff fold for LID subcatchments (subcatch.c:746-751)
        if dt_runoff > 0.0 {
            for i in 0..self.ctx_.n_subcatches() {
                if self.owns_spatial_source(i) {
                    continue;
                }
                let ui = i as usize;
                if ui >= self.ctx_.subcatches.total_lid_area_ft2.len() || self.ctx_.subcatches.total_lid_area_ft2[ui] <= 0.0 {
                    continue;
                }
                let v_lid_out = self.lid_qsurf_cfs_[ui] * dt_runoff;
                let v_outflow = self.runoff_.soa().outflow_vol[ui] - self.lid_vlidin_vol_[ui] + v_lid_out;
                self.ctx_.subcatches.runoff[ui] = v_outflow / dt_runoff;
            }
        }

        // A6b'. Groundwater, after the LID units
        if !self.ctx_.options.ignore_groundwater {
            self.assemble_gw_coupling(dt_runoff);
            self.step_groundwater(dt_runoff)?;
        }

        // A6c. Runoff mass balance
        self.accumulate_runoff_mass_balance(dt_runoff);

        // A6c'. The quality half
        if !self.ctx_.options.ignore_quality {
            crate::quality::msx_lid::publish_units(&mut self.ctx_, dt_runoff);
            self.book_washoff_loads(dt_runoff)?;
        }
        Ok(())
    }

    /// stepRunoff Phase 2: interpolate runoff and GW to the node inflow arrays
    /// (legacy addWetWeatherInflows / addGroundwaterInflows / addLidDrainInflows).
    fn interpolate_runoff_inflows(&mut self) {
        let span = self.new_runoff_ms_ - self.old_runoff_ms_;
        let mut f = if span > 0.0 { (self.ctx_.elapsed_ms - self.old_runoff_ms_) / span } else { 1.0 };
        if f < 0.0 {
            f = 0.0;
        }
        if f > 1.0 {
            f = 1.0;
        }
        self.ctx_.runoff_interp_f = f;

        self.ctx_.nodes.runoff_inflow.fill(0.0);
        self.ctx_.nodes.gw_inflow.fill(0.0);

        let un_sub = self.ctx_.n_subcatches() as usize;
        self.wet_q_interp_.clear();
        self.wet_q_interp_.resize(un_sub, 0.0);
        self.gw_q_interp_.clear();
        self.gw_q_interp_.resize(un_sub, 0.0);
        if self.old_gw_rate_.len() != un_sub {
            self.old_gw_rate_.resize(un_sub, 0.0);
        }
        self.gw_q_node_.clear();
        self.gw_q_node_.resize(un_sub, -1);
        self.lid_drain_q_interp_.clear();
        self.lid_drain_q_interp_.resize(self.ctx_.n_nodes() as usize, 0.0);
        self.lid_drain_q_units_.clear();
        for k in 0..self.lid_.usage_order().len() {
            let (t, u) = self.lid_.usage_order()[k];
            let g = self.lid_.group(t);
            let uu = u as usize;
            let sc = g.subcatch_idx[uu];
            if self.owns_spatial_source(sc) {
                continue;
            }
            if sc < 0 || sc >= self.ctx_.n_subcatches() {
                continue;
            }
            let usc = sc as usize;
            if !(self.ctx_.subcatches.area[usc] > 0.0) || usc >= self.ctx_.subcatches.total_lid_area_ft2.len() || !(self.ctx_.subcatches.total_lid_area_ft2[usc] > 0.0) {
                continue;
            }
            let mut dn = g.drain_node[uu];
            let dsc = g.drain_subcatch[uu];
            if dn < 0 && (dsc < 0 || dsc == sc) {
                dn = self.ctx_.subcatches.outlet_node[usc];
            }
            if dn < 0 || dn >= self.ctx_.n_nodes() {
                continue;
            }
            let q_old = g.old_drain_flow[uu] * g.area[uu];
            let q_new = g.drain_flow[uu] * g.area[uu];
            let q = (1.0 - f) * q_old + f * q_new;
            self.lid_drain_q_interp_[dn as usize] += q;
            self.lid_drain_q_units_.push((dn, q));
        }

        let ucf_land = ucf::ucf(ucf::LANDAREA, &self.ctx_.options);
        for i in 0..self.ctx_.n_subcatches() {
            if self.owns_spatial_source(i) {
                continue;
            }
            let ui = i as usize;
            if self.ctx_.subcatches.area[ui] <= 0.0 {
                continue;
            }
            let q = (1.0 - f) * self.ctx_.subcatches.old_runoff[ui] + f * self.ctx_.subcatches.runoff[ui];

            let out_node = self.ctx_.subcatches.outlet_node[ui];
            if out_node >= 0 && out_node < self.ctx_.n_nodes() {
                self.ctx_.nodes.runoff_inflow[out_node as usize] += q;
                self.wet_q_interp_[ui] = q;
            }

            let mut gw_node = self.ctx_.subcatches.gw_node[ui];
            if gw_node < 0 {
                gw_node = self.ctx_.subcatches.outlet_node[ui];
            }
            if gw_node >= 0 && gw_node < self.ctx_.n_nodes() {
                let area_ft2 = self.ctx_.subcatches.area[ui] / ucf_land;
                let gw_rate_new = if ui < self.groundwater_.state().gw_flow.len() { self.groundwater_.state().gw_flow[ui] } else { 0.0 };
                let gw_rate_old = if ui < self.old_gw_rate_.len() { self.old_gw_rate_[ui] } else { 0.0 };
                let gw_q = ((1.0 - f) * gw_rate_old + f * gw_rate_new) * area_ft2;
                if gw_q.abs() >= crate::core::constants::FLOW_TOL {
                    self.ctx_.nodes.gw_inflow[gw_node as usize] += gw_q;
                    self.gw_q_interp_[ui] = gw_q;
                    self.gw_q_node_[ui] = gw_node;
                }
            }
        }
    }

    /// SWMMEngine::computeRunoffTimestep (legacy runoff_getTimeStep).
    pub(crate) fn compute_runoff_timestep(&self, abs_time: f64, is_raining: bool, has_runoff: bool, has_snow: bool) -> f64 {
        let mut max_step = self.ctx_.options.dry_step as i64;

        {
            let to_evap = datetime::time_diff(self.ctx_.climate_state.next_evap_date, abs_time);
            if to_evap > 0 && to_evap < max_step {
                max_step = to_evap;
            }
        }

        for g in 0..self.ctx_.n_gages() {
            let next_rain_date = gage::gage_next_rain_date(&self.ctx_, g, abs_time);
            let secs_to_change = datetime::time_diff(next_rain_date, abs_time);
            if secs_to_change > 0 && secs_to_change < max_step {
                max_step = secs_to_change;
            }
        }

        let mut time_step = if is_raining || has_snow || has_runoff || self.lid_any_wet() { self.ctx_.options.wet_step as i64 } else { self.ctx_.options.dry_step as i64 };

        if time_step > max_step {
            time_step = max_step;
        }
        time_step as f64
    }

    /// LIDSolver::anyWet(ordinarySourceOwners()).
    fn lid_any_wet(&self) -> bool {
        let owners = self.ordinary_source_owners();
        for g in &self.lid_.groups_ {
            for u in 0..g.count as usize {
                if (owners.is_none() || g.subcatch_idx[u] < 0 || owners.unwrap()[g.subcatch_idx[u] as usize] < 0) && g.is_wet[u] != 0 {
                    return true;
                }
            }
        }
        false
    }

    /// SWMMEngine::subcatchSnowDepth: the subcatchment's area-weighted snow depth (ft).
    pub(crate) fn subcatch_snow_depth(&self, idx: i32) -> f64 {
        let us = idx as usize;
        if self.ctx_.subcatches.snowpack[us] < 0 {
            return 0.0;
        }
        let soa = self.snow_.state();
        let fi = self.ctx_.subcatches.frac_imperv[us];
        let sn = if us < soa.snn.len() { soa.snn[us] } else { 0.0 };
        let f_area = [sn * fi, (1.0 - sn) * fi, 1.0 - fi];
        let base = idx * snow::N_SUBAREAS;
        let mut sd = 0.0;
        for (k, fa) in f_area.iter().enumerate() {
            let uk = (base + k as i32) as usize;
            if uk < soa.wsnow.len() {
                sd += soa.wsnow[uk] * fa;
            }
        }
        sd
    }

    /// stepRunoff A3: the per-subcatchment rain/snow split, plowing, melt, the
    /// snow-modified net precipitation and the snow depth snapshot.
    fn snow_split_and_melt(&mut self, dt_runoff: f64) -> Result<(), Terminated> {
        let un_sc = self.ctx_.n_subcatches() as usize;
        self.snow_rain_.clear();
        self.snow_rain_.resize(un_sc, 0.0);
        self.snow_snow_.clear();
        self.snow_snow_.resize(un_sc, 0.0);
        let ucf_rain_us = ucf::UCF[ucf::RAINFALL as usize][0];
        for i in 0..self.ctx_.n_subcatches() {
            if self.owns_spatial_source(i) {
                continue;
            }
            let ui = i as usize;
            if self.ctx_.subcatches.snowpack[ui] < 0 {
                continue;
            }
            let p = gage::split_precip(&self.ctx_, ui);
            self.snow_rain_[ui] = p.rainfall;
            if ui < self.ctx_.forcing.subcatch_rainfall_mode.len() && self.ctx_.forcing.subcatch_rainfall_mode[ui] != ForcingMode::None {
                let mut rain_inhr = p.rainfall * ucf_rain_us;
                rain_inhr = self.ctx_.forcing.effective_rainfall(ui, rain_inhr);
                self.snow_rain_[ui] = rain_inhr / ucf_rain_us;
            }
            self.snow_snow_[ui] = self.ctx_.forcing.effective_snowfall(ui, p.snowfall);
        }

        {
            let precip_age = self.ctx_.water_age_config.global_age[WaterAgeSource::Rainfall as usize];
            let snow_soa = self.snow_.state_mut();
            snow_soa.track_age = self.ctx_.options.water_age;
            snow_soa.precip_age = precip_age;
        }

        // Accumulation + plowing BEFORE melt (legacy runoff.c:254), then melt.
        self.snow_.plow_snow(&self.ctx_, dt_runoff, &self.snow_snow_);
        let cs = self.ctx_.climate_state;
        self.snow_.execute(&self.ctx_, dt_runoff, cs.temperature, cs.wind_speed, &self.snow_rain_, &self.snow_snow_, cs.gamma, cs.ea);
        // The EQUILIBRIUM surface temperature kernel is not translated: refused here, before
        // advanceSubcatchments changes any state.
        stub_if(self.ctx_.snowpacks.runtime.active && self.ctx_.snowpacks.config.parameters.temperature == crate::hydrology::ueb_types::TemperatureMethod::Equilibrium, "ueb_reference::equilibriumStep (UEB EQUILIBRIUM surface temperature)")?;
        let snow_date = self.ctx_.options.start_date + self.old_runoff_time_ / 86400.;
        match crate::hydrology::ueb_subcatchments::advance_subcatchments(&mut self.ctx_, self.snow_.state_mut(), snow_date, self.old_runoff_time_, dt_runoff, &self.snow_rain_, &self.snow_snow_) {
            Ok(error) => {
                if !error.is_empty() {
                    return Err(thrown(&error));
                }
            }
            Err(what) => return Err(thrown(what.as_bytes())),
        }

        // A3a (Gap #20): per-subcatch snow-modified net precip
        for i in 0..self.ctx_.n_subcatches() {
            if self.owns_spatial_source(i) {
                continue;
            }
            let ui = i as usize;
            if self.ctx_.subcatches.snowpack[ui] < 0 {
                continue;
            }
            let soa = self.snow_.state();
            let rainfall_ft = self.snow_rain_[ui];
            let plow_idx = (i * snow::N_SUBAREAS + snow::SNOW_PLOWABLE) as usize;
            let imperv_idx = (i * snow::N_SUBAREAS + snow::SNOW_IMPERV) as usize;
            let perv_idx = (i * snow::N_SUBAREAS + snow::SNOW_PERV) as usize;

            let f_plow = soa.fArea[plow_idx];
            let f_imperv = soa.fArea[imperv_idx];
            let f_total_i = self.ctx_.subcatches.frac_imperv[ui];

            let mut net_plow = soa.imelt[plow_idx] + rainfall_ft * (1.0 - soa.asc[plow_idx]);
            let mut net_imperv = soa.imelt[imperv_idx] + rainfall_ft * (1.0 - soa.asc[imperv_idx]);
            let mut net_perv = soa.imelt[perv_idx] + rainfall_ft * (1.0 - soa.asc[perv_idx]);
            if soa.external[ui] {
                net_plow = soa.external_net[plow_idx];
                net_imperv = soa.external_net[imperv_idx];
                net_perv = soa.external_net[perv_idx];
            }

            self.ctx_.subcatches.snow_net_imperv[ui] = if f_total_i > 0.0 { (net_plow * f_plow + net_imperv * f_imperv) / f_total_i } else { net_plow };
            self.ctx_.subcatches.snow_net_perv[ui] = net_perv;

            self.ctx_.subcatches.snow_melt_imperv[ui] = if f_total_i > 0.0 { (soa.imelt[plow_idx] * f_plow + soa.imelt[imperv_idx] * f_imperv) / f_total_i } else { soa.imelt[plow_idx] };
            self.ctx_.subcatches.snow_melt_perv[ui] = soa.imelt[perv_idx];

            if soa.track_age {
                let m_plow = soa.imelt[plow_idx] * f_plow;
                let m_imp = soa.imelt[imperv_idx] * f_imperv;
                let m_sum = m_plow + m_imp;
                self.ctx_.subcatches.snow_melt_age_imperv[ui] = if m_sum > 0.0 { (soa.out_age[plow_idx] * m_plow + soa.out_age[imperv_idx] * m_imp) / m_sum } else { soa.out_age[plow_idx] };
                self.ctx_.subcatches.snow_melt_age_perv[ui] = soa.out_age[perv_idx];
            }
        }

        // Legacy newSnowDepth (subcatch.c:816)
        for i in 0..self.ctx_.n_subcatches() {
            if self.owns_spatial_source(i) {
                continue;
            }
            let d = self.subcatch_snow_depth(i);
            self.ctx_.subcatches.snow_depth[i as usize] = d;
        }
        Ok(())
    }

    /// SWMMEngine::saveRunoffIfaceStep (SAVE RUNOFF).
    fn save_runoff_iface_step(&mut self, _dt: f64) -> Result<(), Terminated> {
        match self.runoff_iface_file_.as_ref() {
            None => Ok(()),
            Some(f) => stub_if(f.is_open() && f.writing_, "RunoffInterfaceFile::saveResults (SAVE RUNOFF)"),
        }
    }

    /// stepRunoff A4b': the street sweeping decisions of this step (legacy surfqual.c).
    fn plan_sweeping(&mut self, abs_time: f64) {
        let np = self.ctx_.n_pollutants();
        let nlu = self.ctx_.n_landuses();
        let nsc = self.ctx_.n_subcatches();
        self.ctx_.subcatches.sweep_now.clear();
        self.ctx_.subcatches.sweep_now.resize((cmax(nsc, 0) as usize) * (cmax(nlu, 0) as usize), 0);
        if (np > 0 || self.ctx_.reactions.surface.active()) && nlu > 0 && !self.ctx_.options.ignore_quality {
            let sweep_doy = datetime::day_of_year(abs_time);
            let ss = self.ctx_.options.sweep_start;
            let se = self.ctx_.options.sweep_end;
            let in_season = if ss <= se { sweep_doy >= ss && sweep_doy <= se } else { sweep_doy >= ss || sweep_doy <= se };
            if self.ctx_.subcatches.sweep_last_date.len() != self.ctx_.subcatches.sweep_last_swept.len() {
                let anchor = if self.old_runoff_ms_ <= 0.0 { self.ctx_.options.start_date } else { abs_time };
                let n = self.ctx_.subcatches.sweep_last_swept.len();
                self.ctx_.subcatches.sweep_last_date.resize(n, 0.0);
                for k in 0..n {
                    self.ctx_.subcatches.sweep_last_date[k] = anchor - self.ctx_.subcatches.sweep_last_swept[k];
                }
            }
            if in_season {
                let nlu_sz = self.ctx_.subcatches.sweep_last_swept.len() as i32 / cmax(nsc, 1);
                const MIN_RUNOFF_FTS: f64 = 2.31481e-8;
                for i in 0..nsc {
                    if self.owns_spatial_source(i) {
                        continue;
                    }
                    let ui = i as usize;
                    if ui < self.ctx_.subcatches.rainfall.len() && self.ctx_.subcatches.rainfall[ui] > MIN_RUNOFF_FTS {
                        continue;
                    }
                    let sp = self.ctx_.subcatches.snowpack[ui];
                    if sp >= 0 {
                        const MIN_TOTAL_DEPTH: f64 = 0.004167;
                        let sn = self.snow_.state();
                        let ueb = ui < sn.external.len() && sn.external[ui];
                        let skip = if ueb { self.ctx_.subcatches.snow_net_imperv[ui] > 0.0 } else { sn.wsnow[ui * snow::N_SUBAREAS as usize + snow::SNOW_PLOWABLE as usize] > MIN_TOTAL_DEPTH };
                        if skip {
                            continue;
                        }
                    }
                    for lu in 0..nlu {
                        let ulu = lu as usize;
                        let interval = self.ctx_.landuses.sweep_interval[ulu];
                        if interval <= 0.0 {
                            continue;
                        }
                        let cov_idx = ui * nlu as usize + ulu;
                        let frac = if cov_idx < self.ctx_.subcatches.coverage.len() { self.ctx_.subcatches.coverage[cov_idx] / 100.0 } else { 0.0 };
                        if frac <= 0.0 {
                            continue;
                        }
                        let sw_idx = ui * (cmax(nlu_sz, nlu) as usize) + ulu;
                        if sw_idx >= self.ctx_.subcatches.sweep_last_swept.len() {
                            continue;
                        }
                        if abs_time - self.ctx_.subcatches.sweep_last_date[sw_idx] < interval {
                            continue;
                        }
                        self.ctx_.subcatches.sweep_last_date[sw_idx] = abs_time;
                        self.ctx_.subcatches.sweep_now[cov_idx] = 1;
                    }
                }
            }
            let n = self.ctx_.subcatches.sweep_last_swept.len();
            for k in 0..n {
                self.ctx_.subcatches.sweep_last_swept[k] = abs_time - self.ctx_.subcatches.sweep_last_date[k];
            }
        }
    }

    /// stepRunoff A6a: each LID unit's inflow from non-LID subarea runoff and rainfall
    /// (legacy lid_getRunoff), in [LID_USAGE] order.
    fn lid_inflows(&mut self, dt_runoff: f64) -> Result<(), Terminated> {
        let ns = self.ctx_.n_subcatches() as usize;
        self.lid_vlidin_vol_.clear();
        self.lid_vlidin_vol_.resize(ns, 0.0);
        self.lid_qsurf_cfs_.clear();
        self.lid_qsurf_cfs_.resize(ns, 0.0);
        for k in 0..self.lid_.usage_order().len() {
            let (t, u) = self.lid_.usage_order()[k];
            let uu = u as usize;
            let sc = self.lid_.group(t).subcatch_idx[uu];
            if self.owns_spatial_source(sc) {
                continue;
            }
            if sc < 0 || sc >= self.ctx_.n_subcatches() {
                self.lid_.group_mut(t).inflow[uu] = 0.0;
                continue;
            }
            let usc = sc as usize;
            let rain = self.ctx_.subcatches.rainfall[usc];
            let q_imperv = self.runoff_.soa().imperv_runoff_cfs[usc];
            let q_perv = self.runoff_.soa().perv_runoff_cfs[usc];
            let g = self.lid_.group_mut(t);
            let lid_area = g.area[uu];
            let captured_cfs = q_imperv * g.from_imperv[uu] + q_perv * g.from_perv[uu];
            let q_from_sc = if lid_area > 0.0 { captured_cfs / lid_area } else { 0.0 };
            self.lid_vlidin_vol_[usc] += q_from_sc * lid_area * dt_runoff;
            if g.type_ == crate::hydrology::lid::LIDType::RainBarrel && g.stor_covered[uu] != 0 {
                g.inflow[uu] = q_from_sc;
                if usc < self.ctx_.subcatches.total_lid_area_ft2.len() {
                    self.ctx_.subcatches.lid_return_to_perv_cfs[usc] += rain * self.ctx_.subcatches.total_lid_area_ft2[usc];
                }
            } else {
                g.inflow[uu] = q_from_sc + rain;
            }
            g.inflow[uu] += self.ctx_.subcatches.runon_rate[usc];
            g.subcatch_rain[uu] = rain;
            let q_runon_unit = self.ctx_.subcatches.runon_rate[usc] * lid_area;
            let (fi, fp) = (g.from_imperv[uu], g.from_perv[uu]);
            crate::transport::water_age_lid::set_lid_inflow_age(&mut self.ctx_, t, u, sc, rain, q_imperv * fi, q_perv * fp, q_runon_unit, lid_area);
            crate::transport::heat_lid::set_lid_inflow_temperature(&mut self.ctx_, t, u, sc, rain, q_imperv * fi, q_perv * fp, q_runon_unit, lid_area);
        }
        Ok(())
    }

    /// stepRunoff: legacy findNativeInfil for the LID units (native infiltration, its
    /// aquifer ceiling, the infiltration factors, the runoff clock).
    fn lid_native_infil(&mut self, dt_runoff: f64) {
        let ns_l = self.ctx_.n_subcatches() as usize;
        let mut native = vec![0.0f64; ns_l];
        let mut max_native = vec![1.0e10f64; ns_l];
        let mut factors = vec![1.0f64; ns_l];
        let rec = self.ctx_.climate_state.recovery_factor;
        for i in 0..ns_l as i32 {
            if self.owns_spatial_source(i) {
                continue;
            }
            let ui = i as usize;
            let lid_area = if ui < self.ctx_.subcatches.total_lid_area_ft2.len() { self.ctx_.subcatches.total_lid_area_ft2[ui] } else { 0.0 };
            if lid_area <= 0.0 {
                continue;
            }
            let non_lid = self.runoff_.soa().area[ui];
            if non_lid > 0.0 && self.ctx_.subcatches.frac_imperv[ui] < 1.0 {
                native[ui] = self.runoff_.soa().infil_vol[ui] / non_lid / dt_runoff;
            } else {
                native[ui] = self.runoff_.native_infil_full_lid(&self.ctx_, i, dt_runoff, rec);
            }
            if !self.ctx_.options.ignore_groundwater && self.ctx_.subcatches.gw_aquifer[ui] >= 0 && self.ctx_.subcatches.gw_max_infil_vol[ui] < 1.0e30 {
                max_native[ui] = self.ctx_.subcatches.gw_max_infil_vol[ui] / dt_runoff;
            }
            factors[ui] = self.runoff_.infil_factor_used(i);
        }
        self.lid_.set_native_infil(native, max_native);
        self.lid_.set_infil_factors(factors, rec);
        self.lid_.set_runoff_time(self.old_runoff_time_);
    }

    /// stepRunoff A6b: the LID units' outputs back to the subcatchments and nodes, in
    /// [LID_USAGE] order (legacy evalLidUnit, lid_addDrainRunon, lid_addDrainInflow and
    /// lid_addDrainLoads).
    fn lid_outputs(&mut self, dt_runoff: f64) -> Result<(), Terminated> {
        let ucf_land = ucf::ucf(ucf::LANDAREA, &self.ctx_.options);
        for k in 0..self.lid_.usage_order().len() {
            let (t, u) = self.lid_.usage_order()[k];
            let uu = u as usize;
            let sc = self.lid_.group(t).subcatch_idx[uu];
            if self.owns_spatial_source(sc) {
                continue;
            }
            if sc < 0 || sc >= self.ctx_.n_subcatches() {
                continue;
            }
            let usc = sc as usize;
            let lid_area = self.lid_.group(t).area[uu];

            let lid_area_sc = if usc < self.ctx_.subcatches.total_lid_area_ft2.len() { self.ctx_.subcatches.total_lid_area_ft2[usc] } else { 0.0 };
            let sub_area_ft2 = self.ctx_.subcatches.area[usc] / ucf_land;
            let can_return = self.lid_.group(t).to_perv[uu] != 0 && sub_area_ft2 > lid_area_sc;
            {
                let g = self.lid_.group(t);
                if !can_return {
                    self.lid_qsurf_cfs_[usc] += g.surface_runoff[uu] * lid_area;
                } else {
                    self.ctx_.subcatches.lid_return_to_perv_cfs[usc] += g.surface_runoff[uu] * lid_area;
                    let dn0 = g.drain_node[uu];
                    let dsc0 = g.drain_subcatch[uu];
                    let to_outlet = (dn0 < 0 && dsc0 < 0) || (dn0 == self.ctx_.subcatches.outlet_node[usc] && dsc0 == self.ctx_.subcatches.outlet_subcatch[usc]);
                    if to_outlet {
                        self.ctx_.subcatches.lid_return_to_perv_cfs[usc] += g.drain_flow[uu] * lid_area;
                        self.lid_.group_mut(t).drain_flow[uu] = 0.0;
                    }
                }
            }
            {
                let gm = self.lid_.group_mut(t);
                gm.drain_open[uu] = if gm.drain_flow[uu] > 0.0 { 1 } else { 0 };
            }
            let g = self.lid_.group(t);
            self.ctx_.subcatches.lid_drain_flow[usc] += g.drain_flow[uu] * lid_area;
            let mut dn = g.drain_node[uu];
            let mut dsc = g.drain_subcatch[uu];
            if dn < 0 && (dsc < 0 || dsc == sc) {
                dn = self.ctx_.subcatches.outlet_node[usc];
                dsc = self.ctx_.subcatches.outlet_subcatch[usc];
            }
            if dn >= 0 {
                let un = dn as usize;
                if un < self.ctx_.nodes.lid_drain_inflow.len() {
                    self.ctx_.nodes.lid_drain_inflow[un] += g.drain_flow[uu] * lid_area;
                }
            } else if dsc >= 0 && dsc != sc {
                let utsc = dsc as usize;
                if utsc < self.ctx_.subcatches.lid_drain_runon_cfs.len() {
                    let q_dr = g.drain_flow[uu] * lid_area;
                    self.ctx_.subcatches.lid_drain_runon_cfs[utsc] += q_dr;
                    let lst = &self.ctx_.lid_layer_state;
                    let ut = t as usize;
                    // A4: its age, so assembleRunon hands the rate and the q*age together
                    if self.ctx_.options.water_age && lst.active() && utsc < self.ctx_.water_age_state.subcatch_lid_drain_age_cfs.len() && ut + 1 < lst.group_offset.len() {
                        let fl = lst.group_offset[ut] + u;
                        if fl >= 0 && fl < lst.n_units {
                            let v = lst.drain_value[fl as usize * lst.n_species as usize + LidSpecies::Age as usize];
                            self.ctx_.water_age_state.subcatch_lid_drain_age_cfs[utsc] += q_dr * v;
                        }
                    }
                    // H5b: the same seam for temperature
                    let lst = &self.ctx_.lid_layer_state;
                    if self.ctx_.options.heat_transport && lst.active() && ut + 1 < lst.group_offset.len() {
                        let fl = lst.group_offset[ut] + u;
                        if fl >= 0 && fl < lst.n_units && utsc < self.ctx_.heat_state.subcatch_lid_drain_temp_cfs.len() {
                            let v = lst.drain_value[fl as usize * lst.n_species as usize + LidSpecies::Temperature as usize];
                            self.ctx_.heat_state.subcatch_lid_drain_temp_cfs[utsc] += q_dr * v;
                        }
                    }
                }
            }

            // Gap #26: LID drain quality routing
            let np_ctx = self.ctx_.n_pollutants();
            let np_lid = g.n_pollutants;
            let np_use = cmin(np_ctx, np_lid);
            let species = self.ctx_.reactions.surface.active();
            if (np_use > 0 || self.ctx_.options.water_age || self.ctx_.options.heat_transport || species) && g.drain_flow[uu] > 0.0 && !self.ctx_.options.ignore_quality {
                let drain_cfs = g.drain_flow[uu] * lid_area;
                let mut dest_node = g.drain_node[uu];
                if dest_node < 0 {
                    let target_sc = if g.drain_subcatch[uu] >= 0 { g.drain_subcatch[uu] } else { sc };
                    if target_sc >= 0 && target_sc < self.ctx_.n_subcatches() {
                        dest_node = self.ctx_.subcatches.outlet_node[target_sc as usize];
                    }
                }
                let mut runon_sc = -1;
                {
                    let mut dn2 = g.drain_node[uu];
                    let mut dsc2 = g.drain_subcatch[uu];
                    if dn2 < 0 && (dsc2 < 0 || dsc2 == sc) {
                        dn2 = self.ctx_.subcatches.outlet_node[usc];
                        dsc2 = self.ctx_.subcatches.outlet_subcatch[usc];
                    }
                    if dn2 < 0 && dsc2 >= 0 && dsc2 != sc && dsc2 < self.ctx_.n_subcatches() {
                        runon_sc = dsc2;
                    }
                }
                let mass_ucf_d = ucf::ucf(ucf::MASS, &self.ctx_.options);
                let mcf_of = |c: &crate::core::context::SimulationContext, p: usize| -> f64 {
                    let mut mcf_d = mass_ucf_d;
                    if p < c.pollutants.units.len() {
                        match c.pollutants.units[p] {
                            MassUnits::UgPerL => mcf_d = mass_ucf_d / 1000.0,
                            MassUnits::CountsPerL => mcf_d = 1.0,
                            _ => {}
                        }
                    }
                    mcf_d
                };
                if runon_sc >= 0 {
                    let urs = runon_sc as usize;
                    let is_runoff_load_r = g.drain_node[uu] >= 0 || g.drain_subcatch[uu] == sc;
                    for p in 0..np_use as usize {
                        let rmvl_idx = uu * np_lid as usize + p;
                        let rmvl = if rmvl_idx < g.drain_rmvl.len() { g.drain_rmvl[rmvl_idx] } else { 0.0 };
                        let src_idx = usc * np_ctx as usize + p;
                        let c_src = if src_idx < self.ctx_.subcatches.conc.len() { self.ctx_.subcatches.conc[src_idx] } else { 0.0 };
                        let ri = urs * np_ctx as usize + p;
                        if ri < self.ctx_.subcatches.lid_drain_runon_qual.len() {
                            self.ctx_.subcatches.lid_drain_runon_qual[ri] += (drain_cfs * c_src * 28.317) * (1.0 - rmvl);
                        }
                        if g.to_perv[uu] == 0 {
                            let mcf_d = mcf_of(&self.ctx_, p);
                            let w = drain_cfs * c_src * dt_runoff * 28.317 * mcf_d;
                            if w > 0.0 {
                                self.ctx_.mass_balance.qual_bmp_removal[p] += rmvl * w;
                                if is_runoff_load_r {
                                    self.ctx_.mass_balance.qual_runoff_load[p] += w * (1.0 - rmvl);
                                }
                            }
                        }
                    }
                    let source: Option<Vec<f64>> = crate::quality::msx_lid::drain_values(&self.ctx_, t, u, g.control_idx[uu]).map(|v| v.to_vec());
                    crate::transport::msx_surface::add_lid_drain_runon(&mut self.ctx_, g.control_idx[uu], sc, runon_sc, drain_cfs, g.to_perv[uu] != 0, is_runoff_load_r, dt_runoff, source.as_deref());
                    dest_node = -1;
                }
                if dest_node >= 0 && dest_node < self.ctx_.n_nodes() {
                    let udn = dest_node as usize;
                    if udn < self.ctx_.nodes.lid_drain_qual_vol.len() {
                        self.ctx_.nodes.lid_drain_qual_vol[udn] += drain_cfs;
                    }
                    // A4: the drain's age as the usual q*age rate (storage layer)
                    if self.ctx_.options.water_age {
                        let lst = &self.ctx_.lid_layer_state;
                        let uts = t as usize;
                        if lst.active() && uts + 1 < lst.group_offset.len() {
                            let flat = lst.group_offset[uts] + u;
                            let ai = flat as usize * lst.n_species as usize + LidSpecies::Age as usize;
                            if flat >= 0 && flat < lst.n_units && udn < self.ctx_.water_age_state.node_lid_drain_age_vol_in.len() {
                                let v = lst.drain_value[ai];
                                self.ctx_.water_age_state.node_lid_drain_age_vol_in[udn] += drain_cfs * v;
                            }
                        }
                    }
                    // H5b: and its temperature, as q*T
                    if self.ctx_.options.heat_transport {
                        let lst = &self.ctx_.lid_layer_state;
                        let uts = t as usize;
                        if lst.active() && uts + 1 < lst.group_offset.len() {
                            let flat = lst.group_offset[uts] + u;
                            let ti = flat as usize * lst.n_species as usize + LidSpecies::Temperature as usize;
                            if flat >= 0 && flat < lst.n_units && udn < self.ctx_.heat_state.node_lid_drain_temp_vol_in.len() {
                                let v = lst.drain_value[ti];
                                self.ctx_.heat_state.node_lid_drain_temp_vol_in[udn] += drain_cfs * v;
                            }
                        }
                    }
                    for p in 0..np_use as usize {
                        let rmvl_idx = uu * np_lid as usize + p;
                        let rmvl = if rmvl_idx < g.drain_rmvl.len() { g.drain_rmvl[rmvl_idx] } else { 0.0 };
                        let src_idx = usc * np_ctx as usize + p;
                        let c_src = if src_idx < self.ctx_.subcatches.conc.len() { self.ctx_.subcatches.conc[src_idx] } else { 0.0 };
                        let nd_idx = udn * np_ctx as usize + p;
                        if nd_idx < self.ctx_.nodes.lid_drain_qual_load.len() {
                            self.ctx_.nodes.lid_drain_qual_load[nd_idx] += drain_cfs * c_src * (1.0 - rmvl);
                        }
                        if g.to_perv[uu] == 0 {
                            let mcf_d = mcf_of(&self.ctx_, p);
                            let w = drain_cfs * c_src * dt_runoff * 28.317 * mcf_d;
                            if w > 0.0 {
                                if p < self.ctx_.mass_balance.qual_bmp_removal.len() {
                                    self.ctx_.mass_balance.qual_bmp_removal[p] += rmvl * w;
                                }
                                let is_runoff_load = g.drain_node[uu] >= 0 || g.drain_subcatch[uu] == sc;
                                if is_runoff_load && p < self.ctx_.mass_balance.qual_runoff_load.len() {
                                    self.ctx_.mass_balance.qual_runoff_load[p] += w * (1.0 - rmvl);
                                }
                            }
                        }
                    }
                    // Multispecies Phase 2 (S3): the species ride the same drain water to the
                    // same node.
                    let source: Option<Vec<f64>> = crate::quality::msx_lid::drain_values(&self.ctx_, t, u, g.control_idx[uu]).map(|v| v.to_vec());
                    crate::transport::msx_surface::add_lid_drain(&mut self.ctx_, g.control_idx[uu], sc, dest_node, drain_cfs, g.to_perv[uu] != 0, g.drain_node[uu] >= 0 || g.drain_subcatch[uu] == sc, dt_runoff, source.as_deref());
                }
            }
        }
        Ok(())
    }

    /// SWMMEngine::assembleRunon: run-on from outfalls, upstream subcatchments and LID
    /// drains into runon_inflow / runon_rate (legacy subcatch_getRunon).
    pub(crate) fn assemble_runon(&mut self) -> Result<(), Terminated> {
        self.ctx_.subcatches.runon_inflow.fill(0.0);
        self.ctx_.subcatches.runon_rate.fill(0.0);
        self.ctx_.subcatches.runon_qual_rate.fill(0.0);
        self.ctx_.reactions.surface.runon_qual_rate.fill(0.0); // msxsurf::clearRunon
        {
            let ws = &mut self.ctx_.water_age_state;
            ws.subcatch_runon_age_vol_in.fill(0.0);
            let hs = &mut self.ctx_.heat_state;
            hs.subcatch_runon_temp_vol_in.fill(0.0);
            hs.subcatch_runon_temp_rate.fill(0.0);
        }
        let nsc = self.ctx_.n_subcatches();
        let land2ft2 = 1.0 / ucf::ucf(ucf::LANDAREA, &self.ctx_.options);

        // legacy subcatch_addRunonFlow
        fn add_runon(e: &mut SwmmEngine, k: i32, q: f64, land2ft2: f64) {
            let uk = k as usize;
            let full_ft2 = e.ctx_.subcatches.area[uk] * land2ft2;
            if full_ft2 <= 0.0 {
                return;
            }
            let nonlid = e.runoff_.soa().area[uk];
            e.ctx_.subcatches.runon_rate[uk] += if nonlid > 0.0 { q / nonlid } else { q / full_ft2 };
            e.ctx_.subcatches.runon_inflow[uk] += q;
        }

        // Gap #28: outfall runon FIRST
        let dt_prev = self.prev_runoff_step_sec_;
        if dt_prev > 0.0 {
            for i in 0..nsc {
                if self.owns_spatial_source(i) {
                    continue;
                }
                let ui = i as usize;
                let vol = self.ctx_.subcatches.outfall_runon_vol[ui];
                if vol > 0.0 {
                    add_runon(self, i, vol / dt_prev, land2ft2);
                    {
                        let ws = &mut self.ctx_.water_age_state;
                        if self.ctx_.options.water_age && ui < ws.subcatch_outfall_age_vol.len() && ui < ws.subcatch_runon_age_vol_in.len() {
                            ws.subcatch_runon_age_vol_in[ui] += ws.subcatch_outfall_age_vol[ui] / dt_prev;
                            ws.subcatch_outfall_age_vol[ui] = 0.0;
                        }
                    }
                    if self.ctx_.options.heat_transport && ui < self.ctx_.heat_state.subcatch_outfall_temp_vol.len() {
                        let tv = self.ctx_.heat_state.subcatch_outfall_temp_vol[ui];
                        crate::transport::heat_watershed::add_runon_temperature_at(&mut self.ctx_, i, vol / dt_prev, tv / vol);
                        self.ctx_.heat_state.subcatch_outfall_temp_vol[ui] = 0.0;
                    }
                    self.ctx_.mass_balance.runoff_runon += vol;
                    self.ctx_.subcatches.outfall_runon_vol[ui] = 0.0;
                }
            }
        }

        // Subcatch-to-subcatch, in subcatchment order
        for i in 0..nsc {
            if self.owns_spatial_source(i) {
                continue;
            }
            let ui = i as usize;
            if self.ctx_.subcatches.area[ui] == 0.0 {
                continue;
            }
            let out_sc = self.ctx_.subcatches.outlet_subcatch[ui];
            if out_sc >= 0 && out_sc < nsc && out_sc != i {
                let q = self.ctx_.subcatches.runoff[ui];
                add_runon(self, out_sc, q, land2ft2);
                {
                    let np_r = self.ctx_.n_pollutants() as usize;
                    for p in 0..np_r {
                        let src = ui * np_r + p;
                        let dst = out_sc as usize * np_r + p;
                        if src < self.ctx_.subcatches.conc_old.len() && dst < self.ctx_.subcatches.runon_qual_rate.len() {
                            self.ctx_.subcatches.runon_qual_rate[dst] += self.ctx_.subcatches.runoff[ui] * self.ctx_.subcatches.conc_old[src] * 28.317;
                        }
                    }
                }
                msx_add_runon(&mut self.ctx_, i, out_sc, q);
                crate::transport::water_age_watershed::add_runon_age(&mut self.ctx_, i, out_sc, q);
                let q_sc = self.ctx_.subcatches.runoff[ui];
                crate::transport::heat_watershed::add_runon_temperature(&mut self.ctx_, i, out_sc, q_sc);
            }
        }

        // Gap #25: LID drain runon
        for i in 0..nsc {
            if self.owns_spatial_source(i) {
                continue;
            }
            let ui = i as usize;
            let q = self.ctx_.subcatches.lid_drain_runon_cfs[ui];
            {
                let np_r = self.ctx_.n_pollutants() as usize;
                for p in 0..np_r {
                    let k = ui * np_r + p;
                    if k < self.ctx_.subcatches.lid_drain_runon_qual.len() && k < self.ctx_.subcatches.runon_qual_rate.len() {
                        self.ctx_.subcatches.runon_qual_rate[k] += self.ctx_.subcatches.lid_drain_runon_qual[k];
                        self.ctx_.subcatches.lid_drain_runon_qual[k] = 0.0;
                    }
                }
                msx_consume_lid_drain_runon(&mut self.ctx_, i);
            }
            if q > 0.0 {
                add_runon(self, i, q, land2ft2);
                self.ctx_.subcatches.lid_drain_runon_cfs[ui] = 0.0;
                {
                    let ws = &mut self.ctx_.water_age_state;
                    if self.ctx_.options.water_age && ui < ws.subcatch_lid_drain_age_cfs.len() && ui < ws.subcatch_runon_age_vol_in.len() {
                        ws.subcatch_runon_age_vol_in[ui] += ws.subcatch_lid_drain_age_cfs[ui];
                        ws.subcatch_lid_drain_age_cfs[ui] = 0.0;
                    }
                }
                if self.ctx_.options.heat_transport && ui < self.ctx_.heat_state.subcatch_lid_drain_temp_cfs.len() {
                    let tq = self.ctx_.heat_state.subcatch_lid_drain_temp_cfs[ui];
                    crate::transport::heat_watershed::add_runon_temperature_at(&mut self.ctx_, i, q, tq / q);
                    self.ctx_.heat_state.subcatch_lid_drain_temp_cfs[ui] = 0.0;
                }
            }
        }
        Ok(())
    }

    /// SWMMEngine::assembleGWCoupling: each aquifer's surface water head and available
    /// node flow from the routing state.
    pub(crate) fn assemble_gw_coupling(&mut self, dt_runoff: f64) {
        let ns = self.ctx_.n_subcatches();
        for i in 0..ns {
            let ui = i as usize;
            self.ctx_.subcatches.gw_sw_head[ui] = 0.0;
            self.ctx_.subcatches.gw_node_avail_flow[ui] = 0.0;

            let aq_idx = self.ctx_.subcatches.gw_aquifer[ui];
            if aq_idx < 0 {
                continue;
            }
            let mut gw_node = self.ctx_.subcatches.gw_node[ui];
            if gw_node < 0 {
                gw_node = self.ctx_.subcatches.outlet_node[ui];
            }
            if gw_node >= 0 && gw_node < self.ctx_.n_nodes() {
                let un = gw_node as usize;
                let bottom_elev = self.groundwater_.state().bottom_elev[ui];
                let fixed_depth = self.ctx_.subcatches.gw_tw[ui] / ucf::ucf(ucf::LENGTH, &self.ctx_.options);
                if fixed_depth > 0.0 {
                    self.ctx_.subcatches.gw_sw_head[ui] = fixed_depth + self.ctx_.nodes.invert_elev[un] - bottom_elev;
                } else {
                    self.ctx_.subcatches.gw_sw_head[ui] = self.ctx_.nodes.depth[un] + self.ctx_.nodes.invert_elev[un] - bottom_elev;
                }

                let area_ft2 = self.ctx_.subcatches.area[ui] / ucf::ucf(ucf::LANDAREA, &self.ctx_.options);
                if area_ft2 > 0.0 && dt_runoff > 0.0 {
                    self.ctx_.subcatches.gw_node_avail_flow[ui] = (self.ctx_.nodes.inflow[un] + self.ctx_.nodes.volume[un] / dt_runoff) / area_ft2;
                }
            }
        }
    }

    /// SWMMEngine::stepGroundwater: the aquifer inputs as legacy subcatch_getRunoff hands
    /// them to gwater_getGroundwater, the GW solver, then the flows and statistics.
    pub(crate) fn step_groundwater(&mut self, dt_runoff: f64) -> Result<(), Terminated> {
        let ns = self.ctx_.n_subcatches();
        let uns = ns as usize;
        let ucf_landarea = ucf::ucf(ucf::LANDAREA, &self.ctx_.options);

        self.gw_frac_perv_.clear();
        self.gw_frac_perv_.resize(uns, 0.0);
        self.gw_perv_evap_.clear();
        self.gw_perv_evap_.resize(uns, 0.0);
        self.gw_infil_rate_.clear();
        self.gw_infil_rate_.resize(uns, 0.0);
        {
            self.gw_lid_perv_area_.clear();
            self.gw_lid_perv_area_.resize(uns, 0.0);
            self.gw_lid_infil_vol_.clear();
            self.gw_lid_infil_vol_.resize(uns, 0.0);
            self.gw_lid_perv_evap_vol_.clear();
            self.gw_lid_perv_evap_vol_.resize(uns, 0.0);
            for t in 0..self.lid_.num_groups() {
                let g = self.lid_.group(t);
                for u in 0..g.count as usize {
                    let sc = g.subcatch_idx[u];
                    if sc < 0 || sc >= ns {
                        continue;
                    }
                    let usc = sc as usize;
                    let pervious = g.stor_thick[u] == 0.0 || g.stor_ksat[u] > 0.0;
                    self.gw_lid_infil_vol_[usc] += g.infil_loss[u] * g.area[u];
                    if pervious {
                        self.gw_lid_perv_area_[usc] += g.area[u];
                        self.gw_lid_perv_evap_vol_[usc] += g.evap_loss[u] * g.area[u];
                    }
                }
            }
            let rsoa = self.runoff_.soa();
            for ui in 0..uns {
                let total_area = self.ctx_.subcatches.area[ui] / ucf_landarea;
                if total_area <= 0.0 {
                    continue;
                }
                let mut frac_perv = 1.0 - self.ctx_.subcatches.frac_imperv[ui];
                let lid_area = if ui < self.ctx_.subcatches.total_lid_area_ft2.len() { self.ctx_.subcatches.total_lid_area_ft2[ui] } else { 0.0 };
                if lid_area > 0.0 {
                    frac_perv = (frac_perv * (total_area - lid_area) + self.gw_lid_perv_area_[ui]) / total_area;
                    frac_perv = cmin(frac_perv, 1.0);
                }
                self.gw_frac_perv_[ui] = frac_perv;
                let v_pevap = (if ui < rsoa.perv_evap_vol.len() { rsoa.perv_evap_vol[ui] } else { 0.0 }) + self.gw_lid_perv_evap_vol_[ui];
                let v_infil = (if ui < rsoa.infil_vol.len() { rsoa.infil_vol[ui] } else { 0.0 }) + self.gw_lid_infil_vol_[ui];
                self.gw_perv_evap_[ui] = v_pevap / total_area / dt_runoff;
                self.gw_infil_rate_[ui] = v_infil / total_area / dt_runoff;
            }
        }

        // U3 (track I-b): 2D infiltration whose destination is SUBCATCH_AQUIFER recharges the
        // containing subcatchment's aquifer (m3 since the last drain, as ft/s over the full
        // subcatchment area for this runoff step).
        #[cfg(feature = "twod")]
        let mut use_2d_recharge = false;
        #[cfg(feature = "twod")]
        if self.surface_router_.is_active() && dt_runoff > 0.0 {
            let mut vols = std::mem::take(&mut self.gw_2d_recharge_vol_);
            self.surface_router_.drain_subcatch_recharge(&mut vols);
            self.gw_2d_recharge_vol_ = vols;
            if !self.gw_2d_recharge_vol_.is_empty() {
                self.gw_infil_with_2d_.clone_from(&self.gw_infil_rate_);
                const K_M3_TO_FT3: f64 = 1.0 / (0.3048 * 0.3048 * 0.3048);
                let mut i = 0usize;
                while (i as i32) < ns && i < self.gw_2d_recharge_vol_.len() {
                    let vol_m3 = self.gw_2d_recharge_vol_[i];
                    if vol_m3 > 0.0 {
                        let area_ft2 = self.ctx_.subcatches.area[i] / ucf_landarea;
                        if area_ft2 > 0.0 {
                            self.gw_infil_with_2d_[i] += (vol_m3 * K_M3_TO_FT3) / (area_ft2 * dt_runoff);
                            self.ctx_.mass_balance.gw_infil_2d_recharge += vol_m3 * K_M3_TO_FT3;
                            self.ctx_.mass_balance_2d.infil_to_aquifer += vol_m3;
                        }
                    }
                    i += 1;
                }
                use_2d_recharge = true;
            }
        }

        // Legacy dates the monthly aquifer evaporation pattern on NewRunoffTime.
        let gw_evap_month = datetime::month_of_year(datetime::add_seconds(self.ctx_.options.start_date, (self.new_runoff_ms_ + 1.0) / 1000.0));
        let sw_head = self.ctx_.subcatches.gw_sw_head.clone();
        let evap = self.ctx_.climate_state.evap_rate;
        #[cfg(not(feature = "twod"))]
        let gw_infil = &self.gw_infil_rate_;
        #[cfg(feature = "twod")]
        let gw_infil = if use_2d_recharge { &self.gw_infil_with_2d_ } else { &self.gw_infil_rate_ };
        if let Err(w) = self.groundwater_.execute(&mut self.ctx_, dt_runoff, evap, gw_infil, &sw_head, &self.gw_frac_perv_, &self.gw_perv_evap_, gw_evap_month) {
            return Err(stub(w));
        }

        for ui in 0..uns {
            let gw_rate = self.groundwater_.state().gw_flow[ui];
            let area_ft2 = self.ctx_.subcatches.area[ui] / ucf::ucf(ucf::LANDAREA, &self.ctx_.options);
            self.ctx_.subcatches.gw_flow[ui] = gw_rate * area_ft2;
            if self.ctx_.subcatches.gw_aquifer[ui] >= 0 {
                self.ctx_.subcatches.gw_max_infil_vol[ui] = self.groundwater_.state().max_infil_vol[ui];
            }
        }

        let gws = self.groundwater_.state();
        let sc = &mut self.ctx_.subcatches;
        for ui in 0..uns {
            if sc.gw_aquifer[ui] < 0 {
                continue;
            }
            let area_ft2 = sc.area[ui] / ucf::ucf(ucf::LANDAREA, &self.ctx_.options);
            sc.stat_gw_infil_vol[ui] += self.gw_infil_rate_[ui] * area_ft2 * dt_runoff;
            sc.stat_gw_upper_evap_vol[ui] += gws.upper_evap[ui] * area_ft2 * dt_runoff;
            sc.stat_gw_lower_evap_vol[ui] += gws.lower_evap[ui] * area_ft2 * dt_runoff;
            sc.stat_gw_deep_perc_vol[ui] += gws.deep_loss[ui] * area_ft2 * dt_runoff;
            let flow_cfs = sc.gw_flow[ui];
            sc.stat_gw_flow_vol[ui] += flow_cfs * dt_runoff;
            if flow_cfs > sc.stat_gw_max_flow[ui] {
                sc.stat_gw_max_flow[ui] = flow_cfs;
            }
            let water_table = gws.lower_depth[ui] + gws.bottom_elev[ui];
            sc.stat_gw_sum_theta[ui] += gws.theta[ui] * dt_runoff;
            sc.stat_gw_sum_depth[ui] += water_table * dt_runoff;
            sc.stat_gw_final_theta[ui] = gws.theta[ui];
            sc.stat_gw_final_depth[ui] = water_table;
            sc.stat_gw_time[ui] += dt_runoff;
            sc.stat_gw_steps[ui] += 1;
        }
        Ok(())
    }

    /// SWMMEngine::stepSurfaceQuality: buildup and washoff (surfacespecies::step over the
    /// pollutant family).
    pub(crate) fn step_surface_quality(&mut self, dt_runoff: f64) -> Result<(), Terminated> {
        if self.ctx_.n_pollutants() > 0 {
            let mcf = sskernel::pollutant_surface_mcf(&self.ctx_);
            sskernel::step_pollutants(&mut self.ctx_, &self.landuse_solver_, &mut self.surface_quality_.buildup, &mcf, dt_runoff, self.old_runoff_time_);
        }
        Ok(())
    }

    /// SWMMEngine::bookWashoffLoads: the washoff mass bookings against the post-LID volumes.
    fn book_washoff_loads(&mut self, dt_runoff: f64) -> Result<(), Terminated> {
        let np = self.ctx_.n_pollutants();
        let ns = self.ctx_.n_subcatches();
        let species = self.ctx_.reactions.surface.active();
        if (np <= 0 && !species) || ns <= 0 || self.ctx_.options.ignore_quality {
            return Ok(());
        }
        let mut lid_infil_vol = vec![0.0f64; ns as usize];
        for t in 0..self.lid_.num_groups() {
            let g = self.lid_.group(t);
            for u in 0..g.count as usize {
                let sc = g.subcatch_idx[u];
                if sc < 0 || sc >= ns {
                    continue;
                }
                lid_infil_vol[sc as usize] += g.infil_loss[u] * g.area[u];
            }
        }
        if np > 0 {
            let mcf = sskernel::pollutant_surface_mcf(&self.ctx_);
            sskernel::book_post_lid_pollutants(&mut self.ctx_, &mcf, dt_runoff, &lid_infil_vol);
        }
        // Multispecies Phase 2: the species' post-LID bookings, same volumes.
        if species {
            crate::transport::msx_surface::book_post_lid(&mut self.ctx_, dt_runoff, &lid_infil_vol);
        }
        Ok(())
    }

    /// SWMMEngine::accumulateRunoffMassBalance: the runoff ledger for one runoff step.
    pub(crate) fn accumulate_runoff_mass_balance(&mut self, dt_runoff: f64) {
        let rain_to_ftsec = 1.0 / ucf::ucf(ucf::RAINFALL, &self.ctx_.options);
        let _ = rain_to_ftsec;
        let landarea_to_ft2 = 1.0 / ucf::ucf(ucf::LANDAREA, &self.ctx_.options);
        let ucf_land = ucf::ucf(ucf::LANDAREA, &self.ctx_.options);

        for i in 0..self.ctx_.n_subcatches() {
            if self.owns_spatial_source(i) {
                continue;
            }
            let ui = i as usize;
            let c = &mut self.ctx_;
            let area_ft2 = c.subcatches.area[ui] * landarea_to_ft2;
            let nonlid_area_ft2 = self.runoff_.soa().area[ui];

            let rain_ftsec = c.subcatches.rainfall[ui];
            c.mass_balance.runoff_rainfall += rain_ftsec * area_ft2 * dt_runoff;
            c.mass_balance.runoff_evap += c.subcatches.evap_loss[ui] * nonlid_area_ft2 * dt_runoff;
            c.mass_balance.runoff_infil += c.subcatches.infil_loss[ui] * nonlid_area_ft2 * dt_runoff;

            let sheds_to_node = c.subcatches.outlet_node[ui] >= 0;
            let sheds_to_self = c.subcatches.outlet_subcatch[ui] == i;
            if sheds_to_node || sheds_to_self {
                c.mass_balance.runoff_runoff += c.subcatches.runoff[ui] * dt_runoff;
            }

            if c.subcatches.runoff[ui] > c.subcatches.stat_max_runoff[ui] {
                c.subcatches.stat_max_runoff[ui] = c.subcatches.runoff[ui];
            }
            c.subcatches.stat_precip_vol[ui] += rain_ftsec * (c.subcatches.area[ui] / ucf_land) * dt_runoff;
            c.subcatches.stat_runon_vol[ui] += c.subcatches.runon_rate[ui] * (c.subcatches.area[ui] / ucf_land) * dt_runoff;
            c.subcatches.stat_runoff_vol[ui] += c.subcatches.runoff[ui] * dt_runoff;
        }

        for t in 0..self.lid_.num_groups() {
            let g = self.lid_.group(t);
            for u in 0..g.count as usize {
                if self.owns_spatial_source(g.subcatch_idx[u]) {
                    continue;
                }
                self.ctx_.mass_balance.runoff_evap += g.evap_loss[u] * g.area[u];
                self.ctx_.mass_balance.runoff_infil += g.infil_loss[u] * g.area[u];
            }
        }

        let mut lid_drain_cfs = 0.0;
        for &q in &self.ctx_.nodes.lid_drain_inflow {
            lid_drain_cfs += q;
        }
        self.ctx_.mass_balance.runoff_lid_drain += lid_drain_cfs * dt_runoff;
    }
}

/// msxsurf::rollLidDrain (quality/MsxSurfaceQuality.cpp), active or not.
fn msx_roll_lid_drain(ctx: &mut crate::core::context::SimulationContext) {
    let s = &mut ctx.reactions.surface;
    s.lid_drain_load_old = s.lid_drain_load.clone();
    s.lid_drain_load.fill(0.0);
}

/// msxsurf::addRunon (quality/MsxSurfaceQuality.cpp).
fn msx_add_runon(ctx: &mut crate::core::context::SimulationContext, donor: i32, receiver: i32, q: f64) {
    let s = &mut ctx.reactions.surface;
    if !s.active() {
        return;
    }
    for m in 0..s.n_species {
        let src = s.sidx(donor, m);
        let dst = s.sidx(receiver, m);
        if src < s.washoff_conc_old.len() && dst < s.runon_qual_rate.len() {
            s.runon_qual_rate[dst] += q * s.washoff_conc_old[src] * 28.317;
        }
    }
}

/// msxsurf::consumeLidDrainRunon (quality/MsxSurfaceQuality.cpp).
fn msx_consume_lid_drain_runon(ctx: &mut crate::core::context::SimulationContext, receiver: i32) {
    let s = &mut ctx.reactions.surface;
    if !s.active() {
        return;
    }
    for m in 0..s.n_species {
        let k = s.sidx(receiver, m);
        if k < s.lid_drain_runon_load.len() && k < s.runon_qual_rate.len() {
            s.runon_qual_rate[k] += s.lid_drain_runon_load[k];
            s.lid_drain_runon_load[k] = 0.0;
        }
    }
}

/// Gate only (cargo feature `gate`): the members the runoff step reads that neither the
/// context loader nor the hydrology loader fills, as gate/diff/ctxstep's runoff mode writes
/// them (SIDE.eng). Not part of the translation.
#[cfg(feature = "gate")]
#[derive(Debug, Clone, Default)]
pub struct GateRunoffSide {
    pub old_runoff_time_: f64,
    pub new_runoff_time_: f64,
    pub old_runoff_ms_: f64,
    pub prev_runoff_step_sec_: f64,
    pub new_runoff_ms_: f64,
    pub has_runoff_: bool,
    pub has_snow_: bool,
    pub wet_q_interp_: Vec<f64>,
    pub lid_drain_q_interp_: Vec<f64>,
    pub lid_drain_q_units_: Vec<(i32, f64)>,
    pub lid_vlidin_vol_: Vec<f64>,
    pub lid_qsurf_cfs_: Vec<f64>,
    pub buildup: Vec<f64>,
}

#[cfg(feature = "gate")]
impl SwmmEngine {
    pub fn gate_install_runoff_side(&mut self, s: GateRunoffSide) {
        self.old_runoff_time_ = s.old_runoff_time_;
        self.new_runoff_time_ = s.new_runoff_time_;
        self.old_runoff_ms_ = s.old_runoff_ms_;
        self.prev_runoff_step_sec_ = s.prev_runoff_step_sec_;
        self.new_runoff_ms_ = s.new_runoff_ms_;
        self.has_runoff_ = s.has_runoff_;
        self.has_snow_ = s.has_snow_;
        self.wet_q_interp_ = s.wet_q_interp_;
        self.lid_drain_q_interp_ = s.lid_drain_q_interp_;
        self.lid_drain_q_units_ = s.lid_drain_q_units_;
        self.lid_vlidin_vol_ = s.lid_vlidin_vol_;
        self.lid_qsurf_cfs_ = s.lid_qsurf_cfs_;
        self.surface_quality_.buildup = s.buildup;
    }

    /// One call of the runoff gate: what step() does before stepRunoff
    /// (save_lat_qual_state, resetStepMassBalance, applyForcings), then stepRunoff.
    pub fn gate_runoff_call(&mut self, dt: f64) -> Result<(), Terminated> {
        self.ctx_.save_lat_qual_state();
        self.reset_step_mass_balance();
        self.apply_forcings(dt)?;
        self.step_runoff(dt)
    }

    pub fn gate_surface_buildup(&self) -> Vec<f64> {
        self.surface_quality_.buildup.clone()
    }

    /// The routing clock advanced by hand between runoff gate calls.
    pub fn gate_advance_clock(&mut self, dt: f64) {
        let total = self.ctx_.options.total_duration_ms();
        self.ctx_.old_elapsed_ms = self.ctx_.elapsed_ms;
        self.ctx_.elapsed_ms += 1000.0 * dt;
        if self.ctx_.elapsed_ms > total {
            self.ctx_.elapsed_ms = total;
        }
    }
}
