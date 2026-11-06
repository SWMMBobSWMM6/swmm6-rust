// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/SWMMEngine.hpp and SWMMEngine.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The engine lifecycle manager, so far: construction, open() (1D build: OPENSWMM_HAS_2D
// off) and the private helpers open() calls (validate_project, set_error, emit_warning,
// push_report_warning, write_open_failure_report). initialize/start/step/end/report/close
// and the solver members they drive come with the solver phases.
//
// Differences of form:
//   * The members exist that open() and the init legs translated so far touch:
//     initHydraulics (core/engine_init_hydraulics.rs: router_, controls_, inflow_, rdii_,
//     exfil_, inlet_, hydstruct_, iface_, next_event_, between_events_), initHydrology
//     (core/engine_init_hydrology.rs: runoff_, climate_file_, snow_, last_melt_doy_,
//     groundwater_, lid_, gage_states_ and the gw_*/snow_* work vectors), and initQuality,
//     initGeometry and initMassBalance (core/engine_init_quality.rs: hotstart_deferred_,
//     quality_, landuse_solver_, surface_quality_). They are declared in upstream's order.
//     The other members (ARD, LAGRANGIAN, interface files, checkpoints, run-time state) are
//     added as their solvers are translated. ctx_ and the solver members are pub(crate) so
//     the init methods in those sibling files can borrow them field by field.
//   * ctx_.state_accessors: upstream installs closures over the engine's solvers
//     (FV dry state, channel thermal/species regions, snow, infiltration and groundwater
//     state). They read and write solver state no open() result depends on, and the
//     solvers are not translated yet, so they are left unset here. NOT YET TRANSLATED.
//   * write_open_failure_report writes the .rpt of a failed open through
//     DefaultReportPlugin::prepare/finalize, which take the context const. The .rpt writer
//     is not translated yet, so no report file is written. NOT YET TRANSLATED.
//   * perf:: timers (reset_load, reset_fv, ScopedTimer, begin_process_model) only time
//     the run; they are left out.
//   * emit_warning's callback is a boxed closure; the C callback's engine pointer and user
//     data are the closure's captures.

use crate::core::cmath::{cmin, cvt_i32};
use crate::core::context::{EngineState, SimulationContext};
use crate::core::error_codes::{format_error, format_warning, ERR_RAIN_GAGE_INTERVAL, ERR_SUBCATCH_OUTLET, WARN_DRY_STEP_INCREASED, WARN_GAGE_INTERVAL, WARN_ROUTING_STEP_REDUCED, WARN_WET_STEP_REDUCED};
use crate::core::fs_path::FsPath;
use crate::core::options::{Bytes, QualitySolverKind};
use crate::core::plugins::{DefaultInputPlugin, DefaultOutputPlugin, DefaultReportPlugin, DefaultStateIOPlugin, PluginFactory};
use crate::data::gage::RainSource;
use crate::data::species_registry::SpeciesKind;
use crate::hydrology::climate_file::ClimateFileReader;
use crate::hydrology::gage::GageState;
use crate::hydrology::groundwater::GWSolver;
use crate::hydrology::lid::LIDSolver;
use crate::hydrology::rdii::RDIISolver;
use crate::hydrology::runoff::RunoffSolver;
use crate::hydrology::snow::SnowSolver;
use crate::quality::landuse::{LanduseSolver, SurfaceQualitySoA};
use crate::quality::quality_routing::QualitySolver;

// Error codes returned across the C ABI (the public SWMM_ErrorCode values).
pub const SWMM_OK: i32 = 0;
pub const SWMM_ERR_MEMORY: i32 = 1;
pub const SWMM_ERR_FILE_NOT_FOUND: i32 = 2;
pub const SWMM_ERR_PARSE: i32 = 5;
pub const SWMM_ERR_WRONG_STATE: i32 = 6;
pub const SWMM_ERR_BADPARAM: i32 = 9;
pub const SWMM_ERR_PLUGIN: i32 = 10;
pub const SWMM_ERR_IO: i32 = 11;
// core/ErrorCodes.hpp CffiErrorCode values the lifecycle returns
pub const CFFI_ERR_LIFECYCLE: i32 = 6;
pub const CFFI_ERR_BADPARAM: i32 = 9;
pub const CFFI_ERR_HOTSTART: i32 = 12;

/// The warning callback of EngineCallbacks (code, message).
pub type WarningCallback = Box<dyn FnMut(i32, &[u8])>;

#[derive(Default)]
pub struct EngineCallbacks {
    pub on_warning: Option<WarningCallback>,
    /// on_progress (fraction complete, current date); the C callback's engine pointer and
    /// user data are the closure's captures
    pub on_progress: Option<Box<dyn FnMut(f64, f64)>>,
    /// on_step_begin (current date, dt)
    pub on_step_begin: Option<Box<dyn FnMut(f64, f64)>>,
    /// on_step_end (current date, dt)
    pub on_step_end: Option<Box<dyn FnMut(f64, f64)>>,
}

/// SWMMEngine::AvgAccumulator: float32 display-unit sums of node and link results between
/// report boundaries (legacy RptFlags.averages).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct AvgAccumulator {
    pub node_depth: Vec<f32>,
    pub node_head: Vec<f32>,
    pub node_volume: Vec<f32>,
    pub node_lat_inflow: Vec<f32>,
    pub node_total_inflow: Vec<f32>,
    pub node_overflow: Vec<f32>,
    pub link_flow: Vec<f32>,
    pub link_depth: Vec<f32>,
    pub link_velocity: Vec<f32>,
    pub link_volume: Vec<f32>,
    pub link_capacity: Vec<f32>,
    /// Number of routing steps accumulated
    pub n_steps: i32,
}

impl AvgAccumulator {
    /// Allocate zeroed report-sum arrays and reset the sample count.
    pub fn resize(&mut self, n_nodes: i32, n_links: i32) {
        let un = n_nodes as usize;
        let ul = n_links as usize;
        for v in [&mut self.node_depth, &mut self.node_head, &mut self.node_volume, &mut self.node_lat_inflow, &mut self.node_total_inflow, &mut self.node_overflow] {
            v.clear();
            v.resize(un, 0.0);
        }
        for v in [&mut self.link_flow, &mut self.link_depth, &mut self.link_velocity, &mut self.link_volume, &mut self.link_capacity] {
            v.clear();
            v.resize(ul, 0.0);
        }
        self.n_steps = 0;
    }

    /// Zero the report-sum arrays and the sample count without resizing.
    pub fn reset(&mut self) {
        for v in [
            &mut self.node_depth,
            &mut self.node_head,
            &mut self.node_volume,
            &mut self.node_lat_inflow,
            &mut self.node_total_inflow,
            &mut self.node_overflow,
            &mut self.link_flow,
            &mut self.link_depth,
            &mut self.link_velocity,
            &mut self.link_volume,
            &mut self.link_capacity,
        ] {
            v.iter_mut().for_each(|x| *x = 0.0);
        }
        self.n_steps = 0;
    }
}

pub struct SwmmEngine {
    pub(crate) ctx_: SimulationContext,
    /// permissive open (see set_lenient_open)
    lenient_open_: bool,
    /// parse-only handle cannot initialize a run
    pub(crate) metadata_only_: bool,
    /// initialize(): a native USE HOTSTART waits for the end of initialize()
    pub(crate) hotstart_deferred_: bool,
    pub(crate) plugins_: PluginFactory,

    // Computational modules (SWMMEngine.hpp order)
    /// Hydraulic routing (owns XSectGroups, the DW/KW solvers and the FV state)
    pub(crate) router_: crate::hydraulics::routing::Router,
    /// Subcatchment runoff (batch nonlinear reservoir)
    pub(crate) runoff_: RunoffSolver,
    /// Climate file reader (temp/evap/wind from file)
    pub(crate) climate_file_: ClimateFileReader,
    /// Snowmelt (batch over subcatch x subareas)
    pub(crate) snow_: SnowSolver,
    /// Day the seasonal melt coeffs were last set (legacy climate.c LastDay)
    pub(crate) last_melt_doy_: i32,
    /// Groundwater (batch ODE per subcatchment)
    pub(crate) groundwater_: GWSolver,
    /// LID (batch by type group)
    pub(crate) lid_: LIDSolver,
    /// Quality routing (batch link-load + mixing)
    pub(crate) quality_: QualitySolver,
    /// Buildup/washoff computation
    pub(crate) landuse_solver_: LanduseSolver,
    /// Per-subcatch surface quality state
    pub(crate) surface_quality_: SurfaceQualitySoA,
    /// Control rule evaluation
    pub(crate) controls_: crate::controls::controls::ControlEngine,
    /// External + DWF inflows
    pub(crate) inflow_: crate::hydrology::inflow::InflowSolver,
    /// RDII (unit hydrograph convolution); initHydraulics initializes it
    pub(crate) rdii_: RDIISolver,
    /// Storage node exfiltration
    pub(crate) exfil_: crate::hydraulics::exfiltration::ExfilSolver,
    /// Street inlet capture
    pub(crate) inlet_: crate::hydraulics::inlet::InletSolver,
    /// Per-subcatch pervious fraction for GW evap
    pub(crate) gw_frac_perv_: Vec<f64>,
    /// Per-subcatch pervious evap rate (ft/sec)
    pub(crate) gw_perv_evap_: Vec<f64>,
    pub(crate) gw_lid_perv_area_: Vec<f64>,
    pub(crate) gw_lid_infil_vol_: Vec<f64>,
    pub(crate) gw_lid_perv_evap_vol_: Vec<f64>,
    /// U3: the 2D surface's per-subcatchment infiltration recharge (m3)
    pub(crate) gw_2d_recharge_vol_: Vec<f64>,
    pub(crate) gw_infil_with_2d_: Vec<f64>,
    /// (Vinfil + VlidInfil) / area / dt, legacy's `infil` (ft/s)
    pub(crate) gw_infil_rate_: Vec<f64>,
    /// Per-subcatch rainfall into snow step (ft/sec)
    pub(crate) snow_rain_: Vec<f64>,
    /// Per-subcatch snowfall into snow step (ft/sec)
    pub(crate) snow_snow_: Vec<f64>,
    /// Pumps, orifices, weirs, outlets
    pub(crate) hydstruct_: crate::hydraulics::hyd_structures::StructureSolver,
    /// Routing interface file I/O
    pub(crate) iface_: crate::core::interface_file::InterfaceManager,
    /// Index of next event in ctx_.events
    pub(crate) next_event_: i32,
    /// legacy BetweenEvents, as the last stepRouting evaluated it
    pub(crate) between_events_: bool,
    /// Per-gage state
    pub(crate) gage_states_: Vec<GageState>,

    /// Report file path
    rpt_path_: Bytes,
    /// Binary output file path
    out_path_: Bytes,
    /// per-subcatch interpolated GW flow (cfs)
    pub(crate) gw_q_interp_: Vec<f64>,
    /// receiving node for gw_q_interp_ (-1 = skip)
    pub(crate) gw_q_node_: Vec<i32>,
    pub(crate) old_gw_rate_: Vec<f64>,
    pub(crate) callbacks_: EngineCallbacks,
    /// index of the built-in .out writer in plugins_.output_plugins() (owned there)
    binary_output_plugin_: Option<usize>,

    // ---- the run lifecycle's members (initialize, start, step: core/engine_initialize.rs,
    // engine_start.rs, engine_step.rs), in upstream's relative order. Declared together
    // here, after the others, rather than interleaved at their upstream positions. They are
    // pub (upstream's are private) so the gate's state dumper can read them.
    /// start() reached prepare_all (report file may be open)
    pub plugins_prepare_attempted_: bool,
    /// per SAVE HOTSTART row: its dated save was written this run
    pub hotstart_save_done_: Vec<u8>,
    /// USE HOTSTART runoff section, applied in init_modules()
    pub hotstart_runoff_: Vec<f64>,
    /// USE HOTSTART species surface block, applied in init_modules()
    pub hotstart_species_: crate::core::hot_start_manager::LegacySpeciesTrailer,
    /// optional runoff interface file (legacy "Frunoff"), opened by start()
    pub runoff_iface_file_: Option<crate::hydrology::runoff_interface::RunoffInterfaceFile>,
    /// RDII interface file (legacy "Frdii"), opened by start()
    pub rdii_iface_file_: crate::hydrology::rdii_interface::RdiiInterfaceFile,
    pub step_active_: bool,
    pub runtime_pipeline_: bool,
    pub checkpoint_enabled_: bool,
    pub checkpoint_ready_: bool,
    pub checkpoint_output_: bool,
    pub checkpoint_model_: Bytes,
    pub pending_hotstart_checkpoint_: Bytes,
    pub hotstart_checkpoint_error_: Bytes,
    pub exchange_target_: f64,
    /// Previous runoff boundary (seconds from start)
    pub old_runoff_time_: f64,
    /// Next runoff boundary (seconds from start)
    pub new_runoff_time_: f64,
    /// legacy OldRunoffTime (msec)
    pub old_runoff_ms_: f64,
    /// span of the runoff step just completed (legacy oldRunoffStep)
    pub prev_runoff_step_sec_: f64,
    /// legacy NewRunoffTime (msec)
    pub new_runoff_ms_: f64,
    /// legacy NewRuleTime (msec, routing.c:58)
    pub new_rule_time_ms_: f64,
    pub prev_lid_infil_vol_: f64,
    pub prev_lid_evap_vol_: f64,
    /// per-subcatch interpolated runoff (cfs)
    pub wet_q_interp_: Vec<f64>,
    /// per-node interpolated LID drain inflow this routing step (cfs)
    pub lid_drain_q_interp_: Vec<f64>,
    /// the same drain flows per unit as (node, q), in legacy's order
    pub lid_drain_q_units_: Vec<(i32, f64)>,
    /// per-subcatch VlidIn this runoff step (ft3)
    pub lid_vlidin_vol_: Vec<f64>,
    /// per-subcatch LID surface outflow to outlet (cfs)
    pub lid_qsurf_cfs_: Vec<f64>,
    /// Prev step generated runoff (legacy HasRunoff)
    pub has_runoff_: bool,
    /// Prev step had snow cover (legacy HasSnow)
    pub has_snow_: bool,
    /// Whether to save binary results
    pub save_results_: i32,
    /// Legacy DoRouting analog (swmm5.c:748), set by start()
    pub do_routing_: bool,
    /// Averaging accumulator (only used when rpt_averages == true)
    pub avg_: AvgAccumulator,
    /// xsp_cache_: the reporting-path XSectParams cache. Upstream's entries borrow the
    /// context's transect tables; here the cache keeps its size only, and its users build
    /// the (verbatim-copy, so bit-identical) params from the link store when they need them.
    pub xsp_cache_len_: usize,
    /// generation the cache was built at
    pub xsp_cache_gen_: u64,
    /// True when the last stepRouting() skipped routing as steady state.
    pub last_step_steady_: bool,
    /// Previous routing step's system flow RATES (cfs) for the eleven routing continuity terms
    pub routing_prev_rates_: [f64; 11],
    /// lazy EULERIAN_ARD init after Router::init (read and set by stepRouting)
    pub ard_init_attempted_: bool,
}

impl Default for SwmmEngine {
    fn default() -> Self {
        Self::new()
    }
}

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut m = Vec::new();
    for p in parts {
        m.extend_from_slice(p);
    }
    m
}

impl SwmmEngine {
    /// The constructor. Upstream registers the four process components once per process
    /// here; ProcessComponentRegistry::instance() does that on first use.
    pub fn new() -> Self {
        let mut e = SwmmEngine {
            ctx_: SimulationContext::default(),
            lenient_open_: false,
            metadata_only_: false,
            hotstart_deferred_: false,
            plugins_: PluginFactory::new(),
            router_: Default::default(),
            runoff_: RunoffSolver::default(),
            climate_file_: ClimateFileReader::default(),
            snow_: SnowSolver::default(),
            last_melt_doy_: -1,
            groundwater_: GWSolver::default(),
            lid_: LIDSolver::default(),
            quality_: QualitySolver::default(),
            landuse_solver_: LanduseSolver::default(),
            surface_quality_: SurfaceQualitySoA::default(),
            controls_: Default::default(),
            inflow_: Default::default(),
            rdii_: RDIISolver::default(),
            exfil_: Default::default(),
            inlet_: Default::default(),
            gw_frac_perv_: Vec::new(),
            gw_perv_evap_: Vec::new(),
            gw_lid_perv_area_: Vec::new(),
            gw_lid_infil_vol_: Vec::new(),
            gw_lid_perv_evap_vol_: Vec::new(),
            gw_2d_recharge_vol_: Vec::new(),
            gw_infil_with_2d_: Vec::new(),
            gw_infil_rate_: Vec::new(),
            snow_rain_: Vec::new(),
            snow_snow_: Vec::new(),
            hydstruct_: Default::default(),
            iface_: Default::default(),
            next_event_: 0,
            between_events_: false,
            gage_states_: Vec::new(),
            rpt_path_: Vec::new(),
            out_path_: Vec::new(),
            gw_q_interp_: Vec::new(),
            gw_q_node_: Vec::new(),
            old_gw_rate_: Vec::new(),
            callbacks_: EngineCallbacks::default(),
            binary_output_plugin_: None,
            plugins_prepare_attempted_: false,
            hotstart_save_done_: Vec::new(),
            hotstart_runoff_: Vec::new(),
            hotstart_species_: Default::default(),
            runoff_iface_file_: None,
            rdii_iface_file_: Default::default(),
            step_active_: false,
            runtime_pipeline_: false,
            checkpoint_enabled_: false,
            checkpoint_ready_: false,
            checkpoint_output_: false,
            checkpoint_model_: Vec::new(),
            pending_hotstart_checkpoint_: Vec::new(),
            hotstart_checkpoint_error_: Vec::new(),
            exchange_target_: -1.0,
            old_runoff_time_: 0.0,
            new_runoff_time_: 0.0,
            old_runoff_ms_: 0.0,
            prev_runoff_step_sec_: 0.0,
            new_runoff_ms_: 0.0,
            new_rule_time_ms_: 0.0,
            prev_lid_infil_vol_: 0.0,
            prev_lid_evap_vol_: 0.0,
            wet_q_interp_: Vec::new(),
            lid_drain_q_interp_: Vec::new(),
            lid_drain_q_units_: Vec::new(),
            lid_vlidin_vol_: Vec::new(),
            lid_qsurf_cfs_: Vec::new(),
            has_runoff_: false,
            has_snow_: false,
            save_results_: 0,
            do_routing_: true,
            avg_: AvgAccumulator::default(),
            xsp_cache_len_: 0,
            xsp_cache_gen_: !0u64,
            last_step_steady_: false,
            routing_prev_rates_: [0.0; 11],
            ard_init_attempted_: false,
        };
        e.ctx_.state = EngineState::Created;
        e
    }

    pub fn context(&self) -> &SimulationContext {
        &self.ctx_
    }
    pub fn context_mut(&mut self) -> &mut SimulationContext {
        &mut self.ctx_
    }
    pub fn plugin_factory(&self) -> &PluginFactory {
        &self.plugins_
    }
    pub fn set_lenient_open(&mut self, on: bool) {
        self.lenient_open_ = on;
    }
    pub fn lenient_open(&self) -> bool {
        self.lenient_open_
    }
    pub fn metadata_only(&self) -> bool {
        self.metadata_only_
    }
    pub fn rpt_path(&self) -> &[u8] {
        &self.rpt_path_
    }
    pub fn out_path(&self) -> &[u8] {
        &self.out_path_
    }
    pub fn binary_output_plugin(&self) -> Option<usize> {
        self.binary_output_plugin_
    }
    pub fn last_error_message(&self) -> &[u8] {
        &self.ctx_.error_message
    }
    pub fn set_warning_callback(&mut self, cb: Option<WarningCallback>) {
        self.callbacks_.on_warning = cb;
    }

    // Read access to the init legs' members, for the gates (no accessor upstream). The
    // hydrology solvers' accessors are in core/engine_init_hydrology.rs.
    pub fn router(&self) -> &crate::hydraulics::routing::Router {
        &self.router_
    }
    pub fn control_engine(&self) -> &crate::controls::controls::ControlEngine {
        &self.controls_
    }
    pub fn inflow_solver(&self) -> &crate::hydrology::inflow::InflowSolver {
        &self.inflow_
    }
    pub fn exfil_solver(&self) -> &crate::hydraulics::exfiltration::ExfilSolver {
        &self.exfil_
    }
    pub fn inlet_solver(&self) -> &crate::hydraulics::inlet::InletSolver {
        &self.inlet_
    }
    pub fn structure_solver(&self) -> &crate::hydraulics::hyd_structures::StructureSolver {
        &self.hydstruct_
    }
    pub fn interface_manager(&self) -> &crate::core::interface_file::InterfaceManager {
        &self.iface_
    }
    pub fn next_event(&self) -> i32 {
        self.next_event_
    }
    pub fn between_events(&self) -> bool {
        self.between_events_
    }
    pub fn hotstart_deferred(&self) -> bool {
        self.hotstart_deferred_
    }
    /// quality_ (no accessor upstream; read by the gate)
    pub fn quality(&self) -> &QualitySolver {
        &self.quality_
    }
    /// landuse_solver_ (no accessor upstream; read by the gate)
    pub fn landuse_solver(&self) -> &LanduseSolver {
        &self.landuse_solver_
    }
    /// surfaceQuality()
    pub fn surface_quality(&self) -> &SurfaceQualitySoA {
        &self.surface_quality_
    }
    /// The three init legs of core/engine_init_quality.rs, for the gate (examples/initqual.rs).
    pub fn gate_init_quality(&mut self) -> Result<(), crate::core::engine_init_hydraulics::Terminated> {
        self.init_quality()
    }
    pub fn gate_init_geometry(&mut self) {
        self.init_geometry();
    }
    pub fn gate_init_mass_balance(&mut self) {
        self.init_mass_balance();
    }

    // ========================================================================
    // open()
    // ========================================================================

    /// `swmm_engine_open` / `swmm_engine_open_metadata`. A null C string is `None`.
    pub fn open(&mut self, inp_path: Option<&[u8]>, rpt_path: Option<&[u8]>, out_path: Option<&[u8]>, input_plugin_lib: Option<&[u8]>, metadata_only: bool) -> i32 {
        if self.ctx_.state != EngineState::Created && self.ctx_.state != EngineState::Closed {
            self.set_error(SWMM_ERR_WRONG_STATE, b"swmm_engine_open: engine is not in CREATED or CLOSED state");
            return SWMM_ERR_WRONG_STATE;
        }

        // Reset context for a fresh run
        self.ctx_.reset();
        self.metadata_only_ = metadata_only;

        // Stamp the report wall clock before any parsing work (std::time).
        self.ctx_.wall_start = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs() as i64).unwrap_or(0);

        self.rpt_path_ = rpt_path.unwrap_or(b"").to_vec();
        self.out_path_ = out_path.unwrap_or(b"").to_vec();
        let inp: &[u8] = inp_path.unwrap_or(b"");

        // Resolve input plugin: path, id:version, or fall back to default
        if self.plugins_.input_plugins().is_empty() {
            if let Some(lib) = input_plugin_lib {
                if !lib.is_empty() {
                    let mut warned: Vec<Bytes> = Vec::new();
                    let info = self.plugins_.find_component(lib, Some(&mut |m: &[u8]| warned.push(m.to_vec())));
                    for m in &warned {
                        self.emit_warning(crate::core::engine::SWMM_ERR_PLUGIN, m);
                    }
                    // A built-in info's create_input_plugin() returns null: nothing added.
                    let _ = info;
                }
            }
            // Fall back to built-in .inp reader
            if self.plugins_.input_plugins().is_empty() {
                let mut ip = DefaultInputPlugin::new();
                ip.initialize(&[]);
                self.plugins_.add_input_plugin(ip);
            }
        }

        if self.plugins_.input_plugins_mut()[0].read(inp, &mut self.ctx_) != 0 {
            return if self.ctx_.error_code != 0 { self.ctx_.error_code } else { SWMM_ERR_PARSE };
        }

        // Species registry: pollutants occupy the first slots.
        self.ctx_.species_registry.clear();
        for p in 0..self.ctx_.n_pollutants() {
            let name = self.ctx_.pollutant_names.name_of(p).to_vec();
            self.ctx_.species_registry.add(name, SpeciesKind::Pollutant, Vec::new());
        }

        // Resolve [PROCESS_COMPONENTS] registrations, then the embedded fallback, the ARD
        // transport rows and the per-element heat attributes.
        {
            let mut base_dir: Bytes = Vec::new();
            if !inp.is_empty() {
                base_dir = FsPath::new(inp).parent_path().native().to_vec();
            }

            let mut errs: Vec<Bytes> = Vec::new();
            if !self.ctx_.process_component_specs.is_empty() {
                errs = crate::core::process_components::resolve_process_components(&mut self.ctx_, &base_dir);
            }

            let mut reactions_registered = false;
            for spec in &self.ctx_.process_component_specs {
                if spec.id == b"org.hydrocouple.openswmm.reactions" {
                    reactions_registered = true;
                }
            }
            crate::transport::reactions_component::apply_embedded_reaction_sections(&mut self.ctx_, reactions_registered, &mut errs);

            crate::transport::ard_config::resolve_ard_transport_rows(&mut self.ctx_, &mut errs);

            {
                let pe = crate::transport::heat_overrides::resolve_heat_overrides(&mut self.ctx_);
                for e in pe {
                    errs.push(e);
                }
            }

            if !errs.is_empty() {
                for e in &errs {
                    self.ctx_.errors.push(e.clone());
                }
                if !self.lenient_open_ || self.metadata_only_ {
                    self.ctx_.error_code = SWMM_ERR_PARSE;
                    self.ctx_.error_message = errs[0].clone();
                    self.write_open_failure_report();
                    return SWMM_ERR_PARSE;
                }
            }

            crate::transport::reaction_legacy_binding::warn_if_legacy_binding_bypassed(&mut self.ctx_);
            crate::transport::ard_config::warn_if_fv_dispersion_key_ignored(&mut self.ctx_);

            if self.ctx_.heat_config.sediment_exchange && !self.ctx_.heat_config.sediment.has_ground_temp {
                self.ctx_.warnings.push(
                    "[SEDIMENT_EXCHANGE] GROUND_TEMPERATURE is not set; the default 12 degC is in force. For a buried conduit this term is often the largest one in the energy balance, so the default is unlikely to be what you meant."
                        .as_bytes()
                        .to_vec(),
                );
            }

            if self.ctx_.options.quality_solver != QualitySolverKind::Lagrangian && !self.ctx_.options.ignore_quality {
                if self.ctx_.options.quality_step > 0.0 {
                    self.ctx_.warnings.push(
                        "[OPTIONS] QUALITY_STEP is set but QUALITY_SOLVER is not LAGRANGIAN \u{2014} only the LARD engine substeps on this key; it has no effect this simulation."
                            .as_bytes()
                            .to_vec(),
                    );
                }
                if self.ctx_.options.max_segments_per_link != 100 {
                    self.ctx_.warnings.push(
                        "[OPTIONS] MAX_SEGMENTS_PER_LINK is set but QUALITY_SOLVER is not LAGRANGIAN \u{2014} only the LARD engine holds segments; it has no effect this simulation."
                            .as_bytes()
                            .to_vec(),
                    );
                }
                if self.ctx_.options.lard_rwpt {
                    self.ctx_.warnings.push(
                        "[OPTIONS] DISPERSION RWPT is set but QUALITY_SOLVER is not LAGRANGIAN \u{2014} RWPT is the LARD engine's dispersion; the ARD engine reads dispersion from the transport.ard component instead. No dispersion is applied from this key this simulation."
                            .as_bytes()
                            .to_vec(),
                    );
                }
            }
            // D-NS1 (X6): a configured negative quality inflow is EXTRACTION.
            for i in 0..self.ctx_.ext_inflows.count() {
                let ui = i as usize;
                if self.ctx_.ext_inflows.inflow_type[ui] == b"FLOW" {
                    continue;
                }
                if self.ctx_.ext_inflows.baseline[ui] >= 0.0 {
                    continue;
                }
                let m = cat(&[
                    b"[INFLOWS] node '",
                    &self.ctx_.ext_inflows.node_name[ui],
                    b"' pollutant '",
                    &self.ctx_.ext_inflows.constituent[ui],
                    "': the baseline is negative \u{2014} this row EXTRACTS mass (D-NS1). Extraction clamps per step to the mass the element holds.".as_bytes(),
                ]);
                self.ctx_.warnings.push(m);
            }

            // A2b: the REPORTED species list.
            self.ctx_.reported_species_names.clear();
            for p in 0..self.ctx_.n_pollutants() {
                let name = self.ctx_.pollutant_names.name_of(p).to_vec();
                self.ctx_.reported_species_names.push(name);
            }

            if self.ctx_.options.water_age {
                self.ctx_.species_registry.add(b"__WATER_AGE__".to_vec(), SpeciesKind::ReservedAge, b"hours".to_vec());
                self.ctx_.reported_species_names.push(b"__WATER_AGE__".to_vec());
                if self.ctx_.options.ignore_quality {
                    self.ctx_.warnings.push(
                        "[OPTIONS] WATER_AGE ON but IGNORE_QUALITY is YES \u{2014} the quality stage does not run, so no age is tracked this simulation."
                            .as_bytes()
                            .to_vec(),
                    );
                }
            }

            if self.ctx_.options.heat_transport {
                self.ctx_.species_registry.add(b"__TEMPERATURE__".to_vec(), SpeciesKind::ReservedTemperature, b"degC".to_vec());
                self.ctx_.reported_species_names.push(b"__TEMPERATURE__".to_vec());
                if self.ctx_.options.ignore_quality {
                    self.ctx_.warnings.push(
                        "[OPTIONS] HEAT_TRANSPORT ON but IGNORE_QUALITY is YES \u{2014} the quality stage does not run, so no temperature is tracked this simulation."
                            .as_bytes()
                            .to_vec(),
                    );
                }
            }

            // E1: reactions-component species report trailing the reserved rows.
            self.ctx_.refresh_reported_species();
        }

        // Warn about unknown/skipped sections (legacy warns once per pass, twice in all).
        if !self.ctx_.unknown_section_headers.is_empty() {
            let mut msgs: Vec<Bytes> = Vec::new();
            for (tok, line) in &self.ctx_.unknown_section_headers {
                msgs.push(cat(&[b"\n  WARNING: Unknown section '", tok, b"' at line ", line.to_string().as_bytes(), b" will be skipped."]));
            }
            let tail = std::mem::take(&mut self.ctx_.warnings);
            self.ctx_.warnings = msgs.clone();
            self.ctx_.warnings.extend(tail);
            for m in &msgs {
                self.emit_warning(100, m);
            }
            for m in &msgs {
                self.push_report_warning(m, 100);
            }
        } else {
            for tag in self.plugins_.input_plugins()[0].skipped_sections() {
                let m = cat(&[b"WARNING: Unknown section '[", &tag, b"]' will be skipped."]);
                self.push_report_warning(&m, 100);
            }
        }

        // Resolve cross-references (forward refs, final array sizing, head init)
        crate::input::post_parse_resolver::resolve_cross_references(&mut self.ctx_);
        // BW-MSX: bind the [BUILDUP]/[WASHOFF]/[LOADINGS] rows that name a species.
        for w in crate::transport::msx_surface::resolve(&mut self.ctx_) {
            self.ctx_.warnings.push(w);
        }

        // Project-level sanity checks (legacy project_validate), only without input errors.
        if self.ctx_.errors.is_empty() {
            self.validate_project();
        }

        // Post-parse validation errors are fatal (strict open).
        if !self.ctx_.errors.is_empty() {
            self.write_open_failure_report();
            if !self.lenient_open_ || self.metadata_only_ {
                let first = self.ctx_.errors[0].clone();
                self.set_error(SWMM_ERR_PARSE, &first);
                return SWMM_ERR_PARSE;
            }
            // Lenient open: keep the errors and fall through to OPENED.
        }

        // Metadata discovery stops before runtime plugin factories and output registration.
        if self.metadata_only_ {
            self.ctx_.state = EngineState::Opened;
            return SWMM_OK;
        }

        // Phase 4: load plugins listed in [PLUGINS]
        if !self.ctx_.plugin_specs.is_empty() {
            let mut warned: Vec<Bytes> = Vec::new();
            let specs = self.ctx_.plugin_specs.clone();
            self.plugins_.load_plugins(&specs, &mut |m: &[u8]| warned.push(m.to_vec()), Some(&self.ctx_));
            for m in &warned {
                self.emit_warning(SWMM_ERR_PLUGIN, m);
            }
        }

        // Explicit per-role writer choices are independent of the plugin catalogue.
        if !self.selected_writer(false) || !self.selected_writer(true) {
            self.set_error(SWMM_ERR_PLUGIN, b"The selected output/report writer or destination is unavailable");
            return SWMM_ERR_PLUGIN;
        }
        // Native companions (IO_OUTPUT_WRITER / IO_REPORT_WRITER naming a plugin): reached
        // only after selected_writer accepted a plugin, which needs a loaded external library
        // (never, here). Upstream then renames out_path_/rpt_path_ to the .inp's .out/.rpt.

        // Inject built-in output/report plugins when paths are configured.
        if !self.out_path_.is_empty() {
            let mut op = DefaultOutputPlugin::new(self.out_path_.clone());
            op.initialize(&[]);
            op.validate(&self.ctx_);
            self.plugins_.add_output_plugin(op);
            self.binary_output_plugin_ = Some(self.plugins_.output_plugins().len() - 1);
        }
        if !self.rpt_path_.is_empty() {
            let mut rp = DefaultReportPlugin::new(self.rpt_path_.clone());
            rp.initialize(&[]);
            rp.validate(&self.ctx_);
            self.plugins_.add_report_plugin(rp);
        }

        // ctx_.state_accessors: NOT YET TRANSLATED (see the file header).

        // Inject the default state-IO plugin if none resolved during [PLUGINS] processing.
        {
            let has_default_state_io = !self.plugins_.state_io_plugins().is_empty();
            if !has_default_state_io {
                let mut sp = DefaultStateIOPlugin::new();
                sp.initialize(&[]);
                sp.validate(&self.ctx_);
                self.plugins_.add_state_io_plugin(sp);
            }
        }

        self.ctx_.state = EngineState::Opened;
        SWMM_OK
    }

    /// open()'s selectedWriter lambda. An empty or BUILTIN choice keeps the built-in
    /// writer. Any other id must name a component that creates a writer: a built-in info's
    /// create_*_plugin() returns null, an unknown id is not found, and a library path is
    /// never loaded here, so every such choice is refused, as upstream refuses it when the
    /// library does not load. NOT TRANSLATED: a loadable external writer library.
    fn selected_writer(&mut self, report: bool) -> bool {
        let key: &[u8] = if report { b"IO_REPORT_WRITER" } else { b"IO_OUTPUT_WRITER" };
        let id = self.ctx_.options.ext_options.get(key).cloned().unwrap_or_default();
        if id.is_empty() || id == b"BUILTIN" {
            return true;
        }
        let info = self.plugins_.find_component(&id, None);
        let info = match info {
            Some(i) => i,
            None => return false,
        };
        if if report { !info.has_report } else { !info.has_output } {
            return false;
        }
        // The destination is computed (and a native deck guarded against overwrite) before
        // the writer is created; a built-in info creates no writer.
        false
    }

    // ========================================================================
    // Private helpers
    // ========================================================================

    pub fn validate_project(&mut self) {
        // ERROR 108: an outlet name matching both a node and a subcatchment.
        for s in 0..self.ctx_.n_subcatches() {
            let us = s as usize;
            if us >= self.ctx_.subcatches.outlet_name.len() {
                continue;
            }
            let name = &self.ctx_.subcatches.outlet_name[us];
            if !name.is_empty() && self.ctx_.node_names.find(name) >= 0 && self.ctx_.subcatch_names.find(name) >= 0 {
                let e = format_error(ERR_SUBCATCH_OUTLET, self.ctx_.subcatch_names.name_of(s));
                self.ctx_.errors.push(e);
            }
        }

        // Legacy gage_validate checks only a used TIMESERIES gage that is not a co-gage.
        let ng = self.ctx_.n_gages();
        for g in 0..ng {
            let ug = g as usize;
            if !crate::hydrology::gage::gage_is_used(&self.ctx_, g) {
                continue;
            }
            if self.ctx_.gages.source[ug] != RainSource::Timeseries {
                continue;
            }
            if self.ctx_.gages.co_gage_index[ug] >= 0 {
                continue;
            }

            // ERROR 159 / WARNING 09: recording interval against the series' smallest spacing.
            let ts = self.ctx_.gages.ts_index[ug];
            if ts >= 0 && ts < self.ctx_.tables.tables.len() as i32 {
                let tx = &self.ctx_.tables.tables[ts as usize].x;
                if tx.len() > 1 {
                    let mut dx_min = tx[1] - tx[0];
                    for k in 2..tx.len() {
                        dx_min = cmin(dx_min, tx[k] - tx[k - 1]);
                    }
                    let ts_interval = cvt_i32((dx_min * 86400.0 + 0.5).floor());
                    if ts_interval > 0 && self.ctx_.gages.interval_sec[ug] > ts_interval {
                        let e = format_error(ERR_RAIN_GAGE_INTERVAL, self.ctx_.gage_names.name_of(g));
                        self.ctx_.errors.push(e);
                    }
                    if self.ctx_.gages.interval_sec[ug] < ts_interval {
                        let w = format_warning(WARN_GAGE_INTERVAL, self.ctx_.gage_names.name_of(g));
                        self.ctx_.warnings.push(w);
                    }
                }
            }

            let interval = self.ctx_.gages.interval_sec[ug] as f64;
            if interval > 0.0 && self.ctx_.options.wet_step > interval {
                self.ctx_.options.wet_step = interval;
                let w = format_warning(WARN_WET_STEP_REDUCED, self.ctx_.gage_names.name_of(g));
                self.ctx_.warnings.push(w);
            }
        }

        // WARNING 06: dry-weather step increased to the wet-weather step.
        let opt = &mut self.ctx_.options;
        if opt.dry_step < opt.wet_step {
            opt.dry_step = opt.wet_step;
            self.ctx_.warnings.push(format_warning(WARN_DRY_STEP_INCREASED, b""));
        }

        // WARNING 07: routing step reduced to the wet-weather step; the authored step is kept.
        let opt = &mut self.ctx_.options;
        opt.routing_step_authored = opt.routing_step;
        if opt.routing_step > opt.wet_step {
            opt.routing_step = opt.wet_step;
            self.ctx_.warnings.push(format_warning(WARN_ROUTING_STEP_REDUCED, b""));
        }
    }

    /// Writes the .rpt of a failed open. NOT YET TRANSLATED (the .rpt writer): no file is
    /// written. The context is not changed either way (prepare/finalize take it const).
    fn write_open_failure_report(&mut self) {
        if self.rpt_path_.is_empty() || self.ctx_.options.rpt_disabled {
            return;
        }
    }

    /// `message` is a C string upstream (callers pass `.c_str()`), so it ends at a NUL.
    pub fn set_error(&mut self, code: i32, message: &[u8]) {
        let message = crate::core::cstdlib::c_str(message);
        self.ctx_.error_code = code;
        self.ctx_.error_message = message.to_vec();
        self.ctx_.state = EngineState::ErrorState;
        self.emit_warning(code, message);
    }

    pub fn emit_warning(&mut self, code: i32, message: &[u8]) {
        if let Some(cb) = self.callbacks_.on_warning.as_mut() {
            cb(code, crate::core::cstdlib::c_str(message));
        }
    }

    pub fn push_report_warning(&mut self, message: &[u8], code: i32) {
        self.ctx_.warnings.push(message.to_vec());
        self.emit_warning(code, message);
    }
}

/// Gate only (cargo feature `gate`): the members gate/diff's engine dump covers, writable, for
/// the context loader in examples/common/ctxload.rs. Not part of the translation.
#[cfg(feature = "gate")]
pub struct GateMembers<'a> {
    pub lenient_open_: &'a mut bool,
    pub metadata_only_: &'a mut bool,
    pub rpt_path_: &'a mut Bytes,
    pub out_path_: &'a mut Bytes,
    pub plugins_: &'a mut PluginFactory,
    pub binary_output_plugin_: &'a mut Option<usize>,
}

/// Gate only (cargo feature `gate`): the hydrology solvers initQuality and initMassBalance
/// read, writable, so examples/initqual.rs can install upstream's post-initialize state
/// into them (see gate/diff/initqual.cpp, massbal). Not part of the translation.
#[cfg(feature = "gate")]
pub struct GateHydrology<'a> {
    pub runoff_: &'a mut RunoffSolver,
    pub snow_: &'a mut SnowSolver,
    pub groundwater_: &'a mut GWSolver,
    pub lid_: &'a mut LIDSolver,
    /// the rest of what the hydrology dump covers (examples/common/hydroload.rs)
    pub rdii_: &'a mut RDIISolver,
    pub gage_states_: &'a mut Vec<GageState>,
    pub last_melt_doy_: &'a mut i32,
    /// gw_frac_perv_ ... old_gw_rate_, in the hydrology dump's order
    pub vectors: [&'a mut Vec<f64>; 12],
    pub gw_q_node_: &'a mut Vec<i32>,
    pub climate_file_: &'a mut ClimateFileReader,
}

/// Gate only (cargo feature `gate`): the solver members rtdump prints, writable, for the
/// routing-step cross test's loader (examples/common/ctxload.rs load_rt). Not part of the
/// translation.
#[cfg(feature = "gate")]
pub struct GateSolvers<'a> {
    pub router_: &'a mut crate::hydraulics::routing::Router,
    pub controls_: &'a mut crate::controls::controls::ControlEngine,
    pub inflow_: &'a mut crate::hydrology::inflow::InflowSolver,
    pub rdii_: &'a mut RDIISolver,
    pub exfil_: &'a mut crate::hydraulics::exfiltration::ExfilSolver,
    pub inlet_: &'a mut crate::hydraulics::inlet::InletSolver,
    pub hydstruct_: &'a mut crate::hydraulics::hyd_structures::StructureSolver,
    pub iface_: &'a mut crate::core::interface_file::InterfaceManager,
    pub next_event_: &'a mut i32,
    pub between_events_: &'a mut bool,
}

/// Gate only (cargo feature `gate`): the engine members the routing-step cross test's
/// stepdump prints (examples/routestep.rs). Not part of the translation.
#[cfg(feature = "gate")]
pub struct GateStepState<'a> {
    pub last_step_steady_: bool,
    pub new_rule_time_ms_: f64,
    pub checkpoint_enabled_: bool,
    pub ard_init_attempted_: bool,
    pub rdii_iface_open_: bool,
    pub rdii_iface_writing_: bool,
    pub wet_q_interp_: &'a [f64],
    pub lid_drain_q_units_: &'a [(i32, f64)],
    pub gw_q_interp_: &'a [f64],
    pub gw_q_node_: &'a [i32],
}
/// Writable twin of GateStepState.
#[cfg(feature = "gate")]
pub struct GateStepStateMut<'a> {
    pub last_step_steady_: &'a mut bool,
    pub new_rule_time_ms_: &'a mut f64,
    pub checkpoint_enabled_: &'a mut bool,
    pub ard_init_attempted_: &'a mut bool,
    pub wet_q_interp_: &'a mut Vec<f64>,
    pub lid_drain_q_units_: &'a mut Vec<(i32, f64)>,
    pub gw_q_interp_: &'a mut Vec<f64>,
    pub gw_q_node_: &'a mut Vec<i32>,
}

/// Gate only (cargo feature `gate`): the quality members qualdump prints, writable.
#[cfg(feature = "gate")]
pub struct GateQuality<'a> {
    pub hotstart_deferred_: &'a mut bool,
    pub quality_: &'a mut QualitySolver,
    pub landuse_solver_: &'a mut LanduseSolver,
    pub surface_quality_: &'a mut SurfaceQualitySoA,
}

#[cfg(feature = "gate")]
impl SwmmEngine {
    pub fn gate_solvers(&mut self) -> GateSolvers<'_> {
        GateSolvers {
            router_: &mut self.router_,
            controls_: &mut self.controls_,
            inflow_: &mut self.inflow_,
            rdii_: &mut self.rdii_,
            exfil_: &mut self.exfil_,
            inlet_: &mut self.inlet_,
            hydstruct_: &mut self.hydstruct_,
            iface_: &mut self.iface_,
            next_event_: &mut self.next_event_,
            between_events_: &mut self.between_events_,
        }
    }
    pub fn gate_step_state(&self) -> GateStepState<'_> {
        GateStepState {
            last_step_steady_: self.last_step_steady_,
            new_rule_time_ms_: self.new_rule_time_ms_,
            checkpoint_enabled_: self.checkpoint_enabled_,
            ard_init_attempted_: self.ard_init_attempted_,
            rdii_iface_open_: self.rdii_iface_file_.is_open(),
            rdii_iface_writing_: self.rdii_iface_file_.is_writing(),
            wet_q_interp_: &self.wet_q_interp_,
            lid_drain_q_units_: &self.lid_drain_q_units_,
            gw_q_interp_: &self.gw_q_interp_,
            gw_q_node_: &self.gw_q_node_,
        }
    }
    pub fn gate_step_state_mut(&mut self) -> GateStepStateMut<'_> {
        GateStepStateMut {
            last_step_steady_: &mut self.last_step_steady_,
            new_rule_time_ms_: &mut self.new_rule_time_ms_,
            checkpoint_enabled_: &mut self.checkpoint_enabled_,
            ard_init_attempted_: &mut self.ard_init_attempted_,
            wet_q_interp_: &mut self.wet_q_interp_,
            lid_drain_q_units_: &mut self.lid_drain_q_units_,
            gw_q_interp_: &mut self.gw_q_interp_,
            gw_q_node_: &mut self.gw_q_node_,
        }
    }
    pub fn gate_quality(&mut self) -> GateQuality<'_> {
        GateQuality {
            hotstart_deferred_: &mut self.hotstart_deferred_,
            quality_: &mut self.quality_,
            landuse_solver_: &mut self.landuse_solver_,
            surface_quality_: &mut self.surface_quality_,
        }
    }
    pub fn gate_hydrology(&mut self) -> GateHydrology<'_> {
        GateHydrology {
            runoff_: &mut self.runoff_,
            snow_: &mut self.snow_,
            groundwater_: &mut self.groundwater_,
            lid_: &mut self.lid_,
            rdii_: &mut self.rdii_,
            gage_states_: &mut self.gage_states_,
            last_melt_doy_: &mut self.last_melt_doy_,
            vectors: [
                &mut self.gw_frac_perv_,
                &mut self.gw_perv_evap_,
                &mut self.gw_lid_perv_area_,
                &mut self.gw_lid_infil_vol_,
                &mut self.gw_lid_perv_evap_vol_,
                &mut self.gw_2d_recharge_vol_,
                &mut self.gw_infil_with_2d_,
                &mut self.gw_infil_rate_,
                &mut self.snow_rain_,
                &mut self.snow_snow_,
                &mut self.gw_q_interp_,
                &mut self.old_gw_rate_,
            ],
            gw_q_node_: &mut self.gw_q_node_,
            climate_file_: &mut self.climate_file_,
        }
    }
    pub fn gate_members(&mut self) -> GateMembers<'_> {
        GateMembers {
            lenient_open_: &mut self.lenient_open_,
            metadata_only_: &mut self.metadata_only_,
            rpt_path_: &mut self.rpt_path_,
            out_path_: &mut self.out_path_,
            plugins_: &mut self.plugins_,
            binary_output_plugin_: &mut self.binary_output_plugin_,
        }
    }
}
