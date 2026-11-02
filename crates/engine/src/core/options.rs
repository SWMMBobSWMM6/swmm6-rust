// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/SimulationOptions.hpp,
// src/engine/hydraulics/fv/FvOptions.hpp and src/engine/core/FilePathPair.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Defaults are upstream's member initializers, value for value. Strings are byte strings.

use std::collections::HashMap;

pub type Bytes = Vec<u8>;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i32)]
pub enum FlowUnits {
    #[default]
    Cfs = 0,
    Gpm = 1,
    Mgd = 2,
    Cms = 3,
    Lps = 4,
    Mld = 5,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i32)]
pub enum RoutingModel {
    Steady = 0,
    Kinwave = 1,
    #[default]
    Dynwave = 2,
    Fv = 3,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i32)]
pub enum QualitySolverKind {
    #[default]
    Legacy = 0,
    EulerianArd = 1,
    Lagrangian = 2,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i32)]
pub enum InfiltrationModel {
    #[default]
    Horton = 0,
    ModHorton = 1,
    GreenAmpt = 2,
    ModGreenAmpt = 3,
    CurveNumber = 4,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i32)]
pub enum RunoffModel {
    None = 0,
    #[default]
    NlPond = 1,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i32)]
pub enum NodeContinuity {
    #[default]
    Explicit = 0,
    SemiImplicit = 1,
}

/// An authored file path and, once resolved, its absolute form. Equality compares the
/// authored form only, as upstream's operator== does.
#[derive(Debug, Clone, Default)]
pub struct FilePathPair {
    pub absolute: Bytes,
    pub original: Bytes,
}

impl FilePathPair {
    pub fn set(&mut self, s: &[u8]) {
        self.original = s.to_vec();
        self.absolute.clear();
    }
    pub fn is_empty(&self) -> bool {
        self.original.is_empty()
    }
}

impl PartialEq for FilePathPair {
    fn eq(&self, o: &Self) -> bool {
        self.original == o.original
    }
}

pub mod fv {
    pub const K_MAX_LTS_TIERS: i32 = 8;

    #[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
    #[repr(i32)]
    pub enum RiemannSolver {
        Hll = 0,
        #[default]
        Hllc = 1,
    }
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
    #[repr(i32)]
    pub enum Limiter {
        #[default]
        Minmod = 0,
        Vanleer = 1,
        Superbee = 2,
    }
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
    #[repr(i32)]
    pub enum ScalarScheme {
        Upwind = 0,
        #[default]
        Muscl = 1,
        QuickestUltimate = 2,
    }
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
    #[repr(i32)]
    pub enum TimeIntegration {
        #[default]
        Euler = 0,
        Rk2 = 1,
    }
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
    #[repr(i32)]
    pub enum StructureCoupling {
        #[default]
        Substep = 0,
        RoutingStep = 1,
    }
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
    #[repr(i32)]
    pub enum Backend {
        Cpu = 0,
        #[default]
        Auto = 1,
        Omp = 2,
        Cuda = 3,
        Hip = 4,
        Sycl = 5,
    }

    #[derive(Debug, Clone)]
    pub struct FvOptions {
        pub cell_length: f64,
        pub min_cells: i32,
        pub pressure_closure: i32,
        pub unsteady_friction: i32,
        pub uf_k3: f64,
        pub cfl: f64,
        pub riemann: RiemannSolver,
        pub order: i32,
        pub limiter: Limiter,
        pub scalar_scheme: ScalarScheme,
        pub time_integration: TimeIntegration,
        pub slot_celerity: f64,
        pub pressurized_implicit: bool,
        pub structure_coupling: StructureCoupling,
        pub backend: Backend,
        pub min_parallel_cells: i64,
        pub compaction: bool,
        pub lts: bool,
        pub node_feedback_dt: bool,
        pub lts_max_tiers: i32,
        pub cfl_census_interval: i32,
        pub dispersion: f64,
    }

    impl Default for FvOptions {
        fn default() -> Self {
            FvOptions {
                cell_length: 0.0,
                min_cells: 4,
                pressure_closure: 0,
                unsteady_friction: 0,
                uf_k3: 0.015,
                cfl: 0.5,
                riemann: RiemannSolver::Hllc,
                order: 1,
                limiter: Limiter::Minmod,
                scalar_scheme: ScalarScheme::Muscl,
                time_integration: TimeIntegration::Euler,
                slot_celerity: 100.0,
                pressurized_implicit: false,
                structure_coupling: StructureCoupling::Substep,
                backend: Backend::Auto,
                min_parallel_cells: 20000,
                compaction: true,
                lts: true,
                node_feedback_dt: false,
                lts_max_tiers: 6,
                cfl_census_interval: 1,
                dispersion: 0.0,
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct SimulationOptions {
    // ---- dates and steps (dates are OLE days since 12/30/1899; steps seconds)
    pub start_date: f64,
    pub end_date: f64,
    pub total_duration_ms: f64,
    pub report_start: f64,
    pub routing_step: f64,
    pub routing_step_authored: f64,
    pub min_routing_step: f64,
    pub dry_step: f64,
    pub wet_step: f64,
    pub report_step: f64,
    pub dry_days: f64,
    pub variable_step: f64,
    pub lengthening_step: f64,
    pub sweep_start: i32,
    pub sweep_end: i32,
    // ---- models
    pub flow_units: FlowUnits,
    pub routing_model: RoutingModel,
    pub quality_solver: QualitySolverKind,
    pub quality_step: f64,
    pub max_segments_per_link: i32,
    pub lard_rwpt: bool,
    pub rwpt_seed: i32,
    pub water_age: bool,
    pub outfall_backflow_zero: bool,
    pub heat_transport: bool,
    pub humidity: [f64; 12],
    pub humidity_type: i32,
    pub humidity_var: i32,
    pub humidity_ts_name: Bytes,
    pub water_density: f64,
    pub water_specific_heat: f64,
    pub wind_func_coeff_a: f64,
    pub wind_func_coeff_b: f64,
    pub pressure_ratio: f64,
    pub fv: fv::FvOptions,
    pub infiltration: InfiltrationModel,
    pub runoff_model: RunoffModel,
    // ---- dynamic wave
    pub max_trials: i32,
    pub surcharge_method: i32,
    pub tpa_celerity: f64,
    pub dps_target_celerity: f64,
    pub dps_alpha: f64,
    pub dps_decay_time: f64,
    pub unsteady_friction: i32,
    pub report_signed_heads: i32,
    pub uf_k3: f64,
    pub node_continuity: NodeContinuity,
    pub virtual_junction_momentum: i32,
    pub anderson_accel: bool,
    pub normal_flow_ltd: i32,
    pub force_main_eqn: i32,
    pub inertial_damping: i32,
    pub link_offsets: i32,
    pub min_slope: f64,
    pub min_surf_area: f64,
    pub head_tol: f64,
    pub sys_flow_tol: f64,
    pub lat_flow_tol: f64,
    pub skip_steady_state: bool,
    pub write_absolute_paths: bool,
    // ---- process switches
    pub quality_routing: bool,
    pub allow_ponding: bool,
    pub ignore_rainfall: bool,
    pub ignore_snow_melt: bool,
    pub ignore_groundwater: bool,
    pub ignore_rdii: bool,
    pub ignore_routing: bool,
    pub ignore_quality: bool,
    pub ignore_2d: bool,
    pub rule_step: f64,
    pub num_threads: i32,
    pub crs: Bytes,
    pub ext_options: HashMap<Bytes, Bytes>,
    // ---- evaporation, temperature, wind, snow, ADC
    pub evap_type: i32,
    pub evap_values: [f64; 12],
    pub evap_ts_name: Bytes,
    pub evap_recovery_pat: Bytes,
    pub evap_dry_only: bool,
    pub pan_coeff: [f64; 12],
    pub temp_source: i32,
    pub temp_ts_name: Bytes,
    pub temp_file: FilePathPair,
    pub temp_file_start: f64,
    pub temp_units: i32,
    pub wind_type: i32,
    pub wind_speed: [f64; 12],
    pub snow_divt: f64,
    pub snow_ati_wt: f64,
    pub snow_nrg_ratio: f64,
    pub snow_lat: f64,
    pub snow_dtlong: f64,
    pub snow_elev: f64,
    pub snow_min_melt: f64,
    pub snow_max_melt: f64,
    pub adc_imperv: [f64; 10],
    pub adc_perv: [f64; 10],
    // ---- report
    pub rpt_disabled: bool,
    pub rpt_subcatchments: i32,
    pub rpt_nodes: i32,
    pub rpt_links: i32,
    pub rpt_input: bool,
    pub rpt_continuity: bool,
    pub rpt_flowstats: bool,
    pub rpt_controls: bool,
    pub rpt_averages: bool,
    pub rpt_link_steps: bool,
    pub rpt_subcatch_names: Vec<Bytes>,
    pub rpt_node_names: Vec<Bytes>,
    pub rpt_link_names: Vec<Bytes>,
}

impl Default for SimulationOptions {
    fn default() -> Self {
        SimulationOptions {
            start_date: 37987.0, // 1 Jan 2004 (the legacy default)
            end_date: 0.0,
            total_duration_ms: -1.0,
            report_start: 0.0,
            routing_step: 20.0,
            routing_step_authored: 0.0,
            min_routing_step: 0.5,
            dry_step: 3600.0,
            wet_step: 300.0,
            report_step: 900.0,
            dry_days: 0.0,
            variable_step: 0.75,
            lengthening_step: 0.0,
            sweep_start: 1,
            sweep_end: 365,
            flow_units: FlowUnits::Cfs,
            routing_model: RoutingModel::Dynwave,
            quality_solver: QualitySolverKind::Legacy,
            quality_step: 0.0,
            max_segments_per_link: 100,
            lard_rwpt: false,
            rwpt_seed: 0,
            water_age: false,
            outfall_backflow_zero: false,
            heat_transport: false,
            humidity: [50.0; 12],
            humidity_type: 0,
            humidity_var: 0,
            humidity_ts_name: Vec::new(),
            water_density: 1000.0,
            water_specific_heat: 4184.0,
            wind_func_coeff_a: 1.505e-8,
            wind_func_coeff_b: 1.6e-8,
            pressure_ratio: 1.0,
            fv: fv::FvOptions::default(),
            infiltration: InfiltrationModel::Horton,
            runoff_model: RunoffModel::NlPond,
            max_trials: 8,
            surcharge_method: 0,
            tpa_celerity: 100.0,
            dps_target_celerity: 25.0,
            dps_alpha: 3.0,
            dps_decay_time: 0.5,
            unsteady_friction: 0,
            report_signed_heads: 0,
            uf_k3: 0.015,
            node_continuity: NodeContinuity::Explicit,
            virtual_junction_momentum: 0,
            anderson_accel: false,
            normal_flow_ltd: 2,  // BOTH (legacy default)
            force_main_eqn: 0,
            inertial_damping: 1, // PARTIAL (legacy default)
            link_offsets: 0,     // DEPTH_OFFSET (legacy default)
            min_slope: 0.0,
            min_surf_area: 0.0,
            head_tol: 0.0,
            sys_flow_tol: 0.05,
            lat_flow_tol: 0.05,
            skip_steady_state: false,
            write_absolute_paths: false,
            quality_routing: false,
            allow_ponding: false,
            ignore_rainfall: false,
            ignore_snow_melt: false,
            ignore_groundwater: false,
            ignore_rdii: false,
            ignore_routing: false,
            ignore_quality: false,
            ignore_2d: false,
            rule_step: 0.0,
            num_threads: 1,
            crs: Vec::new(),
            ext_options: HashMap::new(),
            evap_type: 0,
            evap_values: [0.0; 12],
            evap_ts_name: Vec::new(),
            evap_recovery_pat: Vec::new(),
            evap_dry_only: false,
            pan_coeff: [1.0; 12],
            temp_source: 0,
            temp_ts_name: Vec::new(),
            temp_file: FilePathPair::default(),
            temp_file_start: 0.0,
            temp_units: -1,
            wind_type: 0,
            wind_speed: [0.0; 12],
            snow_divt: 34.0,
            snow_ati_wt: 0.5,
            snow_nrg_ratio: 0.6,
            snow_lat: 0.0,
            snow_dtlong: 0.0,
            snow_elev: 0.0,
            snow_min_melt: 0.0,
            snow_max_melt: 0.0,
            adc_imperv: [1.0; 10],
            adc_perv: [1.0; 10],
            rpt_disabled: false,
            rpt_subcatchments: 1,
            rpt_nodes: 1,
            rpt_links: 1,
            rpt_input: false,
            rpt_continuity: true,
            rpt_flowstats: true,
            rpt_controls: false,
            rpt_averages: false,
            rpt_link_steps: false,
            rpt_subcatch_names: Vec::new(),
            rpt_node_names: Vec::new(),
            rpt_link_names: Vec::new(),
        }
    }
}

impl SimulationOptions {
    /// The routing step a link-length check uses: the authored one if there was one.
    pub fn link_validate_routing_step(&self) -> f64 {
        if self.routing_step_authored > 0.0 {
            self.routing_step_authored
        } else {
            self.routing_step
        }
    }

    pub fn total_duration_ms(&self) -> f64 {
        if self.total_duration_ms >= 0.0 {
            return self.total_duration_ms;
        }
        ((self.end_date - self.start_date) * 86400.0).floor() * 1000.0
    }
}
