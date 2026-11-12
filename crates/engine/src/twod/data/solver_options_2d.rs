// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/data/SolverOptions2D.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The [2D_OPTIONS] store. Fields keep upstream's names, order and defaults; the
// `int8_t` tri-state fields are i8.

#![allow(clippy::all)]

use crate::core::options::{Bytes, FilePathPair};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i8)]
pub enum CellClosure2D {
    #[default]
    Flat = 0,
    Vfr = 1,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i8)]
pub enum FaceDepth2D {
    #[default]
    Mean = 0,
    VfrFace = 1,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i8)]
pub enum RainfallMode {
    #[default]
    NaturalNeighbour = 0,
    System = 1,
    None = 2,
    NearestNeighbour = 3,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i8)]
pub enum CouplingLaw2D {
    #[default]
    Orifice = 0,
    WeirOrifice = 1,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i8)]
pub enum CouplingHead2D {
    #[default]
    Cell = 0,
    Vertex = 1,
    StencilWetMean = 2,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i8)]
pub enum CouplingPondMode2D {
    #[default]
    Footprint = 0,
    Chamber = 1,
    Tracking = 2,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i8)]
pub enum Backend2D {
    Cpu = 0,
    #[default]
    Auto = 1,
    Omp = 2,
    Cuda = 3,
    Hip = 4,
    Sycl = 5,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i8)]
pub enum Momentum2D {
    #[default]
    LocalInertial = 0,
    FullSwe = 1,
    DiffusiveWave = 2,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i8)]
pub enum OutputPrecision2D {
    #[default]
    Float32 = 0,
    Float64 = 1,
}

/// `report2d::Var`: the REPORT_2D_VARIABLES bit groups.
pub mod report2d {
    pub const DEPTH: u32 = 1 << 0;
    pub const VELOCITY: u32 = 1 << 1;
    pub const EDGE_FLUX: u32 = 1 << 2;
    pub const NODE_HEAD: u32 = 1 << 3;
    pub const SPECIES: u32 = 1 << 4;
    pub const RAINFALL: u32 = 1 << 5;
    pub const INFILTRATION: u32 = 1 << 6;
    pub const COUPLING: u32 = 1 << 7;
    pub const GRADIENTS: u32 = 1 << 8;
    pub const CONTINUITY: u32 = 1 << 9;
    pub const ENVELOPES: u32 = 1 << 10;
    pub const BUILDUP: u32 = 1 << 11;
    pub const GROUNDWATER: u32 = 1 << 12;
    pub const GW_DETAILED: u32 = 1 << 13;
    pub const HEAT: u32 = 1 << 14;
    pub const TURBULENCE: u32 = 1 << 15;
    pub const ALL_MASK: u32 = (1u32 << 16) - 1u32;
    pub const DEFAULT_MASK: u32 =
        DEPTH | VELOCITY | EDGE_FLUX | NODE_HEAD | SPECIES | RAINFALL | INFILTRATION | ENVELOPES | BUILDUP | GROUNDWATER;
    pub const MINIMAL_MASK: u32 = DEPTH | NODE_HEAD | ENVELOPES;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i8)]
pub enum Coriolis2D {
    #[default]
    Off = 0,
    Latitude = 1,
    ConstantF = 2,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i8)]
pub enum Turbulence2D {
    #[default]
    None = 0,
    Smagorinsky = 1,
    MixingLength = 2,
}

/// turbulenceName()
pub fn turbulence_name(t: Turbulence2D) -> &'static [u8] {
    match t {
        Turbulence2D::Smagorinsky => b"SMAGORINSKY",
        Turbulence2D::MixingLength => b"MIXING_LENGTH",
        Turbulence2D::None => b"NONE",
    }
}

#[derive(Debug, Clone)]
pub struct SolverOptions2D {
    pub coupled_heat: bool,
    pub heat_model_authored: bool,
    pub heat_surface: i32,
    pub heat_radiative: i32,
    pub heat_surface_authored: bool,
    pub heat_radiative_authored: bool,
    pub bed_heat_model: i32,
    pub bed_heat_authored: bool,
    pub water_heat_capacity: f64,
    pub momentum: Momentum2D,
    pub coriolis: Coriolis2D,
    pub turbulence: Turbulence2D,
    pub coriolis_latitude: f64,
    pub coriolis_f: f64,
    pub smagorinsky_cs: f64,
    pub eddy_viscosity_background: f64,
    pub mixing_length_coefficient: f64,
    pub coriolis_latitude_set: bool,
    pub coriolis_f_set: bool,
    pub wall_drag_enabled: bool,
    pub front_rebuild: i32,
    pub reconstruction_order: i32,
    pub max_timestep: f64,
    pub dry_depth: f64,
    pub limiter_epsilon: f64,
    pub flux_dh_eps: f64,
    pub coupling_sync: f64,
    pub coupling_cd: f64,
    pub coupling_law: CouplingLaw2D,
    pub coupling_head: CouplingHead2D,
    pub coupling_coeffs: i32,
    pub coupling_pond_mode: CouplingPondMode2D,
    pub coupling_pond_relax: f64,
    pub coupling_rim_band: f64,
    pub coupling_equalize: i8,
    pub report_2d: bool,
    pub output_precision: OutputPrecision2D,
    pub output_compression: i32,
    pub report_2d_vars: u32,
    pub report_2d_species: Vec<Bytes>,
    pub report_2d_step: f64,
    pub infiltration: i8,
    pub infil_step: f64,
    pub infil_default_method: Bytes,
    pub infil_destination: Bytes,
    pub evaporation: i8,
    pub transport_pollutants: bool,
    pub transport_msx: bool,
    pub transport_age: bool,
    pub transport_temperature: bool,
    pub groundwater: i8,
    pub gw_et: Bytes,
    pub rainfall_mode: RainfallMode,
    pub cell_closure: CellClosure2D,
    pub face_reconstruction: FaceDepth2D,
    pub vfr_min_wet_frac: f64,
    pub theta: f64,
    pub theta_reference_step: f64,
    pub cfl_number: f64,
    pub h_move: f64,
    pub lts_tiers: i32,
    pub froude_max: f64,
    pub advection: bool,
    pub exchange_beta: f64,
    pub dispersion: f64,
    pub exchange_relax: f64,
    pub coupling_area_auto: bool,
    pub coupling_in_flooding: bool,
    pub backend: Backend2D,
    pub auto_backend_cpu: bool,
    pub mesh_file: FilePathPair,
    pub output_file: FilePathPair,
    pub len_1d_to_2d: f64,
    pub len_2d_to_1d: f64,
    pub vol_1d_to_2d: f64,
    pub flow_1d_to_2d: f64,
    pub flow_2d_to_1d: f64,
    pub num_threads: i32,
    pub requested_threads: i32,
    pub mesh_units_si: bool,
    pub mesh_scaled_to_si: bool,
    pub mesh_to_si_factor: f64,
    pub pending_rows_drained: bool,
    pub bc_flow_to_si_applied: f64,
}

impl Default for SolverOptions2D {
    fn default() -> Self {
        SolverOptions2D {
            coupled_heat: false,
            heat_model_authored: false,
            heat_surface: -1,
            heat_radiative: -1,
            heat_surface_authored: false,
            heat_radiative_authored: false,
            bed_heat_model: 0,
            bed_heat_authored: false,
            water_heat_capacity: 4186000.0,
            momentum: Momentum2D::LocalInertial,
            coriolis: Coriolis2D::Off,
            turbulence: Turbulence2D::None,
            coriolis_latitude: 0.0,
            coriolis_f: 0.0,
            smagorinsky_cs: 0.10,
            eddy_viscosity_background: 0.0,
            mixing_length_coefficient: 0.5,
            coriolis_latitude_set: false,
            coriolis_f_set: false,
            wall_drag_enabled: false,
            front_rebuild: -1,
            reconstruction_order: 1,
            max_timestep: 10.0,
            dry_depth: 0.001,
            limiter_epsilon: 1.0e-6,
            flux_dh_eps: 0.004,
            coupling_sync: 0.0,
            coupling_cd: 0.65,
            coupling_law: CouplingLaw2D::Orifice,
            coupling_head: CouplingHead2D::Cell,
            coupling_coeffs: 0,
            coupling_pond_mode: CouplingPondMode2D::Footprint,
            coupling_pond_relax: 0.0,
            coupling_rim_band: -1.0,
            coupling_equalize: -1,
            report_2d: true,
            output_precision: OutputPrecision2D::Float32,
            output_compression: 4,
            report_2d_vars: report2d::DEFAULT_MASK,
            report_2d_species: Vec::new(),
            report_2d_step: 0.0,
            infiltration: -1,
            infil_step: 0.0,
            infil_default_method: Bytes::new(),
            infil_destination: Bytes::new(),
            evaporation: 1,
            transport_pollutants: true,
            transport_msx: true,
            transport_age: true,
            transport_temperature: true,
            groundwater: -1,
            gw_et: Bytes::new(),
            rainfall_mode: RainfallMode::NaturalNeighbour,
            cell_closure: CellClosure2D::Flat,
            face_reconstruction: FaceDepth2D::Mean,
            vfr_min_wet_frac: 0.01,
            theta: 0.8,
            theta_reference_step: 0.0,
            cfl_number: 0.7,
            h_move: 0.003,
            lts_tiers: 4,
            froude_max: 1.5,
            advection: false,
            exchange_beta: 0.8,
            dispersion: 0.0,
            exchange_relax: 1.0,
            coupling_area_auto: false,
            coupling_in_flooding: false,
            backend: Backend2D::Auto,
            auto_backend_cpu: false,
            mesh_file: FilePathPair::default(),
            output_file: FilePathPair::default(),
            len_1d_to_2d: 1.0,
            len_2d_to_1d: 1.0,
            vol_1d_to_2d: 1.0,
            flow_1d_to_2d: 1.0,
            flow_2d_to_1d: 1.0,
            num_threads: 1,
            requested_threads: 0,
            mesh_units_si: false,
            mesh_scaled_to_si: false,
            mesh_to_si_factor: 1.0,
            pending_rows_drained: false,
            bc_flow_to_si_applied: 1.0,
        }
    }
}

impl SolverOptions2D {
    /// momentumPhysics()
    pub fn momentum_physics(&self) -> bool {
        self.wall_drag_enabled || self.coriolis != Coriolis2D::Off || self.turbulence != Turbulence2D::None
    }

    /// coriolisParameter()
    pub fn coriolis_parameter(&self) -> f64 {
        if self.coriolis == Coriolis2D::ConstantF {
            self.coriolis_f
        } else if self.coriolis == Coriolis2D::Latitude {
            1.4584231711061317e-4 * (self.coriolis_latitude * 0.017453292519943295).sin() // LIBM
        } else {
            0.0
        }
    }
}
