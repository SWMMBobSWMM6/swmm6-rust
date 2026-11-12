// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/heat/MeshHeat.hpp and MeshHeat.cpp
// (initializeMeshHeat, configureGroundLidOwners) (HydroCouple/openswmm.engine, branch swmm6_rel,
// commit 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version
// 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// MeshHeatState and the two initialize-time entry points, as far as phase 2a reaches them.
// Heat on the mesh needs a surface temperature row, and SurfaceRouter2D::initialize refuses any
// transported surface row first (phase 2b), so initializeMeshHeat always takes its early
// return here (after clear()), and configureGroundLidOwners always finds an inactive owner.
// The early-return tests are upstream's; the rest of each function is phase 2b and returns the
// crate's stub. The context and aquifer_ground pointers are flags (set or not).

#![allow(clippy::all)]

use crate::core::context::SimulationContext;
use crate::core::engine_init_hydraulics::Terminated;
use crate::core::engine_step::stub;
use crate::core::options::Bytes;
use crate::twod::data::mesh_data::MeshData;
use crate::twod::data::solver_options_2d::SolverOptions2D;
use crate::twod::data::surface_state_data::SurfaceStateData;
use crate::twod::heat::ground_heat::{GroundHeatState, GroundLidThermalOwner};

#[derive(Debug, Clone)]
pub struct MeshHeatState {
    pub context: bool,
    pub enabled: bool,
    pub surface: bool,
    pub radiative: bool,
    pub bed: bool,
    pub surface_inherits: bool,
    pub radiative_inherits: bool,
    pub ground_column: bool,
    pub standalone_ground: GroundHeatState,
    pub aquifer_ground: bool,
    pub aquifer_owned: Vec<u8>,
    pub channel_footprint: Vec<f64>,
    pub ground_now: f64,
    pub bed_energy: Vec<f64>,
    pub initial_bed_energy: Vec<f64>,
    pub bed_capacity: Vec<f64>,
    pub reported_temperature: Vec<f64>,
    pub atmospheric_j: Vec<f64>,
    pub radiative_j: Vec<f64>,
    pub sensible_j: Vec<f64>,
    pub ground_j: Vec<f64>,
    pub water_bed_j: Vec<f64>,
    pub water_j: Vec<f64>,
    pub latent_j: Vec<f64>,
    pub evaporated_m3: Vec<f64>,
    pub failures: Vec<i32>,
    pub ground_column_seconds: Vec<f64>,
    pub hyporheic_eligible: Vec<u8>,
    pub bed_species: i32,
    pub bed_age_row: i32,
    pub bed_temp_row: i32,
    pub bed_pore_volume: Vec<f64>,
    pub bed_carried_mass: Vec<f64>,
    pub bed_initial_mass: Vec<f64>,
    pub bed_exchange_mass: Vec<f64>,
    pub bed_generated_age: Vec<f64>,
}

impl Default for MeshHeatState {
    fn default() -> Self {
        MeshHeatState {
            context: false,
            enabled: false,
            surface: false,
            radiative: false,
            bed: false,
            surface_inherits: true,
            radiative_inherits: true,
            ground_column: false,
            standalone_ground: GroundHeatState::default(),
            aquifer_ground: false,
            aquifer_owned: Vec::new(),
            channel_footprint: Vec::new(),
            ground_now: 0.0,
            bed_energy: Vec::new(),
            initial_bed_energy: Vec::new(),
            bed_capacity: Vec::new(),
            reported_temperature: Vec::new(),
            atmospheric_j: Vec::new(),
            radiative_j: Vec::new(),
            sensible_j: Vec::new(),
            ground_j: Vec::new(),
            water_bed_j: Vec::new(),
            water_j: Vec::new(),
            latent_j: Vec::new(),
            evaporated_m3: Vec::new(),
            failures: Vec::new(),
            ground_column_seconds: Vec::new(),
            hyporheic_eligible: Vec::new(),
            bed_species: 0,
            bed_age_row: -1,
            bed_temp_row: -1,
            bed_pore_volume: Vec::new(),
            bed_carried_mass: Vec::new(),
            bed_initial_mass: Vec::new(),
            bed_exchange_mass: Vec::new(),
            bed_generated_age: Vec::new(),
        }
    }
}

impl MeshHeatState {
    pub fn clear(&mut self) {
        *self = MeshHeatState::default();
    }
}

/// initializeMeshHeat(ctx, mesh, state, options, aquifer, coverage): the 2a part is the reset
/// and the early return; past it (coupled heat with a surface temperature row) is phase 2b.
pub fn initialize_mesh_heat(ctx: &SimulationContext, _mesh: &MeshData, state: &mut SurfaceStateData, options: &SolverOptions2D) -> Result<Bytes, Terminated> {
    state.heat.clear();
    if !options.coupled_heat || !ctx.options.heat_transport || state.transport.temp_row < 0 {
        return Ok(Bytes::new());
    }
    Err(stub("2D heat on the mesh (initializeMeshHeat with a surface temperature row)"))
}

/// configureGroundLidOwners(ctx, mesh, owner, beds): sizes the owner's LID footprint table and
/// returns for an inactive owner (the only owner phase 2a has).
pub fn configure_ground_lid_owners(_ctx: &SimulationContext, mesh: &MeshData, owner: &mut GroundHeatState, n_beds: usize) -> Result<Bytes, Terminated> {
    owner.lid_owners.clear();
    owner.lid_owners.resize(mesh.n_cells() as usize, Vec::<GroundLidThermalOwner>::new());
    if !owner.active() {
        return Ok(Bytes::new());
    }
    let _ = n_beds;
    Err(stub("2D ground LID heat owners (configureGroundLidOwners with an active ground owner)"))
}
