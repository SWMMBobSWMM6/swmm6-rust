// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/BedZoneData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The bed / hyporheic transient-storage zone (plan H6b): its configuration
// ([SEDIMENT_EXCHANGE]) with upstream's defaults value for value, and the runtime bed
// state per link and per ARD cell. resizeCells is resize_cells. clear() replaces the
// whole struct with its default, as `*this = BedZoneState{}` does.

/// `[SEDIMENT_EXCHANGE]` — the bed zone's material and geometry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SedimentConfig {
    /// Explicit routed-channel finite mesh contact; default prescribed ground.
    pub mesh_ground: bool,
    /// Preserve an explicitly authored PRESCRIBED choice.
    pub ground_source_configured: bool,
    /// Pore fraction of the new mesh bed; legacy 1D capacity is unchanged.
    pub mesh_porosity: f64,
    /// Bed thermal diffusivity alpha_sed, m²/s.
    pub thermal_diffusivity: f64,
    /// Bed effective solute diffusivity D_sed, m²/s (one value for every species).
    pub solute_diffusivity: f64,
    /// Bed layer thickness Y_hts, m.
    pub bed_thickness: f64,
    /// Depth from the bed to the deep-ground boundary Y_gr, m.
    pub ground_depth: f64,
    /// Deep-ground temperature T_gr, °C. Constant spelling.
    pub ground_temp: f64,
    /// Set by the parser when the deck stated `GROUND_TEMPERATURE` (warn if not).
    pub has_ground_temp: bool,
    /// `GROUND_TEMPERATURE TIMESERIES <name>` — index into ctx.tables, -1 for constant.
    pub ground_ts_index: i32,
    /// Hyporheic exchange velocity v_hyp, m/s. `Q_hts = v_hyp * A_bed`.
    pub hyporheic_velocity: f64,
    /// Bed bulk density rho_sed, kg/m³.
    pub sed_density: f64,
    /// Bed specific heat c_sed, J/kg/K.
    pub sed_specific_heat: f64,
    /// Initial bed temperature, °C.
    pub initial_temp: f64,
    pub has_initial_temp: bool,
}

impl Default for SedimentConfig {
    fn default() -> Self {
        SedimentConfig {
            mesh_ground: false,
            ground_source_configured: false,
            mesh_porosity: 0.30,
            thermal_diffusivity: 1.0e-6,
            solute_diffusivity: 1.0e-9,
            bed_thickness: 0.20,
            ground_depth: 2.0,
            ground_temp: 12.0,
            has_ground_temp: false,
            ground_ts_index: -1,
            hyporheic_velocity: 0.0,
            sed_density: 1670.0,
            sed_specific_heat: 1807.0,
            initial_temp: 12.0,
            has_initial_temp: false,
        }
    }
}

/// Runtime bed state, one entry per LINK (LEGACY/LARD) and one per CELL (EULERIAN_ARD).
#[derive(Debug, Clone, Default)]
pub struct BedZoneState {
    /// [link], °C
    pub link_temp: Vec<f64>,
    /// [link * n_species + s], species-minor
    pub link_conc: Vec<f64>,
    pub n_species: i32,
    pub seeded: bool,

    /// [cell], °C
    pub cell_temp: Vec<f64>,
    /// SPECIES-MAJOR [s * n_cells + c]
    pub cell_conc: Vec<f64>,
    pub cell_n_species: i32,
    pub cells_seeded: bool,
}

// The member initializers (n_species = cell_n_species = 0, seeded = cells_seeded =
// false) are Rust's defaults, so Default is derived.

impl BedZoneState {
    /// Replace cell-resolved bed state with initial temperature and zero solutes.
    pub fn resize_cells(&mut self, n_cells: i32, n_spec: i32, t_init: f64) {
        let nc = (if n_cells > 0 { n_cells } else { 0 }) as usize;
        let ns = (if n_spec > 0 { n_spec } else { 0 }) as usize;
        self.cell_n_species = ns as i32;
        self.cell_temp.clear();
        self.cell_temp.resize(nc, t_init);
        self.cell_conc.clear();
        self.cell_conc.resize(ns * nc, 0.0);
        self.cells_seeded = true;
    }

    /// Replace link bed state and request subsequent concentration seeding.
    pub fn resize(&mut self, n_links: i32, n_spec: i32, t_init: f64) {
        let nl = (if n_links > 0 { n_links } else { 0 }) as usize;
        let ns = (if n_spec > 0 { n_spec } else { 0 }) as usize;
        self.n_species = ns as i32;
        self.link_temp.clear();
        self.link_temp.resize(nl, t_init);
        self.link_conc.clear();
        self.link_conc.resize(nl * ns, 0.0);
        self.seeded = false;
    }

    /// Check the link-temperature array length only.
    pub fn sized(&self, n_links: i32) -> bool {
        self.link_temp.len() == (if n_links > 0 { n_links } else { 0 }) as usize
    }

    pub fn clear(&mut self) {
        *self = BedZoneState::default();
    }
}
