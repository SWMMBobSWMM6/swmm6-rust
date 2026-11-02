// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/quality/MsxLidReactions.hpp, MsxLidReactions.cpp and
// MsxLidUnits.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Multispecies reactions in LID storage and treatment layers (namespace msxlid), init part:
// initialize (LID storage-node layer species stores, with enabledNode) and initializeUnits
// (conventional LID unit species pools, with the file-local volumes()).
//
// Differences of form:
//   * `double v[layers]` / `double (&water)[layers]` are `[f64; LAYERS]`.
//   * `std::vector<bool> untracked` is Vec<bool>.
//
// Run-time functions, translated with the quality phase (not here):
//   MsxLidReactions.cpp: temperature (file-local), the scratch structs, heldMass,
//     mobileWater, prepareNodes, outlet, receive, reactNodes, reactPool,
//     reactedSubcatchment.
//   MsxLidUnits.cpp: routeUnits, publishUnits, drainValues.

use crate::core::cmath::cmax;
use crate::core::context::SimulationContext;
use crate::core::unit_conversion::{ucf, LANDAREA};
use crate::hydrology::lid::{completed_unit_storage, LIDGroupSoA, LIDSolver, LIDType};
use crate::transport::lid_layer_common::{build_offsets, donors_for, layer_volumes, K_ABSENT, K_NL};

const LAYERS: usize = 4;

/// enabledNode: a storage row with LID cells whose control opted into [REACTION_LIDS].
fn enabled_node(c: &SimulationContext, row: i32) -> bool {
    let s = &c.node_subtypes.storages;
    row >= 0
        && row < s.count()
        && !s.lid_state[row as usize].cells.is_empty()
        && c.reactions.lid_enabled(s.lid[row as usize].control)
}

/// msxlid::initialize: ensure the opted-in LID-node layer species stores exist.
pub fn initialize(c: &mut SimulationContext) {
    if c.reactions.lid_controls.is_empty() {
        return;
    }
    let ns = c.reactions.n_species();
    let count = c.node_subtypes.storages.count();
    c.reactions.lid.node_mass.resize(count as usize, Vec::new());
    c.reactions.lid.node_outlet.resize(count as usize, Vec::new());
    c.reactions.lid.node_mobile_old.resize(count as usize, 0.);
    for row in 0..count {
        if enabled_node(c, row) {
            let ur = row as usize;
            let s = &c.node_subtypes.storages;
            let state = &s.lid_state[ur];
            let node = s.node_idx[ur];
            let r = &mut c.reactions;
            if r.lid.node_mass[ur].len() == state.cells.len() * (ns as usize) {
                continue;
            }
            r.lid.node_mobile_old[ur] = c.nodes.volume[node as usize];
            r.lid.node_mass[ur].resize(state.cells.len() * (ns as usize), 0.0);
            for i in 0..state.cells.len() {
                let water = if state.richards && i < state.richards_water.len() {
                    state.richards_water[i]
                } else {
                    state.cells[i].theta * state.cells[i].geometric_volume
                };
                for sp in 0..ns {
                    let v = water * r.msx_node_conc[(node * ns + sp) as usize];
                    r.lid.node_mass[ur][i * (ns as usize) + sp as usize] = v;
                }
            }
        }
    }
}

/// volumes(g, u, water): conventional LID active-layer water volumes (ft3) of one unit,
/// surface/pavement/soil/storage.
fn volumes(g: &LIDGroupSoA, u: i32, water: &mut [f64; LAYERS]) {
    let mut v = [0.0f64; K_NL];
    layer_volumes(g, u as usize, &mut v);
    water.copy_from_slice(&v[..LAYERS]);
    for value in water.iter_mut() {
        *value *= g.area[u as usize];
    }
    let d = donors_for(g.type_, g.soil_thick[u as usize] > 0.0);
    let donors = [d.surface, d.pavement, d.soil, d.storage];
    for k in 0..LAYERS {
        if donors[k] == K_ABSENT {
            water[k] = 0.0;
        }
    }
    if g.type_ == LIDType::VegSwale {
        water[0] = completed_unit_storage(g, u);
    }
}

/// msxlid::initializeUnits: conventional reactive LID species pools and the full-coverage
/// ownership mask.
pub fn initialize_units(c: &mut SimulationContext, solver: &LIDSolver) {
    if c.reactions.lid_controls.is_empty() {
        return;
    }
    let ns = c.reactions.n_species();
    let nsc = c.n_subcatches();
    {
        let s = &mut c.reactions.lid;
        s.group_offsets = build_offsets(solver);
        let units = *s.group_offsets.last().unwrap();
        s.mass.clear();
        s.mass.resize((units as usize) * LAYERS * (ns as usize), 0.0);
        s.old_water.clear();
        s.old_water.resize((units as usize) * LAYERS, 0.0);
        s.surface_conc.clear();
        s.surface_conc.resize((units as usize) * (ns as usize), 0.0);
        s.drain_conc.clear();
        s.drain_conc.resize((units as usize) * (ns as usize), 0.0);
        s.surface_export.clear();
        s.surface_export.resize((nsc as usize) * (ns as usize), 0.0);
        s.surface_carry.clear();
        s.surface_carry.resize((nsc as usize) * (ns as usize), 0.0);
        s.subcatchments.clear();
        s.subcatchments.resize(nsc as usize, 0);
    }
    let mut area: Vec<f64> = vec![0.; nsc as usize];
    let mut untracked: Vec<bool> = vec![false; nsc as usize];
    for t in 0..solver.num_groups() {
        let g = solver.group(t);
        for u in 0..g.count {
            let uu = u as usize;
            let sc = g.subcatch_idx[uu];
            if sc < 0 || sc >= nsc {
                continue;
            }
            if !c.reactions.lid_enabled(g.control_idx[uu]) {
                untracked[sc as usize] = true;
                continue;
            }
            if c.reactions.species_is_wall.iter().any(|&value| value != 0) {
                c.errors.push(b"Reactive conventional LIDs require a BULK-only reaction system.".to_vec());
            }
            c.reactions.lid.subcatchments[sc as usize] = 1;
            area[sc as usize] += g.area[uu];
            if g.type_ == LIDType::RainBarrel && g.stor_covered[uu] != 0 {
                c.errors.push(
                    b"Reactive rain barrels must be uncovered; covered-barrel rain returns need a separate soil receipt.".to_vec(),
                );
            }
            let flat = c.reactions.lid.group_offsets[t as usize] + u;
            let mut v = [0.0f64; LAYERS];
            volumes(g, u, &mut v);
            for k in 0..LAYERS {
                let r = &mut c.reactions;
                r.lid.old_water[(flat as usize) * LAYERS + k] = v[k];
                for sp in 0..ns as usize {
                    let mass = r.init_global[sp] * v[k];
                    r.lid.mass[((flat as usize) * LAYERS + k) * (ns as usize) + sp] = mass;
                    r.surface.led_init_buildup[sp] += mass * 28.317 * r.surface.mcf[sp];
                }
            }
        }
    }
    for sc in 0..nsc as usize {
        if c.reactions.lid.subcatchments[sc] != 0 {
            let full = c.subcatches.area[sc] / ucf(LANDAREA, &c.options);
            if untracked[sc] || (area[sc] - full).abs() > 0.001 * cmax(1., full) {
                c.errors.push(
                    b"Reactive conventional LIDs require full subcatchment coverage and every unit to opt into [REACTION_LIDS]. Partial capture needs a separate donor-mass adapter.".to_vec(),
                );
            }
        }
    }
}
