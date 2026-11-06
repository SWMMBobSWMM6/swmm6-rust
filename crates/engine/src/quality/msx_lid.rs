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
// The runoff-step part (6 November 2026): reactPool, reactedSubcatchment, routeUnits,
// publishUnits and drainValues. PoolScratchLease's thread_local PoolScratch is per-call
// locals (forms and spans are assigned, values and before written, and the workspace
// re-initialized before any read).
//
// NOT YET TRANSLATED (the routing step's node LIDs): MsxLidReactions.cpp's temperature
// (file-local), heldMass, mobileWater, prepareNodes, outlet, receive, reactNodes; except
// the read helpers translated with the quality stage.

use crate::core::cmath::cmax;
use crate::core::context::SimulationContext;
use crate::core::unit_conversion::{ucf, LANDAREA};
use crate::hydrology::lid::{completed_unit_storage, LIDGroupSoA, LIDSolver, LIDType};
use crate::core::engine_step_routing::StepError;
use crate::data::lid_layer_species::{LidLayer, LidSpecies};
use crate::data::reaction::ReactionExprForm;
use crate::data::reaction_tokens::{RxExprSpan, RxHydVar};
use crate::transport::lid_layer_common::{build_offsets, donors_for, layer_volumes, K_ABSENT, K_NL, K_STOR, K_SURF};
use crate::transport::reaction_legacy_step::{integrator_step, note_negative_clamp, Ws};

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

/// msxlid::drainValues (MsxLidUnits.cpp): the LID unit's drain concentrations for its species,
/// or None when the control is not a reactions LID or the unit has no row.
pub fn drain_values(c: &SimulationContext, type_: i32, unit: i32, control: i32) -> Option<&[f64]> {
    let r = &c.reactions;
    if !r.lid_enabled(control) || type_ < 0 || type_ + 1 >= r.lid.group_offsets.len() as i32 {
        return None;
    }
    let k = (r.lid.group_offsets[type_ as usize] + unit) as usize * r.n_species() as usize;
    if k + r.n_species() as usize <= r.lid.drain_conc.len() {
        Some(&r.lid.drain_conc[k..k + r.n_species() as usize])
    } else {
        None
    }
}

/// msxlid::reactPool: react one LID layer pool's species under its [REACTIONS] LID laws.
/// `reacted` and `clamped` are reassigned (zero) first, as upstream.
#[allow(clippy::too_many_arguments)]
pub fn react_pool(c: &mut SimulationContext, control: i32, layer: i32, volume: f64, mass: &mut [f64], dt: f64, temp: f64, reacted: &mut Vec<f64>, clamped: &mut Vec<f64>) -> Result<bool, StepError> {
    let ns = c.reactions.n_species();
    let uns = ns as usize;
    reacted.clear();
    reacted.resize(uns, 0.0);
    clamped.clear();
    clamped.resize(uns, 0.0);
    if !(volume > crate::quality::quality_routing::ZERO_VOLUME) || !(dt > 0.0) {
        return Ok(true);
    }
    let mut forms = vec![ReactionExprForm::None; uns];
    let mut spans = vec![RxExprSpan::default(); uns];
    let mut any = false;
    for select in [0, layer] {
        for law in c.reactions.lid_laws.iter() {
            if law.control == control && law.layer == select {
                forms[law.species as usize] = law.form;
                spans[law.species as usize] = law.compiled.clone();
                any = true;
            }
        }
    }
    if !any {
        return Ok(true);
    }
    let mut values = vec![0.0f64; uns];
    let mut before = vec![0.0f64; uns];
    for s in 0..uns {
        values[s] = mass[s] / volume;
        before[s] = values[s];
    }
    let mut variables = [0.0f64; RxHydVar::Count as usize];
    variables[RxHydVar::Temp as usize] = temp;
    let mut work = Ws::new(&c.reactions);
    let result = {
        let r = &c.reactions;
        integrator_step(r, &forms, &spans, dt, &mut values, &mut variables, &mut work, None, None, &r.coef_value)?
    };
    let mut valid = result.ok;
    for s in 0..uns {
        valid = valid && values[s].is_finite() && (values[s] >= 0.0 || forms[s] == ReactionExprForm::Rate);
    }
    if !valid {
        if !c.reactions.warned_react_failure {
            c.reactions.warned_react_failure = true;
            c.warnings.push(b"LID reaction solve failed; its unreacted inventory was retained.".to_vec());
        }
        return Ok(false);
    }
    for s in 0..uns {
        reacted[s] = (before[s] - values[s]) * volume;
        clamped[s] = cmax(0., -values[s]) * volume;
        mass[s] = cmax(0., values[s]) * volume;
        if clamped[s] > 0.0 {
            note_negative_clamp(c, b"a LID layer");
        }
    }
    Ok(true)
}

/// msxlid::reactedSubcatchment.
pub fn reacted_subcatchment(c: &SimulationContext, sc: i32) -> bool {
    sc >= 0 && sc < c.reactions.lid.subcatchments.len() as i32 && c.reactions.lid.subcatchments[sc as usize] != 0
}

/// msxlid::routeUnits: carry, react and drain the species of every reactions LID unit's
/// layers over one runoff step.
pub fn route_units(c: &mut SimulationContext, solver: &LIDSolver, dt: f64) -> Result<(), StepError> {
    if c.reactions.lid_controls.is_empty() || c.reactions.lid.group_offsets.is_empty() || !(dt > 0.0) {
        return Ok(());
    }
    let ns = c.reactions.n_species();
    let uns = ns as usize;
    c.reactions.lid.surface_export.fill(0.);
    let mut footprint = vec![0.0f64; c.n_subcatches().max(0) as usize];
    for t in 0..solver.num_groups() {
        let g = solver.group(t);
        for u in 0..g.count as usize {
            if c.reactions.lid_enabled(g.control_idx[u]) && g.subcatch_idx[u] >= 0 {
                footprint[g.subcatch_idx[u] as usize] += g.area[u];
            }
        }
    }
    for t in 0..solver.num_groups() {
        let g = solver.group(t);
        for uu in 0..g.count as usize {
            let u = uu as i32;
            if !c.reactions.lid_enabled(g.control_idx[uu]) {
                continue;
            }
            let sc = g.subcatch_idx[uu];
            let flat = c.reactions.lid.group_offsets[t as usize] + u;
            let control = g.control_idx[uu];
            if !reacted_subcatchment(c, sc) {
                continue;
            }
            let usc = sc as usize;
            let area = g.area[uu];
            let rain = cmax(0., c.subcatches.rainfall[usc]) * area * dt;
            let fraction = if footprint[usc] > 0.0 { area / footprint[usc] } else { 0. };
            let d = donors_for(g.type_, g.soil_thick[uu] > 0.0);
            let donor = [d.surface, d.pavement, d.soil, d.storage];
            let rates = [g.in_surf[uu], g.in_pave[uu], g.in_soil[uu], g.in_stor[uu]];
            let evap = [g.evap_surf[uu], g.evap_pave[uu], g.evap_soil[uu], g.evap_stor[uu]];
            let mut water = [0.0f64; LAYERS];
            volumes(g, u, &mut water);
            let mut incoming = [0.0f64; LAYERS];
            let mut internal = [0.0f64; LAYERS];
            let mut arrivals = vec![0.0f64; LAYERS * uns];
            let mut first: i32 = -1;
            let mut last: i32 = -1;
            for k in 0..LAYERS {
                if donor[k] != K_ABSENT {
                    if first < 0 {
                        first = k as i32;
                    }
                    last = k as i32;
                    incoming[k] = cmax(0., rates[k]) * area * dt;
                    if donor[k] >= 0 {
                        internal[donor[k] as usize] += incoming[k];
                    }
                }
            }
            if first < 0 {
                continue;
            }
            let ufirst = first as usize;
            incoming[ufirst] = cmax(0., g.inflow[uu]) * area * dt;
            for sp in 0..uns {
                let cr = if sp < c.reactions.rain_conc.len() { c.reactions.rain_conc[sp] } else { 0. };
                let index = usc * uns + sp;
                let runon = if index < c.reactions.surface.runon_qual_rate.len() { c.reactions.surface.runon_qual_rate[index] * dt * fraction / 28.317 } else { 0. };
                arrivals[ufirst * uns + sp] = cr * rain + runon;
                c.reactions.surface.led_wet_deposition[sp] += cr * rain * 28.317 * c.reactions.surface.mcf[sp];
            }
            let surface = cmax(0., g.surface_runoff[uu]) * area * dt;
            let drain = cmax(0., g.drain_flow[uu]) * area * dt;
            let infiltrated = cmax(0., g.infil_loss[uu]) * area;
            let surf_layer = if donor[K_SURF as usize] == K_ABSENT { first } else { K_SURF };
            let drain_layer = if donor[K_STOR as usize] == K_ABSENT { last } else { K_STOR };
            for k in 0..LAYERS {
                if donor[k] == K_ABSENT {
                    continue;
                }
                let ki = k as i32;
                let row = (flat as usize * LAYERS + k) * uns;
                let mut mass: Vec<f64> = c.reactions.lid.mass[row..row + uns].to_vec();
                for sp in 0..uns {
                    mass[sp] += arrivals[k * uns + sp];
                }
                let ow = flat as usize * LAYERS + k;
                let pool = c.reactions.lid.old_water[ow] + incoming[k];
                let mut reacted: Vec<f64> = Vec::new();
                let mut clamped: Vec<f64> = Vec::new();
                if cmax(c.reactions.lid.old_water[ow], water[k]) > 0.0353147 {
                    let mut temp = c.reactions.default_temp_c;
                    {
                        let thermal = &c.lid_layer_state;
                        if c.options.heat_transport && flat < thermal.n_units && thermal.n_species > LidSpecies::Temperature as i32 {
                            let layer = match k {
                                0 => LidLayer::Surface,
                                1 => LidLayer::Pavement,
                                2 => LidLayer::Soil,
                                _ => LidLayer::Storage,
                            };
                            let value = thermal.value[thermal.layer_index(flat, layer, LidSpecies::Temperature as i32)];
                            if value.is_finite() {
                                temp = value;
                            }
                        }
                    }
                    react_pool(c, control, ki + 1, pool, &mut mass, dt, temp, &mut reacted, &mut clamped)?;
                    for sp in 0..uns {
                        let factor = 28.317 * c.reactions.surface.mcf[sp];
                        c.reactions.surface.led_reacted[sp] += reacted[sp] * factor;
                        c.reactions.surface.led_reaction_clamped[sp] += clamped[sp] * factor;
                    }
                }
                let vs = if ki == surf_layer { surface } else { 0. };
                let vd = if ki == drain_layer { drain } else { 0. };
                let vi = if ki == last { infiltrated } else { 0. };
                let outgoing = internal[k] + vs + vd + vi;
                let mixed = cmax(cmax(pool - cmax(0., evap[k]) * area * dt, outgoing + water[k]), outgoing);
                for sp in 0..uns {
                    let concentration = if mixed > 0.0 { mass[sp] / mixed } else { 0. };
                    for dest in 0..LAYERS {
                        if donor[dest] == ki {
                            arrivals[dest * uns + sp] += concentration * incoming[dest];
                        }
                    }
                    c.reactions.lid.surface_export[usc * uns + sp] += concentration * vs;
                    if ki == surf_layer {
                        c.reactions.lid.surface_conc[flat as usize * uns + sp] = if surface > 0.0 { concentration } else { 0. };
                    }
                    if ki == drain_layer {
                        c.reactions.lid.drain_conc[flat as usize * uns + sp] = if drain > 0.0 { concentration } else { 0. };
                    }
                    c.reactions.surface.led_infil_loss[sp] += concentration * vi * 28.317 * c.reactions.surface.mcf[sp];
                    mass[sp] = cmax(0., mass[sp] - concentration * outgoing);
                }
                c.reactions.lid.mass[row..row + uns].copy_from_slice(&mass);
                c.reactions.lid.old_water[ow] = water[k];
            }
        }
    }
    Ok(())
}

/// msxlid::publishUnits: the reacted subcatchments' washoff concentration from the LID units'
/// surface export (carried while there is no runoff).
pub fn publish_units(c: &mut SimulationContext, dt: f64) {
    let ns = c.reactions.n_species();
    for sc in 0..c.n_subcatches() {
        if reacted_subcatchment(c, sc) {
            for sp in 0..ns {
                let k = sc as usize * ns as usize + sp as usize;
                let mass = c.reactions.lid.surface_export[k] + c.reactions.lid.surface_carry[k];
                let water = c.subcatches.runoff[sc as usize] * dt;
                c.reactions.surface.washoff_conc[k] = if water > 0.0 { mass / water } else { 0. };
                c.reactions.lid.surface_carry[k] = if water > 0.0 { 0. } else { mass };
            }
        }
    }
}

/// msxlid::prepareNodes, the part without a reactions-LID storage node: initialize the
/// node stores, arm the species routing ledger (`fv_retained`: the FV transport's retained
/// species mass, as route_legacy_msx takes it) and size msx_ext_mass_in. A deck with an
/// opted-in LID storage node row (enabledNode) is refused before anything changes: the node
/// half (the mobile/held split, the layer transfers, outlet, receive) is not translated.
pub fn prepare_nodes(c: &mut SimulationContext, dt: f64, fv_retained: &[f64]) -> Result<(), StepError> {
    if c.reactions.lid_controls.is_empty() || !(dt > 0.0) {
        return Ok(());
    }
    for row in 0..c.node_subtypes.storages.count() {
        if enabled_node(c, row) {
            return Err(StepError::NotTranslated("reaction species LID storage nodes (msxlid::prepareNodes)"));
        }
    }
    initialize(c);
    crate::transport::reaction_legacy_step::arm_msx_routing_ledger(c, fv_retained);
    let want = c.n_nodes() as usize * c.reactions.n_species() as usize;
    c.reactions.msx_ext_mass_in.resize(want, 0.0);
    Ok(())
}

/// msxlid::reactNodes: a no-op without an opted-in LID storage node row (prepare_nodes
/// refuses a deck that has one).
pub fn react_nodes(c: &mut SimulationContext, _dt: f64) -> Result<(), StepError> {
    if c.reactions.lid_controls.is_empty() {
        return Ok(());
    }
    for row in 0..c.node_subtypes.storages.count() {
        if enabled_node(c, row) {
            return Err(StepError::NotTranslated("reaction species LID storage nodes (msxlid::reactNodes)"));
        }
    }
    Ok(())
}

/// Whether any storage row is an opted-in reactions-LID node (enabledNode): the part of
/// prepareNodes / reactNodes not translated yet.
pub fn has_enabled_node(c: &SimulationContext) -> bool {
    (0..c.node_subtypes.storages.count()).any(|row| enabled_node(c, row))
}
