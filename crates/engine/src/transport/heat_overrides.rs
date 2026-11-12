// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/transport/components/HeatFluxModules/
// HeatOverrides.hpp and HeatOverrides.cpp (HydroCouple/openswmm.engine, branch swmm6_rel,
// commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache License,
// Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Per-element heat attribute overrides (PE2): resolve the parsed rows into dense
// per-link / per-node tables (called from SWMMEngine::open after every component has
// applied), and the lookups that fall back to the global configuration.
//
// 1D build: the OPENSWMM_HAS_2D branch (CELL2D / MESH_TAG rows) is the #else branch,
// which refuses mesh overrides. With the `twod` feature (12 November 2026) the 2D branch is
// translated: the rows land on the router's mesh, reached through ctx.twod_io.mesh.

use crate::core::context::SimulationContext;
use crate::core::options::Bytes;
use crate::data::bed_zone::SedimentConfig;
use crate::data::heat::RadiativeConfig;
use crate::data::heat_override::{is_sediment_attr, HeatAttr, HeatElemKind, HeatElement, HeatScope};

pub fn patch_radiative(c: &mut RadiativeConfig, a: HeatAttr, v: f64) {
    match a {
        HeatAttr::Albedo => c.albedo = v,
        HeatAttr::ShadeFactor => c.shade_factor = v,
        HeatAttr::SkyView => c.sky_view = v,
        HeatAttr::EmissWater => c.emiss_water = v,
        HeatAttr::EmissLandcover => c.emiss_landcover = v,
        HeatAttr::LandcoverTemp => c.landcover_temp = v,
        _ => {} // sediment attrs never reach here (is_sediment_attr)
    }
}

pub fn patch_sediment(c: &mut SedimentConfig, a: HeatAttr, v: f64) {
    match a {
        HeatAttr::SedThermalDiffusivity => c.thermal_diffusivity = v,
        HeatAttr::SedSoluteDiffusivity => c.solute_diffusivity = v,
        HeatAttr::SedBedThickness => c.bed_thickness = v,
        HeatAttr::SedGroundDepth => c.ground_depth = v,
        HeatAttr::SedGroundTemp => {
            c.ground_temp = v;
            // A per-element ground temperature is by definition stated.
            c.has_ground_temp = true;
        }
        HeatAttr::SedHyporheicVelocity => c.hyporheic_velocity = v,
        HeatAttr::SedDensity => c.sed_density = v,
        HeatAttr::SedSpecificHeat => c.sed_specific_heat = v,
        HeatAttr::SedMeshPorosity => c.mesh_porosity = v,
        HeatAttr::SedGroundSource => c.mesh_ground = v != 0.0,
        _ => {}
    }
}

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut m = Vec::new();
    for p in parts {
        m.extend_from_slice(p);
    }
    m
}

pub fn resolve_heat_overrides(ctx: &mut SimulationContext) -> Vec<Bytes> {
    let mut errors: Vec<Bytes> = Vec::new();
    ctx.heat_config.overrides.resolved = true;
    if ctx.heat_config.overrides.rows.is_empty() {
        return errors; // nothing sized, nothing allocated
    }

    let nl = ctx.n_links();
    let nn = ctx.n_nodes();

    let mut want_rad_link = false;
    let mut want_rad_node = false;
    let mut want_sed = false;
    let mut want_rad_cell = false;
    let mut want_sed_cell = false;
    for r in &ctx.heat_config.overrides.rows {
        if r.scope == HeatScope::Global {
            continue; // already in the base
        }
        if r.scope == HeatScope::Cell2d || r.scope == HeatScope::MeshTag {
            if is_sediment_attr(r.attr) {
                want_sed_cell = true;
            } else {
                want_rad_cell = true;
            }
            continue;
        }
        if is_sediment_attr(r.attr) {
            want_sed = true;
        } else if r.scope == HeatScope::Node {
            want_rad_node = true;
        } else {
            // TAG rows can land on either kind, so a TAG row arms both.
            want_rad_link = true;
            if r.scope == HeatScope::Tag {
                want_rad_node = true;
            }
        }
    }

    let radiative = ctx.heat_config.radiative;
    let sediment = ctx.heat_config.sediment;
    let ov = &mut ctx.heat_config.overrides;
    if want_rad_link && nl > 0 {
        ov.rad_link.clear();
        ov.rad_link.resize(nl as usize, radiative);
    }
    if want_rad_node && nn > 0 {
        ov.rad_node.clear();
        ov.rad_node.resize(nn as usize, radiative);
    }
    if want_sed && nl > 0 {
        ov.sed_link.clear();
        ov.sed_link.resize(nl as usize, sediment);
    }

    #[cfg(not(feature = "twod"))]
    if want_rad_cell || want_sed_cell {
        errors.push(b"model.heat: mesh overrides require a 2D-enabled engine.".to_vec());
    }
    #[cfg(feature = "twod")]
    if want_rad_cell || want_sed_cell {
        let mesh = ctx.twod_io.mesh.clone();
        let mesh = mesh.as_ref().map(|m| m.borrow());
        let nc = mesh.as_ref().map_or(0, |m| m.n_cells());
        let ov = &mut ctx.heat_config.overrides;
        if want_rad_cell {
            ov.rad_cell.clear();
            ov.rad_cell.resize(nc as usize, radiative);
        }
        if want_sed_cell {
            ov.sed_cell.clear();
            ov.sed_cell.resize(nc as usize, sediment);
        }
        for scope in [HeatScope::MeshTag, HeatScope::Cell2d] {
            for ri in 0..ov.rows.len() {
                if ov.rows[ri].scope != scope {
                    continue;
                }
                let name = ov.rows[ri].name.clone();
                let attr = ov.rows[ri].attr;
                let value = ov.rows[ri].value;
                let mut cell = -1;
                if scope == HeatScope::Cell2d {
                    let (parsed, used) = crate::core::charconv::from_chars_i32(&name);
                    match parsed {
                        Ok(c) if used == name.len() && c >= 1 && c <= nc => cell = c - 1,
                        _ => {
                            let mut m = b"model.heat: CELL2D '".to_vec();
                            m.extend_from_slice(&name);
                            m.extend_from_slice(b"' is not a valid 1-based mesh cell.");
                            errors.push(m);
                            continue;
                        }
                    }
                }
                let mut hit = false;
                for i in 0..nc {
                    let skip = if scope == HeatScope::Cell2d {
                        i != cell
                    } else {
                        let m = mesh.as_ref().unwrap();
                        i as usize >= m.tri_tag.len() || m.tri_tag[i as usize] != name
                    };
                    if skip {
                        continue;
                    }
                    hit = true;
                    if is_sediment_attr(attr) {
                        patch_sediment(&mut ov.sed_cell[i as usize], attr, value);
                    } else {
                        patch_radiative(&mut ov.rad_cell[i as usize], attr, value);
                    }
                }
                if !hit {
                    let mut m = b"model.heat: MESH_TAG '".to_vec();
                    m.extend_from_slice(&name);
                    m.extend_from_slice(b"' matches no mesh region.");
                    errors.push(m);
                }
            }
        }
    }

    // Two passes, in specificity order (D-PE3): TAG first, then elements.
    let rows = ctx.heat_config.overrides.rows.clone();
    for pass in 0..2 {
        for r in &rows {
            if r.scope == HeatScope::Global {
                continue;
            }
            let is_elem = r.scope == HeatScope::Node || r.scope == HeatScope::Link;
            if pass == 0 && r.scope != HeatScope::Tag {
                continue;
            }
            if pass == 1 && !is_elem {
                continue;
            }

            if r.scope == HeatScope::Tag {
                let mut hit = false;
                for j in 0..nl {
                    let uj = j as usize;
                    if uj >= ctx.links.tags.len() || ctx.links.tags[uj] != r.name {
                        continue;
                    }
                    hit = true;
                    let ov = &mut ctx.heat_config.overrides;
                    if is_sediment_attr(r.attr) {
                        if uj < ov.sed_link.len() {
                            patch_sediment(&mut ov.sed_link[uj], r.attr, r.value);
                        }
                    } else if uj < ov.rad_link.len() {
                        patch_radiative(&mut ov.rad_link[uj], r.attr, r.value);
                    }
                }
                if !is_sediment_attr(r.attr) {
                    for i in 0..nn {
                        let ui = i as usize;
                        if ui >= ctx.nodes.tags.len() || ctx.nodes.tags[ui] != r.name {
                            continue;
                        }
                        hit = true;
                        let ov = &mut ctx.heat_config.overrides;
                        if ui < ov.rad_node.len() {
                            patch_radiative(&mut ov.rad_node[ui], r.attr, r.value);
                        }
                    }
                }
                if !hit {
                    errors.push(cat(&[
                        b"model.heat: TAG '",
                        &r.name,
                        b"' matches no link or node, so this override would silently do nothing. Check [TAGS].",
                    ]));
                }
                continue;
            }

            if r.scope == HeatScope::Link {
                let j = ctx.link_names.find(&r.name);
                if j < 0 {
                    errors.push(cat(&[b"model.heat: unknown link '", &r.name, b"' in a per-element heat override."]));
                    continue;
                }
                let uj = j as usize;
                let ov = &mut ctx.heat_config.overrides;
                if is_sediment_attr(r.attr) {
                    if uj < ov.sed_link.len() {
                        patch_sediment(&mut ov.sed_link[uj], r.attr, r.value);
                    }
                } else if uj < ov.rad_link.len() {
                    patch_radiative(&mut ov.rad_link[uj], r.attr, r.value);
                }
                continue;
            }

            // NODE scope: only radiative rows reach here.
            let i = ctx.node_names.find(&r.name);
            if i < 0 {
                errors.push(cat(&[b"model.heat: unknown node '", &r.name, b"' in a per-element heat override."]));
                continue;
            }
            let ui = i as usize;
            let ov = &mut ctx.heat_config.overrides;
            if ui < ov.rad_node.len() {
                patch_radiative(&mut ov.rad_node[ui], r.attr, r.value);
            }
        }
    }
    errors
}

pub fn radiative_for<'a>(ctx: &'a SimulationContext, e: &HeatElement) -> &'a RadiativeConfig {
    let ov = &ctx.heat_config.overrides;
    if e.index >= 0 {
        let ui = e.index as usize;
        if e.kind == HeatElemKind::Link && ui < ov.rad_link.len() {
            return &ov.rad_link[ui];
        }
        if e.kind == HeatElemKind::Node && ui < ov.rad_node.len() {
            return &ov.rad_node[ui];
        }
        if e.kind == HeatElemKind::Cell2d && ui < ov.rad_cell.len() {
            return &ov.rad_cell[ui];
        }
    }
    &ctx.heat_config.radiative
}

pub fn sediment_for<'a>(ctx: &'a SimulationContext, e: &HeatElement) -> &'a SedimentConfig {
    let ov = &ctx.heat_config.overrides;
    if e.kind == HeatElemKind::Cell2d && e.index >= 0 && (e.index as usize) < ov.sed_cell.len() {
        return &ov.sed_cell[e.index as usize];
    }
    if e.kind == HeatElemKind::Link && e.index >= 0 {
        let ui = e.index as usize;
        if ui < ov.sed_link.len() {
            return &ov.sed_link[ui];
        }
    }
    &ctx.heat_config.sediment
}
