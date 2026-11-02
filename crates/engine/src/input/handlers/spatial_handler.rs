// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/handlers/SpatialHandler.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// [MAP], [VERTICES], [POLYGONS], [SYMBOLS] and [TAGS]. Differences of form only:
//   * `std::transform(..., ::toupper)` is ASCII upper-casing (the C locale).
//   * [VERTICES] and [POLYGONS] hoist a reused token vector in upstream; here each row
//     gets its own. Which tokenizer a row goes through is kept exactly: the view
//     tokenizer unless the line holds a '"', the owned tokenizer otherwise (they differ
//     on empty CSV fields and quotes, so the choice can change a result).

#![allow(clippy::all)]

use crate::core::context::SimulationContext;
use crate::core::options::Bytes;
use crate::input::parse_utils::to_double;
use crate::input::tokenizer;

/// [MAP]: DIMENSIONS x1 y1 x2 y2 sets the extents; UNITS stores an upper-cased token.
pub fn handle_map(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.is_empty() {
            continue;
        }

        // Normalize keyword to uppercase for comparison
        let mut keyword = tok[0].clone();
        keyword.make_ascii_uppercase();

        if keyword == b"DIMENSIONS" && tok.len() >= 5 {
            ctx.spatial.map_x1 = to_double(&tok[1], 0.0);
            ctx.spatial.map_y1 = to_double(&tok[2], 0.0);
            ctx.spatial.map_x2 = to_double(&tok[3], 0.0);
            ctx.spatial.map_y2 = to_double(&tok[4], 0.0);
        } else if keyword == b"UNITS" && tok.len() >= 2 {
            ctx.spatial.map_units = tok[1].clone();
            ctx.spatial.map_units.make_ascii_uppercase();
        }
    }
}

/// The row split of [VERTICES] and [POLYGONS]: the view tokenizer for a line with no
/// double quote, the owned tokenizer otherwise.
fn split_row(line: &[u8]) -> Vec<Bytes> {
    if !line.contains(&b'"') {
        let mut tok: Vec<&[u8]> = Vec::new();
        tokenizer::tokenize_views_into(line, &mut tok);
        tok.iter().map(|t| t.to_vec()).collect()
    } else {
        tokenizer::tokenize(line)
    }
}

/// [VERTICES]: rebuild every link's intermediate vertices from the section, in order.
pub fn handle_vertices(ctx: &mut SimulationContext, lines: &[Bytes]) {
    let nl = ctx.link_names.size() as usize;
    if ctx.spatial.link_vertices_x.len() < nl {
        ctx.spatial.link_vertices_x.resize(nl, Vec::new());
    }
    if ctx.spatial.link_vertices_y.len() < nl {
        ctx.spatial.link_vertices_y.resize(nl, Vec::new());
    }

    // Rebuild from section contents to preserve deterministic per-link order.
    for xs in ctx.spatial.link_vertices_x.iter_mut() {
        xs.clear();
    }
    for ys in ctx.spatial.link_vertices_y.iter_mut() {
        ys.clear();
    }

    for line in lines {
        let tok = split_row(line);
        if tok.len() < 3 {
            continue;
        }

        let idx = ctx.link_names.find(&tok[0]);
        if idx < 0 {
            continue;
        }

        let uidx = idx as usize;
        if ctx.spatial.link_vertices_x.len() <= uidx {
            ctx.spatial.link_vertices_x.resize(uidx + 1, Vec::new());
            ctx.spatial.link_vertices_y.resize(uidx + 1, Vec::new());
        }

        ctx.spatial.link_vertices_x[uidx].push(to_double(&tok[1], 0.0));
        ctx.spatial.link_vertices_y[uidx].push(to_double(&tok[2], 0.0));
    }
}

/// [POLYGONS]: rebuild every subcatchment's outline from the section, in order.
pub fn handle_polygons(ctx: &mut SimulationContext, lines: &[Bytes]) {
    let ns = ctx.subcatch_names.size() as usize;
    if ctx.spatial.subcatch_polygon_x.len() < ns {
        ctx.spatial.subcatch_polygon_x.resize(ns, Vec::new());
    }
    if ctx.spatial.subcatch_polygon_y.len() < ns {
        ctx.spatial.subcatch_polygon_y.resize(ns, Vec::new());
    }

    // Rebuild from section contents to preserve deterministic per-subcatch order.
    for xs in ctx.spatial.subcatch_polygon_x.iter_mut() {
        xs.clear();
    }
    for ys in ctx.spatial.subcatch_polygon_y.iter_mut() {
        ys.clear();
    }

    for line in lines {
        let tok = split_row(line);
        if tok.len() < 3 {
            continue;
        }

        let idx = ctx.subcatch_names.find(&tok[0]);
        if idx < 0 {
            continue;
        }

        let uidx = idx as usize;
        if ctx.spatial.subcatch_polygon_x.len() <= uidx {
            ctx.spatial.subcatch_polygon_x.resize(uidx + 1, Vec::new());
            ctx.spatial.subcatch_polygon_y.resize(uidx + 1, Vec::new());
        }

        ctx.spatial.subcatch_polygon_x[uidx].push(to_double(&tok[1], 0.0));
        ctx.spatial.subcatch_polygon_y[uidx].push(to_double(&tok[2], 0.0));
    }
}

/// [SYMBOLS]: rain-gage symbol coordinates; the last row for a gage wins.
pub fn handle_symbols(ctx: &mut SimulationContext, lines: &[Bytes]) {
    let ng = ctx.gage_names.size() as usize;
    if ctx.spatial.gage_x.len() < ng {
        ctx.spatial.gage_x.resize(ng, 0.0);
    }
    if ctx.spatial.gage_y.len() < ng {
        ctx.spatial.gage_y.resize(ng, 0.0);
    }

    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 3 {
            continue;
        }

        let idx = ctx.gage_names.find(&tok[0]);
        if idx < 0 {
            continue;
        }

        let uidx = idx as usize;
        if ctx.spatial.gage_x.len() <= uidx {
            ctx.spatial.gage_x.resize(uidx + 1, 0.0);
            ctx.spatial.gage_y.resize(uidx + 1, 0.0);
        }

        ctx.spatial.gage_x[uidx] = to_double(&tok[1], 0.0);
        ctx.spatial.gage_y[uidx] = to_double(&tok[2], 0.0);
    }
}

/// [TAGS]: `NODE|LINK|SUBCATCH name tag`, stored per object index.
pub fn handle_tags(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 3 {
            continue;
        }

        // Normalize object type to uppercase
        let mut obj_type = tok[0].clone();
        obj_type.make_ascii_uppercase();

        let name: &[u8] = &tok[1];
        let tag: &Bytes = &tok[2];

        // Resolve name → index against the SoA stores. Tags are stored
        // per-NodeData/LinkData/SubcatchData index, not name-keyed, so
        // they survive a subsequent swmm_*_rename.
        if obj_type == b"NODE" {
            let idx = ctx.node_names.find(name);
            if idx >= 0 {
                let u = idx as usize;
                if u >= ctx.nodes.tags.len() {
                    ctx.nodes.tags.resize(u + 1, Vec::new());
                }
                ctx.nodes.tags[u] = tag.clone();
            }
        } else if obj_type == b"LINK" {
            let idx = ctx.link_names.find(name);
            if idx >= 0 {
                let u = idx as usize;
                if u >= ctx.links.tags.len() {
                    ctx.links.tags.resize(u + 1, Vec::new());
                }
                ctx.links.tags[u] = tag.clone();
            }
        } else if obj_type == b"SUBCATCH" {
            let idx = ctx.subcatch_names.find(name);
            if idx >= 0 {
                let u = idx as usize;
                if u >= ctx.subcatches.tags.len() {
                    ctx.subcatches.tags.resize(u + 1, Vec::new());
                }
                ctx.subcatches.tags[u] = tag.clone();
            }
        }
    }
}
