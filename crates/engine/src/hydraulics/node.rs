// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/Node.hpp and Node.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Node volume (namespace node), the part SWMMEngine::initGeometry reaches: getVolume with its
// storage-geometry accessor (storageGeom, isGeometric); and getSurfArea, which the FV mesh
// builder reaches from SWMMEngine::initHydraulics.
//
// Differences of form:
//   * `TableData* tables` (null allowed) is Option<&TableData>; `const NodeSubtypes* subs` is
//     Option<&NodeSubtypes>. table_getStorageVolume takes its Table non-const upstream but
//     does not write it.
//   * `ucf::Ucf[q][unit_sys]` is `UCF[q][unit_sys]`, indexed with unit_sys as given (0 or 1).
//
// getDepth (the inverse of getVolume) is at the end of the file (5 November 2026, for
// SWMMEngine::initialize).
//
// Not translated here (run-time; translated with the hydraulics phase): getPondedArea, getMaxOutflow, getOverflow, both storeVolume overloads,
// computeHeads, computeVolumes, computeOverflows and the batch kernels.

use crate::core::cmath::{cclamp, cmax};
use crate::core::constants::MIN_SURFAREA;
use crate::core::unit_conversion::{LENGTH, UCF, VOLUME};
use crate::data::node::{NodeData, NodeType, StorageShape};
use crate::data::node_subtypes::NodeSubtypes;
use crate::data::storage_geometry::storage_shape_is_geometric;
use crate::data::table::{table_get_storage_volume, table_lookup_ex, TableData};

/// Storage geometry accessor result (relational side-table).
#[derive(Debug, Clone, Copy)]
struct StorageGeom {
    curve: i32,
    a: f64,
    b: f64,
    c: f64,
    shape: StorageShape,
}

/// storageGeom: the storage row's (curve, a, b, c, shape), or the resize defaults
/// (-1 / 0 / 0 / 0 / FUNCTIONAL) without a side table or a storage row.
#[inline]
fn storage_geom(_nodes: &NodeData, subs: Option<&NodeSubtypes>, ui: usize) -> StorageGeom {
    if let Some(subs) = subs {
        let r = subs.storage_row(ui as i32);
        if r >= 0 {
            let ur = r as usize;
            return StorageGeom {
                curve: subs.storages.curve[ur],
                a: subs.storages.a[ur],
                b: subs.storages.b[ur],
                c: subs.storages.c[ur],
                shape: subs.storages.shape[ur],
            };
        }
    }
    StorageGeom { curve: -1, a: 0.0, b: 0.0, c: 0.0, shape: StorageShape::Functional }
}

/// isGeometric: CYLINDRICAL, CONICAL, PARABOLOID or PYRAMIDAL.
#[inline]
fn is_geometric(g: &StorageGeom) -> bool {
    storage_shape_is_geometric(g.shape)
}

/// node::storeVolume (Node.cpp): the water the continuity equation attributes to a node
/// (legacy's convention: 0 for a junction below its rim; storage and outfalls their volume).
pub fn store_volume_at(nodes: &NodeData, idx: i32, depth: f64, volume: f64) -> f64 {
    let ui = idx as usize;
    let v = cmax(volume, 0.0);
    if nodes.type_[ui] == NodeType::Storage || nodes.type_[ui] == NodeType::Outfall {
        return v;
    }
    let fd = nodes.full_depth[ui];
    if !(fd > 0.0) {
        return v;
    }
    let mut store = nodes.full_volume[ui] * (crate::core::cmath::cmin(depth, fd) / fd);
    if depth > fd {
        store += v;
    }
    store
}

/// node::storeVolume(nodes, idx): at the node's current depth and volume.
pub fn store_volume(nodes: &NodeData, idx: i32) -> f64 {
    let ui = idx as usize;
    store_volume_at(nodes, idx, nodes.depth[ui], nodes.volume[ui])
}

/// node::getVolume(nodes, idx, depth, tables, unit_sys, subs): volume (ft3) at `depth`.
pub fn get_volume(nodes: &NodeData, idx: i32, depth: f64, tables: Option<&TableData>, unit_sys: i32, subs: Option<&NodeSubtypes>) -> f64 {
    if depth <= 0.0 {
        return 0.0;
    }
    let ui = idx as usize;
    if nodes.storage_area[ui] > 0.0 {
        return nodes.storage_area[ui] * cclamp(depth, 0.0, nodes.full_depth[ui]);
    }

    if nodes.type_[ui] == NodeType::Storage {
        let lr = match subs {
            Some(s) => s.storage_row(idx),
            None => -1,
        };
        if lr >= 0 && !subs.unwrap().storages.lid_state[lr as usize].cells.is_empty() {
            let state = &subs.unwrap().storages.lid_state[lr as usize];
            if state.richards {
                let surface = &state.cells[0];
                return cclamp(depth - surface.bottom, 0., surface.top - surface.bottom) * surface.area * surface.porosity;
            }
            let mut volume = 0.0;
            for cell in &subs.unwrap().storages.lid_state[lr as usize].cells {
                let fraction = cclamp((depth - cell.bottom) / (cell.top - cell.bottom), 0.0, 1.0);
                volume += fraction * cell.geometric_volume * cmax(0.0, cell.porosity - cell.theta);
            }
            return volume;
        }
        // Clamp at fullDepth -> fullVolume (legacy node.c 909-910)
        if depth >= nodes.full_depth[ui] && nodes.full_volume[ui] > 0.0 {
            return nodes.full_volume[ui];
        }

        // Geometry fetched lazily, only after the clamp early-out.
        let g = storage_geom(nodes, subs, ui);

        if g.curve >= 0 {
            // Tabulated: trapezoidal integration of the area curve
            let ci = g.curve as usize;
            if let Some(tables) = tables {
                if ci < tables.tables.len() {
                    let ucf_len = UCF[LENGTH as usize][unit_sys as usize];
                    let ucf_vol = UCF[VOLUME as usize][unit_sys as usize];
                    let vol_disp = table_get_storage_volume(&tables.tables[ci], depth * ucf_len);
                    return vol_disp / ucf_vol;
                }
            }
            return 0.0;
        }
        // Geometric shapes: V = d*(c + d*(a/2 + d*b/3)), Horner form as legacy.
        if is_geometric(&g) {
            let d = depth * UCF[LENGTH as usize][unit_sys as usize];
            let v = d * (g.c + d * (g.a / 2.0 + d * g.b / 3.0));
            return v / UCF[VOLUME as usize][unit_sys as usize];
        }

        // Functional: V = a0*d + a1/(a2+1)*d^(a2+1), user units per call (legacy).
        let d = depth * UCF[LENGTH as usize][unit_sys as usize];
        let n = g.b + 1.0;
        let v = (g.c * d) + g.a / n * d.powf(n); // LIBM
        return v / UCF[VOLUME as usize][unit_sys as usize];
    }

    // JUNCTION / OUTFALL / DIVIDER: linear V = fullVolume * (d / fullDepth)
    let fd = nodes.full_depth[ui];
    if fd <= 0.0 {
        return 0.0;
    }
    let full_vol = if nodes.full_volume[ui] > 0.0 { nodes.full_volume[ui] } else { MIN_SURFAREA * fd };
    full_vol * (depth / fd)
}

/// getSurfArea: surface area (ft2) of node `idx` at `depth` (ft); 0 for a non-storage node.
pub fn get_surf_area(nodes: &NodeData, idx: i32, depth: f64, tables: Option<&TableData>, unit_sys: i32, subs: Option<&NodeSubtypes>) -> f64 {
    let ui = idx as usize;
    let us = unit_sys as usize;
    if nodes.storage_area[ui] > 0.0 {
        return if depth > nodes.full_depth[ui] { 0.0 } else { nodes.storage_area[ui] };
    }

    if nodes.type_[ui] == NodeType::Storage {
        let lr = match subs {
            Some(s) => s.storage_row(idx),
            None => -1,
        };
        if lr >= 0 && !subs.unwrap().storages.lid_state[lr as usize].cells.is_empty() {
            let st = &subs.unwrap().storages.lid_state[lr as usize];
            let cells = &st.cells;
            if st.richards {
                return cells[0].area * cells[0].porosity;
            }
            for c in cells {
                if depth >= c.bottom && depth < c.top {
                    return c.area * cmax(0.0, c.porosity - c.theta);
                }
            }
            let c = if depth <= 0.0 { cells.last().unwrap() } else { &cells[0] };
            return c.area * cmax(0.0, c.porosity - c.theta);
        }

        let g = storage_geom(nodes, subs, ui);
        // RAW storage-curve area (no MIN_SURFAREA clamp here; legacy storage_getSurfArea).
        if g.curve >= 0 {
            let ci = g.curve as usize;
            if let Some(t) = tables {
                if ci < t.tables.len() {
                    let ucf_len = UCF[LENGTH as usize][us];
                    // PARITY: extrapolating table_lookupEx, two successive divisions.
                    let area = table_lookup_ex(&t.tables[ci], depth * ucf_len);
                    return area / ucf_len / ucf_len;
                }
            }
            return 0.0;
        }
        let ucf_len = UCF[LENGTH as usize][us];

        // Geometric shapes: quadratic area = a0 + d*(a1 + d*a2), Horner form (node.c:993-999).
        if is_geometric(&g) {
            let d = depth * ucf_len;
            let area = g.c + d * (g.a + d * g.b);
            return area / ucf_len / ucf_len;
        }

        // Functional: area = a0 + a1 * d^a2 (user-unit coefficients, converted per call).
        let area = g.c + g.a * (depth * ucf_len).powf(g.b); // LIBM
        return area / ucf_len / ucf_len;
    }

    // Non-storage nodes: legacy node_getSurfArea returns 0.
    0.0
}

// ============================================================================
// Per-element: getDepth (inverse of getVolume). Translated 5 November 2026 for
// SWMMEngine::initialize (the LID-node seeding loop); the run-time callers use it too.
// ============================================================================

/// node::getDepth(nodes, idx, volume, tables, unit_sys, subs): depth (ft) holding `volume`.
pub fn get_depth(nodes: &NodeData, idx: i32, volume: f64, tables: Option<&TableData>, unit_sys: i32, subs: Option<&NodeSubtypes>) -> f64 {
    let mut volume = volume;
    if volume < 0.0 {
        return 0.0;
    }
    let ui = idx as usize;
    if nodes.storage_area[ui] > 0.0 {
        return cclamp(volume / nodes.storage_area[ui], 0.0, nodes.full_depth[ui]);
    }

    if nodes.type_[ui] == NodeType::Storage {
        let lr = match subs {
            Some(s) => s.storage_row(idx),
            None => -1,
        };
        if lr >= 0 && !subs.unwrap().storages.lid_state[lr as usize].cells.is_empty() {
            let st = &subs.unwrap().storages.lid_state[lr as usize];
            let cells = &st.cells;
            if st.richards {
                let surface = &cells[0];
                return crate::core::cmath::cmin(surface.top, surface.bottom + volume / (surface.area * surface.porosity));
            }
            for c in cells.iter().rev() {
                let capacity = c.geometric_volume * cmax(0.0, c.porosity - c.theta);
                if capacity > 0.0 && volume < capacity {
                    return c.bottom + volume / capacity * (c.top - c.bottom);
                }
                volume -= capacity;
            }
            return nodes.full_depth[ui];
        }
        if volume <= 0.0 {
            return 0.0;
        }
        let fd = nodes.full_depth[ui];
        let fv = nodes.full_volume[ui];
        if fv > 0.0 && volume >= fv {
            return fd;
        }

        // Geometry fetched lazily, only after the clamp early-out (as legacy did).
        let g = storage_geom(nodes, subs, ui);

        if g.curve >= 0 {
            // Tabulated: quadratic solve per interval (legacy table_getStorageDepth).
            let ci = g.curve as usize;
            if let Some(t) = tables {
                if ci < t.tables.len() {
                    let ucf_len = UCF[LENGTH as usize][unit_sys as usize];
                    let ucf_vol = UCF[VOLUME as usize][unit_sys as usize];
                    let vol_disp = volume * ucf_vol;
                    let d_disp = crate::data::table::table_get_storage_depth(&t.tables[ci], vol_disp);
                    return d_disp / ucf_len;
                }
            }
            if fv > 0.0 {
                return fd * (volume / fv);
            }
            return 0.0;
        }

        // Functional: V = a0*d + a1/(a2+1) * d^(a2+1), in USER units (legacy storage_getDepth).
        let a0 = g.c;
        let a1 = g.a;
        let a2 = g.b;
        let ucf_len = UCF[LENGTH as usize][unit_sys as usize];
        let ucf_vol = UCF[VOLUME as usize][unit_sys as usize];
        let v = volume * ucf_vol;
        let mut d: f64;

        let newton_solve = |seed: f64| -> f64 {
            let mut d = seed;
            crate::math::find_root::newton(0.0, fd * ucf_len, &mut d, 0.001, |y: f64, f: &mut f64, df: &mut f64| {
                let y_ft = y / ucf_len;
                *f = get_volume(nodes, idx, y_ft, tables, unit_sys, subs) * ucf_vol - v;
                *df = get_surf_area(nodes, idx, y_ft, tables, unit_sys, subs) * ucf_vol / ucf_len;
            });
            d
        };

        if is_geometric(&g) {
            match g.shape {
                StorageShape::Cylindrical => {
                    d = v / a0;
                }
                StorageShape::Paraboloid => {
                    d = (2.0 * v / a1).sqrt();
                }
                _ => {
                    d = newton_solve(v / a0);
                }
            }
        } else if a2 == 0.0 {
            d = v / (a0 + a1);
        } else if a0 == 0.0 {
            let e = 1.0 / (a2 + 1.0);
            d = (v / (a1 * e)).powf(e); // LIBM
        } else if a2 == 1.0 && a1 > 0.0 {
            d = ((a0 * a0 + 2. * a1 * v).sqrt() - a0) / a1;
        } else {
            d = newton_solve(v / (a0 + a1));
        }
        d /= ucf_len;
        return crate::core::cmath::cmin(d, fd);
    }

    // JUNCTION / OUTFALL / DIVIDER: V = MIN_SURFAREA * d -> d = V / MIN_SURFAREA
    volume / MIN_SURFAREA
}
