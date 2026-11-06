// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/Node.cpp (HydroCouple/openswmm.engine,
// branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache
// License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The node functions the routing step calls that hydraulics/node.rs (getVolume, getSurfArea)
// does not hold: getDepth (with a copy of node.rs's file-local storageGeom), getPondedArea,
// getMaxOutflow, getOverflow.
//
// MERGE NOTE: branch rs6-init adds a getDepth to hydraulics/node.rs for SWMMEngine::initialize.
// Both translate Node.cpp's getDepth line for line; keep one.

use crate::core::cmath::{cclamp, cmax, cmin};
use crate::core::constants::MIN_SURFAREA;
use crate::core::unit_conversion::{LENGTH, UCF, VOLUME};
use crate::data::node::{NodeData, NodeType, StorageShape};
use crate::data::node_subtypes::NodeSubtypes;
use crate::data::storage_geometry::storage_shape_is_geometric;
use crate::data::table::{table_get_storage_depth, TableData};
use crate::hydraulics::node::{get_surf_area, get_volume};
use crate::math::find_root;

/// node::getDepth: the depth (ft) of node `idx` holding `volume` (ft3).
pub fn get_depth(nodes: &NodeData, idx: i32, mut volume: f64, tables: Option<&TableData>, unit_sys: i32, subs: Option<&NodeSubtypes>) -> f64 {
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
        if lr >= 0 {
            let st = &subs.unwrap().storages.lid_state[lr as usize];
            if !st.cells.is_empty() {
                let cells = &st.cells;
                if st.richards {
                    let surface = &cells[0];
                    return cmin(surface.top, surface.bottom + volume / (surface.area * surface.porosity));
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
        }
        if volume <= 0.0 {
            return 0.0;
        }
        let fd = nodes.full_depth[ui];
        let fv = nodes.full_volume[ui];
        if fv > 0.0 && volume >= fv {
            return fd;
        }

        // storageGeom
        let (curve, ga, gb, gc, shape) = match subs {
            Some(s) if s.storage_row(idx) >= 0 => {
                let ur = s.storage_row(idx) as usize;
                (s.storages.curve[ur], s.storages.a[ur], s.storages.b[ur], s.storages.c[ur], s.storages.shape[ur])
            }
            _ => (-1, 0.0, 0.0, 0.0, StorageShape::Functional),
        };

        let us = unit_sys as usize;
        if curve >= 0 {
            let ci = curve as usize;
            if let Some(t) = tables {
                if ci < t.tables.len() {
                    let ucf_len = UCF[LENGTH as usize][us];
                    let ucf_vol = UCF[VOLUME as usize][us];
                    let vol_disp = volume * ucf_vol;
                    let d_disp = table_get_storage_depth(&t.tables[ci], vol_disp);
                    return d_disp / ucf_len;
                }
            }
            if fv > 0.0 {
                return fd * (volume / fv);
            }
            return 0.0;
        }

        let a0 = gc;
        let a1 = ga;
        let a2 = gb;
        let ucf_len = UCF[LENGTH as usize][us];
        let ucf_vol = UCF[VOLUME as usize][us];
        let v = volume * ucf_vol;
        let mut d;

        let newton_solve = |seed: f64| -> f64 {
            let mut d = seed;
            find_root::newton(0.0, fd * ucf_len, &mut d, 0.001, |y: f64, f: &mut f64, df: &mut f64| {
                let y_ft = y / ucf_len;
                *f = get_volume(nodes, idx, y_ft, tables, unit_sys, subs) * ucf_vol - v;
                *df = get_surf_area(nodes, idx, y_ft, tables, unit_sys, subs) * ucf_vol / ucf_len;
            });
            d
        };

        if storage_shape_is_geometric(shape) {
            d = match shape {
                StorageShape::Cylindrical => v / a0,
                StorageShape::Paraboloid => (2.0 * v / a1).sqrt(),
                _ => newton_solve(v / a0),
            };
        } else if a2 == 0.0 {
            d = v / (a0 + a1);
        } else if a0 == 0.0 {
            let e = 1.0 / (a2 + 1.0);
            d = (v / (a1 * e)).powf(e); // LIBM pow
        } else if a2 == 1.0 && a1 > 0.0 {
            d = ((a0 * a0 + 2. * a1 * v).sqrt() - a0) / a1;
        } else {
            d = newton_solve(v / (a0 + a1));
        }
        d /= ucf_len;
        return cmin(d, fd);
    }

    volume / MIN_SURFAREA
}

/// node::getPondedArea
pub fn get_ponded_area(nodes: &NodeData, idx: i32, depth: f64, tables: Option<&TableData>, unit_sys: i32, subs: Option<&NodeSubtypes>) -> f64 {
    let ui = idx as usize;
    if depth <= nodes.full_depth[ui] || nodes.ponded_area[ui] == 0.0 {
        return get_surf_area(nodes, idx, depth, tables, unit_sys, subs);
    }
    let mut a = nodes.ponded_area[ui];
    if a <= 0.0 {
        a = get_surf_area(nodes, idx, nodes.full_depth[ui], tables, unit_sys, subs);
    }
    a
}

/// node::getMaxOutflow: caps only nodes with a legacy full volume.
pub fn get_max_outflow(nodes: &NodeData, idx: i32, mut q: f64, dt: f64) -> f64 {
    let ui = idx as usize;
    if nodes.rpt_full_volume[ui] > 0.0 {
        let q_max = nodes.inflow[ui] + nodes.old_volume[ui] / dt;
        q = cmin(q, q_max);
    }
    cmax(0.0, q)
}

/// node::getOverflow
pub fn get_overflow(new_volume: f64, full_volume: f64, dt: f64) -> f64 {
    if new_volume > full_volume && dt > 0.0 {
        return (new_volume - full_volume) / dt;
    }
    0.0
}
