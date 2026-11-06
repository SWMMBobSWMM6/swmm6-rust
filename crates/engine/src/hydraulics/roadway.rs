// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/Roadway.cpp and Roadway.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// ROADWAY weir flow (legacy roadway.c, FHWA HDS-5 coefficients).

use crate::core::cmath::cmax;

pub const SURFACE_PAVED: i32 = 1;

const CR_LOW_PAVED: [[f64; 2]; 4] = [[0.0, 2.85], [0.2, 2.95], [0.7, 3.03], [4.0, 3.05]];
const CR_LOW_GRAVEL: [[f64; 2]; 8] = [[0.0, 2.5], [0.5, 2.7], [1.0, 2.8], [1.5, 2.9], [2.0, 2.98], [2.5, 3.02], [3.0, 3.03], [4.0, 3.05]];
const CR_HIGH_PAVED: [[f64; 2]; 2] = [[0.15, 3.05], [0.25, 3.10]];
const CR_HIGH_GRAVEL: [[f64; 2]; 2] = [[0.15, 2.95], [0.30, 3.10]];
const KT_PAVED: [[f64; 2]; 9] = [[0.8, 1.0], [0.85, 0.98], [0.90, 0.92], [0.93, 0.85], [0.95, 0.80], [0.97, 0.70], [0.98, 0.60], [0.99, 0.50], [1.00, 0.40]];
const KT_GRAVEL: [[f64; 2]; 12] = [[0.75, 1.00], [0.80, 0.985], [0.83, 0.97], [0.86, 0.93], [0.89, 0.90], [0.90, 0.87], [0.92, 0.80], [0.94, 0.70], [0.96, 0.60], [0.98, 0.50], [0.99, 0.40], [1.00, 0.24]];

/// legacy getY: linear interpolation in a two-column table, clamped at its ends.
fn get_y(x: f64, table: &[[f64; 2]]) -> f64 {
    let n = table.len();
    if x <= table[0][0] {
        return table[0][1];
    }
    if x >= table[n - 1][0] {
        return table[n - 1][1];
    }
    for i in 1..n {
        if x <= table[i][0] {
            let x1 = table[i - 1][0];
            let dx = table[i][0] - x1;
            let y1 = table[i - 1][1];
            let dy = table[i][1] - y1;
            return y1 + (x - x1) * dy / dx;
        }
    }
    table[n - 1][1]
}

/// legacy getCd
fn get_cd(h_wr: f64, ht: f64, road_width: f64, road_surf: i32) -> f64 {
    let mut k_t = 1.0;
    if h_wr <= 0.0 {
        return 0.0;
    }
    let h_l = h_wr / road_width;
    let c_r = if h_l <= 0.15 {
        if road_surf == SURFACE_PAVED {
            get_y(h_wr, &CR_LOW_PAVED)
        } else {
            get_y(h_wr, &CR_LOW_GRAVEL)
        }
    } else if road_surf == SURFACE_PAVED {
        get_y(h_l, &CR_HIGH_PAVED)
    } else {
        get_y(h_l, &CR_HIGH_GRAVEL)
    };
    if ht > 0.0 {
        let ht_h = ht / h_wr;
        k_t = if road_surf == SURFACE_PAVED { get_y(ht_h, &KT_PAVED) } else { get_y(ht_h, &KT_GRAVEL) };
    }
    c_r * k_t
}

/// roadway::getInflow: the signed flow over a roadway weir.
#[allow(clippy::too_many_arguments)]
pub fn get_inflow(dir: f64, h_road: f64, h1: f64, h2: f64, cd_user: f64, si: bool, road_width: f64, road_surf: i32, length: f64, dqdh: &mut f64, depth: &mut f64) -> f64 {
    const FUDGE: f64 = 0.0001;
    let mut q = 0.0;
    *dqdh = 0.0;
    let mut c_d = cd_user;
    if si {
        c_d /= 0.552;
    }
    let use_variable_cd = road_width > 0.0 && road_surf >= 1;
    let h_wr = h1 - h_road;
    let ht = h2 - h_road;
    if h_wr > FUDGE {
        if use_variable_cd {
            c_d = get_cd(h_wr, ht, road_width, road_surf);
        }
        q = c_d * length * h_wr.powf(1.5); // LIBM pow
        *dqdh = 1.5 * q / h_wr;
    }
    *depth = cmax(h1 - h_road, 0.0);
    dir * q
}
