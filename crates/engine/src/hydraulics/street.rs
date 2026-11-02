// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/Street.hpp and
// src/engine/hydraulics/Street.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit
// 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Street cross-section geometry (gutter + road crown + backing), a literal port of legacy
// transect_createStreetTransect. Differences in form: buildTransect is build_transect, and
// the point-appending lambda is a nested function; Default is written by hand because the
// roughness starts at 0.013 and sides at 2.

use crate::core::cmath::cmax;
use crate::hydraulics::transect::build_tables;
use crate::hydraulics::transect_types::TransectData;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StreetParams {
    /// Distance curb to crown (ft)
    pub width: f64,
    /// Curb elevation (ft)
    pub curb_height: f64,
    /// Transverse slope (fraction)
    pub slope: f64,
    /// Manning's n (road surface)
    pub roughness: f64,
    /// Depressed gutter depth (ft)
    pub gutter_depression: f64,
    /// Depressed gutter width (ft)
    pub gutter_width: f64,
    /// 1=half street, 2=full street
    pub sides: i32,
    /// Backing width (ft)
    pub back_width: f64,
    /// Backing slope (fraction)
    pub back_slope: f64,
    /// Backing Manning's n
    pub back_roughness: f64,
}

impl Default for StreetParams {
    fn default() -> Self {
        StreetParams {
            width: 0.0,
            curb_height: 0.0,
            slope: 0.0,
            roughness: 0.013,
            gutter_depression: 0.0,
            gutter_width: 0.0,
            sides: 2,
            back_width: 0.0,
            back_slope: 0.0,
            back_roughness: 0.0,
        }
    }
}

/// Build a transect from street parameters: its station/elevation points, Manning's n and
/// bank stations, then the tables (build_tables without end walls).
pub fn build_transect(sp: &StreetParams, td: &mut TransectData) {
    // Literal port of legacy transect_createStreetTransect (transect.c:625).
    //
    //   Point 0 = top of backing    Point 3 = bottom of depressed gutter
    //   Point 1 = top of curb       Point 4 = top of depressed gutter (crown)
    //   Point 2 = bottom of curb    Point 5 = street crown wall
    //
    //   y3 = Hdep + Sx*Wg   — the gutter's outer edge is the depression PLUS
    //                         the road cross-slope over the gutter width,
    //   y4 = y3 + Sx*(W-Wg) — and the crown builds on y3, not on Hdep.
    td.stations.clear();
    td.elevations.clear();

    let w1 = sp.back_width; // backing
    let w2 = sp.gutter_width; // depressed gutter
    let w3 = sp.width; // curb to crown
    let w4 = w3 - w2; // road, less the gutter
    let y3 = sp.gutter_depression + sp.slope * w2;
    let y1 = sp.curb_height + sp.gutter_depression;
    let y4 = y3 + sp.slope * w4;
    // The backing's top is raised to the crown elevation when the crown is
    // higher (legacy `ymax = MAX(ymax, y4)` before Elev[0] is assigned).
    let ymax = cmax(sp.back_slope * w1 + y1, y4);

    fn pt(td: &mut TransectData, x: f64, y: f64) {
        td.stations.push(x);
        td.elevations.push(y);
    }

    pt(td, 0.0, ymax); // 0: top of backing
    pt(td, w1, y1); // 1: top of curb
    pt(td, w1, 0.0); // 2: bottom of curb (the thalweg)
    pt(td, w1 + w2, y3); // 3: outer edge of the depressed gutter
    pt(td, w1 + w3, y4); // 4: crown

    if sp.sides == 1 {
        // A half street is closed by a vertical wall at the crown.
        pt(td, w1 + w3, ymax);
    } else {
        // The right side mirrors the left.
        pt(td, w1 + w3 + w4, y3);
        pt(td, w1 + w3 + w4 + w2, 0.0);
        pt(td, w1 + w3 + w4 + w2, y1);
        pt(td, w1 + w3 + w4 + w2 + w1, ymax);
    }

    // Manning's n and the bank stations (legacy: a street with no backing is
    // ALL main channel, keyed on backWidth — not on whether nBack was given).
    td.n_channel = sp.roughness;
    if sp.back_width == 0.0 {
        td.n_left = sp.roughness;
        td.n_right = sp.roughness;
        td.x_left_bank = td.stations[0];
        td.x_right_bank = td.stations[td.stations.len() - 1];
    } else {
        td.n_left = sp.back_roughness;
        td.n_right = sp.back_roughness;
        td.x_left_bank = td.stations[1];
        td.x_right_bank = if sp.sides == 2 {
            td.stations[td.stations.len() - 2]
        } else {
            td.stations[td.stations.len() - 1]
        };
    }

    build_tables(td, false);
}
