// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/Constants.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Literals exactly as upstream spells them; computed ones keep upstream's expression
// (5.0 / 3.0 and 1.0 / 12.0 are evaluated in double, as constexpr evaluates them).
#![allow(clippy::approx_constant)]

pub const GRAVITY: f64 = 32.2;
pub const SQRT_GRAVITY: f64 = 5.67450438786;
pub const INV_SQRT_GRAVITY: f64 = 0.17622692584;
pub const PI: f64 = 3.141592654;
pub const PHI: f64 = 1.486;
pub const MEXP: f64 = 5.0 / 3.0;
pub const VISCOS: f64 = 1.1e-5;
pub const SEC_PER_DAY: f64 = 86400.0;
pub const FT_PER_IN: f64 = 1.0 / 12.0;
pub const IN_PER_FT: f64 = 12.0;
pub const FUDGE: f64 = 0.0001;
pub const MIN_SURFAREA: f64 = 12.566;
pub const FLOW_TOL: f64 = 0.00001;
pub const MIN_DELTA_Z: f64 = 0.001;
pub const XTOL: f64 = 0.001;
pub const TINY: f64 = 1.0e-6;
pub const MISSING: f64 = -1.0e10;
pub const OMEGA: f64 = 0.5;
pub const DEFAULT_HEAD_TOL: f64 = 0.005;
pub const DEFAULT_MAX_TRIALS: i32 = 8;
pub const MAX_VELOCITY: f64 = 50.0;
pub const MIN_TIMESTEP: f64 = 0.001;
pub const EXTRAN_CROWN_CUTOFF: f64 = 0.96;
pub const SLOT_CROWN_CUTOFF: f64 = 0.985257;
pub const SLOT_WIDTH_FACTOR: f64 = 0.001;

pub mod ode {
    pub const MAX_STEPS: i32 = 10000;
    pub const ODE_TINY: f64 = 1.0e-30;
    pub const SAFETY: f64 = 0.9;
    pub const PGROW: f64 = -0.2;
    pub const PSHRNK: f64 = -0.25;
    pub const ERRCON: f64 = 1.89e-4;
    pub const GWTOL: f64 = 0.0001;
    pub const ODETOL: f64 = 0.0001;
}

pub const MIN_LINK_FLOW: f64 = 0.001;
pub const RECT_ALFMAX: f64 = 0.97;
pub const TRAP_ALFMAX: f64 = 0.98;
pub const N_TRANSECT_TBL: i32 = 51;
pub const MAX_SHAPES: i32 = 26;
pub const N_FLOW_CLASSES: i32 = 7;
pub const MAX_PAST_RAIN: i32 = 48;
pub const N_TIME_BINS: i32 = 5;
pub const DATE_DELTA: i32 = 693594;
