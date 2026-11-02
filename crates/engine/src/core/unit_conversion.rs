// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/UnitConversion.hpp and
// src/engine/core/UnitConversion.cpp (HydroCouple/openswmm.engine, branch swmm6_rel,
// commit 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the Apache License,
// Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Upstream's namespace `openswmm::ucf` is this module. Differences in form, never in
// result:
//
//   * `enum Quantity` is unscoped and `UCF` takes it as `int`, so `Quantity` is an alias
//     for i32 and RAINFALL..FLOW are `pub const` items with upstream's values.
//   * Names follow Rust case: `Ucf`/`Qcf`/`Ucf_inv`/`Qcf_inv`/`FlowUnitWords` are
//     `UCF`/`QCF`/`UCF_INV`/`QCF_INV`/`FLOW_UNIT_WORDS`; `UCF()`/`UCF_inv()`/
//     `getUnitSystem()` are `ucf()`/`ucf_inv()`/`get_unit_system()`.
//   * Every literal is spelled as upstream spells it. The reciprocal tables are built by
//     `const fn`, as upstream's are by `constexpr`: one IEEE division `1.0 / x` each,
//     correctly rounded in both compilers.
//   * Label words (`const char*` to static literals) are `&'static [u8]`.
//   * `SimulationOptions::flow_units` is the closed `FlowUnits` enum, so the `fu` range
//     checks below can never fail here; they are kept as upstream writes them.

use crate::core::options::SimulationOptions;

// ============================================================================
// Quantity codes — matching legacy enums.h
// ============================================================================

/// Upstream `enum Quantity`.
pub type Quantity = i32;

/// Divide in/hr (US) or mm/hr (SI) by this → ft/sec
pub const RAINFALL: Quantity = 0;
/// Divide in (US) or mm (SI) by this → ft
pub const RAINDEPTH: Quantity = 1;
/// Divide in/day (US) or mm/day (SI) by this → ft/sec
pub const EVAPRATE: Quantity = 2;
/// Divide ft (US) or m (SI) by this → ft
pub const LENGTH: Quantity = 3;
/// Divide ac (US) or ha (SI) by this → ft²
pub const LANDAREA: Quantity = 4;
/// Divide ft³ (US) or m³ (SI) by this → ft³
pub const VOLUME: Quantity = 5;
/// Divide mph (US) or km/hr (SI) by this → mph
pub const WINDSPEED: Quantity = 6;
/// Legacy scale entry only; absolute temperatures require offset conversion
pub const TEMPERATURE: Quantity = 7;
/// Divide lb (US) or kg (SI) by this → mg
pub const MASS: Quantity = 8;
/// Divide cfs/ac (US) or cms/ha (SI) by this → ft/sec
pub const GWFLOW: Quantity = 9;
/// Divide display flow units by this → cfs
pub const FLOW: Quantity = 10;

// ============================================================================
// Conversion factor tables — identical to legacy Ucf[10][2] and Qcf[6]
// ============================================================================

/// Ucf[quantity][unit_system] — multiply project units by 1/Ucf to get internal
/// Or equivalently: internal = project / UCF; project = internal * UCF
pub const UCF: [[f64; 2]; 10] = [
    //  US           SI
    [43200.0, 1097280.0],    // RAINFALL: in/hr, mm/hr → ft/sec
    [12.0, 304.8],           // RAINDEPTH: in, mm → ft
    [1036800.0, 26334720.0], // EVAPRATE: in/day, mm/day → ft/sec
    [1.0, 0.3048],           // LENGTH: ft, m → ft
    [2.2956e-5, 0.92903e-5], // LANDAREA: ac, ha → ft²
    [1.0, 0.02832],          // VOLUME: ft³, m³ → ft³
    [1.0, 1.608],            // WINDSPEED: mph, km/hr → mph
    [1.0, 1.8],              // TEMPERATURE: °F, °C → °F
    [2.203e-6, 1.0e-6],      // MASS: lb, kg → mg
    [43560.0, 3048.0],       // GWFLOW: cfs/ac, cms/ha → ft/sec
];

/// Qcf[flow_units] — conversion from display flow to cfs
pub const QCF: [f64; 6] = [
    1.0,     // CFS → cfs
    448.831, // GPM → cfs
    0.64632, // MGD → cfs
    0.02832, // CMS → cfs
    28.317,  // LPS → cfs
    2.4466,  // MLD → cfs
];

// ============================================================================
// Precomputed reciprocal tables — multiply instead of divide
// ============================================================================
//
// Display → internal conversion is `internal = display / Ucf` everywhere. To
// let hot/peripheral conversion sites multiply by a hoisted loop-invariant
// instead of dividing, expose the compile-time reciprocals. These are derived
// from Ucf/Qcf so they can never drift out of sync.
//
// Convention:  internal = display * Ucf_inv[q][us]   (display → internal)
//              display  = internal * Ucf[q][us]       (internal → display)

pub mod detail {
    use super::{QCF, UCF};

    /// Compile-time reciprocals of every quantity/system conversion: entry 1/Ucf[q][system].
    pub const fn make_ucf_inv() -> [[f64; 2]; 10] {
        let mut r = [[0.0f64; 2]; 10];
        let mut q = 0;
        while q < 10 {
            let mut u = 0;
            while u < 2 {
                r[q][u] = 1.0 / UCF[q][u];
                u += 1;
            }
            q += 1;
        }
        r
    }

    /// Compile-time reciprocals of flow conversion factors: entry 1/Qcf[flow_units].
    pub const fn make_qcf_inv() -> [f64; 6] {
        let mut r = [0.0f64; 6];
        let mut f = 0;
        while f < 6 {
            r[f] = 1.0 / QCF[f];
            f += 1;
        }
        r
    }
}

/// Reciprocal of Ucf[quantity][unit_system] — multiply display by this → internal.
pub const UCF_INV: [[f64; 2]; 10] = detail::make_ucf_inv();

/// Reciprocal of Qcf[flow_units] — multiply display flow by this → cfs.
pub const QCF_INV: [f64; 6] = detail::make_qcf_inv();

// ============================================================================
// Convenience constants — common fixed conversions not unit-system dependent
// ============================================================================

/// 1 acre = 43560 ft²
pub const ACRES_TO_FT2: f64 = 43560.0;
/// ft³ → million gallons
pub const FT3_TO_MGAL: f64 = 7.48052e-6;
/// seconds per day
pub const SEC_PER_DAY: f64 = 86400.0;

// ============================================================================
// UCF function — matching legacy UCF()
// ============================================================================

/// Determine unit system (0=US, 1=SI) from flow units: CFS/GPM/MGD → US (0),
/// CMS/LPS/MLD → SI (1).
pub fn get_unit_system(flow_units: i32) -> i32 {
    // CFS=0, GPM=1, MGD=2 → US (0)
    // CMS=3, LPS=4, MLD=5 → SI (1)
    if flow_units >= 3 {
        1
    } else {
        0
    }
}

/// Unit conversion factor for a quantity: Ucf[quantity][unit_system] for a non-flow
/// quantity, Qcf[flow_units] for FLOW. internal = display / ucf; display = internal * ucf.
#[allow(clippy::manual_range_contains)]
pub fn ucf(quantity: i32, opts: &SimulationOptions) -> f64 {
    let fu = opts.flow_units as i32;
    if quantity < FLOW {
        let us = get_unit_system(fu);
        if quantity >= 0 && quantity <= 9 {
            return UCF[quantity as usize][us as usize];
        }
        return 1.0;
    }
    // FLOW
    if fu >= 0 && fu <= 5 {
        return QCF[fu as usize];
    }
    1.0
}

/// Reciprocal of `ucf`, from the precomputed tables: internal = display * ucf_inv.
#[allow(clippy::manual_range_contains)]
pub fn ucf_inv(quantity: i32, opts: &SimulationOptions) -> f64 {
    let fu = opts.flow_units as i32;
    if quantity < FLOW {
        let us = get_unit_system(fu);
        if quantity >= 0 && quantity <= 9 {
            return UCF_INV[quantity as usize][us as usize];
        }
        return 1.0;
    }
    // FLOW
    if fu >= 0 && fu <= 5 {
        return QCF_INV[fu as usize];
    }
    1.0
}

/// Legacy FlowUnitWords — index by FlowUnits enum (CFS,GPM,MGD,CMS,LPS,MLD).
pub const FLOW_UNIT_WORDS: [&[u8]; 6] = [b"CFS", b"GPM", b"MGD", b"CMS", b"LPS", b"MLD"];

// ============================================================================
// DisplayUnits — single source of truth for internal → display conversion
// ============================================================================
//
// Built once from SimulationOptions via DisplayUnits::from(). Carries the
// internal→display multipliers (display = internal * factor) plus the matching
// unit-label words, so every output / report consumer shares one definition
// instead of re-deriving Ucf / Qcf inline. Mirrors legacy statsrpt.c / report.c
// UCF()-based output and the VolUnitsWords / PondingUnitsWords tables.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DisplayUnits {
    /// 0 = US, 1 = SI
    pub unit_system: i32,
    /// FlowUnits enum value
    pub flow_units: i32,

    // --- internal → display multipliers (display = internal * factor) ---
    /// cfs  → display flow (Qcf)
    pub flow: f64,
    /// ft   → ft | m
    pub length: f64,
    /// ft³  → ft³ | m³
    pub volume: f64,
    /// ft/s → in/hr | mm/hr
    pub rainfall: f64,
    /// ft   → in | mm (also runoff-depth columns)
    pub raindepth: f64,
    /// ft/s → in/day | mm/day
    pub evaprate: f64,
    /// ft²  → ac | ha
    pub landarea: f64,
    /// ft³  → 10^6 gal | 10^6 ltr
    pub mvol: f64,
    /// ft³  → acre-ft | hectare-m
    pub landvol: f64,

    // --- unit-label words (legacy statsrpt.c / report.c) ---
    /// CFS | CMS | ...
    pub flow_word: &'static [u8],
    /// depth / HGL / surcharge column
    pub length_word: &'static [u8],
    /// rainfall / runoff depth (in | mm)
    pub depth_word: &'static [u8],
    /// link velocity (ft/sec | m/sec)
    pub vel_word: &'static [u8],
    /// 10^6 gal | 10^6 ltr
    pub mvol_word: &'static [u8],
    /// acre-feet | hectare-m
    pub landvol_word: &'static [u8],
    /// 1000 ft3 | 1000 m3
    pub storage_vol_word: &'static [u8],
}

impl Default for DisplayUnits {
    fn default() -> Self {
        DisplayUnits {
            unit_system: 0,
            flow_units: 0,
            flow: 1.0,
            length: 1.0,
            volume: 1.0,
            rainfall: 1.0,
            raindepth: 1.0,
            evaprate: 1.0,
            landarea: 1.0,
            mvol: 1.0,
            landvol: 1.0,
            flow_word: b"CFS",
            length_word: b"Feet",
            depth_word: b"in",
            vel_word: b"ft/sec",
            mvol_word: b"10^6 gal",
            landvol_word: b"acre-feet",
            storage_vol_word: b"1000 ft3",
        }
    }
}

impl DisplayUnits {
    /// Build display factors and static labels from the options' flow units (US/SI follows
    /// flow_units alone); an out-of-range flow_units falls back to CFS/US.
    #[allow(clippy::manual_range_contains)]
    pub fn from(opts: &SimulationOptions) -> DisplayUnits {
        let mut d = DisplayUnits::default();
        let mut fu = opts.flow_units as i32;
        if fu < 0 || fu > 5 {
            fu = 0;
        }
        let us = get_unit_system(fu);
        let si = us == 1;

        d.unit_system = us;
        d.flow_units = fu;

        // internal → display multipliers
        d.flow = QCF[fu as usize];
        d.length = UCF[LENGTH as usize][us as usize];
        d.volume = UCF[VOLUME as usize][us as usize];
        d.rainfall = UCF[RAINFALL as usize][us as usize];
        d.raindepth = UCF[RAINDEPTH as usize][us as usize];
        d.evaprate = UCF[EVAPRATE as usize][us as usize];
        d.landarea = UCF[LANDAREA as usize][us as usize];
        // ft³ → 10^6 gal | 10^6 ltr  (legacy statsrpt.c Vcf: 7.48/1e6 | 28.3168/1e6)
        d.mvol = if si { 28.3168 / 1.0e6 } else { FT3_TO_MGAL };
        // ft³ → acre-ft | hectare-m  (legacy report.c x*UCF(LENGTH)*UCF(LANDAREA))
        d.landvol = if si { 0.0283168 / 1.0e4 } else { 1.0 / ACRES_TO_FT2 };

        // unit-label words
        d.flow_word = FLOW_UNIT_WORDS[fu as usize];
        d.length_word = if si { b"Meters" } else { b"Feet" };
        d.depth_word = if si { b"mm" } else { b"in" };
        d.vel_word = if si { b"m/sec" } else { b"ft/sec" };
        d.mvol_word = if si { b"10^6 ltr" } else { b"10^6 gal" };
        d.landvol_word = if si { b"hectare-m" } else { b"acre-feet" };
        d.storage_vol_word = if si { b"1000 m3" } else { b"1000 ft3" };
        d
    }

    /// Convert an internal Fahrenheit temperature to this bundle's display scale: tF for
    /// unit_system == 0, otherwise (5/9)*(tF-32) in deg C. No finite/range validation.
    pub fn temperature(&self, t_f: f64) -> f64 {
        if self.unit_system == 0 {
            t_f
        } else {
            (5.0 / 9.0) * (t_f - 32.0)
        }
    }
}
