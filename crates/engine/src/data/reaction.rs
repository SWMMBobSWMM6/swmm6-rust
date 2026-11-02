// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/ReactionData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Multispecies reaction system data (EPANET-MSX conventions), SoA, hot/cold split.
// Rust has no nested type declarations, so ReactionData::X is ReactionDataX
// (ParameterOverride, WallArea, LidLaw, LidState, SourcePattern, Source, SubsurfaceLaw,
// SurfaceLaw, MsxRoutingLedger). ParameterOverride, WallArea and Source leave some members
// without initializers in upstream; their Default here is zero, which is what upstream's
// value-initialized `T{}` gives (a default-initialized one is indeterminate in C++).
// coefficients() returns a slice of coef_value or of the caller's scratch in place of a
// raw pointer. clear() keeps `surface`, as upstream's move-aside / reset / move-back does.

use crate::core::options::Bytes;
use crate::data::msx_surface::MsxSurfaceData;
use crate::data::reaction_tokens::{RxExprSpan, RxToken};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum ReactionSolverKind {
    Eul = 0,
    Rk5 = 1,
    Ros2 = 2,
    Bdf2 = 3,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum ReactionCoupling {
    None = 0,
    Full = 1,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum ReactionRateUnits {
    Sec = 0,
    Min = 1,
    Hr = 2,
    Day = 3,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum ReactionAreaUnits {
    Ft2 = 0,
    M2 = 1,
    Cm2 = 2,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum ReactionExprForm {
    None = 0,
    Rate = 1,
    Equil = 2,
    Formula = 3,
}

/// ReactionData::ParameterOverride — a sparse authored per-element PARAMETER override.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ReactionDataParameterOverride {
    pub link: bool,
    pub element: i32,
    pub coefficient: i32,
    pub value: f64,
}

/// ReactionData::WallArea — declared AREA_UNITS.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ReactionDataWallArea {
    pub link: bool,
    pub element: i32,
    pub area: f64,
}

/// ReactionData::LidLaw
#[derive(Debug, Clone, PartialEq)]
pub struct ReactionDataLidLaw {
    pub control: i32,
    pub layer: i32,
    pub species: i32,
    pub form: ReactionExprForm,
    pub expression: Bytes,
    pub compiled: RxExprSpan,
}

impl Default for ReactionDataLidLaw {
    fn default() -> Self {
        ReactionDataLidLaw {
            control: -1,
            layer: 0,
            species: -1,
            form: ReactionExprForm::None,
            expression: Vec::new(),
            compiled: RxExprSpan::default(),
        }
    }
}

/// ReactionData::LidState
#[derive(Debug, Clone, Default)]
pub struct ReactionDataLidState {
    pub node_mass: Vec<Vec<f64>>,
    pub node_outlet: Vec<Vec<f64>>,
    pub node_mobile_old: Vec<f64>,
    pub group_offsets: Vec<i32>,
    /// conventional units, conc*ft3
    pub mass: Vec<f64>,
    pub old_water: Vec<f64>,
    pub surface_conc: Vec<f64>,
    pub drain_conc: Vec<f64>,
    /// [subcatchment*species], conc*ft3
    pub surface_export: Vec<f64>,
    pub surface_carry: Vec<f64>,
    pub subcatchments: Vec<u8>,
}

/// ReactionData::SourcePattern
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ReactionDataSourcePattern {
    pub name: Bytes,
    pub factors: Vec<f64>,
}

/// ReactionData::Source — booster: 1 FLOWPACED, 2 SETPOINT.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ReactionDataSource {
    pub mass: bool,
    pub node: i32,
    pub species: i32,
    pub pattern: i32,
    pub value: f64,
    pub booster: i32,
}

/// ReactionData::SubsurfaceLaw — scope: global/tag/cell; zone: sat/unsat/layer.
#[derive(Debug, Clone, PartialEq)]
pub struct ReactionDataSubsurfaceLaw {
    pub scope: i32,
    pub cell: i32,
    pub zone: i32,
    pub layer: i32,
    pub species: i32,
    pub tag: Bytes,
    pub expression: Bytes,
    pub form: ReactionExprForm,
    pub compiled: RxExprSpan,
}

impl Default for ReactionDataSubsurfaceLaw {
    fn default() -> Self {
        ReactionDataSubsurfaceLaw {
            scope: 0,
            cell: -1,
            zone: 0,
            layer: -1,
            species: -1,
            tag: Vec::new(),
            expression: Vec::new(),
            form: ReactionExprForm::None,
            compiled: RxExprSpan::default(),
        }
    }
}

/// ReactionData::SurfaceLaw
#[derive(Debug, Clone, PartialEq)]
pub struct ReactionDataSurfaceLaw {
    pub subcatchment: i32,
    pub landuse: i32,
    pub species: i32,
    pub form: ReactionExprForm,
    pub expression: Bytes,
    pub compiled: RxExprSpan,
}

impl Default for ReactionDataSurfaceLaw {
    fn default() -> Self {
        ReactionDataSurfaceLaw {
            subcatchment: -1,
            landuse: -1,
            species: -1,
            form: ReactionExprForm::None,
            expression: Vec::new(),
            compiled: RxExprSpan::default(),
        }
    }
}

/// ReactionData::MsxRoutingLedger — the reaction-species routing ledger under LEGACY,
/// per species, in conc-units × ft³.
#[derive(Debug, Clone, Default)]
pub struct ReactionDataMsxRoutingLedger {
    pub armed: bool,
    pub init: Vec<f64>,
    pub inflow: Vec<f64>,
    pub outflow: Vec<f64>,
    pub flooding: Vec<f64>,
    pub seepage: Vec<f64>,
    pub reacted: Vec<f64>,
    pub clamped: Vec<f64>,
    pub final_dry: Vec<f64>,
}

impl ReactionDataMsxRoutingLedger {
    /// Zero and activate the ledger (armed is set before the columns are sized).
    pub fn arm(&mut self, ns: usize) {
        self.armed = true;
        for v in [
            &mut self.init,
            &mut self.inflow,
            &mut self.outflow,
            &mut self.flooding,
            &mut self.seepage,
            &mut self.reacted,
            &mut self.clamped,
            &mut self.final_dry,
        ] {
            v.clear();
            v.resize(ns, 0.0);
        }
    }
}

#[derive(Debug, Clone)]
pub struct ReactionData {
    // ---- [REACTION_OPTIONS]
    /// Default integrator (RK5; see upstream's measured table for why).
    pub solver: ReactionSolverKind,
    pub coupling: ReactionCoupling,
    pub rate_units: ReactionRateUnits,
    pub area_units: ReactionAreaUnits,
    /// WALL densities use a fixed reference contact area
    pub wall_geometry: bool,
    /// 0 ⇒ follow QUALITY_STEP
    pub timestep: f64,
    pub atol: f64,
    pub rtol: f64,
    /// TEMP's value when no heat-transport temperature exists for the element.
    pub default_temp_c: f64,

    // ---- [REACTION_SPECIES]
    /// first MSX index in the registry
    pub registry_base: i32,
    /// cold
    pub species_name: Vec<Bytes>,
    /// 0 = BULK, 1 = WALL
    pub species_is_wall: Vec<u8>,
    /// cold
    pub species_units: Vec<Bytes>,
    /// hot (0 ⇒ global)
    pub species_atol: Vec<f64>,
    /// hot (0 ⇒ global)
    pub species_rtol: Vec<f64>,

    // ---- [REACTION_COEFFICIENTS]
    /// cold
    pub coef_name: Vec<Bytes>,
    /// 1 = PARAMETER (overridable), 0 = CONSTANT
    pub coef_is_param: Vec<u8>,
    /// hot (global values)
    pub coef_value: Vec<f64>,

    pub parameter_overrides: Vec<ReactionDataParameterOverride>,
    pub wall_areas: Vec<ReactionDataWallArea>,
    pub lid_controls: Vec<i32>,
    pub lid_laws: Vec<ReactionDataLidLaw>,
    pub lid: ReactionDataLidState,

    /// actual cumulative injection, concentration*ft3
    pub booster_mass: Vec<f64>,
    /// seconds, independent of RATE_UNITS
    pub pattern_timestep: f64,
    pub source_patterns: Vec<ReactionDataSourcePattern>,
    pub sources: Vec<ReactionDataSource>,

    /// Native result selection. Empty means all species (legacy default).
    pub report_enabled: Vec<u8>,
    pub reported_species: Vec<i32>,

    // ---- [REACTION_TERMS]
    /// cold
    pub term_name: Vec<Bytes>,
    /// cold (compiled in R2)
    pub term_expr_src: Vec<Bytes>,

    pub subsurface_laws: Vec<ReactionDataSubsurfaceLaw>,

    // ---- [REACTION_PIPES] / [REACTION_TANKS] — dense per species
    /// size n_species
    pub pipe_form: Vec<ReactionExprForm>,
    /// size n_species (cold)
    pub pipe_expr_src: Vec<Bytes>,
    pub tank_form: Vec<ReactionExprForm>,
    pub tank_expr_src: Vec<Bytes>,

    pub wet_laws: Vec<ReactionDataSurfaceLaw>,
    pub buildup_laws: Vec<ReactionDataSurfaceLaw>,

    // Wet, non-LID subcatchment pools have their own reactor context.
    pub subcatch_form: Vec<ReactionExprForm>,
    pub subcatch_expr_src: Vec<Bytes>,
    pub subcatch_expr: Vec<RxExprSpan>,

    // ---- [REACTION_QUALITY] GLOBAL initial values (R1)
    /// size n_species
    pub init_global: Vec<f64>,
    /// [REACTION_RAINFALL]: concentration of each species in rainfall. Size n_species.
    pub rain_conc: Vec<f64>,

    // ---- [REACTION_QUALITY] NODE/LINK per-element overrides (E-B)
    /// 0 = NODE row, 1 = LINK row
    pub init_elem_is_link: Vec<u8>,
    /// resolved node/link index
    pub init_elem_idx: Vec<i32>,
    /// 0..n_species()-1
    pub init_elem_species: Vec<i32>,
    pub init_elem_value: Vec<f64>,

    // ---- Compiled bytecode (R2)
    /// hot: one contiguous pool
    pub token_pool: Vec<RxToken>,
    /// per term, in-order evaluation
    pub term_expr: Vec<RxExprSpan>,
    /// per species (len 0 ⇒ none)
    pub pipe_expr: Vec<RxExprSpan>,
    /// per species
    pub tank_expr: Vec<RxExprSpan>,
    /// R2 compile pass succeeded
    pub compiled: bool,
    pub warned_react_failure: bool,
    /// review A2: once-per-run note
    pub warned_negative_clamp: bool,

    // ---- R4: MSX species element state under QUALITY_SOLVER LEGACY
    /// [element * n_species + s]
    pub msx_node_conc: Vec<f64>,
    pub msx_link_conc: Vec<f64>,

    // ---- U2: external MSX loads at the node seam
    /// [node * n_species + m], conc-units × ft3/s
    pub msx_ext_mass_in: Vec<f64>,
    /// species MASS handed to a node by the 2D surface drain and the aquifer bed
    pub msx_coupling_queue: Vec<f64>,
    pub msx_ledger: ReactionDataMsxRoutingLedger,

    /// BW-MSX: surface buildup / washoff / sweeping of MSX species.
    pub surface: MsxSurfaceData,

    /// a reactions component applied
    pub configured: bool,
}

impl Default for ReactionData {
    fn default() -> Self {
        ReactionData {
            solver: ReactionSolverKind::Rk5,
            coupling: ReactionCoupling::None,
            rate_units: ReactionRateUnits::Hr,
            area_units: ReactionAreaUnits::Ft2,
            wall_geometry: false,
            timestep: 0.0,
            atol: 1.0e-6,
            rtol: 1.0e-4,
            default_temp_c: 20.0,
            registry_base: -1,
            species_name: Vec::new(),
            species_is_wall: Vec::new(),
            species_units: Vec::new(),
            species_atol: Vec::new(),
            species_rtol: Vec::new(),
            coef_name: Vec::new(),
            coef_is_param: Vec::new(),
            coef_value: Vec::new(),
            parameter_overrides: Vec::new(),
            wall_areas: Vec::new(),
            lid_controls: Vec::new(),
            lid_laws: Vec::new(),
            lid: ReactionDataLidState::default(),
            booster_mass: Vec::new(),
            pattern_timestep: 3600.0,
            source_patterns: Vec::new(),
            sources: Vec::new(),
            report_enabled: Vec::new(),
            reported_species: Vec::new(),
            term_name: Vec::new(),
            term_expr_src: Vec::new(),
            subsurface_laws: Vec::new(),
            pipe_form: Vec::new(),
            pipe_expr_src: Vec::new(),
            tank_form: Vec::new(),
            tank_expr_src: Vec::new(),
            wet_laws: Vec::new(),
            buildup_laws: Vec::new(),
            subcatch_form: Vec::new(),
            subcatch_expr_src: Vec::new(),
            subcatch_expr: Vec::new(),
            init_global: Vec::new(),
            rain_conc: Vec::new(),
            init_elem_is_link: Vec::new(),
            init_elem_idx: Vec::new(),
            init_elem_species: Vec::new(),
            init_elem_value: Vec::new(),
            token_pool: Vec::new(),
            term_expr: Vec::new(),
            pipe_expr: Vec::new(),
            tank_expr: Vec::new(),
            compiled: false,
            warned_react_failure: false,
            warned_negative_clamp: false,
            msx_node_conc: Vec::new(),
            msx_link_conc: Vec::new(),
            msx_ext_mass_in: Vec::new(),
            msx_coupling_queue: Vec::new(),
            msx_ledger: ReactionDataMsxRoutingLedger::default(),
            surface: MsxSurfaceData::default(),
            configured: false,
        }
    }
}

impl ReactionData {
    /// True on the first lid_controls entry equal to `control`; linear scan.
    pub fn lid_enabled(&self, control: i32) -> bool {
        for &c in &self.lid_controls {
            if c == control {
                return true;
            }
        }
        false
    }

    /// Global coefficients, or (when overrides exist and element >= 0) a copy in
    /// `scratch` with every matching override applied in row order (later rows win).
    pub fn coefficients<'a>(
        &'a self,
        link: bool,
        element: i32,
        scratch: &'a mut Vec<f64>,
    ) -> &'a [f64] {
        if self.parameter_overrides.is_empty() || element < 0 {
            return &self.coef_value;
        }
        scratch.clone_from(&self.coef_value);
        for row in &self.parameter_overrides {
            if row.link == link && row.element == element {
                scratch[row.coefficient as usize] = row.value;
            }
        }
        scratch
    }

    /// True if any global subcatch_form is non-NONE or any sparse wet law names sc.
    pub fn wet_enabled(&self, sc: i32) -> bool {
        for &f in &self.subcatch_form {
            if f != ReactionExprForm::None {
                return true;
            }
        }
        for l in &self.wet_laws {
            if l.subcatchment == sc {
                return true;
            }
        }
        false
    }

    /// True for a non-NONE global form or any wet_laws row.
    pub fn has_subcatch_reactions(&self) -> bool {
        for &form in &self.subcatch_form {
            if form != ReactionExprForm::None {
                return true;
            }
        }
        !self.wet_laws.is_empty()
    }

    pub fn n_species(&self) -> i32 {
        self.species_name.len() as i32
    }

    /// First byte-exact species name, or -1.
    pub fn find_species(&self, n: &[u8]) -> i32 {
        for i in 0..self.species_name.len() {
            if self.species_name[i] == n {
                return i as i32;
            }
        }
        -1
    }

    /// First byte-exact coefficient name, or -1.
    pub fn find_coef(&self, n: &[u8]) -> i32 {
        for i in 0..self.coef_name.len() {
            if self.coef_name[i] == n {
                return i as i32;
            }
        }
        -1
    }

    /// First byte-exact term name, or -1.
    pub fn find_term(&self, n: &[u8]) -> i32 {
        for i in 0..self.term_name.len() {
            if self.term_name[i] == n {
                return i as i32;
            }
        }
        -1
    }

    /// Reset everything except `surface`, which keeps its authored rows.
    pub fn clear(&mut self) {
        let keep = std::mem::take(&mut self.surface);
        *self = ReactionData::default();
        self.surface = keep;
    }
}
