// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/LidNodeData.hpp, with the types it needs
// from src/engine/hydrology/RichardsColumn.hpp (Material, Options, Transfer, Report) and
// src/engine/hydrology/Infiltration.hpp (GreenAmptState)
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson (first pass by tools/soa2rs.py, finished by hand). Not an
// OpenSWMM release.
//
// Plain value types for LID nodes: layers, cells, the per-node state and config. The
// three included hydrology types (richards::Material/Options/Report with Transfer, and
// GreenAmptState) are translated here as data only, because no hydrology module exists
// yet; the functions declared beside them in those headers are not. Member initializers
// carry into hand-written Default impls. Fields keep upstream's names, including the
// capitalized ones (Ks, IMDmax, ...), so non_snake_case is allowed where they occur.
// Upstream's two initializer-less structs (LidWaterTransfer, LidPortTransfer) default to
// zero here, as value-initialization (`T{}`) gives in C++.

#![allow(clippy::all)]
use crate::core::options::Bytes;

/// The types this file needs from hydrology/RichardsColumn.hpp (namespace
/// openswmm::richards). All kernel quantities are SI.
pub mod richards {
    use crate::core::options::Bytes;

    #[derive(Debug, Clone, Copy)]
    pub struct Material {
        pub theta_r: f64,
        pub alpha: f64,
        pub n: f64,
        pub l: f64,
        pub specific_storage: f64,
    }

    impl Default for Material {
        fn default() -> Self {
            Material {
                theta_r: 0.0,
                alpha: 0.0,
                n: 0.0,
                l: 0.5,
                specific_storage: 0.0,
            }
        }
    }

    #[derive(Debug, Clone, Copy)]
    pub struct Options {
        pub enabled: bool,
        pub cells_per_layer: i32,
        pub atol: f64,
        pub rtol: f64,
        pub max_step: f64,
    }

    impl Default for Options {
        fn default() -> Self {
            Options {
                enabled: false,
                cells_per_layer: 8,
                atol: 1.0e-7,
                rtol: 1.0e-5,
                max_step: 30.0,
            }
        }
    }

    /// Negative indices: -2 surface supply, -3 evaporation, -4 bottom outflow,
    /// -5 bottom inflow.
    #[derive(Debug, Clone, Copy)]
    pub struct Transfer {
        pub from: i32,
        pub to: i32,
        pub volume: f64,
    }

    impl Default for Transfer {
        fn default() -> Self {
            Transfer {
                from: -1,
                to: -1,
                volume: 0.0,
            }
        }
    }

    #[derive(Debug, Clone)]
    pub struct Report {
        pub ok: bool,
        pub error: Bytes,
        pub accepted: i32,
        pub rejected: i32,
        pub rhs: i32,
        pub newton: i32,
        pub min_step: f64,
        pub bottom_volume: f64,
        pub evaporation: f64,
        pub balance: f64,
        pub bottom_outflow: f64,
        pub bottom_inflow: f64,
        pub inflow: f64,
        pub bottom_out_limited: i32,
        pub bottom_in_limited: i32,
        pub transfers: Vec<Transfer>,
    }

    impl Default for Report {
        fn default() -> Self {
            Report {
                ok: true,
                error: Bytes::new(),
                accepted: 0,
                rejected: 0,
                rhs: 0,
                newton: 0,
                min_step: 0.0,
                bottom_volume: 0.0,
                evaporation: 0.0,
                balance: 0.0,
                bottom_outflow: 0.0,
                bottom_inflow: 0.0,
                inflow: 0.0,
                bottom_out_limited: 0,
                bottom_in_limited: 0,
                transfers: Vec::new(),
            }
        }
    }
}

/// Green-Ampt state, from hydrology/Infiltration.hpp. Every initializer is zero or
/// false, so Default is derived.
#[allow(non_snake_case)]
#[derive(Debug, Clone, Copy, Default)]
pub struct GreenAmptState {
    /// Capillary suction head (ft)
    pub S: f64,
    /// Saturated hydraulic conductivity (ft/sec)
    pub Ks: f64,
    /// Maximum initial moisture deficit (0-1)
    pub IMDmax: f64,
    /// Current moisture deficit (0-1)
    pub IMD: f64,
    /// Cumulative infiltration (ft)
    pub F: f64,
    /// Upper zone saturation volume (ft)
    pub Fu: f64,
    /// Max upper zone saturation (ft)
    pub Fumax: f64,
    /// Upper zone depth (ft)
    pub Lu: f64,
    /// Inter-event timer (sec, counts down)
    pub T: f64,
    /// True when surface is saturated
    pub saturated: bool,
}

// Authored, ordered layers. Thickness/suction are in rain-depth units (in/mm),
// conductivity in in/hr or mm/hr. BOTTOM is a boundary, not a numbered layer.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum LidNodeLayerKind {
    Surface = 0,
    Media = 1,
    Aggregate = 2,
    Bottom = 3,
}

// removal = 0.0 and decay = 0.0 are Rust's defaults, so Default is derived.
#[derive(Debug, Clone, Default)]
pub struct LidLayerTreatment {
    pub pollutant: Bytes,
    /// fraction of outgoing mass
    pub removal: f64,
    /// first-order rate, per day
    pub decay: f64,
    /// optional R= or C=, applied after fixed removal
    pub expression: Bytes,
}

#[derive(Debug, Clone)]
pub struct LidNodeLayer {
    pub kind: LidNodeLayerKind,
    pub params: [f64; 7],
    pub treatment: Vec<LidLayerTreatment>,
    /// explicit inverse-metre parameters
    pub retention: richards::Material,
    /// control options, stored on first layer only
    pub flow: richards::Options,
}

impl Default for LidNodeLayer {
    fn default() -> Self {
        LidNodeLayer {
            kind: LidNodeLayerKind::Media,
            params: [0.0; 7],
            treatment: Vec::new(),
            retention: richards::Material::default(),
            flow: richards::Options::default(),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct LidNodeCell {
    /// feet above invert
    pub bottom: f64,
    pub top: f64,
    pub geometric_volume: f64,
    pub area: f64,
    pub porosity: f64,
    pub theta: f64,
    pub field_capacity: f64,
    pub wilting_point: f64,
    pub conductivity: f64,
    pub conductivity_slope: f64,
    pub suction: f64,
    pub kind: LidNodeLayerKind,
    pub layer: i32,
}

impl Default for LidNodeCell {
    fn default() -> Self {
        LidNodeCell {
            bottom: 0.0,
            top: 0.0,
            geometric_volume: 0.0,
            area: 0.0,
            porosity: 1.0,
            theta: 0.0,
            field_capacity: 0.0,
            wilting_point: 0.0,
            conductivity: 0.0,
            conductivity_slope: 0.0,
            suction: 0.0,
            kind: LidNodeLayerKind::Aggregate,
            layer: 0,
        }
    }
}

/// -1 mobile, -2 lateral, -3 evaporation. Upstream gives no initializers.
#[derive(Debug, Clone, Copy, Default)]
pub struct LidWaterTransfer {
    pub from: i32,
    pub to: i32,
    pub volume: f64,
}

/// positive into held water. Upstream gives no initializers.
#[derive(Debug, Clone, Copy, Default)]
pub struct LidPortTransfer {
    pub link: i32,
    pub cell: i32,
    pub volume: f64,
}

#[derive(Debug, Clone)]
pub struct LidNodeState {
    pub richards: bool,
    /// ft3, full porous stores; surface row is unused
    pub richards_water: Vec<f64>,
    /// ft, negative means suction
    pub richards_pressure: Vec<f64>,
    /// SI diagnostics for last accepted interval
    pub richards_report: richards::Report,
    /// ft3, outside-mesh system loss
    pub richards_bottom_loss: f64,
    /// internal owner clock; R5 joint restart still guarded
    pub richards_completed_time: f64,
    /// accepted hydraulic interval; transient until joint restart qualifies
    pub richards_mesh_start: f64,
    /// accepted interval awaits one hydraulic port finish
    pub richards_mesh_pending_finish: bool,
    /// final native routing trial awaits finish; no replay/reset
    pub richards_mesh_routed: bool,
    /// One modified Green-Ampt history for the media directly below SURFACE.
    /// Other layer/port wetting uses the physical cell and mobile inventories.
    pub surface_infil: GreenAmptState,
    pub infiltration_cell: i32,
    /// last reconciled mobile water-table depth
    pub infiltration_head: f64,

    /// `[cell * pollutant + pollutant]`, retained stores
    pub quality_mass: Vec<f64>,
    pub quality_old_water: Vec<f64>,
    pub quality_transfers: Vec<LidWaterTransfer>,
    pub quality_ports: Vec<LidPortTransfer>,
    /// `[link * pollutant + pollutant]`, NaN means mobile
    pub quality_outlet_conc: Vec<f64>,
    /// canonical J per physical cell; never a dry-temperature display value
    pub thermal_energy: Vec<f64>,
    /// cumulative signed soil-contact receipt, positive into each LID cell
    pub thermal_ground_contact_j: Vec<f64>,
    pub quality_old_mobile: f64,
    /// accepted port balance in the current hydraulic iteration
    pub mobile_delta: f64,
    /// accepted link exchange volumes, reset each hydraulic iteration
    pub port_delta: Vec<f64>,
    /// top to bottom; MEDIA subdivided dynamically
    pub cells: Vec<LidNodeCell>,
    pub held_volume: f64,
    pub evap_volume: f64,
    pub captured_flow: f64,
    pub treated_volume: f64,
}

impl Default for LidNodeState {
    fn default() -> Self {
        LidNodeState {
            richards: false,
            richards_water: Vec::new(),
            richards_pressure: Vec::new(),
            richards_report: richards::Report::default(),
            richards_bottom_loss: 0.0,
            richards_completed_time: 0.0,
            richards_mesh_start: 0.0,
            richards_mesh_pending_finish: false,
            richards_mesh_routed: false,
            surface_infil: GreenAmptState::default(),
            infiltration_cell: -1,
            infiltration_head: 0.0,
            quality_mass: Vec::new(),
            quality_old_water: Vec::new(),
            quality_transfers: Vec::new(),
            quality_ports: Vec::new(),
            quality_outlet_conc: Vec::new(),
            thermal_energy: Vec::new(),
            thermal_ground_contact_j: Vec::new(),
            quality_old_mobile: 0.0,
            mobile_delta: 0.0,
            port_delta: Vec::new(),
            cells: Vec::new(),
            held_volume: 0.0,
            evap_volume: 0.0,
            captured_flow: 0.0,
            treated_volume: 0.0,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum LidNodeBottom {
    Current = 0,
    Sealed = 1,
    MeshAquifer = 2,
}

#[derive(Debug, Clone, Copy)]
pub struct LidNodeConfig {
    pub control: i32,
    /// percent
    pub initial_saturation: f64,
    pub bottom: LidNodeBottom,
}

impl Default for LidNodeConfig {
    fn default() -> Self {
        LidNodeConfig {
            control: -1,
            initial_saturation: 0.0,
            bottom: LidNodeBottom::Current,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct LidNodeOutlet {
    pub link: i32,
    /// one-based, top to bottom
    pub layer: i32,
    pub top: bool,
}

impl Default for LidNodeOutlet {
    fn default() -> Self {
        LidNodeOutlet {
            link: -1,
            layer: 1,
            top: false,
        }
    }
}
