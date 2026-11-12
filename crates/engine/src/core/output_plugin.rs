// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/plugins/DefaultOutputPlugin.cpp and
// include/openswmm/plugin_sdk/SimulationSnapshot.hpp (HydroCouple/openswmm.engine, branch
// swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache
// License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The SWMM 5.x binary .out writer: the header (writeHeader, from prepare()), one record per
// report period (update), the closing records (finalize); and the 1D part of the snapshot
// the engine hands it.
//
// Differences of form, never in result:
//   * SimulationSnapshot keeps the fields DefaultOutputPlugin reads and the engine fills in
//     a 1D build, and, with the `twod` feature, the 2D surface block SWMMEngine::
//     fillSurfaceSnapshot fills (SurfaceSnapshot, 12 November 2026, phase 2a wave 5; the
//     subsurface (gw2d_*) and thermal blocks belong to phase 2b and are left out, gw2d_active
//     stays false). No writer of this build reads the surface block: upstream's reader is the
//     HDF5 2D writer, which the reference build leaves out (OPENSWMM_2D_HDF5_OUTPUT=OFF). The
//     name-table pointers (read by no writer of this build) are left out.
//   * The FILE* is a CFile (core/cprintf.rs); fwrite of an int/float/double is its native
//     (little-endian) bytes; ftell is the count of bytes written to the fresh file.
//   * The managed checkpoint members (checkpoint_fresh_, readPrefix and the rest) serve a
//     checkpoint session, which a CLI run never opens; they are not translated.

use crate::core::context::SimulationContext;
use crate::core::cprintf::CFile;
use crate::core::plugins::{DefaultOutputPlugin, PluginState};
use crate::data::link::LinkType;
use crate::data::pollutant::MassUnits;

/// SubcatchSnapshot::snow_process_count
pub const SNOW_PROCESS_COUNT: usize = 10;
/// SWMM_OUT_SUBCATCH_SNOW_SWE (include/openswmm/engine/openswmm_output.h)
pub const SWMM_OUT_SUBCATCH_SNOW_SWE: i32 = 65536;

#[derive(Debug, Clone, Default)]
pub struct NodeSnapshot {
    pub depth: Vec<f64>,
    pub head: Vec<f64>,
    pub volume: Vec<f64>,
    pub lateral_inflow: Vec<f64>,
    pub total_inflow: Vec<f64>,
    pub overflow: Vec<f64>,
}

#[derive(Debug, Clone, Default)]
pub struct LinkSnapshot {
    pub flow: Vec<f64>,
    pub depth: Vec<f64>,
    pub velocity: Vec<f64>,
    pub volume: Vec<f64>,
    pub capacity: Vec<f64>,
}

#[derive(Debug, Clone, Default)]
pub struct SubcatchSnapshot {
    /// [subcatch * SNOW_PROCESS_COUNT + field]
    pub snow_process: Vec<f64>,
    pub rainfall: Vec<f64>,
    pub snow_depth: Vec<f64>,
    pub evap: Vec<f64>,
    pub infil: Vec<f64>,
    pub runoff: Vec<f64>,
    pub gw_flow: Vec<f64>,
    pub gw_elev: Vec<f64>,
    pub soil_moist: Vec<f64>,
}

/// SnowResults.hpp meshSnowFields.size(): the per-cell mesh snow fields of the snapshot.
#[cfg(feature = "twod")]
pub const MESH_SNOW_FIELDS: usize = 23;

/// The 2D surface block of SimulationSnapshot (its surface_* members, SI-native, in
/// upstream's order; `snow` is surface_snow, `species_*` surface_species_*, and so on).
/// species_names / species_units (borrowed name tables) are left out: they are set only with
/// surface transport rows (phase 2b), and nothing of this build reads them.
#[cfg(feature = "twod")]
#[derive(Debug, Clone)]
pub struct SurfaceSnapshot {
    pub tri_count: i32,
    pub vert_count: i32,
    pub quad_count: i32,
    pub edge_stride: i32,
    pub cell_nv: Vec<u8>,
    pub depth: Vec<f64>,
    pub head: Vec<f64>,
    pub grad_hx: Vec<f64>,
    pub grad_hy: Vec<f64>,
    pub grad_hx_lim: Vec<f64>,
    pub grad_hy_lim: Vec<f64>,
    pub snow: [Vec<f64>; MESH_SNOW_FIELDS],
    pub rainfall: Vec<f64>,
    pub coupling_flux: Vec<f64>,
    pub net_source: Vec<f64>,
    pub infil_rate: Vec<f64>,
    pub infil_cum: Vec<f64>,
    pub rain_cum: Vec<f64>,
    pub species_conc: Vec<f64>,
    pub species_count: i32,
    pub buildup: Vec<f64>,
    pub buildup_count: i32,
    pub gw2d_active: bool,
    pub edge_flux: Vec<f64>,
    pub vert_head: Vec<f64>,
    pub vert_depth: Vec<f64>,
    pub face_vx: Vec<f64>,
    pub face_vy: Vec<f64>,
    pub continuity_err: Vec<f64>,
    pub stat_max_depth: Vec<f64>,
    pub stat_max_velocity: Vec<f64>,
    pub stat_max_cont_err: Vec<f64>,
}

#[cfg(feature = "twod")]
impl Default for SurfaceSnapshot {
    fn default() -> Self {
        SurfaceSnapshot {
            tri_count: 0,
            vert_count: 0,
            quad_count: 0,
            edge_stride: 3,
            cell_nv: Vec::new(),
            depth: Vec::new(),
            head: Vec::new(),
            grad_hx: Vec::new(),
            grad_hy: Vec::new(),
            grad_hx_lim: Vec::new(),
            grad_hy_lim: Vec::new(),
            snow: Default::default(),
            rainfall: Vec::new(),
            coupling_flux: Vec::new(),
            net_source: Vec::new(),
            infil_rate: Vec::new(),
            infil_cum: Vec::new(),
            rain_cum: Vec::new(),
            species_conc: Vec::new(),
            species_count: 0,
            buildup: Vec::new(),
            buildup_count: 0,
            gw2d_active: false,
            edge_flux: Vec::new(),
            vert_head: Vec::new(),
            vert_depth: Vec::new(),
            face_vx: Vec::new(),
            face_vy: Vec::new(),
            continuity_err: Vec::new(),
            stat_max_depth: Vec::new(),
            stat_max_velocity: Vec::new(),
            stat_max_cont_err: Vec::new(),
        }
    }
}

/// The 1D part of SimulationSnapshot, and with `twod` its 2D surface block.
#[derive(Debug, Clone)]
pub struct SimulationSnapshot {
    pub sim_time: f64,
    pub node_count: i32,
    pub link_count: i32,
    pub subcatch_count: i32,
    pub gage_count: i32,
    pub pollut_count: i32,
    pub nodes: NodeSnapshot,
    pub links: LinkSnapshot,
    pub subcatch: SubcatchSnapshot,
    pub sys_temperature: f64,
    pub sys_rainfall: f64,
    pub sys_snow_depth: f64,
    pub sys_evap: f64,
    pub sys_infil: f64,
    pub sys_runoff: f64,
    pub sys_dw_inflow: f64,
    pub sys_gw_inflow: f64,
    pub sys_ii_inflow: f64,
    pub sys_ext_inflow: f64,
    pub sys_flooding: f64,
    pub sys_outflow: f64,
    pub sys_storage: f64,
    pub sys_pet: f64,
    pub has_subcatchments: bool,
    pub node_quality: Vec<f64>,
    pub link_quality: Vec<f64>,
    pub subcatch_quality: Vec<f64>,
    pub flow_units_code: i32,
    /// The surface_* members (SWMMEngine::fillSurfaceSnapshot).
    #[cfg(feature = "twod")]
    pub surface: SurfaceSnapshot,
}

impl Default for SimulationSnapshot {
    fn default() -> Self {
        SimulationSnapshot {
            sim_time: 0.0,
            node_count: 0,
            link_count: 0,
            subcatch_count: 0,
            gage_count: 0,
            pollut_count: 0,
            nodes: NodeSnapshot::default(),
            links: LinkSnapshot::default(),
            subcatch: SubcatchSnapshot::default(),
            sys_temperature: 0.0,
            sys_rainfall: 0.0,
            sys_snow_depth: 0.0,
            sys_evap: 0.0,
            sys_infil: 0.0,
            sys_runoff: 0.0,
            sys_dw_inflow: 0.0,
            sys_gw_inflow: 0.0,
            sys_ii_inflow: 0.0,
            sys_ext_inflow: 0.0,
            sys_flooding: 0.0,
            sys_outflow: 0.0,
            sys_storage: 0.0,
            sys_pet: 0.0,
            has_subcatchments: true,
            node_quality: Vec::new(),
            link_quality: Vec::new(),
            subcatch_quality: Vec::new(),
            flow_units_code: 0,
            #[cfg(feature = "twod")]
            surface: SurfaceSnapshot::default(),
        }
    }
}

const MAGIC_NUMBER: i32 = 516114522;
const VERSION: i32 = 60000; // SWMM 6.0.000
const MAX_SYS_RESULTS: usize = 15;

/// msxUnitCode: the .out unit code of a reported reactions-species column (p is past the
/// pollutants and the reserved rows); 0 for anything that is not MG, UG or COUNT.
fn msx_unit_code(ctx: &SimulationContext, p: i32) -> i32 {
    let rx = &ctx.reactions;
    let mut m = p - ctx.n_pollutants() - if ctx.options.water_age { 1 } else { 0 } - if ctx.options.heat_transport { 1 } else { 0 };
    if m < 0 || m as usize >= rx.reported_species.len() {
        return 0;
    }
    m = rx.reported_species[m as usize];
    let mut u = rx.species_units[m as usize].clone();
    for c in u.iter_mut() {
        *c = c.to_ascii_uppercase();
    }
    if u.len() > 2 && &u[u.len() - 2..] == b"/L" {
        u.truncate(u.len() - 2);
    }
    if u == b"MG" {
        return MassUnits::MgPerL as i32;
    }
    if u == b"UG" {
        return MassUnits::UgPerL as i32;
    }
    if u == b"COUNT" {
        return MassUnits::CountsPerL as i32;
    }
    0
}

fn write_int4(f: &mut CFile, v: i32) {
    f.write_bytes(&v.to_ne_bytes());
}
fn write_real4(f: &mut CFile, v: f32) {
    f.write_bytes(&v.to_ne_bytes());
}
fn write_real8(f: &mut CFile, v: f64) {
    f.write_bytes(&v.to_ne_bytes());
}
/// writeID: the name up to its first NUL (strlen), after its length.
fn write_id(f: &mut CFile, id: &[u8]) {
    let id = crate::core::cstdlib::c_str(id);
    write_int4(f, id.len() as i32);
    f.write_bytes(id);
}

/// `uj < v.size() ? v[uj] : 0.0`
fn at(v: &[f64], i: usize) -> f64 {
    if i < v.len() {
        v[i]
    } else {
        0.0
    }
}

impl DefaultOutputPlugin {
    /// update(snapshot): one report period.
    pub fn update(&mut self, snapshot: &SimulationSnapshot) -> i32 {
        let mut f = match self.out_file_.take() {
            Some(f) => f,
            None => return -1,
        };
        self.set_state(PluginState::Updating);

        write_real8(&mut f, snapshot.sim_time);

        let np = self.n_polluts_;
        for j in 0..snapshot.subcatch_count {
            let uj = j as usize;
            if uj < self.subcatch_rpt_flag_.len() && self.subcatch_rpt_flag_[uj] == 0 {
                continue;
            }
            let s = &snapshot.subcatch;
            write_real4(&mut f, at(&s.rainfall, uj) as f32);
            write_real4(&mut f, at(&s.snow_depth, uj) as f32);
            write_real4(&mut f, at(&s.evap, uj) as f32);
            write_real4(&mut f, at(&s.infil, uj) as f32);
            write_real4(&mut f, at(&s.runoff, uj) as f32);
            write_real4(&mut f, at(&s.gw_flow, uj) as f32);
            write_real4(&mut f, at(&s.gw_elev, uj) as f32);
            write_real4(&mut f, at(&s.soil_moist, uj) as f32);
            for p in 0..np {
                let qi = uj * np as usize + p as usize;
                write_real4(&mut f, at(&snapshot.subcatch_quality, qi) as f32);
            }
            if self.snow_process_ {
                for field in 0..SNOW_PROCESS_COUNT {
                    let i = uj * SNOW_PROCESS_COUNT + field;
                    let v = if i < s.snow_process.len() { s.snow_process[i] } else { f64::NAN };
                    write_real4(&mut f, v as f32);
                }
            }
        }

        for j in 0..snapshot.node_count {
            let uj = j as usize;
            if uj < self.node_rpt_flag_.len() && self.node_rpt_flag_[uj] == 0 {
                continue;
            }
            let n = &snapshot.nodes;
            write_real4(&mut f, at(&n.depth, uj) as f32);
            write_real4(&mut f, at(&n.head, uj) as f32);
            write_real4(&mut f, at(&n.volume, uj) as f32);
            write_real4(&mut f, at(&n.lateral_inflow, uj) as f32);
            write_real4(&mut f, at(&n.total_inflow, uj) as f32);
            write_real4(&mut f, at(&n.overflow, uj) as f32);
            for p in 0..np {
                let qi = uj * np as usize + p as usize;
                write_real4(&mut f, at(&snapshot.node_quality, qi) as f32);
            }
        }

        for j in 0..snapshot.link_count {
            let uj = j as usize;
            if uj < self.link_rpt_flag_.len() && self.link_rpt_flag_[uj] == 0 {
                continue;
            }
            let l = &snapshot.links;
            write_real4(&mut f, at(&l.flow, uj) as f32);
            write_real4(&mut f, at(&l.depth, uj) as f32);
            write_real4(&mut f, at(&l.velocity, uj) as f32);
            write_real4(&mut f, at(&l.volume, uj) as f32);
            write_real4(&mut f, at(&l.capacity, uj) as f32);
            for p in 0..np {
                let qi = uj * np as usize + p as usize;
                write_real4(&mut f, at(&snapshot.link_quality, qi) as f32);
            }
        }

        let mut sys = [0.0f32; MAX_SYS_RESULTS];
        sys[0] = snapshot.sys_temperature as f32;
        sys[1] = snapshot.sys_rainfall as f32;
        sys[2] = snapshot.sys_snow_depth as f32;
        sys[3] = snapshot.sys_infil as f32;
        sys[4] = snapshot.sys_runoff as f32;
        sys[5] = snapshot.sys_dw_inflow as f32;
        sys[6] = snapshot.sys_gw_inflow as f32;
        sys[7] = snapshot.sys_ii_inflow as f32;
        sys[8] = snapshot.sys_ext_inflow as f32;
        sys[9] = sys[4] + sys[5] + sys[6] + sys[7] + sys[8];
        sys[10] = snapshot.sys_flooding as f32;
        sys[11] = snapshot.sys_outflow as f32;
        sys[12] = snapshot.sys_storage as f32;
        sys[13] = snapshot.sys_evap as f32;
        sys[14] = snapshot.sys_pet as f32;
        for v in sys {
            write_real4(&mut f, v);
        }
        f.flush();
        self.out_file_ = Some(f);

        self.count_period();
        0
    }

    /// finalize(ctx): the closing records, then fclose.
    pub fn finalize(&mut self, _ctx: &SimulationContext) -> i32 {
        let mut f = match self.out_file_.take() {
            Some(f) => f,
            None => {
                self.set_state(PluginState::Finalized);
                return 0;
            }
        };
        write_int4(&mut f, self.id_start_pos_ as i32);
        write_int4(&mut f, self.input_start_pos_ as i32);
        write_int4(&mut f, self.output_start_pos_ as i32);
        write_int4(&mut f, self.n_periods());
        write_int4(&mut f, 0); // error code
        write_int4(&mut f, MAGIC_NUMBER);
        drop(f); // fclose
        self.set_state(PluginState::Finalized);
        0
    }

    /// writeHeader(ctx): the header, matching legacy output_open().
    pub(crate) fn write_header(&mut self, ctx: &SimulationContext) {
        let mut f = match self.out_file_.take() {
            Some(f) => f,
            None => return,
        };
        self.subcatch_rpt_flag_ = ctx.subcatches.rpt_flag.clone();
        self.node_rpt_flag_ = ctx.nodes.rpt_flag.clone();
        self.link_rpt_flag_ = ctx.links.rpt_flag.clone();

        self.n_subcatch_ = self.subcatch_rpt_flag_.iter().filter(|&&x| x != 0).count() as i32;
        self.n_nodes_ = self.node_rpt_flag_.iter().filter(|&&x| x != 0).count() as i32;
        self.n_links_ = self.link_rpt_flag_.iter().filter(|&&x| x != 0).count() as i32;
        self.n_polluts_ = if ctx.options.ignore_quality { 0 } else { ctx.n_reported_species() };

        self.snow_process_ = crate::hydrology::native_snow::has_ueb_snow(ctx);
        self.n_subcatch_vars_ = 8 + self.n_polluts_ + if self.snow_process_ { SNOW_PROCESS_COUNT as i32 } else { 0 };
        self.n_node_vars_ = 6 + self.n_polluts_;
        self.n_link_vars_ = 5 + self.n_polluts_;

        write_int4(&mut f, MAGIC_NUMBER);
        write_int4(&mut f, VERSION);
        write_int4(&mut f, self.flow_units_code_);
        write_int4(&mut f, self.n_subcatch_);
        write_int4(&mut f, self.n_nodes_);
        write_int4(&mut f, self.n_links_);
        write_int4(&mut f, self.n_polluts_);

        self.id_start_pos_ = f.tell();
        for j in 0..ctx.n_subcatches() {
            let uj = j as usize;
            if uj < self.subcatch_rpt_flag_.len() && self.subcatch_rpt_flag_[uj] != 0 {
                write_id(&mut f, ctx.subcatch_names.name_of(j));
            }
        }
        for j in 0..ctx.n_nodes() {
            let uj = j as usize;
            if uj < self.node_rpt_flag_.len() && self.node_rpt_flag_[uj] != 0 {
                write_id(&mut f, ctx.node_names.name_of(j));
            }
        }
        for j in 0..ctx.n_links() {
            let uj = j as usize;
            if uj < self.link_rpt_flag_.len() && self.link_rpt_flag_[uj] != 0 {
                write_id(&mut f, ctx.link_names.name_of(j));
            }
        }
        for p in 0..self.n_polluts_ {
            write_id(&mut f, &ctx.reported_species_names[p as usize]);
        }

        for p in 0..self.n_polluts_ {
            let is_reserved = p >= ctx.n_pollutants();
            write_int4(&mut f, if is_reserved { msx_unit_code(ctx, p) } else { ctx.pollutants.units[p as usize] as i32 });
        }

        self.input_start_pos_ = f.tell();

        write_int4(&mut f, 1);
        write_int4(&mut f, 1); // INPUT_AREA code
        for j in 0..ctx.n_subcatches() {
            let uj = j as usize;
            if uj >= self.subcatch_rpt_flag_.len() || self.subcatch_rpt_flag_[uj] == 0 {
                continue;
            }
            write_real4(&mut f, ctx.subcatches.area[uj] as f32);
        }

        write_int4(&mut f, 3);
        write_int4(&mut f, 0);
        write_int4(&mut f, 2);
        write_int4(&mut f, 4);
        for j in 0..ctx.n_nodes() {
            let uj = j as usize;
            if uj >= self.node_rpt_flag_.len() || self.node_rpt_flag_[uj] == 0 {
                continue;
            }
            write_int4(&mut f, ctx.nodes.type_[uj] as i32);
            write_real4(&mut f, (ctx.nodes.invert_elev[uj] * self.ucf_length_) as f32);
            write_real4(&mut f, (ctx.nodes.full_depth[uj] * self.ucf_length_) as f32);
        }

        write_int4(&mut f, 5);
        write_int4(&mut f, 0);
        write_int4(&mut f, 6);
        write_int4(&mut f, 7);
        write_int4(&mut f, 4);
        write_int4(&mut f, 8);
        for j in 0..ctx.n_links() {
            let uj = j as usize;
            if uj >= self.link_rpt_flag_.len() || self.link_rpt_flag_[uj] == 0 {
                continue;
            }
            let lt = ctx.links.type_[uj];
            write_int4(&mut f, lt as i32);
            if lt == LinkType::Pump {
                write_real4(&mut f, 0.0);
                write_real4(&mut f, 0.0);
                write_real4(&mut f, 0.0);
                write_real4(&mut f, 0.0);
            } else {
                let mut off1 = ctx.links.offset1[uj] * self.ucf_length_;
                let mut off2 = ctx.links.offset2[uj] * self.ucf_length_;
                if ctx.links.direction[uj] < 0 {
                    std::mem::swap(&mut off1, &mut off2);
                }
                write_real4(&mut f, off1 as f32);
                write_real4(&mut f, off2 as f32);
                if lt == LinkType::Outlet {
                    write_real4(&mut f, 0.0);
                } else {
                    write_real4(&mut f, (ctx.links.xsect_y_full[uj] * self.ucf_length_) as f32);
                }
                if lt == LinkType::Conduit {
                    let cr = ctx.link_subtypes.conduit_row(j);
                    let len = if cr >= 0 { ctx.link_subtypes.conduits.length[cr as usize] } else { 0.0 };
                    write_real4(&mut f, (len * self.ucf_length_) as f32);
                } else {
                    write_real4(&mut f, 0.0);
                }
            }
        }

        write_int4(&mut f, self.n_subcatch_vars_);
        for v in 0..8 + self.n_polluts_ {
            write_int4(&mut f, v);
        }
        if self.snow_process_ {
            for v in 0..SNOW_PROCESS_COUNT as i32 {
                write_int4(&mut f, SWMM_OUT_SUBCATCH_SNOW_SWE + v);
            }
        }
        write_int4(&mut f, self.n_node_vars_);
        for v in 0..self.n_node_vars_ {
            write_int4(&mut f, v);
        }
        write_int4(&mut f, self.n_link_vars_);
        for v in 0..self.n_link_vars_ {
            write_int4(&mut f, v);
        }
        write_int4(&mut f, MAX_SYS_RESULTS as i32);
        for v in 0..MAX_SYS_RESULTS as i32 {
            write_int4(&mut f, v);
        }

        write_real8(&mut f, ctx.options.start_date);
        write_int4(&mut f, ctx.options.report_step as i32);

        self.output_start_pos_ = f.tell();
        self.out_file_ = Some(f);
    }
}
