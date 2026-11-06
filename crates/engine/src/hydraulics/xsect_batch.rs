// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/XSectBatch.hpp,
// src/engine/hydraulics/XSectBatch.cpp and the scalar semantics of simd::multiply in
// src/engine/math/SIMD.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit
// 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The cross-section data types (XSectShape, XSectParams), the shape-grouped batch
// geometry (ShapeGroup, XSectGroups) and the batch kernels (upstream namespace
// xsect_batch). The per-element functions this header declares are defined in
// XSection.cpp and live in hydraulics::xsection.
//
// Differences in form, never in result:
//   * Raw `const double*` spans plus a count become slices plus the same `i32` count;
//     outputs are `&mut [f64]`. A negative count runs no iteration, as upstream's loops.
//   * XSectParams' borrowed table pointers are `Option<&'a [f64]>` and
//     `Option<&'a LocateLut>`, so the struct carries the lifetime of the tables it
//     borrows (upstream: "their backing arrays must remain live/unmoved").
//   * ShapeGroup's per-link transect table pointers (`std::vector<const double*>`) are
//     `Vec<Option<Arc<[f64]>>>`: attach_transect_tables copies each attached table out
//     of ctx.transect_tables once, instead of pointing into it, so the groups do not
//     borrow the context the solver mutates. Upstream requires those tables to stay
//     unchanged for the run, so the copies hold the same values.
//   * `mutable` members (the three scratch buffers, the packed bypass-mask views, the
//     mask flag) are RefCell/Cell, so the compute methods and set_bypass_mask keep
//     upstream's const (`&self`) signatures. RefCell is not Sync: the *Team variants run
//     one worker slice per call, and a caller runs the slices in turn; each element's
//     result does not depend on the slicing (upstream's own contract), so the outputs
//     are those of the serial pass.
//   * simd::multiply (SIMD.hpp) is kept scalar: dst[i] = a[i] * b[i]. The reference
//     build compiles its AVX2 path (-mavx2 -mfma), which is _mm256_mul_pd plus a scalar
//     tail; with -ffp-contract=off that is an element-wise IEEE multiply, the same
//     bits as the scalar product.
//   * Only the default build (SWMM_XSECT_FAST_LOOKUP off) is translated: norm_param is
//     y_full, and the fast-lookup branches of the kernels are omitted.
//   * `void build(const SimulationContext&)` is declared upstream but never defined, so
//     it is not translated. `static_cast<int>(x / delta)` in perlink_tabulated_pair is
//     cvt_i32 (x86 INT_MIN on a non-finite x, where upstream then reads out of bounds;
//     Rust panics on that index instead).
//   * XSectGroups' trailing-underscore member names drop the underscore.
#![allow(clippy::too_many_arguments, clippy::needless_range_loop)]

use std::cell::{Cell, RefCell};
use std::sync::Arc;

use crate::core::cmath::{cmin, cvt_i32};
use crate::core::context::SimulationContext;
use crate::hydraulics::transect_types::N_TRANSECT_TBL;
use crate::hydraulics::xsect_lookup::{norm_lookup, norm_x, LocateLut};
use crate::hydraulics::xsect_tables as xt;
use crate::hydraulics::{xsect_kernels, xsection};

// ============================================================================
// Cross-section shape codes (matches legacy enums.h XsectType)
// ============================================================================

#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub enum XSectShape {
    #[default]
    Dummy = 0,
    Circular = 1,
    FilledCircular = 2,
    RectClosed = 3,
    RectOpen = 4,
    Trapezoidal = 5,
    Triangular = 6,
    Parabolic = 7,
    Powerfunc = 8,
    RectTriang = 9,
    RectRound = 10,
    ModBasket = 11,
    HorizEllipse = 12,
    VertEllipse = 13,
    Arch = 14,
    Eggshaped = 15,
    Horseshoe = 16,
    Gothic = 17,
    Catenary = 18,
    Semielliptical = 19,
    Baskethandle = 20,
    Semicircular = 21,
    Irregular = 22,
    Custom = 23,
    ForceMain = 24,
    StreetXsect = 25,
}

impl XSectShape {
    /// The shape with this kernel code, or None for a code upstream's enum does not name
    /// (a `switch` on such a code takes its `default` arm).
    pub fn from_code(t: i32) -> Option<XSectShape> {
        use XSectShape::*;
        Some(match t {
            0 => Dummy,
            1 => Circular,
            2 => FilledCircular,
            3 => RectClosed,
            4 => RectOpen,
            5 => Trapezoidal,
            6 => Triangular,
            7 => Parabolic,
            8 => Powerfunc,
            9 => RectTriang,
            10 => RectRound,
            11 => ModBasket,
            12 => HorizEllipse,
            13 => VertEllipse,
            14 => Arch,
            15 => Eggshaped,
            16 => Horseshoe,
            17 => Gothic,
            18 => Catenary,
            19 => Semielliptical,
            20 => Baskethandle,
            21 => Semicircular,
            22 => Irregular,
            23 => Custom,
            24 => ForceMain,
            25 => StreetXsect,
            _ => return None,
        })
    }
}

// ============================================================================
// Cross-section parameter struct (mirrors legacy TXsect)
// ============================================================================

/// Native cross-section geometry plus borrowed optional normalized transect tables.
/// Lengths in ft, areas in ft2.
#[derive(Debug, Clone, Copy)]
pub struct XSectParams<'a> {
    pub type_: i32,
    pub culvert_code: i32,
    pub transect: i32,

    /// Full depth (ft)
    pub y_full: f64,
    /// Width at widest point (ft)
    pub w_max: f64,
    /// Depth at widest point (ft)
    pub yw_max: f64,
    /// Area when full (ft2)
    pub a_full: f64,
    /// Hydraulic radius when full (ft)
    pub r_full: f64,
    /// Section factor when full (normally ft^(8/3); FORCE_MAIN uses R^0.63)
    pub s_full: f64,
    /// Section factor at max flow (same convention as s_full)
    pub s_max: f64,

    /// Depth of bottom section / fill depth
    pub y_bot: f64,
    /// Area of bottom section
    pub a_bot: f64,
    /// Slope of bottom section / exponent
    pub s_bot: f64,
    /// Radius of bottom section / coefficient
    pub r_bot: f64,

    /// Normalized area vs y/y_full (IRREGULAR / CUSTOM / STREET_XSECT only)
    pub area_tbl: Option<&'a [f64]>,
    /// Normalized hyd-radius vs y/y_full
    pub hrad_tbl: Option<&'a [f64]>,
    /// Normalized width vs y/y_full
    pub width_tbl: Option<&'a [f64]>,
    /// Entry count of the tables above
    pub transect_tbl_size: i32,

    /// Bucket map over `area_tbl` for the getYofA inversion; None whenever area_tbl is.
    pub area_lut: Option<&'a LocateLut>,
}

// By hand: `transect` starts at -1.
impl Default for XSectParams<'_> {
    fn default() -> Self {
        XSectParams {
            type_: 0,
            culvert_code: 0,
            transect: -1,
            y_full: 0.0,
            w_max: 0.0,
            yw_max: 0.0,
            a_full: 0.0,
            r_full: 0.0,
            s_full: 0.0,
            s_max: 0.0,
            y_bot: 0.0,
            a_bot: 0.0,
            s_bot: 0.0,
            r_bot: 0.0,
            area_tbl: None,
            hrad_tbl: None,
            width_tbl: None,
            transect_tbl_size: 0,
            area_lut: None,
        }
    }
}

/// A per-link borrowed table (upstream `const double*`, null = None).
pub type TablePtr = Option<Arc<[f64]>>;

// ============================================================================
// Shape group — contiguous SoA for all links sharing one shape type
// ============================================================================

/// SoA parameter block for all links of one cross-section shape.
#[derive(Debug, Clone, Default)]
pub struct ShapeGroup {
    pub shape: XSectShape,
    pub count: i32,

    /// link_idx[i] = index in SimulationContext
    pub link_idx: Vec<i32>,

    /// Full depth (ft)
    pub y_full: Vec<f64>,
    /// Full area (ft2)
    pub a_full: Vec<f64>,
    /// Hyd. radius at full (ft)
    pub r_full: Vec<f64>,
    /// Section factor at full
    pub s_full: Vec<f64>,
    /// Max width (ft)
    pub w_max: Vec<f64>,

    /// 1.0 / y_full (or 0 if y_full <= 0)
    pub inv_y_full: Vec<f64>,

    pub y_bot: Vec<f64>,
    pub a_bot: Vec<f64>,
    pub s_bot: Vec<f64>,
    pub r_bot: Vec<f64>,

    /// Per-link area table (IRREGULAR / CUSTOM / STREET_XSECT)
    pub area_tables: Vec<TablePtr>,
    /// Per-link hyd-rad table
    pub hrad_tables: Vec<TablePtr>,
    /// Per-link width table
    pub width_tables: Vec<TablePtr>,
    /// Table size (same for all)
    pub transect_tbl_size: i32,

    /// Gather buffer for depths
    pub buf_d: RefCell<Vec<f64>>,
    /// Scatter buffer for results
    pub buf_r: RefCell<Vec<f64>>,
    /// Second scatter buffer (for fused ops)
    pub buf_r2: RefCell<Vec<f64>>,
}

impl ShapeGroup {
    /// Set count and resize link indices, scalar geometry and the three scratch vectors,
    /// keeping existing prefixes. The per-link table vectors are not touched.
    pub fn resize(&mut self, n: i32) {
        self.count = n;
        let un = n as usize;
        self.link_idx.resize(un, 0);
        self.y_full.resize(un, 0.0);
        self.inv_y_full.resize(un, 0.0);
        self.a_full.resize(un, 0.0);
        self.r_full.resize(un, 0.0);
        self.s_full.resize(un, 0.0);
        self.w_max.resize(un, 0.0);
        self.y_bot.resize(un, 0.0);
        self.a_bot.resize(un, 0.0);
        self.s_bot.resize(un, 0.0);
        self.r_bot.resize(un, 0.0);
        // Pre-allocate working buffers for hot loop
        self.buf_d.get_mut().resize(un, 0.0);
        self.buf_r.get_mut().resize(un, 0.0);
        self.buf_r2.get_mut().resize(un, 0.0);
    }
}

// ============================================================================
// XSectGroups — the shape-grouped index over all links
// ============================================================================

/// Shape-grouped cross-section manager for batch computation.
#[derive(Debug, Default)]
pub struct XSectGroups {
    groups: Vec<ShapeGroup>,

    // Per-iteration bypass-mask state (see set_bypass_mask). packed_groups[gi] is a
    // packed view of groups[gi] holding only the active links. packed_count[gi]:
    // K_MASK_FULL_GROUP -> no link in the group is bypassed; 0..count -> number of
    // active links in the packed view.
    packed_groups: RefCell<Vec<ShapeGroup>>,
    packed_count: RefCell<Vec<i32>>,
    mask_active: Cell<bool>,
}

impl XSectGroups {
    const K_MASK_FULL_GROUP: i32 = -1;

    /// Number of non-empty shape groups.
    pub fn num_groups(&self) -> i32 {
        self.groups.len() as i32
    }

    /// One indexed shape group (unchecked upstream; panics here when out of range).
    pub fn group(&self, i: i32) -> &ShapeGroup {
        &self.groups[i as usize]
    }

    /// Gate only (cargo feature `gate`): one group, writable, for the routing-step cross
    /// test's loader. Not part of the translation.
    #[cfg(feature = "gate")]
    pub fn gate_group_mut(&mut self, i: i32) -> &mut ShapeGroup {
        &mut self.groups[i as usize]
    }

    /// The first group with the requested kernel shape, or None.
    pub fn find_group(&self, shape: XSectShape) -> Option<&ShapeGroup> {
        self.groups.iter().find(|g| g.shape == shape)
    }

    /// The original, packed or no group under the current bypass mask.
    fn masked_group<'s>(&'s self, gi: usize, packed: &'s [ShapeGroup], pc: &[i32]) -> Option<&'s ShapeGroup> {
        let mut gp = &self.groups[gi];
        if self.mask_active.get() {
            let na = pc[gi];
            if na == 0 {
                return None;
            }
            if na != Self::K_MASK_FULL_GROUP {
                gp = &packed[gi];
            }
        }
        Some(gp)
    }

    // ========================================================================
    // XSectGroups::build (from XSectParams array)
    // ========================================================================

    /// Replace the groups from per-link geometry: one group per shape code 1..25 present,
    /// in code order, rows in link order. DUMMY and unknown codes are skipped.
    pub fn build(&mut self, params: &[XSectParams], n_links: i32) {
        self.groups.clear();
        // Invalidate the bypass-mask packed mirrors; they are lazily re-sized on the next
        // set_bypass_mask() call against the new group layout.
        self.packed_groups.get_mut().clear();
        self.packed_count.get_mut().clear();
        self.mask_active.set(false);

        // Count links per shape (skip DUMMY=0: non-conduit links that need no geometry)
        const MAX_SHAPES: i32 = 26;
        let mut shape_count = [0i32; MAX_SHAPES as usize];
        for i in 0..n_links {
            let t = params[i as usize].type_;
            if t > 0 && t < MAX_SHAPES {
                shape_count[t as usize] += 1;
            }
        }

        // Create groups for non-empty shapes (DUMMY excluded)
        let mut cursor = [0i32; MAX_SHAPES as usize];
        let mut group_map = [-1i32; MAX_SHAPES as usize];

        for s in 1..MAX_SHAPES {
            // start at 1, skip DUMMY
            if shape_count[s as usize] == 0 {
                continue;
            }
            group_map[s as usize] = self.groups.len() as i32;
            self.groups.push(ShapeGroup::default());
            let g = self.groups.last_mut().unwrap();
            g.shape = XSectShape::from_code(s).unwrap();
            g.resize(shape_count[s as usize]);
        }

        // Fill groups (skip DUMMY links)
        for i in 0..n_links {
            let pi = &params[i as usize];
            let t = pi.type_;
            if t <= 0 || t >= MAX_SHAPES {
                continue;
            }
            let gi = group_map[t as usize];
            if gi < 0 {
                continue;
            }
            let g = &mut self.groups[gi as usize];
            let c = cursor[t as usize];
            cursor[t as usize] += 1;
            let uc = c as usize;
            g.link_idx[uc] = i;
            g.y_full[uc] = pi.y_full;
            g.inv_y_full[uc] = if pi.y_full > 0.0 { 1.0 / pi.y_full } else { 0.0 };
            g.a_full[uc] = pi.a_full;
            g.r_full[uc] = pi.r_full;
            g.s_full[uc] = pi.s_full;
            g.w_max[uc] = pi.w_max;
            g.y_bot[uc] = pi.y_bot;
            g.a_bot[uc] = pi.a_bot;
            g.s_bot[uc] = pi.s_bot;
            g.r_bot[uc] = pi.r_bot;
        }
    }

    /// Attach the profile tables of irregular, custom and street groups from
    /// ctx.transect_tables (by ctx.links.xsect_curve), update their full dimensions and
    /// stable-sort each such group by transect index. Invalid indices keep their
    /// existing (or None) table slots.
    pub fn attach_transect_tables(&mut self, ctx: &SimulationContext) {
        // The packed mirrors must include the per-link table pointers; force a lazy
        // re-init so they pick up the arrays attached below.
        self.packed_groups.get_mut().clear();
        self.packed_count.get_mut().clear();
        self.mask_active.set(false);
        for g in self.groups.iter_mut() {
            if (g.shape != XSectShape::Irregular
                && g.shape != XSectShape::Custom
                && g.shape != XSectShape::StreetXsect)
                || g.count == 0
            {
                continue;
            }

            let uc = g.count as usize;
            g.area_tables.resize(uc, None);
            g.hrad_tables.resize(uc, None);
            g.width_tables.resize(uc, None);
            g.transect_tbl_size = N_TRANSECT_TBL;

            for k in 0..g.count {
                let uk = k as usize;
                let link_j = g.link_idx[uk];
                let uj = link_j as usize;

                let ci = ctx.links.xsect_curve[uj];
                if ci >= 0 && (ci as usize) < ctx.transect_tables.len() {
                    let td = &ctx.transect_tables[ci as usize];
                    g.area_tables[uk] = Some(Arc::from(&td.area_tbl[..]));
                    g.hrad_tables[uk] = Some(Arc::from(&td.hrad_tbl[..]));
                    g.width_tables[uk] = Some(Arc::from(&td.width_tbl[..]));
                    // Update full-depth properties from transect
                    g.y_full[uk] = td.y_full;
                    g.inv_y_full[uk] = if td.y_full > 0.0 { 1.0 / td.y_full } else { 0.0 };
                    g.a_full[uk] = td.a_full;
                    g.r_full[uk] = td.r_full;
                    g.w_max[uk] = td.w_max;
                }
            }

            sort_group_by_transect(g, ctx);
        }
    }

    // ========================================================================
    // Batch compute — results scattered to global link arrays
    // ========================================================================

    /// Grouped-link area (ft2). Ungrouped entries of `areas` are not written.
    pub fn compute_areas(&self, depths: &[f64], areas: &mut [f64], _n_links: i32) {
        for g in &self.groups {
            if g.count == 0 {
                continue;
            }
            let mut local_d = g.buf_d.borrow_mut();
            let mut local_a = g.buf_r.borrow_mut();
            gather_depths_all(g, depths, &mut local_d);
            apply_area_kernel_all(g, &local_d, &mut local_a);
            scatter_results_all(g, &local_a, areas);
        }
    }

    /// Grouped-link hydraulic radius (ft).
    pub fn compute_hyd_rad(&self, depths: &[f64], hydrad: &mut [f64], _n_links: i32) {
        for g in &self.groups {
            if g.count == 0 {
                continue;
            }
            let mut local_d = g.buf_d.borrow_mut();
            let mut local_r = g.buf_r.borrow_mut();
            gather_depths_all(g, depths, &mut local_d);
            apply_hydrad_kernel_all(g, &local_d, &mut local_r);
            scatter_results_all(g, &local_r, hydrad);
        }
    }

    /// Grouped-link area and hydraulic radius from one gather per group.
    pub fn compute_area_and_hyd_rad(
        &self,
        depths: &[f64],
        areas: &mut [f64],
        hydrad: &mut [f64],
        _n_links: i32,
    ) {
        for g in &self.groups {
            if g.count == 0 {
                continue;
            }
            let mut local_d = g.buf_d.borrow_mut();
            let mut local_a = g.buf_r.borrow_mut(); // area results
            let mut local_h = g.buf_r2.borrow_mut(); // hydrad results
            gather_depths_all(g, depths, &mut local_d); // single gather

            apply_area_hydrad_kernel_all(g, &local_d, &mut local_a, &mut local_h);

            // Two scatters (same group, different output arrays)
            scatter_results_all(g, &local_a, areas);
            scatter_results_all(g, &local_h, hydrad);
        }
    }

    /// Grouped-link top width (ft).
    pub fn compute_widths(&self, depths: &[f64], widths: &mut [f64], _n_links: i32) {
        for g in &self.groups {
            if g.count == 0 {
                continue;
            }
            let mut local_d = g.buf_d.borrow_mut();
            let mut local_w = g.buf_r.borrow_mut();
            gather_depths_all(g, depths, &mut local_d);
            apply_width_kernel_all(g, &local_d, &mut local_w);

            scatter_results_all(g, &local_w, widths);
        }
    }

    /// Fused d1 -> (a1, hrad1), d2 -> a2, dm -> (am, hrad_mid) on one worker's slice of
    /// every group; honors the bypass mask.
    pub fn compute_area_hyd_rad_triple_team(
        &self,
        d1: &[f64],
        d2: &[f64],
        dm: &[f64],
        a1: &mut [f64],
        a2: &mut [f64],
        am: &mut [f64],
        hrad1: &mut [f64],
        hrad_mid: &mut [f64],
        _n_links: i32,
        tid: i32,
        mut nthreads: i32,
    ) {
        if nthreads < 1 {
            nthreads = 1;
        }
        let packed = self.packed_groups.borrow();
        let pc = self.packed_count.borrow();
        for gi in 0..self.groups.len() {
            let g = match self.masked_group(gi, &packed, &pc) {
                Some(g) if g.count != 0 => g,
                _ => continue,
            };

            let (lo, hi) = team_slice(g.count, tid, nthreads);
            if lo >= hi {
                continue;
            }
            let n = hi - lo;

            let mut ld = g.buf_d.borrow_mut();
            let mut la = g.buf_r.borrow_mut();
            let mut lh = g.buf_r2.borrow_mut();

            // d1 → (a1, hrad1)
            gather_depths(g, d1, &mut ld, lo, n);
            apply_area_hydrad_kernel(g, &ld, &mut la, &mut lh, lo, n);
            scatter_results(g, &la, a1, lo, n);
            scatter_results(g, &lh, hrad1, lo, n);

            // d2 → a2
            gather_depths(g, d2, &mut ld, lo, n);
            apply_area_kernel(g, &ld, &mut la, lo, n);
            scatter_results(g, &la, a2, lo, n);

            // dm → (am, hrad_mid)
            gather_depths(g, dm, &mut ld, lo, n);
            apply_area_hydrad_kernel(g, &ld, &mut la, &mut lh, lo, n);
            scatter_results(g, &la, am, lo, n);
            scatter_results(g, &lh, hrad_mid, lo, n);
        }
    }

    /// Serial form of compute_area_hyd_rad_triple_team (worker 0 of 1).
    pub fn compute_area_hyd_rad_triple(
        &self,
        d1: &[f64],
        d2: &[f64],
        dm: &[f64],
        a1: &mut [f64],
        a2: &mut [f64],
        am: &mut [f64],
        hrad1: &mut [f64],
        hrad_mid: &mut [f64],
        n_links: i32,
    ) {
        self.compute_area_hyd_rad_triple_team(d1, d2, dm, a1, a2, am, hrad1, hrad_mid, n_links, 0, 1);
    }

    /// Restrict the triple kernels to non-bypassed links (1 = bypassed) for the current
    /// Picard iteration; None clears the mask.
    pub fn set_bypass_mask(&self, bypassed_by_link: Option<&[u8]>) {
        self.mask_active.set(bypassed_by_link.is_some());
        let mask = match bypassed_by_link {
            Some(m) => m,
            None => return,
        };

        let mut packed = self.packed_groups.borrow_mut();
        let mut pc = self.packed_count.borrow_mut();

        // Lazy (re)build of the packed mirrors, once per group layout.
        if packed.len() != self.groups.len() {
            packed.clear();
            packed.resize_with(self.groups.len(), ShapeGroup::default);
            pc.clear();
            pc.resize(self.groups.len(), Self::K_MASK_FULL_GROUP);
            for gi in 0..self.groups.len() {
                let g = &self.groups[gi];
                let p = &mut packed[gi];
                p.shape = g.shape;
                p.resize(g.count);
                if !g.area_tables.is_empty() {
                    let uc = g.count as usize;
                    p.area_tables.resize(uc, None);
                    p.hrad_tables.resize(uc, None);
                    p.width_tables.resize(uc, None);
                    p.transect_tbl_size = g.transect_tbl_size;
                }
            }
        }

        for gi in 0..self.groups.len() {
            let g = &self.groups[gi];
            let p = &mut packed[gi];
            if g.count == 0 {
                pc[gi] = Self::K_MASK_FULL_GROUP;
                continue;
            }

            // Cheap pre-count so the common all-active case pays one byte-scan.
            let mut nb = 0i32;
            for k in 0..g.count {
                nb += if mask[g.link_idx[k as usize] as usize] != 0 {
                    1
                } else {
                    0
                };
            }
            if nb == 0 {
                pc[gi] = Self::K_MASK_FULL_GROUP;
                continue;
            }
            if nb == g.count {
                pc[gi] = 0;
                continue;
            }

            let has_tbl = !g.area_tables.is_empty();
            let mut n = 0i32;
            for k in 0..g.count {
                let uk = k as usize;
                let li = g.link_idx[uk];
                if mask[li as usize] != 0 {
                    continue;
                }
                let un = n as usize;
                p.link_idx[un] = li;
                p.y_full[un] = g.y_full[uk];
                p.inv_y_full[un] = g.inv_y_full[uk];
                p.a_full[un] = g.a_full[uk];
                p.r_full[un] = g.r_full[uk];
                p.s_full[un] = g.s_full[uk];
                p.w_max[un] = g.w_max[uk];
                p.y_bot[un] = g.y_bot[uk];
                p.a_bot[un] = g.a_bot[uk];
                p.s_bot[un] = g.s_bot[uk];
                p.r_bot[un] = g.r_bot[uk];
                if has_tbl {
                    p.area_tables[un] = g.area_tables[uk].clone();
                    p.hrad_tables[un] = g.hrad_tables[uk].clone();
                    p.width_tables[un] = g.width_tables[uk].clone();
                }
                n += 1;
            }
            p.count = n; // the kernels and gather/scatter read the view's count
            pc[gi] = n;
        }
    }

    /// Fused d1 -> w1, d2 -> w2, dm -> wm on one worker's slice; honors the bypass mask.
    pub fn compute_widths_triple_team(
        &self,
        d1: &[f64],
        d2: &[f64],
        dm: &[f64],
        w1: &mut [f64],
        w2: &mut [f64],
        wm: &mut [f64],
        _n_links: i32,
        tid: i32,
        mut nthreads: i32,
    ) {
        if nthreads < 1 {
            nthreads = 1;
        }
        let packed = self.packed_groups.borrow();
        let pc = self.packed_count.borrow();
        for gi in 0..self.groups.len() {
            let g = match self.masked_group(gi, &packed, &pc) {
                Some(g) if g.count != 0 => g,
                _ => continue,
            };

            let (lo, hi) = team_slice(g.count, tid, nthreads);
            if lo >= hi {
                continue;
            }
            let n = hi - lo;

            let mut ld = g.buf_d.borrow_mut();
            let mut lw = g.buf_r.borrow_mut();

            gather_depths(g, d1, &mut ld, lo, n);
            apply_width_kernel(g, &ld, &mut lw, lo, n);
            scatter_results(g, &lw, w1, lo, n);

            gather_depths(g, d2, &mut ld, lo, n);
            apply_width_kernel(g, &ld, &mut lw, lo, n);
            scatter_results(g, &lw, w2, lo, n);

            gather_depths(g, dm, &mut ld, lo, n);
            apply_width_kernel(g, &ld, &mut lw, lo, n);
            scatter_results(g, &lw, wm, lo, n);
        }
    }

    /// Serial form of compute_widths_triple_team (worker 0 of 1).
    pub fn compute_widths_triple(
        &self,
        d1: &[f64],
        d2: &[f64],
        dm: &[f64],
        w1: &mut [f64],
        w2: &mut [f64],
        wm: &mut [f64],
        n_links: i32,
    ) {
        self.compute_widths_triple_team(d1, d2, dm, w1, w2, wm, n_links, 0, 1);
    }

    /// Section factor per grouped link through the scalar params_at fallback (which omits
    /// s_max and the table pointers).
    pub fn compute_section_factors(&self, areas: &[f64], sfact: &mut [f64], _n_links: i32) {
        for g in &self.groups {
            for k in 0..g.count {
                let uk = k as usize;
                let li = g.link_idx[uk] as usize;
                sfact[li] = xsection::get_s_of_a(&params_at(g, k), areas[li]);
            }
        }
    }

    /// Depth per grouped link through the scalar params_at fallback.
    pub fn compute_depths_from_area(&self, areas: &[f64], depths: &mut [f64], _n_links: i32) {
        for g in &self.groups {
            for k in 0..g.count {
                let uk = k as usize;
                let li = g.link_idx[uk] as usize;
                depths[li] = xsection::get_y_of_a(&params_at(g, k), areas[li]);
            }
        }
    }
}

/// Stable-sort a profile-backed group by model transect index, permuting every column
/// whose length equals count (scratch buffers untouched).
fn sort_group_by_transect(g: &mut ShapeGroup, ctx: &SimulationContext) {
    let count = g.count;
    if count < 2 || g.area_tables.is_empty() {
        return;
    }

    let mut order: Vec<i32> = (0..count).collect();
    {
        let link_idx = &g.link_idx;
        let curve_of = |k: i32| -> i32 {
            let uj = link_idx[k as usize] as usize;
            ctx.links.xsect_curve[uj]
        };
        // std::stable_sort with `curve_of(lhs) < curve_of(rhs)`: Rust's sort_by_key is
        // stable, so the permutation is the same.
        order.sort_by_key(|&k| curve_of(k));
    }

    let mut already_sorted = true;
    for k in 0..count {
        if order[k as usize] != k {
            already_sorted = false;
            break;
        }
    }
    if already_sorted {
        return;
    }

    fn permute<T: Clone>(order: &[i32], count: i32, vec: &mut Vec<T>) {
        if vec.len() as i32 != count {
            return;
        }
        let mut tmp: Vec<T> = Vec::with_capacity(vec.len());
        for k in 0..count {
            tmp.push(vec[order[k as usize] as usize].clone());
        }
        std::mem::swap(vec, &mut tmp);
    }
    permute(&order, count, &mut g.link_idx);
    permute(&order, count, &mut g.y_full);
    permute(&order, count, &mut g.inv_y_full);
    permute(&order, count, &mut g.a_full);
    permute(&order, count, &mut g.r_full);
    permute(&order, count, &mut g.s_full);
    permute(&order, count, &mut g.w_max);
    permute(&order, count, &mut g.y_bot);
    permute(&order, count, &mut g.a_bot);
    permute(&order, count, &mut g.s_bot);
    permute(&order, count, &mut g.r_bot);
    permute(&order, count, &mut g.area_tables);
    permute(&order, count, &mut g.hrad_tables);
    permute(&order, count, &mut g.width_tables);
}

// ============================================================================
// Batch kernels (upstream namespace xsect_batch)
// ============================================================================

/// SIMD.hpp simd::multiply, scalar form: dst[i] = a[i] * b[i] over n doubles. The
/// reference build's AVX2 path (_mm256_mul_pd + scalar tail, -ffp-contract=off) is an
/// element-wise IEEE multiply, so this gives the same bits.
#[inline]
fn simd_multiply(a: &[f64], b: &[f64], dst: &mut [f64], n: usize) {
    for i in 0..n {
        dst[i] = a[i] * b[i];
    }
}

/// a_full * A_Circ(depth/y_full).
pub fn area_circular(depth: &[f64], nrm: &[f64], a_full: &[f64], area: &mut [f64], count: i32) {
    let table: &[f64] = &xt::A_Circ;
    let n_items = xt::N_A_Circ;
    for k in 0..count {
        let k = k as usize;
        area[k] = a_full[k] * norm_lookup(depth[k], nrm[k], table, n_items);
    }
}

/// depth * w_max through simd::multiply.
pub fn area_rect(depth: &[f64], w_max: &[f64], area: &mut [f64], count: i32) {
    simd_multiply(depth, w_max, area, count as usize);
}

pub fn area_trapezoidal(depth: &[f64], y_bot: &[f64], s_bot: &[f64], area: &mut [f64], count: i32) {
    for k in 0..count {
        let k = k as usize;
        area[k] = xsect_kernels::shape::trapez_a_of_y(depth[k], y_bot[k], s_bot[k]);
    }
}

pub fn area_triangular(depth: &[f64], s_bot: &[f64], area: &mut [f64], count: i32) {
    for k in 0..count {
        let k = k as usize;
        area[k] = xsect_kernels::shape::triang_a_of_y(depth[k], s_bot[k]);
    }
}

pub fn area_parabolic(depth: &[f64], r_bot: &[f64], area: &mut [f64], count: i32) {
    for k in 0..count {
        let k = k as usize;
        area[k] = xsect_kernels::shape::parab_a_of_y(depth[k], r_bot[k]);
    }
}

pub fn area_powerfunc(depth: &[f64], s_bot: &[f64], r_bot: &[f64], area: &mut [f64], count: i32) {
    for k in 0..count {
        let k = k as usize;
        area[k] = xsect_kernels::shape::powerfunc_a_of_y(depth[k], s_bot[k], r_bot[k]);
    }
}

/// a_full * table(depth/y_full) for one shared table.
pub fn area_tabulated(
    depth: &[f64],
    nrm: &[f64],
    a_full: &[f64],
    table: &[f64],
    table_size: i32,
    area: &mut [f64],
    count: i32,
) {
    for k in 0..count {
        let k = k as usize;
        area[k] = a_full[k] * norm_lookup(depth[k], nrm[k], table, table_size);
    }
}

/// a_full * invLookup(depth/y_full) for one shared depth-versus-area table.
pub fn area_inv_tabulated(
    depth: &[f64],
    nrm: &[f64],
    a_full: &[f64],
    table: &[f64],
    table_size: i32,
    lut: Option<&LocateLut>,
    area: &mut [f64],
    count: i32,
) {
    for k in 0..count {
        let k = k as usize;
        area[k] = a_full[k] * xsection::inv_lookup(norm_x(depth[k], nrm[k]), table, table_size, lut);
    }
}

/// Fused per-link area + hyd-rad lookup on one shared normalized depth grid;
/// bit-identical to two perlink_tabulated passes.
pub fn perlink_tabulated_pair<T: AsRef<[f64]>>(
    depth: &[f64],
    nrm: &[f64],
    scale_a: &[f64],
    scale_b: &[f64],
    tables_a: &[Option<T>],
    tables_b: &[Option<T>],
    table_size: i32,
    out_a: &mut [f64],
    out_b: &mut [f64],
    count: i32,
) {
    let n = table_size;
    let delta = 1.0 / ((n - 1) as f64);
    let dd = delta * delta;

    for k in 0..count {
        let k = k as usize;
        let ta = tables_a[k].as_ref().map(|t| t.as_ref());
        let tb = tables_b[k].as_ref().map(|t| t.as_ref());
        let (ta, tb) = match (ta, tb) {
            (Some(ta), Some(tb)) => (ta, tb),
            _ => {
                // degenerate link — no fusion to do
                out_a[k] = match ta {
                    Some(t) => scale_a[k] * norm_lookup(depth[k], nrm[k], t, n),
                    None => 0.0,
                };
                out_b[k] = match tb {
                    Some(t) => scale_b[k] * norm_lookup(depth[k], nrm[k], t, n),
                    None => 0.0,
                };
                continue;
            }
        };

        let x = norm_x(depth[k], nrm[k]);
        let i = cvt_i32(x / delta);
        if i >= n - 1 {
            // lookup_exact's early-out (no clamp)
            out_a[k] = scale_a[k] * ta[(n - 1) as usize];
            out_b[k] = scale_b[k] * tb[(n - 1) as usize];
            continue;
        }

        let iu = i as usize;
        let x0 = (i as f64) * delta;
        let mut a = ta[iu] + (x - x0) * (ta[iu + 1] - ta[iu]) / delta;
        let mut b = tb[iu] + (x - x0) * (tb[iu + 1] - tb[iu]) / delta;
        if i < 2 {
            // quadratic refinement (legacy lookup)
            let x1 = ((i as f64) + 1.0) * delta;
            let a2 = a + (x - x0) * (x - x1) / dd * (ta[iu] / 2.0 - ta[iu + 1] + ta[iu + 2] / 2.0);
            let b2 = b + (x - x0) * (x - x1) / dd * (tb[iu] / 2.0 - tb[iu + 1] + tb[iu + 2] / 2.0);
            if a2 > 0.0 {
                a = a2;
            }
            if b2 > 0.0 {
                b = b2;
            }
        }
        if a < 0.0 {
            a = 0.0;
        }
        if b < 0.0 {
            b = 0.0;
        }
        out_a[k] = scale_a[k] * a;
        out_b[k] = scale_b[k] * b;
    }
}

/// Per-link tabulated lookup; a None table writes 0.
pub fn perlink_tabulated<T: AsRef<[f64]>>(
    depth: &[f64],
    nrm: &[f64],
    scale: &[f64],
    tables: &[Option<T>],
    table_size: i32,
    result: &mut [f64],
    count: i32,
) {
    for k in 0..count {
        let k = k as usize;
        match tables[k].as_ref() {
            None => {
                result[k] = 0.0;
                continue;
            }
            Some(t) => {
                result[k] = scale[k] * norm_lookup(depth[k], nrm[k], t.as_ref(), table_size);
            }
        }
    }
}

// ---- Hydraulic radius ----

/// r_full * R_Circ(depth/y_full); R_Circ[0] (0.01) at zero depth, as legacy.
pub fn hydrad_circular(depth: &[f64], nrm: &[f64], r_full: &[f64], hydrad: &mut [f64], count: i32) {
    let table: &[f64] = &xt::R_Circ;
    let n_items = xt::N_R_Circ;
    for k in 0..count {
        let k = k as usize;
        hydrad[k] = r_full[k] * norm_lookup(depth[k], nrm[k], table, n_items);
    }
}

/// Fused circular area + hydraulic radius sharing one index; bit-identical to two
/// lookup_exact calls.
pub fn area_hydrad_circular(
    depth: &[f64],
    nrm: &[f64],
    a_full: &[f64],
    r_full: &[f64],
    area: &mut [f64],
    hydrad: &mut [f64],
    count: i32,
) {
    let aa: &[f64] = &xt::A_Circ;
    let rr: &[f64] = &xt::R_Circ;
    let n = xt::N_A_Circ; // == N_R_Circ
    let delta = 1.0 / ((n - 1) as f64);
    let dd = delta * delta;

    for k in 0..count {
        let k = k as usize;
        let yn = norm_x(depth[k], nrm[k]); // (nrm>0 && d>0) ? d/nrm : 0
        let i = cvt_i32(yn / delta);
        if i >= n - 1 {
            area[k] = a_full[k] * aa[(n - 1) as usize];
            hydrad[k] = r_full[k] * rr[(n - 1) as usize];
            continue;
        }
        let iu = i as usize;
        let x0 = (i as f64) * delta;
        let mut a = aa[iu] + (yn - x0) * (aa[iu + 1] - aa[iu]) / delta;
        let mut r = rr[iu] + (yn - x0) * (rr[iu + 1] - rr[iu]) / delta;
        if i < 2 {
            // quadratic refinement (legacy lookup)
            let x1 = ((i as f64) + 1.0) * delta;
            let a2 = a + (yn - x0) * (yn - x1) / dd * (aa[iu] / 2.0 - aa[iu + 1] + aa[iu + 2] / 2.0);
            let r2 = r + (yn - x0) * (yn - x1) / dd * (rr[iu] / 2.0 - rr[iu + 1] + rr[iu + 2] / 2.0);
            if a2 > 0.0 {
                a = a2;
            }
            if r2 > 0.0 {
                r = r2;
            }
        }
        if a < 0.0 {
            a = 0.0;
        }
        if r < 0.0 {
            r = 0.0;
        }
        area[k] = a_full[k] * a;
        hydrad[k] = r_full[k] * r;
    }
}

pub fn hydrad_trapezoidal(
    depth: &[f64],
    y_bot: &[f64],
    s_bot: &[f64],
    r_bot: &[f64],
    hydrad: &mut [f64],
    count: i32,
) {
    for k in 0..count {
        let k = k as usize;
        hydrad[k] = xsect_kernels::shape::trapez_r_of_y(depth[k], y_bot[k], s_bot[k], r_bot[k]);
    }
}

pub fn hydrad_triangular(depth: &[f64], s_bot: &[f64], r_bot: &[f64], hydrad: &mut [f64], count: i32) {
    for k in 0..count {
        let k = k as usize;
        hydrad[k] = xsect_kernels::shape::triang_r_of_y(depth[k], s_bot[k], r_bot[k]);
    }
}

/// R = (w*d) / (w + 2*d). Upstream keeps this kernel as it was; it has no callers.
pub fn hydrad_rect(depth: &[f64], w_max: &[f64], hydrad: &mut [f64], count: i32) {
    for k in 0..count {
        let k = k as usize;
        let d = depth[k];
        let w = w_max[k];
        hydrad[k] = (w * d) / (w + 2.0 * d);
    }
}

/// RECT_CLOSED radius with the near-full crown correction (not declared in upstream's
/// header; external linkage in XSectBatch.cpp).
pub fn hydrad_rect_closed(depth: &[f64], w_max: &[f64], a_full: &[f64], hydrad: &mut [f64], count: i32) {
    for k in 0..count {
        let k = k as usize;
        hydrad[k] = xsect_kernels::shape::rect_closed_r_of_a(
            xsect_kernels::shape::rect_a_of_y(depth[k], w_max[k]),
            w_max[k],
            a_full[k],
        );
    }
}

/// RECT_OPEN radius honoring the s_bot "sides removed" term (not declared in upstream's
/// header; external linkage in XSectBatch.cpp).
pub fn hydrad_rect_open(depth: &[f64], w_max: &[f64], s_bot: &[f64], hydrad: &mut [f64], count: i32) {
    for k in 0..count {
        let k = k as usize;
        hydrad[k] = xsect_kernels::shape::rect_open_r_of_a(
            xsect_kernels::shape::rect_a_of_y(depth[k], w_max[k]),
            w_max[k],
            s_bot[k],
        );
    }
}

/// RECT_CLOSED top width: 0 at the exact crown (not declared in upstream's header;
/// external linkage in XSectBatch.cpp).
pub fn width_rect_closed(depth: &[f64], nrm: &[f64], w_max: &[f64], width: &mut [f64], count: i32) {
    for k in 0..count {
        let k = k as usize;
        width[k] = xsect_kernels::shape::rect_closed_w_of_y(norm_x(depth[k], nrm[k]), w_max[k]);
    }
}

/// r_full * table(depth/y_full); the table's first value at zero depth.
pub fn hydrad_tabulated(
    depth: &[f64],
    nrm: &[f64],
    r_full: &[f64],
    table: &[f64],
    table_size: i32,
    hydrad: &mut [f64],
    count: i32,
) {
    for k in 0..count {
        let k = k as usize;
        hydrad[k] = r_full[k] * norm_lookup(depth[k], nrm[k], table, table_size);
    }
}

// ---- Top width ----

pub fn width_circular(depth: &[f64], nrm: &[f64], w_max: &[f64], width: &mut [f64], count: i32) {
    let table: &[f64] = &xt::W_Circ;
    let n_items = xt::N_W_Circ;
    for k in 0..count {
        let k = k as usize;
        width[k] = w_max[k] * norm_lookup(depth[k], nrm[k], table, n_items);
    }
}

pub fn width_trapezoidal(depth: &[f64], y_bot: &[f64], s_bot: &[f64], width: &mut [f64], count: i32) {
    for k in 0..count {
        let k = k as usize;
        width[k] = xsect_kernels::shape::trapez_w_of_y(depth[k], y_bot[k], s_bot[k]);
    }
}

pub fn width_triangular(depth: &[f64], s_bot: &[f64], width: &mut [f64], count: i32) {
    for k in 0..count {
        let k = k as usize;
        width[k] = xsect_kernels::shape::triang_w_of_y(depth[k], s_bot[k]);
    }
}

pub fn width_rect(w_max: &[f64], width: &mut [f64], count: i32) {
    for k in 0..count {
        let k = k as usize;
        width[k] = w_max[k];
    }
}

pub fn width_tabulated(
    depth: &[f64],
    nrm: &[f64],
    w_max: &[f64],
    table: &[f64],
    table_size: i32,
    width: &mut [f64],
    count: i32,
) {
    for k in 0..count {
        let k = k as usize;
        width[k] = w_max[k] * norm_lookup(depth[k], nrm[k], table, table_size);
    }
}

// ============================================================================
// Helpers: gather depths for a group, compute via kernel, scatter results
// ============================================================================

/// The group normalization array of the default build (y_full).
#[inline]
fn norm_param(g: &ShapeGroup) -> &[f64] {
    &g.y_full
}

fn gather_depths(g: &ShapeGroup, global_depths: &[f64], local_depths: &mut [f64], lo: i32, n: i32) {
    let hi = lo + n;
    for k in lo..hi {
        let k = k as usize;
        local_depths[k] = global_depths[g.link_idx[k] as usize];
    }
}

fn gather_depths_all(g: &ShapeGroup, global_depths: &[f64], local_depths: &mut [f64]) {
    gather_depths(g, global_depths, local_depths, 0, g.count);
}

fn scatter_results(g: &ShapeGroup, local_results: &[f64], global_results: &mut [f64], lo: i32, n: i32) {
    let hi = lo + n;
    for k in lo..hi {
        let k = k as usize;
        global_results[g.link_idx[k] as usize] = local_results[k];
    }
}

fn scatter_results_all(g: &ShapeGroup, local_results: &[f64], global_results: &mut [f64]) {
    scatter_results(g, local_results, global_results, 0, g.count);
}

/// The stored scalar subset of one group element (no s_max, yw_max, tables or LUT).
#[inline]
fn params_at(g: &ShapeGroup, k: i32) -> XSectParams<'static> {
    let uk = k as usize;
    XSectParams {
        type_: g.shape as i32,
        y_full: g.y_full[uk],
        a_full: g.a_full[uk],
        r_full: g.r_full[uk],
        s_full: g.s_full[uk],
        w_max: g.w_max[uk],
        y_bot: g.y_bot[uk],
        a_bot: g.a_bot[uk],
        s_bot: g.s_bot[uk],
        r_bot: g.r_bot[uk],
        ..XSectParams::default()
    }
}

/// A built-in table and its size (upstream TableRef; None is {nullptr, 0}).
type TableRef = Option<(&'static [f64], i32)>;

fn area_table_for(shape: XSectShape) -> TableRef {
    match shape {
        XSectShape::Eggshaped => Some((&xt::A_Egg, xt::N_A_Egg)),
        XSectShape::Horseshoe => Some((&xt::A_Horseshoe, xt::N_A_Horseshoe)),
        XSectShape::Baskethandle => Some((&xt::A_Baskethandle, xt::N_A_Baskethandle)),
        XSectShape::HorizEllipse => Some((&xt::A_HorizEllipse, xt::N_A_HorizEllipse)),
        XSectShape::VertEllipse => Some((&xt::A_VertEllipse, xt::N_A_VertEllipse)),
        XSectShape::Arch => Some((&xt::A_Arch, xt::N_A_Arch)),
        _ => None,
    }
}

fn area_inv_table_for(shape: XSectShape) -> TableRef {
    match shape {
        XSectShape::Gothic => Some((&xt::Y_Gothic, xt::N_Y_Gothic)),
        XSectShape::Catenary => Some((&xt::Y_Catenary, xt::N_Y_Catenary)),
        XSectShape::Semielliptical => Some((&xt::Y_SemiEllip, xt::N_Y_SemiEllip)),
        XSectShape::Semicircular => Some((&xt::Y_SemiCirc, xt::N_Y_SemiCirc)),
        _ => None,
    }
}

fn area_inv_lut_for(shape: XSectShape) -> Option<&'static LocateLut> {
    use crate::hydraulics::xsect_kernels::LutId;
    match shape {
        XSectShape::Gothic => xsection::host_tables().lut(LutId::YGothic),
        XSectShape::Catenary => xsection::host_tables().lut(LutId::YCatenary),
        XSectShape::Semielliptical => xsection::host_tables().lut(LutId::YSemiEllip),
        XSectShape::Semicircular => xsection::host_tables().lut(LutId::YSemiCirc),
        _ => None,
    }
}

fn hydrad_table_for(shape: XSectShape) -> TableRef {
    match shape {
        XSectShape::Eggshaped => Some((&xt::R_Egg, xt::N_R_Egg)),
        XSectShape::Horseshoe => Some((&xt::R_Horseshoe, xt::N_R_Horseshoe)),
        XSectShape::Baskethandle => Some((&xt::R_Baskethandle, xt::N_R_Baskethandle)),
        XSectShape::HorizEllipse => Some((&xt::R_HorizEllipse, xt::N_R_HorizEllipse)),
        XSectShape::VertEllipse => Some((&xt::R_VertEllipse, xt::N_R_VertEllipse)),
        XSectShape::Arch => Some((&xt::R_Arch, xt::N_R_Arch)),
        _ => None,
    }
}

fn width_table_for(shape: XSectShape) -> TableRef {
    match shape {
        XSectShape::Eggshaped => Some((&xt::W_Egg, xt::N_W_Egg)),
        XSectShape::Horseshoe => Some((&xt::W_Horseshoe, xt::N_W_Horseshoe)),
        XSectShape::Gothic => Some((&xt::W_Gothic, xt::N_W_Gothic)),
        XSectShape::Catenary => Some((&xt::W_Catenary, xt::N_W_Catenary)),
        XSectShape::Semielliptical => Some((&xt::W_SemiEllip, xt::N_W_SemiEllip)),
        XSectShape::Baskethandle => Some((&xt::W_BasketHandle, xt::N_W_BasketHandle)),
        XSectShape::Semicircular => Some((&xt::W_SemiCirc, xt::N_W_SemiCirc)),
        XSectShape::HorizEllipse => Some((&xt::W_HorizEllipse, xt::N_W_HorizEllipse)),
        XSectShape::VertEllipse => Some((&xt::W_VertEllipse, xt::N_W_VertEllipse)),
        XSectShape::Arch => Some((&xt::W_Arch, xt::N_W_Arch)),
        _ => None,
    }
}

// Range form [lo, lo+n): the buffers are the full group-local buffers; parameter arrays
// and buffers are offset by `lo`. Per-element results are position-independent.

/// Group area (ft2) over one slice; per-link branches with an empty table vector leave
/// the outputs unchanged.
fn apply_area_kernel(g: &ShapeGroup, local_d: &[f64], local_a: &mut [f64], lo: i32, n: i32) {
    let ulo = lo as usize;
    let ld = &local_d[ulo..];
    match g.shape {
        XSectShape::Circular | XSectShape::ForceMain => {
            area_circular(
                ld,
                &norm_param(g)[ulo..],
                &g.a_full[ulo..],
                &mut local_a[ulo..],
                n,
            );
        }
        XSectShape::RectClosed | XSectShape::RectOpen => {
            area_rect(ld, &g.w_max[ulo..], &mut local_a[ulo..], n);
        }
        XSectShape::Trapezoidal => {
            area_trapezoidal(ld, &g.y_bot[ulo..], &g.s_bot[ulo..], &mut local_a[ulo..], n);
        }
        XSectShape::Triangular => {
            area_triangular(ld, &g.s_bot[ulo..], &mut local_a[ulo..], n);
        }
        XSectShape::Parabolic => {
            area_parabolic(ld, &g.r_bot[ulo..], &mut local_a[ulo..], n);
        }
        XSectShape::Powerfunc => {
            area_powerfunc(ld, &g.s_bot[ulo..], &g.r_bot[ulo..], &mut local_a[ulo..], n);
        }
        XSectShape::Irregular | XSectShape::Custom | XSectShape::StreetXsect => {
            if !g.area_tables.is_empty() {
                perlink_tabulated(
                    ld,
                    &norm_param(g)[ulo..],
                    &g.a_full[ulo..],
                    &g.area_tables[ulo..],
                    g.transect_tbl_size,
                    &mut local_a[ulo..],
                    n,
                );
            }
        }
        _ => {
            if let Some((tbl, size)) = area_table_for(g.shape) {
                area_tabulated(
                    ld,
                    &norm_param(g)[ulo..],
                    &g.a_full[ulo..],
                    tbl,
                    size,
                    &mut local_a[ulo..],
                    n,
                );
            } else if let Some((inv, size)) = area_inv_table_for(g.shape) {
                area_inv_tabulated(
                    ld,
                    &norm_param(g)[ulo..],
                    &g.a_full[ulo..],
                    inv,
                    size,
                    area_inv_lut_for(g.shape),
                    &mut local_a[ulo..],
                    n,
                );
            } else {
                for k in lo..lo + n {
                    let uk = k as usize;
                    local_a[uk] = xsection::get_a_of_y(&params_at(g, k), local_d[uk]);
                }
            }
        }
    }
}

fn apply_area_kernel_all(g: &ShapeGroup, local_d: &[f64], local_a: &mut [f64]) {
    apply_area_kernel(g, local_d, local_a, 0, g.count);
}

/// Group hydraulic radius (ft) over one slice.
fn apply_hydrad_kernel(g: &ShapeGroup, local_d: &[f64], local_h: &mut [f64], lo: i32, n: i32) {
    let ulo = lo as usize;
    let ld = &local_d[ulo..];
    match g.shape {
        XSectShape::Circular | XSectShape::ForceMain => {
            hydrad_circular(
                ld,
                &norm_param(g)[ulo..],
                &g.r_full[ulo..],
                &mut local_h[ulo..],
                n,
            );
        }
        XSectShape::RectClosed => {
            hydrad_rect_closed(ld, &g.w_max[ulo..], &g.a_full[ulo..], &mut local_h[ulo..], n);
        }
        XSectShape::RectOpen => {
            hydrad_rect_open(ld, &g.w_max[ulo..], &g.s_bot[ulo..], &mut local_h[ulo..], n);
        }
        XSectShape::Trapezoidal => {
            hydrad_trapezoidal(
                ld,
                &g.y_bot[ulo..],
                &g.s_bot[ulo..],
                &g.r_bot[ulo..],
                &mut local_h[ulo..],
                n,
            );
        }
        XSectShape::Triangular => {
            hydrad_triangular(ld, &g.s_bot[ulo..], &g.r_bot[ulo..], &mut local_h[ulo..], n);
        }
        XSectShape::Irregular | XSectShape::Custom | XSectShape::StreetXsect => {
            if !g.hrad_tables.is_empty() {
                perlink_tabulated(
                    ld,
                    &norm_param(g)[ulo..],
                    &g.r_full[ulo..],
                    &g.hrad_tables[ulo..],
                    g.transect_tbl_size,
                    &mut local_h[ulo..],
                    n,
                );
            }
        }
        _ => {
            if let Some((tbl, size)) = hydrad_table_for(g.shape) {
                hydrad_tabulated(
                    ld,
                    &norm_param(g)[ulo..],
                    &g.r_full[ulo..],
                    tbl,
                    size,
                    &mut local_h[ulo..],
                    n,
                );
            } else {
                for k in lo..lo + n {
                    let uk = k as usize;
                    local_h[uk] = xsection::get_r_of_y(&params_at(g, k), local_d[uk]);
                }
            }
        }
    }
}

fn apply_hydrad_kernel_all(g: &ShapeGroup, local_d: &[f64], local_h: &mut [f64]) {
    apply_hydrad_kernel(g, local_d, local_h, 0, g.count);
}

/// Group top width (ft) over one slice.
fn apply_width_kernel(g: &ShapeGroup, local_d: &[f64], local_w: &mut [f64], lo: i32, n: i32) {
    let ulo = lo as usize;
    let ld = &local_d[ulo..];
    match g.shape {
        XSectShape::Circular | XSectShape::ForceMain => {
            width_circular(ld, &norm_param(g)[ulo..], &g.w_max[ulo..], &mut local_w[ulo..], n);
        }
        XSectShape::RectOpen => {
            width_rect(&g.w_max[ulo..], &mut local_w[ulo..], n);
        }
        XSectShape::RectClosed => {
            width_rect_closed(ld, &norm_param(g)[ulo..], &g.w_max[ulo..], &mut local_w[ulo..], n);
        }
        XSectShape::Trapezoidal => {
            width_trapezoidal(ld, &g.y_bot[ulo..], &g.s_bot[ulo..], &mut local_w[ulo..], n);
        }
        XSectShape::Triangular => {
            width_triangular(ld, &g.s_bot[ulo..], &mut local_w[ulo..], n);
        }
        XSectShape::Irregular | XSectShape::Custom | XSectShape::StreetXsect => {
            if !g.width_tables.is_empty() {
                perlink_tabulated(
                    ld,
                    &norm_param(g)[ulo..],
                    &g.w_max[ulo..],
                    &g.width_tables[ulo..],
                    g.transect_tbl_size,
                    &mut local_w[ulo..],
                    n,
                );
            }
        }
        _ => {
            if let Some((tbl, size)) = width_table_for(g.shape) {
                width_tabulated(
                    ld,
                    &norm_param(g)[ulo..],
                    &g.w_max[ulo..],
                    tbl,
                    size,
                    &mut local_w[ulo..],
                    n,
                );
            } else {
                for k in lo..lo + n {
                    let uk = k as usize;
                    local_w[uk] = xsection::get_w_of_y(&params_at(g, k), local_d[uk]);
                }
            }
        }
    }
}

fn apply_width_kernel_all(g: &ShapeGroup, local_d: &[f64], local_w: &mut [f64]) {
    apply_width_kernel(g, local_d, local_w, 0, g.count);
}

/// Group area (ft2) and hydraulic radius (ft) over one slice: the fused circular and
/// per-link pair kernels share one index, other shapes run the two dispatchers.
fn apply_area_hydrad_kernel(
    g: &ShapeGroup,
    local_d: &[f64],
    local_a: &mut [f64],
    local_h: &mut [f64],
    lo: i32,
    n: i32,
) {
    let ulo = lo as usize;
    if g.shape == XSectShape::Circular || g.shape == XSectShape::ForceMain {
        area_hydrad_circular(
            &local_d[ulo..],
            &norm_param(g)[ulo..],
            &g.a_full[ulo..],
            &g.r_full[ulo..],
            &mut local_a[ulo..],
            &mut local_h[ulo..],
            n,
        );
        return;
    }
    // Transect-backed shapes: area and hyd-rad tables share one depth grid (plan A2).
    if (g.shape == XSectShape::Irregular
        || g.shape == XSectShape::Custom
        || g.shape == XSectShape::StreetXsect)
        && !g.area_tables.is_empty()
        && !g.hrad_tables.is_empty()
    {
        perlink_tabulated_pair(
            &local_d[ulo..],
            &norm_param(g)[ulo..],
            &g.a_full[ulo..],
            &g.r_full[ulo..],
            &g.area_tables[ulo..],
            &g.hrad_tables[ulo..],
            g.transect_tbl_size,
            &mut local_a[ulo..],
            &mut local_h[ulo..],
            n,
        );
        return;
    }
    apply_area_kernel(g, local_d, local_a, lo, n);
    apply_hydrad_kernel(g, local_d, local_h, lo, n);
}

fn apply_area_hydrad_kernel_all(g: &ShapeGroup, local_d: &[f64], local_a: &mut [f64], local_h: &mut [f64]) {
    apply_area_hydrad_kernel(g, local_d, local_a, local_h, 0, g.count);
}

/// Team-split slice: worker `tid` of `nthreads` owns [lo, hi) of a count-element group,
/// in 16-element-aligned chunks.
#[inline]
fn team_slice(count: i32, tid: i32, nthreads: i32) -> (i32, i32) {
    let mut chunk = (count + nthreads - 1) / nthreads;
    chunk = (chunk + 15) & !15;
    let lo = cmin(count, tid * chunk);
    let hi = cmin(count, lo + chunk);
    (lo, hi)
}
