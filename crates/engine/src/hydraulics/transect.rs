// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/Transect.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Irregular transect geometry, numerically identical to legacy transect.c and shape.c: the
// three functions Transect.hpp declares (the TransectData record is in transect_types).
//
// Differences in form, never in result:
//   * Names are snake_case (buildTables -> build_tables, buildCustomTables ->
//     build_custom_tables, buildFromStore -> build_from_store).
//   * build_custom_tables takes the curve as two slices plus upstream's int count; it reads
//     exactly the entries upstream reads (indices below n_pts).
//   * getFlow, a lambda capturing the station arrays, is a closure over the same values.
//   * The SWMM_TRACE_TRANSECT dump writes the same bytes: the file name is built as
//     snprintf("%s.ref.%s") builds it (the name up to its first NUL, the whole truncated
//     to 511 bytes) and `%a` is glibc's hex-float form (fmt_a below).

#![allow(
    clippy::manual_memcpy,
    clippy::collapsible_if,
    clippy::collapsible_else_if,
    clippy::manual_range_contains,
    clippy::assign_op_pattern
)]

use std::io::Write as _;

use crate::core::cmath::{cmax, cmin};
use crate::core::file_io::path_from_bytes;
use crate::data::infra::TransectStore;
use crate::hydraulics::transect_types::{TransectData, N_TRANSECT_TBL};
use crate::hydraulics::xsect_lookup::build_invlookup_lut;

const PHI: f64 = 1.486;

// Faithful 1:1 port of legacy transect.c (transect_validate + createTables +
// getGeometry + getSliceGeom + getFlow). AREA and WIDTH are the same geometric
// integrals as a simple area/perimeter sweep, but the hydraulic radius is NOT
// area/perimeter — legacy back-solves an *effective* hyd-radius from the
// divided-channel composite-roughness conveyance, so that (used with the main
// channel n) it reproduces the total normal flow.
/// Build tabulated geometry (area, hydraulic radius, width at N_TRANSECT_TBL equal depth
/// increments, normalized) from station-elevation data. `add_end_walls` closes the section
/// with vertical walls as legacy transect_validate does (false for STREET sections).
pub fn build_tables(td: &mut TransectData, add_end_walls: bool) {
    let n_in = td.stations.len() as i32;
    if n_in < 2 || td.elevations.len() as i32 != n_in {
        return;
    }

    let mut n_channel = td.n_channel;
    if n_channel <= 0.0 {
        return;
    }
    let n_left = if td.n_left > 0.0 { td.n_left } else { n_channel };
    let n_right = if td.n_right > 0.0 { td.n_right } else { n_channel };

    let l_factor = if td.length_factor > 0.0 { td.length_factor } else { 1.0 };
    n_channel *= l_factor.sqrt();

    let x_left_bank = td.x_left_bank;
    let x_right_bank = td.x_right_bank;

    let mut ymin = td.elevations[0];
    let mut ymax = td.elevations[0];
    for i in 1..n_in as usize {
        ymin = cmin(ymin, td.elevations[i]);
        ymax = cmax(ymax, td.elevations[i]);
    }
    td.y_full = ymax - ymin;
    if td.y_full <= 0.0 {
        return;
    }

    // add vertical end-walls reaching full height (legacy transect_validate)
    let n = n_in;
    // N+3, not N+2: getFlow() below reads X[k+1] at k == Nsta == N+1 when the
    // right-bank branch is taken, which is one past the last written element;
    // the trailing element is 0.0, like legacy's zero-initialised slack.
    let mut xs = vec![0.0f64; n as usize + 3];
    let mut ys = vec![0.0f64; n as usize + 3];
    let nsta: i32;
    if add_end_walls {
        xs[0] = td.stations[0];
        ys[0] = ymax;
        for i in 0..n as usize {
            xs[i + 1] = td.stations[i];
            ys[i + 1] = td.elevations[i];
        }
        xs[(n + 1) as usize] = td.stations[(n - 1) as usize];
        ys[(n + 1) as usize] = ymax;
        nsta = n + 1;
    } else {
        // A street's points already include its end walls, so they are laid
        // out exactly as legacy's Station[0..Nstations] with Nstations = N - 1.
        for i in 0..n as usize {
            xs[i] = td.stations[i];
            ys[i] = td.elevations[i];
        }
        nsta = n - 1;
    }

    let get_flow = |k: i32, a: f64, wp: f64, find_flow: bool| -> f64 {
        let mut find_flow = find_flow;
        let ku = k as usize;
        if !find_flow {
            if k == nsta - 1 {
                find_flow = true;
            } else if xs[ku] == x_left_bank {
                if n_left != n_channel && xs[ku] != xs[ku - 1] {
                    find_flow = true;
                }
            } else if xs[ku] == x_right_bank {
                if n_right != n_channel && xs[ku] != xs[ku + 1] {
                    find_flow = true;
                }
            }
        }
        if find_flow {
            let mut nn = n_channel;
            if xs[ku - 1] < x_left_bank {
                nn = n_left;
            }
            if xs[ku] > x_right_bank {
                nn = n_right;
            }
            return PHI / nn * a * (a / wp).powf(2.0 / 3.0); // LIBM
        }
        0.0
    };

    let dy = td.y_full / (N_TRANSECT_TBL - 1) as f64;

    td.area_tbl[0] = 0.0;
    td.hrad_tbl[0] = 0.0;
    td.width_tbl[0] = 0.0;

    // PARITY: legacy createTables (transect.c:281-289) ACCUMULATES the depth
    // `y = ymin; y += dy` each row rather than recomputing `ymin + idx*dy`.
    let mut y = ymin;
    for idx in 1..N_TRANSECT_TBL as usize {
        y += dy;
        let mut wp_sum = 0.0;
        let mut a_sum = 0.0;
        let mut q_sum = 0.0;
        let mut area_t = 0.0;
        let mut width_t = 0.0;

        for k in 1..=nsta {
            let ku = k as usize;
            let ek = ys[ku];
            let ekm = ys[ku - 1];
            let yhi = cmax(ekm, ek);
            let ylo = cmin(ekm, ek);
            if ylo >= y {
                continue;
            }

            let width = (xs[ku] - xs[ku - 1]).abs();
            let mut w = width;
            let mut wp = (width * width + (yhi - ylo) * (yhi - ylo)).sqrt();
            let mut a = 0.0;
            if y > yhi {
                a = width * ((y - yhi) + (y - ylo)) / 2.0;
            } else if yhi > ylo {
                let ratio = (y - ylo) / (yhi - ylo);
                a = width * (yhi - ylo) / 2.0 * ratio * ratio;
                w *= ratio;
                wp *= ratio;
            }

            wp_sum += wp;
            a_sum += a;
            area_t += a;
            width_t += w;

            let find_flow = ek >= y;
            let q = get_flow(k, a_sum, wp_sum, find_flow);
            if q > 0.0 {
                q_sum += q;
                a_sum = 0.0;
                wp_sum = 0.0;
            }
        }

        td.area_tbl[idx] = area_t;
        td.width_tbl[idx] = width_t;
        if area_t == 0.0 {
            td.hrad_tbl[idx] = td.hrad_tbl[idx - 1];
        } else {
            td.hrad_tbl[idx] = (q_sum * n_channel / PHI / area_t).powf(1.5); // LIBM
        }
    }

    // PARITY: legacy setMaxSectionFactor (transect.c:596-618) — the physical
    // maximum section factor and the area at which it occurs, computed on the
    // UNNORMALIZED tables BEFORE the normalization below; stored ABSOLUTE.
    td.s_max = 0.0;
    td.a_max = 0.0;
    for i in 1..N_TRANSECT_TBL as usize {
        let sf = td.area_tbl[i] * td.hrad_tbl[i].powf(2.0 / 3.0); // LIBM
        if sf > td.s_max {
            td.s_max = sf;
            td.a_max = td.area_tbl[i];
        }
    }

    let n_last = (N_TRANSECT_TBL - 1) as usize;
    td.a_full = td.area_tbl[n_last];
    td.r_full = td.hrad_tbl[n_last];
    td.w_max = td.width_tbl[n_last];

    for i in 1..=n_last {
        if td.a_full > 0.0 {
            td.area_tbl[i] /= td.a_full;
        }
        if td.r_full > 0.0 {
            td.hrad_tbl[i] /= td.r_full;
        }
        if td.w_max > 0.0 {
            td.width_tbl[i] /= td.w_max;
        }
    }
    // width at zero depth = width at first increment (legacy createTables:309)
    td.width_tbl[0] = td.width_tbl[1];

    // Bucket map for the getYofA inversion (plan XSECT_LOOKUP_ACCEL A1) — a
    // search accelerator only.
    build_invlookup_lut(&mut td.area_lut, Some(&td.area_tbl[..]), N_TRANSECT_TBL);

    // --- A3 transect-table parity dump (env-gated), mirrors legacy transect.c
    if let Some(p) = std::env::var_os("SWMM_TRACE_TRANSECT") {
        let p = p.as_encoded_bytes();
        if !p.is_empty() {
            // snprintf(fn, 512, "%s.ref.%s", p, td.name.c_str())
            let name = match td.name.iter().position(|&c| c == 0) {
                Some(z) => &td.name[..z],
                None => &td.name[..],
            };
            let mut fname: Vec<u8> = Vec::with_capacity(p.len() + 5 + name.len());
            fname.extend_from_slice(p);
            fname.extend_from_slice(b".ref.");
            fname.extend_from_slice(name);
            fname.truncate(511);
            if let Ok(mut tf) = std::fs::File::create(path_from_bytes(&fname)) {
                let mut s = String::new();
                s.push_str("yFull=");
                s.push_str(&fmt_a(td.y_full));
                s.push_str(" aFull=");
                s.push_str(&fmt_a(td.a_full));
                s.push_str(" rFull=");
                s.push_str(&fmt_a(td.r_full));
                s.push_str(" wMax=");
                s.push_str(&fmt_a(td.w_max));
                s.push('\n');
                for q in 0..N_TRANSECT_TBL as usize {
                    s.push_str(&q.to_string());
                    s.push(',');
                    s.push_str(&fmt_a(td.area_tbl[q]));
                    s.push(',');
                    s.push_str(&fmt_a(td.hrad_tbl[q]));
                    s.push(',');
                    s.push_str(&fmt_a(td.width_tbl[q]));
                    s.push('\n');
                }
                let _ = tf.write_all(s.as_bytes());
            }
        }
    }
}

/// glibc printf's `%a` for a double: "0x1.8p+1", "0x0p+0", "0x0.0000000000001p-1022" for
/// a subnormal, "inf"/"nan" with a '-' when the sign bit is set.
fn fmt_a(v: f64) -> String {
    let bits = v.to_bits();
    let neg = (bits >> 63) != 0;
    let exp = ((bits >> 52) & 0x7ff) as i32;
    let mant = bits & ((1u64 << 52) - 1);
    let mut s = String::new();
    if neg {
        s.push('-');
    }
    if v.is_nan() {
        s.push_str("nan");
        return s;
    }
    if v.is_infinite() {
        s.push_str("inf");
        return s;
    }
    let (lead, e) = if exp == 0 {
        ('0', if mant == 0 { 0 } else { -1022 })
    } else {
        ('1', exp - 1023)
    };
    s.push_str("0x");
    s.push(lead);
    if mant != 0 {
        let hex = format!("{:013x}", mant);
        s.push('.');
        s.push_str(hex.trim_end_matches('0'));
    }
    s.push('p');
    s.push(if e < 0 { '-' } else { '+' });
    s.push_str(&e.unsigned_abs().to_string());
    s
}

// ============================================================================
// PARITY helpers — op-for-op transliterations of legacy src/legacy/engine/
// shape.c. Every expression keeps legacy grouping and operand order so the
// 51-entry CUSTOM geometry tables are bit-identical to legacy TShape tables.
// ============================================================================

// shape.c:307-324 getWidth. NOTE grouping: (y-y1)/(y2-y1) evaluated FIRST,
// then multiplied by (w2-w1).
#[inline]
fn shape_get_width(y: f64, y1: f64, y2: f64, w1: f64, w2: f64) -> f64 {
    if y2 == y1 {
        return w2;
    }
    w1 + (y - y1) / (y2 - y1) * (w2 - w1)
}

// shape.c:328-349 getArea — trapezoid as (wMin + (wMax-wMin)/2)*(y-y1).
#[inline]
fn shape_get_area(y: f64, w: f64, y1: f64, w1: f64) -> f64 {
    let w_min;
    let w_max;
    if w > w1 {
        w_min = w1;
        w_max = w;
    } else {
        w_min = w;
        w_max = w1;
    }
    (w_min + (w_max - w_min) / 2.0) * (y - y1)
}

// shape.c:353-369 getPerim
#[inline]
fn shape_get_perim(y: f64, w: f64, y1: f64, w1: f64) -> f64 {
    let dy = y - y1;
    let dw = (w - w1).abs() / 2.0;
    2.0 * (dy * dy + dw * dw).sqrt()
}

// shape.c:246-303 getNextInterval. Atotal/Ptotal are the legacy file-scope
// statics passed by reference; the cursor ci replaces table_getNextEntry.
// yLast/wLast are by value, exactly like legacy (caller re-syncs from y1/w1).
#[inline]
#[allow(clippy::too_many_arguments)]
fn shape_get_next_interval(
    curve_x: &[f64],
    curve_y: &[f64],
    n_pts: i32,
    ci: &mut i32,
    y: f64,
    y_last: f64,
    w_last: f64,
    y1: &mut f64,
    y2: &mut f64,
    w1: &mut f64,
    w2: &mut f64,
    w_max: &mut f64,
    a_total: &mut f64,
    p_total: &mut f64,
) -> bool {
    let mut y_last = y_last;
    let mut w_last = w_last;
    while y > *y2 {
        // shape.c:272-277 — book the sub-interval up to the curve boundary
        // using the CURVE's endpoint width w2 (no interpolation).
        if *y2 > y_last {
            *a_total += shape_get_area(*y2, *w2, y_last, w_last);
            *p_total += shape_get_perim(*y2, *w2, y_last, w_last);
            y_last = *y2;
            w_last = *w2;
        }
        // shape.c:280-285 — advance to the next curve table interval
        *y1 = *y2;
        *w1 = *w2;
        if *ci >= n_pts {
            // table_getNextEntry failed
            *y2 = 1.0; // shape.c:284
            return true;
        }
        *y2 = curve_x[*ci as usize];
        *w2 = curve_y[*ci as usize];
        *ci += 1;
        if *w2 > *w_max {
            *w_max = *w2; // shape.c:288-290
        }
        if *y2 < *y1 || *w2 < 0.0 {
            return false; // shape.c:293-295
        }
        if *y2 > 1.0 {
            *y2 = 1.0; // shape.c:297-299
        }
    }
    true
}

/// Build tabulated geometry from a CUSTOM shape curve (y/yFull vs w/wMax), as legacy
/// shape_validate plus xsect_setCustomXsectParams' scaling.
pub fn build_custom_tables(
    td: &mut TransectData,
    y_full: f64,
    curve_x: &[f64],
    curve_y: &[f64],
    n_pts: i32,
) {
    if n_pts < 1 || y_full <= 0.0 {
        return;
    }

    // PARITY: op-for-op transliteration of legacy shape.c shape_validate
    //   = computeShapeTables (shape.c:64) + getSmax (shape.c:189)
    //   + normalizeShapeTables (shape.c:216),
    // then scaled to physical units exactly as xsect_setCustomXsectParams
    // does (xsect.c:681-686).

    td.y_full = y_full;

    // --- get first entry of user's shape curve (shape.c:81-89)
    let mut y1 = curve_x[0];
    let mut w1 = curve_y[0];
    if y1 < 0.0 || y1 >= 1.0 || w1 < 0.0 {
        return;
    }
    let mut w_max = w1;

    let mut y2: f64;
    let mut w2: f64;
    let mut ci: i32 = 1; // curve cursor == legacy table iterator position

    // --- if first entry not at zero ht. then add an initial entry (shape.c:92-97)
    if y1 != 0.0 {
        y2 = y1;
        w2 = w1;
        y1 = 0.0;
        w1 = 0.0;
    }
    // --- otherwise get next entry in the user's shape curve (shape.c:99-115)
    else {
        if ci >= n_pts {
            return;
        }
        y2 = curve_x[ci as usize];
        w2 = curve_y[ci as usize];
        ci += 1;
        if y2 < y1 || w2 < 0.0 {
            return;
        }
        if y2 > 1.0 {
            y2 = 1.0;
        }
        if w2 > w_max {
            w_max = w2;
        }
    }

    // --- number of entries & interval size in geom. tables (shape.c:118-120)
    let n = (N_TRANSECT_TBL - 1) as usize; // == legacy N_SHAPE_TBL - 1
    let dy = 1.0 / n as f64;

    // --- initialize geometry tables (shape.c:123-127)
    td.area_tbl[0] = 0.0;
    td.hrad_tbl[0] = 0.0;
    td.width_tbl[0] = w1;
    let mut p_total = w1;
    let mut a_total = 0.0;

    // --- fill in rest of geometry tables (shape.c:130-174)
    let mut y = 0.0;
    let mut w = w1;
    for i in 1..=n {
        let mut y_last = y;
        let mut w_last = w;
        y = y + dy;

        // do not allow height to exceed 1.0 (shape.c:140-142, TINY = 1e-6)
        if (y - 1.0).abs() < 1.0e-6 {
            y = 1.0;
        }

        // if height exceeds current curve interval, move to next one
        // (shape.c:146-154 — an if, not a loop; getNextInterval loops inside)
        if y > y2 {
            if !shape_get_next_interval(
                curve_x,
                curve_y,
                n_pts,
                &mut ci,
                y,
                y_last,
                w_last,
                &mut y1,
                &mut y2,
                &mut w1,
                &mut w2,
                &mut w_max,
                &mut a_total,
                &mut p_total,
            ) {
                return;
            }
            y_last = y1;
            w_last = w1;
        }

        // top width, area & perimeter of current interval (shape.c:157-159)
        w = shape_get_width(y, y1, y2, w1, w2);
        a_total += shape_get_area(y, w, y_last, w_last);
        p_total += shape_get_perim(y, w, y_last, w_last);

        // add top width to total perimeter if at top of shape
        // (shape.c:162-164 — legacy adds w2, NOT the interpolated w)
        if y == 1.0 {
            p_total += w2;
        }

        // update table values (shape.c:167-173)
        td.width_tbl[i] = w;
        td.area_tbl[i] = a_total;
        if p_total > 0.0 {
            td.hrad_tbl[i] = a_total / p_total;
        } else {
            td.hrad_tbl[i] = 0.0;
        }
    }

    // --- unit-height full-depth properties (shape.c:177-178)
    let a_full = td.area_tbl[n];
    let r_full = td.hrad_tbl[n];

    // --- max. section factor & its area (getSmax, shape.c:189-212), computed
    //     on the UNNORMALIZED tables before normalization, like legacy.
    let mut s_max = 0.0;
    let mut a_max = 0.0;
    for i in 1..=n {
        let sf = td.area_tbl[i] * td.hrad_tbl[i].powf(2.0 / 3.0); // LIBM
        if sf > s_max {
            s_max = sf;
            a_max = td.area_tbl[i];
        }
    }
    td.s_max = s_max; // unit-height values; caller scales per xsect.c:685-686
    td.a_max = a_max;

    // --- normalize tables (normalizeShapeTables, shape.c:229-239)
    if a_full == 0.0 || r_full == 0.0 || w_max == 0.0 {
        return;
    }
    for i in 0..=n {
        td.area_tbl[i] /= a_full;
        td.hrad_tbl[i] /= r_full;
        td.width_tbl[i] /= w_max;
    }

    // --- scale to physical units (xsect_setCustomXsectParams, xsect.c:681-683)
    td.w_max = w_max * y_full;
    td.a_full = a_full * y_full * y_full;
    td.r_full = r_full * y_full;

    build_invlookup_lut(&mut td.area_lut, Some(&td.area_tbl[..]), N_TRANSECT_TBL);
}

/// Copy one [TRANSECTS] store row into `td`, apply the legacy parse-time transforms to the
/// copy, then build_tables(); false when the channel Manning's n is invalid.
pub fn build_from_store(ts: &TransectStore, index: i32, ucf_length: f64, td: &mut TransectData) -> bool {
    let ut = index as usize;
    td.name = ts.names[ut].clone();
    td.n_left = ts.n_left[ut];
    td.n_right = ts.n_right[ut];
    td.n_channel = ts.n_channel[ut];
    if td.n_channel <= 0.0 {
        return false;
    }
    td.stations = ts.stations[ut].clone();
    td.elevations = ts.elevations[ut].clone();
    td.length_factor = ts.length_factor[ut];
    // PARITY: legacy transect setParams/addStation (transect.c:360-410)
    // transform the raw [TRANSECTS] GR data BEFORE building the tables:
    //   Station = x * Xfactor / UCF(LENGTH)          (mult, then divide)
    //   Elev    = (y + Yfactor) / UCF(LENGTH),  Yfactor = x9 / UCF
    //   Xbank   = (xbank / UCF) * Xfactor            (divide, then mult)
    let mut x_factor = ts.x_factor[ut]; // parse-time 0→1 default
    if x_factor == 0.0 {
        x_factor = 1.0;
    }
    let y_factor = ts.y_factor[ut] / ucf_length; // x9 / UCF
    for s in td.stations.iter_mut() {
        *s = *s * x_factor / ucf_length;
    }
    for e in td.elevations.iter_mut() {
        *e = (*e + y_factor) / ucf_length;
    }
    td.x_left_bank = (ts.x_left_bank[ut] / ucf_length) * x_factor;
    td.x_right_bank = (ts.x_right_bank[ut] / ucf_length) * x_factor;
    build_tables(td, true);
    true
}

#[cfg(test)]
mod tests {
    use super::fmt_a;
    // Expected strings recorded from glibc printf("%a") on the reference machine.
    #[test]
    fn fmt_a_matches_glibc() {
        let cases: [(u64, &str); 20] = [
            (0x0000000000000000, "0x0p+0"),
            (0x8000000000000000, "-0x0p+0"),
            (0x3ff0000000000000, "0x1p+0"),
            (0xbff0000000000000, "-0x1p+0"),
            (0x3fe0000000000000, "0x1p-1"),
            (0x4008000000000000, "0x1.8p+1"),
            (0x3fb999999999999a, "0x1.999999999999ap-4"),
            (0x7e37e43c8800759c, "0x1.7e43c8800759cp+996"),
            (0x81a56e1fc2f8f359, "-0x1.56e1fc2f8f359p-997"),
            (0x0000000000000001, "0x0.0000000000001p-1022"),
            (0x0008000000000000, "0x0.8p-1022"),
            (0x0010000000000000, "0x1p-1022"),
            (0x000fffffffffffff, "0x0.fffffffffffffp-1022"),
            (0x7ff0000000000000, "inf"),
            (0xfff0000000000000, "-inf"),
            (0x7ff8000000000000, "nan"),
            (0xfff8000000000000, "-nan"),
            (0x40fe240c9fbe76c9, "0x1.e240c9fbe76c9p+16"),
            (0x3fd5555555555555, "0x1.5555555555555p-2"),
            (0x7fefffffffffffff, "0x1.fffffffffffffp+1023"),
        ];
        for (bits, s) in cases {
            assert_eq!(fmt_a(f64::from_bits(bits)), s);
        }
    }
}
