// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/surface/SolarGeometry.cpp and
// SolarGeometry.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright
// 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 6 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// The exact interval integral of the solar-vector cosine on horizontal ground and on an
// inclined plane, between the roots of each harmonic. std::invalid_argument is Err(what).
// std::remainder has no Rust counterpart; `ieee_remainder` computes it exactly (fmod, then the
// half-way correction to the even quotient), as IEEE 754 defines it.

const PI: f64 = 3.14159265358979323846;
const RAD: f64 = PI / 180.0;

/// SolarIntegral
#[derive(Debug, Clone, Copy, Default)]
pub struct SolarIntegral {
    pub cosine: f64,
    pub horizontal: f64,
    pub hri: f64,
}

#[derive(Clone, Copy)]
struct Harmonic {
    a: f64,
    b: f64,
    c: f64,
}

impl Harmonic {
    fn value(&self, x: f64) -> f64 {
        self.a + self.b * x.cos() + self.c * x.sin() // LIBM
    }
    fn primitive(&self, x: f64) -> f64 {
        self.a * x + self.b * x.sin() - self.c * x.cos() // LIBM
    }
    fn roots(&self, lo: f64, hi: f64, cuts: &mut Vec<f64>) {
        let r = self.b.hypot(self.c); // LIBM
        if r == 0.0 || self.a.abs() >= r {
            return;
        }
        let offset = self.c.atan2(self.b); // LIBM
        let angle = (-self.a / r).acos(); // LIBM
        for base in [offset - angle, offset + angle] {
            let first = crate::core::cmath::cvt_i32(((lo - base) / (2.0 * PI)).ceil());
            let last = crate::core::cmath::cvt_i32(((hi - base) / (2.0 * PI)).floor());
            let mut k = first;
            while k <= last {
                let x = base + (k * 2) as f64 * PI;
                if x > lo && x < hi {
                    cuts.push(x);
                }
                k += 1;
            }
        }
    }
}

/// std::remainder(x, y) for finite x and nonzero finite y.
pub fn ieee_remainder(x: f64, y: f64) -> f64 {
    let ay = y.abs();
    let mut r = x % y;
    let ar = r.abs();
    let half = ay * 0.5;
    let odd = {
        let q = ((x - r) / y).abs();
        q % 2.0 == 1.0
    };
    if ar > half || (ar == half && odd) {
        r -= if r > 0.0 { ay } else { -ay };
    }
    r
}

/// solarIntegral (degrees; aspect clockwise from north; hour angle 0 at solar noon).
pub fn solar_integral(latitude: f64, declination: f64, hour_start: f64, duration_hours: f64, slope: f64, aspect: f64) -> Result<SolarIntegral, &'static str> {
    for v in [latitude, declination, hour_start, duration_hours, slope, aspect] {
        if !v.is_finite() {
            return Err("Nonfinite solar geometry.");
        }
    }
    if latitude.abs() > 90.0 || declination.abs() > 90.0 || duration_hours <= 0.0 || duration_hours > 24.0 || slope < 0.0 || slope >= 90.0 {
        return Err("Invalid solar geometry range.");
    }
    // LIBM (sin, cos)
    let lat = latitude * RAD;
    let dec = declination * RAD;
    let tilt = slope * RAD;
    let azi = aspect * RAD;
    let north = tilt.sin() * azi.cos();
    let east = tilt.sin() * azi.sin();
    let up = tilt.cos();
    let horizontal = Harmonic { a: lat.sin() * dec.sin(), b: lat.cos() * dec.cos(), c: 0.0 };
    let plane = Harmonic { a: dec.sin() * (north * lat.cos() + up * lat.sin()), b: dec.cos() * (up * lat.cos() - north * lat.sin()), c: -east * dec.cos() };
    let lo = ieee_remainder(hour_start, 360.0) * RAD;
    let hi = lo + duration_hours * PI / 12.0;
    let mut cuts = vec![lo, hi];
    horizontal.roots(lo, hi, &mut cuts);
    plane.roots(lo, hi, &mut cuts);
    cuts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut result = SolarIntegral::default();
    for i in 1..cuts.len() {
        let a = cuts[i - 1];
        let b = cuts[i];
        let mid = 0.5 * (a + b);
        if horizontal.value(mid) <= 0.0 {
            continue;
        }
        result.horizontal += horizontal.primitive(b) - horizontal.primitive(a);
        if plane.value(mid) > 0.0 {
            result.cosine += plane.primitive(b) - plane.primitive(a);
        }
    }
    result.cosine = crate::core::cmath::cclamp(result.cosine / (hi - lo), 0.0, 1.0);
    result.horizontal = crate::core::cmath::cclamp(result.horizontal / (hi - lo), 0.0, 1.0);
    result.hri = result.cosine / up;
    Ok(result)
}
