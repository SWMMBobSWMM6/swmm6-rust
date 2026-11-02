// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/DateTime.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Dates are f64 days since 30 Dec 1899 (the OLE / Delphi convention, as in SWMM 5).
// `static_cast<int>(std::floor(x))` is `x.floor() as i32`: the same for every in-range
// value; out of range it is undefined in C++ and saturates here.

pub const DATE_DELTA: i32 = 693594;
pub const SECS_PER_DAY: f64 = 86400.0;
pub const ONE_SECOND: f64 = 1.1574074e-5;

pub type DateTime = f64;

const DAYS_PER_MONTH: [[i32; 12]; 2] = [
    [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
    [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
];

/// (n / d, n - d * (n / d)), or (0, 0) when d is 0. C++ integer division truncates, as `/` does.
#[inline]
pub fn div_mod(n: i32, d: i32) -> (i32, i32) {
    if d == 0 {
        (0, 0)
    } else {
        let r = n / d;
        (r, n - d * r)
    }
}

#[inline]
pub fn is_leap_year(year: i32) -> bool {
    (year % 4 == 0) && ((year % 100 != 0) || (year % 400 == 0))
}

pub fn encode_date(year: i32, month: i32, mut day: i32) -> DateTime {
    let i = if is_leap_year(year) { 1 } else { 0 };
    if (1..=9999).contains(&year) && (1..=12).contains(&month) && day >= 1 && day <= DAYS_PER_MONTH[i][(month - 1) as usize] {
        for j in 0..(month - 1) as usize {
            day += DAYS_PER_MONTH[i][j];
        }
        let y = year - 1;
        return (y * 365 + y / 4 - y / 100 + y / 400 + day - DATE_DELTA) as f64;
    }
    -DATE_DELTA as f64
}

pub fn encode_time(hour: i32, minute: i32, second: i32) -> DateTime {
    if hour >= 0 && minute >= 0 && second >= 0 {
        let s = hour * 3600 + minute * 60 + second;
        return s as f64 / SECS_PER_DAY;
    }
    0.0
}

/// (year, month, day).
pub fn decode_date(date: DateTime) -> (i32, i32, i32) {
    const D1: i32 = 365;
    const D4: i32 = 1461;
    const D100: i32 = 36524;
    const D400: i32 = 146097;
    let mut t = (date.floor() as i32).wrapping_add(DATE_DELTA);
    if t <= 0 {
        return (0, 1, 1);
    }
    t -= 1;
    let mut y = 1;
    while t >= D400 {
        t -= D400;
        y += 400;
    }
    let (mut i, mut d) = div_mod(t, D100);
    if i == 4 {
        i -= 1;
        d += D100;
    }
    y += i * 100;
    let (i2, d2) = div_mod(d, D4);
    y += i2 * 4;
    let (mut i3, mut d3) = div_mod(d2, D1);
    if i3 == 4 {
        i3 -= 1;
        d3 += D1;
    }
    y += i3;
    let k = if is_leap_year(y) { 1 } else { 0 };
    let mut m = 1;
    loop {
        let len = DAYS_PER_MONTH[k][(m - 1) as usize];
        if d3 < len {
            break;
        }
        d3 -= len;
        m += 1;
    }
    (y, m, d3 + 1)
}

/// (hours, minutes, seconds), rounded to the nearest second, never past 23:59:59.
pub fn decode_time(time: DateTime) -> (i32, i32, i32) {
    let frac_day = (time - time.floor()) * SECS_PER_DAY;
    let mut secs = (frac_day + 0.5).floor() as i32;
    if secs >= 86400 {
        secs = 86399;
    }
    let (mins, s) = div_mod(secs, 60);
    let (mut h, m) = div_mod(mins, 60);
    if h > 23 {
        h = 0;
    }
    (h, m, s)
}

pub fn add_seconds(date: DateTime, seconds: f64) -> DateTime {
    let d = date.floor();
    let (h, m, s) = decode_time(date);
    d + (3600.0 * h as f64 + 60.0 * m as f64 + s as f64 + seconds) / SECS_PER_DAY
}

pub fn time_diff(date1: DateTime, date2: DateTime) -> i64 {
    let d1 = date1.floor();
    let d2 = date2.floor();
    let (h, m, s) = decode_time(date1);
    let s1 = 3600i64 * h as i64 + 60i64 * m as i64 + s as i64;
    let (h, m, s) = decode_time(date2);
    let s2 = 3600i64 * h as i64 + 60i64 * m as i64 + s as i64;
    let mut secs = ((d1 - d2) * SECS_PER_DAY + 0.5).floor() as i64;
    secs += s1 - s2;
    secs
}

pub fn month_of_year(date: DateTime) -> i32 {
    decode_date(date).1
}

pub fn day_of_year(date: DateTime) -> i32 {
    let (y, _, _) = decode_date(date);
    let start_of_year = encode_date(y, 1, 1);
    (date - start_of_year).floor() as i32 + 1
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn round_trip() {
        assert_eq!(encode_date(2004, 1, 1), 37987.0);
        assert_eq!(decode_date(37987.0), (2004, 1, 1));
        assert_eq!(decode_date(encode_date(2000, 2, 29)), (2000, 2, 29));
        assert_eq!(decode_time(encode_time(13, 45, 30)), (13, 45, 30));
        assert_eq!(encode_date(2001, 2, 29), -693594.0);
    }
}
