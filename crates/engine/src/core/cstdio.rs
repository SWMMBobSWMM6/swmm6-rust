// SPDX-License-Identifier: Apache-2.0
//
// The C and C++ stdio pieces upstream's input code calls, as glibc 2.39 / libstdc++
// (GCC 13) perform them in the C locale:
//   * printf_g: `%.<prec>g` (snprintf/printf) on a double. The digits come from Rust's
//     exact formatting, which rounds the exact binary value half-to-even as glibc's
//     printf_fp does; the %e/%f choice, trailing-zero removal, two-digit signed exponent,
//     "-0", "inf" and "-nan" are glibc's.
//   * IfStream: `std::ifstream in(path)` read with `std::getline(in, line)`. The path is a
//     C string (it ends at the first NUL). A read error sets badbit, as the exception
//     libstdc++'s underflow throws does inside getline; a line cut short by it is not
//     returned.
//   * to_string_long: `std::to_string(long)`.
// Written for swmm6-rust, 2 November 2026, by Claude (Anthropic) for Robert Dickinson.

use crate::core::cstdlib::c_str;
use crate::core::file_io::path_from_bytes;
use crate::core::options::Bytes;
use std::io::Read;

/// `snprintf(buf, n, "%.<prec>g", v)` for a buffer large enough.
pub fn printf_g(v: f64, prec: i32) -> Bytes {
    let neg = v.is_sign_negative();
    let mut out = Bytes::new();
    if v.is_nan() {
        if neg {
            out.push(b'-');
        }
        out.extend_from_slice(b"nan");
        return out;
    }
    if neg {
        out.push(b'-');
    }
    let a = v.abs();
    if a.is_infinite() {
        out.extend_from_slice(b"inf");
        return out;
    }
    let p: i32 = if prec == 0 { 1 } else { prec };
    // Exponent of the value rounded to p significant digits.
    let e_str = format!("{:.*e}", (p - 1) as usize, a);
    let epos = e_str.find('e').unwrap();
    let x: i32 = e_str[epos + 1..].parse().unwrap();
    let body: String = if p > x && x >= -4 {
        format!("{:.*}", (p - 1 - x) as usize, a)
    } else {
        let mant = &e_str[..epos];
        let sign = if x < 0 { '-' } else { '+' };
        format!("{}e{}{:02}", mant, sign, x.abs())
    };
    // Remove trailing zeros of the fraction, then a bare decimal point.
    let (num, exp) = match body.find('e') {
        Some(k) => (&body[..k], &body[k..]),
        None => (&body[..], ""),
    };
    let mut num = num.to_string();
    if num.contains('.') {
        while num.ends_with('0') {
            num.pop();
        }
        if num.ends_with('.') {
            num.pop();
        }
    }
    out.extend_from_slice(num.as_bytes());
    out.extend_from_slice(exp.as_bytes());
    out
}

/// `std::to_string` of an integer.
pub fn to_string_long(v: i64) -> Bytes {
    v.to_string().into_bytes()
}

/// An input file stream read whole on open and handed out line by line.
pub struct IfStream {
    data: Bytes,
    pos: usize,
    /// a read failed after `data`
    bad: bool,
}

impl IfStream {
    /// `std::ifstream in(path)`; None where `!in` (is_open() false).
    pub fn open(path: &[u8]) -> Option<IfStream> {
        let mut f = std::fs::File::open(path_from_bytes(c_str(path))).ok()?;
        let mut data = Bytes::new();
        let mut chunk = [0u8; 65536];
        let mut bad = false;
        loop {
            match f.read(&mut chunk) {
                Ok(0) => break,
                Ok(k) => data.extend_from_slice(&chunk[..k]),
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
                Err(_) => {
                    bad = true;
                    break;
                }
            }
        }
        Some(IfStream { data, pos: 0, bad })
    }

    /// `std::getline(in, line)` as a loop condition: false once nothing is left.
    pub fn getline(&mut self, line: &mut Bytes) -> bool {
        if self.pos >= self.data.len() {
            return false;
        }
        let rest = &self.data[self.pos..];
        match rest.iter().position(|&c| c == b'\n') {
            Some(k) => {
                *line = rest[..k].to_vec();
                self.pos += k + 1;
                true
            }
            None => {
                if self.bad {
                    // the read error interrupted this line: getline fails
                    self.pos = self.data.len();
                    return false;
                }
                *line = rest.to_vec();
                self.pos = self.data.len();
                true
            }
        }
    }

    /// `in.bad()`, meaningful once getline has returned false.
    pub fn bad(&self) -> bool {
        self.bad
    }
}

/// Differential test of data::table::sscanf (and so cstdlib::strtof) against glibc: set
/// SSCANF_GATE_IN to the input lines and SSCANF_GATE_OUT to where the Rust side writes
/// what gate/diff/scan.c prints for them. Skipped when the variables are unset.
#[cfg(test)]
mod sscanf_gate {
    use crate::data::table::{sscanf, ScanArg::*};
    use std::fmt::Write;

    fn s(b: &[u8]) -> String {
        String::from_utf8_lossy(b).into_owned()
    }

    #[test]
    fn sscanf_against_glibc() {
        let (inp, outp) = match (std::env::var("SSCANF_GATE_IN"), std::env::var("SSCANF_GATE_OUT")) {
            (Ok(a), Ok(b)) => (a, b),
            _ => return,
        };
        let data = std::fs::read(inp).unwrap();
        let mut o = String::new();
        for raw in data.split(|&c| c == b'\n') {
            let line: Vec<u8> = raw.iter().map(|&c| if c == b'~' { b'\n' } else { c }).collect();
            let l = &line[..];
            let fl = |f: f32| format!(" f{:08x}", f.to_bits());
            {
                let (mut m, mut d, mut y, mut s1, mut s2) = (7u32, 7u32, 7u32, b'#', b'#');
                let r = sscanf(l, b"%u%c%u%c%u", &mut [UInt(&mut m), Char(&mut s1), UInt(&mut d), Char(&mut s2), UInt(&mut y)]);
                writeln!(o, "A {} {} {} {} {} {}", r, m, s1 as i8, d, s2 as i8, y).unwrap();
            }
            {
                let (mut mon, mut d, mut y, mut s1, mut s2) = (b"?".to_vec(), 7u32, 7u32, b'#', b'#');
                let r = sscanf(l, b"%3[A-Za-z]%c%u%c%u", &mut [Str(&mut mon), Char(&mut s1), UInt(&mut d), Char(&mut s2), UInt(&mut y)]);
                writeln!(o, "B {} {} {} {} {} {}", r, s(&mon), s1 as i8, d, s2 as i8, y).unwrap();
            }
            {
                let (mut h, mut m, mut sec) = (7, 7, 7);
                let r = sscanf(l, b"%d:%d:%d", &mut [Int(&mut h), Int(&mut m), Int(&mut sec)]);
                writeln!(o, "C {} {} {} {}", r, h, m, sec).unwrap();
            }
            {
                let (mut a, mut b, mut e) = (7i64, 7, b"?".to_vec());
                let r = sscanf(l, b"%6ld %2d %4s", &mut [Long(&mut a), Int(&mut b), Str(&mut e)]);
                writeln!(o, "D {} {} {} {}", r, a, b, s(&e)).unwrap();
            }
            {
                let (mut a, mut b, mut e) = (7i64, 7, b"?".to_vec());
                let r = sscanf(l, b"%6ld,%2d,%4s", &mut [Long(&mut a), Int(&mut b), Str(&mut e)]);
                writeln!(o, "E {} {} {} {}", r, a, b, s(&e)).unwrap();
            }
            {
                let (mut b, mut y, mut e, mut rc) = (7, 7, b"?".to_vec(), b"?".to_vec());
                let r = sscanf(l, b"%2d %4s %2s %4d", &mut [Int(&mut b), Str(&mut e), Str(&mut rc), Int(&mut y)]);
                writeln!(o, "F {} {} {} {} {}", r, b, s(&e), s(&rc), y).unwrap();
            }
            {
                let (mut b, mut y, mut e, mut rc) = (7, 7, b"?".to_vec(), b"?".to_vec());
                let r = sscanf(l, b"%2d,%4s,%2s,%4d", &mut [Int(&mut b), Str(&mut e), Str(&mut rc), Int(&mut y)]);
                writeln!(o, "G {} {} {} {} {}", r, b, s(&e), s(&rc), y).unwrap();
            }
            {
                let (mut a, mut b, mut e, mut rc) = (7i64, 7, b"?".to_vec(), b"?".to_vec());
                let r = sscanf(l, b"%3s%6ld%2d%4s", &mut [Str(&mut rc), Long(&mut a), Int(&mut b), Str(&mut e)]);
                writeln!(o, "H {} {} {} {} {}", r, s(&rc), a, b, s(&e)).unwrap();
            }
            {
                let (mut a, mut c) = (7i64, b"?".to_vec());
                let r = sscanf(l, b"%5s%6ld", &mut [Str(&mut c), Long(&mut a)]);
                writeln!(o, "I {} {} {}", r, s(&c), a).unwrap();
            }
            {
                let (mut x, mut c) = (99.99f32, b'#');
                let r = sscanf(l, b"%f %c", &mut [Float(&mut x), Char(&mut c)]);
                writeln!(o, "J {} {}{}", r, c as i8, fl(x)).unwrap();
            }
            {
                let (mut v, mut c) = (7i64, b'#');
                let r = sscanf(l, b"%ld %c", &mut [Long(&mut v), Char(&mut c)]);
                writeln!(o, "K {} {} {}", r, v, c as i8).unwrap();
            }
            {
                let (mut a, mut b, mut c, mut d) = (7, 7, 7, 7);
                let r = sscanf(l, b"%4d%2d%4d%3d", &mut [Int(&mut a), Int(&mut b), Int(&mut c), Int(&mut d)]);
                writeln!(o, "L {} {} {} {} {}", r, a, b, c, d).unwrap();
            }
            for (tag, f) in [("M", &b"%4d %2d %2d"[..]), ("N", &b"%4d,%2d,%2d"[..]), ("O", &b"%4d%2d%2d"[..])] {
                let (mut a, mut b, mut c) = (7, 7, 7);
                let r = sscanf(l, f, &mut [Int(&mut a), Int(&mut b), Int(&mut c)]);
                writeln!(o, "{} {} {} {} {}", tag, r, a, b, c).unwrap();
            }
            for (tag, f) in [("P", &b"%2d%2d%6ld%c%c"[..]), ("Q", &b" %2d%2d %6ld %c %c"[..]), ("R", &b",%2d%2d,%6ld,%c,%c"[..])] {
                let (mut h, mut m, mut v, mut f1, mut f2) = (7, 7, 7i64, b'#', b'#');
                let r = sscanf(l, f, &mut [Int(&mut h), Int(&mut m), Long(&mut v), Char(&mut f1), Char(&mut f2)]);
                writeln!(o, "{} {} {} {} {} {} {}", tag, r, h, m, v, f1 as i8, f2 as i8).unwrap();
            }
            {
                let (mut h, mut m) = (7, 7);
                let r = sscanf(l, b" %2d:%2d", &mut [Int(&mut h), Int(&mut m)]);
                writeln!(o, "S {} {} {}", r, h, m).unwrap();
            }
            {
                let (mut t, mut a, mut b, mut c, mut d, mut e, mut v) = (b"?".to_vec(), 7, 7, 7, 7, 7, -1.5f32);
                let r = sscanf(
                    l,
                    b"%63s %d %d %d %d %d %f",
                    &mut [Str(&mut t), Int(&mut a), Int(&mut b), Int(&mut c), Int(&mut d), Int(&mut e), Float(&mut v)],
                );
                writeln!(o, "T {} {} {} {} {} {} {}{}", r, s(&t), a, b, c, d, e, fl(v)).unwrap();
            }
            {
                let (mut m, mut d, mut y, mut e) = (7, 7, 7, -9);
                let r = sscanf(l, b"%d/%d/%d%n", &mut [Int(&mut m), Int(&mut d), Int(&mut y), Count(&mut e)]);
                writeln!(o, "U {} {} {} {} {}", r, m, d, y, e).unwrap();
            }
            {
                let mut st = b"?".to_vec();
                let r = sscanf(l, b"%3[^a-c-]", &mut [Str(&mut st)]);
                writeln!(o, "V {} {}", r, s(&st)).unwrap();
            }
            {
                let mut st = b"?".to_vec();
                let r = sscanf(l, b"%[]a-]", &mut [Str(&mut st)]);
                writeln!(o, "W {} {}", r, s(&st)).unwrap();
            }
            {
                let mut x = 5.0f32;
                let r = sscanf(l, b"%3f", &mut [Float(&mut x)]);
                writeln!(o, "X {}{}", r, fl(x)).unwrap();
            }
        }
        std::fs::write(outp, o).unwrap();
    }
}
