// SPDX-License-Identifier: Apache-2.0
//
// C's printf family as glibc 2.39 performs it in the C locale, and a C stdio output stream,
// for the .out and .rpt writers (plugins/DefaultOutputPlugin.cpp, DefaultReportPlugin.cpp)
// and the engine's own messages.
//   * format(fmt, args): sprintf. Flags - + space # 0, width and precision (also `*`), the
//     length modifiers h l ll L q j z t, and the conversions d i u x X o c s f F e E g G p %.
//     The digits of %f/%e/%g come from Rust's exact formatting, which rounds the exact
//     binary value half to even as glibc's printf_fp does; the %e/%f choice of %g,
//     trailing-zero removal, the two-digit signed exponent, "-0", "inf" and "-nan" are
//     glibc's. Ported from swmm5-rust's cio.rs (the same project family), where it is gated
//     byte for byte against EPA SWMM 5.2.4's reports.
//   * CFile: a FILE* opened for writing ("w", "w+b"), buffered, with ftell; write errors
//     are ignored, as an unchecked fprintf ignores them.
// Written for swmm6-rust, 5 November 2026, by Claude (Anthropic) for Robert Dickinson.

use crate::core::options::Bytes;
use std::io::Write;

/// One argument of a variadic C call, after C's default promotions (float becomes
/// double, char and short become int).
#[derive(Clone, Copy, Debug)]
pub enum Arg<'a> {
    I(i64),
    U(u64),
    F(f64),
    S(&'a [u8]),
}

/// Conversion of a Rust value into the argument C would have received.
pub trait ToArg {
    fn to_arg(&self) -> Arg<'_>;
}
impl ToArg for i32 {
    fn to_arg(&self) -> Arg<'_> {
        Arg::I(*self as i64)
    }
}
impl ToArg for i64 {
    fn to_arg(&self) -> Arg<'_> {
        Arg::I(*self)
    }
}
impl ToArg for u32 {
    fn to_arg(&self) -> Arg<'_> {
        Arg::U(*self as u64)
    }
}
impl ToArg for u64 {
    fn to_arg(&self) -> Arg<'_> {
        Arg::U(*self)
    }
}
impl ToArg for usize {
    fn to_arg(&self) -> Arg<'_> {
        Arg::U(*self as u64)
    }
}
/// A C `char` passed to `%c` (or `%d`).
impl ToArg for u8 {
    fn to_arg(&self) -> Arg<'_> {
        Arg::I(*self as i64)
    }
}
impl ToArg for f64 {
    fn to_arg(&self) -> Arg<'_> {
        Arg::F(*self)
    }
}
impl ToArg for f32 {
    fn to_arg(&self) -> Arg<'_> {
        Arg::F(*self as f64)
    }
}
impl ToArg for [u8] {
    fn to_arg(&self) -> Arg<'_> {
        Arg::S(self)
    }
}
impl ToArg for Vec<u8> {
    fn to_arg(&self) -> Arg<'_> {
        Arg::S(self.as_slice())
    }
}
impl<const N: usize> ToArg for [u8; N] {
    fn to_arg(&self) -> Arg<'_> {
        Arg::S(&self[..])
    }
}
impl ToArg for str {
    fn to_arg(&self) -> Arg<'_> {
        Arg::S(self.as_bytes())
    }
}
impl<T: ToArg + ?Sized> ToArg for &T {
    fn to_arg(&self) -> Arg<'_> {
        (**self).to_arg()
    }
}
impl<T: ToArg + ?Sized> ToArg for &mut T {
    fn to_arg(&self) -> Arg<'_> {
        (**self).to_arg()
    }
}

/// `sprintf` into a new byte vector: `cformat!(b"%10.3f", x)`.
#[macro_export]
macro_rules! cformat {
    ($fmt:expr $(, $a:expr)* $(,)?) => {
        $crate::core::cprintf::format($fmt, &[$($crate::core::cprintf::ToArg::to_arg(&$a)),*])
    };
}

/// `fprintf(f, ...)` to a CFile: `cfprintf!(f, b"%d", n)`.
#[macro_export]
macro_rules! cfprintf {
    ($f:expr, $fmt:expr $(, $a:expr)* $(,)?) => {
        $f.write_bytes(&$crate::core::cprintf::format($fmt, &[$($crate::core::cprintf::ToArg::to_arg(&$a)),*]))
    };
}

/// The bytes of a C string: up to the first NUL.
fn cstr(s: &[u8]) -> &[u8] {
    match s.iter().position(|&c| c == 0) {
        Some(k) => &s[..k],
        None => s,
    }
}

struct Spec {
    left: bool,
    plus: bool,
    space: bool,
    alt: bool,
    zero: bool,
    width: usize,
    prec: Option<usize>,
}

/// C `sprintf`. A conversion whose argument type does not match (undefined behaviour in C)
/// panics: a loud failure in the gate rather than a quietly different report.
pub fn format(fmt: &[u8], args: &[Arg]) -> Bytes {
    let fmt = cstr(fmt);
    let mut out = Vec::with_capacity(fmt.len() + 16);
    let mut ai = 0usize;
    let mut i = 0usize;
    let next = |ai: &mut usize| -> Arg {
        let a = *args.get(*ai).unwrap_or_else(|| panic!("printf: missing argument {} for {:?}", *ai, String::from_utf8_lossy(fmt)));
        *ai += 1;
        a
    };
    while i < fmt.len() {
        let c = fmt[i];
        if c != b'%' {
            out.push(c);
            i += 1;
            continue;
        }
        i += 1;
        if i >= fmt.len() {
            break;
        }
        if fmt[i] == b'%' {
            out.push(b'%');
            i += 1;
            continue;
        }
        let mut sp = Spec { left: false, plus: false, space: false, alt: false, zero: false, width: 0, prec: None };
        while i < fmt.len() {
            match fmt[i] {
                b'-' => sp.left = true,
                b'+' => sp.plus = true,
                b' ' => sp.space = true,
                b'#' => sp.alt = true,
                b'0' => sp.zero = true,
                _ => break,
            }
            i += 1;
        }
        if i < fmt.len() && fmt[i] == b'*' {
            let w = match next(&mut ai) {
                Arg::I(v) => v,
                Arg::U(v) => v as i64,
                a => panic!("printf: * width needs int, got {:?}", a),
            };
            if w < 0 {
                sp.left = true;
                sp.width = (-w) as usize;
            } else {
                sp.width = w as usize;
            }
            i += 1;
        } else {
            while i < fmt.len() && fmt[i].is_ascii_digit() {
                sp.width = sp.width * 10 + (fmt[i] - b'0') as usize;
                i += 1;
            }
        }
        if i < fmt.len() && fmt[i] == b'.' {
            i += 1;
            if i < fmt.len() && fmt[i] == b'*' {
                let p = match next(&mut ai) {
                    Arg::I(v) => v,
                    Arg::U(v) => v as i64,
                    a => panic!("printf: * precision needs int, got {:?}", a),
                };
                sp.prec = if p < 0 { None } else { Some(p as usize) };
                i += 1;
            } else {
                let mut p = 0usize;
                while i < fmt.len() && fmt[i].is_ascii_digit() {
                    p = p * 10 + (fmt[i] - b'0') as usize;
                    i += 1;
                }
                sp.prec = Some(p);
            }
        }
        // length modifiers: without one, C reads an int (32 bits)
        let mut long = false;
        while i < fmt.len() && matches!(fmt[i], b'h' | b'l' | b'L' | b'q' | b'j' | b'z' | b't') {
            if matches!(fmt[i], b'l' | b'L' | b'q' | b'j' | b'z' | b't') {
                long = true;
            }
            i += 1;
        }
        if i >= fmt.len() {
            break;
        }
        let conv = fmt[i];
        i += 1;
        match conv {
            b'd' | b'i' => {
                let v = match next(&mut ai) {
                    Arg::I(v) => v as i128,
                    Arg::U(v) => v as i128,
                    a => panic!("printf: %d needs an integer, got {:?} in {:?}", a, String::from_utf8_lossy(fmt)),
                };
                let v = if long { v as i64 as i128 } else { v as i32 as i128 };
                emit_int(&mut out, &sp, v, 10, false);
            }
            b'u' | b'x' | b'X' | b'o' => {
                let v = match next(&mut ai) {
                    Arg::I(v) => v as u64 as i128,
                    Arg::U(v) => v as i128,
                    a => panic!("printf: %u needs an integer, got {:?}", a),
                };
                let v = if long { v as u64 as i128 } else { v as u32 as i128 };
                let base = match conv {
                    b'o' => 8,
                    b'u' => 10,
                    _ => 16,
                };
                emit_int(&mut out, &sp, v, base, conv == b'X');
            }
            b'c' => {
                let ch = match next(&mut ai) {
                    Arg::I(v) => v as u8,
                    Arg::U(v) => v as u8,
                    a => panic!("printf: %c needs a char, got {:?}", a),
                };
                pad(&mut out, &sp, &[], &[ch], false);
            }
            b's' => {
                let s = match next(&mut ai) {
                    Arg::S(s) => cstr(s),
                    a => panic!("printf: %s needs a string, got {:?} in {:?}", a, String::from_utf8_lossy(fmt)),
                };
                let s = match sp.prec {
                    Some(p) if p < s.len() => &s[..p],
                    _ => s,
                };
                pad(&mut out, &sp, &[], s, false);
            }
            b'f' | b'F' | b'e' | b'E' | b'g' | b'G' => {
                let x = match next(&mut ai) {
                    Arg::F(x) => x,
                    a => panic!("printf: %{} needs a double, got {:?} in {:?}", conv as char, a, String::from_utf8_lossy(fmt)),
                };
                emit_float(&mut out, &sp, x, conv);
            }
            b'p' => {
                let _ = next(&mut ai);
                out.extend_from_slice(b"(nil)");
            }
            _ => {
                // an unknown conversion: glibc prints it literally
                out.push(b'%');
                out.push(conv);
            }
        }
    }
    out
}

fn pad(out: &mut Vec<u8>, sp: &Spec, sign: &[u8], body: &[u8], zero_ok: bool) {
    let len = sign.len() + body.len();
    let fill = sp.width.saturating_sub(len);
    if sp.left {
        out.extend_from_slice(sign);
        out.extend_from_slice(body);
        out.extend(std::iter::repeat_n(b' ', fill));
    } else if sp.zero && zero_ok {
        out.extend_from_slice(sign);
        out.extend(std::iter::repeat_n(b'0', fill));
        out.extend_from_slice(body);
    } else {
        out.extend(std::iter::repeat_n(b' ', fill));
        out.extend_from_slice(sign);
        out.extend_from_slice(body);
    }
}

fn emit_int(out: &mut Vec<u8>, sp: &Spec, v: i128, base: u32, upper: bool) {
    let neg = v < 0;
    let mag = v.unsigned_abs();
    let mut digits: Vec<u8> = if mag == 0 && sp.prec == Some(0) {
        Vec::new()
    } else {
        let s = match base {
            8 => format!("{:o}", mag),
            16 if upper => format!("{:X}", mag),
            16 => format!("{:x}", mag),
            _ => format!("{}", mag),
        };
        s.into_bytes()
    };
    if let Some(p) = sp.prec {
        while digits.len() < p {
            digits.insert(0, b'0');
        }
    }
    let sign: &[u8] = if neg {
        b"-"
    } else if sp.plus && base == 10 {
        b"+"
    } else if sp.space && base == 10 {
        b" "
    } else {
        b""
    };
    // the 0 flag is ignored when a precision is given
    pad(out, sp, sign, &digits, sp.prec.is_none());
}

fn emit_float(out: &mut Vec<u8>, sp: &Spec, x: f64, conv: u8) {
    let upper = conv.is_ascii_uppercase();
    let neg = x.is_sign_negative();
    let sign: &[u8] = if neg {
        b"-"
    } else if sp.plus {
        b"+"
    } else if sp.space {
        b" "
    } else {
        b""
    };
    if !x.is_finite() {
        let body: &[u8] = match (x.is_nan(), upper) {
            (true, false) => b"nan",
            (true, true) => b"NAN",
            (false, false) => b"inf",
            (false, true) => b"INF",
        };
        pad(out, sp, sign, body, false);
        return;
    }
    let a = x.abs();
    let prec = sp.prec.unwrap_or(6);
    let body: Vec<u8> = match conv {
        b'f' | b'F' => fixed(a, prec, sp.alt),
        b'e' | b'E' => expo(a, prec, sp.alt, upper),
        _ => general(a, prec, sp.alt, upper),
    };
    pad(out, sp, sign, &body, true);
}

/// `%.{prec}f` of a non-negative finite value. Rust's fixed-precision formatting is exact
/// and rounds ties to even, as glibc does.
fn fixed(a: f64, prec: usize, alt: bool) -> Vec<u8> {
    let mut s = format!("{:.*}", prec, a).into_bytes();
    if alt && prec == 0 {
        s.push(b'.');
    }
    s
}

/// Mantissa digits and decimal exponent of `%.{prec}e`.
fn expo_parts(a: f64, prec: usize) -> (Vec<u8>, i32) {
    let s = format!("{:.*e}", prec, a);
    let k = s.find('e').unwrap();
    let e: i32 = s[k + 1..].parse().unwrap();
    (s[..k].as_bytes().to_vec(), e)
}

fn exp_suffix(e: i32, upper: bool) -> Vec<u8> {
    let mut v = vec![if upper { b'E' } else { b'e' }];
    v.push(if e < 0 { b'-' } else { b'+' });
    let m = e.unsigned_abs();
    if m < 10 {
        v.push(b'0');
    }
    v.extend_from_slice(m.to_string().as_bytes());
    v
}

fn expo(a: f64, prec: usize, alt: bool, upper: bool) -> Vec<u8> {
    let (mut m, e) = expo_parts(a, prec);
    if alt && prec == 0 {
        m.push(b'.');
    }
    m.extend(exp_suffix(e, upper));
    m
}

fn general(a: f64, prec: usize, alt: bool, upper: bool) -> Vec<u8> {
    let p = if prec == 0 { 1 } else { prec };
    let x_exp = if a == 0.0 { 0 } else { expo_parts(a, p - 1).1 };
    let mut body;
    if (x_exp as i64) < p as i64 && x_exp >= -4 {
        body = fixed(a, (p as i64 - 1 - x_exp as i64) as usize, alt);
        if !alt {
            strip_zeros(&mut body);
        }
    } else {
        let (mut m, e) = expo_parts(a, p - 1);
        if !alt {
            strip_zeros(&mut m);
        } else if p - 1 == 0 {
            m.push(b'.');
        }
        m.extend(exp_suffix(e, upper));
        body = m;
    }
    body
}

fn strip_zeros(s: &mut Vec<u8>) {
    if s.contains(&b'.') {
        while s.last() == Some(&b'0') {
            s.pop();
        }
        if s.last() == Some(&b'.') {
            s.pop();
        }
    }
}

/// `std::to_string(double)`: `%f`.
pub fn to_string_f64(v: f64) -> Bytes {
    format(b"%f", &[Arg::F(v)])
}

// ============================================================================
// CFile: a C output stream
// ============================================================================

/// A FILE* opened for writing. Buffered (std::io::BufWriter); `tell` is ftell (the stream
/// is opened truncated and only written, so the position is the count of bytes written).
pub struct CFile {
    w: std::io::BufWriter<std::fs::File>,
    pos: u64,
}

impl CFile {
    /// fopen(path, mode) for mode "w" (write) or "w+b" (read and write), creating or
    /// truncating; None where fopen returns null.
    pub fn create(path: &[u8], read_too: bool) -> Option<CFile> {
        let p = crate::core::file_io::path_from_bytes(crate::core::cstdlib::c_str(path));
        let f = std::fs::OpenOptions::new().read(read_too).write(true).create(true).truncate(true).open(p).ok()?;
        Some(CFile { w: std::io::BufWriter::with_capacity(1 << 16, f), pos: 0 })
    }
    pub fn write_bytes(&mut self, b: &[u8]) {
        let _ = self.w.write_all(b);
        self.pos += b.len() as u64;
    }
    pub fn tell(&self) -> i64 {
        self.pos as i64
    }
    pub fn flush(&mut self) {
        let _ = self.w.flush();
    }
}
