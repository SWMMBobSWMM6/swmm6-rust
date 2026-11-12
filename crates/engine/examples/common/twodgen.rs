// Shared input generator and printer for the 2D kernel differentials (examples/twod_*.rs), the
// twin of gate/diff/twodgen.hpp: the same xorshift64* stream and the same exact operations, so
// both sides see identical input bits. Every output line carries its inputs and results as hex
// bits.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
#![allow(dead_code)]
use std::fmt::Write as _;
use std::io::Write;

/// An opaque copy (the twin of the C++ side's asm barrier): keeps rustc from constant-folding
/// a kernel on grid inputs.
#[inline(always)]
pub fn opq(x: f64) -> f64 {
    std::hint::black_box(x)
}

pub const INF: f64 = f64::INFINITY;
pub const QNAN: f64 = f64::from_bits(0x7ff8000000000000);
pub const NQNAN: f64 = f64::from_bits(0xfff8000000000000);
pub const DMIN: f64 = f64::from_bits(1);
pub const NMIN: f64 = f64::MIN_POSITIVE;
pub const DMAX: f64 = f64::MAX;

/// Any NaN becomes +NaN (the twin of the C++ side's canon(): every derived input goes through it,
/// so both sides see the same input bits even where two NaNs of opposite sign meet).
pub fn canon(x: f64) -> f64 {
    if x.is_nan() {
        QNAN
    } else {
        x
    }
}

pub struct Gen(pub u64);
impl Gen {
    pub fn nx(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545f4914f6cdd1d)
    }
    pub fn below(&mut self, n: u64) -> u64 {
        self.nx() % n
    }
    pub fn u01(&mut self) -> f64 {
        (self.nx() >> 11) as f64 * f64::from_bits((1023u64 - 53) << 52)
    }
    pub fn uni(&mut self, lo: f64, hi: f64) -> f64 {
        canon(lo + (hi - lo) * self.u01())
    }
    pub fn ebits(&mut self, lo: i32, hi: i32) -> f64 {
        let r = self.nx();
        let e = lo + (r % (hi - lo + 1) as u64) as i32;
        f64::from_bits(((r >> 63) << 63) | ((e as u64) << 52) | (self.nx() & ((1u64 << 52) - 1)))
    }
    pub fn pos_ebits(&mut self, lo: i32, hi: i32) -> f64 {
        self.ebits(lo, hi).abs()
    }
    pub fn pick(&mut self, v: &[f64]) -> f64 {
        v[self.below(v.len() as u64) as usize]
    }
    pub fn nudge(&mut self, mut x: f64) -> f64 {
        let k = self.below(9) as i32 - 4;
        for _ in 0..k.max(0) {
            x = x.next_up();
        }
        for _ in 0..(-k).max(0) {
            x = x.next_down();
        }
        x
    }
    pub fn nonfinite(&mut self) -> f64 {
        match self.below(4) {
            0 => INF,
            1 => -INF,
            2 => QNAN,
            _ => NQNAN,
        }
    }
    pub fn draw(&mut self, sp: &[f64], lo: f64, hi: f64) -> f64 {
        if self.below(64) == 0 {
            return self.nonfinite();
        }
        match self.below(10) {
            0..=2 => self.pick(sp),
            3 | 4 => {
                let v = self.pick(sp);
                self.nudge(v)
            }
            5..=7 => self.uni(lo, hi),
            8 => self.ebits(1023 - 45, 1023 + 12),
            _ => self.ebits(1, 2046),
        }
    }
}

pub struct Out {
    pub buf: String,
    pub lines: i64,
}
impl Default for Out {
    fn default() -> Self {
        Self::new()
    }
}
impl Out {
    pub fn new() -> Self {
        Out { buf: String::new(), lines: 0 }
    }
    pub fn tag(&mut self, t: &str) {
        self.buf.push_str(t);
        self.buf.push(' ');
    }
    pub fn hx(&mut self, v: f64) {
        let _ = write!(self.buf, "{:016x} ", v.to_bits());
    }
    pub fn hi(&mut self, v: i64) {
        let _ = write!(self.buf, "{} ", v);
    }
    /// inputs | results
    pub fn sep(&mut self) {
        self.buf.push_str("| ");
    }
    pub fn nl(&mut self) {
        self.buf.push('\n');
        self.lines += 1;
        if self.buf.len() > (1 << 22) {
            self.flush();
        }
    }
    pub fn flush(&mut self) {
        let _ = std::io::stdout().lock().write_all(self.buf.as_bytes());
        self.buf.clear();
    }
}

/// The cartesian product of index ranges: calls f(idx) for every combination (first index
/// fastest), as the C++ side's twodgen::grid.
pub fn grid(sizes: &[usize], mut f: impl FnMut(&[usize])) {
    if sizes.iter().any(|&s| s == 0) {
        return;
    }
    let mut idx = vec![0usize; sizes.len()];
    loop {
        f(&idx);
        let mut k = 0;
        while k < sizes.len() {
            idx[k] += 1;
            if idx[k] < sizes[k] {
                break;
            }
            idx[k] = 0;
            k += 1;
        }
        if k == sizes.len() {
            return;
        }
    }
}
