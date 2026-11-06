// The two C library hazards of the FV step against g++ and glibc (twin of
// gate/diff/ext80dot.cpp): core::ext80::dot_ext80 (NativeGmres.cpp's long double `dot`) and
// core::cmath::glibc_cbrt (std::cbrt), over the same xorshift64* inputs, as hex bits.
use std::io::Write;
use swmm6_rust::core::cmath::glibc_cbrt;
use swmm6_rust::core::ext80::dot_ext80;

struct Rng(u64);
impl Rng {
    fn nx(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545f4914f6cdd1d)
    }
    fn mk(&mut self, lo: i32, hi: i32) -> f64 {
        let r = self.nx();
        let e = lo + (r % (hi - lo + 1) as u64) as i32;
        let b = ((r >> 63) << 63) | ((e as u64) << 52) | (self.nx() & ((1u64 << 52) - 1));
        f64::from_bits(b)
    }
}

fn main() {
    let mut g = Rng(0x9e3779b97f4a7c15);
    let mut o = Vec::new();
    for _t in 0..300000 {
        let kind = g.nx() % 6;
        let n = 1 + (g.nx() % 64) as usize;
        let mut a = vec![0.0f64; n];
        let mut b = vec![0.0f64; n];
        for i in 0..n {
            match kind {
                0 => {
                    a[i] = g.mk(1023 - 30, 1023 + 30);
                    b[i] = g.mk(1023 - 30, 1023 + 30);
                }
                1 => {
                    a[i] = g.mk(1, 2046);
                    b[i] = g.mk(1023 - 200, 1023 + 200);
                }
                2 => {
                    a[i] = g.mk(1023 - 5, 1023 + 5);
                    b[i] = if i % 2 == 1 { -a[i - 1] * (1.0 + (g.nx() % 8) as f64 * 2f64.powi(-50)) } else { g.mk(1023 - 5, 1023 + 5) };
                }
                3 => {
                    a[i] = f64::from_bits(g.nx() & ((1u64 << 52) - 1));
                    b[i] = g.mk(1023 + 500, 1023 + 1000);
                }
                4 => {
                    a[i] = if g.nx() % 4 == 0 { 0.0 } else { g.mk(1023 - 60, 1023 + 60) };
                    b[i] = if g.nx() % 4 == 0 { -0.0 } else { g.mk(1023 - 60, 1023 + 60) };
                }
                _ => {
                    a[i] = g.mk(1, 400);
                    b[i] = g.mk(1, 400);
                }
            }
        }
        writeln!(o, "{:016x}", dot_ext80(&a, &b).to_bits()).unwrap();
    }
    for t in 0..300000 {
        let mut x = if t % 3 == 0 { g.mk(1, 2046) } else { g.mk(1023 - 40, 1023 + 40) };
        if t % 97 == 0 {
            x = f64::from_bits(g.nx() & ((1u64 << 52) - 1));
        }
        writeln!(o, "{:016x}", glibc_cbrt(x).to_bits()).unwrap();
    }
    std::io::stdout().write_all(&o).unwrap();
}
