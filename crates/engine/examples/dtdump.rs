// Rust side of gate/diff/dtdump.cpp.
use std::fmt::Write as _;
use swmm6_rust::core::datetime::*;
fn main() {
    let mut o = String::new();
    let hx = |o: &mut String, v: f64| { write!(o, "{:016x} ", v.to_bits()).unwrap(); };
    let mut y = -2;
    while y <= 10002 {
        for m in 0..=13 { let mut d = -1; while d <= 32 { hx(&mut o, encode_date(y, m, d)); d += 3; } }
        y += if y < 1905 || y > 2105 { 97 } else { 1 };
    }
    o.push('\n');
    for h in -1..=25 { let mut mi = -1; while mi <= 61 { let mut s = -1; while s <= 61 { hx(&mut o, encode_time(h, mi, s)); s += 13; } mi += 7; } }
    o.push('\n');
    let mut st: u64 = 88172645463325252;
    for _ in 0..200000 {
        st ^= st << 13; st ^= st >> 7; st ^= st << 17;
        let x = (st % 4000000000) as f64 / 1000.0 - 700000.0;
        let (y, m, d) = decode_date(x); let (h, mi, s) = decode_time(x);
        write!(o, "{} {} {} {} {} {} {} {} ", y, m, d, h, mi, s, month_of_year(x), day_of_year(x)).unwrap();
        hx(&mut o, add_seconds(x, (st % 100000) as f64 / 7.0));
        writeln!(o, "{}", time_diff(x, x - (st % 977) as f64 / 13.0)).unwrap();
    }
    print!("{o}");
}
