// Rust side of gate/diff/printfg.c: core::cstdio::printf_g on the same doubles.
use std::io::BufRead;
use swmm6_rust::core::cstdio::printf_g;
fn main() {
    for line in std::io::stdin().lock().lines() {
        let u = u64::from_str_radix(line.unwrap().trim(), 16).unwrap();
        let d = f64::from_bits(u);
        let s = |p| String::from_utf8(printf_g(d, p)).unwrap();
        println!("{}|{}|{}|{}|{}", s(10), s(3), s(1), s(6), s(17));
    }
}
