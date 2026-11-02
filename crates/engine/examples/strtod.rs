// Rust side of gate/diff/strtod.c.
use std::io::Write;
use swmm6_rust::core::cstdlib::{atoi, strtod};
fn main() {
    let data = std::fs::read(std::env::args().nth(1).unwrap()).unwrap();
    let mut lines: Vec<&[u8]> = data.split(|&c| c == b'\n').collect();
    if data.last() == Some(&b'\n') { lines.pop(); }
    let mut o = Vec::new();
    for l in lines { let (v, used) = strtod(l); writeln!(o, "{} {:016x} {}", used, v.to_bits(), atoi(l)).unwrap(); }
    std::io::stdout().write_all(&o).unwrap();
}
