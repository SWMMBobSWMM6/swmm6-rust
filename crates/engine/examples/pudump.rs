// Rust side of gate/diff/pudump.cpp.
use std::io::Write;
use swmm6_rust::data::name_index::NameIndex;
use swmm6_rust::input::parse_utils::*;
fn main() {
    let data = std::fs::read(std::env::args().nth(1).unwrap()).unwrap();
    let mut lines: Vec<&[u8]> = data.split(|&c| c == b'\n').collect();
    if data.last() == Some(&b'\n') { lines.pop(); }
    let mut o: Vec<u8> = Vec::new(); let mut ix = NameIndex::new(); let mut errs = Vec::new();
    let hx = |o: &mut Vec<u8>, v: f64| { write!(o, "{:016x} ", v.to_bits()).unwrap(); };
    for t in lines {
        let mut s = -7.0; let ok = parse_double_strict(t, &mut s);
        hx(&mut o, to_double(t, -1.0)); write!(o, "{} ", ok as i32).unwrap(); hx(&mut o, s); write!(o, "{} ", to_int(t, -99)).unwrap();
        hx(&mut o, parse_date(t)); hx(&mut o, parse_time_seconds(t)); hx(&mut o, parse_time_day_fraction(t)); hx(&mut o, parse_datetime(t, t));
        writeln!(o, "{}", add_unique(&mut ix, t, &mut errs)).unwrap();
    }
    for e in errs { o.extend_from_slice(&e); o.push(b'\n'); }
    std::io::stdout().write_all(&o).unwrap();
}
