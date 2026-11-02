// Rust side of gate/diff/tabdump.cpp.
use std::io::{Read, Write};
use swmm6_rust::data::table::*;
fn main() {
    let mut b = Vec::new();
    std::fs::File::open(std::env::args().nth(1).unwrap()).unwrap().read_to_end(&mut b).unwrap();
    let mut o = 0usize;
    let i32_ = |o: &mut usize| { let v = i32::from_le_bytes(b[*o..*o+4].try_into().unwrap()); *o += 4; v };
    let nt = i32_(&mut o);
    let mut out = String::new();
    let hx = |s: &mut String, v: f64| { s.push_str(&format!("{:016x} ", v.to_bits())); };
    for _ in 0..nt {
        let n = i32_(&mut o) as usize;
        let rd = |o: &mut usize, k: usize| -> Vec<f64> { let v = (0..k).map(|j| f64::from_le_bytes(b[*o+8*j..*o+8*j+8].try_into().unwrap())).collect(); *o += 8*k; v };
        let x = rd(&mut o, n); let y = rd(&mut o, n);
        let nq = i32_(&mut o) as usize; let q = rd(&mut o, nq);
        let mk = || { let mut t = Table::new(b"T", TableType::CurveStorage); t.x = x.clone(); t.y = y.clone(); t };
        let (tb, mut a, mut bb, mut c) = (mk(), mk(), mk(), mk());
        for &xq in &q {
            hx(&mut out, table_lookup_cursor(&mut a, xq)); hx(&mut out, table_tseries_lookup_cursor(&mut bb, xq)); hx(&mut out, table_step_cursor(&mut c, xq));
            hx(&mut out, table_get_storage_volume(&tb, xq)); hx(&mut out, table_get_storage_depth(&tb, xq)); hx(&mut out, table_inverse_lookup(&tb, xq));
            hx(&mut out, table_lookup_ex(&tb, xq)); hx(&mut out, table_interval_lookup(&tb, xq)); hx(&mut out, table_get_slope(&tb, xq)); hx(&mut out, table_get_max_y(&tb));
            out.push_str(&format!("{} {} {}\n", a.cursor.index, bb.cursor.index, c.cursor.index));
        }
    }
    std::io::stdout().write_all(out.as_bytes()).unwrap();
}
