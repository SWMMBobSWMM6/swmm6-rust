// Rust side of gate/diff/stgeo.cpp.
use swmm6_rust::data::node::StorageShape;
use swmm6_rust::data::storage_geometry::*;
fn shape(i: u64) -> StorageShape { [StorageShape::Tabular, StorageShape::Functional, StorageShape::Cylindrical, StorageShape::Conical, StorageShape::Paraboloid, StorageShape::Pyramidal][i as usize] }
fn main() {
    let mut st: u64 = 12345; let mut o = String::new();
    let nx = |st: &mut u64| { *st ^= *st << 13; *st ^= *st >> 7; *st ^= *st << 17; };
    for _ in 0..60000 {
        nx(&mut st); let sh = st % 6; let mut p = [0.0f64; 3];
        for j in 0..3 { nx(&mut st); let m = st % 10; p[j] = if m == 0 { 0.0 } else if m == 1 { -1.5 } else { (st % 100000) as f64 / 97.0 }; }
        let (mut a, mut b, mut c) = (-9.0, -9.0, -9.0);
        let ok = storage_shape_coeffs(shape(sh), p[0], p[1], p[2], &mut a, &mut b, &mut c);
        o += &format!("{} {} {} {} {:016x} {:016x} {:016x} {}\n", sh, storage_shape_params_valid(shape(sh), p[0], p[1], p[2]) as i32, ok as i32, storage_shape_is_geometric(shape(sh)) as i32, a.to_bits(), b.to_bits(), c.to_bits(), String::from_utf8_lossy(storage_shape_keyword(shape(sh))));
    }
    for k in ["TABULAR", "functional", "Cylindrical", "CONICAL", "PARABOLOID", "PYRAMIDAL", "x", ""] {
        let mut out = StorageShape::Tabular; let r = storage_shape_from_keyword(k.as_bytes(), &mut out);
        o += &format!("{} {} {}\n", k, r as i32, out as i32);
    }
    for v in -2..9 { o += &format!("{} ", storage_shape_is_valid_code(v) as i32); }
    o.push('\n'); print!("{o}");
}
