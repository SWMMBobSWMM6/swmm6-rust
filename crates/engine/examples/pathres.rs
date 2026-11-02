// Rust side of gate/diff/pathres.cpp: core::path_resolver and core::fs_path on the same lines.
use std::io::{BufRead, Write};
use swmm6_rust::core::fs_path::FsPath;
use swmm6_rust::core::path_resolver as io;
fn p(o: &mut Vec<u8>, l: &str, s: &[u8]) { write!(o, "{} {}:", l, s.len()).unwrap(); o.extend_from_slice(s); o.push(b'\n'); }
fn main() {
    let stdin = std::io::stdin();
    let mut o = Vec::new();
    for line in stdin.lock().split(b'\n') {
        let line = line.unwrap();
        let (a, b) = match line.iter().position(|&c| c == b'\t') { Some(t) => (line[..t].to_vec(), line[t + 1..].to_vec()), None => (line.clone(), Vec::new()) };
        p(&mut o, "parentDir", &io::parent_dir(&a));
        p(&mut o, "normalise", &io::normalise_separators(&a));
        p(&mut o, "isAbs", if io::is_absolute_path(&a) { b"1" } else { b"0" });
        p(&mut o, "resolve", &io::resolve_relative(&a, &b));
        let r = io::make_relative(&a, &b, 16);
        p(&mut o, "mkrel.path", &r.path); p(&mut o, "mkrel.class", (r.classification as i32).to_string().as_bytes());
        p(&mut o, "mkrel.up", r.up_levels.to_string().as_bytes()); p(&mut o, "mkrel.warn", &r.warning);
        let r2 = io::make_relative(&a, &b, 1);
        p(&mut o, "mkrel1.path", &r2.path); p(&mut o, "mkrel1.warn", &r2.warning);
        let pa = FsPath::new(&a); let pb = FsPath::new(&b);
        p(&mut o, "norm", &pa.lexically_normal().generic_string());
        p(&mut o, "normnative", pa.lexically_normal().native());
        p(&mut o, "gen", &pa.generic_string());
        p(&mut o, "join", pa.join(&pb).native());
        p(&mut o, "joinnorm", pa.join(&pb).lexically_normal().native());
        p(&mut o, "rel", pa.lexically_relative(&pb).native());
        p(&mut o, "prox", pa.lexically_proximate(&pb).native());
        p(&mut o, "parent", pa.parent_path().native());
        let mut rf = pa.clone(); rf.remove_filename(); p(&mut o, "rmfn", rf.native());
        p(&mut o, "fname", pa.filename().native());
        p(&mut o, "hasfn", if pa.has_filename() { b"1" } else { b"0" });
        p(&mut o, "hasrel", if pa.has_relative_path() { b"1" } else { b"0" });
    }
    std::io::stdout().write_all(&o).unwrap();
}
