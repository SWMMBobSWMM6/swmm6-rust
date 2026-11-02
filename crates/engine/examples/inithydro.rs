// Rust side of gate/diff/inithydro.cpp: create the translated engine, open the deck as the
// harness does (report and output "/dev/null", no input plugin library), print open_rc,
// run init_hydrology when the open succeeded, then dump the context, the engine members
// open() sets and the hydrology members init_hydrology fills.
// Usage: inithydro deck.inp [open|lenient|rdii]   (run from the deck's directory, as
// gate/ctxgate.sh does); "rdii" first runs rdii_.init, as initHydraulics does.
#[path = "common/ctxdump.rs"]
mod ctxdump;
#[path = "common/hydrodump.rs"]
mod hydrodump;
use std::io::Write;
use swmm6_rust::core::engine::SwmmEngine;
fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mode = if args.len() > 2 { args[2].as_str() } else { "open" };
    let mut e = SwmmEngine::new();
    if mode == "lenient" {
        e.set_lenient_open(true);
    }
    let rc = e.open(Some(args[1].as_bytes()), Some(b"/dev/null"), Some(b"/dev/null"), None, false);
    let mut out = std::io::stdout().lock();
    writeln!(out, "open_rc {}", rc).unwrap();
    if rc == 0 {
        if mode == "rdii" {
            e.init_rdii_for_gate();
            writeln!(out, "initrdii").unwrap();
        }
        e.init_hydrology_for_gate();
        writeln!(out, "inithydro").unwrap();
    }
    out.write_all(&ctxdump::dump(e.context())).unwrap();
    out.write_all(&ctxdump::engdump(&e)).unwrap();
    out.write_all(&hydrodump::hydrodump(&mut e)).unwrap();
}
