// Rust side of gate/diff/ctxopen.cpp: create the translated engine, open the deck as the
// harness does (report and output "/dev/null", no input plugin library), print open_rc,
// then dump the context and the engine members open() sets.
// Usage: ctxopen deck.inp [open|lenient|metadata]
#[path = "common/ctxdump.rs"]
mod ctxdump;
use std::io::Write;
use swmm6_rust::core::engine::SwmmEngine;
fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mode = if args.len() > 2 { args[2].as_str() } else { "open" };
    let mut e = SwmmEngine::new();
    if mode == "lenient" {
        e.set_lenient_open(true);
    }
    let rc = if mode == "metadata" {
        // swmm_engine_open_metadata refuses a null or empty path with SWMM_ERR_BADPARAM
        if args[1].is_empty() { 9 } else { e.open(Some(args[1].as_bytes()), None, None, None, true) }
    } else {
        e.open(Some(args[1].as_bytes()), Some(b"/dev/null"), Some(b"/dev/null"), None, false)
    };
    let mut out = std::io::stdout().lock();
    writeln!(out, "open_rc {}", rc).unwrap();
    if rc == 0 && mode == "init" {
        panic!("initialize is not translated yet");
    }
    out.write_all(&ctxdump::dump(e.context())).unwrap();
    out.write_all(&ctxdump::engdump(&e)).unwrap();
}
