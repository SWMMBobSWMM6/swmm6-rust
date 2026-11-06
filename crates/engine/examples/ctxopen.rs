// Rust side of gate/diff/ctxopen.cpp: create the translated engine, open the deck as the
// harness does (report and output "/dev/null", no input plugin library), print open_rc,
// then dump the context and the engine members open() sets.
// Usage: ctxopen deck.inp [open|lenient|metadata|init]
// "init" also runs SwmmEngine::initialize (core/engine_initialize.rs) and prints init_rc;
// where upstream terminates inside initialize, it prints libstdc++'s terminate message and
// aborts, as the C++ harness's process does.
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
    // buffered and written at the end, as the C++ harness's stdout to a pipe is: where
    // upstream terminates, its buffered lines are lost with the process
    let mut out: Vec<u8> = Vec::new();
    writeln!(out, "open_rc {}", rc).unwrap();
    if rc == 0 && mode == "init" {
        match e.initialize() {
            Ok(irc) => writeln!(out, "init_rc {}", irc).unwrap(),
            Err(t) => {
                eprintln!("terminate called after throwing an instance of '{}'\n  what():  {}", t.exception, t.what);
                std::process::abort();
            }
        }
    }
    out.extend_from_slice(&ctxdump::dump(e.context()));
    out.extend_from_slice(&ctxdump::engdump(&e));
    std::io::stdout().lock().write_all(&out).unwrap();
}
