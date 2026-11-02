// Rust side of gate/diff/inithyd.cpp: open the deck as the harness does (report and
// output "/dev/null", no input plugin library), run the translated initHydraulics, then
// dump the context, the engine members open() sets and the state initHydraulics builds.
// Where upstream's initHydraulics would end the process (std::terminate from an exception
// escaping the noexcept function, finding U1), this harness aborts too, after printing
// nothing, as the reference harness's buffered stdout is lost on abort.
// Usage: inithyd deck.inp
#[path = "common/ctxdump.rs"]
mod ctxdump;
#[path = "common/rtdump.rs"]
mod rtdump;
use std::io::Write;
use swmm6_rust::core::engine::SwmmEngine;
fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        std::process::exit(2);
    }
    let mut e = SwmmEngine::new();
    let rc = e.open(Some(args[1].as_bytes()), Some(b"/dev/null"), Some(b"/dev/null"), None, false);
    let mut out: Vec<u8> = Vec::new();
    writeln!(out, "open_rc {}", rc).unwrap();
    if rc == 0 {
        if let Err(t) = e.gate_init_hydraulics() {
            eprintln!("terminate called after throwing an instance of '{}'\n  what():  {}", t.exception, t.what);
            std::process::abort();
        }
    }
    out.extend_from_slice(&ctxdump::dump(e.context()));
    out.extend_from_slice(&ctxdump::engdump(&e));
    out.extend_from_slice(&rtdump::rtdump(&e));
    std::io::stdout().lock().write_all(&out).unwrap();
}
