// Context loader round trip: read a dump as gate/diff/ctxopen (C++) or examples/ctxopen
// (Rust) prints it, load it into a fresh SwmmEngine (context and the engine members the
// dump covers), and dump that again. The output must equal the input byte for byte.
// Build with --features gate. Usage: ctxload dump.txt
#[path = "common/ctxdump.rs"]
mod ctxdump;
#[path = "common/ctxload.rs"]
mod ctxload;
use std::io::Write;
use swmm6_rust::core::engine::SwmmEngine;
fn main() {
    let args: Vec<String> = std::env::args().collect();
    let text = std::fs::read(&args[1]).expect("ctxload: cannot read the dump");
    let mut r = ctxload::Rd::new(&text);
    let mut out = std::io::stdout().lock();
    for h in ["open_rc", "init_rc"] {
        if let Some(line) = r.header_line(h) {
            out.write_all(line).unwrap();
        }
    }
    let mut e = SwmmEngine::new();
    ctxload::load(&mut r, e.context_mut());
    ctxload::load_engine(&mut r, &mut e);
    r.finish();
    out.write_all(&ctxdump::dump(e.context())).unwrap();
    out.write_all(&ctxdump::engdump(&e)).unwrap();
}
