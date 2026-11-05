// Rust side of gate/diff/ctxstep.cpp, the per-step state gate: open the deck in the
// translated engine as the harness does (report and output "/dev/null", no input plugin
// library), then initialize, start(1) and step as the CLI's run loop does, printing the same
// lines and dumps as the C++ side.
// Usage: ctxstep deck.inp init|start|step [N [hash|full K]]   (see gate/diff/ctxstep.cpp)
//
// Where the translation calls code that is not translated yet (a stub file, or a branch
// marked as a stub), the engine returns Err(Terminated) with exception
// swmm6_rust::core::engine_step::STUB_EXCEPTION, or (the Router stubs, which return no
// Result) panics with "STUB: ...". The harness then prints what it has so far and a last line
// "stub <phase> <what>" (phase: init, start or step k) and exits 0, so ctxgate.sh can report
// how far the deck got. Any other Err(Terminated) is where upstream terminates: the harness
// prints libstdc++'s terminate message to stderr and aborts, as the C++ side does.
// Build: cargo build --release --features gate --example ctxstep
#[path = "common/ctxdump.rs"]
mod ctxdump;
#[path = "common/hydrodump.rs"]
mod hydrodump;
#[path = "common/lifedump.rs"]
mod lifedump;
#[path = "common/rtdump.rs"]
mod rtdump;
use std::io::Write;
use swmm6_rust::core::engine::SwmmEngine;
use swmm6_rust::core::engine_init_hydraulics::Terminated;
use swmm6_rust::core::engine_step::STUB_EXCEPTION;

fn full_dump(e: &mut SwmmEngine) -> Vec<u8> {
    let mut o = ctxdump::dump(e.context());
    o.extend_from_slice(&ctxdump::engdump(e));
    o.extend_from_slice(&rtdump::rtdump(e));
    o.extend_from_slice(&hydrodump::hydrodump(e));
    o.extend_from_slice(&ctxdump::qualdump(e));
    o.extend_from_slice(&lifedump::lifedump(e));
    o
}

fn fnv1a(b: &[u8]) -> u64 {
    let mut h: u64 = 1469598103934665603;
    for &c in b {
        h ^= c as u64;
        h = h.wrapping_mul(1099511628211);
    }
    h
}

fn dump_or_hash(e: &mut SwmmEngine, hash: bool, what: &str, out: &mut Vec<u8>) {
    let d = full_dump(e);
    if hash {
        writeln!(out, "hash {} {:016x}", what, fnv1a(&d)).unwrap();
    } else {
        out.extend_from_slice(&d);
    }
}

/// Upstream's process ends here (std::terminate); a stub ends the harness's run instead.
fn finish(r: Result<i32, Terminated>, phase: &str, out: &mut Vec<u8>) -> Option<i32> {
    match r {
        Ok(rc) => Some(rc),
        Err(t) if t.exception == STUB_EXCEPTION => {
            writeln!(out, "stub {} {}", phase, t.what).unwrap();
            std::io::stdout().lock().write_all(out).unwrap();
            std::process::exit(0);
        }
        Err(t) => {
            eprintln!("terminate called after throwing an instance of '{}'\n  what():  {}", t.exception, t.what);
            std::process::abort();
        }
    }
}

/// Run `f`, turning a "STUB: ..." panic (the Router stubs) into Err(stub).
fn guarded<F: FnOnce() -> Result<i32, Terminated>>(f: F) -> Result<i32, Terminated> {
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(f)) {
        Ok(r) => r,
        Err(p) => {
            let msg = if let Some(s) = p.downcast_ref::<&str>() {
                s.to_string()
            } else if let Some(s) = p.downcast_ref::<String>() {
                s.clone()
            } else {
                String::from("panic")
            };
            match msg.strip_prefix("STUB: ") {
                Some(w) => Err(Terminated { exception: STUB_EXCEPTION, what: Box::leak(w.to_string().into_boxed_str()) }),
                None => {
                    eprintln!("panic: {}", msg);
                    std::process::exit(101);
                }
            }
        }
    }
}

fn main() {
    // the Router stubs panic; keep the default hook's message off stderr for those
    std::panic::set_hook(Box::new(|info| {
        let s = info.to_string();
        if !s.contains("STUB: ") {
            eprintln!("{}", s);
        }
    }));
    let args: Vec<String> = std::env::args().collect();
    let mode = if args.len() > 2 { args[2].as_str() } else { "init" };
    let nsteps: i32 = if args.len() > 3 { args[3].parse().unwrap_or(1) } else { 1 };
    let how = if args.len() > 4 { args[4].as_str() } else { "" };
    let full_at: i32 = if how == "full" && args.len() > 5 { args[5].parse().unwrap_or(-1) } else { -1 };
    let hashing = how == "hash" || full_at >= 0;
    let mut out: Vec<u8> = Vec::new();
    let mut e = SwmmEngine::new();
    let mut rc = e.open(Some(args[1].as_bytes()), Some(b"/dev/null"), Some(b"/dev/null"), None, false);
    writeln!(out, "open_rc {}", rc).unwrap();
    if rc == 0 {
        let r = guarded(|| e.initialize());
        rc = finish(r, "init", &mut out).unwrap();
        writeln!(out, "init_rc {}", rc).unwrap();
    }
    if mode == "init" {
        out.extend_from_slice(&full_dump(&mut e));
        std::io::stdout().lock().write_all(&out).unwrap();
        return;
    }
    if rc == 0 {
        let r = guarded(|| e.start(1));
        rc = finish(r, "start", &mut out).unwrap();
        writeln!(out, "start_rc {}", rc).unwrap();
    }
    if mode == "start" {
        out.extend_from_slice(&full_dump(&mut e));
        std::io::stdout().lock().write_all(&out).unwrap();
        return;
    }
    dump_or_hash(&mut e, hashing && full_at != 0, "start", &mut out);
    let mut k = 1;
    while rc == 0 && k <= nsteps {
        let mut elapsed = -1.0f64;
        let r = guarded(|| e.step(&mut elapsed));
        rc = finish(r, &format!("step {}", k), &mut out).unwrap();
        writeln!(out, "step {} rc {} elapsed {:016x} state {}", k, rc, elapsed.to_bits(), e.context().state as i32).unwrap();
        dump_or_hash(&mut e, hashing && k != full_at, &k.to_string(), &mut out);
        if elapsed == 0.0 {
            break;
        }
        k += 1;
    }
    std::io::stdout().lock().write_all(&out).unwrap();
}
