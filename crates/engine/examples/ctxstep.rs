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
// Mode "helper K NAME SIDE" (see gate/diff/ctxstep.cpp): after the translation's own open,
// initialize and start, load upstream's state after K steps (SIDE.ctx, SIDE.eng) and call the
// translation of the named step() helper on it.
// Build: cargo build --release --features gate --example ctxstep
#[path = "common/ctxdump.rs"]
mod ctxdump;
#[path = "common/ctxload.rs"]
mod ctxload;
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

/// SIDE.eng as gate/diff/ctxstep's helper mode writes it: the "helper_at k rc r" line, k,
/// and the members.
fn read_side(path: &str) -> (String, i32, swmm6_rust::core::engine_step::GateStepSide) {
    let text = std::fs::read_to_string(path).expect("ctxstep: cannot read the side file");
    let first = text.lines().next().unwrap_or("").to_string();
    let mut it = text.split_ascii_whitespace();
    let mut next = || it.next().expect("side: short").to_string();
    assert_eq!(next(), "helper_at");
    let k: i32 = next().parse().unwrap();
    assert_eq!(next(), "rc");
    let _rc = next();
    let int = |name: &str, next: &mut dyn FnMut() -> String| -> i64 {
        assert_eq!(next(), name);
        next().parse().unwrap()
    };
    let bits = |x: String| f64::from_bits(u64::from_str_radix(&x, 16).unwrap());
    let vec = |name: &str, next: &mut dyn FnMut() -> String| -> Vec<f64> {
        assert_eq!(next(), name);
        let n: usize = next().parse().unwrap();
        (0..n).map(|_| bits(next())).collect()
    };
    let mut s = swmm6_rust::core::engine_step::GateStepSide::default();
    s.next_event_ = int("next_event_", &mut next) as i32;
    s.fv_accepted = int("fv_accepted", &mut next) != 0;
    for (k, name) in ["area", "imperv_pct", "frac_imperv0", "frac_imperv1", "depth_imperv0", "depth_imperv1", "depth_perv"].iter().enumerate() {
        s.runoff[k] = vec(name, &mut next);
    }
    s.snow_enabled = vec("snow_enabled", &mut next).iter().map(|&x| x as u8).collect();
    s.snow_wsnow = vec("snow_wsnow", &mut next);
    s.snow_fw = vec("snow_fw", &mut next);
    s.snow_farea = vec("snow_fArea", &mut next);
    assert_eq!(next(), "snow_removed");
    s.snow_removed = bits(next());
    s.buildup = vec("buildup", &mut next);
    s.lid_drain_q_interp_ = vec("lid_drain_q_interp_", &mut next);
    s.wet_q_interp_ = vec("wet_q_interp_", &mut next);
    s.gw_q_interp_ = vec("gw_q_interp_", &mut next);
    s.gw_q_node_ = vec("gw_q_node_", &mut next).iter().map(|&x| x as i32).collect();
    let un = vec("lid_units_node", &mut next);
    let uq = vec("lid_units_q", &mut next);
    s.lid_drain_q_units_ = un.iter().zip(uq.iter()).map(|(&n, &q)| (n as i32, q)).collect();
    let ng = int("lid_groups", &mut next);
    for _ in 0..ng {
        let count = int("lidg_count", &mut next) as i32;
        let sd = vec("lidg_surf_depth", &mut next);
        let sm = vec("lidg_soil_moist", &mut next);
        let st = vec("lidg_stor_depth", &mut next);
        let pd = vec("lidg_pave_depth", &mut next);
        s.lid_groups.push((count, [sd, sm, st, pd]));
    }
    assert_eq!(next(), "end");
    (first, k, s)
}

const HELPERS: [&str; 9] = ["stats", "massbal", "finalstore", "finalqual", "nodetotals", "diverged", "events", "lateral", "forcings"];

/// Load upstream's state after K steps (SIDE.ctx, SIDE.eng) into `e` and run one helper;
/// the side file's first line, the helper's own lines, and whether the deck was skipped.
fn load_and_run(e: &mut SwmmEngine, name: &str, side: &str, out: &mut Vec<u8>) -> (String, bool) {
    let ctx_text = std::fs::read(format!("{}.ctx", side)).expect("ctxstep: cannot read SIDE.ctx");
    let mut r = ctxload::Rd::new(&ctx_text);
    ctxload::load(&mut r, e.context_mut());
    r.finish();
    let (first, k, side_members) = read_side(&format!("{}.eng", side));
    e.gate_install_step_side(side_members);
    let c = e.context();
    let dt = if k > 0 { (c.elapsed_ms - c.old_elapsed_ms) / 1000.0 } else { c.options.routing_step };
    let mut lines = Vec::new();
    let r = guarded(|| {
        e.gate_step_helper(name, dt).map(|b| {
            lines = b;
            0
        })
    });
    if let Err(t) = r {
        finish(Err(t), &format!("helper {}", name), out);
    }
    let skipped = lines.ends_with(b"skipped (FV)\n") || lines.ends_with(b"skipped (inlets)\n");
    out.extend_from_slice(&lines);
    (first, skipped)
}

/// The helper mode: one helper with the full context dump after it, or "all" with hashes.
fn helper_mode(e: &mut SwmmEngine, args: &[String], out: &mut Vec<u8>) {
    let name = args[4].as_str();
    let side = args[5].as_str();
    if name != "all" {
        let mut lines = Vec::new();
        let (first, skipped) = load_and_run(e, name, side, &mut lines);
        writeln!(out, "{}", first).unwrap();
        out.extend_from_slice(&lines);
        if !skipped {
            out.extend_from_slice(&ctxdump::dump(e.context()));
        }
        return;
    }
    let mut first_line = None;
    let mut body = Vec::new();
    for n in HELPERS {
        // a fresh engine through the translation's own open, initialize and start each time
        let mut f = SwmmEngine::new();
        let _ = f.open(Some(args[1].as_bytes()), Some(b"/dev/null"), Some(b"/dev/null"), None, false);
        let _ = guarded(|| f.initialize());
        let _ = guarded(|| f.start(1));
        writeln!(body, "helper {}", n).unwrap();
        let (first, skipped) = load_and_run(&mut f, n, side, &mut body);
        first_line.get_or_insert(first);
        if !skipped {
            writeln!(body, "hash {} {:016x}", n, fnv1a(&ctxdump::dump(f.context()))).unwrap();
        }
    }
    let _ = e;
    writeln!(out, "{}", first_line.unwrap_or_default()).unwrap();
    out.extend_from_slice(&body);
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
    if mode == "helper" {
        if rc == 0 {
            helper_mode(&mut e, &args, &mut out);
        }
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
