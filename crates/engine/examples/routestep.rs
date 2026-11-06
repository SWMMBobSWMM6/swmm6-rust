// The routing-step cross test, translation side (twin of gate/diff/routestep.cpp).
//
// Reads the stream gate/diff/routestep printed for a deck: for every
//   BEFORE c <dt bits> <dump> AFTER c <dump>
// block it opens the deck in a fresh translated engine, runs the translated init legs
// (initHydraulics, initQuality) so every solver has its structure, loads upstream's BEFORE dump
// over the context and the solver members (examples/common/ctxload.rs: load, load_rt,
// load_qual, and the stepdump members by hand), runs the translated SWMMEngine::stepRouting
// with upstream's dt, and prints "AFTER c" and the same dumps. The gate compares that text with
// upstream's AFTER blocks (gate/ctxgate.sh route).
//
// Usage: routestep deck.inp stream.txt [step|load]
//        routestep deck.inp stream|streamfresh     (upstream's output on stdin; see stream())
//   step  (default) as above
//   load  print the loaded state without stepping (the loader's round trip: the output must
//         equal upstream's BEFORE blocks)
// A translated stepRouting that reaches a path not translated yet prints "AFTER c" and
// "nottranslated <what>"; an exception upstream would throw prints "throw <what>", as the
// C++ harness does.
#[path = "common/ctxdump.rs"]
mod ctxdump;
#[path = "common/ctxload.rs"]
mod ctxload;
#[path = "common/rtdump.rs"]
mod rtdump;
use std::io::Write;
use swmm6_rust::core::context::SimulationContext;
use swmm6_rust::core::engine::SwmmEngine;
use swmm6_rust::core::engine_step_routing::{STUB_EXCEPTION, THROWN_RUNTIME_ERROR};

fn pf(o: &mut Vec<u8>, l: &str, v: f64) {
    writeln!(o, "{} {:016x}", l, v.to_bits()).unwrap();
}
fn pi(o: &mut Vec<u8>, l: &str, v: i64) {
    writeln!(o, "{} {}", l, v).unwrap();
}
fn vf(o: &mut Vec<u8>, l: &str, v: &[f64]) {
    writeln!(o, "{} n={}", l, v.len()).unwrap();
    for x in v {
        write!(o, " {:016x}", x.to_bits()).unwrap();
    }
    if !v.is_empty() {
        o.push(b'\n');
    }
}
fn vi(o: &mut Vec<u8>, l: &str, v: &[i64]) {
    writeln!(o, "{} n={}", l, v.len()).unwrap();
    for x in v {
        write!(o, " {}", x).unwrap();
    }
    if !v.is_empty() {
        o.push(b'\n');
    }
}

/// Twin of routestep.cpp's stepdump.
fn stepdump(e: &SwmmEngine) -> Vec<u8> {
    let mut o = Vec::new();
    let s = e.gate_step_state();
    pi(&mut o, "engine.last_step_steady_", s.last_step_steady_ as i64);
    pf(&mut o, "engine.new_rule_time_ms_", s.new_rule_time_ms_);
    pi(&mut o, "engine.checkpoint_enabled_", s.checkpoint_enabled_ as i64);
    pi(&mut o, "engine.ard_init_attempted_", s.ard_init_attempted_ as i64);
    pi(&mut o, "engine.rdii_iface_file_.open", s.rdii_iface_open_ as i64);
    pi(&mut o, "engine.rdii_iface_file_.writing", s.rdii_iface_writing_ as i64);
    vf(&mut o, "engine.wet_q_interp_", s.wet_q_interp_);
    pi(&mut o, "engine.lid_drain_q_units_.n", s.lid_drain_q_units_.len() as i64);
    for u in s.lid_drain_q_units_.iter() {
        pi(&mut o, "engine.lid_drain_q_units_[].first", u.0 as i64);
        pf(&mut o, "engine.lid_drain_q_units_[].second", u.1);
    }
    vf(&mut o, "engine.gw_q_interp_", s.gw_q_interp_);
    vi(&mut o, "engine.gw_q_node_", &s.gw_q_node_.iter().map(|x| *x as i64).collect::<Vec<i64>>());
    o
}

fn load_stepdump(r: &mut ctxload::Rd, e: &mut SwmmEngine) {
    let s = e.gate_step_state_mut();
    *s.last_step_steady_ = r.pi("engine.last_step_steady_") != 0;
    *s.new_rule_time_ms_ = r.pf("engine.new_rule_time_ms_");
    *s.checkpoint_enabled_ = r.pi("engine.checkpoint_enabled_") != 0;
    *s.ard_init_attempted_ = r.pi("engine.ard_init_attempted_") != 0;
    *s.rdii_iface_open_ = r.pi("engine.rdii_iface_file_.open") != 0;
    *s.rdii_iface_writing_ = r.pi("engine.rdii_iface_file_.writing") != 0;
    *s.wet_q_interp_ = r.vf("engine.wet_q_interp_");
    let n = r.count("engine.lid_drain_q_units_.n");
    s.lid_drain_q_units_.clear();
    for _ in 0..n {
        let a = r.pi("engine.lid_drain_q_units_[].first") as i32;
        let b = r.pf("engine.lid_drain_q_units_[].second");
        s.lid_drain_q_units_.push((a, b));
    }
    *s.gw_q_interp_ = r.vf("engine.gw_q_interp_");
    *s.gw_q_node_ = r.vi("engine.gw_q_node_").into_iter().map(|x| x as i32).collect();
}

fn dumpall(e: &SwmmEngine, out: &mut Vec<u8>) {
    out.extend_from_slice(&ctxdump::dump(e.context()));
    out.extend_from_slice(&rtdump::rtdump(e));
    out.extend_from_slice(&ctxdump::qualdump(e));
    out.extend_from_slice(&stepdump(e));
}

/// A fresh translated engine for the deck: open, then the init legs that build the solvers'
/// structure (the loaded dump then overwrites every member it covers).
fn build_engine(deck: &str) -> SwmmEngine {
    let mut e = SwmmEngine::new();
    let rc = e.open(Some(deck.as_bytes()), Some(b"/dev/null"), Some(b"/dev/null"), None, false);
    if rc != 0 {
        panic!("routestep: the translated open failed ({})", rc);
    }
    if e.gate_init_hydraulics().is_err() || e.gate_init_quality().is_err() {
        panic!("routestep: a translated init leg terminated");
    }
    e
}

/// Load one BEFORE dump into `e` (context first reset).
fn load_into(e: &mut SwmmEngine, body: &[u8]) {
    *e.context_mut() = SimulationContext::default();
    let mut r = ctxload::Rd::new(body);
    ctxload::load(&mut r, e.context_mut());
    ctxload::load_rt(&mut r, e);
    ctxload::load_qual(&mut r, e);
    load_stepdump(&mut r, e);
    r.finish();
}

/// Run the translated stepRouting and dump the result (or the line upstream's harness prints
/// for a throw).
fn step_and_dump(e: &mut SwmmEngine, dt: f64, out: &mut Vec<u8>) {
    match e.step_routing(dt) {
        Ok(()) => dumpall(e, out),
        Err(t) if t.exception == THROWN_RUNTIME_ERROR => {
            out.extend_from_slice(format!("throw {}\n", t.what).as_bytes());
        }
        Err(t) if t.exception == STUB_EXCEPTION => {
            out.extend_from_slice(format!("nottranslated {}\n", t.what).as_bytes());
        }
        Err(t) => {
            out.extend_from_slice(format!("terminate {} {}\n", t.exception, t.what).as_bytes());
        }
    }
}

/// The first differing line of two dumps: (label, element index or -1). A vector line's label
/// is the header line above it.
fn first_diff(a: &[u8], b: &[u8]) -> Option<(String, i64)> {
    let la: Vec<&[u8]> = a.split(|&c| c == b'\n').collect();
    let lb: Vec<&[u8]> = b.split(|&c| c == b'\n').collect();
    let mut last: &[u8] = b"?";
    for i in 0..la.len().max(lb.len()) {
        let x: &[u8] = if i < la.len() { la[i] } else { b"<end>" };
        let y: &[u8] = if i < lb.len() { lb[i] } else { b"<end>" };
        if !x.starts_with(b" ") {
            last = x.split(|&c| c == b' ').next().unwrap_or(b"?");
        }
        if x != y {
            if x.starts_with(b" ") && y.starts_with(b" ") {
                let xs: Vec<&[u8]> = x.split(|&c| c == b' ').filter(|t| !t.is_empty()).collect();
                let ys: Vec<&[u8]> = y.split(|&c| c == b' ').filter(|t| !t.is_empty()).collect();
                let k = (0..xs.len().max(ys.len())).find(|&k| k >= xs.len() || k >= ys.len() || xs[k] != ys[k]).map(|k| k as i64).unwrap_or(-1);
                return Some((String::from_utf8_lossy(last).into_owned(), k));
            }
            return Some((String::from_utf8_lossy(last).into_owned(), -1));
        }
    }
    None
}

fn finish(msg: String) -> ! {
    println!("{}", msg);
    std::process::exit(0);
}

/// `stream` mode: read gate/diff/routestep's output from stdin call by call, load each BEFORE
/// block, check the loader's round trip, step, compare with the AFTER block, and print ONE
/// line: identical N | loaddiff c label k | differ c label k | nottranslated c what |
/// nostart rc | nosteps. Nothing is written to disk; the run stops at the first difference.
/// One translated engine is built and reloaded for every call ("fresh" builds one per call).
fn stream(deck: &str, fresh: bool) {
    use std::io::BufRead;
    let stdin = std::io::stdin();
    let mut rd = std::io::BufReader::with_capacity(1 << 20, stdin.lock());
    let mut line: Vec<u8> = Vec::new();
    let mut before: Vec<u8> = Vec::new();
    let mut after: Vec<u8> = Vec::new();
    let mut state = 0; // 0 none, 1 in BEFORE, 2 in AFTER
    let mut call = 0i64;
    let mut dt = 0.0f64;
    let mut n_ok = 0i64;
    let mut eng: Option<SwmmEngine> = None;
    let mut ours: Vec<u8> = Vec::new();
    loop {
        line.clear();
        let n = rd.read_until(b'\n', &mut line).expect("routestep: read error");
        let eof = n == 0;
        let is_before = line.starts_with(b"BEFORE ");
        let is_end = line.starts_with(b"END ") || line.starts_with(b"NOSTART");
        if state == 2 && (eof || is_before || is_end) {
            if fresh || eng.is_none() {
                eng = Some(build_engine(deck));
            }
            let e = eng.as_mut().unwrap();
            load_into(e, &before);
            ours.clear();
            dumpall(e, &mut ours);
            if ours != before {
                let (lab, k) = first_diff(&before, &ours).unwrap_or(("?".into(), -1));
                finish(format!("loaddiff {} {} {}", call, lab, k));
            }
            ours.clear();
            step_and_dump(e, dt, &mut ours);
            if ours.starts_with(b"nottranslated ") {
                let w = String::from_utf8_lossy(&ours[14..ours.len() - 1]).into_owned();
                finish(format!("nottranslated {} {}", call, w));
            }
            if ours != after {
                let (lab, k) = first_diff(&after, &ours).unwrap_or(("?".into(), -1));
                finish(format!("differ {} {} {}", call, lab, k));
            }
            n_ok += 1;
            state = 0;
        }
        if eof {
            break;
        }
        if line.starts_with(b"NOSTART") {
            finish(format!("nostart {}", String::from_utf8_lossy(&line[8..]).trim_end()));
        }
        if is_end {
            break;
        }
        match state {
            0 => {
                if is_before {
                    let t = String::from_utf8_lossy(&line).into_owned();
                    let w: Vec<&str> = t.split_whitespace().collect();
                    call = w[1].parse().unwrap();
                    dt = f64::from_bits(u64::from_str_radix(w[2], 16).unwrap());
                    before.clear();
                    state = 1;
                }
            }
            1 => {
                if line.starts_with(b"AFTER ") && String::from_utf8_lossy(&line).trim_end() == format!("AFTER {}", call) {
                    after.clear();
                    state = 2;
                } else {
                    before.extend_from_slice(&line);
                }
            }
            _ => after.extend_from_slice(&line),
        }
    }
    if n_ok == 0 {
        finish("nosteps".to_string());
    }
    finish(format!("identical {}", n_ok));
}

fn find(hay: &[u8], needle: &[u8], from: usize) -> Option<usize> {
    if needle.len() > hay.len() {
        return None;
    }
    (from..=hay.len() - needle.len()).find(|&i| &hay[i..i + needle.len()] == needle)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        std::process::exit(2);
    }
    if args[2] == "stream" || args[2] == "streamfresh" {
        stream(&args[1], args[2] == "streamfresh");
        return;
    }
    let mode = if args.len() > 3 { args[3].as_str() } else { "step" };
    let text = std::fs::read(&args[2]).expect("routestep: cannot read the stream");
    let mut out: Vec<u8> = Vec::new();
    let mut p = 0usize;
    while p < text.len() {
        let eol = find(&text, b"\n", p).unwrap_or(text.len());
        let line = &text[p..eol];
        if !line.starts_with(b"BEFORE ") {
            // END / NOSTART / a stray line: copied, so both streams carry it
            out.extend_from_slice(&text[p..(eol + 1).min(text.len())]);
            p = eol + 1;
            continue;
        }
        let words: Vec<&[u8]> = line.split(|&c| c == b' ').collect();
        let call = std::str::from_utf8(words[1]).unwrap();
        let dt = f64::from_bits(u64::from_str_radix(std::str::from_utf8(words[2]).unwrap(), 16).unwrap());
        let tag = format!("\nAFTER {}\n", call);
        let body_end = find(&text, tag.as_bytes(), eol).expect("routestep: BEFORE without AFTER") + 1;
        let body = &text[eol + 1..body_end];
        // the AFTER block runs to the next BEFORE or the end
        let after_start = body_end + tag.len() - 1;
        let next = find(&text, b"\nBEFORE ", after_start - 1).map(|k| k + 1).unwrap_or(text.len());
        let next = if text[after_start..].starts_with(b"BEFORE ") { after_start } else { next };

        let mut e = SwmmEngine::new();
        let rc = e.open(Some(args[1].as_bytes()), Some(b"/dev/null"), Some(b"/dev/null"), None, false);
        if rc != 0 {
            panic!("routestep: the translated open failed ({})", rc);
        }
        if e.gate_init_hydraulics().is_err() || e.gate_init_quality().is_err() {
            panic!("routestep: a translated init leg terminated");
        }
        *e.context_mut() = SimulationContext::default();
        let mut r = ctxload::Rd::new(body);
        ctxload::load(&mut r, e.context_mut());
        ctxload::load_rt(&mut r, &mut e);
        ctxload::load_qual(&mut r, &mut e);
        load_stepdump(&mut r, &mut e);
        r.finish();
        if mode == "load" {
            out.extend_from_slice(&text[p..eol + 1]);
            dumpall(&e, &mut out);
            out.extend_from_slice(&text[body_end..next]);
        } else {
            out.extend_from_slice(&text[p..eol + 1]);
            out.extend_from_slice(body);
            out.extend_from_slice(format!("AFTER {}\n", call).as_bytes());
            match e.step_routing(dt) {
                Ok(()) => dumpall(&e, &mut out),
                Err(t) if t.exception == THROWN_RUNTIME_ERROR => {
                    out.extend_from_slice(format!("throw {}\n", t.what).as_bytes());
                }
                Err(t) if t.exception == STUB_EXCEPTION => {
                    out.extend_from_slice(format!("nottranslated {}\n", t.what).as_bytes());
                }
                Err(t) => {
                    out.extend_from_slice(format!("terminate {} {}\n", t.exception, t.what).as_bytes());
                }
            }
        }
        p = next;
    }
    std::io::stdout().lock().write_all(&out).unwrap();
}
