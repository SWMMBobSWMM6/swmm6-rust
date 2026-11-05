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
use swmm6_rust::core::engine_step_routing::StepError;

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
                Err(StepError::Exception(w)) => {
                    out.extend_from_slice(b"throw ");
                    out.extend_from_slice(&w);
                    out.push(b'\n');
                }
                Err(StepError::NotTranslated(w)) => {
                    out.extend_from_slice(format!("nottranslated {}\n", w).as_bytes());
                }
            }
        }
        p = next;
    }
    std::io::stdout().lock().write_all(&out).unwrap();
}
