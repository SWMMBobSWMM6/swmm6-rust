// Rust side of gate/diff/outgate.cpp, the output-side gate (gate (a): the .out and .rpt
// writers without the runoff step).
//
//   outgate deck rpt out start
//       the translation's open, initialize and start as the CLI runs them (failure path
//       included); on success the plugin files are flushed and the process exits without
//       closing anything, as the C++ side's _Exit does. Prints "rc <code>", exits with it.
//   outgate deck rpt out stream UP_RPT UP_OUT
//       reads gate/diff/outgate's "run" stream on stdin. The translation opens, initializes
//       and starts the deck itself (its preamble and header, gated by "start"); then for each
//       "PERIOD p" block it loads upstream's state as it was just before that period's
//       postOutputSnapshot and calls the translated post_output_snapshot, which appends the
//       period's record to its .out; at "ENDSTATE" it loads the state end() starts from and
//       runs end, report and close as the CLI does. Then it compares: the exit code; the
//       .rpt with upstream's (UP_RPT) without the three wall-clock lines; the .out header,
//       each period record written against upstream's record of the same period in UP_OUT,
//       and the closing records. One result line: "identical P" (P periods compared),
//       "differ <what> ..." with the first differing period, byte or line, "nottranslated
//       <what>", "nostart <rc>" (upstream itself stops before stepping; then the files are
//       compared as in "start").
// Build: cargo build --release --features gate --example outgate
#[path = "common/ctxload.rs"]
mod ctxload;
use std::io::{BufRead, Write};
use swmm6_rust::core::context::SimulationContext;
use swmm6_rust::core::engine::SwmmEngine;
use swmm6_rust::core::engine_init_hydraulics::Terminated;
use swmm6_rust::core::engine_step::STUB_EXCEPTION;

fn finish(line: String) -> ! {
    let _ = writeln!(std::io::stdout(), "{}", line);
    let _ = std::io::stdout().flush();
    std::process::exit(0);
}

/// An engine call: its code; a stub ends the gate's run; any other Err is upstream's
/// std::terminate.
fn call(r: Result<i32, Terminated>, phase: &str) -> i32 {
    match r {
        Ok(rc) => rc,
        Err(t) if t.exception == STUB_EXCEPTION => finish(format!("nottranslated {} {}", phase, t.what)),
        Err(t) => {
            eprintln!("terminate called after throwing an instance of '{}'\n  what():  {}", t.exception, t.what);
            std::process::abort();
        }
    }
}

/// open, initialize, start as src/cli/main.cpp runs them; Ok(engine) after a successful
/// start, Err(rc) after the CLI's failure path.
fn cli_start(deck: &[u8], rpt: &[u8], out: &[u8]) -> Result<SwmmEngine, i32> {
    let mut e = SwmmEngine::new();
    let err = e.open(Some(deck), Some(rpt), Some(out), None, false);
    if err != 0 {
        drop(e);
        return Err(err);
    }
    let err = call(e.initialize(), "initialize");
    if err != 0 {
        e.close();
        drop(e);
        return Err(err);
    }
    let err = call(e.start(if out.is_empty() { 0 } else { 1 }), "start");
    if err != 0 {
        e.close();
        drop(e);
        return Err(err);
    }
    Ok(e)
}

fn load_outside(r: &mut ctxload::Rd, e: &mut SwmmEngine) {
    e.old_runoff_ms_ = r.pf("engine.old_runoff_ms_");
    e.new_runoff_ms_ = r.pf("engine.new_runoff_ms_");
    *e.gate_old_gw_rate_mut() = r.vf("engine.old_gw_rate_");
    e.save_results_ = r.pi("engine.save_results_") as i32;
    e.do_routing_ = r.pi("engine.do_routing_") != 0;
    e.hotstart_save_done_ = r.vi("engine.hotstart_save_done_").into_iter().map(|x| x as u8).collect();
    {
        let f = |v: Vec<f64>| v.into_iter().map(|x| x as f32).collect::<Vec<f32>>();
        let a = &mut e.avg_;
        a.node_depth = f(r.vf("engine.avg_.node_depth"));
        a.node_head = f(r.vf("engine.avg_.node_head"));
        a.node_volume = f(r.vf("engine.avg_.node_volume"));
        a.node_lat_inflow = f(r.vf("engine.avg_.node_lat_inflow"));
        a.node_total_inflow = f(r.vf("engine.avg_.node_total_inflow"));
        a.node_overflow = f(r.vf("engine.avg_.node_overflow"));
        a.link_flow = f(r.vf("engine.avg_.link_flow"));
        a.link_depth = f(r.vf("engine.avg_.link_depth"));
        a.link_velocity = f(r.vf("engine.avg_.link_velocity"));
        a.link_volume = f(r.vf("engine.avg_.link_volume"));
        a.link_capacity = f(r.vf("engine.avg_.link_capacity"));
        a.n_steps = r.pi("engine.avg_.n_steps") as i32;
    }
    let h = e.gate_hydrology();
    h.runoff_.soa_.area = r.vf("runoff_.soa_.area");
    {
        let g = &mut h.groundwater_.soa_;
        g.gw_flow = r.vf("groundwater_.gw_flow");
        g.bottom_elev = r.vf("groundwater_.bottom_elev");
        g.lower_depth = r.vf("groundwater_.lower_depth");
        g.theta = r.vf("groundwater_.theta");
        g.upper_evap = r.vf("groundwater_.upper_evap");
        g.lower_evap = r.vf("groundwater_.lower_evap");
        g.total_depth = r.vf("groundwater_.total_depth");
        g.porosity = r.vf("groundwater_.porosity");
    }
    {
        let sn = &mut h.snow_.soa_;
        sn.n_subcatch = r.pi("snow_.n_subcatch") as i32;
        sn.external = r.vi("snow_.external").into_iter().map(|x| x != 0).collect();
        sn.external_net = r.vf("snow_.external_net");
        sn.imelt = r.vf("snow_.imelt");
        sn.fArea = r.vf("snow_.fArea");
        sn.enabled = r.vi("snow_.enabled").into_iter().map(|x| x as u8).collect();
    }
    let ng = r.pi("lid_.groups.n") as usize;
    if ng != h.lid_.groups_.len() {
        finish(format!("loaddiff lid_.groups.n {} against {}", ng, h.lid_.groups_.len()));
    }
    for t in 0..ng {
        let g = &mut h.lid_.groups_[t];
        g.count = r.pi("lid_.group.count") as i32;
        g.subcatch_idx = r.vi("lid_.group.subcatch_idx").into_iter().map(|x| x as i32).collect();
        g.control_idx = r.vi("lid_.group.control_idx").into_iter().map(|x| x as i32).collect();
        g.area = r.vf("lid_.group.area");
        g.evap_loss = r.vf("lid_.group.evap_loss");
        g.infil_loss = r.vf("lid_.group.infil_loss");
        g.wb_inflow = r.vf("lid_.group.wb_inflow");
        g.wb_evap = r.vf("lid_.group.wb_evap");
        g.wb_infil = r.vf("lid_.group.wb_infil");
        g.wb_surf_flow = r.vf("lid_.group.wb_surf_flow");
        g.wb_drain_flow = r.vf("lid_.group.wb_drain_flow");
        g.wb_init_vol = r.vf("lid_.group.wb_init_vol");
        g.wb_final_vol = r.vf("lid_.group.wb_final_vol");
    }
    let np = r.pi("engine.out_plugin.n_periods_");
    if np >= 0 {
        for p in e.gate_plugins_mut().output_plugins_mut().iter_mut() {
            p.gate_set_n_periods(np as i32);
        }
    }
}

fn load_state(e: &mut SwmmEngine, body: &[u8]) {
    *e.context_mut() = SimulationContext::default();
    let mut r = ctxload::Rd::new(body);
    ctxload::load(&mut r, e.context_mut());
    ctxload::load_rt(&mut r, e);
    ctxload::load_qual(&mut r, e);
    load_outside(&mut r, e);
    r.finish();
}

fn out_pos(e: &mut SwmmEngine) -> i64 {
    e.gate_plugins_mut().output_plugins_mut().first().and_then(|p| p.out_file_.as_ref().map(|f| f.tell())).unwrap_or(-1)
}

fn rpt_lines(b: &[u8]) -> Vec<&[u8]> {
    let ts = |l: &[u8]| -> bool {
        let s = String::from_utf8_lossy(l);
        s.contains("Analysis begun") || s.contains("Analysis ended") || s.contains("Total elapsed time")
    };
    b.split(|&c| c == b'\n').filter(|l| !ts(l)).collect()
}

/// The first differing line of two reports without their clock lines, or None.
fn rpt_diff(a: &[u8], b: &[u8]) -> Option<String> {
    let (la, lb) = (rpt_lines(a), rpt_lines(b));
    for (i, (x, y)) in la.iter().zip(lb.iter()).enumerate() {
        if x != y {
            return Some(format!("line {}: upstream {:?} translation {:?}", i + 1, String::from_utf8_lossy(&x[..x.len().min(100)]), String::from_utf8_lossy(&y[..y.len().min(100)])));
        }
    }
    if la.len() != lb.len() {
        return Some(format!("lines: upstream {} translation {}", la.len(), lb.len()));
    }
    None
}

fn first_byte_diff(a: &[u8], b: &[u8]) -> Option<usize> {
    let n = a.len().min(b.len());
    (0..n).find(|&i| a[i] != b[i]).or(if a.len() != b.len() { Some(n) } else { None })
}

fn read(p: &str) -> Vec<u8> {
    std::fs::read(p).unwrap_or_default()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 5 {
        std::process::exit(2);
    }
    let (deck, rpt, out) = (args[1].as_bytes(), args[2].as_bytes(), args[3].as_bytes());
    if args[4] == "start" {
        match cli_start(deck, rpt, out) {
            Ok(mut e) => {
                e.gate_flush_plugin_files();
                println!("rc 0");
                let _ = std::io::stdout().flush();
                std::process::exit(0);
            }
            Err(rc) => {
                println!("rc {}", rc);
                std::process::exit(rc);
            }
        }
    }
    // stream
    let (up_rpt, up_out) = (&args[5], &args[6]);
    let stdin = std::io::stdin();
    let mut rd = std::io::BufReader::with_capacity(1 << 20, stdin.lock());
    let mut line: Vec<u8> = Vec::new();
    let mut body: Vec<u8> = Vec::new();
    // what the current block is: None, Some(period) for PERIOD, Some(0) for ENDSTATE
    let mut block: Option<i64> = None;
    let mut eng: Option<SwmmEngine> = None;
    let mut records: Vec<(i64, i64, i64)> = Vec::new();
    let mut ended = false;
    let mut rust_rc = 0;
    let mut up_rc: Option<i32> = None;
    loop {
        line.clear();
        let n = rd.read_until(b'\n', &mut line).expect("outgate: read error");
        let eof = n == 0;
        let header = line.starts_with(b"PERIOD ") || line.starts_with(b"ENDSTATE ") || line.starts_with(b"DONE ") || line.starts_with(b"rc ");
        if (eof || header) && block.is_some() {
            let e = eng.as_mut().unwrap();
            load_state(e, &body);
            let k = block.unwrap();
            if k > 0 {
                let a = out_pos(e);
                if let Err(t) = e.gate_post_output_snapshot() {
                    if t.exception == STUB_EXCEPTION {
                        finish(format!("nottranslated period {} {}", k, t.what));
                    }
                    finish(format!("terminate period {} {} {}", k, t.exception, t.what));
                }
                let b = out_pos(e);
                records.push((k, a, b));
            } else {
                // the CLI's end, report and close
                let mut result = rust_rc;
                let mut fin = |s: i32| {
                    if s != 0 && result == 0 {
                        result = s;
                    }
                };
                fin(call(e.end(), "end"));
                fin(call(e.report(), "report"));
                fin(e.close());
                rust_rc = result;
                ended = true;
            }
            body.clear();
            block = None;
        }
        if eof {
            break;
        }
        if line.starts_with(b"rc ") {
            // upstream stopped before stepping
            let rc: i32 = String::from_utf8_lossy(&line[3..]).trim().parse().unwrap();
            let r = match cli_start(deck, rpt, out) {
                Ok(e) => {
                    drop(e);
                    0
                }
                Err(rc) => rc,
            };
            if r != rc {
                finish(format!("differ rc upstream {} translation {}", rc, r));
            }
            if let Some(d) = rpt_diff(&read(up_rpt), &read(&args[2])) {
                finish(format!("differ rpt nostart {}", d));
            }
            let (uo, ro) = (read(up_out), read(&args[3]));
            if let Some(k) = first_byte_diff(&uo, &ro) {
                finish(format!("differ out nostart byte {} (sizes {} {})", k, uo.len(), ro.len()));
            }
            finish(format!("nostart {}", rc));
        }
        if line.starts_with(b"PERIOD ") || line.starts_with(b"ENDSTATE ") {
            if eng.is_none() {
                match cli_start(deck, rpt, out) {
                    Ok(e) => eng = Some(e),
                    Err(rc) => finish(format!("differ the translation's start failed ({}) where upstream's ran", rc)),
                }
            }
            let t = String::from_utf8_lossy(&line).into_owned();
            let w: Vec<&str> = t.split_whitespace().collect();
            if w[0] == "PERIOD" {
                block = Some(w[1].parse().unwrap());
            } else {
                rust_rc = w[1].parse().unwrap();
                block = Some(0);
            }
            continue;
        }
        if line.starts_with(b"DONE ") {
            up_rc = Some(String::from_utf8_lossy(&line[5..]).trim().parse().unwrap());
            continue;
        }
        if block.is_some() {
            body.extend_from_slice(&line);
        }
    }
    drop(eng);
    if !ended || up_rc.is_none() {
        finish("noresult".to_string());
    }
    let up_rc = up_rc.unwrap();
    if up_rc != rust_rc {
        finish(format!("differ rc upstream {} translation {}", up_rc, rust_rc));
    }
    // the .out: header, each period written, closing records
    let (uo, ro) = (read(up_out), read(&args[3]));
    if uo.len() >= 24 {
        let ints = |b: &[u8], k: usize| i32::from_ne_bytes(b[b.len() - 24 + 4 * k..b.len() - 20 + 4 * k].try_into().unwrap()) as i64;
        let out_start = ints(&uo, 2) as usize;
        let nper = ints(&uo, 3);
        let rec = if nper > 0 { (uo.len() - out_start - 24) / nper as usize } else { 0 };
        if ro.len() < out_start + 24 || uo[..out_start] != ro[..out_start] {
            let k = first_byte_diff(&uo[..out_start.min(uo.len())], &ro[..out_start.min(ro.len())]).unwrap_or(out_start);
            finish(format!("differ out header byte {}", k));
        }
        for &(p, a, b) in &records {
            let (a, b) = (a as usize, b as usize);
            if b - a != rec {
                finish(format!("differ out period {} record size {} against {}", p, b - a, rec));
            }
            let u0 = out_start + (p as usize - 1) * rec;
            if let Some(k) = first_byte_diff(&uo[u0..u0 + rec], &ro[a..b]) {
                finish(format!("differ out period {} byte {} of {}", p, k, rec));
            }
        }
        if uo[uo.len() - 24..] != ro[ro.len() - 24..] {
            finish(format!("differ out closing records {:?} {:?}", &uo[uo.len() - 24..], &ro[ro.len() - 24..]));
        }
        let written: usize = records.iter().map(|r| (r.2 - r.1) as usize).sum();
        if ro.len() != out_start + written + 24 {
            finish(format!("differ out size {} against {}", ro.len(), out_start + written + 24));
        }
    } else if !ro.is_empty() || !uo.is_empty() {
        if let Some(k) = first_byte_diff(&uo, &ro) {
            finish(format!("differ out byte {} (sizes {} {})", k, uo.len(), ro.len()));
        }
    }
    if let Some(d) = rpt_diff(&read(up_rpt), &read(&args[2])) {
        finish(format!("differ rpt {}", d));
    }
    finish(format!("identical {}", records.len()));
}
