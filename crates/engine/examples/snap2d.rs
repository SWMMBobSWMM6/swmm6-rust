// Rust side of gate/diff/snap2d.cpp, the 2D output-side differential (phase 2a wave 5).
//
//   snap2d deck rpt out
// Runs the deck as the CLI does (open, initialize, start(save = 1), the step loop, end,
// report, close) and prints, in the C++ side's format: at every report period the 2D surface
// block fill_surface_snapshot filled ("SNAP k t <bits>", the counts, one line per array with
// its length and FNV-1a 64 hash; every value with SNAP2D_FULL=1), after end() "END <rc>" and
// ctx.mass_balance_2d as hex bits, then "DONE <exit code>". Where the translation reaches code
// not translated yet it prints "STUB <what>" and stops (exit 3); where upstream would terminate
// it aborts.
// Build: cargo build --release --features "twod gate" --example snap2d (into target2d).
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
use std::io::Write;
use swmm6_rust::core::engine::{SwmmEngine, SWMM_OK};
use swmm6_rust::core::engine_init_hydraulics::Terminated;
use swmm6_rust::core::engine_output_twod::gate_set_snapshot_hook;
use swmm6_rust::core::engine_step::STUB_EXCEPTION;
use swmm6_rust::core::output_plugin::SimulationSnapshot;

fn fnv(b: &[u8]) -> u64 {
    let mut h: u64 = 1469598103934665603;
    for &x in b {
        h ^= x as u64;
        h = h.wrapping_mul(1099511628211);
    }
    h
}

fn vd(o: &mut String, full: bool, name: &str, v: &[f64]) {
    let mut b = Vec::with_capacity(v.len() * 8);
    for x in v {
        b.extend_from_slice(&x.to_bits().to_le_bytes());
    }
    o.push_str(&format!("  {} n {} h {:016x}\n", name, v.len(), fnv(&b)));
    if full {
        for (i, x) in v.iter().enumerate() {
            o.push_str(&format!("    {}[{}] {:016x}\n", name, i, x.to_bits()));
        }
    }
}

fn out(s: &str) {
    let _ = std::io::stdout().lock().write_all(s.as_bytes());
}

fn call(r: Result<i32, Terminated>) -> i32 {
    match r {
        Ok(rc) => rc,
        Err(t) if t.exception == STUB_EXCEPTION => {
            out(&format!("STUB {}\n", t.what));
            let _ = std::io::stdout().lock().flush();
            std::process::exit(3);
        }
        Err(t) => {
            let _ = std::io::stdout().lock().flush();
            eprintln!("terminate called after throwing an instance of '{}'\n  what():  {}", t.exception, t.what);
            std::process::abort();
        }
    }
}

fn main() {
    use std::os::unix::ffi::OsStringExt;
    let args: Vec<Vec<u8>> = std::env::args_os().map(|a| a.into_vec()).collect();
    if args.len() < 4 {
        std::process::exit(2);
    }
    std::panic::set_hook(Box::new(|info| {
        let msg = if let Some(s) = info.payload().downcast_ref::<&str>() {
            s.to_string()
        } else if let Some(s) = info.payload().downcast_ref::<String>() {
            s.clone()
        } else {
            String::new()
        };
        let _ = std::io::stdout().lock().flush();
        if msg.starts_with("STUB") {
            println!("STUB {}", msg);
            let _ = std::io::stdout().lock().flush();
            std::process::exit(3);
        }
        eprintln!("snap2d: panic: {}", msg);
        std::process::exit(101);
    }));
    let full = std::env::var_os("SNAP2D_FULL").is_some();
    let mut k = 0;
    gate_set_snapshot_hook(Some(Box::new(move |snap: &SimulationSnapshot| {
        k += 1;
        let s = &snap.surface;
        let mut o = format!("SNAP {} t {:016x}\n", k, snap.sim_time.to_bits());
        o.push_str(&format!(
            "  counts {} {} {} {} species {} buildup {} gw {}\n",
            s.tri_count, s.vert_count, s.quad_count, s.edge_stride, s.species_count, s.buildup_count, s.gw2d_active as i32
        ));
        o.push_str(&format!("  cell_nv n {} h {:016x}\n", s.cell_nv.len(), fnv(&s.cell_nv)));
        vd(&mut o, full, "depth", &s.depth);
        vd(&mut o, full, "head", &s.head);
        vd(&mut o, full, "grad_hx", &s.grad_hx);
        vd(&mut o, full, "grad_hy", &s.grad_hy);
        vd(&mut o, full, "grad_hx_lim", &s.grad_hx_lim);
        vd(&mut o, full, "grad_hy_lim", &s.grad_hy_lim);
        for (i, f) in s.snow.iter().enumerate() {
            o.push_str(&format!("  snow[{}] n {}\n", i, f.len()));
        }
        vd(&mut o, full, "rainfall", &s.rainfall);
        vd(&mut o, full, "coupling_flux", &s.coupling_flux);
        vd(&mut o, full, "net_source", &s.net_source);
        vd(&mut o, full, "infil_rate", &s.infil_rate);
        vd(&mut o, full, "infil_cum", &s.infil_cum);
        vd(&mut o, full, "rain_cum", &s.rain_cum);
        vd(&mut o, full, "species_conc", &s.species_conc);
        vd(&mut o, full, "buildup", &s.buildup);
        vd(&mut o, full, "edge_flux", &s.edge_flux);
        vd(&mut o, full, "vert_head", &s.vert_head);
        vd(&mut o, full, "vert_depth", &s.vert_depth);
        vd(&mut o, full, "face_vx", &s.face_vx);
        vd(&mut o, full, "face_vy", &s.face_vy);
        vd(&mut o, full, "continuity_err", &s.continuity_err);
        vd(&mut o, full, "stat_max_depth", &s.stat_max_depth);
        vd(&mut o, full, "stat_max_velocity", &s.stat_max_velocity);
        vd(&mut o, full, "stat_max_cont_err", &s.stat_max_cont_err);
        out(&o);
    })));

    let (inp, rpt, outp) = (&args[1], &args[2], &args[3]);
    let mut engine = SwmmEngine::new();
    let err = engine.open(Some(inp), Some(rpt), Some(outp), None, false);
    if err != SWMM_OK {
        drop(engine);
        out(&format!("DONE {}\n", err));
        return;
    }
    let err = call(engine.initialize());
    if err != SWMM_OK {
        engine.close();
        drop(engine);
        out(&format!("DONE {}\n", err));
        return;
    }
    let err = call(engine.start(if outp.is_empty() { 0 } else { 1 }));
    if err != SWMM_OK {
        engine.close();
        drop(engine);
        out(&format!("DONE {}\n", err));
        return;
    }
    let mut elapsed = 0.0f64;
    let mut err;
    loop {
        err = call(engine.step(&mut elapsed));
        if err != SWMM_OK || elapsed <= 0.0 {
            break;
        }
    }
    let mut result = err;
    let rc_end = call(engine.end());
    if rc_end != SWMM_OK && result == SWMM_OK {
        result = rc_end;
    }
    out(&format!("END {}\n", rc_end));
    {
        let m = &engine.context().mass_balance_2d;
        let d: [(&str, f64); 21] = [
            ("init_storage", m.init_storage),
            ("final_storage", m.final_storage),
            ("external_in", m.external_in),
            ("external_out", m.external_out),
            ("rainfall_in", m.rainfall_in),
            ("coupling_1d_to_2d_in", m.coupling_1d_to_2d_in),
            ("coupling_2d_to_1d_out", m.coupling_2d_to_1d_out),
            ("outfall_in", m.outfall_in),
            ("outfall_out", m.outfall_out),
            ("boundary_in", m.boundary_in),
            ("boundary_out", m.boundary_out),
            ("evap_out", m.evap_out),
            ("infil_out", m.infil_out),
            ("infil_to_aquifer", m.infil_to_aquifer),
            ("infil_aquifer_pending", m.infil_aquifer_pending),
            ("aquifer_in", m.aquifer_in),
            ("solver_avg_h", m.solver_avg_h),
            ("solver_last_h", m.solver_last_h),
            ("solver_active_min", m.solver_active_min),
            ("solver_active_mean", m.solver_active_mean),
            ("solver_active_max", m.solver_active_max),
        ];
        let mut o = String::new();
        for (k, v) in d {
            o.push_str(&format!("  mb2.{} {:016x}\n", k, v.to_bits()));
        }
        o.push_str(&format!("  mb2.active {}\n", m.active as i32));
        o.push_str(&format!("  mb2.solver_nsteps {}\n", m.solver_nsteps));
        o.push_str(&format!("  mb2.solver_nrhs {}\n", m.solver_nrhs));
        o.push_str(&format!("  mb2.solver_n_tiers {}\n", m.solver_n_tiers));
        for k in 0..8 {
            o.push_str(&format!("  mb2.solver_tier_cells[{}] {}\n", k, m.solver_tier_cells[k]));
        }
        out(&o);
    }
    let s = call(engine.report());
    if s != SWMM_OK && result == SWMM_OK {
        result = s;
    }
    let s = engine.close();
    if s != SWMM_OK && result == SWMM_OK {
        result = s;
    }
    drop(engine);
    out(&format!("DONE {}\n", result));
}
