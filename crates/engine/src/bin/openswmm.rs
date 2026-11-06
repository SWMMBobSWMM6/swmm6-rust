// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/cli/main.cpp (HydroCouple/openswmm.engine, branch
// swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the Apache
// License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The command-line program:
//   openswmm <input.inp> <report.rpt> [output.out]
//   openswmm --help | -h
//   openswmm --version | -v
// It runs the engine through the same sequence as upstream's CLI (open, initialize, start, the
// step loop, end, report, close), prints the same lines and returns the same exit codes.
//
// Differences of form:
//   * The C API calls (swmm_engine_create, _open, ...) are the engine's methods; the handle's
//     destruction is the engine going out of scope.
//   * Where the engine returns Err(Terminated) (upstream's process would end in
//     std::terminate), the program prints libstdc++'s terminate message and aborts, as
//     upstream's does (exit status 134). A call into code the translation does not have yet
//     (Err with STUB_EXCEPTION, or a "STUB:" panic) ends the program with a message on
//     stderr and exit status 3 (upstream has no such case).

use std::io::Write;
use swmm6_rust::core::engine::{SwmmEngine, SWMM_OK};
use swmm6_rust::core::engine_init_hydraulics::Terminated;
use swmm6_rust::core::engine_step::STUB_EXCEPTION;
use swmm6_rust::core::report_plugin::OPENSWMM_VERSION_FULL;

fn out(b: &[u8]) {
    let _ = std::io::stdout().lock().write_all(b);
}
fn flush() {
    let _ = std::io::stdout().lock().flush();
}
fn v() -> String {
    String::from_utf8_lossy(OPENSWMM_VERSION_FULL).into_owned()
}

fn print_help() {
    out(b"\n");
    out(format!("Open Source Storm Water Management Model {}\n", v()).as_bytes());
    out(b"===================================================\n\n");
    out(b"USAGE:\n");
    out(b"  openswmm <input.inp> <report.rpt> [output.out]\n\n");
    out(b"OPTIONS:\n");
    out(b"  --help, -h       Show this help message\n");
    out(b"  --version, -v    Show version number\n\n");
    out(b"DESCRIPTION:\n");
    out(b"  Runs a SWMM simulation using the new data-oriented engine.\n");
    out(b"  The engine reads a standard SWMM .inp file and produces\n");
    out(b"  a report (.rpt) and optional binary output (.out) file.\n\n");
    out(b"  For the legacy EPA SWMM 5.x engine, use: openswmm-legacy\n\n");
}

fn print_version() {
    out(format!("openswmm.engine {}\n", v()).as_bytes());
}

/// An engine call: its return code, or the end of the process where upstream's ends.
fn call(r: Result<i32, Terminated>) -> i32 {
    match r {
        Ok(rc) => rc,
        Err(t) if t.exception == STUB_EXCEPTION => {
            flush();
            eprintln!("openswmm (swmm6-rust): not translated yet: {}", t.what);
            std::process::exit(3);
        }
        Err(t) => {
            flush();
            eprintln!("terminate called after throwing an instance of '{}'\n  what():  {}", t.exception, t.what);
            std::process::abort();
        }
    }
}

fn lossy(b: &[u8]) -> String {
    String::from_utf8_lossy(b).into_owned()
}

fn run(args: &[Vec<u8>]) -> i32 {
    let argc = args.len();
    if argc < 2 {
        out(b"\nError: not enough arguments. Use --help for usage.\n\n");
        return 1;
    }
    if args[1] == b"--help" || args[1] == b"-h" {
        print_help();
        return 0;
    }
    if args[1] == b"--version" || args[1] == b"-v" {
        print_version();
        return 0;
    }
    if argc < 3 {
        out(b"\nError: need at least input and report file paths.\n");
        out(b"Usage: openswmm <input.inp> <report.rpt> [output.out]\n\n");
        return 1;
    }
    let inp_file = &args[1];
    let rpt_file = &args[2];
    let empty: Vec<u8> = Vec::new();
    let out_file = if argc > 3 { &args[3] } else { &empty };

    let mut m = format!("\n... OpenSWMM Engine {}\n... Input: ", v()).into_bytes();
    m.extend_from_slice(inp_file);
    m.extend_from_slice(b"\n... Report: ");
    m.extend_from_slice(rpt_file);
    m.extend_from_slice(b"\n... Output: ");
    m.extend_from_slice(out_file);
    m.push(b'\n');
    out(&m);
    flush();

    let start = std::time::SystemTime::now();

    // ---- Create engine ----
    let mut engine = SwmmEngine::new();

    // ---- Open input file ----
    let err = engine.open(Some(inp_file), Some(rpt_file), Some(out_file), None, false);
    if err != SWMM_OK {
        out(format!("Error opening input file: {}\n", lossy(engine.last_error_message())).as_bytes());
        drop(engine);
        return err;
    }

    // ---- Initialize ----
    let err = call(engine.initialize());
    if err != SWMM_OK {
        out(format!("Error initializing: {}\n", lossy(engine.last_error_message())).as_bytes());
        engine.close();
        drop(engine);
        return err;
    }

    // ---- Start simulation ----
    let save = if !out_file.is_empty() { 1 } else { 0 };
    let err = call(engine.start(save));
    if err != SWMM_OK {
        out(format!("Error starting simulation: {}\n", lossy(engine.last_error_message())).as_bytes());
        engine.close();
        drop(engine);
        return err;
    }

    // ---- Step loop ----
    let mut elapsed = 0.0f64;
    let mut step_count: i64 = 0;
    let mut err;
    loop {
        err = call(engine.step(&mut elapsed));
        if err != SWMM_OK {
            out(format!("\nError at step {}: {}\n", step_count, lossy(engine.last_error_message())).as_bytes());
            break;
        }
        if elapsed <= 0.0 {
            break;
        }
        step_count += 1;
        if step_count % 1000 == 0 {
            out(&swmm6_rust::cformat!(b"\r... %ld steps, elapsed = %.4f days", step_count, elapsed));
            flush();
        }
    }
    out(format!("\n... {} steps completed.\n", step_count).as_bytes());

    // ---- End and report ----
    let mut result = err;
    let mut finalize_status = |engine: &SwmmEngine, status: i32, phase: &str| {
        if status != SWMM_OK {
            out(format!("Error {}: {}\n", phase, lossy(engine.last_error_message())).as_bytes());
            if result == SWMM_OK {
                result = status;
            }
        }
    };
    let s = call(engine.end());
    finalize_status(&engine, s, "ending simulation");
    let s = call(engine.report());
    finalize_status(&engine, s, "writing report");

    // ---- Close and destroy ----
    let s = engine.close();
    finalize_status(&engine, s, "closing simulation");
    drop(engine);

    // ---- Summary ----
    let run_time = start.elapsed().map(|d| d.as_secs()).unwrap_or(0) as f64;
    out(&swmm6_rust::cformat!(b"... OpenSWMM Engine %s in %.2f seconds.\n\n", if result == SWMM_OK { b"completed" as &[u8] } else { b"failed" }, run_time));
    result
}

fn main() {
    use std::os::unix::ffi::OsStringExt;
    let args: Vec<Vec<u8>> = std::env::args_os().map(|a| a.into_vec()).collect();
    // A "STUB:" panic (code not translated yet) ends with exit status 3, as a stub Err does.
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
            eprintln!("openswmm (swmm6-rust): not translated yet: {}", msg);
            std::process::exit(3);
        }
        eprintln!("openswmm (swmm6-rust): panic: {} ({})", msg, info.location().map(|l| l.to_string()).unwrap_or_default());
        std::process::exit(101);
    }));
    let rc = run(&args);
    flush();
    std::process::exit(rc);
}
