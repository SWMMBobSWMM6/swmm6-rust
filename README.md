# swmm6-rust

The OpenSWMM 6 engine's 1D path, translated from C++ into safe Rust and checked byte for byte against OpenSWMM 6 itself.

The goal is a Rust engine that gives the same `.out` file, bit for bit, as the C++ engine on the same deck. It is not a rewrite or a new solver. Every floating-point operation keeps upstream's order, so the same inputs give the same bits. Where upstream does something odd, the translation does the same odd thing, and the oddity is written down in [UPSTREAM_NOTES.md](UPSTREAM_NOTES.md).

## Status

Phase 1a (1D) is under way. Done, and identical to upstream on every gate deck:

* reading the `.inp` (every section handler), cross-reference resolution, and `open()` in strict, lenient and metadata modes;
* the initialize legs: `initHydraulics`, `initHydrology`, `initQuality`, `initGeometry`, `initMassBalance`;
* the leaf numerics: cross-section geometry and its tables, transects, streets, force mains, the C library functions upstream depends on (`strtod`, `from_chars`, `sscanf`, `%.10g`, `std::sort` and the `unordered_map`/`map` orders);
* the run lifecycle: `initialize()` as a whole and `start()` (identical on all 932 gate decks), and `step()` with its own bookkeeping (statistics, routing mass balance, final storage), gated helper by helper;
* the runoff (hydrology) step, `stepRunoff`: rain gages, climate, the nonlinear reservoir, all infiltration methods, degree-day snow, groundwater, RDII, the LID units, buildup and washoff, the runoff mass balance (identical on every gate deck that reaches it, loaded from upstream's own state mid-run; the SWMM6-only UEB snow, heat and reaction-species hooks are still stubs);
* the hydraulic routing step (`stepRouting`) for STEADY, KINWAVE and DYNWAVE routing: the dynamic-wave solver (EXTRAN and SLOT surcharge, virtual junctions), kinematic wave and steady flow, nodes, outfalls, dividers, storage and its exfiltration, pumps, orifices, weirs, outlets, culverts, force mains, roadways, street and channel inlets, node LID columns, control rules and the inflows. It is checked call by call against upstream's own `stepRouting` (`ctxgate.sh route`, PLAN.md);
* the FV router (FLOW_ROUTING FV) on its CPU explicit backend: the HLL/HLLC finite-volume step with MUSCL-Hancock reconstruction, local time stepping, the algebraic junction solves, structures, virtual junctions, Preissmann-slot surcharge and storage nodes, the Router's forcing and publish, the FV divergence check and run statistics, and FV water quality transport (the conservative complete-mix solve with its GMRES, whose `long double` dot product is emulated exactly in `core/ext80.rs`), checked by `ctxgate.sh route` and end to end on every FV gate deck;
* the routing step's water quality under the LEGACY quality solver: pollutant loads from every source, mixing in nodes and links, decay, evaporation and seepage, `[TREATMENT]` expressions, the storage-node LID layer treatment, the SWMM6 water-age and temperature mirrors, and reaction species (transport and kinetics with the EUL and RK5 solvers), also checked by `ctxgate.sh route`.

* the output side and the end of the run: the `.out` writer (header, every report period, the closing records), the `.rpt` writer (preamble, warnings and errors, every continuity block and summary table, the timing lines, the report of a run that fails in open, initialize or start), the report-period snapshot with its interpolation and time-step averaging, SAVE HOTSTART in the legacy format, `end()`, `report()`, `close()`, and the command-line program `openswmm` (`src/bin/openswmm.rs`, upstream's `src/cli/main.cpp`). Checked by `ctxgate.sh out`, which loads upstream's own mid-run and end-of-run state into the translated writers (PLAN.md).

Not yet translated: the EULERIAN_ARD and LAGRANGIAN quality engines, the heat flux modules and the other SWMM6-only hooks PLAN.md lists, and the FV solver options no gate deck sets (FV_PRESSURIZED_IMPLICIT with its implicit and pressurized solvers, TPA, RK2, unsteady friction); `step()` calls the missing pieces through stub files, and `openswmm` stops at one with exit status 3. A deck that needs none of them runs end to end through `openswmm`, and on every such gate deck (904 of the 932: all 480 corpus decks and 424 of the 452 extra decks) the `.out` is byte-identical to the C++ engine's, the `.rpt` identical but its three clock lines, and the exit code the same. [PLAN.md](PLAN.md) has the phases, the end-to-end counts and what comes next. 2D comes after 1D.

## How it is checked

The reference is OpenSWMM 6 itself, built here from the source below. Four kinds of gate:

| Gate | What must hold |
|---|---|
| G1 | `.out` byte-identical to the C++ engine on every gate deck |
| G2 | `.rpt` identical except the wall-clock lines |
| G3 | State after every routing step bit-identical |
| G4 | Leaf functions identical to upstream's own code, compiled from upstream source |

Until the time step exists, the state gates stand in for G1–G3. Each one opens a deck in both engines, dumps every member of the engine's context and solvers (doubles as hex bits), and compares the dumps and the exit codes. The deck set is 688 decks: 236 from the [1729-SWMM5-Models-2030](https://github.com/SWMMBobSWMM6/1729-SWMM5-Models-2030) corpus, 336 of upstream's own test decks, and 116 synthetic decks in `gate/decks`. As of 2 November 2026 every mode is identical on every deck. The table is in PLAN.md.

Where upstream crashes (an exception escaping a `noexcept` function), the Rust function returns `Err(Terminated)` at the same point with the same exception type and message. The library never aborts; only the gate harness does, so the gate can compare exit codes.

## Layout

| Path | What |
|---|---|
| `crates/engine` | the library (`swmm6_rust`), `#![forbid(unsafe_code)]`, and the command-line program `openswmm` (`src/bin`) |
| `crates/engine/examples` | the Rust side of each gate harness |
| `gate/` | the C++ harnesses that link the reference library, the deck set, `ctxgate.sh` |
| `tools/` | generators for the state dumpers and loader, and a first-pass translator for structure-of-arrays headers |
| `census/` | OpenSWMM 6 against EPA SWMM 5.2.4 on a corpus of SWMM5 models |
| `TRANSLATING.md` | the translation rules |
| `UPSTREAM_NOTES.md` | 96 upstream oddities found along the way, each quoting the code |

`cargo build --release` builds the library. The gates also need the reference build of OpenSWMM 6 and the deck corpora; `gate/build.sh` and `gate/ctxgate.sh` say where they expect them.

## License and attribution

This project is a translation of the OpenSWMM 6 engine.

OpenSWMM Engine: Copyright 2026 HydroCouple Developers, licensed under the Apache License, Version 2.0. See [LICENSE](LICENSE) and [NOTICE](NOTICE).
Source: https://github.com/HydroCouple/openswmm.engine (branch `swmm6_rel`, commit `88329d6df2f8dd63fe7e7c601954ec2bde3573e6`, 6.0.0-alpha.4).
Lead developer Caleb Buahin; developer and documentation Corinne Wiesner-Friedman.

Every translated file is derived from an upstream file, and its header names that file, keeps its copyright line, and says it was modified (translated into Rust).

Portions derive from the USEPA Storm Water Management Model, which is in the public domain under 17 USC § 105. The USEPA does not endorse this product.

swmm6-rust is not an OpenSWMM or HydroCouple release and is not endorsed by HydroCouple. The translation, the gate harnesses and the tools are distributed under the same Apache License, Version 2.0.

## AI use

Claude (Anthropic's AI model, in Claude Code) wrote the translation, the gate harnesses, the generators and these documents. Robert Dickinson set the goal, chose the approach (1D first, byte for byte against OpenSWMM 6), and directed the work. Nothing here is accepted on the strength of the translation looking right: each piece counts as done only when its gate shows it identical to the reference build of upstream.
