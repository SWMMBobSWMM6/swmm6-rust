# swmm6-rust

The OpenSWMM 6 engine's 1D path, translated from C++ into safe Rust and checked byte for byte against OpenSWMM 6 itself.

The goal is a Rust engine that gives the same `.out` file, bit for bit, as the C++ engine on the same deck. It is not a rewrite or a new solver. Every floating-point operation keeps upstream's order, so the same inputs give the same bits. Where upstream does something odd, the translation does the same odd thing, and the oddity is written down in [UPSTREAM_NOTES.md](UPSTREAM_NOTES.md).

## Status

Phase 1 (1D) is done for every path a gate deck reaches. Done, and identical to upstream on every gate deck:

* reading the `.inp` (every section handler), cross-reference resolution, and `open()` in strict, lenient and metadata modes;
* the initialize legs: `initHydraulics`, `initHydrology`, `initQuality`, `initGeometry`, `initMassBalance`;
* the leaf numerics: cross-section geometry and its tables, transects, streets, force mains, the C library functions upstream depends on (`strtod`, `from_chars`, `sscanf`, `%.10g`, `std::sort` and the `unordered_map`/`map` orders);
* the run lifecycle: `initialize()` as a whole and `start()` (identical on all 932 gate decks), and `step()` with its own bookkeeping (statistics, routing mass balance, final storage), gated helper by helper;
* the runoff (hydrology) step, `stepRunoff`: rain gages, climate, the nonlinear reservoir, all infiltration methods, degree-day snow, groundwater, RDII, the LID units, buildup and washoff, the runoff mass balance, and the SWMM6-only runoff hooks: UEB snow on subcatchments (with its modified force-restore point model), the heat and water-age hooks and LID layer columns, and the reaction species' buildup, washoff, ponded-pool and LID-unit kinetics (identical on every gate deck that reaches it, loaded from upstream's own state mid-run);
* the hydraulic routing step (`stepRouting`) for STEADY, KINWAVE and DYNWAVE routing: the dynamic-wave solver (EXTRAN and SLOT surcharge, virtual junctions), kinematic wave and steady flow, nodes, outfalls, dividers, storage and its exfiltration, pumps, orifices, weirs, outlets, culverts, force mains, roadways, street and channel inlets, node LID columns (with the Richards-equation columns: the adaptive BDF1 soil-column integrator and its ports), control rules and the inflows. It is checked call by call against upstream's own `stepRouting` (`ctxgate.sh route`, PLAN.md);
* the FV router (FLOW_ROUTING FV) on its CPU explicit backend: the HLL/HLLC finite-volume step with MUSCL-Hancock reconstruction, local time stepping, the algebraic junction solves, structures, virtual junctions, Preissmann-slot surcharge and storage nodes, the Router's forcing and publish, the FV divergence check and run statistics, and FV water quality transport (the conservative complete-mix solve with its GMRES, whose `long double` dot product is emulated exactly in `core/ext80.rs`), checked by `ctxgate.sh route` and end to end on every FV gate deck;
* the routing step's water quality under the LEGACY quality solver: pollutant loads from every source, mixing in nodes and links, decay, evaporation and seepage, `[TREATMENT]` expressions, the storage-node LID layer treatment, the SWMM6 water-age and temperature mirrors, and reaction species (transport and kinetics with the EUL and RK5 solvers, also in the layers and mobile water of LID storage nodes), also checked by `ctxgate.sh route`;
* the SWMM6 quality engines EULERIAN_ARD (the finite-volume species transport on the network mesh, with flux limiting, dispersion, node stores and reactions) and LAGRANGIAN (segment transport with random-walk dispersion), and the heat flux modules (surface and radiative exchange, the solar forcing, the bed), under every quality engine, checked by `ctxgate.sh route` and `step`.

* the output side and the end of the run: the `.out` writer (header, every report period, the closing records), the `.rpt` writer (preamble, warnings and errors, every continuity block and summary table, the timing lines, the report of a run that fails in open, initialize or start), the report-period snapshot with its interpolation and time-step averaging, SAVE HOTSTART in the legacy format, `end()`, `report()`, `close()`, and the command-line program `openswmm` (`src/bin/openswmm.rs`, upstream's `src/cli/main.cpp`). Checked by `ctxgate.sh out`, which loads upstream's own mid-run and end-of-run state into the translated writers (PLAN.md).

Every gate deck now runs end to end through `openswmm`: on all 932 (the 480 corpus decks and the 452 extra decks) the `.out` is byte-identical to the C++ engine's, the `.rpt` identical but its three clock lines, and the exit code the same. Still not translated, because no gate deck reaches them: the UEB EQUILIBRIUM surface-temperature model, the reaction solvers ROS2 and BDF2, EQUIL species and booster sources, the FV solver options FV_PRESSURIZED_IMPLICIT (with its implicit and pressurized solvers), TPA, RK2, unsteady friction and species on the FV mesh, reaction species, heat and bed exchange under FV quality transport, DYNAMIC_SLOT and TPA surcharge and physical junction areas in the dynamic-wave solver, control-rule expression premises, the runoff interface file (USE and SAVE RUNOFF), a native (OPENSWMM_HS) hotstart file and engine checkpointing. A deck that needs one of them stops where upstream would first act, before anything changes, and `openswmm` exits with status 3; [PLAN.md](PLAN.md) lists them with the reason for each, the phases, the end-to-end counts and what comes next.

Phase 2a (2D surface flow and 1D-2D coupling, under the `twod` cargo feature) is under way. With `twod`, `open()` reads every 2D section (the mesh, its options and rows, external mesh files, 2D infiltration, surface quality, groundwater transport and aquifer rows) and leaves the context and the 2D router identical to OpenSWMM 6's 2D build on 1,223 gate decks in strict, lenient and metadata modes; the 2D solver's leaf numerics match upstream's kernels in differential tests (all but the sign of a few NaNs from non-finite inputs); and `initialize()` and `start()` build the mesh topology, the rainfall weights, the 1D-2D coupling points, the solver and the 2D infiltration identically on every deck that does not need a later feature (1,125 of 1,223 gate decks under the 2D reference, 228 of them with an active mesh). A deck that needs a phase-2b feature (a 2D aquifer, surface quality or species transport, snow or heat on the mesh) stops where upstream would first act. The 2D time step on a mesh that exchanges no water with the 1D network is translated too: the explicit local-inertial / diffusive-wave / full shallow-water marcher with its local time stepping, wet/dry fronts, every boundary type, rain, 2D infiltration and evaporation, the optional momentum physics, and the router's batches and ledgers; every such deck stays identical to the 2D reference after every step to the end of the run (298 decks, 182 of them synthetic, in `gate/decks/twods`). The run-time 1D-2D exchange (the junction exchange laws, the outfall tailwater and outfall water exchange, the coupled ponds, the delivery queue) is translated, gives upstream's bits on every exchange call of every coupled gate deck, and is wired into the step: every gate deck whose mesh drains into or spills from 1D junctions and outfalls (122 decks, among them upstream's coupled example) stays identical to the 2D reference after every step to the end of its run. The 2D output is next.

## How it is checked

The reference is OpenSWMM 6 itself, built here from the source below. Four kinds of gate:

| Gate | What must hold |
|---|---|
| G1 | `.out` byte-identical to the C++ engine on every gate deck |
| G2 | `.rpt` identical except the wall-clock lines |
| G3 | State after every routing step bit-identical |
| G4 | Leaf functions identical to upstream's own code, compiled from upstream source |

G1 and G2 are checked end to end through `openswmm` on 932 decks: 480 from the [1729-SWMM5-Models-2030](https://github.com/SWMMBobSWMM6/1729-SWMM5-Models-2030) corpus, 336 of upstream's own test decks, and 116 synthetic decks in `gate/decks`. The state gates check what a float32 `.out` can hide (G3): each opens a deck in both engines, dumps every member of the engine's context and solvers (doubles as hex bits) after open, initialize, start, each helper, each runoff step and each routing step, and compares the dumps and the exit codes; the routing-step gate loads upstream's own state before every call. Synthetic decks added since the 932 were fixed (`gate/decks/fv`, `heat`, `richards`) reach paths the 932 do not, and are checked the same way. The tables are in PLAN.md.

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
| `UPSTREAM_NOTES.md` | upstream oddities found along the way, each quoting the code |

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
