# swmm6-rust: plan and gates

The OpenSWMM 6 engine's 1D path, translated into safe Rust and held byte for byte against OpenSWMM 6 itself. Phase 1 is 1D; 2D comes later.

## The reference

* Upstream is HydroCouple/openswmm.engine, branch `swmm6_rel`, commit `88329d6df2f8dd63fe7e7c601954ec2bde3573e6` (6.0.0-alpha.4, 31 Oct 2026). It is Apache-2.0 with a NOTICE (see `NOTICE`, `LICENSE`).
* The reference build is Release, 1D only: `-DOPENSWMM_BUILD_2D=OFF -DOPENSWMM_BUILD_GPU_PLUGIN=OFF -DOPENSWMM_WITH_HDF5=OFF -DOPENSWMM_WITH_GEOPACKAGE=OFF`, built with g++ 13.3.0 on Ubuntu 24.04.
  * The fast kernels are off (the shipping default, documented as bit-exact).
  * OpenMP is on.
* Measured:
  * All 27 non-mesh decks of upstream's parity corpus run, and each gives the same `.out` on 1 thread as on 2, and on repeated runs.
  * The 2 decks that need a 2D mesh stop with ERROR 145, as they should without 2D.
* Upstream states that THREADS defaults to 1 and that the dynamic-wave solver is bit-identical at any thread count. No floating-point reduction crosses threads on the 1D path.

## Gates

| Gate | What must hold | Tool |
|---|---|---|
| G1 | `.out` byte-identical to native `openswmm` on every gate deck | `cmp` |
| G2 | `.rpt` identical except the lines that print wall-clock time | diff with those lines masked |
| G3 | State after every routing step bit-identical: subcatchment runoff, node depth and volume, link flow and depth | `gate/steptrace` (runs OpenSWMM 6 through its own C API) + `gate/trcmp.py` |
| G4 | Leaf functions identical to upstream's own code, compiled from upstream source, over every line or value in the corpus | `gate/diff/*` |

G3 points at the first step and object where a translation departs. G4 catches a departure before any engine runs.

**State gates** (`gate/ctxgate.sh <mode>`; build both sides first with `gate/build.sh`). Each runs the 688 decks (236 corpus, 336 upstream test decks, 116 synthetic in `gate/decks`) through upstream (the reference library, `gate/diff/*`) and the translation (`target/release/examples/*`), dumps the state and compares stdout and exit code. Every mode is identical on every deck (2 Nov):

| Mode | Upstream does | Dump | Decks identical |
|---|---|---|---|
| `ALL`, `ALL resolve` | InputReader + handlers; then resolve_cross_references | context | 688, 688 |
| `open`, `open lenient`, `open metadata` | swmm_engine_open (lenient; open_metadata) | context + open() members | 688, 688, 688 |
| `inithyd` | open, initHydraulics | + Router and the hydraulic solvers (`rtdump`) | 688 |
| `inithydro`, `inithydro lenient`, `inithydro rdii` | open, initHydrology (rdii: rdii_.init first) | + hydrology solvers (`hydrodump`) | 688, 688, 688 |
| `initqual quality`, `geometry`, `massbal` | open, initQuality / initGeometry / initMassBalance | + quality members (`qualdump`) | 688, 688, 688 (massbal skipped on 6 on both sides: upstream's own initialize stops first; two of the 6 are the U57 terminate decks, where both sides abort with exit 134 and empty stdout) |
| `loadrt open`, `lenient`, `metadata` | dump → `examples/ctxload` → dump again, from the translation's dump and from upstream's | context + open() members | 688 rust and 688 cross, each |
| `loadrt init` | the same from upstream's dump after swmm_engine_initialize (cross only) | context + open() members | 685 (no upstream dump on 3, where upstream terminates inside initialize: IO-41 (U1) and the two U57 decks) |
| `init`, `start` | swmm_engine_open, swmm_engine_initialize (and swmm_engine_start(1)) in `gate/diff/ctxstep` against `examples/ctxstep` | everything: context, open() members, `rtdump`, `hydrodump`, `qualdump` and the lifecycle members (`lifedump`) | 932 and 932 on decks932 (480 + 452); 20 and 20 on decks20 |
| `open init` | ctxopen after swmm_engine_initialize | context + open() members | 932 (decks932) |
| `loadrt init` | the loader round trip after initialize, both sides | context + open() members | rust 930, cross 930 (decks932; no dump on the 2 U57 decks, where both sides abort) |
| `step [N]` | the same as `start`, then up to N swmm_engine_step calls; each step's dump hashed | the same | N = 3: decks480 49 identical + 431 stopped at a stub; decks452 124 + 328; decks20 4 + 16; 0 differ (see "Run lifecycle" below) |
| `helpers [K]` | upstream runs K real steps, then calls each of step()'s helpers directly; the translation loads that state and calls its translation | context after each helper (hashed) | K = 10: 480 (decks480), 452 (decks452), 20 (decks20) |

**Routing-step cross test** (`gate/ctxgate.sh route [K1 [K0 [hyd]]]`, harness `gate/diff/routestep` and `examples/routestep`). Upstream runs a deck through its own C API; the harness interposes on SWMMEngine::stepRouting (exported with `-Wl,-E`, the library's call goes through its PLT) and, for calls K0..K1, dumps the whole state (context, rtdump, qualdump and the engine members stepRouting uses) before and after the library's stepRouting. The translation loads each BEFORE dump into its own engine (built by the translated open and init legs), checks the loader's round trip, runs `step_routing` with upstream's dt and compares with the AFTER dump. Nothing is written to disk (the dumps are piped). `hyd` sets ignore_quality on both sides for the call. A deck's [FILES] SAVE files are backed up and restored around each run. The loader does not reopen interface files, so a deck with an open RDII or routing interface file shows as loaddiff here; `ctxgate.sh step` covers those once step() reaches stepRouting.

**Output-side gate** (`gate/ctxgate.sh out start | out run [K0 [K1 [STRIDE]]]`, harness `gate/diff/outgate` and `examples/outgate`; the writers without the runoff step). `out start`: both sides run open, initialize and start(save = 1) as the CLI does and stop there, flushing their files (a failure takes the CLI's failure path: close, destroy); the two `.rpt` files without the three clock lines, the two `.out` files and the exit codes must be identical. That is the preamble and the `.out` header of every deck that starts, and the whole report of every deck that stops in open, initialize or start. `out run`: upstream runs the deck to the end through its C API; the harness defines SWMMEngine::postOutputSnapshot (the library's step() calls it through its PLT, as routestep does stepRouting) and prints the whole state (context, rtdump, qualdump and the engine members the writers read: the runoff clock, old_gw_rate_, the averaging accumulators, the groundwater, snow and LID solver members, the period count) before each report period it writes, and again before end(). The translation loads each state into its own started engine and runs the translated post_output_snapshot, then end, report and close; it compares every period record it wrote with upstream's record of the same period, the `.out` header and closing records, the `.rpt` without its clock lines and the exit code. STRIDE -n dumps about n periods spread over the run (the end of the run is always dumped). Results: identical P (periods compared), differ (what, first differing period, byte or line), nottranslated, nostart (upstream stops before stepping; the files are then compared as in `out start`), noresult, panic. `gate/withsaves.sh DECK CMD...` runs a command on a corpus deck with its SAVE files backed up and restored; `gate/openswmm_saves.sh` is the CLI for `gate.py` with the same protection.

The dumpers are generated twice, C++ and Rust, by `tools/ctxdump.py` (`gate/diff/{ctxdump,rtdump,hydrodump}.hpp`, `examples/common/{ctxdump,rtdump,hydrodump}.rs`); the loader `examples/common/ctxload.rs` by `tools/ctxload.py`, which checks itself against ctxdump.py's output. The context dump covers every `SimulationContext` member that is data. Left out: `wall_start` (the wall clock at open), `state_accessors` (closures over the engine's solvers), `twod_io` and `spatial.coupling_points` (pointers; null or empty in a 1D build) and `Table::file_handle` (a `FILE*`). Inside the engine dumps: `lid_.rpt_files_` (open streams, printed as a count), the FV quality mixers' `gmres_` workspaces, and the RDII UH parameters and name map, which are printed by hand.

The step-trace harness was checked against the CLI on `orif_legacy.inp`: running the engine through the C API writes the same `.out` as the `openswmm` command.

**Gate decks:**
* the 25 one-dimensional decks of upstream's parity corpus (the 27, less the two `twod_*` decks, whose 2D sections a 1D build skips);
* the SWMM5 corpus decks that OpenSWMM 6 runs without error (see Findings).

## Findings so far

1. **An OpenSWMM 6 `.out` differs from an EPA 5.2.4 `.out` in seven header bytes on every file:**
   * the version word (60000 against 52004);
   * the input-property codes. OpenSWMM's enum adds OUTPUT_INVERT, OUTPUT_MAX_DEPTH and OUTPUT_OFFSET, so node codes are 0,2,4 rather than 0,2,3, and link codes are 0,6,7,4,8 rather than 0,4,4,3,5.
   
   The values after the codes are the same quantities in the same order. `gate/outnorm.py` rewrites the seven bytes; after it, `orif_legacy.inp` is byte-identical between the two engines.
2. **OpenSWMM 6 is not EPA 5.2.4 with a new header.** The SWMM5 corpus census ran 221 models through both engines, with the `.out` header normalized:

   | Outcome | Models |
   |---|---|
   | both run, results byte-identical | 64 |
   | both run, results differ | 112 |
   | only OpenSWMM 6 runs (EPA 5.2.4 stops with an error) | 22 |
   | only EPA 5.2.4 runs | 4 |
   | neither runs | 19 |

   Of the 4 models only EPA runs, 3 are inputs OpenSWMM 6 now rejects on purpose (code-review decks NUM-35, CON-12, API-18). The fourth is a crash (finding 4).
   
   198 models run in OpenSWMM 6 and are gate candidates. OpenSWMM's legacy baseline is 5.2.4 plus its own "Build 5.3.0" fixes in 23 legacy files (dynwave, xsect, infil, node, routing, lid and others), and the new engine is held to that baseline, not to EPA. So the reference is OpenSWMM 6 itself, and swmm5-rust cannot stand in for it.
3. **The `.inp` tokenizer is not SWMM5's.** A comma separates fields on a line with no white space, or where it touches white space; between two other characters it is part of a name. A quote groups text inside a token. Empty CSV fields are kept. Numbers are read with `std::from_chars`, which rejects a leading `+`, accepts `nan(...)`, and treats `2e-324` as out of range (so not a number) while `2.5e-324` is the smallest subnormal.
4. **OpenSWMM 6 aborts on `IO-41_point.inp`** (`[CONTROLS] EXPRESSION E1 = 2 + .`), with "terminate called after throwing an instance of 'std::invalid_argument' what(): stod". In `src/engine/math/MathExpr.cpp`, the scanner takes a lone `.` as the start of a number and converts it unguarded:

   ```cpp
           if (std::isdigit(c) || c == '.') {
   ...
               t.value = std::stod(s.substr(start, i - start));
   ```

   `ControlEngine::addExpression` (`Controls.cpp:738`) calls `mathexpr::parse` with no `try`, so the exception ends the process. EPA 5.2.4 runs the deck. The translation must reproduce upstream's behavior up to this point, and must not abort.

## Translation hazards found so far

* **Iterating an `unordered_map`.** Wherever upstream loops over a `std::unordered_map` and the order reaches output, the line order is libstdc++'s internal hash-table order. Two cases so far:
  * `InpWriter.cpp` writes `[USER_FLAG_VALUES]` by looping over `ctx.user_flags.all_values()`. The `.inp` writer is off the run path, so this is not urgent.
  * `SectionRegistry::registered_tags` iterates two maps; whether anything uses that order has not been checked yet.
  
  Each such loop on the run path must either be shown not to reach a result, or reproduce libstdc++'s `std::hash<std::string>` (`_Hash_bytes`, seed 0xc70f6907) and its bucket growth.
* **NaN sign when both operands are NaN.** x86 `mulsd` returns the first operand's NaN. g++ keeps source order; LLVM may commute a multiply. When both operands of `a * b` are NaN with different signs, the result's sign bit can differ. The geometry differential showed exactly this: 47 of 3,222,094 values, all from degenerate (NaN) dimensions. Safe Rust cannot pin the operand order. It matters only when a NaN would reach the output, which a valid model does not produce.
* **Errors as exceptions.** Upstream throws (`std::invalid_argument`, `length_error`) where the translation returns `Err` or treats a negative count as empty. This only matters on paths upstream itself fails on. Where an exception escapes a `noexcept` function upstream and the process terminates (U1/U65 in initHydraulics, U57 in initQuality), the translated function returns `Err(Terminated { exception, what })` at the same point, with upstream's exception type and `what()`. The library never aborts; the gate harnesses (`examples/inithyd.rs`, `examples/initqual.rs`) print libstdc++'s terminate message and abort, so the gate compares exit code 134 and the empty stdout on both sides (decks `gate/decks/initqual/initq_treat_dot.inp` and `initq_treat_huge.inp`).

## Size

Counted by the module survey: lines a SWMM5-equivalent deck runs through, excluding the C API, 2D, GeoPackage, checkpoints and the `.inp` writer.

| Module | Upstream C++ on the run path | Counterpart in SWMM 5 |
|---|---|---|
| core (SWMMEngine, SimulationContext, options, dates, units, hotstart, interface files) | ~26k | swmm5.c, routing.c, massbal.c, stats.c, project.c, … |
| data (structure-of-arrays stores, tables, name index) | ~10.9k | objects.h, table.c, hash.c |
| input (reader, tokenizer, handlers, post-parse resolver) | ~14.6k | input.c, keywords.c, project.c |
| plugins (default .inp, .out and .rpt plugins) | ~5.3k | output.c, report.c, statsrpt.c, inputrpt.c |
| hydrology (runoff, infiltration, gage, climate, groundwater, LID, RDII, snow) | ~13.3k | subcatch.c, infil.c, gage.c, climate.c, gwater.c, lid*.c, rdii.c, snow.c |
| hydraulics (DW, KW, cross sections, structures, culvert, force main, street, inlet, transect) | ~27.4k | dynwave.c, dwflow.c, kinwave.c, xsect.c, link.c, node.c, … |
| quality, controls, math | ~9.5k | qualrout.c, treatmnt.c, landuse.c, controls.c, findroot.c, odesolve.c, mathexpr.c |
| **SWMM5-equivalent path** | **~107k** | |
| SWMM6-only 1D: FV router 17.4k, EULERIAN_ARD / LAGRANGIAN quality, water age, heat, UEB snow, Richards LID | ~45k | none |

## Phases

Each phase ends on a gate run. A deck enters the gate set when every section it uses is translated.

| Phase | Translates | Gate at the end |
|---|---|---|
| **1a** foundations | from_chars, tokenizer ✔; SimulationOptions, NameIndex, TableData, node/link/subcatchment stores, SimulationContext; InputReader, SectionRegistry, the handlers and PostParseResolver for the sections hydraulics needs | G4 on the parser; parsed state compared through the C API getters |
| **1b** hydraulics | cross sections (all shapes, transects, custom), DynamicWave, Node/Outfall/Divider/Storage, KinematicWave, Steady, orifice, weir, outlet, pump, culvert, force main, roadway, street, inlet, exfiltration, TopoSort, controls, mass balance, statistics, the `.out` and `.rpt` writers, the CLI | G1–G3 on the hydraulics-only decks (EXTRAN shapes, user1, force main, pumps) |
| **1c** hydrology | gages, rain files, climate, runoff, infiltration (all five), groundwater, degree-day snow, RDII, LID | G1–G3 on every SWMM5-feature deck that runs |
| **1d** quality | pollutants, land uses, buildup and washoff, treatment, the LEGACY quality router | the same, with quality |
| **1e** SWMM 6 only | the FV router (HLL/HLLC, MUSCL, local time stepping, implicit GMRES), DYNAMIC_SLOT, TPA, Anderson acceleration, virtual junctions, EULERIAN_ARD, LAGRANGIAN, water age, heat, UEB snow, Richards LID, the OPENSWMM_HS hotstart | the 25-deck parity corpus |
| 2 | 2D | upstream's 2D bitwise regression |
| 3 | swmm6-rust in the RustGIS plugin and the page, beside swmm5-rust | the plugin against native OpenSWMM 6 |

## Rules for the translation

* Translate the C++ as it stands, file by file, keeping its structure-of-arrays layout (`SimulationContext` is the one context passed to every solver). Do not "improve" arithmetic: the order of floating-point operations is the result.
* Bytes, not text. Names and file contents are byte strings, as `std::string` is.
* The math library is a hazard, as it was in swmm5-rust: Rust's `f64::powf`/`exp` must give glibc's bits on the gate machine. G3 shows any such departure at its first step.
* Each translated file:
  * names the upstream file or files it derives from;
  * keeps their copyright line and license;
  * states that it was modified (translated), when, and by whom.
* No endorsement is implied. This is not "SWMM6", and the crate does not call itself OpenSWMM.

## Done

| Piece | Upstream | Check (identical to upstream compiled from source unless noted) |
|---|---|---|
| Reference build, corpus reproducibility, census | — | 27 parity decks repeatable on 1 and 2 threads; 221-model census |
| Gate tools | — | `gate/steptrace.c` + `trcmp.py` (C API run == CLI run on orif_legacy); `outnorm.py` |
| `core::charconv` | charconv_compat.hpp + std::from_chars | 30 double and 20 integer cases recorded from libstdc++ |
| `core::cmath` | std::min/max/clamp | NaN and signed-zero cases |
| `core::datetime` | DateTime.hpp | 200,002 lines |
| `core::error_codes` | ErrorCodes.hpp/.cpp | 85,456 lines (every code, every formatter) |
| `core::unit_conversion` | UnitConversion.hpp/.cpp | 178 lines (every factor, as bits) |
| `core::options`, `spatial_frame`, `user_flags` | SimulationOptions, FvOptions, FilePathPair, SpatialFrame, UserFlags | compiles; defaults transcribed |
| `input::tokenizer` | Tokenizer.cpp | 584,895 corpus lines |
| `input::parse_utils` | InputParseUtils.hpp | 447,009 lines (every distinct corpus token) |
| `data::table` (in-memory) | TableData.hpp | 134,703 query rows, all ten lookups |
| `data::name_index` | NameIndex.hpp, StringCase.hpp | unit test |
| Node, link, subcatchment, gage, pollutant stores | NodeData, LinkData, SubcatchData, GageData, PollutantData | sizing after every method in a sequence: 3,196 lines |
| Link and node subtypes, inlet usage | LinkSubtypes, NodeSubtypes, InfraData | row bookkeeping through build/retype/erase: 360 + 386 lines |
| Storage geometry | StorageGeometry.hpp | 60,009 lines |
| Inflow, infra, LID node, quality, treatment types, 2D coupling data | InflowData, InfraData, LidNodeData, QualityData, Treatment.hpp, NodeCoupling2DData | compiles; sizing checked by the agents' scripts |
| SWMM6-only data: heat, water age, ARD, bed zone, channel thermal, LID layer species, MSX surface, reactions, species registry, initial quality, forcing | the 13 headers | compiles |

| Section handlers (every section DefaultInputPlugin registers in a 1D build) | input/handlers/*.cpp, DefaultInputPlugin's table | **parse gate: 681 decks (corpus, every upstream test deck, 109 synthetic), the whole context identical** (`gate/ctxgate.sh ALL`) |
| SWMMEngine class and open(): plugins, process components, reactions component (apply + expression compiler), heat, water age, ARD config, surface species resolve, gage setup | SWMMEngine.cpp open() and what it reaches | **open gate: 681 decks identical after swmm_engine_open in strict, lenient and metadata modes: the engine's own context and members, open codes, error text** (`gate/ctxgate.sh open [lenient\|metadata]`) |
| PostParseResolver (cross references, unit conversion, rain/series file loading incl. NWS formats), PathResolver, VirtualJunctionOps (what the resolver reaches), MsxInitialQuality | PostParseResolver.cpp and the rest | **resolve gate: 681 decks identical after resolve_cross_references** (`gate/ctxgate.sh ALL resolve`); path resolution 61,560 pairs; `%.10g` 20,040 doubles; sscanf 30,016 lines × 24 formats vs glibc |
| SWMMEngine::initHydraulics: Router::init (XSectGroups, DW/KW/STEADY setup, FV mesh builder and explicit solver initialize, outfall map), TopoSort, layout checks, thread setup (libgomp model), inflow, RDII, structure, exfiltration, inlet and interface-file solvers, event sort, [CONTROLS] rule parser; libstdc++ `unordered_map` iteration order | SWMMEngine.cpp initHydraulics, Routing, DynamicWave, KinematicWave, fv/*, HydStructures, Outfall, TopoSort, Inflow, RDII, Exfiltration, Inlet, InterfaceFile, Controls, MathExpr | **inithyd gate: 681 decks** (`ctxgate.sh inithyd`); unordered_map order: 200 insertion sequences, `std::hash<int>` and CiHash (`gate/diff/umaporder`) |
| SWMMEngine::initHydrology (+ setNextEvapDate): runoff and infiltration (all methods, InfilBank), gage states, climate file, radiation, degree-day and native UEB snow (lumped, subcatchments, forcing), groundwater with [GWF] expressions, LID (with the heat and water-age LID layers) | SWMMEngine.cpp initHydrology, Runoff, Infiltration, Gage, ClimateFile, Radiation, DegreeDaySnow, Ueb*, Groundwater, LID, HeatLid, WaterAgeLid | **inithydro gate: 681 decks in open, lenient and rdii modes**, incl. 13 synthetic climate-file, UEB, GWF, LID and RDII decks (`gate/decks/inithydro`) |
| SWMMEngine::initQuality, initGeometry, initMassBalance: land uses, buildup/washoff params, the LEGACY quality router's FV transport state, msxlid, msxsurf::initState, ensureMsxState, treatment cycle check, node volume and surface area | SWMMEngine.cpp, Landuse, QualityRouting, FvQualityTransport, ConservativeMixing, MsxLid, MsxSurface, ReactionLegacyBinding, Node (getVolume, getSurfArea), LidNode | **initqual gate: 686 decks (+5 init-leg decks) in quality, geometry and massbal modes** (`ctxgate.sh initqual …`; massbal installs upstream's post-initialize hydrology state into the translation's solvers) |
| SwmmEngine members of all three init legs, one struct in upstream's order | SWMMEngine.hpp | inside the init gates |
| Context loader: reads a dump (either side's) into a fresh SimulationContext and SwmmEngine | — (`tools/ctxload.py` → `examples/common/ctxload.rs`; `gate` feature for the private members) | **loadrt: 681 decks round-trip byte for byte from both sides' dumps after open (strict, lenient, metadata); 680 from upstream's dumps after swmm_engine_initialize** (`ctxgate.sh loadrt …`) |
| Context dump coverage: every SimulationContext member that is data (run clocks, climate, water-age, heat, bed and LID-layer state, forcing, coupling vectors, 2D mass balance, max_* stats, routing stats, VJ and inlet diagnostics, control log, every Table member and the by-name index, unordered maps sorted by key, user-flag index) | tools/ctxdump.py | every state gate above |
| InputReader, SectionRegistry, SectionParser, Constants | the matching files | inside the parse gate |
| C-library input calls (strtod, atof, atoi) | glibc | 486,499 strings identical to glibc |
| External series files | TableData.cpp, MultiColumnSeriesFile | 138 files (82 real, 53 adversarial, 3 error cases), 62,175 lines |
| Transect, street, force main, link, culvert, FindRoot | Transect.cpp, Street.cpp, ForceMain.cpp, HydClosureKernels.hpp, Link.cpp, Culvert.cpp, FindRoot | 7,629,922 values (every corpus transect, street and shape curve; every culvert code): identical but for the NaN-sign hazard (13 lines) |
| Cross-section geometry core | xsect_tables.hpp, XSectKernels.hpp, XSection.cpp, XSectBatch | 3,222,094 values from 338 parameter sets × every function, plus the batch and kernel APIs: identical but for the NaN-sign hazard (47 values) |

| Run lifecycle (5 Nov): SWMMEngine::initialize (state checks; state seeding: MIN_SURFAREA junction storage, LID nodes incl. Richards, initial node volumes, pump wet wells, the legacy DW volume convention, q0 normal depths, initNodeDepths, outfall depths, backwater fill; legacy USE HOTSTART with its runoff section and species trailer; init_modules in upstream's order; forcing.resize; the node inflow seed) and start (routing/runoff/RDII interface files, the RAINFALL warning, the plugin lifecycle and prepare, averaging, NodeInflow seed, SAVE HOTSTART flags) | SWMMEngine.cpp initialize, init_modules, applyLegacyRunoffSection, computeFinalQualityMassBalance, start; HotStartManager (legacy file), LidNode initialize, RichardsColumn theta/storage/pressure, Node getDepth, Outfall setAllOutfallDepths, InterfaceFile open/headers, RunoffInterface and RdiiInterface open, PluginFactory initialize/validate/prepare_all, InfilBank pack/unpack, msxsurf remainingBuildup | **init 932/932 and start 932/932 identical on decks932** (`ctxgate.sh init`, `start`; also `open init` 932, `loadrt init` 930 + 2 aborting on both sides). Reached: on decks480, 431 decks start (7 stop in start with SWMM_ERR_IO, 17 in initialize, 25 in open); on decks452, 328 start |
| SWMMEngine::step and its own helpers: the step-size choice (between events, adaptive, rule grid), applyForcings, resetStepMassBalance, accumulateNodeRoutingTotals, checkRoutingDiverged, updateStatistics, updateRoutingMassBalance (with inlet adjustFloodingTotals, bookMsxSystemOutflow), computeFinalStorage, computeFinalQualityMassBalance, the routing-report split, the clocks (TimestepController), advanceTo; stepRouting's helpers assembleLateralInflows, isBetweenEvents, isInSteadyState | SWMMEngine.cpp, TimestepController.cpp, Inlet.cpp, MsxLegacyTransport.cpp | **helpers gate: 932/932 identical (K = 10 steps of upstream's state)**; step() itself: `ctxgate.sh step 3` stops at the first stub on every deck that starts (Router::getAdaptiveStep on 598, SWMMEngine::stepRunoff on 161 of decks932) with everything before it identical, 0 differ |
| SWMMEngine::stepRouting and what it calls, for STEADY, KINWAVE and DYNWAVE (5 Nov): events and pump hysteresis, control rules (evaluate, setting transitions), inflows (external, DWF, RDII grid and USE RDII file, interface file), MSX coupling queue, the steady-state skip, the DW solver step (Picard loop run as a team of one, EXTRAN and SLOT surcharge, virtual junctions, Anderson, bypass), KW and steady flow, node, outfall, divider, storage and storage exfiltration (Green-Ampt), pumps, orifices, weirs, outlets, culverts, force mains, roadways, street and channel inlets (HEC-22), node LID columns (not Richards), getAdaptiveStep, accumulateLinkStepStats | SWMMEngine.cpp stepRouting; Routing.cpp; DynamicWave.cpp; KinematicWave.cpp; Node.cpp, Outfall.cpp, Divider, Roadway; HydStructures.cpp, Culvert.cpp, ForceMain.cpp; Exfiltration.cpp, Infiltration.cpp (grnampt_getInfil); Inlet.cpp computeAll; Controls.cpp evaluation; Inflow.cpp computeAll; LidNode.cpp run time; InterfaceFile.cpp readInflows | **route gate, decks932, calls 1..300: 665 identical, 0 differ**; the rest: 170 nostart and 2 nosteps (upstream itself does not step), 92 not translated (71 quality routing, 16 FV, 4 species surface washoff, 1 Richards LID column), 3 loaddiff (USE/SAVE RDII: the loader cannot reopen the RDII interface file). Deeper windows: calls 2001..2100 on decks932, 319 identical, 0 differ; calls 5001..5050 on the 887 decks without [FILES] side files, 220 identical, 0 differ. With ignore_quality set on both sides for the call (`route … hyd`), all 75 quality and species-washoff decks are identical over 300 calls; node LID columns: 6 decks (lid_chain_mass_balance, lid_node_resaturation, lid_partition_backwater, lid_partition_free, initq_rxlid_nodes, initq_rxlid_units) identical in hyd mode; virtual junctions: street_inlet_junction, cedar_glen_dw, cedar_glen_dw_extran identical over 300 calls. `ctxgate.sh step 100` on decks932: 0 differ, 172 identical (both sides stop in open, initialize or start), 760 stop at the stepRunoff stub in step 1, so step() does not reach stepRouting yet. Not gated: InterfaceManager::readInflows (every USE INFLOWS deck stops before it on both sides) |

Rust so far: see `git log`; about 59,100 lines in `crates/engine/src`, 47,200 without comments and blank lines (2 Nov).

## Next

Approach for the solver (phase 1b on): the reference library exports every C++ symbol, and the C API handle is the SWMMEngine, so a harness can open, initialize and step upstream and dump its whole context at any point (`gate/diff/ctxopen`). Two uses:
* **State gate.** Dump after open, after initialize and after every step on both sides, and report the first step and member that differ (G3, with every member rather than five).
* **Context loader (built).** `tools/ctxload.py` generates a Rust loader from ctxdump.py's lists that reads a dump into a SimulationContext and SwmmEngine. Load upstream's state at step k into Rust, run one translated step, and compare with upstream's dump at k+1. This tests each solver piece on real states in isolation, without waiting for the whole engine. It loads the context and the open() members now. A step also reads solver state held outside ctx in SWMMEngine members (DW tiles, router state, the hydrology solvers): `rtdump` and `hydrodump` print those, and the loader has to learn to read them back before the per-step test.

**The step's signatures (agreed for the merge).** step() (core/engine_step.rs) calls the solvers by these names; each lives in a file of its own that the agent translating it replaces wholesale. Until then the file is a STUB:

| File | Signature | Stub |
|---|---|---|
| `core/engine_step_routing.rs` | `SwmmEngine::step_routing(&mut self, dt_routing: f64) -> Result<(), Terminated>` | translated (rs6-route); paths not translated yet return `Err(Terminated)` with `exception == STUB_EXCEPTION` and `what` naming the path |
| `hydraulics/routing_step.rs`, `hydraulics/inlet_step.rs` | `Router::get_adaptive_step(&mut self, ctx: &mut SimulationContext, fixed_step: f64, courant: f64) -> f64`; `Router::accumulate_link_step_stats(&mut self, ctx: &mut SimulationContext, dt: f64, in_window: bool)`; `InletSolver::compute_all(&mut self, ctx: &mut SimulationContext, dt: f64)` | translated (rs6-route) |
| `core/engine_step_runoff.rs` | `SwmmEngine::step_runoff(&mut self, dt_routing: f64) -> Result<(), Terminated>` (and step_groundwater, step_surface_quality, which stepRunoff calls) | `Err(stub(...))` |
| `core/engine_output.rs` | `SwmmEngine::post_output_snapshot(&mut self, dt_step: f64) -> Result<(), Terminated>`; `accumulate_avg_results(&mut self) -> Result<(), Terminated>`; `write_due_hotstart_saves(&mut self) -> Result<(), Terminated>`; `DefaultOutputPlugin::write_header(&mut self, ctx: &SimulationContext)` (core/output_plugin.rs); `DefaultReportPlugin::write_preamble(&mut self, f: &mut CFile, ctx: &SimulationContext)` (core/report_plugin.rs) | translated (rs6-output); a NATIVE SAVE HOTSTART row returns `Err(stub)` |

Error convention: `Err(Terminated)` is a termination (an exception escaping a noexcept function) and step() passes it up; `exception == engine_step::STUB_EXCEPTION` marks untranslated code; an exception still propagating to step()'s `catch (const std::exception&)` (from a function that is not noexcept, such as stepRunoff or stepRouting) is an `Err(Terminated)` whose `exception` starts with `engine_step::THROWN_PREFIX` ("throw "), which step() turns into set_error(SWMM_ERR_WRONG_STATE, what) and return 6, as upstream does. stepRouting's helpers assemble_lateral_inflows, is_between_events and is_in_steady_state are in core/engine_step.rs, pub(crate), for the routing step to call.

1. **The solvers behind the stubs** (stepRunoff; stepRouting and the output side are in). Upstream's step() calls stepRunoff on every deck, with or without subcatchments (the runoff clock advances regardless), so until stepRunoff is translated `ctxgate.sh step N` stops there on every deck that starts and never reaches the routing step; the route cross test is what checks stepRouting meanwhile. When a stub file is replaced, `ctxgate.sh step N` goes on to the next stub; with all of them in, it is the per-step state gate (G3 with every member), and with CTXGATE_DIFFDIR set it writes both sides' full dumps at the first differing step.
2. **What stepRouting still refuses** (Err(Terminated) with STUB_EXCEPTION, before any state changes): the quality stage (QualitySolver::execute and assembleExternalLoads, Treatment, the EULERIAN_ARD and LAGRANGIAN engines, InletSolver::adjustQualInflows; 71 decks of decks932 stop here, and all of them are identical in the route gate's hyd mode), the FV router (16 decks), species surface washoff delivery (msxsurf::deliver, 4 decks), Richards node LID columns (LidNodeRichards.cpp; 1 deck), engine checkpointing, and inside the DW solver DYNAMIC_SLOT, TPA, physical junction areas (NodeData::storage_area > 0) and 2D-coupled ponds; control-rule expression premises and checkpoint sessions in the controls. No deck of decks932 reaches DYNAMIC_SLOT, TPA or physical junction areas, so those need synthetic decks.
3. **Not translated in the lifecycle** (each reached only by decks no gate set holds, and each returns Err(stub) where it is met): a native OPENSWMM_HS_V1 USE HOTSTART file (HotStartManager::open/apply, the engine checkpoint), prepareNativeHotstartCheckpointing (SAVE HOTSTART under coupled heat), engine checkpointing in step(), the FV solver's divergence check in checkRoutingDiverged, and a NATIVE SAVE HOTSTART row (HotStartManager::save and the engine checkpoint; a LEGACY or AUTO row of a deck without storage LIDs or heat is written as legacy `.hsf`). applyForcings' ADD rainfall forcing now reads the translated gage::getReportRainfall (hydrology/gage_report.rs).
4. **What the gates do not see yet.** The helpers gate skips assembleLateralInflows on decks with street inlets (its last call, InletSolver::computeAll, is the routing step's; it is translated now, and the route gate covers it) and checkRoutingDiverged under FV (the FV solver's own check). Forcings are set only through the C API, so on gate decks applyForcings has nothing to apply. step()'s own step-size choice (between events, the rule grid, compute_next, IGNORE_ROUTING) and the order of its calls are seen only by `ctxgate.sh step`, which stops at the first stub; they become visible when the stubs are replaced.
5. **Corpus side files.** start() opens a deck's [FILES] SAVE files for writing in place: `Hydrology/save_rtk_long.inp` and `one_sub_one_rdii_node.inp` both SAVE RDII `cho.txt`, which `Use_SeasonalRDIIn.inp` USEs. Any harness that runs start() over the corpus must back these up and put them back after each run, as `ctxgate.sh init|start|step|helpers|route` do (`gate/savefiles.py` lists them), or a killed run leaves `cho.txt` empty and the USE deck then fails start (SWMM_ERR_IO) on both sides. With a full file (the reference CLI's SAVE from save_rtk_long.inp, 1,114,474 bytes) in a scratch copy of the USE deck, start and helpers are identical and the RDII binary reader holds its first record on both sides.
