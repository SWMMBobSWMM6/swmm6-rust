# swmm6-rust: plan and gates

The OpenSWMM 6 engine's 1D path, translated into safe Rust and held byte for byte against OpenSWMM 6 itself. Phase 1 is 1D; 2D comes later.

## The reference

* Upstream is HydroCouple/openswmm.engine, branch `swmm6_rel`, commit `88329d6df2f8dd63fe7e7c601954ec2bde3573e6` (6.0.0-alpha.4, 31 Oct 2026). It is Apache-2.0 with a NOTICE (see `NOTICE`, `LICENSE-APACHE`).
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

Rust so far: see `git log`; about 59,100 lines in `crates/engine/src`, 47,200 without comments and blank lines (2 Nov).

## Next

Approach for the solver (phase 1b on): the reference library exports every C++ symbol, and the C API handle is the SWMMEngine, so a harness can open, initialize and step upstream and dump its whole context at any point (`gate/diff/ctxopen`). Two uses:
* **State gate.** Dump after open, after initialize and after every step on both sides, and report the first step and member that differ (G3, with every member rather than five).
* **Context loader (built).** `tools/ctxload.py` generates a Rust loader from ctxdump.py's lists that reads a dump into a SimulationContext and SwmmEngine. Load upstream's state at step k into Rust, run one translated step, and compare with upstream's dump at k+1. This tests each solver piece on real states in isolation, without waiting for the whole engine. It loads the context and the open() members now. A step also reads solver state held outside ctx in SWMMEngine members (DW tiles, router state, the hydrology solvers): `rtdump` and `hydrodump` print those, and the loader has to learn to read them back before the per-step test.

1. **SWMMEngine::initialize itself.** The four init legs are in; what is left is initialize's own work around them:
   * the state checks and the state-seeding phase before init_modules (SWMMEngine.cpp from about line 1026 to 1569: node and link initial depths, volumes and flows, `lidnode::initialize`, the backwater and area seeding the comments there describe);
   * hotstart (USE HOTSTART: `hotstart_deferred_`, `hotstart_runoff_`, `hotstart_species_`, applyLegacyRunoffSection);
   * init_modules: the legs in upstream's order (initHydraulics, initHydrology, initQuality, the hotstart blocks, initGeometry, initMassBalance), then `ctx_.forcing.resize`.
   
   Gate: `ctxopen <deck> init` on both sides (a `ctxgate.sh open init` mode), which the loader's cross test already reads on 680 decks.
2. **start()**, then **step()**: the per-step state gate and the loader's step k → k+1 test, extended to the solver members.
