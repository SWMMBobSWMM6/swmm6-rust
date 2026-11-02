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
* **Errors as exceptions.** Upstream throws (`std::invalid_argument`, `length_error`) where the translation returns `Err` or treats a negative count as empty. This only matters on paths upstream itself fails on. The one known path is U1, which aborts upstream.

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

| Section handlers (all 1D sections but LID node, 2D coupling and REACTION_*) | input/handlers/*.cpp, DefaultInputPlugin's table | **parse gate: 572 decks (corpus + every upstream test deck), every store of the context identical** (`gate/ctxgate.sh`) |
| InputReader, SectionRegistry, SectionParser, Constants | the matching files | inside the parse gate |
| C-library input calls (strtod, atof, atoi) | glibc | 486,499 strings identical to glibc |

Rust so far: see `git log`; about 25,000 lines in `crates/engine`.

## Next

1. TableData.cpp and MultiColumnSeriesFile (external series files) — in progress; then swap the stub in catchment_handler.rs.
2. Cross-section geometry core (xsect_tables, XSectKernels, XSection, XSectBatch) — in progress.
3. The remaining section readers: LID node sections (hydrology/LidNode.cpp readers), 2D_NODE_COUPLING (core/NodeCoupling2D.cpp), REACTION_* (reactionSectionTags).
4. PostParseResolver (2,500 code lines; needs Transect.cpp, Street, ForceMain, Link, the geometry core, LidNode::validate, VirtualJunctionOps, MsxInitialQuality, PathResolver). Gate: extend ctxread to run upstream's resolver from the reference library after parsing, and the translation's after its parse.
