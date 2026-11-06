# swmm6-rust

The OpenSWMM 6 engine's 1D path, translated from C++ into safe Rust and checked byte for byte against OpenSWMM 6 itself.

The goal is a Rust engine that gives the same `.out` file, bit for bit, as the C++ engine on the same deck. It is not a rewrite or a new solver. Every floating-point operation keeps upstream's order, so the same inputs give the same bits. Where upstream does something odd, the translation does the same odd thing, and the oddity is written down in [UPSTREAM_NOTES.md](UPSTREAM_NOTES.md).

## Status

Phase 1a (1D) is under way. Done, and identical to upstream on every gate deck:

* reading the `.inp` (every section handler), cross-reference resolution, and `open()` in strict, lenient and metadata modes;
* the initialize legs: `initHydraulics`, `initHydrology`, `initQuality`, `initGeometry`, `initMassBalance`;
* the leaf numerics: cross-section geometry and its tables, transects, streets, force mains, the C library functions upstream depends on (`strtod`, `from_chars`, `sscanf`, `%.10g`, `std::sort` and the `unordered_map`/`map` orders);
* the run lifecycle: `initialize()` as a whole and `start()` (identical on all 932 gate decks), and `step()` with its own bookkeeping (statistics, routing mass balance, final storage), gated helper by helper.

Not yet translated: the solvers the step calls (runoff, routing, quality), nor the `.out` and `.rpt` writers; `step()` calls them through stub files. **It does not run a simulation yet.** [PLAN.md](PLAN.md) has the phases and what comes next. 2D comes after 1D.

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
| `crates/engine` | the library (`swmm6_rust`), `#![forbid(unsafe_code)]` |
| `crates/engine/examples` | the Rust side of each gate harness |
| `gate/` | the C++ harnesses that link the reference library, the deck set, `ctxgate.sh` |
| `tools/` | generators for the state dumpers and loader, and a first-pass translator for structure-of-arrays headers |
| `census/` | OpenSWMM 6 against EPA SWMM 5.2.4 on a corpus of SWMM5 models |
| `TRANSLATING.md` | the translation rules |
| `UPSTREAM_NOTES.md` | 77 upstream oddities found along the way, each quoting the code |

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
