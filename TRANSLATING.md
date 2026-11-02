# How a file is translated

These conventions apply to every file in `crates/engine`. `crates/engine/src/data/node.rs` is the worked example of a structure-of-arrays store, and `data/table.rs` of numeric code.

## The result is the arithmetic

* Keep every floating-point operation in upstream's order. `a + b + c` stays `(a + b) + c`; never reassociate, never fuse, never "simplify".
  * Write `x * x` where upstream writes `x * x`, and `pow` where it writes `pow`.
  * Rust does not contract `a * b + c` into an FMA, so leave that expression as it is.
* `std::min`/`std::max`/`std::clamp` become `core::cmath::{cmin, cmax, cclamp}`, never Rust's `f64::min`/`max`/`clamp`, which treat NaN differently.
* `std::fabs` becomes `.abs()`, and `std::sqrt` becomes `.sqrt()`; both are exact.
* Transcendental functions (`pow`, `exp`, `log`, `sin`, …) need glibc's bits. Mark each use with `// LIBM` so the later libm pass can find it.
* `static_cast<int>(double)` becomes `core::cmath::cvt_i32(x)` wherever the value can be NaN, infinite or out of range. The reference build gives INT_MIN for all of those, where Rust's `as i32` saturates (NaN gives 0). For values known to be in range, `as i32` is the same.
* Integer arithmetic keeps C++ widths: `int` is `i32`, `long` is `i64`, `std::size_t` is `usize`. A narrowing cast is `as`, as `static_cast` would be.

## Data

* Strings are bytes. `std::string` becomes `Bytes` (`Vec<u8>`, from `core::options`), and `std::string_view` becomes `&[u8]`. Never decode to UTF-8.
* An enum keeps upstream's `#[repr]` and discriminants. Variants are CamelCase (`OVERFLOW_DIV` becomes `OverflowDiv`).
* A field named after a Rust keyword gets a trailing underscore (`type` becomes `type_`).
* A structure-of-arrays store keeps one `Vec` per attribute, in upstream's field order.
  * Every sizing method (`resize`, `grow_to`, `erase_at`, `resize_quality`, …) sizes and clears exactly the fields upstream's does, with the same default values.
  * `v.assign(n, d)` becomes `v.clear(); v.resize(n, d)`.
  * `reserve` and `shrink_to_fit` are kept, though they cannot change a result.
  * A field upstream leaves unsized in a method stays unsized.
* `std::copy(a.begin(), a.end(), b.begin())` becomes `b[..a.len()].copy_from_slice(&a)`.
* Member initializers carry into `Default`. Derive it only when every initializer equals Rust's default; otherwise write `impl Default` by hand with upstream's values.
* `#![forbid(unsafe_code)]` holds for the whole crate.

## Where upstream throws

* An exception upstream catches becomes an `Err` the caller handles at the catch site.
* An exception that escapes a `noexcept` function (so upstream calls `std::terminate`) becomes `Err(Terminated { exception, what })` from that function, using `core::engine_init_hydraulics::Terminated` with upstream's exception type (`"std::invalid_argument"`) and `what()` (`"stod"`). Never `panic!` and never `std::process::abort` in the library; only a gate harness may abort, to reproduce the exit code.

## Each file's header

```rust
// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/<path> (HydroCouple/openswmm.engine,
// branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under
// the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, <date>, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
```

Copy the copyright line from the upstream file's own header. Add a short note on anything the translation does differently in form, never in result.

## Tools

* `tools/soa2rs.py X.hpp` gives a first pass for a structure-of-arrays header. It translates the per-field statements and leaves every other line as a `// MANUAL:` comment holding the C++ text.
  * Every MANUAL line must be resolved by hand.
  * It flattens `if` blocks, so check each per-field statement's scope against the C++.
* `gate/diff/` holds the upstream-source differentials: a C++ harness compiled from upstream's own file, a Rust example producing the same output, and `cmp`.
