// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/ThreadInfo.hpp and ThreadInfo.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright (c) 2026
// Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Process-wide hardware / OpenMP thread capability query, Linux build with
// SWMM_USE_OPENMP (the reference build).
//
// Differences in form, never in result:
//   * The crate has no OpenMP runtime, so the libgomp state these functions read is
//     modelled, in safe Rust, the way libgomp and glibc compute it on Linux:
//       - logical_cpus (std::thread::hardware_concurrency, i.e. glibc get_nprocs) parses
//         /sys/devices/system/cpu/online, falling back to counting the cpuN lines of
//         /proc/stat, as glibc does;
//       - omp_max_threads (omp_get_max_threads) is the nthreads-var ICV: the value last
//         given to omp_set_num_threads in this process, else the first entry of
//         OMP_NUM_THREADS, else the size of the process's CPU affinity mask
//         (Cpus_allowed_list of /proc/self/status, what libgomp reads through
//         pthread_getaffinity_np);
//       - omp_team_size (what `#pragma omp parallel num_threads(n)` then reports through
//         omp_get_num_threads) is n capped at OMP_THREAD_LIMIT when that is set.
//     OMP_DYNAMIC and nested teams are not modelled (the reference harness sets neither).
//   * perfCores() is the Apple-only P-core query; on Linux it returns 0, and so does this.
//   * snprintf into a fixed buffer becomes format!; every message is far shorter than its
//     buffer, so no truncation applies.
//   * The Kokkos thread count is a process-wide atomic upstream and here.
//   * `const char* what` is `&str`; warnings are byte strings.

use std::sync::atomic::{AtomicI32, Ordering};

use crate::core::options::Bytes;

static G_KOKKOS_OMP_THREADS: AtomicI32 = AtomicI32::new(0);
/// libgomp's nthreads-var after omp_set_num_threads (0: never set in this process).
static G_OMP_NUM_THREADS_SET: AtomicI32 = AtomicI32::new(0);

/// Count the CPUs in a Linux cpu-list string ("0-3,8,10-11").
fn count_cpu_list(s: &str) -> i32 {
    let mut n = 0i32;
    for part in s.trim().split(',') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        if let Some((a, b)) = part.split_once('-') {
            if let (Ok(a), Ok(b)) = (a.trim().parse::<i32>(), b.trim().parse::<i32>()) {
                if b >= a {
                    n += b - a + 1;
                }
            }
        } else if part.parse::<i32>().is_ok() {
            n += 1;
        }
    }
    n
}

/// std::thread::hardware_concurrency() (glibc get_nprocs).
pub fn logical_cpus() -> i32 {
    if let Ok(s) = std::fs::read_to_string("/sys/devices/system/cpu/online") {
        let n = count_cpu_list(&s);
        if n > 0 {
            return n;
        }
    }
    if let Ok(s) = std::fs::read_to_string("/proc/stat") {
        let n = s.lines().filter(|l| l.starts_with("cpu") && l.as_bytes().get(3).is_some_and(|c| c.is_ascii_digit())).count() as i32;
        if n > 0 {
            return n;
        }
    }
    0
}

/// The size of this process's CPU affinity mask (libgomp's default team size).
fn affinity_cpus() -> i32 {
    if let Ok(s) = std::fs::read_to_string("/proc/self/status") {
        for line in s.lines() {
            if let Some(rest) = line.strip_prefix("Cpus_allowed_list:") {
                let n = count_cpu_list(rest);
                if n > 0 {
                    return n;
                }
            }
        }
    }
    let n = logical_cpus();
    if n > 0 {
        n
    } else {
        1
    }
}

/// The first value of OMP_NUM_THREADS, when it parses to a positive count.
fn env_omp_num_threads() -> Option<i32> {
    let v = std::env::var("OMP_NUM_THREADS").ok()?;
    let first = v.split(',').next()?.trim();
    let n: i64 = first.parse().ok()?;
    if n > 0 {
        Some(n.min(i32::MAX as i64) as i32)
    } else {
        None
    }
}

/// omp_set_num_threads(n): libgomp stores n > 0 ? n : 1 in nthreads-var.
pub fn omp_set_num_threads(n: i32) {
    G_OMP_NUM_THREADS_SET.store(if n > 0 { n } else { 1 }, Ordering::Relaxed);
}

/// The team size `#pragma omp parallel num_threads(n)` gets (OMP_THREAD_LIMIT caps it).
pub fn omp_team_size(n: i32) -> i32 {
    if let Ok(v) = std::env::var("OMP_THREAD_LIMIT") {
        if let Ok(l) = v.trim().parse::<i64>() {
            if l > 0 && (n as i64) > l {
                return l as i32;
            }
        }
    }
    n
}

/// ompMaxThreads(): omp_get_max_threads(), at least 1.
pub fn omp_max_threads() -> i32 {
    let set = G_OMP_NUM_THREADS_SET.load(Ordering::Relaxed);
    let n = if set > 0 { set } else { env_omp_num_threads().unwrap_or_else(affinity_cpus) };
    if n > 0 {
        n
    } else {
        1
    }
}

/// ompAvailable(): the reference build has OpenMP.
pub fn omp_available() -> i32 {
    1
}

/// perfCores(): Apple only; 0 elsewhere.
pub fn perf_cores() -> i32 {
    0
}

pub fn kokkos_omp_threads() -> i32 {
    G_KOKKOS_OMP_THREADS.load(Ordering::Relaxed)
}

pub fn set_kokkos_omp_threads(n: i32) {
    G_KOKKOS_OMP_THREADS.store(if n > 0 { n } else { 0 }, Ordering::Relaxed);
}

/// OPENSWMM_HOST_RESERVED_THREADS through atoi, negatives as 0.
pub fn host_reserved_threads() -> i32 {
    match std::env::var_os("OPENSWMM_HOST_RESERVED_THREADS") {
        None => 0,
        Some(s) => {
            let b = s.as_encoded_bytes();
            if b.is_empty() {
                return 0;
            }
            let n = crate::core::cstdlib::atoi(b);
            if n > 0 {
                n
            } else {
                0
            }
        }
    }
}

pub fn is_oversubscribed(threads: i32) -> bool {
    let logical = logical_cpus();
    logical > 0 && threads + host_reserved_threads() > logical
}

pub fn resolve_requested(requested: i32, what: &str, warnings: Option<&mut Vec<Bytes>>) -> i32 {
    let omp_max = omp_max_threads();
    if requested <= 0 {
        return omp_max;
    }

    let nt = requested;
    let warnings = match warnings {
        None => return nt,
        Some(w) => w,
    };

    let logical = logical_cpus();
    let reserved = host_reserved_threads();
    if logical > 0 && nt + reserved > logical {
        if reserved > 0 {
            warnings.push(
                format!(
                    "THREADS = {} plus the {} threads the host application reserves for itself exceed the {} logical processors on this machine ({}). The run is oversubscribed and is likely to be slower; active spin-wait is disabled for this run.",
                    nt, reserved, logical, what
                )
                .into_bytes(),
            );
        } else {
            warnings.push(
                format!(
                    "THREADS = {} exceeds the {} logical processors on this machine ({}). The run will be oversubscribed and is likely to be slower; active spin-wait is disabled for this run.",
                    nt, logical, what
                )
                .into_bytes(),
            );
        }
    } else if nt > omp_max {
        warnings.push(
            format!(
                "THREADS = {} exceeds the OpenMP limit of {} in this process (OMP_NUM_THREADS / OMP_THREAD_LIMIT / CPU affinity) for the {}; requesting {} anyway.",
                nt, omp_max, what, nt
            )
            .into_bytes(),
        );
    }
    // (__APPLE__ performance-core warning: not compiled on Linux.)
    nt
}

/// Minimum conduits of momentum work per dynamic-wave team thread (performance only).
pub const K_MIN_CONDUITS_PER_THREAD: i32 = 100;

pub fn dw_threads(requested: i32, n_conduits: i32, mut warnings: Option<&mut Vec<Bytes>>) -> i32 {
    let mut nt = resolve_requested(requested, "dynamic wave", warnings.as_deref_mut());

    if nt > 1 {
        let mut cap = n_conduits / K_MIN_CONDUITS_PER_THREAD;
        if cap < 1 {
            cap = 1;
        }
        if cap < nt {
            if requested > 0 {
                if let Some(w) = warnings.as_deref_mut() {
                    w.push(
                        format!(
                            "THREADS = {} reduced to {} for dynamic wave: the model has only {} conduits ({} per thread minimum).",
                            requested, cap, n_conduits, K_MIN_CONDUITS_PER_THREAD
                        )
                        .into_bytes(),
                    );
                }
            }
            nt = cap;
        }
    }

    if requested == 0 && nt > 1 {
        let pcores = perf_cores();
        if pcores > 2 && nt > pcores - 2 {
            nt = pcores - 2;
        }
    }
    nt
}

pub fn two_d_threads(requested: i32, n_triangles: i32, mut warnings: Option<&mut Vec<Bytes>>) -> i32 {
    let mut nt = resolve_requested(requested, "2D surface", warnings.as_deref_mut());
    if requested == 0 && nt > 1 {
        let pcores = perf_cores();
        if pcores > 2 && nt > pcores - 2 {
            nt = pcores - 2;
        }
    }
    if n_triangles < 4 * nt {
        if requested > 0 && nt > 1 {
            if let Some(w) = warnings.as_deref_mut() {
                w.push(
                    format!(
                        "THREADS = {} reduced to 1 for the 2D surface: the mesh has only {} triangles (need >= 4 per thread).",
                        requested, n_triangles
                    )
                    .into_bytes(),
                );
            }
        }
        nt = 1;
    }
    nt
}
