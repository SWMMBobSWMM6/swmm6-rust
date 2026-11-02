// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/TopoSort.hpp and TopoSort.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Topological link ordering (Kahn's algorithm over the node-to-outgoing-link CSR), used
// by kinematic-wave and steady routing.
//
// Differences in form, never in result:
//   * The raw node1/node2 pointers plus a count become slices plus the same `i32` count;
//     the nullable divert_link pointer is an Option of a slice.

/// Fill `sorted_links` with links in upstream-to-downstream order; returns how many were
/// ordered (fewer than `n_links` means a cycle).
pub fn sort_links(node1: &[i32], node2: &[i32], n_links: i32, n_nodes: i32, sorted_links: &mut Vec<i32>, divert_link: Option<&[i32]>) -> i32 {
    sorted_links.clear();
    sorted_links.reserve(n_links.max(0) as usize);

    let nn = n_nodes.max(0) as usize;
    let mut degree = vec![0i32; nn];
    let mut in_degree = vec![0i32; nn];

    // Count outgoing links per node
    for j in 0..n_links {
        let n1 = node1[j as usize];
        if n1 >= 0 && n1 < n_nodes {
            degree[n1 as usize] += 1;
        }
    }

    // Build CSR offsets
    let mut start_pos = vec![0i32; (n_nodes + 1).max(0) as usize];
    for i in 0..n_nodes {
        start_pos[(i + 1) as usize] = start_pos[i as usize] + degree[i as usize];
    }

    // Fill adjacency list
    let mut adj_list = vec![0i32; n_links.max(0) as usize];
    let mut cursor = start_pos.clone();
    for j in 0..n_links {
        let n1 = node1[j as usize];
        if n1 >= 0 && n1 < n_nodes {
            let c = cursor[n1 as usize];
            cursor[n1 as usize] += 1;
            adj_list[c as usize] = j;
        }
    }

    // legacy adjustAdjList: at a DIVIDER with exactly two outgoing links, put the
    // NON-diversion link first.
    if let Some(divert_link) = divert_link {
        for i in 0..n_nodes {
            let ui = i as usize;
            if divert_link[ui] < 0 {
                continue;
            }
            if degree[ui] != 2 {
                continue;
            }
            let m = start_pos[ui] as usize;
            if adj_list[m] == divert_link[ui] {
                adj_list.swap(m, m + 1);
            }
        }
    }

    // Compute in-degree
    for j in 0..n_links {
        let n2 = node2[j as usize];
        if n2 >= 0 && n2 < n_nodes {
            in_degree[n2 as usize] += 1;
        }
    }

    // Kahn's algorithm: stack of zero in-degree nodes
    let mut stack: Vec<i32> = Vec::with_capacity(nn);
    for i in 0..n_nodes {
        if in_degree[i as usize] == 0 {
            stack.push(i);
        }
    }

    let mut first = 0usize;
    while first < stack.len() {
        let i1 = stack[first];
        first += 1;
        let ui1 = i1 as usize;

        let k1 = start_pos[ui1];
        let k2 = start_pos[ui1 + 1];

        for k in k1..k2 {
            let j = adj_list[k as usize];
            sorted_links.push(j);

            let i2 = node2[j as usize];
            if i2 >= 0 && i2 < n_nodes {
                let ui2 = i2 as usize;
                in_degree[ui2] -= 1;
                if in_degree[ui2] == 0 {
                    stack.push(i2);
                }
            }
        }
    }

    sorted_links.len() as i32
}
