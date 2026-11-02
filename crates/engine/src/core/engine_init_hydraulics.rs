// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/SWMMEngine.cpp (initHydraulics) and
// SWMMEngine.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// SWMMEngine::initHydraulics(): the router (XSectGroups, DW/KW/STEADY/FV setup, outfall
// map), the layout validations (ERROR 141/145/115/133/135/134), the thread
// configuration, the routing-step histogram, the inflow, RDII, structure, exfiltration,
// inlet and interface-file solvers, event sorting and the [CONTROLS] rule parser.
//
// Differences in form, never in result:
//   * OPENSWMM_HAS_2D is off in the reference build, so the 2D blocks (has_2d_domain,
//     surface_router_.initialize and its thread warnings) are absent here too.
//   * omp_set_num_threads and the `#pragma omp parallel` team-size probe go through
//     core::thread_info's model of libgomp (see that module).
//   * router_.setStructureSolver(&hydstruct_) records the attachment
//     (Router::structures_attached_); the structure solver stays an engine member.
//   * initHydraulics is `noexcept` upstream. The one exception known to escape it is
//     std::stod's, through ControlEngine::parseRuleText -> addExpression ->
//     mathexpr::parse (finding U1, `EXPRESSION E1 = 2 + .`), which ends the process with
//     std::terminate. Here init_hydraulics returns Err(Terminated) at that point instead,
//     with everything before it done exactly as upstream did it.
//   * The routing-event sort is libstdc++'s std::sort (multi_column_series_file::std_sort),
//     so events with equal starts keep upstream's order.
//   * The RDII UH validation iterates the UH name map in libstdc++'s unordered_map order:
//     RDIISolver keeps the names in insertion order, and
//     std_unordered_map::ci_map_iteration_order replays the inserts with CiHash codes, so
//     the error order is upstream's.

use crate::core::cmath::{cmax, cmin};
use crate::core::constants;
use crate::core::std_unordered_map::ci_map_iteration_order;
use crate::core::engine::{SwmmEngine, SWMM_ERR_PARSE};
use crate::core::error_codes::{
    format_error, format_error_detail, ERR_CONTROL_RULE, ERR_DIVIDER, ERR_DUMMY_LINK, ERR_MULTI_OUTLET, ERR_NO_OUTLETS, ERR_OUTFALL, ERR_RDII_AREA, ERR_SLOPE, ERR_UNITHYD_RATIOS, ERR_UNITHYD_TIMES,
};
use crate::core::options::{Bytes, RoutingModel};
use crate::core::thread_info;
use crate::data::link::{LinkType, XsectShape};
use crate::data::node::NodeType;
use crate::hydraulics::routing::RouteModel;

/// Upstream would have ended the process here (std::terminate from a noexcept function).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Terminated {
    /// what() of the escaping exception, as libstdc++'s terminate handler prints it
    pub what: &'static str,
    /// the exception type
    pub exception: &'static str,
}

impl SwmmEngine {
    /// The gate's access to the private initHydraulics (upstream's harness reaches it
    /// through `#define private public`).
    #[doc(hidden)]
    pub fn gate_init_hydraulics(&mut self) -> Result<(), Terminated> {
        self.init_hydraulics()
    }

    /// SWMMEngine::initHydraulics.
    pub(crate) fn init_hydraulics(&mut self) -> Result<(), Terminated> {
        // 1. Router
        let mut rm = RouteModel::Dynwave;
        if self.ctx_.options.routing_model == RoutingModel::Kinwave {
            rm = RouteModel::Kinwave;
        } else if self.ctx_.options.routing_model == RoutingModel::Steady {
            rm = RouteModel::Steady;
        } else if self.ctx_.options.routing_model == RoutingModel::Fv {
            rm = RouteModel::Fv;
        }
        self.router_.init(&mut self.ctx_, rm);
        self.ctx_.routing_stats.lstep_mode = 0;
        if self.ctx_.options.rpt_link_steps && self.ctx_.n_links() > 0 && self.router_.enable_link_step_stats(&self.ctx_) {
            self.ctx_.routing_stats.lstep_mode = if rm == RouteModel::Fv { 2 } else { 1 };
        }

        // Surface what the FV mesh builder found.
        for w in self.router_.fv_warnings().to_vec() {
            self.ctx_.warnings.push(w);
        }
        for e in self.router_.fv_errors().to_vec() {
            self.ctx_.errors.push(e);
        }

        let n_nodes = self.ctx_.nodes.count();
        self.ctx_.node_subtypes.rebuild_index(n_nodes);

        // Gap #83 / #84: outfall connectivity, no outlets, adverse slope, tree layout.
        {
            let n_nl = self.ctx_.n_nodes();
            let n_ll = self.ctx_.n_links();
            let mut n_in = vec![0i32; n_nl as usize];
            let mut n_out = vec![0i32; n_nl as usize];
            for j in 0..n_ll {
                let uj = j as usize;
                let n1 = self.ctx_.links.node1[uj];
                let n2 = self.ctx_.links.node2[uj];
                if n1 >= 0 && n1 < n_nl {
                    n_out[n1 as usize] += 1;
                }
                if n2 >= 0 && n2 < n_nl {
                    n_in[n2 as usize] += 1;
                }
            }
            let mut n_outlets = 0;
            for i in 0..n_nl {
                let ui = i as usize;
                if self.ctx_.nodes.type_[ui] != NodeType::Outfall {
                    continue;
                }
                if rm == RouteModel::Dynwave && n_in[ui] + n_out[ui] > 1 {
                    let e = format_error(ERR_OUTFALL, &self.ctx_.node_names.names()[ui]);
                    self.push_error_and_set(e);
                } else {
                    n_outlets += 1;
                }
            }
            let has_2d_domain = false; // OPENSWMM_HAS_2D off
            if n_outlets == 0 && (rm == RouteModel::Dynwave || rm == RouteModel::Fv) && !has_2d_domain {
                let e = format_error(ERR_NO_OUTLETS, b"");
                self.push_error_and_set(e);
            }
            if rm != RouteModel::Dynwave && rm != RouteModel::Fv {
                for j in 0..n_ll {
                    let uj = j as usize;
                    if self.ctx_.links.type_[uj] != LinkType::Conduit {
                        continue;
                    }
                    if self.ctx_.links.xsect_shape[uj] == XsectShape::Dummy {
                        continue;
                    }
                    let cr = self.ctx_.link_subtypes.conduit_row(j);
                    if cr >= 0 && self.ctx_.link_subtypes.conduits.slope[cr as usize] < 0.0 {
                        let e = format_error(ERR_SLOPE, &self.ctx_.link_names.names()[uj]);
                        self.push_error_and_set(e);
                    }
                }
                for i in 0..n_nl {
                    let ui = i as usize;
                    let deg = n_out[ui];
                    match self.ctx_.nodes.type_[ui] {
                        NodeType::Divider => {
                            if deg > 2 {
                                let e = format_error(ERR_DIVIDER, &self.ctx_.node_names.names()[ui]);
                                self.push_error_and_set(e);
                            }
                        }
                        NodeType::Outfall => {
                            if deg > 0 {
                                let e = format_error(ERR_OUTFALL, &self.ctx_.node_names.names()[ui]);
                                self.push_error_and_set(e);
                            }
                        }
                        NodeType::Storage => {}
                        _ => {
                            if deg > 1 {
                                let e = format_error(ERR_MULTI_OUTLET, &self.ctx_.node_names.names()[ui]);
                                self.push_error_and_set(e);
                            }
                        }
                    }
                }
            }
        }

        // Gap #44: routing loop
        if self.router_.has_cycle() {
            self.set_error(SWMM_ERR_PARSE, b"Routing loop detected: network contains a cycle that cannot be topologically sorted (KW/STEADY routing only).");
        }

        // 1b. OpenMP thread count from THREADS.
        {
            let mut thread_warnings: Vec<Bytes> = Vec::new();
            let nt = thread_info::resolve_requested(self.ctx_.options.num_threads, "OpenMP team", Some(&mut thread_warnings));
            thread_info::omp_set_num_threads(nt);
            if nt > 1 {
                let got = thread_info::omp_team_size(nt);
                if got < nt {
                    thread_warnings.push(format!("THREADS = {} requested but the OpenMP runtime provided a team of {} (OMP_THREAD_LIMIT or runtime policy).", nt, got).into_bytes());
                }
            }
            if rm == RouteModel::Dynwave {
                let n = self.ctx_.options.num_threads;
                self.router_.set_dw_num_threads(n, Some(&mut thread_warnings));
            }
            for w in &thread_warnings {
                self.push_report_warning(w, 0);
            }
        }

        // Routing time-step histogram bins.
        {
            let mut hist_min_step = self.ctx_.options.min_routing_step;
            if rm == RouteModel::Dynwave || rm == RouteModel::Fv {
                hist_min_step = cmax(cmin(hist_min_step, self.ctx_.options.routing_step), constants::MIN_TIMESTEP);
            }
            let rs = self.ctx_.options.routing_step;
            self.ctx_.routing_stats.init_histogram(rs, hist_min_step);
        }

        // 9. Inflow solver
        self.inflow_.init(&mut self.ctx_);

        // 10a. RDII solver + Gap #80 UH validation
        self.rdii_.init(&mut self.ctx_);
        {
            // `for (const auto& [uh_name, uh_i] : rdii_.uhNameIndex())`: the map's order.
            let order: Vec<(Bytes, i32)> = ci_map_iteration_order(self.rdii_.uh_name_index());
            for (uh_name, uh_i) in &order {
                let uh_i = *uh_i;
                if uh_i < 0 || uh_i >= self.rdii_.uh_params.len() as i32 {
                    continue;
                }
                let uh = self.rdii_.uh_params[uh_i as usize];
                let mut time_err = false;
                let mut ratio_err = false;
                let mut m = 0;
                while m < 12 && (!time_err || !ratio_err) {
                    let mut rsum = 0.0;
                    for k in 0..3 {
                        if uh.tBase[m][k] == 0.0 {
                            continue;
                        }
                        if !time_err && uh.tPeak[m][k] < 0.0 {
                            let e = format_error(ERR_UNITHYD_TIMES, uh_name);
                            self.push_error_and_set(e);
                            time_err = true;
                        }
                        if uh.r[m][k] < 0.0 {
                            if !ratio_err {
                                let e = format_error(ERR_UNITHYD_RATIOS, uh_name);
                                self.push_error_and_set(e);
                                ratio_err = true;
                            }
                        } else {
                            rsum += uh.r[m][k];
                        }
                    }
                    if !ratio_err && rsum > 1.01 {
                        let e = format_error(ERR_UNITHYD_RATIOS, uh_name);
                        self.push_error_and_set(e);
                        ratio_err = true;
                    }
                    m += 1;
                }
            }
            // RDII sewer areas
            for i in 0..self.ctx_.rdii_assigns.count() {
                let ui = i as usize;
                if self.ctx_.rdii_assigns.sewer_area[ui] < 0.0 {
                    let ni = self.ctx_.rdii_assigns.node_idx[ui];
                    let nname: Bytes = if ni >= 0 && ni < self.ctx_.n_nodes() { self.ctx_.node_names.names()[ni as usize].clone() } else { b"?".to_vec() };
                    let e = format_error(ERR_RDII_AREA, &nname);
                    self.push_error_and_set(e);
                }
            }
        }

        // 10b. Non-conduit hydraulic structures
        self.hydstruct_.init(&mut self.ctx_);

        // 10b-1. Illegal DUMMY / ideal-pump connections (legacy ERROR 134).
        if rm == RouteModel::Dynwave || rm == RouteModel::Fv {
            let n_ll = self.ctx_.n_links();
            let n_nl = self.ctx_.n_nodes();
            let ctx = &self.ctx_;
            let is_pass_through = |uj: usize| -> bool {
                if ctx.links.type_[uj] == LinkType::Conduit {
                    return ctx.links.xsect_shape[uj] == XsectShape::Dummy;
                }
                if ctx.links.type_[uj] != LinkType::Pump {
                    return false;
                }
                let pr = ctx.link_subtypes.pump_row(uj as i32);
                pr >= 0 && ctx.link_subtypes.pumps.curve_type[pr as usize] == 6
            };

            let mut n_out = vec![0i32; n_nl as usize];
            let mut all_in_pass = vec![0i32; n_nl as usize];
            for j in 0..n_ll {
                let uj = j as usize;
                let n1 = ctx.links.node1[uj];
                let n2 = ctx.links.node2[uj];
                if n1 < 0 || n1 >= n_nl || n2 < 0 || n2 >= n_nl {
                    continue;
                }
                let rev = ctx.links.direction[uj] < 0;

                let mut up = if rev { n2 } else { n1 };
                if ctx.nodes.type_[up as usize] == NodeType::Outfall {
                    up = if rev { n1 } else { n2 };
                }
                n_out[up as usize] += 1;

                let m = &mut all_in_pass[(if rev { n1 } else { n2 }) as usize];
                if !is_pass_through(uj) {
                    *m = -1;
                } else if *m == 0 {
                    *m = 1;
                }
            }

            let mut reported = vec![0u8; n_nl as usize];
            let mut to_report: Vec<i32> = Vec::new();
            let mut report = |node: i32, to_report: &mut Vec<i32>| {
                if node < 0 || node >= n_nl {
                    return;
                }
                let un = node as usize;
                if reported[un] != 0 {
                    return;
                }
                reported[un] = 1;
                to_report.push(node);
            };

            for j in 0..n_ll {
                let uj = j as usize;
                if !is_pass_through(uj) {
                    continue;
                }
                let n1 = ctx.links.node1[uj];
                if n1 < 0 || n1 >= n_nl {
                    continue;
                }
                let u1 = n1 as usize;
                if n_out[u1] > 1 {
                    report(n1, &mut to_report);
                }
                if all_in_pass[u1] > 0 {
                    report(n1, &mut to_report);
                }
                if ctx.links.type_[uj] == LinkType::Conduit && ctx.nodes.type_[u1] == NodeType::Storage {
                    report(n1, &mut to_report);
                }
            }
            // Each report pushes the error and calls set_error, in report order.
            for node in to_report {
                let e = format_error(ERR_DUMMY_LINK, &self.ctx_.node_names.names()[node as usize]);
                self.push_error_and_set(e);
            }
        }

        // KINWAVE/STEADY evaluate structures inline: wire the structure solver.
        self.router_.set_structure_solver_attached(true);

        // 10c. Exfiltration solver
        self.exfil_.init(&mut self.ctx_);

        // 10d. Inlet solver
        self.inlet_.init(&mut self.ctx_);

        // 10d. Interface file manager (files opened later in start())
        self.iface_.init(&mut self.ctx_);

        // 10d-2. Sort routing events chronologically and resolve overlaps.
        if !self.ctx_.events.is_empty() {
            crate::input::multi_column_series_file::std_sort(&mut self.ctx_.events, &|a: &crate::core::context::Event, b: &crate::core::context::Event| a.start < b.start);
            let mut i = 0usize;
            while i + 1 < self.ctx_.events.len() {
                if self.ctx_.events[i].end > self.ctx_.events[i + 1].start {
                    self.ctx_.events[i].end = self.ctx_.events[i + 1].start;
                }
                i += 1;
            }
            self.next_event_ = 0;
            self.between_events_ = true;
        }

        // 10e. Control rule parsing from [CONTROLS] section text.
        if self.ctx_.control_rules.count() > 0 {
            self.controls_.clear_rules();
            for i in 0..self.ctx_.control_rules.rule_text.len() {
                let text = self.ctx_.control_rules.rule_text[i].clone();
                let rc = match self.controls_.parse_rule_text(&text, &self.ctx_) {
                    Ok(rc) => rc,
                    Err(e) => {
                        let (exception, what) = match e {
                            crate::math::math_expr::StodError::InvalidArgument => ("std::invalid_argument", "stod"),
                            crate::math::math_expr::StodError::OutOfRange => ("std::out_of_range", "stod"),
                        };
                        return Err(Terminated { what, exception });
                    }
                };
                if rc < 0 {
                    let pe = self.controls_.last_parse_error().clone();
                    self.ctx_.error_code = 217; // legacy ERR_RULE
                    let mut where_: Bytes = format!("in [CONTROLS] rule block #{}", i + 1).into_bytes();
                    if pe.line > 0 {
                        where_.extend_from_slice(format!(", line {}", pe.line).as_bytes());
                    }
                    if !pe.message.is_empty() {
                        where_.extend_from_slice(b": ");
                        where_.extend_from_slice(&pe.message);
                    }
                    self.ctx_.error_message = format_error_detail(if pe.code > 0 { pe.code } else { ERR_CONTROL_RULE }, &pe.token, &where_);
                    let m = self.ctx_.error_message.clone();
                    self.ctx_.errors.push(m);
                    return Ok(());
                }
            }
        }
        Ok(())
    }

    /// `ctx_.errors.push_back(msg); set_error(SWMM_ERR_PARSE, ctx_.errors.back().c_str());`
    fn push_error_and_set(&mut self, msg: Bytes) {
        self.ctx_.errors.push(msg.clone());
        self.set_error(SWMM_ERR_PARSE, &msg);
    }
}
