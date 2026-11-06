// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/controls/Controls.cpp and Controls.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The control-rule engine's run-time half (controls.rs holds the data and the parser):
// evaluate, batchEvaluateGroup, applyPendingActions, updateActionValue, evaluatePremise,
// compareValues, compareTimes, getVariableValue, computePIDSetting, getNamedVariableValue.
//
// NOT YET TRANSLATED: expression premises (mathexpr::evaluate over the named variables);
// evaluate returns StepError::NotTranslated when a rule holds one.

use crate::controls::controls::{Action, ActionType, CompareOp, ConditionVar, ControlEngine, LogicOp, PIDState, PendingAction, Premise, PremiseSoA};
use crate::core::cmath::cmax;
use crate::core::constants::{self, MISSING};
use crate::core::context::{ControlLogEntry, SimulationContext};
use crate::core::datetime;
use crate::core::engine_step_routing::StepError;
use crate::core::unit_conversion as ucf;
use crate::data::gage::GageData;
use crate::data::link::LinkType;
use crate::data::table::table_lookup_cursor;
use crate::hydraulics::link;
use crate::hydraulics::xsect_batch::XSectParams;

impl ControlEngine {
    /// batchEvaluateGroup: gather and compare one premise group.
    fn batch_evaluate_group(g: &mut PremiseSoA, this: &ControlEngine, premise_results: &mut [u8], ctx: &SimulationContext, current_time: f64, half_step: f64) {
        if g.count == 0 {
            return;
        }
        for i in 0..g.count as usize {
            g.lhs_values[i] = this.get_variable_value(ctx, g.var_type, g.obj_idx[i], current_time, g.param[i]);
        }
        let is_time = matches!(g.var_type, ConditionVar::SimTime | ConditionVar::ClockTime | ConditionVar::LinkTimeopen | ConditionVar::LinkTimeclosed);
        for i in 0..g.count as usize {
            let lhs = g.lhs_values[i];
            let rhs = g.rhs_value[i];
            let op = op_of(g.op[i]);
            let r: u8 = if lhs == MISSING {
                0
            } else if g.rhs_is_variable[i] {
                0
            } else if is_time {
                compare_times(lhs, op, rhs, half_step) as u8
            } else {
                compare_values(lhs, op, rhs) as u8
            };
            g.results[i] = r != 0;
            premise_results[g.flat_idx[i] as usize] = r;
        }
    }

    /// ControlEngine::evaluate: rules at the step's dated start; returns the action count.
    pub fn evaluate(&mut self, ctx: &mut SimulationContext, current_date: f64, dt: f64, rule_time_reached: bool) -> Result<i32, StepError> {
        if !rule_time_reached {
            self.last_action_count_ = 0;
            return Ok(0);
        }
        if self.rules_.iter().any(|r| r.premises.iter().any(|p| p.is_expression)) {
            return Err(StepError::NotTranslated("control-rule expression premises"));
        }

        self.eval_date_ = current_date;
        self.cur_date_ = current_date.floor();
        self.cur_time_ = current_date - self.cur_date_;
        self.elapsed_days_ = current_date - ctx.options.start_date;
        let current_time = current_date;

        let dt_days = dt / constants::SEC_PER_DAY;
        let half_step = dt_days / 2.0;
        self.pending_actions_.clear();

        // Phase 1a
        let mut groups = std::mem::take(&mut self.premise_groups_);
        let mut results = std::mem::take(&mut self.premise_results_);
        for g in groups.iter_mut() {
            Self::batch_evaluate_group(g, self, &mut results, ctx, current_time, half_step);
        }
        self.premise_groups_ = groups;
        self.premise_results_ = results;

        // Phase 1b: variable-RHS and expression premises
        for r in 0..self.rules_.len() {
            let base = self.rule_premise_offset_[r];
            for p in 0..self.rules_[r].premises.len() {
                let prem = self.rules_[r].premises[p];
                if !prem.is_expression && !prem.rhs_is_variable {
                    continue;
                }
                let ok = self.evaluate_premise(ctx, &prem, current_time, half_step);
                self.premise_results_[base as usize + p] = if ok { 1 } else { 0 };
            }
        }

        // Phase 2: per-rule AND/OR reduction
        self.rule_results_.clear();
        self.rule_results_.resize(self.rules_.len(), false);
        for r in 0..self.rules_.len() {
            let base = self.rule_premise_offset_[r];
            let mut result = false;
            for p in 0..self.rules_[r].premises.len() {
                let pval = self.premise_results_[base as usize + p];
                let logic = self.rules_[r].premises[p].logic;
                if p == 0 {
                    result = pval != 0;
                } else if logic == LogicOp::And {
                    if !result {
                        break;
                    }
                    result = pval != 0;
                } else if !result {
                    result = pval != 0;
                }
            }
            self.rule_results_[r] = result;
        }

        // Phase 3: actions of the fired rules (last-authored first)
        for ur in 0..self.rules_.len() {
            let r = ur as i32;
            let actions: Vec<Action> = if self.rule_results_[ur] { self.rules_[ur].then_actions.clone() } else { self.rules_[ur].else_actions.clone() };
            for a0 in actions.iter().rev() {
                let mut a = *a0;
                if a.type_ == ActionType::Pid || a.type_ == ActionType::Curve || a.type_ == ActionType::Timeseries {
                    if let Some(last) = self.rules_[ur].premises.last().copied() {
                        let _ = self.evaluate_premise(ctx, &last, current_time, half_step);
                    }
                }
                self.update_action_value(&mut a, ctx, current_time, dt);
                let pa = PendingAction { link_idx: a.link_idx, value: a.value, priority: self.rules_[ur].priority, rule_idx: r, type_: a.type_ };
                self.pending_actions_.push(pa);
            }
        }

        self.last_action_count_ = self.apply_pending_actions(ctx, current_time);
        Ok(self.last_action_count_)
    }

    /// applyPendingActions: one slot per link, the highest priority wins.
    fn apply_pending_actions(&mut self, ctx: &mut SimulationContext, _current_time: f64) -> i32 {
        for s in self.action_slots_.iter_mut() {
            *s = -1;
        }
        for k in 0..self.pending_actions_.len() as i32 {
            let pa = self.pending_actions_[k as usize];
            if pa.link_idx < 0 || pa.link_idx >= ctx.n_links() {
                continue;
            }
            let mut s = 0usize;
            let mut listed = false;
            while s < self.action_slots_.len() && self.action_slots_[s] >= 0 {
                let cur = self.action_slots_[s];
                if self.pending_actions_[cur as usize].link_idx == pa.link_idx {
                    if pa.priority > self.pending_actions_[cur as usize].priority {
                        self.action_slots_[s] = k;
                    }
                    listed = true;
                    break;
                }
                s += 1;
            }
            if listed {
                continue;
            }
            if s < self.action_slots_.len() {
                self.action_slots_[s] = k;
            } else {
                self.action_slots_.insert(0, k);
            }
        }

        let mut changes = 0;
        for si in 0..self.action_slots_.len() {
            let slot = self.action_slots_[si];
            if slot < 0 {
                break;
            }
            let pa = self.pending_actions_[slot as usize];
            let ul = pa.link_idx as usize;
            if ctx.links.target_setting[ul] != pa.value {
                ctx.links.target_setting[ul] = pa.value;
                if ctx.options.rpt_controls && pa.type_ == ActionType::Numeric {
                    let ri = pa.rule_idx;
                    let known = ri >= 0 && (ri as usize) < self.rules_.len();
                    if known {
                        let names = &mut ctx.control_rule_names;
                        if names.len() < self.rules_.len() {
                            names.resize(self.rules_.len(), Vec::new());
                        }
                        let slot_name = &mut names[ri as usize];
                        if *slot_name != self.rules_[ri as usize].name {
                            *slot_name = self.rules_[ri as usize].name.clone();
                        }
                    }
                    let entry = ControlLogEntry { link_idx: pa.link_idx, rule_idx: if known { ri } else { -1 }, new_setting: pa.value, date: self.eval_date_ };
                    ctx.log_control_action(entry);
                }
                changes += 1;
            }
        }
        changes
    }

    /// updateActionValue: modulated actions (curve, time series, PID).
    fn update_action_value(&mut self, a: &mut Action, ctx: &mut SimulationContext, _current_time: f64, dt: f64) {
        match a.type_ {
            ActionType::Numeric => {}
            ActionType::Curve => {
                if a.curve_idx >= 0 && (a.curve_idx as usize) < ctx.tables.count() {
                    a.value = table_lookup_cursor(&mut ctx.tables.tables[a.curve_idx as usize], self.control_value_);
                }
            }
            ActionType::Timeseries => {
                if a.tseries_idx >= 0 && (a.tseries_idx as usize) < ctx.tables.count() {
                    a.value = table_lookup_cursor(&mut ctx.tables.tables[a.tseries_idx as usize], self.eval_date_);
                }
            }
            ActionType::Pid => {
                if a.pid_idx >= 0 && (a.pid_idx as usize) < self.pid_states_.len() {
                    let set_point = self.set_point_;
                    let control_value = self.control_value_;
                    let pid = &mut self.pid_states_[a.pid_idx as usize];
                    pid.setpoint = set_point;
                    let in_range = a.link_idx >= 0 && a.link_idx < ctx.n_links();
                    let cur_setting = if in_range { ctx.links.target_setting[a.link_idx as usize] } else { 0.0 };
                    let is_pump = if in_range { ctx.links.type_[a.link_idx as usize] == LinkType::Pump } else { false };
                    a.value = compute_pid_setting(pid, control_value, cur_setting, is_pump, dt);
                }
            }
        }
    }

    /// evaluatePremise: the scalar path (sets control_value_ / set_point_).
    fn evaluate_premise(&mut self, ctx: &SimulationContext, p: &Premise, current_time: f64, half_step: f64) -> bool {
        // (expression premises: refused in evaluate)
        let lhs = self.get_variable_value(ctx, p.lhs_var, p.lhs_idx, current_time, p.lhs_param);
        if lhs == MISSING {
            return false;
        }
        let rhs = if p.rhs_is_variable { self.get_variable_value(ctx, p.rhs_var, p.rhs_idx, current_time, p.rhs_param) } else { p.rhs_value };
        if rhs == MISSING {
            return false;
        }
        self.control_value_ = lhs;
        self.set_point_ = rhs;

        let is_time = matches!(p.lhs_var, ConditionVar::SimTime | ConditionVar::ClockTime | ConditionVar::LinkTimeopen | ConditionVar::LinkTimeclosed);
        if is_time {
            let result = compare_times(lhs, p.op, rhs, half_step);
            if p.lhs_var == ConditionVar::LinkTimeopen || p.lhs_var == ConditionVar::LinkTimeclosed {
                self.control_value_ = lhs * 24.0;
            }
            return result;
        }
        compare_values(lhs, p.op, rhs)
    }

    /// getVariableValue: an attribute in display units (MISSING when unavailable).
    fn get_variable_value(&self, ctx: &SimulationContext, var: ConditionVar, idx: i32, _current_time: f64, param: i32) -> f64 {
        let ui = idx as usize;
        let fu = ctx.options.flow_units as i32;
        let us = ucf::get_unit_system(fu) as usize;
        let ucf_len = ucf::UCF[3][us];
        let ucf_flow = ucf::QCF[fu as usize];
        let ucf_vol = ucf::UCF[5][us];
        let node_ok = idx >= 0 && idx < ctx.n_nodes();
        let link_ok = idx >= 0 && idx < ctx.n_links();

        match var {
            ConditionVar::NodeDepth => {
                if node_ok {
                    ctx.nodes.depth[ui] * ucf_len
                } else {
                    MISSING
                }
            }
            ConditionVar::NodeMaxdepth => {
                if node_ok {
                    ctx.nodes.full_depth[ui] * ucf_len
                } else {
                    MISSING
                }
            }
            ConditionVar::NodeHead => {
                if node_ok {
                    ctx.nodes.head[ui] * ucf_len
                } else {
                    MISSING
                }
            }
            ConditionVar::NodeVolume => {
                if node_ok {
                    ctx.nodes.volume[ui] * ucf_vol
                } else {
                    MISSING
                }
            }
            ConditionVar::NodeInflow => {
                if node_ok {
                    ctx.nodes.lat_flow[ui] * ucf_flow
                } else {
                    MISSING
                }
            }
            ConditionVar::LinkFlow => {
                if link_ok {
                    ctx.links.flow[ui] * (ctx.links.direction[ui] as f64) * ucf_flow
                } else {
                    MISSING
                }
            }
            ConditionVar::LinkDepth => {
                if link_ok {
                    ctx.links.depth[ui] * ucf_len
                } else {
                    MISSING
                }
            }
            ConditionVar::LinkSetting => {
                if link_ok {
                    ctx.links.setting[ui]
                } else {
                    MISSING
                }
            }
            ConditionVar::LinkStatus => {
                if !link_ok {
                    return MISSING;
                }
                if ctx.links.type_[ui] != LinkType::Conduit && ctx.links.type_[ui] != LinkType::Pump {
                    return MISSING;
                }
                ctx.links.setting[ui]
            }
            ConditionVar::LinkFullflow => {
                if !link_ok {
                    return MISSING;
                }
                let cr = ctx.link_subtypes.conduit_row(idx);
                (if cr >= 0 { ctx.link_subtypes.conduits.q_full[cr as usize] } else { 0.0 }) * ucf_flow
            }
            ConditionVar::LinkFulldepth => {
                if link_ok {
                    ctx.links.xsect_y_full[ui] * ucf_len
                } else {
                    MISSING
                }
            }
            ConditionVar::LinkLength => {
                if !link_ok {
                    return MISSING;
                }
                let cr = ctx.link_subtypes.conduit_row(idx);
                (if cr >= 0 { ctx.link_subtypes.conduits.length[cr as usize] } else { 0.0 }) * ucf_len
            }
            ConditionVar::LinkSlope => {
                if !link_ok {
                    return MISSING;
                }
                let cr = ctx.link_subtypes.conduit_row(idx);
                if cr >= 0 {
                    ctx.link_subtypes.conduits.slope[cr as usize]
                } else {
                    0.0
                }
            }
            ConditionVar::LinkVelocity => {
                if !link_ok {
                    return MISSING;
                }
                if ctx.links.type_[ui] != LinkType::Conduit {
                    return MISSING;
                }
                let mut xs = XSectParams::default();
                xs.type_ = ctx.links.xsect_batch_shape[ui];
                xs.y_full = ctx.links.xsect_y_full[ui];
                xs.a_full = ctx.links.xsect_a_full[ui];
                xs.w_max = ctx.links.xsect_w_max[ui];
                let cr = ctx.link_subtypes.conduit_row(idx);
                let barrels = if cr >= 0 { ctx.link_subtypes.conduits.barrels[cr as usize] } else { 1 };
                link::get_velocity(&xs, ctx.links.flow[ui], ctx.links.depth[ui], barrels) * ucf_len
            }
            ConditionVar::LinkTimeopen => {
                if !link_ok {
                    return MISSING;
                }
                if ctx.links.setting[ui] <= 0.0 {
                    return MISSING;
                }
                if ui < ctx.links.time_last_set.len() {
                    return self.cur_date_ + self.cur_time_ - ctx.links.time_last_set[ui];
                }
                0.0
            }
            ConditionVar::LinkTimeclosed => {
                if !link_ok {
                    return MISSING;
                }
                if ctx.links.setting[ui] > 0.0 {
                    return MISSING;
                }
                if ui < ctx.links.time_last_set.len() {
                    return self.cur_date_ + self.cur_time_ - ctx.links.time_last_set[ui];
                }
                0.0
            }
            ConditionVar::GageRain => {
                if idx < 0 || idx >= ctx.n_gages() {
                    return MISSING;
                }
                if ctx.gages.st_used[ui] != 0 {
                    ctx.gages.rainfall[ui]
                } else {
                    0.0
                }
            }
            ConditionVar::GageRainPast => {
                if idx < 0 || idx >= ctx.n_gages() {
                    return MISSING;
                }
                if param < 1 || param > GageData::MAXPASTRAIN {
                    return MISSING;
                }
                if ctx.gages.st_used[ui] == 0 {
                    return 0.0;
                }
                let mut total = 0.0;
                let base = ui * GageData::MAXPASTRAIN as usize;
                for h in 0..param as usize {
                    total += ctx.gages.past_rain[base + h];
                }
                total
            }
            ConditionVar::SimTime => self.elapsed_days_,
            ConditionVar::SimDate => self.cur_date_,
            ConditionVar::ClockTime => self.cur_time_,
            ConditionVar::SimDay => {
                let t = (self.cur_date_ as i32).wrapping_add(datetime::DATE_DELTA);
                ((t % 7) + 1) as f64
            }
            ConditionVar::SimMonth => datetime::month_of_year(self.cur_date_) as f64,
            ConditionVar::SimDayofyear => datetime::day_of_year(self.cur_date_) as f64,
        }
    }
}

fn op_of(v: i32) -> CompareOp {
    match v {
        1 => CompareOp::Ne,
        2 => CompareOp::Lt,
        3 => CompareOp::Le,
        4 => CompareOp::Gt,
        5 => CompareOp::Ge,
        _ => CompareOp::Eq,
    }
}

/// compareValues
fn compare_values(lhs: f64, op: CompareOp, rhs: f64) -> bool {
    match op {
        CompareOp::Eq => lhs == rhs,
        CompareOp::Ne => lhs != rhs,
        CompareOp::Lt => lhs < rhs,
        CompareOp::Le => lhs <= rhs,
        CompareOp::Gt => lhs > rhs,
        CompareOp::Ge => lhs >= rhs,
    }
}

/// compareTimes: = and <> within half a step.
fn compare_times(lhs: f64, op: CompareOp, rhs: f64, half_step: f64) -> bool {
    if op == CompareOp::Eq {
        return lhs >= rhs - half_step && lhs < rhs + half_step;
    }
    if op == CompareOp::Ne {
        return lhs < rhs - half_step || lhs >= rhs + half_step;
    }
    compare_values(lhs, op, rhs)
}

/// computePIDSetting (legacy getPIDSetting).
fn compute_pid_setting(pid: &mut PIDState, control_value: f64, current_setting: f64, is_pump: bool, dt: f64) -> f64 {
    if dt <= 0.0 {
        return current_setting;
    }
    let dt_min = (dt / constants::SEC_PER_DAY) * 1440.0;
    const TINY: f64 = 1.0e-6;
    const TOLERANCE: f64 = 0.0001;

    let mut e0 = pid.setpoint - control_value;
    if e0.abs() > TINY {
        if pid.setpoint != 0.0 {
            e0 /= pid.setpoint;
        } else {
            e0 /= control_value;
        }
    }

    if (e0 - pid.e1).abs() < TOLERANCE {
        pid.e2 = 0.0;
        pid.e1 = 0.0;
    }

    let p = e0 - pid.e1;
    let i = if pid.ki != 0.0 { e0 * dt_min / pid.ki } else { 0.0 };
    let d = if dt_min > 0.0 { pid.kd * (e0 - 2.0 * pid.e1 + pid.e2) / dt_min } else { 0.0 };
    let mut update = pid.kp * (p + i + d);
    if update.abs() < TOLERANCE {
        update = 0.0;
    }

    let mut setting = current_setting + update;

    pid.e2 = pid.e1;
    pid.e1 = e0;

    if setting < 0.0 {
        setting = 0.0;
    }
    if !is_pump && setting > 1.0 {
        setting = 1.0;
    }
    let _ = cmax(0.0, 0.0);
    setting
}
