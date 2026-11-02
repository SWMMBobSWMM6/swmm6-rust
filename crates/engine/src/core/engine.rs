// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/SWMMEngine.hpp and SWMMEngine.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The engine lifecycle manager, so far: construction, open() (1D build: OPENSWMM_HAS_2D
// off) and the private helpers open() calls (validate_project, set_error, emit_warning,
// push_report_warning, write_open_failure_report). initialize/start/step/end/report/close
// and the solver members they drive come with the solver phases.
//
// Differences of form:
//   * Only the members open() touches exist yet. The solver members (router, runoff,
//     quality, ...) are added as their solvers are translated.
//   * ctx_.state_accessors: upstream installs closures over the engine's solvers
//     (FV dry state, channel thermal/species regions, snow, infiltration and groundwater
//     state). They read and write solver state no open() result depends on, and the
//     solvers are not translated yet, so they are left unset here. NOT YET TRANSLATED.
//   * write_open_failure_report writes the .rpt of a failed open through
//     DefaultReportPlugin::prepare/finalize, which take the context const. The .rpt writer
//     is not translated yet, so no report file is written. NOT YET TRANSLATED.
//   * perf:: timers (reset_load, reset_fv, ScopedTimer, begin_process_model) only time
//     the run; they are left out.
//   * emit_warning's callback is a boxed closure; the C callback's engine pointer and user
//     data are the closure's captures.

use crate::core::cmath::{cmin, cvt_i32};
use crate::core::context::{EngineState, SimulationContext};
use crate::core::error_codes::{format_error, format_warning, ERR_RAIN_GAGE_INTERVAL, ERR_SUBCATCH_OUTLET, WARN_DRY_STEP_INCREASED, WARN_GAGE_INTERVAL, WARN_ROUTING_STEP_REDUCED, WARN_WET_STEP_REDUCED};
use crate::core::fs_path::FsPath;
use crate::core::options::{Bytes, QualitySolverKind};
use crate::core::plugins::{DefaultInputPlugin, DefaultOutputPlugin, DefaultReportPlugin, DefaultStateIOPlugin, PluginFactory};
use crate::data::gage::RainSource;
use crate::data::species_registry::SpeciesKind;

// Error codes returned across the C ABI (the public SWMM_ErrorCode values).
pub const SWMM_OK: i32 = 0;
pub const SWMM_ERR_MEMORY: i32 = 1;
pub const SWMM_ERR_FILE_NOT_FOUND: i32 = 2;
pub const SWMM_ERR_PARSE: i32 = 5;
pub const SWMM_ERR_WRONG_STATE: i32 = 6;
pub const SWMM_ERR_PLUGIN: i32 = 10;

/// The warning callback of EngineCallbacks (code, message).
pub type WarningCallback = Box<dyn FnMut(i32, &[u8])>;

#[derive(Default)]
pub struct EngineCallbacks {
    pub on_warning: Option<WarningCallback>,
}

pub struct SwmmEngine {
    ctx_: SimulationContext,
    /// permissive open (see set_lenient_open)
    lenient_open_: bool,
    /// parse-only handle cannot initialize a run
    metadata_only_: bool,
    plugins_: PluginFactory,
    callbacks_: EngineCallbacks,
    /// Report file path
    rpt_path_: Bytes,
    /// Binary output file path
    out_path_: Bytes,
    /// index of the built-in .out writer in plugins_.output_plugins() (owned there)
    binary_output_plugin_: Option<usize>,
}

impl Default for SwmmEngine {
    fn default() -> Self {
        Self::new()
    }
}

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut m = Vec::new();
    for p in parts {
        m.extend_from_slice(p);
    }
    m
}

impl SwmmEngine {
    /// The constructor. Upstream registers the four process components once per process
    /// here; ProcessComponentRegistry::instance() does that on first use.
    pub fn new() -> Self {
        let mut e = SwmmEngine {
            ctx_: SimulationContext::default(),
            lenient_open_: false,
            metadata_only_: false,
            plugins_: PluginFactory::new(),
            callbacks_: EngineCallbacks::default(),
            rpt_path_: Vec::new(),
            out_path_: Vec::new(),
            binary_output_plugin_: None,
        };
        e.ctx_.state = EngineState::Created;
        e
    }

    pub fn context(&self) -> &SimulationContext {
        &self.ctx_
    }
    pub fn context_mut(&mut self) -> &mut SimulationContext {
        &mut self.ctx_
    }
    pub fn plugin_factory(&self) -> &PluginFactory {
        &self.plugins_
    }
    pub fn set_lenient_open(&mut self, on: bool) {
        self.lenient_open_ = on;
    }
    pub fn lenient_open(&self) -> bool {
        self.lenient_open_
    }
    pub fn metadata_only(&self) -> bool {
        self.metadata_only_
    }
    pub fn rpt_path(&self) -> &[u8] {
        &self.rpt_path_
    }
    pub fn out_path(&self) -> &[u8] {
        &self.out_path_
    }
    pub fn binary_output_plugin(&self) -> Option<usize> {
        self.binary_output_plugin_
    }
    pub fn last_error_message(&self) -> &[u8] {
        &self.ctx_.error_message
    }
    pub fn set_warning_callback(&mut self, cb: Option<WarningCallback>) {
        self.callbacks_.on_warning = cb;
    }

    // ========================================================================
    // open()
    // ========================================================================

    /// `swmm_engine_open` / `swmm_engine_open_metadata`. A null C string is `None`.
    pub fn open(&mut self, inp_path: Option<&[u8]>, rpt_path: Option<&[u8]>, out_path: Option<&[u8]>, input_plugin_lib: Option<&[u8]>, metadata_only: bool) -> i32 {
        if self.ctx_.state != EngineState::Created && self.ctx_.state != EngineState::Closed {
            self.set_error(SWMM_ERR_WRONG_STATE, b"swmm_engine_open: engine is not in CREATED or CLOSED state");
            return SWMM_ERR_WRONG_STATE;
        }

        // Reset context for a fresh run
        self.ctx_.reset();
        self.metadata_only_ = metadata_only;

        // Stamp the report wall clock before any parsing work (std::time).
        self.ctx_.wall_start = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs() as i64).unwrap_or(0);

        self.rpt_path_ = rpt_path.unwrap_or(b"").to_vec();
        self.out_path_ = out_path.unwrap_or(b"").to_vec();
        let inp: &[u8] = inp_path.unwrap_or(b"");

        // Resolve input plugin: path, id:version, or fall back to default
        if self.plugins_.input_plugins().is_empty() {
            if let Some(lib) = input_plugin_lib {
                if !lib.is_empty() {
                    let mut warned: Vec<Bytes> = Vec::new();
                    let info = self.plugins_.find_component(lib, Some(&mut |m: &[u8]| warned.push(m.to_vec())));
                    for m in &warned {
                        self.emit_warning(crate::core::engine::SWMM_ERR_PLUGIN, m);
                    }
                    // A built-in info's create_input_plugin() returns null: nothing added.
                    let _ = info;
                }
            }
            // Fall back to built-in .inp reader
            if self.plugins_.input_plugins().is_empty() {
                let mut ip = DefaultInputPlugin::new();
                ip.initialize(&[]);
                self.plugins_.add_input_plugin(ip);
            }
        }

        if self.plugins_.input_plugins_mut()[0].read(inp, &mut self.ctx_) != 0 {
            return if self.ctx_.error_code != 0 { self.ctx_.error_code } else { SWMM_ERR_PARSE };
        }

        // Species registry: pollutants occupy the first slots.
        self.ctx_.species_registry.clear();
        for p in 0..self.ctx_.n_pollutants() {
            let name = self.ctx_.pollutant_names.name_of(p).to_vec();
            self.ctx_.species_registry.add(name, SpeciesKind::Pollutant, Vec::new());
        }

        // Resolve [PROCESS_COMPONENTS] registrations, then the embedded fallback, the ARD
        // transport rows and the per-element heat attributes.
        {
            let mut base_dir: Bytes = Vec::new();
            if !inp.is_empty() {
                base_dir = FsPath::new(inp).parent_path().native().to_vec();
            }

            let mut errs: Vec<Bytes> = Vec::new();
            if !self.ctx_.process_component_specs.is_empty() {
                errs = crate::core::process_components::resolve_process_components(&mut self.ctx_, &base_dir);
            }

            let mut reactions_registered = false;
            for spec in &self.ctx_.process_component_specs {
                if spec.id == b"org.hydrocouple.openswmm.reactions" {
                    reactions_registered = true;
                }
            }
            crate::transport::reactions_component::apply_embedded_reaction_sections(&mut self.ctx_, reactions_registered, &mut errs);

            crate::transport::ard_config::resolve_ard_transport_rows(&mut self.ctx_, &mut errs);

            {
                let pe = crate::transport::heat_overrides::resolve_heat_overrides(&mut self.ctx_);
                for e in pe {
                    errs.push(e);
                }
            }

            if !errs.is_empty() {
                for e in &errs {
                    self.ctx_.errors.push(e.clone());
                }
                if !self.lenient_open_ || self.metadata_only_ {
                    self.ctx_.error_code = SWMM_ERR_PARSE;
                    self.ctx_.error_message = errs[0].clone();
                    self.write_open_failure_report();
                    return SWMM_ERR_PARSE;
                }
            }

            crate::transport::reaction_legacy_binding::warn_if_legacy_binding_bypassed(&mut self.ctx_);
            crate::transport::ard_config::warn_if_fv_dispersion_key_ignored(&mut self.ctx_);

            if self.ctx_.heat_config.sediment_exchange && !self.ctx_.heat_config.sediment.has_ground_temp {
                self.ctx_.warnings.push(
                    "[SEDIMENT_EXCHANGE] GROUND_TEMPERATURE is not set; the default 12 degC is in force. For a buried conduit this term is often the largest one in the energy balance, so the default is unlikely to be what you meant."
                        .as_bytes()
                        .to_vec(),
                );
            }

            if self.ctx_.options.quality_solver != QualitySolverKind::Lagrangian && !self.ctx_.options.ignore_quality {
                if self.ctx_.options.quality_step > 0.0 {
                    self.ctx_.warnings.push(
                        "[OPTIONS] QUALITY_STEP is set but QUALITY_SOLVER is not LAGRANGIAN \u{2014} only the LARD engine substeps on this key; it has no effect this simulation."
                            .as_bytes()
                            .to_vec(),
                    );
                }
                if self.ctx_.options.max_segments_per_link != 100 {
                    self.ctx_.warnings.push(
                        "[OPTIONS] MAX_SEGMENTS_PER_LINK is set but QUALITY_SOLVER is not LAGRANGIAN \u{2014} only the LARD engine holds segments; it has no effect this simulation."
                            .as_bytes()
                            .to_vec(),
                    );
                }
                if self.ctx_.options.lard_rwpt {
                    self.ctx_.warnings.push(
                        "[OPTIONS] DISPERSION RWPT is set but QUALITY_SOLVER is not LAGRANGIAN \u{2014} RWPT is the LARD engine's dispersion; the ARD engine reads dispersion from the transport.ard component instead. No dispersion is applied from this key this simulation."
                            .as_bytes()
                            .to_vec(),
                    );
                }
            }
            // D-NS1 (X6): a configured negative quality inflow is EXTRACTION.
            for i in 0..self.ctx_.ext_inflows.count() {
                let ui = i as usize;
                if self.ctx_.ext_inflows.inflow_type[ui] == b"FLOW" {
                    continue;
                }
                if self.ctx_.ext_inflows.baseline[ui] >= 0.0 {
                    continue;
                }
                let m = cat(&[
                    b"[INFLOWS] node '",
                    &self.ctx_.ext_inflows.node_name[ui],
                    b"' pollutant '",
                    &self.ctx_.ext_inflows.constituent[ui],
                    "': the baseline is negative \u{2014} this row EXTRACTS mass (D-NS1). Extraction clamps per step to the mass the element holds.".as_bytes(),
                ]);
                self.ctx_.warnings.push(m);
            }

            // A2b: the REPORTED species list.
            self.ctx_.reported_species_names.clear();
            for p in 0..self.ctx_.n_pollutants() {
                let name = self.ctx_.pollutant_names.name_of(p).to_vec();
                self.ctx_.reported_species_names.push(name);
            }

            if self.ctx_.options.water_age {
                self.ctx_.species_registry.add(b"__WATER_AGE__".to_vec(), SpeciesKind::ReservedAge, b"hours".to_vec());
                self.ctx_.reported_species_names.push(b"__WATER_AGE__".to_vec());
                if self.ctx_.options.ignore_quality {
                    self.ctx_.warnings.push(
                        "[OPTIONS] WATER_AGE ON but IGNORE_QUALITY is YES \u{2014} the quality stage does not run, so no age is tracked this simulation."
                            .as_bytes()
                            .to_vec(),
                    );
                }
            }

            if self.ctx_.options.heat_transport {
                self.ctx_.species_registry.add(b"__TEMPERATURE__".to_vec(), SpeciesKind::ReservedTemperature, b"degC".to_vec());
                self.ctx_.reported_species_names.push(b"__TEMPERATURE__".to_vec());
                if self.ctx_.options.ignore_quality {
                    self.ctx_.warnings.push(
                        "[OPTIONS] HEAT_TRANSPORT ON but IGNORE_QUALITY is YES \u{2014} the quality stage does not run, so no temperature is tracked this simulation."
                            .as_bytes()
                            .to_vec(),
                    );
                }
            }

            // E1: reactions-component species report trailing the reserved rows.
            self.ctx_.refresh_reported_species();
        }

        // Warn about unknown/skipped sections (legacy warns once per pass, twice in all).
        if !self.ctx_.unknown_section_headers.is_empty() {
            let mut msgs: Vec<Bytes> = Vec::new();
            for (tok, line) in &self.ctx_.unknown_section_headers {
                msgs.push(cat(&[b"\n  WARNING: Unknown section '", tok, b"' at line ", line.to_string().as_bytes(), b" will be skipped."]));
            }
            let tail = std::mem::take(&mut self.ctx_.warnings);
            self.ctx_.warnings = msgs.clone();
            self.ctx_.warnings.extend(tail);
            for m in &msgs {
                self.emit_warning(100, m);
            }
            for m in &msgs {
                self.push_report_warning(m, 100);
            }
        } else {
            for tag in self.plugins_.input_plugins()[0].skipped_sections() {
                let m = cat(&[b"WARNING: Unknown section '[", &tag, b"]' will be skipped."]);
                self.push_report_warning(&m, 100);
            }
        }

        // Resolve cross-references (forward refs, final array sizing, head init)
        crate::input::post_parse_resolver::resolve_cross_references(&mut self.ctx_);
        // BW-MSX: bind the [BUILDUP]/[WASHOFF]/[LOADINGS] rows that name a species.
        for w in crate::transport::msx_surface::resolve(&mut self.ctx_) {
            self.ctx_.warnings.push(w);
        }

        // Project-level sanity checks (legacy project_validate), only without input errors.
        if self.ctx_.errors.is_empty() {
            self.validate_project();
        }

        // Post-parse validation errors are fatal (strict open).
        if !self.ctx_.errors.is_empty() {
            self.write_open_failure_report();
            if !self.lenient_open_ || self.metadata_only_ {
                let first = self.ctx_.errors[0].clone();
                self.set_error(SWMM_ERR_PARSE, &first);
                return SWMM_ERR_PARSE;
            }
            // Lenient open: keep the errors and fall through to OPENED.
        }

        // Metadata discovery stops before runtime plugin factories and output registration.
        if self.metadata_only_ {
            self.ctx_.state = EngineState::Opened;
            return SWMM_OK;
        }

        // Phase 4: load plugins listed in [PLUGINS]
        if !self.ctx_.plugin_specs.is_empty() {
            let mut warned: Vec<Bytes> = Vec::new();
            let specs = self.ctx_.plugin_specs.clone();
            self.plugins_.load_plugins(&specs, &mut |m: &[u8]| warned.push(m.to_vec()), Some(&self.ctx_));
            for m in &warned {
                self.emit_warning(SWMM_ERR_PLUGIN, m);
            }
        }

        // Explicit per-role writer choices are independent of the plugin catalogue.
        if !self.selected_writer(false) || !self.selected_writer(true) {
            self.set_error(SWMM_ERR_PLUGIN, b"The selected output/report writer or destination is unavailable");
            return SWMM_ERR_PLUGIN;
        }
        // Native companions (IO_OUTPUT_WRITER / IO_REPORT_WRITER naming a plugin): reached
        // only after selected_writer accepted a plugin, which needs a loaded external library
        // (never, here). Upstream then renames out_path_/rpt_path_ to the .inp's .out/.rpt.

        // Inject built-in output/report plugins when paths are configured.
        if !self.out_path_.is_empty() {
            let mut op = DefaultOutputPlugin::new(self.out_path_.clone());
            op.initialize(&[]);
            op.validate(&self.ctx_);
            self.plugins_.add_output_plugin(op);
            self.binary_output_plugin_ = Some(self.plugins_.output_plugins().len() - 1);
        }
        if !self.rpt_path_.is_empty() {
            let mut rp = DefaultReportPlugin::new(self.rpt_path_.clone());
            rp.initialize(&[]);
            rp.validate(&self.ctx_);
            self.plugins_.add_report_plugin(rp);
        }

        // ctx_.state_accessors: NOT YET TRANSLATED (see the file header).

        // Inject the default state-IO plugin if none resolved during [PLUGINS] processing.
        {
            let has_default_state_io = !self.plugins_.state_io_plugins().is_empty();
            if !has_default_state_io {
                let mut sp = DefaultStateIOPlugin::new();
                sp.initialize(&[]);
                sp.validate(&self.ctx_);
                self.plugins_.add_state_io_plugin(sp);
            }
        }

        self.ctx_.state = EngineState::Opened;
        SWMM_OK
    }

    /// open()'s selectedWriter lambda. An empty or BUILTIN choice keeps the built-in
    /// writer. Any other id must name a component that creates a writer: a built-in info's
    /// create_*_plugin() returns null, an unknown id is not found, and a library path is
    /// never loaded here, so every such choice is refused, as upstream refuses it when the
    /// library does not load. NOT TRANSLATED: a loadable external writer library.
    fn selected_writer(&mut self, report: bool) -> bool {
        let key: &[u8] = if report { b"IO_REPORT_WRITER" } else { b"IO_OUTPUT_WRITER" };
        let id = self.ctx_.options.ext_options.get(key).cloned().unwrap_or_default();
        if id.is_empty() || id == b"BUILTIN" {
            return true;
        }
        let info = self.plugins_.find_component(&id, None);
        let info = match info {
            Some(i) => i,
            None => return false,
        };
        if if report { !info.has_report } else { !info.has_output } {
            return false;
        }
        // The destination is computed (and a native deck guarded against overwrite) before
        // the writer is created; a built-in info creates no writer.
        false
    }

    // ========================================================================
    // Private helpers
    // ========================================================================

    pub fn validate_project(&mut self) {
        // ERROR 108: an outlet name matching both a node and a subcatchment.
        for s in 0..self.ctx_.n_subcatches() {
            let us = s as usize;
            if us >= self.ctx_.subcatches.outlet_name.len() {
                continue;
            }
            let name = &self.ctx_.subcatches.outlet_name[us];
            if !name.is_empty() && self.ctx_.node_names.find(name) >= 0 && self.ctx_.subcatch_names.find(name) >= 0 {
                let e = format_error(ERR_SUBCATCH_OUTLET, self.ctx_.subcatch_names.name_of(s));
                self.ctx_.errors.push(e);
            }
        }

        // Legacy gage_validate checks only a used TIMESERIES gage that is not a co-gage.
        let ng = self.ctx_.n_gages();
        for g in 0..ng {
            let ug = g as usize;
            if !crate::hydrology::gage::gage_is_used(&self.ctx_, g) {
                continue;
            }
            if self.ctx_.gages.source[ug] != RainSource::Timeseries {
                continue;
            }
            if self.ctx_.gages.co_gage_index[ug] >= 0 {
                continue;
            }

            // ERROR 159 / WARNING 09: recording interval against the series' smallest spacing.
            let ts = self.ctx_.gages.ts_index[ug];
            if ts >= 0 && ts < self.ctx_.tables.tables.len() as i32 {
                let tx = &self.ctx_.tables.tables[ts as usize].x;
                if tx.len() > 1 {
                    let mut dx_min = tx[1] - tx[0];
                    for k in 2..tx.len() {
                        dx_min = cmin(dx_min, tx[k] - tx[k - 1]);
                    }
                    let ts_interval = cvt_i32((dx_min * 86400.0 + 0.5).floor());
                    if ts_interval > 0 && self.ctx_.gages.interval_sec[ug] > ts_interval {
                        let e = format_error(ERR_RAIN_GAGE_INTERVAL, self.ctx_.gage_names.name_of(g));
                        self.ctx_.errors.push(e);
                    }
                    if self.ctx_.gages.interval_sec[ug] < ts_interval {
                        let w = format_warning(WARN_GAGE_INTERVAL, self.ctx_.gage_names.name_of(g));
                        self.ctx_.warnings.push(w);
                    }
                }
            }

            let interval = self.ctx_.gages.interval_sec[ug] as f64;
            if interval > 0.0 && self.ctx_.options.wet_step > interval {
                self.ctx_.options.wet_step = interval;
                let w = format_warning(WARN_WET_STEP_REDUCED, self.ctx_.gage_names.name_of(g));
                self.ctx_.warnings.push(w);
            }
        }

        // WARNING 06: dry-weather step increased to the wet-weather step.
        let opt = &mut self.ctx_.options;
        if opt.dry_step < opt.wet_step {
            opt.dry_step = opt.wet_step;
            self.ctx_.warnings.push(format_warning(WARN_DRY_STEP_INCREASED, b""));
        }

        // WARNING 07: routing step reduced to the wet-weather step; the authored step is kept.
        let opt = &mut self.ctx_.options;
        opt.routing_step_authored = opt.routing_step;
        if opt.routing_step > opt.wet_step {
            opt.routing_step = opt.wet_step;
            self.ctx_.warnings.push(format_warning(WARN_ROUTING_STEP_REDUCED, b""));
        }
    }

    /// Writes the .rpt of a failed open. NOT YET TRANSLATED (the .rpt writer): no file is
    /// written. The context is not changed either way (prepare/finalize take it const).
    fn write_open_failure_report(&mut self) {
        if self.rpt_path_.is_empty() || self.ctx_.options.rpt_disabled {
            return;
        }
    }

    /// `message` is a C string upstream (callers pass `.c_str()`), so it ends at a NUL.
    pub fn set_error(&mut self, code: i32, message: &[u8]) {
        let message = crate::core::cstdlib::c_str(message);
        self.ctx_.error_code = code;
        self.ctx_.error_message = message.to_vec();
        self.ctx_.state = EngineState::ErrorState;
        self.emit_warning(code, message);
    }

    pub fn emit_warning(&mut self, code: i32, message: &[u8]) {
        if let Some(cb) = self.callbacks_.on_warning.as_mut() {
            cb(code, crate::core::cstdlib::c_str(message));
        }
    }

    pub fn push_report_warning(&mut self, message: &[u8], code: i32) {
        self.ctx_.warnings.push(message.to_vec());
        self.emit_warning(code, message);
    }
}
