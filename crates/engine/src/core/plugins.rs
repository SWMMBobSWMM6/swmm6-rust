// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/plugins/DefaultInputPlugin.cpp,
// PluginFactory.cpp, BuiltinPluginInfos.hpp, DefaultOutputPlugin.cpp,
// DefaultReportPlugin.cpp and DefaultStateIOPlugin.cpp, and
// include/openswmm/plugin_sdk/PluginState.hpp (HydroCouple/openswmm.engine, branch swmm6_rel,
// commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The part of the plugin layer that SWMMEngine::open reaches in a 1D build:
//   * DefaultInputPlugin: the handler table and read().
//   * PluginFactory: the four built-in component infos (register_builtin_infos),
//     find_component, load_plugins, and the per-role plugin lists open() appends to.
//   * The built-in output, report and state-IO plugins up to validate(); writing .out and
//     .rpt files is not translated here.
//
// Differences of form only:
//   * The plugin lists hold the concrete built-in plugin types (a 1D open creates no other).
//     Upstream's lists hold interface pointers.
//   * A built-in info's create_*_plugin() returns null upstream (IPluginComponentInfo's
//     default); here it is simply absent, with the same effect: nothing is created.
//   * PluginFactory's constructor also runs discover(), which scans the engine library's
//     own directory for plugin shared libraries. Safe Rust cannot load a shared library, so
//     no library is ever loaded: discovery finds nothing, and a library path in [PLUGINS],
//     in an input-plugin argument or in IO_*_WRITER behaves as a library that failed to
//     load (upstream's warning text then carries dlerror(), here a fixed note). NOT
//     TRANSLATED: loading an external plugin library.

use std::collections::BTreeMap;

use crate::core::context::{PluginSpec, SimulationContext};
use crate::core::options::Bytes;
use crate::input::input_reader::InputReader;
use crate::input::section_registry::SectionRegistry;

/// include/openswmm/plugin_sdk/PluginState.hpp
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i32)]
pub enum PluginState {
    #[default]
    Unloaded = 0,
    Loaded = 1,
    Initialized = 2,
    Validated = 3,
    Prepared = 4,
    Updating = 5,
    Finalized = 6,
    Closed = 7,
    Error = -1,
}

// ============================================================================
// DefaultInputPlugin
// ============================================================================

pub struct DefaultInputPlugin {
    registry_: SectionRegistry,
    state_: PluginState,
    last_error_: Bytes,
    skipped_sections_: Vec<Bytes>,
}

impl Default for DefaultInputPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl DefaultInputPlugin {
    /// Construct an UNLOADED plugin and register the built-in handlers.
    pub fn new() -> Self {
        let mut p = DefaultInputPlugin {
            registry_: SectionRegistry::new(),
            state_: PluginState::Unloaded,
            last_error_: Vec::new(),
            skipped_sections_: Vec::new(),
        };
        p.register_builtin_handlers();
        p
    }

    /// register_builtin_handlers(): the same set of tags, handlers and passthrough
    /// captures, kept in input::handlers::register_reference_handlers (the 1D table;
    /// dispatch does not depend on registration order).
    fn register_builtin_handlers(&mut self) {
        crate::input::handlers::register_reference_handlers(&mut self.registry_, &[], true);
    }

    pub fn state(&self) -> PluginState {
        self.state_
    }

    pub fn initialize(&mut self, _init_args: &[Bytes]) -> i32 {
        self.state_ = PluginState::Initialized;
        0
    }

    pub fn validate(&mut self, _ctx: &SimulationContext) -> i32 {
        self.state_ = PluginState::Validated;
        0
    }

    pub fn read(&mut self, path: &[u8], ctx: &mut SimulationContext) -> i32 {
        let mut reader = InputReader::new(&self.registry_);
        let ok = reader.read(path, ctx);

        self.skipped_sections_ = reader.skipped_sections().to_vec();

        if !ok {
            self.last_error_ = ctx.error_message.clone();
            return if ctx.error_code != 0 { ctx.error_code } else { -1 };
        }
        0
    }

    pub fn finalize(&mut self, _ctx: &SimulationContext) -> i32 {
        self.state_ = PluginState::Finalized;
        0
    }

    pub fn last_error_message(&self) -> &[u8] {
        &self.last_error_
    }

    pub fn skipped_sections(&self) -> Vec<Bytes> {
        self.skipped_sections_.clone()
    }

    pub fn registry(&mut self) -> &mut SectionRegistry {
        &mut self.registry_
    }
}

// ============================================================================
// DefaultOutputPlugin / DefaultReportPlugin / DefaultStateIOPlugin, up to validate()
// ============================================================================

pub struct DefaultOutputPlugin {
    out_path_: Bytes,
    state_: PluginState,
    step_count_: i64,
    n_periods_: i32,
}

impl DefaultOutputPlugin {
    pub fn new(out_path: Bytes) -> Self {
        DefaultOutputPlugin { out_path_: out_path, state_: PluginState::Loaded, step_count_: 0, n_periods_: 0 }
    }
    pub fn initialize(&mut self, _init_args: &[Bytes]) -> i32 {
        self.step_count_ = 0;
        self.n_periods_ = 0;
        self.state_ = PluginState::Initialized;
        0
    }
    pub fn validate(&mut self, _ctx: &SimulationContext) -> i32 {
        self.state_ = PluginState::Validated;
        0
    }
    pub fn state(&self) -> PluginState {
        self.state_
    }
    pub fn out_path(&self) -> &[u8] {
        &self.out_path_
    }
}

pub struct DefaultReportPlugin {
    rpt_path_: Bytes,
    state_: PluginState,
}

impl DefaultReportPlugin {
    pub fn new(rpt_path: Bytes) -> Self {
        DefaultReportPlugin { rpt_path_: rpt_path, state_: PluginState::Loaded }
    }
    pub fn initialize(&mut self, _init_args: &[Bytes]) -> i32 {
        self.state_ = PluginState::Initialized;
        0
    }
    pub fn validate(&mut self, _ctx: &SimulationContext) -> i32 {
        self.state_ = PluginState::Validated;
        0
    }
    pub fn state(&self) -> PluginState {
        self.state_
    }
    pub fn rpt_path(&self) -> &[u8] {
        &self.rpt_path_
    }
}

pub struct DefaultStateIOPlugin {
    state_: PluginState,
    last_error_: Bytes,
}

impl Default for DefaultStateIOPlugin {
    fn default() -> Self {
        DefaultStateIOPlugin { state_: PluginState::Unloaded, last_error_: Vec::new() }
    }
}

impl DefaultStateIOPlugin {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn initialize(&mut self, _init_args: &[Bytes]) -> i32 {
        self.state_ = PluginState::Initialized;
        self.last_error_.clear();
        0
    }
    pub fn validate(&mut self, _ctx: &SimulationContext) -> i32 {
        self.state_ = PluginState::Validated;
        0
    }
    pub fn state(&self) -> PluginState {
        self.state_
    }
}

// ============================================================================
// Built-in component infos (BuiltinPluginInfos.hpp)
// ============================================================================

/// The four built-in infos: id, version and capabilities. Their create_*_plugin() is
/// IPluginComponentInfo's default, which returns null.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BuiltinInfo {
    pub id: &'static str,
    pub version: &'static str,
    pub has_input: bool,
    pub has_output: bool,
    pub has_report: bool,
    pub has_state_io: bool,
}

pub const BUILTIN_INFOS: [BuiltinInfo; 4] = [
    BuiltinInfo { id: "org.hydrocouple.openswmm.builtin.input", version: "1.0.0", has_input: true, has_output: false, has_report: false, has_state_io: false },
    BuiltinInfo { id: "org.hydrocouple.openswmm.builtin.output", version: "1.0.0", has_input: false, has_output: true, has_report: false, has_state_io: false },
    BuiltinInfo { id: "org.hydrocouple.openswmm.builtin.report", version: "1.0.0", has_input: false, has_output: false, has_report: true, has_state_io: false },
    BuiltinInfo { id: "org.hydrocouple.openswmm.builtin.state_io", version: "1.0.0", has_input: false, has_output: false, has_report: false, has_state_io: true },
];

// ============================================================================
// PluginFactory
// ============================================================================

pub struct PluginFactory {
    /// "id:version" and "id" -> index into BUILTIN_INFOS (libs_ holds only built-ins here)
    registry_: BTreeMap<Bytes, usize>,
    input_plugins_: Vec<DefaultInputPlugin>,
    output_plugins_: Vec<DefaultOutputPlugin>,
    report_plugins_: Vec<DefaultReportPlugin>,
    state_io_plugins_: Vec<DefaultStateIOPlugin>,
}

impl Default for PluginFactory {
    fn default() -> Self {
        Self::new()
    }
}

impl PluginFactory {
    pub fn new() -> Self {
        let mut f = PluginFactory {
            registry_: BTreeMap::new(),
            input_plugins_: Vec::new(),
            output_plugins_: Vec::new(),
            report_plugins_: Vec::new(),
            state_io_plugins_: Vec::new(),
        };
        f.register_builtin_infos();
        // discover(): finds nothing (no library is loaded; see the file header)
        f
    }

    fn register_builtin_infos(&mut self) {
        for (idx, info) in BUILTIN_INFOS.iter().enumerate() {
            let id = info.id.as_bytes().to_vec();
            let mut key = id.clone();
            key.push(b':');
            key.extend_from_slice(info.version.as_bytes());
            if self.registry_.contains_key(&key) {
                continue; // idempotent
            }
            self.registry_.insert(key, idx);
            self.registry_.entry(id).or_insert(idx);
        }
    }

    /// is_file_path(): a separator, or the POSIX shared-library suffix ".so".
    pub fn is_file_path(s: &[u8]) -> bool {
        if s.is_empty() {
            return false;
        }
        if s.contains(&b'/') || s.contains(&b'\\') {
            return true;
        }
        Self::is_shared_library(s)
    }

    pub fn is_shared_library(filename: &[u8]) -> bool {
        if filename.len() < 3 {
            return false;
        }
        filename.ends_with(b".so")
    }

    /// load_library(): upstream dlopen()s the path. NOT TRANSLATED: no library is ever
    /// loaded here, so this is the dlopen-failure branch.
    fn load_library(&mut self, path: &[u8], warn_cb: &mut dyn FnMut(&[u8])) -> Option<BuiltinInfo> {
        let mut m = b"PluginFactory: failed to load '".to_vec();
        m.extend_from_slice(path);
        m.extend_from_slice(b"': (swmm6-rust does not load plugin libraries)");
        warn_cb(&m);
        None
    }

    pub fn find_component(&mut self, id_or_path: &[u8], warn_cb: Option<&mut dyn FnMut(&[u8])>) -> Option<BuiltinInfo> {
        if id_or_path.is_empty() {
            return None;
        }
        let mut noop = |_: &[u8]| {};
        let has_cb = warn_cb.is_some();
        let cb: &mut dyn FnMut(&[u8]) = match warn_cb {
            Some(c) => c,
            None => &mut noop,
        };
        if Self::is_file_path(id_or_path) {
            return self.load_library(id_or_path, cb);
        }
        if let Some(&i) = self.registry_.get(id_or_path) {
            return Some(BUILTIN_INFOS[i]);
        }
        if has_cb {
            let mut m = b"PluginFactory: plugin '".to_vec();
            m.extend_from_slice(id_or_path);
            m.extend_from_slice(b"' not found in registry or on disk");
            cb(&m);
        }
        None
    }

    /// load_plugins(): a built-in info creates nothing (its create_* returns null).
    pub fn load_plugins(&mut self, specs: &[PluginSpec], warn_cb: &mut dyn FnMut(&[u8]), _io: Option<&SimulationContext>) -> i32 {
        let mut loaded = 0;
        for spec in specs {
            let info = self.find_component(&spec.path, Some(&mut *warn_cb));
            if info.is_none() {
                continue;
            }
            // has_input / has_output / has_report / has_state_io: create_*_plugin() is null
            loaded += 1;
        }
        loaded
    }

    pub fn input_plugins(&self) -> &[DefaultInputPlugin] {
        &self.input_plugins_
    }
    pub fn input_plugins_mut(&mut self) -> &mut Vec<DefaultInputPlugin> {
        &mut self.input_plugins_
    }
    pub fn output_plugins(&self) -> &[DefaultOutputPlugin] {
        &self.output_plugins_
    }
    pub fn report_plugins(&self) -> &[DefaultReportPlugin] {
        &self.report_plugins_
    }
    pub fn state_io_plugins(&self) -> &[DefaultStateIOPlugin] {
        &self.state_io_plugins_
    }

    pub fn add_input_plugin(&mut self, plugin: DefaultInputPlugin) {
        self.input_plugins_.push(plugin);
    }
    pub fn add_output_plugin(&mut self, plugin: DefaultOutputPlugin) {
        self.output_plugins_.push(plugin);
    }
    pub fn add_report_plugin(&mut self, plugin: DefaultReportPlugin) {
        self.report_plugins_.push(plugin);
    }
    pub fn add_state_io_plugin(&mut self, plugin: DefaultStateIOPlugin) {
        self.state_io_plugins_.push(plugin);
    }
}

// Gate only (cargo feature `gate`): set the state the engine dump prints, for the context
// loader in examples/common/ctxload.rs. Not part of the translation.
#[cfg(feature = "gate")]
impl DefaultInputPlugin {
    pub fn gate_set(&mut self, state: PluginState, last_error: Bytes, skipped_sections: Vec<Bytes>) {
        self.state_ = state;
        self.last_error_ = last_error;
        self.skipped_sections_ = skipped_sections;
    }
}
#[cfg(feature = "gate")]
impl DefaultOutputPlugin {
    pub fn gate_set_state(&mut self, state: PluginState) {
        self.state_ = state;
    }
}
#[cfg(feature = "gate")]
impl DefaultReportPlugin {
    pub fn gate_set_state(&mut self, state: PluginState) {
        self.state_ = state;
    }
}
#[cfg(feature = "gate")]
impl DefaultStateIOPlugin {
    pub fn gate_set_state(&mut self, state: PluginState) {
        self.state_ = state;
    }
}
