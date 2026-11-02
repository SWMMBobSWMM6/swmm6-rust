// SPDX-License-Identifier: Apache-2.0
//
// The section handlers (upstream src/engine/input/handlers/*.cpp) and the handler table of
// DefaultInputPlugin::register_builtin_handlers. Each handler module carries its own
// derivation notice.

use crate::input::section_registry::SectionRegistry;

/// Register the handlers named in `on` (or all with `all`), as DefaultInputPlugin does,
/// without the 2D sections. Filled in as handlers are translated.
pub fn register_reference_handlers(reg: &mut SectionRegistry, on: &[String], all: bool) {
    let _ = (reg, on, all);
}
