// The handler table of DefaultInputPlugin::register_builtin_handlers, without the 2D sections.
#include "input/SectionRegistry.hpp"
#include "input/handlers/TitleHandler.hpp"
#include <set>
#include <string>
namespace openswmm::input {
void register_reference_handlers(SectionRegistry& reg, const std::set<std::string>& on, bool all) {
    (void)reg; (void)on; (void)all;   // filled in as the Rust handlers arrive
}
}
