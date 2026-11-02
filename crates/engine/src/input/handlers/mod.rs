// SPDX-License-Identifier: Apache-2.0
//
// The section handlers (upstream src/engine/input/handlers/*.cpp) and the handler table of
// DefaultInputPlugin::register_builtin_handlers. Each handler module carries its own
// derivation notice.

pub mod catchment_handler;
pub mod hydrology_handler;
pub mod infra_handler;
pub mod inflows_handler;
pub mod links_handler;
pub mod nodes_handler;
pub mod quality_handler;
pub mod tables_handler;

use crate::core::context::SimulationContext;
use crate::core::options::Bytes;
use crate::input::section_registry::SectionRegistry;

/// Register the handlers named in `on` (or all with `all`), as DefaultInputPlugin does,
/// without the 2D sections. Filled in as handlers are translated.
pub fn register_reference_handlers(reg: &mut SectionRegistry, on: &[String], all: bool) {
    let mut r = |tag: &str, h: fn(&mut SimulationContext, &[Bytes])| {
        if all || on.iter().any(|s| s == tag) {
            reg.register_builtin(tag.as_bytes(), h);
        }
    };

    // Node sections
    r("JUNCTIONS", nodes_handler::handle_junctions);
    r("VIRTUAL_JUNCTIONS", nodes_handler::handle_virtual_junctions);
    r("INLET_JUNCTIONS", infra_handler::handle_inlet_junctions);
    r("OUTFALLS", nodes_handler::handle_outfalls);
    r("DIVIDERS", nodes_handler::handle_dividers);
    r("STORAGE", nodes_handler::handle_storage);
    r("COORDINATES", nodes_handler::handle_coordinates);

    // Link sections
    r("CONDUITS", links_handler::handle_conduits);
    r("PUMPS", links_handler::handle_pumps);
    r("ORIFICES", links_handler::handle_orifices);
    r("WEIRS", links_handler::handle_weirs);
    r("OUTLETS", links_handler::handle_outlets);
    r("XSECTIONS", links_handler::handle_xsections);
    r("TRANSECTS", links_handler::handle_transects);
    r("LOSSES", links_handler::handle_losses);

    // Inflows / loading
    r("INFLOWS", inflows_handler::handle_inflows);
    r("DWF", inflows_handler::handle_dwf);
    r("RDII", inflows_handler::handle_rdii);
    r("HYDROGRAPHS", inflows_handler::handle_hydrographs);
    r("RDII_DECAY", inflows_handler::handle_rdii_decay);
    r("PATTERNS", inflows_handler::handle_patterns);

    // Tables
    r("CURVES", tables_handler::handle_curves);
    r("TIMESERIES", tables_handler::handle_timeseries);

    // Control / map / reporting
    r("ADJUSTMENTS", infra_handler::handle_adjustments);
    r("EVENTS", infra_handler::handle_events);

    // Streets / inlets
    r("STREETS", infra_handler::handle_streets);
    r("INLETS", infra_handler::handle_inlets);
    r("INLET_USAGE", infra_handler::handle_inlet_usage);

    // Subcatchments, hydrology, quality
    r("SUBCATCHMENTS", catchment_handler::handle_subcatchments);
    r("SUBAREAS", catchment_handler::handle_subareas);
    r("INFILTRATION", catchment_handler::handle_infiltration);
    r("LID_CONTROLS", hydrology_handler::handle_lid_controls);
    r("LID_USAGE", hydrology_handler::handle_lid_usage);
    r("AQUIFERS", hydrology_handler::handle_aquifers);
    r("GROUNDWATER", hydrology_handler::handle_groundwater);
    r("GWF", hydrology_handler::handle_gwf);
    r("RAINGAGES", catchment_handler::handle_raingages);
    r("EVAPORATION", hydrology_handler::handle_evaporation);
    r("TEMPERATURE", hydrology_handler::handle_temperature);
    r("SNOWPACKS", hydrology_handler::handle_snowpacks);
    r("UEB_OPTIONS", hydrology_handler::handle_ueb_options);
    r("2D_SNOW_OPTIONS", hydrology_handler::handle_mesh_snow_options);
    r("2D_SNOW_DEFAULTS", hydrology_handler::handle_mesh_snow_defaults);
    r("2D_SNOW", hydrology_handler::handle_mesh_snow_cells);
    r("LOADINGS", quality_handler::handle_loadings);
    r("POLLUTANTS", quality_handler::handle_pollutants);
    r("LANDUSES", quality_handler::handle_landuses);
    r("COVERAGES", quality_handler::handle_coverages);
    r("BUILDUP", quality_handler::handle_buildup);
    r("WASHOFF", quality_handler::handle_washoff);
    r("TREATMENT", quality_handler::handle_treatment);
    r("INITIAL_QUALITY", quality_handler::handle_initial_quality);
}
