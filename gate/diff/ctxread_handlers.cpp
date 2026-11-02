// The handler table of DefaultInputPlugin::register_builtin_handlers, without the 2D sections.
#include "input/SectionRegistry.hpp"
#include "input/handlers/TitleHandler.hpp"
#include "input/handlers/NodesHandler.hpp"
#include "input/handlers/LinksHandler.hpp"
#include "input/handlers/InfraHandler.hpp"
#include "input/handlers/TablesHandler.hpp"
#include "input/handlers/InflowsHandler.hpp"
#include "input/handlers/CatchmentHandler.hpp"
#include "input/handlers/HydrologyHandler.hpp"
#include "input/handlers/QualityHandler.hpp"
#include "input/handlers/OptionsHandler.hpp"
#include "input/handlers/FilesHandler.hpp"
#include "input/handlers/PluginsHandler.hpp"
#include "input/handlers/ProcessComponentsHandler.hpp"
#include "input/handlers/UserFlagsHandler.hpp"
#include "input/handlers/UserFlagValuesHandler.hpp"
#include "input/handlers/SpatialHandler.hpp"
#include "input/handlers/ControlsHandler.hpp"
#include "core/SimulationContext.hpp"
#include "transport/components/ReactionModule/ReactionsComponent.hpp"
#include <set>
#include <string>
namespace openswmm::input {
void register_reference_handlers(SectionRegistry& reg, const std::set<std::string>& on, bool all) {
    auto r = [&](const char* tag, auto h) { if (all || on.count(tag)) reg.register_builtin(tag, h); };
    // DefaultInputPlugin.cpp's capture lambda (local to that file), copied.
    auto capture = [](const char* tag) {
        return [tag](SimulationContext& ctx, const std::vector<std::string>& lines) {
            auto& dst = ctx.passthrough_sections[tag];
            dst.insert(dst.end(), lines.begin(), lines.end());
        };
    };

    r("TITLE", input::handle_title);
    r("OPTIONS", input::handle_options);
    r("CONTROLS", input::handle_controls);
    r("VERTICES", input::handle_vertices);
    r("POLYGONS", input::handle_polygons);
    r("SYMBOLS", input::handle_symbols);
    r("LABELS", capture("LABELS"));
    r("BACKDROP", capture("BACKDROP"));
    r("MAP", input::handle_map);
    r("TAGS", input::handle_tags);
    r("PROFILE", capture("PROFILE"));
    r("PROFILES", capture("PROFILES"));
    r("REPORT", input::handle_report);
    r("FILES", input::handle_files);
    r("USER_FLAGS", input::handle_user_flags);
    r("USER_FLAG_VALUES", input::handle_user_flag_values);
    r("PLUGINS", input::handle_plugins);
    r("PROCESS_COMPONENTS", input::handle_process_components);

    // Node sections
    r("JUNCTIONS", input::handle_junctions);
    r("VIRTUAL_JUNCTIONS", input::handle_virtual_junctions);
    r("INLET_JUNCTIONS", input::handle_inlet_junctions);
    r("OUTFALLS", input::handle_outfalls);
    r("DIVIDERS", input::handle_dividers);
    r("STORAGE", input::handle_storage);
    r("COORDINATES", input::handle_coordinates);

    // Link sections
    r("CONDUITS", input::handle_conduits);
    r("PUMPS", input::handle_pumps);
    r("ORIFICES", input::handle_orifices);
    r("WEIRS", input::handle_weirs);
    r("OUTLETS", input::handle_outlets);
    r("XSECTIONS", input::handle_xsections);
    r("TRANSECTS", input::handle_transects);
    r("LOSSES", input::handle_losses);

    // Inflows / loading
    r("INFLOWS", input::handle_inflows);
    r("DWF", input::handle_dwf);
    r("RDII", input::handle_rdii);
    r("HYDROGRAPHS", input::handle_hydrographs);
    r("RDII_DECAY", input::handle_rdii_decay);
    r("PATTERNS", input::handle_patterns);

    // Tables
    r("CURVES", input::handle_curves);
    r("TIMESERIES", input::handle_timeseries);

    // Control / map / reporting
    r("ADJUSTMENTS", input::handle_adjustments);
    r("EVENTS", input::handle_events);

    // Streets / inlets (legacy s_STREET, s_INLET, s_INLET_USAGE)
    r("STREETS", input::handle_streets);
    r("INLETS", input::handle_inlets);
    r("INLET_USAGE", input::handle_inlet_usage);

    // Subcatchments, hydrology, quality
    r("SUBCATCHMENTS", input::handle_subcatchments);
    r("SUBAREAS", input::handle_subareas);
    r("INFILTRATION", input::handle_infiltration);
    r("LID_CONTROLS", input::handle_lid_controls);
    r("LID_USAGE", input::handle_lid_usage);
    r("AQUIFERS", input::handle_aquifers);
    r("GROUNDWATER", input::handle_groundwater);
    r("GWF", input::handle_gwf);
    r("RAINGAGES", input::handle_raingages);
    r("EVAPORATION", input::handle_evaporation);
    r("TEMPERATURE", input::handle_temperature);
    r("SNOWPACKS", input::handle_snowpacks);
    r("UEB_OPTIONS", input::handle_ueb_options);
    r("2D_SNOW_OPTIONS", input::handle_mesh_snow_options);
    r("2D_SNOW_DEFAULTS", input::handle_mesh_snow_defaults);
    r("2D_SNOW", input::handle_mesh_snow_cells);
    r("LOADINGS", input::handle_loadings);
    r("POLLUTANTS", input::handle_pollutants);
    r("LANDUSES", input::handle_landuses);
    r("COVERAGES", input::handle_coverages);
    r("BUILDUP", input::handle_buildup);
    r("WASHOFF", input::handle_washoff);
    r("TREATMENT", input::handle_treatment);
    r("INITIAL_QUALITY", input::handle_initial_quality);
    for (const auto& tag : transport::reactionSectionTags())
        if (all || on.count(tag))
            reg.register_builtin(tag, [tag](SimulationContext& ctx, const std::vector<std::string>& lines) {
                ctx.embedded_component_sections.emplace_back(tag, lines);
            });
}
}
