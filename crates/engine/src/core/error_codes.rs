// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/ErrorCodes.hpp and
// src/engine/core/ErrorCodes.cpp (HydroCouple/openswmm.engine, branch swmm6_rel,
// commit 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the Apache License,
// Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Differences in form, never in result:
//
//   * The four `enum X : int` types are passed around as `int`, so each is a type alias
//     for i32 and its enumerators are `pub const` items with upstream's names and values.
//   * Each `std::unordered_map<int, const char*>` table is a static slice of
//     (code, text) pairs in upstream's initializer order, searched by key. No key repeats,
//     so the search finds what `find` finds. Texts are byte strings, byte for byte
//     upstream's literals; adjacent C++ literals are joined with `\`-newline, which keeps
//     the space before the break and drops the indentation after it.
//   * The C++ overloads with a `detail` argument are `format_error_detail` and
//     `format_warning_detail`; the default `name = ""` is passed as `b""`.
//   * `std::snprintf` into `char buf[64]` / `char prefix[32]` is `format!`. No `int`
//     overflows those buffers (the longest results, for INT_MIN, are 44 and 23 bytes), so
//     snprintf never truncates. `%d` is `{}` and `%02d` is `{:02}`; both put the sign
//     inside the width, so -5 gives "-5" either way.

use crate::core::options::Bytes;

// ============================================================================
// Error codes — matching legacy error.h / error.txt numeric values
// ============================================================================

/// Upstream `enum ErrorCode : int`.
pub type ErrorCode = i32;

pub const ERR_NONE: ErrorCode = 0;

// --- System errors (100s) ---
/// memory allocation error
pub const ERR_MEMORY: ErrorCode = 101;
/// cannot solve KW equations for Link %s
pub const ERR_KINWAVE: ErrorCode = 103;
/// cannot open ODE solver
pub const ERR_ODE_SOLVER: ErrorCode = 105;
/// cannot compute a valid time step
pub const ERR_TIMESTEP: ErrorCode = 107;

// --- Object validation errors (108–145) ---
/// ambiguous outlet ID name for Subcatchment %s
pub const ERR_SUBCATCH_OUTLET: ErrorCode = 108;
/// invalid parameter values for Aquifer %s
pub const ERR_AQUIFER_PARAMS: ErrorCode = 109;
/// ground elevation is below water table for Subcatchment %s
pub const ERR_GROUND_ELEV: ErrorCode = 110;
/// invalid length for Conduit %s
pub const ERR_LENGTH: ErrorCode = 111;
/// elevation drop exceeds length for Conduit %s
pub const ERR_ELEV_DROP: ErrorCode = 112;
/// invalid roughness for Conduit %s
pub const ERR_ROUGHNESS: ErrorCode = 113;
/// invalid number of barrels for Conduit %s
pub const ERR_BARRELS: ErrorCode = 114;
/// adverse slope for Conduit %s
pub const ERR_SLOPE: ErrorCode = 115;
/// no cross section defined for Link %s
pub const ERR_NO_XSECT: ErrorCode = 117;
/// invalid cross section for Link %s
pub const ERR_XSECT: ErrorCode = 119;
/// missing or invalid pump curve assigned to Pump %s
pub const ERR_NO_CURVE: ErrorCode = 121;
/// startup depth not higher than shutoff depth for Pump %s
pub const ERR_PUMP_LIMITS: ErrorCode = 122;
/// links form cyclic loops in the drainage system
pub const ERR_LOOP: ErrorCode = 131;
/// Node %s has more than one outlet link
pub const ERR_MULTI_OUTLET: ErrorCode = 133;
/// Node %s has illegal DUMMY link connections
pub const ERR_DUMMY_LINK: ErrorCode = 134;
/// Divider %s does not have two outlet links
pub const ERR_DIVIDER: ErrorCode = 135;
/// Divider %s has invalid diversion link
pub const ERR_DIVIDER_LINK: ErrorCode = 136;
/// Weir Divider %s has invalid parameters
pub const ERR_WEIR_DIVIDER: ErrorCode = 137;
/// Node %s has initial depth greater than maximum depth
pub const ERR_NODE_DEPTH: ErrorCode = 138;
/// Regulator %s is the outlet of a non-storage node
pub const ERR_REGULATOR: ErrorCode = 139;
/// Storage node %s has negative volume at full depth
pub const ERR_STORAGE_VOLUME: ErrorCode = 140;
/// Outfall %s has more than 1 inlet link or an outlet link
pub const ERR_OUTFALL: ErrorCode = 141;
/// Regulator %s has invalid cross-section shape
pub const ERR_REGULATOR_SHAPE: ErrorCode = 143;
/// Drainage system has no acceptable outlet nodes
pub const ERR_NO_OUTLETS: ErrorCode = 145;

// --- Hydrology errors (151–159) ---
/// a Unit Hydrograph in set %s has invalid time base
pub const ERR_UNITHYD_TIMES: ErrorCode = 151;
/// a Unit Hydrograph in set %s has invalid response ratios
pub const ERR_UNITHYD_RATIOS: ErrorCode = 153;
/// invalid sewer area for RDII at node %s
pub const ERR_RDII_AREA: ErrorCode = 155;
/// ambiguous station ID for Rain Gage %s
pub const ERR_RAIN_FILE_CONFLICT: ErrorCode = 156;
/// inconsistent rainfall format for Rain Gage %s
pub const ERR_RAIN_GAGE_FORMAT: ErrorCode = 157;
/// time series for Rain Gage %s is also used by another object
pub const ERR_RAIN_GAGE_TSERIES: ErrorCode = 158;
/// recording interval greater than time series interval for Rain Gage %s
pub const ERR_RAIN_GAGE_INTERVAL: ErrorCode = 159;

// --- Treatment ---
/// cyclic dependency in treatment functions at node %s
pub const ERR_CYCLIC_TREATMENT: ErrorCode = 161;

// --- Curves / timeseries (171–173) ---
/// Curve %s has invalid or out of sequence data
pub const ERR_CURVE_SEQUENCE: ErrorCode = 171;
/// Time Series %s has its data out of sequence
pub const ERR_TIMESERIES_SEQUENCE: ErrorCode = 173;

// --- Snow / LID (181–188) ---
/// invalid Snow Melt Climatology parameters
pub const ERR_SNOWMELT_PARAMS: ErrorCode = 181;
/// invalid parameters for Snow Pack %s
pub const ERR_SNOWPACK_PARAMS: ErrorCode = 182;
/// no type specified for LID %s
pub const ERR_LID_TYPE: ErrorCode = 183;
/// missing layer for LID %s
pub const ERR_LID_LAYER: ErrorCode = 184;
/// invalid parameter value for LID %s
pub const ERR_LID_PARAMS: ErrorCode = 185;
/// LID area exceeds total area for Subcatchment %s
pub const ERR_LID_AREAS: ErrorCode = 187;
/// LID capture area exceeds total impervious area for Subcatchment %s
pub const ERR_LID_CAPTURE_AREA: ErrorCode = 188;

// --- Date / reporting (191–195) ---
/// simulation start date comes after ending date
pub const ERR_START_DATE: ErrorCode = 191;
/// report start date comes after ending date
pub const ERR_REPORT_DATE: ErrorCode = 193;
/// reporting time step or duration is less than routing time step
pub const ERR_REPORT_STEP: ErrorCode = 195;

// --- Input parsing (200–235) ---
/// one or more errors in input file
pub const ERR_INPUT: ErrorCode = 200;
/// too many characters in input line
pub const ERR_INPUT_LINE: ErrorCode = 201;
/// too few items
pub const ERR_ITEMS: ErrorCode = 203;
/// invalid keyword %s
pub const ERR_KEYWORD: ErrorCode = 205;
/// duplicate ID name %s
pub const ERR_DUP_NAME: ErrorCode = 207;
/// undefined object %s
pub const ERR_NAME: ErrorCode = 209;
/// invalid number %s
pub const ERR_NUMBER: ErrorCode = 211;
/// invalid date/time %s
pub const ERR_DATETIME: ErrorCode = 213;
/// control rule clause invalid or out of sequence
pub const ERR_CONTROL_RULE: ErrorCode = 217;
/// data provided for unidentified transect
pub const ERR_TRANSECT_UNKNOWN: ErrorCode = 219;
/// transect station out of sequence
pub const ERR_TRANSECT_SEQUENCE: ErrorCode = 221;
/// Transect %s has too few stations
pub const ERR_TRANSECT_TOO_FEW: ErrorCode = 223;
/// Transect %s has too many stations
pub const ERR_TRANSECT_TOO_MANY: ErrorCode = 225;
/// Transect %s has no Manning's N
pub const ERR_TRANSECT_MANNING: ErrorCode = 227;
/// Transect %s has invalid overbank locations
pub const ERR_TRANSECT_OVERBANK: ErrorCode = 229;
/// Transect %s has no depth
pub const ERR_TRANSECT_NO_DEPTH: ErrorCode = 231;
/// invalid math expression
pub const ERR_MATH_EXPR: ErrorCode = 233;
/// invalid infiltration parameters
pub const ERR_INFILTRATION: ErrorCode = 235;

// --- File I/O (301–363) ---
/// files share same names
pub const ERR_FILE_NAME: ErrorCode = 301;
/// cannot open input file
pub const ERR_INP_FILE: ErrorCode = 303;
/// cannot open report file
pub const ERR_RPT_FILE: ErrorCode = 305;
/// cannot open binary results file
pub const ERR_OUT_FILE: ErrorCode = 307;
/// output will exceed maximum file size
pub const ERR_OUT_SIZE: ErrorCode = 308;
/// error writing to binary results file
pub const ERR_OUT_WRITE: ErrorCode = 309;
/// error reading from binary results file
pub const ERR_OUT_READ: ErrorCode = 311;
/// cannot open scratch rainfall interface file
pub const ERR_RAIN_IFACE_SCRATCH: ErrorCode = 313;
/// cannot open rainfall interface file %s
pub const ERR_RAIN_IFACE: ErrorCode = 315;
/// cannot open rainfall data file %s
pub const ERR_RAIN_FILE_OPEN: ErrorCode = 317;
/// line is out of sequence in rainfall data file %s
pub const ERR_RAIN_FILE_SEQUENCE: ErrorCode = 318;
/// unknown format for rainfall data file %s
pub const ERR_RAIN_FILE_FORMAT: ErrorCode = 319;
/// invalid format for rainfall interface file
pub const ERR_RAIN_IFACE_FORMAT: ErrorCode = 320;
/// no data in rainfall interface file for gage %s
pub const ERR_RAIN_IFACE_GAGE: ErrorCode = 321;
/// cannot open runoff interface file %s
pub const ERR_RUNOFF_IFACE: ErrorCode = 323;
/// incompatible data found in runoff interface file
pub const ERR_RUNOFF_IFACE_COMPAT: ErrorCode = 325;
/// attempting to read beyond end of runoff interface file
pub const ERR_RUNOFF_IFACE_EOF: ErrorCode = 327;
/// error in reading from runoff interface file
pub const ERR_RUNOFF_IFACE_READ: ErrorCode = 329;
/// cannot open hot start interface file %s
pub const ERR_HOTSTART_FILE: ErrorCode = 331;
/// incompatible data found in hot start interface file
pub const ERR_HOTSTART_COMPAT: ErrorCode = 333;
/// error in reading from hot start interface file
pub const ERR_HOTSTART_READ: ErrorCode = 335;
/// no climate file specified for evaporation and/or wind speed
pub const ERR_NO_CLIMATE_FILE: ErrorCode = 336;
/// cannot open climate file %s
pub const ERR_CLIMATE_FILE_OPEN: ErrorCode = 337;
/// error in reading from climate file %s
pub const ERR_CLIMATE_FILE_READ: ErrorCode = 338;
/// attempt to read beyond end of climate file %s
pub const ERR_CLIMATE_FILE_EOF: ErrorCode = 339;
/// cannot open scratch RDII interface file
pub const ERR_RDII_IFACE_SCRATCH: ErrorCode = 341;
/// cannot open RDII interface file %s
pub const ERR_RDII_IFACE: ErrorCode = 343;
/// invalid format for RDII interface file
pub const ERR_RDII_IFACE_FORMAT: ErrorCode = 345;
/// cannot open routing interface file %s
pub const ERR_ROUTING_IFACE: ErrorCode = 351;
/// invalid format for routing interface file %s
pub const ERR_ROUTING_IFACE_FORMAT: ErrorCode = 353;
/// mis-matched names in routing interface file %s
pub const ERR_ROUTING_IFACE_NAMES: ErrorCode = 355;
/// inflows and outflows interface files have same name
pub const ERR_ROUTING_IFACE_SAME: ErrorCode = 357;
/// could not open external file used for Time Series %s
pub const ERR_TABLE_FILE_OPEN: ErrorCode = 361;
/// invalid data in external file used for Time Series %s
pub const ERR_TABLE_FILE_READ: ErrorCode = 363;

// --- System / API errors (500–509) — legacy API error keys ---
/// System exception thrown
pub const ERR_SYSTEM: ErrorCode = 500;
/// project not opened
pub const ERR_API_NOT_OPENED: ErrorCode = 501;
/// simulation not started
pub const ERR_API_NOT_STARTED: ErrorCode = 502;
/// simulation not ended
pub const ERR_API_NOT_ENDED: ErrorCode = 503;
/// invalid object type
pub const ERR_API_OBJECT_TYPE: ErrorCode = 504;
/// invalid object index
pub const ERR_API_OBJECT_INDEX: ErrorCode = 505;
/// invalid object name
pub const ERR_API_OBJECT_NAME: ErrorCode = 506;
/// invalid property type
pub const ERR_API_PROPERTY_TYPE: ErrorCode = 507;
/// invalid property value
pub const ERR_API_PROPERTY_VALUE: ErrorCode = 508;
/// invalid time period
pub const ERR_API_TIME_PERIOD: ErrorCode = 509;

// --- New 6.0 errors (600+) ---
/// Time Series %s has no data
pub const ERR_TIMESERIES_EMPTY: ErrorCode = 601;
/// Time Series %s contains NaN or Inf values
pub const ERR_TIMESERIES_NAN: ErrorCode = 603;
/// column count mismatch in data for %s
pub const ERR_TABLE_COL_MISMATCH: ErrorCode = 605;
/// Rain Gage %s references unknown time series
pub const ERR_GAGE_TSERIES_NOTFOUND: ErrorCode = 607;

// --- Virtual junctions (609–621) — refactored engine only ---
/// Virtual Junction %s must connect exactly two conduits
pub const ERR_VJ_LINK_COUNT: ErrorCode = 609;
/// Virtual Junction %s connects conduits with different cross sections
pub const ERR_VJ_XSECT_MISMATCH: ErrorCode = 611;
/// Virtual Junction %s has a conduit with a nonzero offset
pub const ERR_VJ_OFFSET: ErrorCode = 613;
/// Virtual Junction %s conduit inverts do not agree at the node
pub const ERR_VJ_INVERT_MISMATCH: ErrorCode = 615;
/// Virtual Junction %s cannot be coupled to a 2D surface mesh (point laterals are allowed)
pub const ERR_VJ_LATERAL_INFLOW: ErrorCode = 617;
/// Virtual Junction %s requires DYNWAVE or FV flow routing
pub const ERR_VJ_ROUTING_MODEL: ErrorCode = 619;
/// too many items for Virtual Junction %s
pub const ERR_VJ_EXTRA_TOKENS: ErrorCode = 621;

// --- Inlet junctions (623–635) — refactored engine only ---
/// Inlet Junction %s must connect two STREET conduits (or RECT_OPEN/TRAPEZOIDAL for drop inlets)
pub const ERR_IJ_NOT_STREET: ErrorCode = 623;
/// Inlet Junction %s references an unknown inlet design
pub const ERR_IJ_DESIGN: ErrorCode = 625;
/// Inlet Junction %s has an invalid capture node (missing, itself, or virtual)
pub const ERR_IJ_CAPTURE_NODE: ErrorCode = 627;
/// Inlet Junction %s conduit also carries an [INLET_USAGE] entry
pub const ERR_IJ_USAGE_ON_PAIR: ErrorCode = 629;
/// too many items for Inlet Junction %s
pub const ERR_IJ_EXTRA_TOKENS: ErrorCode = 631;
/// Inlet Junction %s has no inlet design assigned
pub const ERR_IJ_NO_USAGE: ErrorCode = 633;
/// Inlet %s cannot be used with the cross section of %s
pub const ERR_INLET_USAGE_SHAPE: ErrorCode = 635;

// ============================================================================
// CFFI error codes — categorical codes from the C API
// ============================================================================
// These map to the SWMM_ErrorCode enum in openswmm_engine.h and provide
// broad categories. Use cffi_error_get_template() to get descriptions.

/// Upstream `enum CffiErrorCode : int`.
pub type CffiErrorCode = i32;

pub const CFFI_OK: CffiErrorCode = 0;
/// memory allocation error
pub const CFFI_ERR_NOMEM: CffiErrorCode = 1;
/// cannot open input file
pub const CFFI_ERR_INPFILE: CffiErrorCode = 2;
/// cannot open report file
pub const CFFI_ERR_RPTFILE: CffiErrorCode = 3;
/// cannot open output file
pub const CFFI_ERR_OUTFILE: CffiErrorCode = 4;
/// input parsing error
pub const CFFI_ERR_PARSE: CffiErrorCode = 5;
/// invalid engine lifecycle state transition
pub const CFFI_ERR_LIFECYCLE: CffiErrorCode = 6;
/// invalid engine handle
pub const CFFI_ERR_BADHANDLE: CffiErrorCode = 7;
/// invalid object index
pub const CFFI_ERR_BADINDEX: CffiErrorCode = 8;
/// invalid parameter value
pub const CFFI_ERR_BADPARAM: CffiErrorCode = 9;
/// plugin initialization or execution error
pub const CFFI_ERR_PLUGIN: CffiErrorCode = 10;
/// file I/O error
pub const CFFI_ERR_IO: CffiErrorCode = 11;
/// hot start file error
pub const CFFI_ERR_HOTSTART: CffiErrorCode = 12;
/// coordinate reference system error
pub const CFFI_ERR_CRS: CffiErrorCode = 13;
/// numerical solver error
pub const CFFI_ERR_NUMERICAL: CffiErrorCode = 14;
/// internal engine error
pub const CFFI_ERR_INTERNAL: CffiErrorCode = 99;

// ============================================================================
// CFFI warning codes
// ============================================================================

/// Upstream `enum CffiWarnCode : int`.
pub type CffiWarnCode = i32;

pub const CFFI_WARN_NONE: CffiWarnCode = 0;
/// hot start file not found, cold start used
pub const CFFI_WARN_HOTSTART_MISSING: CffiWarnCode = 1;
/// unknown section in input file
pub const CFFI_WARN_UNKNOWN_SECTION: CffiWarnCode = 2;
/// unknown option in input file
pub const CFFI_WARN_UNKNOWN_OPTION: CffiWarnCode = 3;
/// deprecated keyword used
pub const CFFI_WARN_DEPRECATED_KW: CffiWarnCode = 4;
/// plugin initialization warning
pub const CFFI_WARN_PLUGIN_INIT: CffiWarnCode = 5;
/// numerical approximation used
pub const CFFI_WARN_NUMERICAL: CffiWarnCode = 6;
/// time step reduced for stability
pub const CFFI_WARN_STABILITY_LIMIT: CffiWarnCode = 7;

// ============================================================================
// Warning codes — matching legacy text.h WARN01–WARN13
// ============================================================================

/// Upstream `enum WarnCode : int`.
pub type WarnCode = i32;

pub const WARN_NONE: WarnCode = 0;
/// wet weather time step reduced to recording interval for Rain Gage %s
pub const WARN_WET_STEP_REDUCED: WarnCode = 1;
/// maximum depth increased for Node %s
pub const WARN_MAX_DEPTH_INCREASED: WarnCode = 2;
/// negative offset ignored for Link %s
pub const WARN_NEGATIVE_OFFSET: WarnCode = 3;
/// minimum elevation drop used for Conduit %s
pub const WARN_MIN_ELEV_DROP: WarnCode = 4;
/// minimum slope used for Conduit %s
pub const WARN_MIN_SLOPE: WarnCode = 5;
/// dry weather time step increased to the wet weather time step
pub const WARN_DRY_STEP_INCREASED: WarnCode = 6;
/// routing time step reduced to the wet weather time step
pub const WARN_ROUTING_STEP_REDUCED: WarnCode = 7;
/// elevation drop exceeds length for Conduit %s
pub const WARN_ELEV_DROP_EXCEEDS: WarnCode = 8;
/// time series interval greater than recording interval for Rain Gage %s
pub const WARN_GAGE_INTERVAL: WarnCode = 9;
/// crest elevation is below downstream invert for regulator Link %s
pub const WARN_REGULATOR_CREST_LOW: WarnCode = 10;
/// non-matching attributes in Control Rule %s
pub const WARN_CONTROL_RULE_ATTR: WarnCode = 11;
/// inlet removed due to unsupported shape for Conduit %s
pub const WARN_INLET_REMOVED: WarnCode = 12;
/// parameters for same month provided more than once for Unit Hydrograph %s
pub const WARN_UH_PARAMS_REPEATED: WarnCode = 13;

// --- New 6.0 warnings (100+) ---
/// Time Series %s has duplicate x values
pub const WARN_TIMESERIES_DUPLICATE_X: WarnCode = 101;
/// boundary regions overlap for %s
pub const WARN_BOUNDARY_OVERLAP: WarnCode = 102;
/// [FILES] %s is not supported and was ignored
pub const WARN_FILES_SLOT_UNSUPPORTED: WarnCode = 103;
/// [2D_OPTIONS] %s retired with CVODE/ARKODE and was ignored
pub const WARN_2D_OPTION_RETIRED: WarnCode = 104;
/// %s has no effect under FLOW_ROUTING FV
pub const WARN_FV_OPTION_INERT: WarnCode = 105;
/// %s is a dynamic wave option and does not apply under FLOW_ROUTING FV
pub const WARN_DW_OPTION_UNDER_FV: WarnCode = 106;
/// unreadable row(s) skipped in rainfall CSV for Rain Gage %s
pub const WARN_RAIN_CSV_ROWS_SKIPPED: WarnCode = 107;
/// slot width cap (5% of top width) overrides FV_SLOT_CELERITY for %s
pub const WARN_FV_SLOT_CAP: WarnCode = 108;

// ============================================================================
// Error description table — from legacy error.txt
// ============================================================================

/// The legacy/new error description map (upstream `error_table()`).
static ERROR_TABLE: &[(i32, &[u8])] = &[
    // System
    (101, b"memory allocation error."),
    (103, b"cannot solve KW equations for Link %s."),
    (105, b"cannot open ODE solver."),
    (107, b"cannot compute a valid time step."),

    // Object validation
    (108, b"ambiguous outlet ID name for Subcatchment %s."),
    (109, b"invalid parameter values for Aquifer %s."),
    (110, b"ground elevation is below water table for Subcatchment %s."),
    (111, b"invalid length for Conduit %s."),
    (112, b"elevation drop exceeds length for Conduit %s."),
    (113, b"invalid roughness for Conduit %s."),
    (114, b"invalid number of barrels for Conduit %s."),
    (115, b"adverse slope for Conduit %s."),
    (117, b"no cross section defined for Link %s."),
    (119, b"invalid cross section for Link %s."),
    (121, b"missing or invalid pump curve assigned to Pump %s."),
    (122, b"startup depth not higher than shutoff depth for Pump %s."),
    (131, b"the following links form cyclic loops in the drainage system:"),
    (133, b"Node %s has more than one outlet link."),
    (134, b"Node %s has illegal DUMMY link connections."),
    (135, b"Divider %s does not have two outlet links."),
    (136, b"Divider %s has invalid diversion link."),
    (137, b"Weir Divider %s has invalid parameters."),
    (138, b"Node %s has initial depth greater than maximum depth."),
    (139, b"Regulator %s is the outlet of a non-storage node."),
    (140, b"Storage node %s has negative volume at full depth."),
    (141, b"Outfall %s has more than 1 inlet link or an outlet link."),
    (143, b"Regulator %s has invalid cross-section shape."),
    (145, b"Drainage system has no acceptable outlet nodes."),

    // Hydrology
    (151, b"a Unit Hydrograph in set %s has invalid time base."),
    (153, b"a Unit Hydrograph in set %s has invalid response ratios."),
    (155, b"invalid sewer area for RDII at node %s."),
    (156, b"ambiguous station ID for Rain Gage %s."),
    (157, b"inconsistent rainfall format for Rain Gage %s."),
    (158, b"time series for Rain Gage %s is also used by another object."),
    (159, b"recording interval greater than time series interval for Rain Gage %s."),

    // Treatment
    (161, b"cyclic dependency in treatment functions at node %s."),

    // Curves / timeseries
    (171, b"Curve %s has invalid or out of sequence data."),
    (173, b"Time Series %s has its data out of sequence."),

    // Snow / LID
    (181, b"invalid Snow Melt Climatology parameters."),
    (182, b"invalid parameters for Snow Pack %s."),
    (183, b"no type specified for LID %s."),
    (184, b"missing layer for LID %s."),
    (185, b"invalid parameter value for LID %s."),
    (187, b"LID area exceeds total area for Subcatchment %s."),
    (188, b"LID capture area exceeds total impervious area for Subcatchment %s."),

    // Date / reporting
    (191, b"simulation start date comes after ending date."),
    (193, b"report start date comes after ending date."),
    (195, b"reporting time step or duration is less than routing time step."),

    // Input parsing
    (200, b"one or more errors in input file."),
    (201, b"too many characters in input line."),
    (203, b"too few items."),
    (205, b"invalid keyword %s."),
    (207, b"duplicate ID name %s."),
    (209, b"undefined object %s."),
    (211, b"invalid number %s."),
    (213, b"invalid date/time %s."),
    (217, b"control rule clause invalid or out of sequence."),
    (219, b"data provided for unidentified transect."),
    (221, b"transect station out of sequence."),
    (223, b"Transect %s has too few stations."),
    (225, b"Transect %s has too many stations."),
    (227, b"Transect %s has no Manning's N."),
    (229, b"Transect %s has invalid overbank locations."),
    (231, b"Transect %s has no depth."),
    (233, b"invalid math expression."),
    (235, b"invalid infiltration parameters."),

    // File I/O
    (301, b"files share same names."),
    (303, b"cannot open input file."),
    (305, b"cannot open report file."),
    (307, b"cannot open binary results file."),
    (308, b"amount of output produced will exceed maximum file size."),
    (309, b"error writing to binary results file."),
    (311, b"error reading from binary results file."),
    (313, b"cannot open scratch rainfall interface file."),
    (315, b"cannot open rainfall interface file %s."),
    (317, b"cannot open rainfall data file %s."),
    (318, b"the following line is out of sequence in rainfall data file %s."),
    (319, b"unknown format for rainfall data file %s."),
    (320, b"invalid format for rainfall interface file."),
    (321, b"no data in rainfall interface file for gage %s."),
    (323, b"cannot open runoff interface file %s."),
    (325, b"incompatible data found in runoff interface file."),
    (327, b"attempting to read beyond end of runoff interface file."),
    (329, b"error in reading from runoff interface file."),
    (331, b"cannot open hot start interface file %s."),
    (333, b"incompatible data found in hot start interface file."),
    (335, b"error in reading from hot start interface file."),
    (336, b"no climate file specified for evaporation and/or wind speed."),
    (337, b"cannot open climate file %s."),
    (338, b"error in reading from climate file %s."),
    (339, b"attempt to read beyond end of climate file %s."),
    (341, b"cannot open scratch RDII interface file."),
    (343, b"cannot open RDII interface file %s."),
    (345, b"invalid format for RDII interface file."),
    (351, b"cannot open routing interface file %s."),
    (353, b"invalid format for routing interface file %s."),
    (355, b"mis-matched names in routing interface file %s."),
    (357, b"inflows and outflows interface files have same name."),
    (361, b"could not open external file used for Time Series %s."),
    (363, b"invalid data in external file used for Time Series %s."),

    // System / API (500–509)
    (500, b"System exception thrown."),
    (501, b"project not opened."),
    (502, b"simulation not started."),
    (503, b"simulation not ended."),
    (504, b"invalid object type."),
    (505, b"invalid object index."),
    (506, b"invalid object name."),
    (507, b"invalid property type."),
    (508, b"invalid property value."),
    (509, b"invalid time period."),

    // New 6.0
    (601, b"Time Series %s has no data."),
    (603, b"Time Series %s contains NaN or Inf values."),
    (605, b"column count mismatch in data for %s."),
    (607, b"Rain Gage %s references unknown time series."),
    (609, b"Virtual Junction %s must connect exactly two conduits."),
    (611, b"Virtual Junction %s connects conduits with different cross sections."),
    (613, b"Virtual Junction %s has a conduit with a nonzero offset."),
    (615, b"Virtual Junction %s conduit inverts do not agree at the node."),
    (617, b"Virtual Junction %s cannot be coupled to a 2D surface mesh."),
    (619, b"Virtual Junction %s requires DYNWAVE or FV flow routing."),
    (621, b"too many items for Virtual Junction %s."),
    (623, b"Inlet Junction %s must connect two STREET conduits."),
    (625, b"Inlet Junction %s references an unknown inlet design."),
    (627, b"Inlet Junction %s has an invalid capture node."),
    (629, b"Inlet Junction %s has a conduit that also carries an [INLET_USAGE] entry."),
    (631, b"too many items for Inlet Junction %s."),
    (633, b"Inlet Junction %s has no inlet design assigned."),
    (635, b"Inlet %s cannot be used with the cross section of its host."),
];

// ============================================================================
// Warning description table — from legacy text.h
// ============================================================================

/// The legacy/new warning description map (upstream `warning_table()`).
static WARNING_TABLE: &[(i32, &[u8])] = &[
    // Codes 1-13 are legacy's own, and legacy's report_writeWarningMsg
    // prints `"%s %s"` — the text from text.h then the object ID, with
    // NO closing period (report.c:writeWarningMsg, text.h:47-64). v6
    // ended each of these with '.', so every deck that raises one
    // differed from the oracle on punctuation alone: `WARNING 04:
    // minimum elevation drop used for Conduit 31.` against legacy's
    // `... Conduit 31`. Kept verbatim against text.h; the 6.0-only
    // warnings below have no legacy counterpart and keep their own style.
    (1, b"wet weather time step reduced to recording interval for Rain Gage %s"),
    (2, b"maximum depth increased for Node %s"),
    (3, b"negative offset ignored for Link %s"),
    (4, b"minimum elevation drop used for Conduit %s"),
    (5, b"minimum slope used for Conduit %s"),
    (6, b"dry weather time step increased to the wet weather time step"),
    (7, b"routing time step reduced to the wet weather time step"),
    (8, b"elevation drop exceeds length for Conduit %s"),
    (9, b"time series interval greater than recording interval for Rain Gage %s"),
    (10, b"crest elevation is below downstream invert for regulator Link %s"),
    (11, b"non-matching attributes in Control Rule %s"),
    (12, b"inlet removed due to unsupported shape for Conduit %s"),
    (13, b"parameters for same month provided more than once for Unit Hydrograph %s"),

    // New 6.0
    (101, b"Time Series %s has duplicate x values."),
    (102, b"boundary regions overlap for %s."),
    (103, b"[FILES] %s is not supported by this engine and was ignored."),
    (104, b"[2D_OPTIONS] %s was retired with the CVODE/ARKODE 2D solvers \
          and was ignored; the explicit local-inertial marcher is the \
          only 2D integrator."),
    (105, b"%s is accepted but has no effect yet under FLOW_ROUTING FV."),
    (106, b"%s is a dynamic wave option and does not apply under \
          FLOW_ROUTING FV."),
    (107, b"unreadable row(s) were skipped in the rainfall CSV for Rain \
          Gage %s."),
    (108, b"the Preissmann slot width cap (5% of the section top width) \
          overrides the requested FV_SLOT_CELERITY for %s; celerities \
          below the cap-implied value have no effect there."),
];

/// `table.find(code)`: the text stored under `code`, if any.
fn table_find(table: &'static [(i32, &'static [u8])], code: i32) -> Option<&'static [u8]> {
    table.iter().find(|e| e.0 == code).map(|e| e.1)
}

// ============================================================================
// Template lookup
// ============================================================================

/// The description template for an error code, with `%s` for the object name; empty for
/// an unknown code.
pub fn error_get_template(code: i32) -> &'static [u8] {
    match table_find(ERROR_TABLE, code) {
        Some(s) => s,
        None => b"",
    }
}

/// The description template for a warning code; empty for an unknown code.
pub fn warning_get_template(code: i32) -> &'static [u8] {
    match table_find(WARNING_TABLE, code) {
        Some(s) => s,
        None => b"",
    }
}

// ============================================================================
// Formatting helpers
// ============================================================================

/// `tmpl.find("%s")`: the position of the first "%s", if any.
fn find_percent_s(tmpl: &[u8]) -> Option<usize> {
    tmpl.windows(2).position(|w| w == b"%s")
}

/// Replace the first literal "%s" in `tmpl` with `name`, copied exactly; a template with
/// no marker comes back unchanged.
#[allow(clippy::reserve_after_initialization)]
fn substitute(tmpl: &[u8], name: &[u8]) -> Bytes {
    let mut result = Bytes::new();
    result.reserve(tmpl.len() + name.len());
    let pos = find_percent_s(tmpl);
    if let Some(pos) = pos {
        result.extend_from_slice(&tmpl[..pos]);
        result.extend_from_slice(name);
        result.extend_from_slice(&tmpl[pos + 2..]);
    } else {
        result.extend_from_slice(tmpl);
    }
    result
}

/// `"  ERROR <code>: "` and the template with `name` in place of `%s`, or
/// `"  ERROR <code>: unknown error."` for an unknown code.
pub fn format_error(code: i32, name: &[u8]) -> Bytes {
    let tmpl = error_get_template(code);
    if tmpl.is_empty() {
        return format!("  ERROR {}: unknown error.", code).into_bytes();
    }
    let mut msg = format!("  ERROR {}: ", code).into_bytes();
    msg.extend_from_slice(&substitute(tmpl, name));
    msg
}

/// `format_error(code, name)` with `detail` appended after a space, before the final
/// period when there is one; an empty `detail` adds nothing.
pub fn format_error_detail(code: i32, name: &[u8], detail: &[u8]) -> Bytes {
    let mut msg = format_error(code, name);
    if !detail.is_empty() {
        // Insert detail before the trailing period (if any)
        if !msg.is_empty() && msg.last() == Some(&b'.') {
            msg.pop();
            msg.push(b' ');
            msg.extend_from_slice(detail);
            msg.push(b'.');
        } else {
            msg.push(b' ');
            msg.extend_from_slice(detail);
        }
    }
    msg
}

/// `"  WARNING <code, two digits>: "` and the template with `name` in place of `%s`, or
/// `"  WARNING <code>: unknown warning."` for an unknown code.
pub fn format_warning(code: i32, name: &[u8]) -> Bytes {
    let tmpl = warning_get_template(code);
    if tmpl.is_empty() {
        return format!("  WARNING {:02}: unknown warning.", code).into_bytes();
    }
    let mut msg = format!("  WARNING {:02}: ", code).into_bytes();
    msg.extend_from_slice(&substitute(tmpl, name));
    msg
}

/// `format_warning(code, name)` with `detail` appended after a space, before the final
/// period when there is one; an empty `detail` adds nothing.
pub fn format_warning_detail(code: i32, name: &[u8], detail: &[u8]) -> Bytes {
    let mut msg = format_warning(code, name);
    if !detail.is_empty() {
        if !msg.is_empty() && msg.last() == Some(&b'.') {
            msg.pop();
            msg.push(b' ');
            msg.extend_from_slice(detail);
            msg.push(b'.');
        } else {
            msg.push(b' ');
            msg.extend_from_slice(detail);
        }
    }
    msg
}

// ============================================================================
// CFFI / C API description tables
// ============================================================================

/// The C API error description map (upstream `cffi_error_table()`).
static CFFI_ERROR_TABLE: &[(i32, &[u8])] = &[
    (0, b"no error."),
    (1, b"memory allocation error."),
    (2, b"cannot open input file."),
    (3, b"cannot open report file."),
    (4, b"cannot open output file."),
    (5, b"input parsing error."),
    (6, b"invalid engine lifecycle state transition."),
    (7, b"invalid engine handle."),
    (8, b"invalid object index."),
    (9, b"invalid parameter value."),
    (10, b"plugin initialization or execution error."),
    (11, b"file I/O error."),
    (12, b"hot start file error."),
    (13, b"coordinate reference system error."),
    (14, b"numerical solver error."),
    (99, b"internal engine error."),
];

/// The C API warning description map (upstream `cffi_warning_table()`).
static CFFI_WARNING_TABLE: &[(i32, &[u8])] = &[
    (0, b"no warning."),
    (1, b"hot start file not found, cold start used."),
    (2, b"unknown section in input file."),
    (3, b"unknown option in input file."),
    (4, b"deprecated keyword used."),
    (5, b"plugin initialization warning."),
    (6, b"numerical approximation used."),
    (7, b"time step reduced for stability."),
];

/// The description for a CFFI error code (SWMM_ErrorCode); empty for an unknown code.
pub fn cffi_error_get_template(code: i32) -> &'static [u8] {
    match table_find(CFFI_ERROR_TABLE, code) {
        Some(s) => s,
        None => b"",
    }
}

/// The description for a CFFI warning code (SWMM_WarnCode); empty for an unknown code.
pub fn cffi_warning_get_template(code: i32) -> &'static [u8] {
    match table_find(CFFI_WARNING_TABLE, code) {
        Some(s) => s,
        None => b"",
    }
}

/// `"  ERROR <code>: "` and the CFFI description with `name` in place of any `%s`, or
/// `"  ERROR <code>: unknown CFFI error."` for an unknown code.
pub fn format_cffi_error(code: i32, name: &[u8]) -> Bytes {
    let tmpl = cffi_error_get_template(code);
    if tmpl.is_empty() {
        return format!("  ERROR {}: unknown CFFI error.", code).into_bytes();
    }
    let mut msg = format!("  ERROR {}: ", code).into_bytes();
    msg.extend_from_slice(&substitute(tmpl, name));
    msg
}

/// `"  WARNING <code>: "` and the CFFI description with `name` in place of any `%s`, or
/// `"  WARNING <code>: unknown CFFI warning."` for an unknown code.
pub fn format_cffi_warning(code: i32, name: &[u8]) -> Bytes {
    let tmpl = cffi_warning_get_template(code);
    if tmpl.is_empty() {
        return format!("  WARNING {}: unknown CFFI warning.", code).into_bytes();
    }
    let mut msg = format!("  WARNING {}: ", code).into_bytes();
    msg.extend_from_slice(&substitute(tmpl, name));
    msg
}
