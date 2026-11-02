// Rust side of gate/diff/errdump.cpp: the same lines, byte for byte.
use std::io::Write;
use swmm6_rust::core::error_codes::*;

fn line(out: &mut Vec<u8>, tag: &str, code: i32, i: i32, j: i32, s: &[u8]) {
    out.extend_from_slice(format!("{} {} {} {} |", tag, code, i, j).as_bytes());
    out.extend_from_slice(s);
    out.extend_from_slice(b"|\n");
}

fn main() {
    let mut out: Vec<u8> = Vec::new();
    macro_rules! k {
        ($x:ident) => {
            out.extend_from_slice(format!("{} {}\n", stringify!($x), $x).as_bytes())
        };
    }
    k!(ERR_NONE); k!(ERR_MEMORY); k!(ERR_KINWAVE); k!(ERR_ODE_SOLVER);
    k!(ERR_TIMESTEP); k!(ERR_SUBCATCH_OUTLET); k!(ERR_AQUIFER_PARAMS); k!(ERR_GROUND_ELEV);
    k!(ERR_LENGTH); k!(ERR_ELEV_DROP); k!(ERR_ROUGHNESS); k!(ERR_BARRELS);
    k!(ERR_SLOPE); k!(ERR_NO_XSECT); k!(ERR_XSECT); k!(ERR_NO_CURVE);
    k!(ERR_PUMP_LIMITS); k!(ERR_LOOP); k!(ERR_MULTI_OUTLET); k!(ERR_DUMMY_LINK);
    k!(ERR_DIVIDER); k!(ERR_DIVIDER_LINK); k!(ERR_WEIR_DIVIDER); k!(ERR_NODE_DEPTH);
    k!(ERR_REGULATOR); k!(ERR_STORAGE_VOLUME); k!(ERR_OUTFALL); k!(ERR_REGULATOR_SHAPE);
    k!(ERR_NO_OUTLETS); k!(ERR_UNITHYD_TIMES); k!(ERR_UNITHYD_RATIOS); k!(ERR_RDII_AREA);
    k!(ERR_RAIN_FILE_CONFLICT); k!(ERR_RAIN_GAGE_FORMAT); k!(ERR_RAIN_GAGE_TSERIES); k!(ERR_RAIN_GAGE_INTERVAL);
    k!(ERR_CYCLIC_TREATMENT); k!(ERR_CURVE_SEQUENCE); k!(ERR_TIMESERIES_SEQUENCE); k!(ERR_SNOWMELT_PARAMS);
    k!(ERR_SNOWPACK_PARAMS); k!(ERR_LID_TYPE); k!(ERR_LID_LAYER); k!(ERR_LID_PARAMS);
    k!(ERR_LID_AREAS); k!(ERR_LID_CAPTURE_AREA); k!(ERR_START_DATE); k!(ERR_REPORT_DATE);
    k!(ERR_REPORT_STEP); k!(ERR_INPUT); k!(ERR_INPUT_LINE); k!(ERR_ITEMS);
    k!(ERR_KEYWORD); k!(ERR_DUP_NAME); k!(ERR_NAME); k!(ERR_NUMBER);
    k!(ERR_DATETIME); k!(ERR_CONTROL_RULE); k!(ERR_TRANSECT_UNKNOWN); k!(ERR_TRANSECT_SEQUENCE);
    k!(ERR_TRANSECT_TOO_FEW); k!(ERR_TRANSECT_TOO_MANY); k!(ERR_TRANSECT_MANNING); k!(ERR_TRANSECT_OVERBANK);
    k!(ERR_TRANSECT_NO_DEPTH); k!(ERR_MATH_EXPR); k!(ERR_INFILTRATION); k!(ERR_FILE_NAME);
    k!(ERR_INP_FILE); k!(ERR_RPT_FILE); k!(ERR_OUT_FILE); k!(ERR_OUT_SIZE);
    k!(ERR_OUT_WRITE); k!(ERR_OUT_READ); k!(ERR_RAIN_IFACE_SCRATCH); k!(ERR_RAIN_IFACE);
    k!(ERR_RAIN_FILE_OPEN); k!(ERR_RAIN_FILE_SEQUENCE); k!(ERR_RAIN_FILE_FORMAT); k!(ERR_RAIN_IFACE_FORMAT);
    k!(ERR_RAIN_IFACE_GAGE); k!(ERR_RUNOFF_IFACE); k!(ERR_RUNOFF_IFACE_COMPAT); k!(ERR_RUNOFF_IFACE_EOF);
    k!(ERR_RUNOFF_IFACE_READ); k!(ERR_HOTSTART_FILE); k!(ERR_HOTSTART_COMPAT); k!(ERR_HOTSTART_READ);
    k!(ERR_NO_CLIMATE_FILE); k!(ERR_CLIMATE_FILE_OPEN); k!(ERR_CLIMATE_FILE_READ); k!(ERR_CLIMATE_FILE_EOF);
    k!(ERR_RDII_IFACE_SCRATCH); k!(ERR_RDII_IFACE); k!(ERR_RDII_IFACE_FORMAT); k!(ERR_ROUTING_IFACE);
    k!(ERR_ROUTING_IFACE_FORMAT); k!(ERR_ROUTING_IFACE_NAMES); k!(ERR_ROUTING_IFACE_SAME); k!(ERR_TABLE_FILE_OPEN);
    k!(ERR_TABLE_FILE_READ); k!(ERR_SYSTEM); k!(ERR_API_NOT_OPENED); k!(ERR_API_NOT_STARTED);
    k!(ERR_API_NOT_ENDED); k!(ERR_API_OBJECT_TYPE); k!(ERR_API_OBJECT_INDEX); k!(ERR_API_OBJECT_NAME);
    k!(ERR_API_PROPERTY_TYPE); k!(ERR_API_PROPERTY_VALUE); k!(ERR_API_TIME_PERIOD); k!(ERR_TIMESERIES_EMPTY);
    k!(ERR_TIMESERIES_NAN); k!(ERR_TABLE_COL_MISMATCH); k!(ERR_GAGE_TSERIES_NOTFOUND); k!(ERR_VJ_LINK_COUNT);
    k!(ERR_VJ_XSECT_MISMATCH); k!(ERR_VJ_OFFSET); k!(ERR_VJ_INVERT_MISMATCH); k!(ERR_VJ_LATERAL_INFLOW);
    k!(ERR_VJ_ROUTING_MODEL); k!(ERR_VJ_EXTRA_TOKENS); k!(ERR_IJ_NOT_STREET); k!(ERR_IJ_DESIGN);
    k!(ERR_IJ_CAPTURE_NODE); k!(ERR_IJ_USAGE_ON_PAIR); k!(ERR_IJ_EXTRA_TOKENS); k!(ERR_IJ_NO_USAGE);
    k!(ERR_INLET_USAGE_SHAPE); k!(CFFI_OK); k!(CFFI_ERR_NOMEM); k!(CFFI_ERR_INPFILE);
    k!(CFFI_ERR_RPTFILE); k!(CFFI_ERR_OUTFILE); k!(CFFI_ERR_PARSE); k!(CFFI_ERR_LIFECYCLE);
    k!(CFFI_ERR_BADHANDLE); k!(CFFI_ERR_BADINDEX); k!(CFFI_ERR_BADPARAM); k!(CFFI_ERR_PLUGIN);
    k!(CFFI_ERR_IO); k!(CFFI_ERR_HOTSTART); k!(CFFI_ERR_CRS); k!(CFFI_ERR_NUMERICAL);
    k!(CFFI_ERR_INTERNAL); k!(CFFI_WARN_NONE); k!(CFFI_WARN_HOTSTART_MISSING); k!(CFFI_WARN_UNKNOWN_SECTION);
    k!(CFFI_WARN_UNKNOWN_OPTION); k!(CFFI_WARN_DEPRECATED_KW); k!(CFFI_WARN_PLUGIN_INIT); k!(CFFI_WARN_NUMERICAL);
    k!(CFFI_WARN_STABILITY_LIMIT); k!(WARN_NONE); k!(WARN_WET_STEP_REDUCED); k!(WARN_MAX_DEPTH_INCREASED);
    k!(WARN_NEGATIVE_OFFSET); k!(WARN_MIN_ELEV_DROP); k!(WARN_MIN_SLOPE); k!(WARN_DRY_STEP_INCREASED);
    k!(WARN_ROUTING_STEP_REDUCED); k!(WARN_ELEV_DROP_EXCEEDS); k!(WARN_GAGE_INTERVAL); k!(WARN_REGULATOR_CREST_LOW);
    k!(WARN_CONTROL_RULE_ATTR); k!(WARN_INLET_REMOVED); k!(WARN_UH_PARAMS_REPEATED); k!(WARN_TIMESERIES_DUPLICATE_X);
    k!(WARN_BOUNDARY_OVERLAP); k!(WARN_FILES_SLOT_UNSUPPORTED); k!(WARN_2D_OPTION_RETIRED); k!(WARN_FV_OPTION_INERT);
    k!(WARN_DW_OPTION_UNDER_FV); k!(WARN_RAIN_CSV_ROWS_SKIPPED); k!(WARN_FV_SLOT_CAP);

    let mut codes: Vec<i32> = (-3..=640).collect();
    codes.extend_from_slice(&[700, 999, 1000, 12345, -100, i32::MIN, i32::MAX]);
    let longname = vec![b'x'; 200];
    let names: Vec<&[u8]> = vec![
        b"", b"RAIN1", b"Conduit 31", b"%s", b"a%sb%s", b"N\0UL", b"\xff\xfe\xe9", b".",
        &longname,
    ];
    let details: Vec<&[u8]> = vec![b"", b"at 01/01/2024 02:00:00", b".", b"x.", b"d\0e"];
    for &c in &codes {
        out.extend_from_slice(format!("T {} |", c).as_bytes());
        out.extend_from_slice(error_get_template(c));
        out.extend_from_slice(b"|");
        out.extend_from_slice(warning_get_template(c));
        out.extend_from_slice(b"|");
        out.extend_from_slice(cffi_error_get_template(c));
        out.extend_from_slice(b"|");
        out.extend_from_slice(cffi_warning_get_template(c));
        out.extend_from_slice(b"|\n");
        for (i, name) in names.iter().enumerate() {
            let i = i as i32;
            line(&mut out, "FE", c, i, -1, &format_error(c, name));
            line(&mut out, "FW", c, i, -1, &format_warning(c, name));
            line(&mut out, "FCE", c, i, -1, &format_cffi_error(c, name));
            line(&mut out, "FCW", c, i, -1, &format_cffi_warning(c, name));
            for (j, detail) in details.iter().enumerate() {
                let j = j as i32;
                line(&mut out, "FED", c, i, j, &format_error_detail(c, name, detail));
                line(&mut out, "FWD", c, i, j, &format_warning_detail(c, name, detail));
            }
        }
        // default arguments
        line(&mut out, "FE0", c, -1, -1, &format_error(c, b""));
        line(&mut out, "FW0", c, -1, -1, &format_warning(c, b""));
        line(&mut out, "FCE0", c, -1, -1, &format_cffi_error(c, b""));
        line(&mut out, "FCW0", c, -1, -1, &format_cffi_warning(c, b""));
    }
    std::io::stdout().write_all(&out).unwrap();
}
