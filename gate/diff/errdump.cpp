// Differential harness for ErrorCodes.hpp / ErrorCodes.cpp, compiled from upstream's own
// ErrorCodes.cpp. Prints every enumerator's value, then, for every code in -3..640 and a
// few outside it, all four templates and every formatter over a set of names and details
// (empty, plain, "%s"-bearing, embedded NUL, non-UTF-8, long). Bytes are written raw, so
// NUL and high bytes compare exactly. Rust side: crates/engine/examples/errdump.rs.
#include "ErrorCodes.hpp"
#include <climits>
#include <cstdio>
#include <string>
#include <string_view>
#include <vector>
using namespace openswmm;
static void put(std::string_view s){ std::fwrite(s.data(),1,s.size(),stdout); }
static void line(const char* tag,int code,int i,int j,std::string_view s){
  std::printf("%s %d %d %d |",tag,code,i,j); put(s); put("|\n");
}
#define K(x) std::printf("%s %d\n", #x, (int)(x))
int main(){
  K(ERR_NONE); K(ERR_MEMORY); K(ERR_KINWAVE); K(ERR_ODE_SOLVER);
  K(ERR_TIMESTEP); K(ERR_SUBCATCH_OUTLET); K(ERR_AQUIFER_PARAMS); K(ERR_GROUND_ELEV);
  K(ERR_LENGTH); K(ERR_ELEV_DROP); K(ERR_ROUGHNESS); K(ERR_BARRELS);
  K(ERR_SLOPE); K(ERR_NO_XSECT); K(ERR_XSECT); K(ERR_NO_CURVE);
  K(ERR_PUMP_LIMITS); K(ERR_LOOP); K(ERR_MULTI_OUTLET); K(ERR_DUMMY_LINK);
  K(ERR_DIVIDER); K(ERR_DIVIDER_LINK); K(ERR_WEIR_DIVIDER); K(ERR_NODE_DEPTH);
  K(ERR_REGULATOR); K(ERR_STORAGE_VOLUME); K(ERR_OUTFALL); K(ERR_REGULATOR_SHAPE);
  K(ERR_NO_OUTLETS); K(ERR_UNITHYD_TIMES); K(ERR_UNITHYD_RATIOS); K(ERR_RDII_AREA);
  K(ERR_RAIN_FILE_CONFLICT); K(ERR_RAIN_GAGE_FORMAT); K(ERR_RAIN_GAGE_TSERIES); K(ERR_RAIN_GAGE_INTERVAL);
  K(ERR_CYCLIC_TREATMENT); K(ERR_CURVE_SEQUENCE); K(ERR_TIMESERIES_SEQUENCE); K(ERR_SNOWMELT_PARAMS);
  K(ERR_SNOWPACK_PARAMS); K(ERR_LID_TYPE); K(ERR_LID_LAYER); K(ERR_LID_PARAMS);
  K(ERR_LID_AREAS); K(ERR_LID_CAPTURE_AREA); K(ERR_START_DATE); K(ERR_REPORT_DATE);
  K(ERR_REPORT_STEP); K(ERR_INPUT); K(ERR_INPUT_LINE); K(ERR_ITEMS);
  K(ERR_KEYWORD); K(ERR_DUP_NAME); K(ERR_NAME); K(ERR_NUMBER);
  K(ERR_DATETIME); K(ERR_CONTROL_RULE); K(ERR_TRANSECT_UNKNOWN); K(ERR_TRANSECT_SEQUENCE);
  K(ERR_TRANSECT_TOO_FEW); K(ERR_TRANSECT_TOO_MANY); K(ERR_TRANSECT_MANNING); K(ERR_TRANSECT_OVERBANK);
  K(ERR_TRANSECT_NO_DEPTH); K(ERR_MATH_EXPR); K(ERR_INFILTRATION); K(ERR_FILE_NAME);
  K(ERR_INP_FILE); K(ERR_RPT_FILE); K(ERR_OUT_FILE); K(ERR_OUT_SIZE);
  K(ERR_OUT_WRITE); K(ERR_OUT_READ); K(ERR_RAIN_IFACE_SCRATCH); K(ERR_RAIN_IFACE);
  K(ERR_RAIN_FILE_OPEN); K(ERR_RAIN_FILE_SEQUENCE); K(ERR_RAIN_FILE_FORMAT); K(ERR_RAIN_IFACE_FORMAT);
  K(ERR_RAIN_IFACE_GAGE); K(ERR_RUNOFF_IFACE); K(ERR_RUNOFF_IFACE_COMPAT); K(ERR_RUNOFF_IFACE_EOF);
  K(ERR_RUNOFF_IFACE_READ); K(ERR_HOTSTART_FILE); K(ERR_HOTSTART_COMPAT); K(ERR_HOTSTART_READ);
  K(ERR_NO_CLIMATE_FILE); K(ERR_CLIMATE_FILE_OPEN); K(ERR_CLIMATE_FILE_READ); K(ERR_CLIMATE_FILE_EOF);
  K(ERR_RDII_IFACE_SCRATCH); K(ERR_RDII_IFACE); K(ERR_RDII_IFACE_FORMAT); K(ERR_ROUTING_IFACE);
  K(ERR_ROUTING_IFACE_FORMAT); K(ERR_ROUTING_IFACE_NAMES); K(ERR_ROUTING_IFACE_SAME); K(ERR_TABLE_FILE_OPEN);
  K(ERR_TABLE_FILE_READ); K(ERR_SYSTEM); K(ERR_API_NOT_OPENED); K(ERR_API_NOT_STARTED);
  K(ERR_API_NOT_ENDED); K(ERR_API_OBJECT_TYPE); K(ERR_API_OBJECT_INDEX); K(ERR_API_OBJECT_NAME);
  K(ERR_API_PROPERTY_TYPE); K(ERR_API_PROPERTY_VALUE); K(ERR_API_TIME_PERIOD); K(ERR_TIMESERIES_EMPTY);
  K(ERR_TIMESERIES_NAN); K(ERR_TABLE_COL_MISMATCH); K(ERR_GAGE_TSERIES_NOTFOUND); K(ERR_VJ_LINK_COUNT);
  K(ERR_VJ_XSECT_MISMATCH); K(ERR_VJ_OFFSET); K(ERR_VJ_INVERT_MISMATCH); K(ERR_VJ_LATERAL_INFLOW);
  K(ERR_VJ_ROUTING_MODEL); K(ERR_VJ_EXTRA_TOKENS); K(ERR_IJ_NOT_STREET); K(ERR_IJ_DESIGN);
  K(ERR_IJ_CAPTURE_NODE); K(ERR_IJ_USAGE_ON_PAIR); K(ERR_IJ_EXTRA_TOKENS); K(ERR_IJ_NO_USAGE);
  K(ERR_INLET_USAGE_SHAPE); K(CFFI_OK); K(CFFI_ERR_NOMEM); K(CFFI_ERR_INPFILE);
  K(CFFI_ERR_RPTFILE); K(CFFI_ERR_OUTFILE); K(CFFI_ERR_PARSE); K(CFFI_ERR_LIFECYCLE);
  K(CFFI_ERR_BADHANDLE); K(CFFI_ERR_BADINDEX); K(CFFI_ERR_BADPARAM); K(CFFI_ERR_PLUGIN);
  K(CFFI_ERR_IO); K(CFFI_ERR_HOTSTART); K(CFFI_ERR_CRS); K(CFFI_ERR_NUMERICAL);
  K(CFFI_ERR_INTERNAL); K(CFFI_WARN_NONE); K(CFFI_WARN_HOTSTART_MISSING); K(CFFI_WARN_UNKNOWN_SECTION);
  K(CFFI_WARN_UNKNOWN_OPTION); K(CFFI_WARN_DEPRECATED_KW); K(CFFI_WARN_PLUGIN_INIT); K(CFFI_WARN_NUMERICAL);
  K(CFFI_WARN_STABILITY_LIMIT); K(WARN_NONE); K(WARN_WET_STEP_REDUCED); K(WARN_MAX_DEPTH_INCREASED);
  K(WARN_NEGATIVE_OFFSET); K(WARN_MIN_ELEV_DROP); K(WARN_MIN_SLOPE); K(WARN_DRY_STEP_INCREASED);
  K(WARN_ROUTING_STEP_REDUCED); K(WARN_ELEV_DROP_EXCEEDS); K(WARN_GAGE_INTERVAL); K(WARN_REGULATOR_CREST_LOW);
  K(WARN_CONTROL_RULE_ATTR); K(WARN_INLET_REMOVED); K(WARN_UH_PARAMS_REPEATED); K(WARN_TIMESERIES_DUPLICATE_X);
  K(WARN_BOUNDARY_OVERLAP); K(WARN_FILES_SLOT_UNSUPPORTED); K(WARN_2D_OPTION_RETIRED); K(WARN_FV_OPTION_INERT);
  K(WARN_DW_OPTION_UNDER_FV); K(WARN_RAIN_CSV_ROWS_SKIPPED); K(WARN_FV_SLOT_CAP);

  std::vector<int> codes;
  for(int c=-3;c<=640;++c) codes.push_back(c);
  for(int c: {700, 999, 1000, 12345, -100, INT_MIN, INT_MAX}) codes.push_back(c);
  const std::string longname(200,'x');
  const std::vector<std::string_view> names = {
    std::string_view(""), std::string_view("RAIN1"), std::string_view("Conduit 31"),
    std::string_view("%s"), std::string_view("a%sb%s"), std::string_view("N\0UL",4),
    std::string_view("\xff\xfe\xe9"), std::string_view("."), std::string_view(longname)};
  const std::vector<std::string_view> details = {
    std::string_view(""), std::string_view("at 01/01/2024 02:00:00"), std::string_view("."),
    std::string_view("x."), std::string_view("d\0e",3)};
  for(int c: codes){
    std::printf("T %d |",c); put(error_get_template(c)); put("|"); put(warning_get_template(c));
    put("|"); put(cffi_error_get_template(c)); put("|"); put(cffi_warning_get_template(c)); put("|\n");
    for(int i=0;i<(int)names.size();++i){
      line("FE",c,i,-1,format_error(c,names[i]));
      line("FW",c,i,-1,format_warning(c,names[i]));
      line("FCE",c,i,-1,format_cffi_error(c,names[i]));
      line("FCW",c,i,-1,format_cffi_warning(c,names[i]));
      for(int j=0;j<(int)details.size();++j){
        line("FED",c,i,j,format_error(c,names[i],details[j]));
        line("FWD",c,i,j,format_warning(c,names[i],details[j]));
      }
    }
    // default arguments
    line("FE0",c,-1,-1,format_error(c));
    line("FW0",c,-1,-1,format_warning(c));
    line("FCE0",c,-1,-1,format_cffi_error(c));
    line("FCW0",c,-1,-1,format_cffi_warning(c));
  }
}
