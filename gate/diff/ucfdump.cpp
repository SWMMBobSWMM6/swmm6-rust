// Differential harness for UnitConversion.hpp / UnitConversion.cpp, compiled from
// upstream's own UnitConversion.cpp. Prints, as hex bits: the Ucf/Qcf tables and their
// constexpr reciprocals, the fixed constants, UCF() and UCF_inv() for every flow-unit
// system and every quantity (plus out-of-range quantities), getUnitSystem() over a range,
// every field of DisplayUnits (default-constructed and from() for each flow unit), and
// DisplayUnits::temperature() over a set of temperatures. Rust side:
// crates/engine/examples/ucfdump.rs.
#include "UnitConversion.hpp"
#include "SimulationOptions.hpp"
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
using namespace openswmm;
static void hx(double v){ uint64_t u; std::memcpy(&u,&v,8); std::printf(" %016llx",(unsigned long long)u); }
static void du(const char* tag,int fu,const ucf::DisplayUnits& d){
  std::printf("%s %d %d %d",tag,fu,d.unit_system,d.flow_units);
  hx(d.flow); hx(d.length); hx(d.volume); hx(d.rainfall); hx(d.raindepth); hx(d.evaprate);
  hx(d.landarea); hx(d.mvol); hx(d.landvol);
  std::printf(" |%s|%s|%s|%s|%s|%s|%s|\n",d.flow_word,d.length_word,d.depth_word,d.vel_word,
              d.mvol_word,d.landvol_word,d.storage_vol_word);
  const double ts[] = {32.0, 212.0, -40.0, 0.0, -0.0, 98.6, 1e308, -1e308, 5e-324, 33.8, 50.0,
                       std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::quiet_NaN()};
  std::printf("TEMP %d",fu);
  for(double t: ts) hx(d.temperature(t));
  std::printf("\n");
}
int main(){
  for(int q=0;q<10;++q){ std::printf("Ucf %d",q); hx(ucf::Ucf[q][0]); hx(ucf::Ucf[q][1]); std::printf("\n"); }
  for(int q=0;q<10;++q){ std::printf("Ucf_inv %d",q); hx(ucf::Ucf_inv[q][0]); hx(ucf::Ucf_inv[q][1]); std::printf("\n"); }
  for(int f=0;f<6;++f){ std::printf("Qcf %d",f); hx(ucf::Qcf[f]); hx(ucf::Qcf_inv[f]); std::printf(" %s\n",ucf::FlowUnitWords[f]); }
  std::printf("CONST"); hx(ucf::ACRES_TO_FT2); hx(ucf::FT3_TO_MGAL); hx(ucf::SEC_PER_DAY); std::printf("\n");
  std::printf("Q %d %d %d %d %d %d %d %d %d %d %d\n",ucf::RAINFALL,ucf::RAINDEPTH,ucf::EVAPRATE,ucf::LENGTH,
              ucf::LANDAREA,ucf::VOLUME,ucf::WINDSPEED,ucf::TEMPERATURE,ucf::MASS,ucf::GWFLOW,ucf::FLOW);
  for(int fu: {-3,-2,-1,0,1,2,3,4,5,6,7,8,INT_MIN,INT_MAX}) std::printf("US %d %d\n",fu,ucf::getUnitSystem(fu));
  for(int fu=0;fu<6;++fu){
    SimulationOptions opts; opts.flow_units = static_cast<FlowUnits>(fu);
    for(int q: {-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,100,INT_MIN,INT_MAX}){
      std::printf("UCF %d %d",fu,q); hx(ucf::UCF(q,opts)); hx(ucf::UCF_inv(q,opts)); std::printf("\n");
    }
    du("DU",fu,ucf::DisplayUnits::from(opts));
  }
  du("DU0",-1,ucf::DisplayUnits{});
  { SimulationOptions opts; du("DUdef",-1,ucf::DisplayUnits::from(opts)); }
}
