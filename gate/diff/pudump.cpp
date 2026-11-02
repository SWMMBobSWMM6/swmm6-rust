// Differential for input/InputParseUtils.hpp: every parser over every distinct token in the corpus.
#include "InputParseUtils.hpp"
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <string>
using namespace openswmm::input;
static void hx(double v){ uint64_t u; std::memcpy(&u,&v,8); std::printf("%016llx ",(unsigned long long)u); }
int main(int, char** argv){
  std::ifstream in(argv[1], std::ios::binary); std::string t;
  openswmm::NameIndex ix; std::vector<std::string> errs;
  while (std::getline(in, t)) {
    double s = -7.0; bool ok = parse_double_strict(t, s);
    hx(to_double(t, -1.0)); std::printf("%d ", (int)ok); hx(s); std::printf("%d ", to_int(t, -99));
    hx(parse_date(t)); hx(parse_time_seconds(t)); hx(parse_time_day_fraction(t)); hx(parse_datetime(t, t));
    std::printf("%d\n", add_unique(ix, t, errs));
  }
  for (auto& e : errs) std::printf("%s\n", e.c_str());
}
