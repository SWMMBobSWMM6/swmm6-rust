// Differential harness for core/DateTime.hpp: every function over a sweep of dates/times.
#include "DateTime.hpp"
#include <cstdio>
#include <cstdint>
#include <cstring>
using namespace openswmm::datetime;
static void hx(double v){ uint64_t u; std::memcpy(&u,&v,8); std::printf("%016llx ",(unsigned long long)u); }
int main(){
  for (int y=-2; y<=10002; y+= (y<1905||y>2105)? 97 : 1)
    for (int m=0; m<=13; ++m) for (int d=-1; d<=32; d+=3) { hx(encodeDate(y,m,d)); }
  std::printf("\n");
  for (int h=-1; h<=25; ++h) for (int mi=-1; mi<=61; mi+=7) for (int s=-1; s<=61; s+=13) hx(encodeTime(h,mi,s));
  std::printf("\n");
  uint64_t st=88172645463325252ull;
  for (int k=0;k<200000;++k){
    st^=st<<13; st^=st>>7; st^=st<<17;
    double x = (double)(st % 4000000000ull) / 1000.0 - 700000.0;   // -700000 .. 3.3e6 days, ms fractions
    int y,m,d,h,mi,s; decodeDate(x,y,m,d); decodeTime(x,h,mi,s);
    std::printf("%d %d %d %d %d %d %d %d ", y,m,d,h,mi,s, monthOfYear(x), dayOfYear(x));
    hx(addSeconds(x, (double)(st%100000)/7.0)); std::printf("%ld\n", timeDiff(x, x - (double)(st%977)/13.0));
  }
}
