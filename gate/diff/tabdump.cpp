// Differential harness for the TableData.hpp lookups: upstream's inline functions over
// tables and queries read from a binary file. Output: one line per query, every result's
// bits in hex, so a one-ulp difference shows.
#include "TableData.hpp"
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <vector>
using namespace openswmm;
static void hx(double v){ uint64_t u; std::memcpy(&u,&v,8); std::printf("%016llx ",(unsigned long long)u); }
int main(int argc,char**argv){
  FILE*f=std::fopen(argv[1],"rb"); int nt; std::fread(&nt,4,1,f);
  for(int t=0;t<nt;++t){
    int n,nq; std::fread(&n,4,1,f); Table tb; tb.x.resize(n); tb.y.resize(n);
    std::fread(tb.x.data(),8,n,f); std::fread(tb.y.data(),8,n,f); std::fread(&nq,4,1,f);
    std::vector<double> q(nq); std::fread(q.data(),8,nq,f);
    Table a=tb,b=tb,c=tb;   // independent cursors per cursor function
    for(double x:q){
      hx(table_lookup_cursor(a,x)); hx(table_tseries_lookup_cursor(b,x)); hx(table_step_cursor(c,x));
      hx(table_getStorageVolume(tb,x)); hx(table_getStorageDepth(tb,x)); hx(table_inverseLookup(tb,x));
      hx(table_lookupEx(tb,x)); hx(table_intervalLookup(tb,x)); hx(table_getSlope(tb,x)); hx(table_getMaxY(tb));
      std::printf("%d %d %d\n",a.cursor.index,b.cursor.index,c.cursor.index);
    }
  }
}
