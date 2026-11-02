// Differential for data/StorageGeometry.hpp: coefficient derivation over random parameters.
#include "StorageGeometry.hpp"
#include <cstdio>
#include <cstdint>
#include <cstring>
using namespace openswmm;
static void hx(double v){ uint64_t u; std::memcpy(&u,&v,8); std::printf("%016llx ",(unsigned long long)u); }
int main(){ uint64_t st=12345;
 for(int k=0;k<60000;++k){ st^=st<<13; st^=st>>7; st^=st<<17;
  int sh=st%6; double p[3]; for(int j=0;j<3;++j){ st^=st<<13; st^=st>>7; st^=st<<17; int m=st%10; p[j]= m==0?0.0: m==1? -1.5 : (double)(st%100000)/97.0; }
  double a=-9,b=-9,c=-9; bool ok=storage_shape_coeffs((StorageShape)sh,p[0],p[1],p[2],a,b,c);
  std::printf("%d %d %d %d ", sh, (int)storage_shape_params_valid((StorageShape)sh,p[0],p[1],p[2]), (int)ok, (int)storage_shape_is_geometric((StorageShape)sh)); hx(a);hx(b);hx(c); std::printf("%s\n", storage_shape_keyword((StorageShape)sh)); }
 const char* kws[]={"TABULAR","functional","Cylindrical","CONICAL","PARABOLOID","PYRAMIDAL","x",""};
 for(auto k:kws){ StorageShape o=StorageShape::TABULAR; bool r=storage_shape_from_keyword(k,o); std::printf("%s %d %d\n",k,(int)r,(int)o);}
 for(int v=-2;v<9;++v) std::printf("%d ",(int)storage_shape_is_valid_code(v)); std::printf("\n"); }
