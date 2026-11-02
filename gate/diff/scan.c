/* Differential for data::table::sscanf (and cstdlib::strtof) against glibc, with every
   format input/PostParseResolver.cpp uses. Lines on stdin ("~" stands for an embedded
   newline). Rust side: the sscanf_gate test in core/cstdio.rs:
   gcc -O1 scan.c -o scan; ./scan < in.txt > c.txt
   SSCANF_GATE_IN=in.txt (no final newline) SSCANF_GATE_OUT=r.txt cargo test --release --lib sscanf_against_glibc
   cmp c.txt r.txt */
/* glibc sscanf with every format the rain-file readers use; prints ret and all outputs. */
#include <stdio.h>
#include <string.h>
#include <stdint.h>
static void pfl(float f){ uint32_t u; memcpy(&u,&f,4); printf(" f%08x", u); }
int main(int argc, char** argv) {
  char line[8192];
  while (fgets(line, sizeof line, stdin)) {
    size_t n = strlen(line); if (n && line[n-1]=='\n') line[--n]=0;
    for (size_t i=0;i<n;i++) if (line[i]=='~') line[i]='\n';  /* allow embedded newlines */
    { unsigned m=7,d=7,y=7; char s1='#',s2='#'; int r=sscanf(line,"%u%c%u%c%u",&m,&s1,&d,&s2,&y); printf("A %d %u %d %u %d %u\n",r,m,s1,d,s2,y); }
    { char mon[64]="?"; unsigned d=7,y=7; char s1='#',s2='#'; int r=sscanf(line,"%3[A-Za-z]%c%u%c%u",mon,&s1,&d,&s2,&y); printf("B %d %s %d %u %d %u\n",r,mon,s1,d,s2,y); }
    { int h=7,m=7,s=7; int r=sscanf(line,"%d:%d:%d",&h,&m,&s); printf("C %d %d %d %d\n",r,h,m,s); }
    { long a=7; int b=7; char e[64]="?"; int r=sscanf(line,"%6ld %2d %4s",&a,&b,e); printf("D %d %ld %d %s\n",r,a,b,e); }
    { long a=7; int b=7; char e[64]="?"; int r=sscanf(line,"%6ld,%2d,%4s",&a,&b,e); printf("E %d %ld %d %s\n",r,a,b,e); }
    { int b=7,y=7; char e[64]="?", rc[64]="?"; int r=sscanf(line,"%2d %4s %2s %4d",&b,e,rc,&y); printf("F %d %d %s %s %d\n",r,b,e,rc,y); }
    { int b=7,y=7; char e[64]="?", rc[64]="?"; int r=sscanf(line,"%2d,%4s,%2s,%4d",&b,e,rc,&y); printf("G %d %d %s %s %d\n",r,b,e,rc,y); }
    { long a=7; int b=7; char e[64]="?", rc[64]="?"; int r=sscanf(line,"%3s%6ld%2d%4s",rc,&a,&b,e); printf("H %d %s %ld %d %s\n",r,rc,a,b,e); }
    { long a=7; char c[64]="?"; int r=sscanf(line,"%5s%6ld",c,&a); printf("I %d %s %ld\n",r,c,a); }
    { float x=99.99f; char c='#'; int r=sscanf(line,"%f %c",&x,&c); printf("J %d %d",r,c); pfl(x); printf("\n"); }
    { long v=7; char c='#'; int r=sscanf(line,"%ld %c",&v,&c); printf("K %d %ld %d\n",r,v,c); }
    { int a=7,b=7,c=7,d=7; int r=sscanf(line,"%4d%2d%4d%3d",&a,&b,&c,&d); printf("L %d %d %d %d %d\n",r,a,b,c,d); }
    { int a=7,b=7,c=7; int r=sscanf(line,"%4d %2d %2d",&a,&b,&c); printf("M %d %d %d %d\n",r,a,b,c); }
    { int a=7,b=7,c=7; int r=sscanf(line,"%4d,%2d,%2d",&a,&b,&c); printf("N %d %d %d %d\n",r,a,b,c); }
    { int a=7,b=7,c=7; int r=sscanf(line,"%4d%2d%2d",&a,&b,&c); printf("O %d %d %d %d\n",r,a,b,c); }
    { int h=7,m=7; long v=7; char f1='#',f2='#'; int r=sscanf(line,"%2d%2d%6ld%c%c",&h,&m,&v,&f1,&f2); printf("P %d %d %d %ld %d %d\n",r,h,m,v,f1,f2); }
    { int h=7,m=7; long v=7; char f1='#',f2='#'; int r=sscanf(line," %2d%2d %6ld %c %c",&h,&m,&v,&f1,&f2); printf("Q %d %d %d %ld %d %d\n",r,h,m,v,f1,f2); }
    { int h=7,m=7; long v=7; char f1='#',f2='#'; int r=sscanf(line,",%2d%2d,%6ld,%c,%c",&h,&m,&v,&f1,&f2); printf("R %d %d %d %ld %d %d\n",r,h,m,v,f1,f2); }
    { int h=7,m=7; int r=sscanf(line," %2d:%2d",&h,&m); printf("S %d %d %d\n",r,h,m); }
    { char t[64]="?"; int a=7,b=7,c=7,d=7,e=7; float v=-1.5f; int r=sscanf(line,"%63s %d %d %d %d %d %f",t,&a,&b,&c,&d,&e,&v); printf("T %d %s %d %d %d %d %d",r,t,a,b,c,d,e); pfl(v); printf("\n"); }
    { int m=7,d=7,y=7,e=-9; int r=sscanf(line,"%d/%d/%d%n",&m,&d,&y,&e); printf("U %d %d %d %d %d\n",r,m,d,y,e); }
    { char s[64]="?"; int r=sscanf(line,"%3[^a-c-]",s); printf("V %d %s\n",r,s); }
    { char s[64]="?"; int r=sscanf(line,"%[]a-]",s); printf("W %d %s\n",r,s); }
    { float x=5.f; int r=sscanf(line,"%3f",&x); printf("X %d",r); pfl(x); printf("\n"); }
  }
}
