/* Differential for core::cstdio::printf_g: glibc printf %.10g/%.3g/%.1g/%g/%.17g of the
   doubles (hex bits, one per line) on stdin. Rust side: examples/printfg.rs.
   gcc -O2 printfg.c -o printfg; ./printfg < bits.txt > c.txt; printfg (rust) < bits.txt > r.txt */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
int main(){ uint64_t u; char line[64]; while (fgets(line, sizeof line, stdin)) { sscanf(line, "%lx", &u); double d; memcpy(&d,&u,8); printf("%.10g|%.3g|%.1g|%g|%.17g\n", d, d, d, d, d);} }
