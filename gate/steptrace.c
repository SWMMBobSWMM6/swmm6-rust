/* steptrace: run a deck through OpenSWMM 6's engine library and record the state
 * after every routing step, so a translation can be compared step by step and
 * the first step and object where it departs can be found.
 *
 * Trace file: "S6TR" magic, int32 version 1, int32 nsub, nnode, nlink, then per
 * step a float64 elapsed time (days) followed by float64 arrays: subcatchment
 * runoff[nsub], node depth[nnode], node volume[nnode], link flow[nlink],
 * link depth[nlink]. Little-endian, as written by the host.
 *
 * Usage: steptrace deck.inp out.rpt out.out trace.bin
 * Build: cc -O2 steptrace.c -I<openswmm>/include/openswmm/engine -L<lib> -lopenswmm.engine
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "openswmm_engine.h"
#include "openswmm_nodes.h"
#include "openswmm_links.h"
#include "openswmm_subcatchments.h"

static void put(FILE *f, const double *v, int n) { if (n > 0) fwrite(v, sizeof(double), (size_t)n, f); }

int main(int argc, char **argv)
{
    if (argc != 5) { fprintf(stderr, "usage: steptrace deck.inp out.rpt out.out trace.bin\n"); return 2; }
    SWMM_Engine e = swmm_engine_create();
    int rc = swmm_engine_open(e, argv[1], argv[2], argv[3], NULL);
    if (rc != SWMM_OK) { fprintf(stderr, "open failed: %d\n", rc); return 3; }
    rc = swmm_engine_initialize(e);
    if (rc != SWMM_OK) { fprintf(stderr, "initialize failed: %d\n", rc); return 3; }
    rc = swmm_engine_start(e, 1);
    if (rc != SWMM_OK) { fprintf(stderr, "start failed: %d\n", rc); return 3; }
    int ns = swmm_subcatch_count(e), nn = swmm_node_count(e), nl = swmm_link_count(e);
    FILE *f = fopen(argv[4], "wb");
    if (!f) { perror(argv[4]); return 4; }
    int hdr[4] = {1, ns, nn, nl};
    fwrite("S6TR", 1, 4, f); fwrite(hdr, sizeof(int), 4, f);
    int m = ns > nn ? ns : nn; if (nl > m) m = nl; if (m < 1) m = 1;
    double *buf = malloc(sizeof(double) * (size_t)m);
    double elapsed = 0.0;
    long steps = 0;
    while ((rc = swmm_engine_step(e, &elapsed)) == SWMM_OK && elapsed > 0.0) {
        fwrite(&elapsed, sizeof(double), 1, f);
        swmm_subcatch_get_runoff_bulk(e, buf, ns); put(f, buf, ns);
        swmm_node_get_depths_bulk(e, buf, nn);     put(f, buf, nn);
        swmm_node_get_volumes_bulk(e, buf, nn);    put(f, buf, nn);
        swmm_link_get_flows_bulk(e, buf, nl);      put(f, buf, nl);
        swmm_link_get_depths_bulk(e, buf, nl);     put(f, buf, nl);
        steps++;
    }
    fclose(f);
    swmm_engine_end(e);
    swmm_engine_report(e);
    swmm_engine_close(e);
    swmm_engine_destroy(e);
    fprintf(stderr, "%ld steps, last rc %d\n", steps, rc);
    return rc == SWMM_OK ? 0 : 1;
}
