// Where does a coupled 1D-2D run's 1D flow-routing continuity error come from? Runs a deck on
// the 2D reference (the C API handle IS the SWMMEngine) and after every step compares the 1D
// ledger (cumulative booked inflow + initial storage - booked outflow) with the volume the 1D
// network actually holds (the same sum computeFinalStorage makes: reportedNodeVolume over the
// nodes plus links.volume), printing the running discrepancy and each coupled node's exchange.
// A diagnostic, not a gate: it reads upstream's state and changes nothing.
// Build: as gate/build.sh's 2D harnesses (build2d, -DOPENSWMM_HAS_2D=1), e.g.
//   g++ -std=c++20 -O2 -ffp-contract=off -mavx2 -mfma -DNDEBUG -DSWMM_OS_LINUX -DSWMM_USE_OPENMP
//       -fopenmp -DOPENSWMM_HAS_2D=1 -I$U/src/engine -I$U/include -I$B2/include
//       gate/diff/couplebal.cpp -o gate/diff/couplebal $B2/src/engine/libopenswmm.engine.so
//       -Wl,-rpath,$B2/src/engine
// Usage: couplebal deck.inp [every]   (print every Nth step; default 1, and every change)
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
#define private public
#define protected public
#include "core/SWMMEngine.hpp"
#undef private
#undef protected
#include <openswmm/engine/openswmm_engine.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace openswmm;

static double stored(SWMMEngine& e) {
    auto& c = e.ctx_;
    double v = 0.0;
    for (int j = 0; j < c.n_nodes(); ++j) v += e.reportedNodeVolume(j);
    for (int j = 0; j < c.n_links(); ++j) v += c.links.volume[static_cast<std::size_t>(j)];
    return v;
}
static double raw_node_volume(SWMMEngine& e) {
    auto& c = e.ctx_;
    double v = 0.0;
    for (int j = 0; j < c.n_nodes(); ++j) v += c.nodes.volume[static_cast<std::size_t>(j)];
    return v;
}

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    const int every = argc > 2 ? std::atoi(argv[2]) : 1;
    // optional: argv[3] = a file to write each node's coupling_inflow per step as an
    // [INFLOWS]-ready time series (cfs), for replaying the exchange into a 1D-only deck
    FILE* ts = argc > 3 ? std::fopen(argv[3], "w") : nullptr;
    SWMM_Engine h = swmm_engine_create();
    int rc = swmm_engine_open(h, argv[1], "/dev/null", "/dev/null", nullptr);
    if (rc == 0) rc = swmm_engine_initialize(h);
    if (rc == 0) rc = swmm_engine_start(h, 0);
    if (rc != 0) { std::printf("rc %d\n", rc); return 1; }
    auto& e = *reinterpret_cast<SWMMEngine*>(h);
    auto& c = e.ctx_;
    auto& mb = c.mass_balance;
    std::vector<int> coupled;
    for (int j = 0; j < c.n_nodes(); ++j) {
        // a node the 2D router exchanges with: anything that ever carries a coupling rate is
        // printed once it does; start with every node and filter at print time
        coupled.push_back(j);
    }
    std::printf("step time(s) dt  booked_in booked_out stored(rep) stored(raw) disc(rep) disc(raw) | per node: name coupling_inflow lat_flow depth volume\n");
    double elapsed = 0.0, last_disc = 0.0;
    int step = 0;
    const double v0_rep = stored(e), v0_raw = raw_node_volume(e);
    for (;;) {
        double t0 = c.current_time;
        rc = swmm_engine_step(h, &elapsed);
        if (rc != 0 || elapsed <= 0.0) break;
        ++step;
        const double in = mb.routing_dry_weather + mb.routing_wet_weather + mb.routing_gw_inflow +
                          mb.routing_link_gw_inflow + mb.routing_rdii + mb.routing_external;
        const double out = mb.routing_flooding + mb.routing_coupling_out + mb.routing_outflow +
                           mb.routing_evap_loss + mb.routing_seep_loss;
        const double s_rep = stored(e);
        double s_raw = raw_node_volume(e);
        for (int j = 0; j < c.n_links(); ++j) s_raw += c.links.volume[static_cast<std::size_t>(j)];
        const double disc_rep = (mb.routing_init_storage + in - out) - s_rep;
        const double disc_raw = (v0_raw + in - out) - s_raw;
        const bool jump = std::fabs(disc_rep - last_disc) > 1e-6;
        if (step % every == 0 || jump) {
            std::printf("%d %.3f %.3f %.6f %.6f %.6f %.6f %.6f %.6f |", step, c.current_time * 86400.0,
                        (c.current_time - t0) * 86400.0, in, out, s_rep, s_raw, disc_rep, disc_raw);
            for (int j : coupled) {
                auto uj = static_cast<std::size_t>(j);
                if (c.nodes.coupling_inflow[uj] != 0.0 || c.nodes.volume[uj] != 0.0)
                    std::printf(" %s %.6f %.6f %.6f %.6f", c.node_names.name_of(j).c_str(),
                                c.nodes.coupling_inflow[uj], c.nodes.lat_flow[uj], c.nodes.depth[uj],
                                c.nodes.volume[uj]);
            }
            std::printf("%s\n", jump ? "  <-" : "");
        }
        last_disc = disc_rep;
        if (ts) {
            for (int j : coupled) {
                auto uj = static_cast<std::size_t>(j);
                std::fprintf(ts, "%s %d %.17g\n", c.node_names.name_of(j).c_str(), step,
                             c.nodes.coupling_inflow[uj]);
            }
        }
    }
    if (ts) std::fclose(ts);
    swmm_engine_end(h);
    std::printf("final: init %.6f in %.6f out %.6f final(computeFinalStorage) %.6f error %.4f%%  v0_rep %.6f\n",
                mb.routing_init_storage,
                mb.routing_dry_weather + mb.routing_wet_weather + mb.routing_gw_inflow + mb.routing_link_gw_inflow +
                    mb.routing_rdii + mb.routing_external,
                mb.routing_flooding + mb.routing_coupling_out + mb.routing_outflow + mb.routing_evap_loss +
                    mb.routing_seep_loss,
                mb.routing_final_storage, mb.routing_error() * 100.0, v0_rep);
    swmm_engine_close(h);
    swmm_engine_destroy(h);
    return 0;
}
