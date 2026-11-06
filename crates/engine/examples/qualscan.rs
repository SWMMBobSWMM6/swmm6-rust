// Which parts of the routing step's quality stage a deck reaches: open the deck in the
// translation, run the init legs that size the solvers, and print one line of the flags that
// select upstream's quality paths (SWMMEngine::stepRouting B4b/B5 and QualitySolver::execute).
// Usage: qualscan deck.inp
use swmm6_rust::core::engine::SwmmEngine;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut e = SwmmEngine::new();
    let rc = e.open(Some(args[1].as_bytes()), Some(b"/dev/null"), Some(b"/dev/null"), None, false);
    if rc != 0 {
        println!("open {} {}", rc, args[1]);
        return;
    }
    let _ = e.gate_init_hydraulics();
    let _ = e.gate_init_quality();
    let ninlet = e.gate_solvers().inlet_.soa_.count;
    let c = e.context();
    let lid_rows = c.node_subtypes.storages.lid.iter().filter(|l| l.control >= 0).count();
    println!(
        "np={} route={:?} qs={:?} age={} heat={} bed={} rx={} nsp={} rxlid={} surf={} treat={} nodelid={} inlets={} lqforce={} ardvol={} coupling={} {}",
        c.n_pollutants(),
        c.options.routing_model,
        c.options.quality_solver,
        c.options.water_age as i32,
        c.options.heat_transport as i32,
        c.heat_config.sediment_exchange as i32,
        swmm6_rust::transport::reaction_legacy_binding::legacy_reactions_active(c) as i32,
        c.reactions.n_species(),
        c.reactions.lid_controls.len(),
        c.reactions.surface.active() as i32,
        c.treatment.has_any() as i32,
        lid_rows,
        ninlet,
        c.forcing.link_quality_mode.len(),
        swmm6_rust::transport::ard_config::ard_boundaries_need_external_volumes(c) as i32,
        c.nodes.coupling_inflow.len(),
        args[1]
    );
}
