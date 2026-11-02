// Rust side of gate/diff/msxiq.cpp.
#[path = "common/ctxdump.rs"]
mod ctxdump;
use std::io::Write;
use swmm6_rust::core::context::SimulationContext;
use swmm6_rust::input::input_reader::InputReader;
use swmm6_rust::input::section_registry::SectionRegistry;
fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut reg = SectionRegistry::new();
    swmm6_rust::input::handlers::register_reference_handlers(&mut reg, &[], true);
    let mut ctx = SimulationContext::default();
    let mut rd = InputReader::new(&reg);
    rd.read(args[1].as_bytes(), &mut ctx);
    ctx.reactions.configured = true;
    ctx.reactions.species_name = vec![b"CL2".to_vec(), b"NH3".to_vec()];
    ctx.reactions.init_elem_is_link = vec![0, 1];
    ctx.reactions.init_elem_idx = vec![0, 0];
    ctx.reactions.init_elem_species = vec![1, 0];
    ctx.reactions.init_elem_value = vec![0.5, 9.0];
    swmm6_rust::input::post_parse_resolver::resolve_cross_references(&mut ctx);
    std::io::stdout().write_all(&ctxdump::dump(&ctx)).unwrap();
}
