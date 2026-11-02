// Rust side of gate/diff/ctxread.cpp: parse with the translated reader and handlers,
// optionally resolve (argv[3] == "resolve"), dump.
#[path = "common/ctxdump.rs"]
mod ctxdump;
use std::io::Write;
use swmm6_rust::core::context::SimulationContext;
use swmm6_rust::input::input_reader::InputReader;
use swmm6_rust::input::section_registry::SectionRegistry;
fn main() {
    let args: Vec<String> = std::env::args().collect();
    let all = args.len() > 2 && args[2] == "ALL";
    let on: Vec<String> = if args.len() > 2 && !all { args[2].split(',').filter(|s| !s.is_empty()).map(String::from).collect() } else { Vec::new() };
    let resolve = args.len() > 3 && args[3] == "resolve";
    let mut reg = SectionRegistry::new();
    swmm6_rust::input::handlers::register_reference_handlers(&mut reg, &on, all);
    let mut ctx = SimulationContext::default();
    let mut rd = InputReader::new(&reg);
    rd.read(args[1].as_bytes(), &mut ctx);
    if resolve {
        swmm6_rust::input::post_parse_resolver::resolve_cross_references(&mut ctx);
    }
    std::io::stdout().write_all(&ctxdump::dump(&ctx)).unwrap();
}
