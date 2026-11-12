// Rust twin of gate/diff/twod_exch.cpp, the run-time 1D-2D exchange differential (phase 2a
// wave 4b). Opens, initializes and starts the deck in the translation (twod), as the C++ side
// does in upstream's 2D build, then reads the C++ side's records on stdin: for each one it
// loads the record's inputs into its own engine (the router's state, the context's node and
// outfall arrays, the outfall exchange's ledgers), calls the translation of the recorded
// function (and the same variants), formats the outputs as the C++ side does and compares
// them with upstream's, bit for bit. The mesh, the coupling points, the options and the
// outfall exchange's donor structure are the translation's own (the start2d gate holds them
// identical to upstream's).
// Output: one line, "identical N (Q a G b ...)" or "differ N of M (...)", then up to
// TWOD_EXCH_SHOW (default 3) differing records with both sides' outputs.
// Usage: twod_exch deck.inp < upstream_records
// Build: cargo build --release --features "twod gate" --example twod_exch
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
use std::collections::BTreeMap;
use std::io::{BufRead, Write};

use swmm6_rust::core::engine::SwmmEngine;
use swmm6_rust::twod::coupling::node_coupling::{
    accumulate_outfall_discharge_step, compute_node_coupling_dqdh1d, compute_node_coupling_q, inject_accumulated_exchange, scatter_coupling_to_ydot, update_outfall_boundaries, CouplingAttachments, CouplingPoint,
};

/// One upstream record: kind, header, inputs by name, and the output text.
struct Record {
    kind: String,
    header: String,
    ins: Vec<(String, Vec<f64>)>,
    out: String,
}

impl Record {
    fn get(&self, name: &str) -> &[f64] {
        for (n, v) in &self.ins {
            if n == name {
                return v;
            }
        }
        panic!("record {} has no input {}", self.header, name)
    }
    fn has(&self, name: &str) -> bool {
        self.ins.iter().any(|(n, _)| n == name)
    }
    fn sc(&self, name: &str) -> f64 {
        self.get(name)[0]
    }
}

/// The output text, in the C++ side's format.
#[derive(Default)]
struct Out {
    s: String,
}

impl Out {
    fn arr(&mut self, name: &str, v: &[f64]) {
        self.s.push_str(&format!("o {} {}", name, v.len()));
        for x in v {
            self.s.push_str(&format!(" {:016x}", x.to_bits()));
        }
        self.s.push('\n');
    }
    fn sc(&mut self, name: &str, v: f64) {
        self.arr(name, &[v]);
    }
    fn thrown(&mut self, what: &str) {
        self.s.push_str("t ");
        self.s.push_str(what);
        self.s.push('\n');
    }
}

fn parse_hex(t: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(t, 16).expect("bad hex"))
}

fn read_record(lines: &mut impl Iterator<Item = std::io::Result<String>>, header: String) -> Record {
    let kind = header.split_whitespace().nth(1).unwrap_or("").to_string();
    let mut r = Record { kind, header, ins: Vec::new(), out: String::new() };
    for l in lines.by_ref() {
        let l = l.expect("read");
        if l == "." {
            break;
        }
        if l.starts_with("o ") || l.starts_with("t ") {
            r.out.push_str(&l);
            r.out.push('\n');
            continue;
        }
        let mut it = l.split_whitespace();
        assert_eq!(it.next(), Some("i"), "bad line {}", l);
        let name = it.next().unwrap().to_string();
        let n: usize = it.next().unwrap().parse().unwrap();
        let v: Vec<f64> = it.map(parse_hex).collect();
        assert_eq!(v.len(), n, "bad count in {}", l);
        r.ins.push((name, v));
    }
    r
}

fn ints(v: &[f64]) -> Vec<i32> {
    v.iter().map(|&x| x as i32).collect()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let show: usize = std::env::var("TWOD_EXCH_SHOW").ok().and_then(|s| s.parse().ok()).unwrap_or(3);
    let mut e = SwmmEngine::new();
    let mut rc = e.open(Some(args[1].as_bytes()), Some(b"/dev/null"), Some(b"/dev/null"), None, false);
    // A path the translation refuses (STUB_EXCEPTION) ends the comparison as "stub".
    let stubbed = |r: Result<i32, swmm6_rust::core::engine_init_hydraulics::Terminated>, phase: &str| -> i32 {
        match r {
            Ok(rc) => rc,
            Err(t) => {
                let mut out = std::io::stdout().lock();
                if t.exception == swmm6_rust::core::engine_step::STUB_EXCEPTION {
                    writeln!(out, "stub {} {}", phase, t.what).unwrap();
                } else {
                    writeln!(out, "terminated {} {} {}", phase, t.exception, t.what).unwrap();
                }
                std::process::exit(0);
            }
        }
    };
    if rc == 0 {
        rc = stubbed(e.initialize(), "init");
    }
    if rc == 0 {
        rc = stubbed(e.start(0), "start");
    }
    let stdin = std::io::stdin();
    let mut lines = stdin.lock().lines();
    let mut first = String::new();
    for l in lines.by_ref() {
        let l = l.expect("read");
        if l.starts_with("START ") {
            first = l;
            break;
        }
    }
    let r = e.surface_router_2d();
    let mine = format!("START {} active {} points {} junctions {}", rc, if r.active_ { 1 } else { 0 }, r.coupling_points_.len(), r.node_coupling_points_.len());
    let mut out = std::io::stdout().lock();
    if first.is_empty() {
        writeln!(out, "noresult").unwrap();
        return;
    }
    if first != mine {
        writeln!(out, "differ start: upstream '{}', translation '{}'", first, mine).unwrap();
        return;
    }
    if rc != 0 || !r.active_ {
        writeln!(out, "nostart {}", first).unwrap();
        return;
    }

    let mut counts: BTreeMap<String, (u64, u64)> = BTreeMap::new();
    let mut cov: BTreeMap<String, u64> = BTreeMap::new();
    let mut shown = 0usize;
    let mut differ_lines: Vec<String> = Vec::new();
    let mut end = String::new();
    while let Some(l) = lines.next() {
        let l = l.expect("read");
        if l.starts_with("END ") {
            end = l;
            break;
        }
        if !l.starts_with("R ") {
            continue;
        }
        let rec = read_record(&mut lines, l);
        coverage(&e, &rec, &mut cov);
        let mine = run(&mut e, &rec);
        let c = counts.entry(rec.kind.clone()).or_insert((0, 0));
        c.0 += 1;
        if mine.s != rec.out {
            c.1 += 1;
            if shown < show {
                shown += 1;
                differ_lines.push(format!("{}\n-- upstream\n{}-- translation\n{}", rec.header, rec.out, mine.s));
            }
        }
    }
    let total: u64 = counts.values().map(|c| c.0).sum();
    let bad: u64 = counts.values().map(|c| c.1).sum();
    let detail: Vec<String> = counts.iter().map(|(k, c)| if c.1 == 0 { format!("{} {}", k, c.0) } else { format!("{} {}/{}", k, c.1, c.0) }).collect();
    if end.is_empty() {
        writeln!(out, "noend {} records ({})", total, detail.join(" ")).unwrap();
    } else if bad == 0 {
        writeln!(out, "identical {} ({}) {}", total, detail.join(" "), end).unwrap();
    } else {
        writeln!(out, "differ {} of {} ({})", bad, total, detail.join(" ")).unwrap();
    }
    let c: Vec<String> = cov.iter().map(|(k, v)| format!("{}={}", k, v)).collect();
    writeln!(out, "cov {}", c.join(" ")).unwrap();
    for d in differ_lines {
        writeln!(out, "{}", d).unwrap();
    }
}

/// Which paths a record reaches (read from upstream's record and the translation's own
/// coupling points, before the record is run), for the gate's coverage line.
fn coverage(e: &SwmmEngine, rec: &Record, cov: &mut BTreeMap<String, u64>) {
    let mut hit = |k: &str| *cov.entry(k.to_string()).or_insert(0) += 1;
    let r = e.surface_router_2d();
    let outv = |name: &str| -> Vec<f64> {
        for l in rec.out.lines() {
            let mut it = l.split_whitespace();
            if it.next() == Some("o") && it.next() == Some(name) {
                return it.skip(1).map(parse_hex).collect();
            }
        }
        Vec::new()
    };
    match rec.kind.as_str() {
        "Q" | "G" => {
            let cp = &r.node_coupling_points_[rec.sc("k") as usize];
            let q = if rec.kind == "Q" { outv("q")[0] } else { outv("g")[0] };
            let p = if rec.kind == "Q" { "Q." } else { "G." };
            hit(&format!("{}{}", p, if q > 0.0 { "pos" } else if q < 0.0 { "neg" } else { "zero" }));
            hit(&format!("{}{}", p, if cp.law == 1 { "weir" } else { "orifice" }));
            let headcells = rec.sc("hc") != 0.0 && cp.head_end > cp.head_begin;
            let closure = r.options_.borrow().cell_closure as i32;
            hit(&format!("{}{}", p, if headcells { "headcells" } else if cp.vertex_idx < 0 { "cellhead" } else if closure == 1 { "vfrvertex" } else { "flatvertex" }));
            if rec.kind == "Q" {
                let v = outv("qv");
                if cp.equalize && rec.sc("pa") != 0.0 && rec.sc("win") > 0.0 && q != v[5] {
                    hit("Q.eqbind");
                }
                if v[6] != 0.0 {
                    hit("Q.lawflip_nonzero");
                }
                if v[0] != q {
                    hit("Q.prov_differs");
                }
                let w = outv("qw");
                if w[0] != -12345.0 {
                    hit(if closure == 1 { "Q.vertexvariant_vfr" } else { "Q.vertexvariant_flat" });
                    if w[0] != 0.0 {
                        hit("Q.vertexvariant_nonzero");
                    }
                }
            }
        }
        "CC" => {
            if !outv("pairs").is_empty() {
                hit("CC.nonempty");
            }
        }
        "AP" => {
            if rec.out.starts_with("t ") {
                hit("AP.throw");
            } else if outv("changed")[0] != 0.0 {
                hit("AP.changed");
                let acc0 = rec.get("acc");
                let acc1 = outv("acc");
                if acc1.iter().zip(acc0).any(|(a, b)| a < b) {
                    hit("AP.withdraw");
                }
                if acc1.iter().zip(acc0).any(|(a, b)| a > b) {
                    hit("AP.discharge");
                }
                if r.outfall_exchange_.points_.iter().any(|p| !p.spread.is_empty()) {
                    hit("AP.spread");
                }
            }
        }
        "PR" => {
            if outv("total").iter().any(|&t| t > 0.0) {
                hit("PR.reserved");
            }
        }
        "OB" => {
            if outv("r2d").iter().any(|&x| x > 0.0) {
                hit("OB.ramp");
            }
            if r.coupling_points_.iter().any(|cp| cp.is_outfall && cp.tw_relax > 0.0) {
                hit("OB.relax");
            }
            if r.coupling_points_.iter().any(|cp| cp.is_outfall && cp.vertex_idx >= 0) {
                hit("OB.vertex");
            }
            if rec.sc("hc") != 0.0 && r.coupling_points_.iter().any(|cp| cp.is_outfall && cp.head_end > cp.head_begin) {
                hit("OB.headcells");
            }
        }
        "TW" => {
            if outv("eng") != rec.get("eng") {
                hit("TW.engaged");
            }
            if outv("sub") != rec.get("sub") {
                hit("TW.submerged");
            }
        }
        "TP" => {
            if outv("parea") != rec.get("parea") {
                hit("TP.moved");
            }
            if outv("ndepth") != rec.get("ndepth") {
                hit("TP.rederived");
            }
        }
        "EQ" => {
            if outv("pa")[0] != 0.0 {
                hit("EQ.attached");
            }
        }
        "DW" => {
            let pa = rec.get("parea");
            let cn = rec.get("cnode");
            if !pa.is_empty() {
                hit("DW.ponds");
                let d = rec.get("ndepth");
                let full = &e.context().nodes.full_depth;
                if (0..pa.len()).any(|i| cn[i] != 0.0 && d[i] - full[i] > 0.0) {
                    hit("DW.pond_above_rim");
                }
            }
        }
        "CA" => {
            if outv("queue") != rec.get("queue") {
                hit("CA.queued");
            }
        }
        "AD" => {
            if outv("clamped").iter().any(|&c| c > 0.0) {
                hit("AD.clamped");
            }
        }
        _ => {}
    }
}

/// Load the state inputs every junction-exchange record carries.
fn load_state(e: &mut SwmmEngine, rec: &Record) {
    let (rt, ctx) = e.gate_surface_router_2d_and_context();
    let st = &mut rt.state_;
    st.depth = rec.get("depth").to_vec();
    st.head = rec.get("head").to_vec();
    st.volume = rec.get("vol").to_vec();
    st.coupling_window_s = rec.sc("win");
    st.coupling_head_cells = rec.sc("hc") != 0.0;
    if st.coupling_head_cells {
        rt.coupling_head_cells_ = ints(rec.get("hcells"));
    }
    st.node_pond_area_m2 = rec.sc("pa") != 0.0;
    if st.node_pond_area_m2 {
        rt.node_pond_area_m2_ = rec.get("pond").to_vec();
    }
    if rec.has("nhead") {
        ctx.nodes.head = rec.get("nhead").to_vec();
        ctx.nodes.depth = rec.get("ndepth").to_vec();
        ctx.nodes.invert_elev = rec.get("ninv").to_vec();
        ctx.nodes.full_depth = rec.get("nfull").to_vec();
    }
}

fn reserved_load(e: &mut SwmmEngine, res: &[f64], tot: &[f64]) {
    let (rt, _) = e.gate_surface_router_2d_and_context();
    let mut at = 0usize;
    assert_eq!(tot.len(), rt.outfall_exchange_.points_.len(), "outfall point count");
    for (pi, p) in rt.outfall_exchange_.points_.iter_mut().enumerate() {
        let n = p.cells.len();
        p.reserved = res[at..at + n].to_vec();
        at += n;
        p.total = tot[pi];
    }
    assert_eq!(at, res.len(), "outfall reservation count");
}

fn reserved_out(e: &mut SwmmEngine, o: &mut Out) {
    let (rt, _) = e.gate_surface_router_2d_and_context();
    let mut res = Vec::new();
    let mut tot = Vec::new();
    for p in &rt.outfall_exchange_.points_ {
        res.extend_from_slice(&p.reserved);
        tot.push(p.total);
    }
    o.arr("reserved", &res);
    o.arr("total", &tot);
}

/// The k-th junction point as buildCouplingPoints made it (gate/diff/twod_exch.cpp
/// original_point): its vertex, if it had one.
fn original_vertex(e: &SwmmEngine, k: usize) -> i32 {
    let r = e.surface_router_2d();
    match r.coupling_points_.iter().filter(|p| !p.is_outfall).nth(k) {
        Some(p) => p.vertex_idx,
        None => -1,
    }
}

fn flip_law(cp: &CouplingPoint) -> CouplingPoint {
    let mut c = cp.clone();
    c.law = if cp.law == 1 { 0 } else { 1 };
    c
}

/// Run the translation of one record and format its outputs.
fn run(e: &mut SwmmEngine, rec: &Record) -> Out {
    let mut o = Out::default();
    match rec.kind.as_str() {
        "Q" => {
            load_state(e, rec);
            let k = rec.sc("k") as usize;
            let ov = original_vertex(e, k);
            let hoff = rec.sc("hoff");
            let prov: Option<Vec<f64>> = if rec.sc("prov") != 0.0 { Some(rec.get("provv").to_vec()) } else { None };
            let (rt, ctx) = e.gate_surface_router_2d_and_context();
            let mesh = rt.mesh_.borrow();
            let opts = rt.options_.borrow();
            let st = &rt.state_;
            let att = CouplingAttachments::of(st, &rt.coupling_head_cells_, &rt.node_pond_area_m2_);
            let cp = &rt.node_coupling_points_[k];
            let pv = prov.as_deref();
            let q = compute_node_coupling_q(cp, &mesh, st, &att, &ctx.nodes, &opts, pv, hoff);
            o.sc("q", q);
            let quarter: Vec<f64> = st.volume.iter().map(|v| v * 0.25).collect();
            let zero = vec![0.0; st.volume.len()];
            let mut eq = cp.clone();
            eq.equalize = !cp.equalize;
            let law = flip_law(cp);
            let v = [
                compute_node_coupling_q(cp, &mesh, st, &att, &ctx.nodes, &opts, Some(&st.volume), hoff),
                compute_node_coupling_q(cp, &mesh, st, &att, &ctx.nodes, &opts, Some(&quarter), hoff),
                compute_node_coupling_q(cp, &mesh, st, &att, &ctx.nodes, &opts, Some(&zero), hoff),
                compute_node_coupling_q(cp, &mesh, st, &att, &ctx.nodes, &opts, pv, hoff + 0.05),
                compute_node_coupling_q(cp, &mesh, st, &att, &ctx.nodes, &opts, pv, hoff - 0.05),
                compute_node_coupling_q(&eq, &mesh, st, &att, &ctx.nodes, &opts, pv, hoff),
                compute_node_coupling_q(&law, &mesh, st, &att, &ctx.nodes, &opts, pv, hoff),
            ];
            o.arr("qv", &v);
            let mut w = [-12345.0; 4];
            if ov >= 0 {
                let mut vx = cp.clone();
                vx.vertex_idx = ov;
                let mut vl = vx.clone();
                vl.law = law.law;
                w[0] = compute_node_coupling_q(&vx, &mesh, st, &att, &ctx.nodes, &opts, pv, hoff);
                w[1] = compute_node_coupling_q(&vx, &mesh, st, &att, &ctx.nodes, &opts, Some(&st.volume), hoff);
                w[2] = compute_node_coupling_q(&vl, &mesh, st, &att, &ctx.nodes, &opts, pv, hoff);
                w[3] = compute_node_coupling_q(&vx, &mesh, st, &att, &ctx.nodes, &opts, Some(&quarter), hoff);
            }
            o.arr("qw", &w);
        }
        "G" => {
            load_state(e, rec);
            let k = rec.sc("k") as usize;
            let ov = original_vertex(e, k);
            let (rt, ctx) = e.gate_surface_router_2d_and_context();
            let mesh = rt.mesh_.borrow();
            let opts = rt.options_.borrow();
            let st = &rt.state_;
            let att = CouplingAttachments::of(st, &rt.coupling_head_cells_, &rt.node_pond_area_m2_);
            let cp = &rt.node_coupling_points_[k];
            o.sc("g", compute_node_coupling_dqdh1d(cp, &mesh, st, &att, &ctx.nodes, &opts));
            o.sc("gv", compute_node_coupling_dqdh1d(&flip_law(cp), &mesh, st, &att, &ctx.nodes, &opts));
            let mut w = [-12345.0; 2];
            if ov >= 0 {
                let mut vx = cp.clone();
                vx.vertex_idx = ov;
                w[0] = compute_node_coupling_dqdh1d(&vx, &mesh, st, &att, &ctx.nodes, &opts);
                vx.law = flip_law(cp).law;
                w[1] = compute_node_coupling_dqdh1d(&vx, &mesh, st, &att, &ctx.nodes, &opts);
            }
            o.arr("gw", &w);
        }
        "CC" => {
            load_state(e, rec);
            let (rt, ctx) = e.gate_surface_router_2d_and_context();
            let mut out: Vec<(i32, f64)> = Vec::new();
            rt.compute_coupling_conductances(ctx, &mut out);
            let pairs: Vec<f64> = out.iter().flat_map(|&(n, g)| [n as f64, g]).collect();
            o.arr("pairs", &pairs);
        }
        "OB" => {
            load_state(e, rec);
            let (rt, ctx) = e.gate_surface_router_2d_and_context();
            ctx.node_subtypes.outfalls.head_2d = rec.get("h2d").to_vec();
            ctx.node_subtypes.outfalls.ramp_2d = rec.get("r2d").to_vec();
            assert_eq!(rec.sc("all"), 1.0, "updateOutfallBoundaries on another point list");
            let mesh = rt.mesh_.borrow();
            let opts = rt.options_.borrow();
            let att = CouplingAttachments::of(&rt.state_, &rt.coupling_head_cells_, &rt.node_pond_area_m2_);
            update_outfall_boundaries(&rt.coupling_points_, &mesh, &rt.state_, &att, ctx, &opts, rec.sc("dt"));
            o.arr("h2d", &ctx.node_subtypes.outfalls.head_2d);
            o.arr("r2d", &ctx.node_subtypes.outfalls.ramp_2d);
        }
        "TW" => {
            let (rt, ctx) = e.gate_surface_router_2d_and_context();
            ctx.node_subtypes.outfalls.head_2d = rec.get("h2d").to_vec();
            ctx.node_subtypes.outfalls.ramp_2d = rec.get("r2d").to_vec();
            let x = &mut rt.outfall_exchange_;
            x.tw_engaged_s_ = rec.get("eng").to_vec();
            x.tw_submerged_s_ = rec.get("sub").to_vec();
            x.tw_max_ = rec.get("max").to_vec();
            x.record_tailwater(ctx, rec.sc("dt"));
            o.arr("eng", &x.tw_engaged_s_);
            o.arr("sub", &x.tw_submerged_s_);
            o.arr("max", &x.tw_max_);
        }
        "PR" => {
            reserved_load(e, rec.get("reserved"), rec.get("total"));
            let (rt, ctx) = e.gate_surface_router_2d_and_context();
            rt.state_.volume = rec.get("vol").to_vec();
            ctx.surface_outfall_link_limit = rec.get("limit").to_vec();
            let x = &mut rt.outfall_exchange_;
            x.prepared_dt_ = rec.sc("pdt");
            x.prepared_ = rec.sc("prep") != 0.0;
            x.node_reserved_ = rec.get("nres").to_vec();
            let mesh = rt.mesh_.borrow();
            let opts = rt.options_.borrow();
            x.prepare(&mesh, &rt.state_, &opts, ctx, rec.sc("dt"));
            o.sc("pdt", x.prepared_dt_);
            o.sc("prep", if x.prepared_ { 1.0 } else { 0.0 });
            o.arr("nres", &x.node_reserved_);
            o.arr("limit", &ctx.surface_outfall_link_limit);
            drop(mesh);
            drop(opts);
            reserved_out(e, &mut o);
        }
        "AP" => {
            reserved_load(e, rec.get("reserved"), rec.get("total"));
            let (rt, ctx) = e.gate_surface_router_2d_and_context();
            assert_eq!(rec.sc("ns"), rt.state_.transport.n_species as f64, "transported rows");
            let st = &mut rt.state_;
            st.volume = rec.get("vol").to_vec();
            st.depth = rec.get("depth").to_vec();
            st.head = rec.get("head").to_vec();
            st.coupling_applied = rec.get("capp").to_vec();
            st.transport.cell_runoff_vol = rec.get("runoff").to_vec();
            ctx.nodes.inflow = rec.get("inflow").to_vec();
            ctx.nodes.outflow = rec.get("outflow").to_vec();
            ctx.mass_balance_2d.outfall_in = rec.sc("mbin");
            ctx.mass_balance_2d.outfall_out = rec.sc("mbout");
            ctx.mass_balance_2d.final_storage = rec.sc("mbfs");
            let x = &mut rt.outfall_exchange_;
            x.node_reserved_ = rec.get("nres").to_vec();
            x.pending_ = rec.get("pend").to_vec();
            x.prepared_ = rec.sc("prep") != 0.0;
            let mut acc = rec.get("acc").to_vec();
            let conc: Option<Vec<f64>> = if rec.sc("nc") != 0.0 { Some(rec.get("conc").to_vec()) } else { None };
            let mut cacc: Option<Vec<f64>> = if rec.sc("ca") != 0.0 { Some(rec.get("cacc").to_vec()) } else { None };
            let mesh = rt.mesh_.borrow();
            let opts = rt.options_.borrow();
            let r = x.apply(&mesh, st, &opts, ctx, rec.sc("dt"), &mut acc, conc.as_deref(), cacc.as_deref_mut());
            let outs = |o: &mut Out| {
                o.arr("vol", &st.volume);
                o.arr("depth", &st.depth);
                o.arr("head", &st.head);
                o.arr("capp", &st.coupling_applied);
                o.arr("pend", &x.pending_);
                o.arr("acc", &acc);
                o.sc("mbin", ctx.mass_balance_2d.outfall_in);
                o.sc("mbout", ctx.mass_balance_2d.outfall_out);
                o.sc("mbfs", ctx.mass_balance_2d.final_storage);
                o.sc("prep", if x.prepared_ { 1.0 } else { 0.0 });
                o.arr("runoff", &st.transport.cell_runoff_vol);
                if let Some(c) = &cacc {
                    o.arr("cacc", c);
                }
            };
            match r {
                Ok(changed) => {
                    o.sc("changed", if changed { 1.0 } else { 0.0 });
                    outs(&mut o);
                }
                Err(t) => {
                    o.thrown(t.what);
                    outs(&mut o);
                }
            }
            o.sc("recv", x.recv_.len() as f64);
        }
        "TP" => {
            let (rt, ctx) = e.gate_surface_router_2d_and_context();
            rt.state_.depth = rec.get("depth").to_vec();
            ctx.coupled_pond_area = rec.get("parea").to_vec();
            ctx.coupled_pond_band = rec.get("pband").to_vec();
            ctx.nodes.depth = rec.get("ndepth").to_vec();
            ctx.nodes.head = rec.get("nhead").to_vec();
            rt.update_tracking_ponds(ctx, rec.sc("dt"));
            o.arr("parea", &ctx.coupled_pond_area);
            o.arr("ndepth", &ctx.nodes.depth);
            o.arr("nhead", &ctx.nodes.head);
        }
        "EQ" => {
            let (rt, ctx) = e.gate_surface_router_2d_and_context();
            ctx.nodes.depth = rec.get("ndepth").to_vec();
            ctx.coupled_pond_area = rec.get("parea").to_vec();
            ctx.coupled_pond_band = rec.get("pband").to_vec();
            rt.node_pond_area_m2_ = rec.get("pond").to_vec();
            rt.publish_equalize_areas(ctx, rec.sc("win"));
            let pa = if rt.state_.node_pond_area_m2 { 1.0 } else { 0.0 };
            o.sc("pa", pa);
            o.sc("pa_target", pa);
            o.sc("win", rt.state_.coupling_window_s);
            o.arr("pond", &rt.node_pond_area_m2_);
        }
        "SY" => {
            load_state(e, rec);
            let k = rec.sc("k") as usize;
            let ov = original_vertex(e, k);
            let (rt, ctx) = e.gate_surface_router_2d_and_context();
            let mesh = rt.mesh_.borrow();
            let opts = rt.options_.borrow();
            let st = &rt.state_;
            let att = CouplingAttachments::of(st, &rt.coupling_head_cells_, &rt.node_pond_area_m2_);
            let cp = &rt.node_coupling_points_[k];
            let qs = [compute_node_coupling_q(cp, &mesh, st, &att, &ctx.nodes, &opts, None, 0.0), 1.0e-3, -1.0e-3];
            o.arr("qs", &qs);
            for q in qs {
                let mut ydot = vec![0.0; st.volume.len()];
                scatter_coupling_to_ydot(&mesh, st, cp, q, &mut ydot);
                o.arr("ydot", &ydot);
            }
            if ov >= 0 {
                let mut vx = cp.clone();
                vx.vertex_idx = ov;
                for q in qs {
                    let mut ydot = vec![0.0; st.volume.len()];
                    scatter_coupling_to_ydot(&mesh, st, &vx, q, &mut ydot);
                    o.arr("ydotv", &ydot);
                }
            }
        }
        "AD" => {
            let (rt, ctx) = e.gate_surface_router_2d_and_context();
            rt.state_.volume = rec.get("vol").to_vec();
            rt.state_.head = rec.get("head").to_vec();
            ctx.nodes.inflow = rec.get("inflow").to_vec();
            ctx.nodes.outflow = rec.get("outflow").to_vec();
            let mesh = rt.mesh_.borrow();
            let opts = rt.options_.borrow();
            let all = &rt.coupling_points_;
            let dt = rec.sc("dt");
            for pass in 0..2 {
                let mut accum = vec![0.0; all.len()];
                let mut budget = rt.state_.volume.clone();
                let mut sample = vec![-7.0; all.len()];
                if pass == 1 {
                    for b in budget.iter_mut() {
                        *b *= 0.25;
                    }
                }
                let clamped = accumulate_outfall_discharge_step(all, &mesh, &rt.state_, ctx, &opts, dt, &mut accum, &mut budget, if pass == 0 { Some(&mut sample[..]) } else { None });
                o.sc("clamped", clamped as f64);
                o.arr("accum", &accum);
                o.arr("budget", &budget);
                o.arr("sample", &sample);
                for sign in [1.0, -1.0] {
                    let mut cp = rt.state_.clone();
                    cp.coupling_flux.iter_mut().for_each(|f| *f = 0.0);
                    let mut acc2 = accum.clone();
                    for (k, a) in acc2.iter_mut().enumerate() {
                        if *a == 0.0 {
                            *a = 1.0e-3 * (k + 1) as f64;
                        }
                    }
                    inject_accumulated_exchange(all, &mesh, &mut cp, &acc2, dt, sign, None);
                    o.arr("flux", &cp.coupling_flux);
                }
            }
        }
        "CA" => {
            let (rt, ctx) = e.gate_surface_router_2d_and_context();
            assert_eq!(rec.sc("common"), if rt.common_stage_ { 1.0 } else { 0.0 }, "common source stage");
            assert_eq!(rec.sc("ns"), rt.state_.transport.n_species as f64, "transported rows");
            ctx.nodes.coupling_queue = rec.get("queue").to_vec();
            ctx.nodes.coupling_volume = rec.get("cvol").to_vec();
            ctx.coupling_delivery_remaining = rec.sc("cdr");
            let exch = rec.get("exch").to_vec();
            rt.co_advance_junction_ledger(ctx, &exch);
            rt.co_advance_junction_queue(ctx, &exch, rec.sc("dt"));
            o.arr("queue", &ctx.nodes.coupling_queue);
            o.arr("cvol", &ctx.nodes.coupling_volume);
            o.sc("cdr", ctx.coupling_delivery_remaining);
        }
        "DW" => {
            let (rt, ctx) = e.gate_router_and_context();
            ctx.options.allow_ponding = rec.sc("pond") != 0.0;
            ctx.coupled_node = rec.get("cnode").iter().map(|&x| x as u8).collect();
            ctx.coupled_pond_area = rec.get("parea").to_vec();
            ctx.coupled_pond_band = rec.get("pband").to_vec();
            ctx.nodes.depth = rec.get("ndepth").to_vec();
            ctx.nodes.lat_flow = rec.get("lat").to_vec();
            ctx.nodes.losses = rec.get("loss").to_vec();
            let (surf, dqdh, conv) = rt.dw_solver_.gate_init_node_states(ctx);
            o.arr("surf", &surf);
            o.arr("dqdh", &dqdh);
            o.arr("conv", &conv);
            o.arr("inflow", &ctx.nodes.inflow);
            o.arr("outflow", &ctx.nodes.outflow);
        }
        "TQ" => {
            let (rt, ctx) = e.gate_surface_router_2d_and_context();
            let mut v = Vec::new();
            for n in 0..ctx.n_nodes() {
                let (mut a, mut b, mut c) = (-1.0, -2.0, -3.0);
                let ok = rt.outfall_exchange_.tailwater(n, &mut a, &mut b, &mut c);
                v.extend_from_slice(&[if ok { 1.0 } else { 0.0 }, a, b, c]);
            }
            o.arr("tw", &v);
        }
        k => panic!("unknown record kind {}", k),
    }
    o
}
