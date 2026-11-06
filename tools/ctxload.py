#!/usr/bin/env python3
"""Generate the Rust context loader, crates/engine/examples/common/ctxload.rs: it reads the
text that gate/diff/ctxdump.hpp (C++) or examples/common/ctxdump.rs (Rust) prints and fills a
fresh SimulationContext and SwmmEngine, so that dumping them again gives the same bytes.

The loader is generated from tools/ctxdump.py's own field lists and kind table, walked in the
dumper's order. Two checks keep the twins from drifting:
  * at generation, the walk re-derives the Rust dumper with ctxdump.py's own emit/dump_struct
    and hand-written lines and compares it, line for line, with what ctxdump.py writes; for
    each generated segment it also compares the sequence of labels the dumper prints with the
    sequence the loader reads;
  * at load time, every line's label is checked against the expected one, and the input must
    end exactly where the dump ends.
Members the dumper skips ("skip <label>") are consumed and left at their defaults; the list is
printed. Usage: tools/ctxload.py   (CTXDUMP_ROOT as for ctxdump.py)
"""
import io, os, re, sys, tempfile, contextlib

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ctxdump as D  # noqa: E402

ROOT = os.environ.get('CTXDUMP_ROOT', '/root/work/swmm6-rust')
SRC = f'{ROOT}/crates/engine/src'
OUT = f'{ROOT}/crates/engine/examples/common/ctxload.rs'

# Element structs of vectors that are not Default in the translation: Rust constructor.
RUST_NEW = {'Table': 'swmm6_rust::data::table::Table::new(b"", swmm6_rust::data::table::TableType::Timeseries)'}

CNT = [0]
def fresh(p):
    CNT[0] += 1
    return f'{p}{CNT[0]}'

SKIPPED = []

def skip(lb, L):
    SKIPPED.append(lb)
    L.append(f'    r.sk("{lb}");')

def load(p, label, k, L):
    """Loader statements filling place `p` with one value of kind k (mirror of ctxdump.emit)."""
    if k[0] == 'f':
        L.append(f'    {p} = FromF64::from_f64(r.pf("{label}"));')
    elif k[0] == 'i':
        L.append(f'    {p} = FromI64::from_i64(r.pi("{label}"));')
    elif k[0] == 's':
        if k[1] == 'FilePathPair':
            L.append(f'    ({p}).original = r.ps("{label}");')
            L.append(f'    ({p}).absolute = r.ps("{label}.abs");')
        else:
            L.append(f'    {p} = r.ps("{label}");')
    elif k[0] == 'v' and k[1][0] in 'fis':
        ik = k[1]
        if ik[0] == 'i':
            L.append(f'    set_vi(&mut {p}, r.vi("{label}"));')
        elif ik[0] == 'f':
            L.append(f'    set_vf(&mut {p}, r.vf("{label}"));')
        elif ik[1] == 'FilePathPair':
            L.append(f'    {{ let o = r.vs("{label}"); let a = r.vs("{label}.abs"); set_fpp(&mut {p}, o, a); }}')
        else:
            L.append(f'    SeqSet::seq_set(&mut {p}, r.vs("{label}"));')
    elif k[0] == 'va':
        L.append(f'    set_va(&mut {p}, r.vf("{label}"));')
    elif k[0] in ('v', 'arr'):
        e = fresh('e')
        L.append(f'    {{ let n = r.count("{label}.n"); SeqLen::seq_len(&mut {p}, n); for {e} in ({p}).iter_mut() {{')
        load(f'(*{e})', label + '[]', k[1], L)
        L.append('    } }')
    elif k[0] == 'opt':
        e = fresh('e')
        L.append(f'    if !r.at_none("{label}") {{ let mut {e} = opt_new(&{p});')
        load(e, label, k[1], L)
        L.append(f'    {p} = Some({e}); }} else {{ r.pn("{label}"); {p} = None; }}')
    elif k[0] == 'map':
        e = fresh('e')
        L.append(f'    {{ let n = r.count("{label}.n"); MapLoad::map_clear(&mut {p}); for _ in 0..n {{')
        L.append(f'    let mut {e}k = map_k(&{p});')
        load(f'{e}k', label + '.key', k[1], L)
        L.append(f'    let mut {e}v = map_v(&{p});')
        load(f'{e}v', label + '.val', k[2], L)
        L.append(f'    MapLoad::map_put(&mut {p}, {e}k, {e}v); }} }}')
    elif k[0] == 'umap':
        e = fresh('e')
        L.append(f'    {{ let n = r.count("{label}.n"); MapLoad::map_clear(&mut {p}); for _ in 0..n {{')
        L.append(f'    let mut {e}k = map_k(&{p});')
        load(f'{e}k', label + '.key', k[1], L)
        L.append(f'    let mut {e}v = map_v(&{p});')
        load(f'{e}v', label + '.val', k[2], L)
        L.append(f'    MapLoad::map_put(&mut {p}, {e}k, {e}v); }} }}')
    elif k[0] == 'pair':
        load(f'{p}.0', label + '.first', k[1], L)
        load(f'{p}.1', label + '.second', k[2], L)
    elif k[0] == 'struct':
        load_struct(k[1], k[2], p, label, L)
    else:
        raise SystemExit(f'ctxload: no loader for kind {k} ({label})')

def load_struct(header, struct, path_r, label, L):
    """Mirror of ctxdump.dump_struct."""
    for name, ty, arr in D.fields(header, struct):
        k = D.kind(ty)
        pr, lb = f'{path_r}.{D.rname(name)}', f'{label}.{name}'
        if arr:
            if k and k[0] in 'fi':
                L.append(f'    for q in 0..({pr}).len() {{')
                load(f'{pr}[q]', lb + '[]', k, L)
                L.append('    }')
            elif k and k[0] == 'struct':
                q = fresh('q')
                L.append(f'    for {q} in 0..({pr}).len() {{')
                load_struct(k[1], k[2], f'{pr}[{q}]', lb + '[]', L)
                L.append('    }')
            else:
                skip(lb, L)
            continue
        em = re.fullmatch(r'std::vector<(\w+)>', ty.replace('const ', '').strip())
        if k is None and em and em.group(1) in D.ELEM_STRUCTS:
            new = f'seq_len_with(&mut {pr}, n, || {RUST_NEW[em.group(1)]})' if em.group(1) in RUST_NEW else f'SeqLen::seq_len(&mut {pr}, n)'
            L.append(f'    {{ let n = r.count("{lb}.n"); {new}; for el in ({pr}).iter_mut() {{')
            load_struct(D.ELEM_STRUCTS[em.group(1)], em.group(1), 'el', lb + '[]', L)
            L.append('    } }')
            continue
        if k is None:
            SKIPPED.append(f'{lb} ({ty})')
            L.append(f'    r.sk("{lb}");')
            continue
        load(pr, lb, k, L)

def load_ctx_member(n, t, a, L):
    """Mirror of ctxdump.ctx_member."""
    k = D.kind(t)
    if a:
        q = fresh('q')
        L.append(f'    for {q} in 0..(c.{n}).len() {{')
        if k[0] == 'struct':
            load_struct(k[1], k[2], f'c.{n}[{q}]', n + '[]', L)
        else:
            load(f'c.{n}[{q}]', n + '[]', k, L)
        L.append('    }')
    else:
        load(f'c.{n}', n, k, L)

# ---------------------------------------------------------------- label-sequence check
R_CALL = re.compile(r'\b(pf|pi|ps|vf|vi|vs|sk|pn)\(o, "([^"]*)"')
L_CALL = re.compile(r'\br\.(pf|pi|ps|vf|vi|vs|sk|pn|count)\("([^"]*)"\)')

def labels_r(lines):
    return [m.groups() for ln in lines for m in R_CALL.finditer(ln)]

def labels_l(lines):
    return [('pi' if f == 'count' else f, l) for ln in lines for f, l in L_CALL.findall(ln)]

def check(seg, R, L):
    a, b = labels_r(R), labels_l(L)
    if a != b:
        for i, (x, y) in enumerate(zip(a + [None] * len(b), b + [None] * len(a))):
            if x != y:
                raise SystemExit(f'ctxload: segment {seg}: dumper prints {x}, loader reads {y} (item {i})')

# ---------------------------------------------------------------- enums of the translation
def block(src, i):
    """Index just past the brace that closes the one opened before src[i]."""
    dep = 1
    while dep:
        dep += {'{': 1, '}': -1}.get(src[i], 0); i += 1
    return i

def scan(src, path, out):
    """Collect pub fieldless enums at `path`, descending into inline `pub mod x { ... }`."""
    for m in re.finditer(r'^\s*pub mod (\w+)\s*\{', src, re.M):
        j = block(src, m.end())
        scan(src[m.end():j - 1], path + '::' + m.group(1), out)
        src = src[:m.start()] + ' ' * (j - m.start()) + src[j:]
    for m in re.finditer(r'^\s*mod (\w+)\s*\{', src, re.M):  # private inline module
        j = block(src, m.end())
        src = src[:m.start()] + ' ' * (j - m.start()) + src[j:]
    for m in re.finditer(r'^\s*pub enum (\w+)\s*\{', src, re.M):
        body = src[m.end():block(src, m.end()) - 1]
        body = re.sub(r'#\[[^\]]*\]', '', body)
        names = []
        for v in (v.strip() for v in D.split_top(body)):
            if not v:
                continue
            vm = re.fullmatch(r'([A-Za-z_]\w*)\s*(?:=\s*[^,]+)?', v)
            if not vm:
                names = None
                break
            names.append(vm.group(1))
        if names:
            out.append((path + '::' + m.group(1), names))

def rust_enums():
    """Every pub fieldless enum in a pub module of crates/engine/src: (rust path, [variants])."""
    out = []
    private_files = set()
    for d, _, fs in os.walk(SRC):
        for f in fs:
            if f.endswith('.rs'):
                src = open(os.path.join(d, f)).read()
                for m in re.finditer(r'^\s*mod (\w+);', src, re.M):
                    base = os.path.join(d, os.path.splitext(f)[0] if f not in ('mod.rs', 'lib.rs') else '', m.group(1))
                    private_files.add(os.path.normpath(base + '.rs'))
    for d, _, fs in sorted(os.walk(SRC)):
        for f in sorted(fs):
            if not f.endswith('.rs'):
                continue
            fp = os.path.normpath(os.path.join(d, f))
            if fp in private_files:
                continue
            rel = os.path.relpath(fp, SRC)[:-3].split(os.sep)
            if rel[-1] in ('mod', 'lib'):
                rel = rel[:-1]
            scan(D.strip(open(fp).read()), 'swmm6_rust::' + '::'.join(rel), out)
    return out

# ---------------------------------------------------------------- the walk (ctxdump.main order)
def walk():
    """Returns (segments of the Rust dumper re-derived, loader lines) for dump() and engdump()."""
    R, L = [], []
    def hand(rlines, llines):
        R.extend(rlines); L.extend(llines)
    def gen(seg, fn):
        r0, l0 = len(R), len(L)
        fn()
        check(seg, R[r0:], L[l0:])
    C = []  # C++ side of ctxdump.emit, discarded
    for n in D.NAMES:
        hand([f'    vs(o, "{n}", c.{n}.names());'], [f'    load_names(&mut c.{n}, r.vs("{n}"));'])
    ctx_fields = {n: (t, a) for n, t, a in D.fields('core/SimulationContext.hpp', 'SimulationContext')}
    for n in D.CTX_SCALARS:
        t, a = ctx_fields[n]; k = D.kind(t)
        def one(n=n, a=a, k=k):
            if a:
                R.append('    for q in 0..12 {'); D.emit(f'c.{n}[q]', f'c.{n}[q]', n + '[]', k, C, R); R.append('    }')
                L.append('    for q in 0..12 {'); load(f'c.{n}[q]', n + '[]', k, L); L.append('    }')
            else:
                D.emit(f'c.{n}', f'c.{n}', n, k, C, R)
                load(f'c.{n}', n, k, L)
        gen(n, one)
    # every Table member and the by-name index (loaded as dumped: folded keys, as the
    # translation stores them)
    gen('tables', lambda: (D.dump_struct('data/TableData.hpp', 'TableData', 'c.tables', 'c.tables', 'tables', C, R),
                           load_struct('data/TableData.hpp', 'TableData', 'c.tables', 'tables', L)))
    hand(['    for h in &c.unknown_section_headers { ps(o, "unknown_header", &h.0); pi(o, "unknown_header.line", h.1); }'],
         ['    c.unknown_section_headers.clear(); while r.peek("unknown_header") { let a = r.ps("unknown_header"); let b = FromI64::from_i64(r.pi("unknown_header.line")); c.unknown_section_headers.push((a, b)); }'])
    hand(['    for (k, v) in &c.passthrough_sections { ps(o, "passthrough", k); vs(o, "passthrough.lines", v); }'],
         ['    c.passthrough_sections.clear(); while r.peek("passthrough") { let k = r.ps("passthrough"); let v = r.vs("passthrough.lines"); c.passthrough_sections.insert(k, v); }'])
    hand(['    vi(o, "node_subtypes.subtype_row", &c.node_subtypes.subtype_row.iter().map(|x| *x as i64).collect::<Vec<i64>>());',
          '    vi(o, "link_subtypes.subtype_row", &c.link_subtypes.subtype_row.iter().map(|x| *x as i64).collect::<Vec<i64>>());'],
         ['    set_vi(&mut c.node_subtypes.subtype_row, r.vi("node_subtypes.subtype_row"));',
          '    set_vi(&mut c.link_subtypes.subtype_row, r.vi("link_subtypes.subtype_row"));'])
    hand(['    pi(o, "events.n", c.events.len() as i64); for e in &c.events { pf(o, "event.start", e.start); pf(o, "event.end", e.end); }'],
         ['    { let n = r.count("events.n"); SeqLen::seq_len(&mut c.events, n); for e in c.events.iter_mut() { e.start = r.pf("event.start"); e.end = r.pf("event.end"); } }'])
    hand(['    pi(o, "pending_link_nodes.n", c.pending_link_nodes.len() as i64); for p in &c.pending_link_nodes { pi(o, "pending_link_node", p.0 as i64); ps(o, "pending_link_node.n1", &(p.1).0); ps(o, "pending_link_node.n2", &(p.1).1); }'],
         ['    { let n = r.count("pending_link_nodes.n"); SeqLen::seq_len(&mut c.pending_link_nodes, n); for p in c.pending_link_nodes.iter_mut() { p.0 = FromI64::from_i64(r.pi("pending_link_node")); (p.1).0 = r.ps("pending_link_node.n1"); (p.1).1 = r.ps("pending_link_node.n2"); } }'])
    hand(['    for p in &c.plugin_specs { ps(o, "plugin.path", &p.path); vs(o, "plugin.init_args", &p.init_args); }'],
         ['    c.plugin_specs.clear(); while r.peek("plugin.path") { let mut p = PluginSpec::default(); p.path = r.ps("plugin.path"); p.init_args = r.vs("plugin.init_args"); c.plugin_specs.push(p); }'])
    hand(['    for p in &c.process_component_specs { ps(o, "pcomp.id", &p.id); ps(o, "pcomp.config_path", &p.config_path); for kv in &p.args { ps(o, "pcomp.arg", &kv.0); ps(o, "pcomp.arg.value", &kv.1); } ps(o, "pcomp.resolved_config_path", &p.resolved_config_path); }'],
         ['    c.process_component_specs.clear(); while r.peek("pcomp.id") { let mut p = ProcessComponentSpec::default(); p.id = r.ps("pcomp.id"); p.config_path = r.ps("pcomp.config_path");',
          '    while r.peek("pcomp.arg") { let k = r.ps("pcomp.arg"); let v = r.ps("pcomp.arg.value"); p.args.push((k, v)); } p.resolved_config_path = r.ps("pcomp.resolved_config_path"); c.process_component_specs.push(p); }'])
    hand(['    for d in c.user_flags.all_defs() { ps(o, "uflag_def.name", &d.name); pi(o, "uflag_def.type", d.type_ as i64); ps(o, "uflag_def.description", &d.description); }'],
         ['    c.user_flags.clear(); while r.peek("uflag_def.name") { let name = r.ps("uflag_def.name"); let type_ = FromI64::from_i64(r.pi("uflag_def.type")); let description = r.ps("uflag_def.description"); c.user_flags.define(UserFlagDef { name, type_, description }); }'])
    hand(['    { let mut e: Vec<_> = c.user_flags.all_values().iter().collect(); e.sort_by(|a, b| a.0.cmp(b.0)); for (k, v) in e { ps(o, "uflag_value", k); match v { UserFlagValue::Bool(b) => { pi(o, "uflag_value.kind", 0); pi(o, "uflag_value.v", *b as i64); } UserFlagValue::Int(i) => { pi(o, "uflag_value.kind", 1); pi(o, "uflag_value.v", *i as i64); } UserFlagValue::Double(d) => { pi(o, "uflag_value.kind", 2); pf(o, "uflag_value.v", *d); } UserFlagValue::String(s) => { pi(o, "uflag_value.kind", 3); ps(o, "uflag_value.v", s); } } } }'],
         ['    while r.peek("uflag_value") { let key = r.ps("uflag_value"); let v = match r.pi("uflag_value.kind") {',
          '    0 => UserFlagValue::Bool(FromI64::from_i64(r.pi("uflag_value.v"))), 1 => UserFlagValue::Int(FromI64::from_i64(r.pi("uflag_value.v"))),',
          '    2 => UserFlagValue::Double(r.pf("uflag_value.v")), 3 => UserFlagValue::String(r.ps("uflag_value.v")), k => r.fail(&format!("user flag value kind {}", k)) };',
          '    let (a, b, f) = split_flag_key(&key).unwrap_or_else(|| r.fail("user flag key without two colons")); c.user_flags.set(a, b, f, v); }'])
    # def_index_ is what define() builds from the definitions: checked, not loaded
    hand(['    { let mut t: Vec<_> = c.user_flags.def_index().iter().collect(); t.sort(); for (k, v) in t { ps(o, "uflag_def_index", k); pi(o, "uflag_def_index.v", *v as i64); } }'],
         ['    { let mut t: Vec<(Bytes, usize)> = Vec::new(); while r.peek("uflag_def_index") { let k = r.ps("uflag_def_index"); let v: usize = FromI64::from_i64(r.pi("uflag_def_index.v")); t.push((k, v)); }',
          '    let mut got: Vec<(Bytes, usize)> = c.user_flags.def_index().iter().map(|(k, v)| (k.clone(), *v)).collect(); got.sort(); if got != t { r.fail("uflag_def_index differs from the index UserFlags::define builds"); } }'])
    hand(['    for h in &c.files.hotstart_saves { ps(o, "hotstart_save.path", &h.path.original); ps(o, "hotstart_save.path.abs", &h.path.absolute); pf(o, "hotstart_save.datetime", h.datetime); pi(o, "hotstart_save.format", h.format as i64); }'],
         ['    c.files.hotstart_saves.clear(); while r.peek("hotstart_save.path") { let mut h = HotstartSaveEntry::default(); h.path.original = r.ps("hotstart_save.path"); h.path.absolute = r.ps("hotstart_save.path.abs");',
          '    h.datetime = r.pf("hotstart_save.datetime"); h.format = FromI64::from_i64(r.pi("hotstart_save.format")); c.files.hotstart_saves.push(h); }'])
    for header, struct, pc, pr, label in D.STORES:
        gen(label, lambda: (D.dump_struct(header, struct, pc, pr, label, C, R), load_struct(header, struct, pr, label, L)))
    def tt():
        R.append('    pi(o, "transect_tables.n", c.transect_tables.len() as i64);'); R.append('    for el in &c.transect_tables {')
        D.dump_struct('hydraulics/Transect.hpp', 'TransectData', 'el', 'el', 'transect_tables[]', C, R)
        R.append('    }')
        L.append('    { let n = r.count("transect_tables.n"); SeqLen::seq_len(&mut c.transect_tables, n); for el in c.transect_tables.iter_mut() {')
        load_struct('hydraulics/Transect.hpp', 'TransectData', 'el', 'transect_tables[]', L)
        L.append('    } }')
    gen('transect_tables', tt)
    hand(['    pi(o, "met.n", c.snowpacks.meteorology.records().len() as i64); for r in c.snowpacks.meteorology.records().iter() { pf(o, "met.date", r.date); for v in r.values.iter() { pf(o, "met.value", *v); } }'],
         ['    { let n = r.count("met.n"); let mut recs: Vec<HourlyMetRecord> = Vec::with_capacity(n); for _ in 0..n { let mut m = HourlyMetRecord::default(); m.date = r.pf("met.date");',
          '    for v in m.values.iter_mut() { *v = r.pf("met.value"); } recs.push(m); } if n > 0 { set_met_records(r, &mut c.snowpacks.meteorology, recs); } }'])
    hand(['    pi(o, "species.n", c.species_registry.count() as i64); pi(o, "species.pollutants", c.species_registry.pollutant_count() as i64); for i in 0..c.species_registry.count() { ps(o, "species.name", c.species_registry.name(i)); pi(o, "species.kind", c.species_registry.kind(i) as i64); ps(o, "species.units", c.species_registry.units(i)); }'],
         ['    { let n = r.count("species.n"); let np = r.pi("species.pollutants"); c.species_registry.clear(); for _ in 0..n { let name = r.ps("species.name"); let kind = FromI64::from_i64(r.pi("species.kind")); let units = r.ps("species.units");',
          '    if c.species_registry.add(name, kind, units) < 0 { r.fail("duplicate species name"); } } if c.species_registry.pollutant_count() as i64 != np { r.fail("species.pollutants differs from the pollutant kinds listed"); } }'])
    for n in D.CTX_MEMBERS:
        t, a = ctx_fields[n]
        gen(n, lambda n=n, t=t, a=a: (D.ctx_member(n, t, a, C, R), load_ctx_member(n, t, a, L)))
    # ---- engdump (ctxdump.engine_dumper order)
    RE, LE = [], []
    for n in D.ENGINE_MEMBERS:
        t = D.member_type(D.ENGINE_HEADER, n); k = D.kind(t)
        r0, l0 = len(RE), len(LE)
        D.emit(f'(e.*get(T_{n}()))', f'e.{n.rstrip("_")}()', f'engine.{n}', k, C, RE)
        load(f'(*m.{n})', f'engine.{n}', k, LE)
        check(f'engine.{n}', RE[r0:], LE[l0:])
    RE += ['    let pf = e.plugin_factory();',
           '    pi(o, "engine.plugins.input.n", pf.input_plugins().len() as i64);',
           '    for p in pf.input_plugins() { pi(o, "engine.plugins.input.state", p.state() as i64); pi(o, "engine.plugins.input.default", 1); vs(o, "engine.plugins.input.skipped", &p.skipped_sections()); ps(o, "engine.plugins.input.last_error", p.last_error_message()); }',
           '    pi(o, "engine.plugins.output.n", pf.output_plugins().len() as i64);',
           '    for p in pf.output_plugins() { pi(o, "engine.plugins.output.state", p.state() as i64); pi(o, "engine.plugins.output.default", 1); ps(o, "engine.plugins.output.path", p.out_path()); }',
           '    pi(o, "engine.plugins.report.n", pf.report_plugins().len() as i64);',
           '    for p in pf.report_plugins() { pi(o, "engine.plugins.report.state", p.state() as i64); pi(o, "engine.plugins.report.default", 1); ps(o, "engine.plugins.report.path", p.rpt_path()); }',
           '    pi(o, "engine.plugins.state_io.n", pf.state_io_plugins().len() as i64);',
           '    for p in pf.state_io_plugins() { pi(o, "engine.plugins.state_io.state", p.state() as i64); pi(o, "engine.plugins.state_io.default", 1); }',
           '    pi(o, "engine.binary_output_plugin_", e.binary_output_plugin().map(|k| k as i64).unwrap_or(-1));']
    LE += ['    let pf = &mut *m.plugins_;',
           '    if !(pf.input_plugins().is_empty() && pf.output_plugins().is_empty() && pf.report_plugins().is_empty() && pf.state_io_plugins().is_empty()) { r.fail("load_engine needs a fresh SwmmEngine"); }',
           '    for _ in 0..r.count("engine.plugins.input.n") { let st = FromI64::from_i64(r.pi("engine.plugins.input.state")); r.builtin("engine.plugins.input.default");',
           '    let sk = r.vs("engine.plugins.input.skipped"); let le = r.ps("engine.plugins.input.last_error"); let mut p = DefaultInputPlugin::new(); p.gate_set(st, le, sk); pf.add_input_plugin(p); }',
           '    for _ in 0..r.count("engine.plugins.output.n") { let st = FromI64::from_i64(r.pi("engine.plugins.output.state")); r.builtin("engine.plugins.output.default");',
           '    let mut p = DefaultOutputPlugin::new(r.ps("engine.plugins.output.path")); p.gate_set_state(st); pf.add_output_plugin(p); }',
           '    for _ in 0..r.count("engine.plugins.report.n") { let st = FromI64::from_i64(r.pi("engine.plugins.report.state")); r.builtin("engine.plugins.report.default");',
           '    let mut p = DefaultReportPlugin::new(r.ps("engine.plugins.report.path")); p.gate_set_state(st); pf.add_report_plugin(p); }',
           '    for _ in 0..r.count("engine.plugins.state_io.n") { let st = FromI64::from_i64(r.pi("engine.plugins.state_io.state")); r.builtin("engine.plugins.state_io.default");',
           '    let mut p = DefaultStateIOPlugin::new(); p.gate_set_state(st); pf.add_state_io_plugin(p); }',
           '    { let k = r.pi("engine.binary_output_plugin_"); *m.binary_output_plugin_ = if k < 0 { None } else { Some(k as usize) }; }']
    return R, L, RE, LE

# ---------------------------------------------------------------- the hydrology loader
HOUT = f'{ROOT}/crates/engine/examples/common/hydroload.rs'

def hydro_walk():
    """(Rust dumper lines re-derived, loader lines) for hydrodump (ctxdump.hydro_dumper order).
    The loader fills the engine through SwmmEngine::gate_hydrology() (`h`)."""
    saved = dict(D.STRUCTS)
    D.STRUCTS.update(D.HYDRO_STRUCTS)
    R, L, C = [], [], []
    try:
        for header, struct, pc, pr, label in D.HYDRO_ENGINE:
            r0, l0 = len(R), len(L)
            D.dump_struct(header, struct, pc, pr, label, C, R)
            load_struct(header, struct, f'(*h.{label})', label, L)
            check(label, R[r0:], L[l0:])
        nstruct = len(R)
        # the hand-written tail of hydro_dumper, mirrored
        L.append('    for u in h.rdii_.uh_params.iter_mut() { for m in 0..12 { for k in 0..3 { u.r[m][k] = r.pf("rdii_.uh_params[].r"); u.tPeak[m][k] = r.pf("rdii_.uh_params[].tPeak"); u.tBase[m][k] = r.pf("rdii_.uh_params[].tBase"); u.iaMax[m][k] = r.pf("rdii_.uh_params[].iaMax"); u.iaRecov[m][k] = r.pf("rdii_.uh_params[].iaRecov"); u.iaInit[m][k] = r.pf("rdii_.uh_params[].iaInit"); } } }')
        # the name map is initHydraulics' (identical on both sides): checked, not loaded
        L.append('    { let mut t: Vec<(Bytes, i32)> = Vec::new(); while r.peek("rdii_.uh_name_to_idx_.key") { let k = r.ps("rdii_.uh_name_to_idx_.key"); let v = r.pi("rdii_.uh_name_to_idx_.val") as i32; t.push((k, v)); } let mut m = h.rdii_.uh_name_to_idx_.clone(); m.sort(); if m != t { r.fail("rdii_.uh_name_to_idx_ differs"); } }')
        # the report files themselves are not reopened: a loaded unit's file is a placeholder
        # that counts as failed, so the translation writes nothing more to it
        L.append('    { let n = r.pi("lid_.rpt_files_.size"); h.lid_.rpt_files_.clear(); for _ in 0..n { let a = r.pi("lid_.rpt_files_.key.first") as i32; let b = r.pi("lid_.rpt_files_.key.second") as i32; let mut f = Box::new(swmm6_rust::hydrology::lid::LidRptFile::default()); f.failed = true; h.lid_.rpt_files_.insert((a, b), f); } }')
        r0, l0 = len(R), len(L)
        D.emit('e.gage_states_', 'e.gage_states()', 'gage_states_', D.kind('std::vector<gage::GageState>'), C, R)
        load('(*h.gage_states_)', 'gage_states_', D.kind('std::vector<gage::GageState>'), L)
        check('gage_states_', R[r0:], L[l0:])
        L.append('    *h.last_melt_doy_ = r.pi("last_melt_doy_") as i32;')
        for k, n in enumerate(D.HYDRO_VECTORS):
            L.append(f'    set_vf(&mut *h.vectors[{k}], r.vf("{n}"));')
        L.append('    set_vi(&mut *h.gw_q_node_, r.vi("gw_q_node_"));')
        L.append('    load_climate_file(r, &mut *h.climate_file_);')
    finally:
        D.STRUCTS.clear(); D.STRUCTS.update(saved)
    return R[:nstruct], R[nstruct:], L

HPRELUDE = r'''// generated by tools/ctxload.py from tools/ctxdump.py's field lists: reads what
// examples/common/hydrodump.rs (or gate/diff/hydrodump.hpp) prints back into a started
// SwmmEngine's hydrology solvers, through SwmmEngine::gate_hydrology() (`gate` feature).
// The RDII unit-hydrograph name map and the LID report files are checked, not loaded; the
// climate file reader is positioned where the dump says it was.
#![allow(dead_code, unused_parens, unused_mut, unused_variables, unused_assignments, unused_imports, clippy::all)]
use super::ctxload::*;
use swmm6_rust::core::engine::SwmmEngine;
use swmm6_rust::core::options::Bytes;
use swmm6_rust::hydrology::climate_file::ClimateFileReader;

fn load_climate_file(r: &mut Rd, cf: &mut ClimateFileReader) {
    let open = r.pi("climate_file_.open");
    let pos = r.pi("climate_file_.pos");
    if (open != 0) != cf.file_.is_some() {
        r.fail("climate_file_.open differs (the reader is not reopened)");
    }
    if let Some(f) = cf.file_.as_mut() {
        if f.fseek(pos) != 0 {
            r.fail("climate_file_: seek failed");
        }
    }
    cf.format_ = FromI64::from_i64(r.pi("climate_file_.format_"));
    cf.unit_system_ = r.pi("climate_file_.unit_system_") as i32;
    for v in 0..4 {
        for d in 0..32 {
            cf.file_data_[v][d] = r.pf("climate_file_.file_data_");
        }
    }
    cf.buf_year_ = r.pi("climate_file_.buf_year_") as i32;
    cf.buf_month_ = r.pi("climate_file_.buf_month_") as i32;
    cf.day_offset_ = r.pf("climate_file_.day_offset_");
    for i in 0..4 {
        cf.field_pos_[i] = r.pi("climate_file_.field_pos_") as i32;
    }
    cf.date_field_pos_ = r.pi("climate_file_.date_field_pos_") as i32;
    cf.wind_type_ = r.pi("climate_file_.wind_type_") as i32;
    cf.temp_units_ = r.pi("climate_file_.temp_units_") as i32;
    cf.saved_line_ = r.ps("climate_file_.saved_line_");
    cf.has_saved_line_ = r.pi("climate_file_.has_saved_line_") != 0;
}
'''
# ---------------------------------------------------------------- router/solver and quality members
# (the routing-step cross test, examples/routestep.rs). The Rust side rebuilds every solver with the
# translated open + init legs, then loads upstream's dump over it, so the members ctxdump.py prints
# through a custom emitter (RT_CUSTOM: tables of the XSectGroups, the FV backend, the RDII UH
# parameters and name map, the control expressions) are consumed, not loaded: they are what the
# init legs build, which the inithyd gate holds identical. The one exception is the XSectGroups
# scratch buffers (buf_d, buf_r, buf_r2), which every step overwrites and which are loaded.
# Checked by the round trip in `ctxgate.sh route` (load, dump again, cmp, before any step).
RT_GATE_PATHS = {'e.router()': 'g.router_', 'e.control_engine()': 'g.controls_', 'e.inflow_solver()': 'g.inflow_',
                 'e.rdii_solver()': 'g.rdii_', 'e.exfil_solver()': 'g.exfil_', 'e.inlet_solver()': 'g.inlet_',
                 'e.structure_solver()': 'g.hydstruct_', 'e.interface_manager()': 'g.iface_'}

def load_rt_struct(header, struct, path_r, label, L):
    """Mirror of ctxdump.rt_dump_struct."""
    for name, ty, arr in D.fields(header, struct):
        pr, lb = f'{path_r}.{D.rname(name)}', f'{label}.{name}'
        if (struct, name) == ('Router', 'groups_'):
            # the group layout (shape, count) must be what the init legs built; the per-link
            # geometry arrays are loaded (upstream's initialize can change them after
            # initHydraulics, e.g. a force main's s_bot), the transect table views are not
            L.append(f'    {{ let n = r.count("{lb}.n"); if n != ({pr}).num_groups() as usize {{ r.fail("XSectGroups layout differs from the init legs"); }} for gi in 0..n {{ let grp = ({pr}).gate_group_mut(gi as i32);')
            L.append(f'    if r.pi("{lb}[].shape") != grp.shape as i64 || r.pi("{lb}[].count") != grp.count as i64 {{ r.fail("XSectGroups group differs from the init legs"); }}')
            L.append(f'    set_vi(&mut grp.link_idx, r.vi("{lb}[].link_idx"));')
            for f in ['y_full', 'a_full', 'r_full', 's_full', 'w_max', 'inv_y_full', 'y_bot', 'a_bot', 's_bot', 'r_bot']:
                L.append(f'    set_vf(&mut grp.{f}, r.vf("{lb}[].{f}"));')
            L.append(f'    while !r.peek("{lb}[].buf_d") {{ r.skip_item(); }}')
            L.append(f'    *grp.buf_d.borrow_mut() = r.vf("{lb}[].buf_d"); *grp.buf_r.borrow_mut() = r.vf("{lb}[].buf_r"); *grp.buf_r2.borrow_mut() = r.vf("{lb}[].buf_r2"); }}')
            L.append(f'    r.sk("{lb}.packed_groups_"); r.sk("{lb}.packed_count_"); r.sk("{lb}.mask_active_"); }}')
            continue
        if (struct, name) in D.RT_CUSTOM:
            L.append(f'    r.skip_prefix("{lb}");')
            continue
        k = D.kind(ty)
        if arr:
            if k and k[0] in 'fi':
                L.append(f'    for q in 0..({pr}).len() {{')
                load(f'{pr}[q]', lb + '[]', k, L)
                L.append('    }')
            else:
                skip(lb, L)
            continue
        t = ty.replace('const ', '').strip()
        em = re.fullmatch(r'std::vector<(?:\w+::)*(\w+)>', t)
        if k is None and em and em.group(1) in D.RT_ELEM_STRUCTS:
            L.append(f'    {{ let n = r.count("{lb}.n"); SeqLen::seq_len(&mut {pr}, n); for el in ({pr}).iter_mut() {{')
            load_rt_struct(D.RT_ELEM_STRUCTS[em.group(1)], em.group(1), '(*el)', lb + '[]', L)
            L.append('    } }')
            continue
        if k is not None and k[0] == 'struct':
            load_rt_struct(k[1], k[2], pr, lb, L)
            continue
        if k is None:
            skip(lb, L)
            continue
        load(pr, lb, k, L)
    for name, flag in D.RT_EXTRA.get(struct, []):
        L.append(f'    r.pi("{label}.{name}");')

def walk_rt():
    saved = dict(D.STRUCTS)
    D.STRUCTS.update(D.RT_STRUCTS)
    L = []
    for header, struct, pc, pr, label in D.RT_MEMBERS:
        load_rt_struct(header, struct, RT_GATE_PATHS[pr], label, L)
    L.append('    *g.next_event_ = FromI64::from_i64(r.pi("engine.next_event_")); *g.between_events_ = FromI64::from_i64(r.pi("engine.between_events_"));')
    D.STRUCTS.clear(); D.STRUCTS.update(saved)
    Q = ['    *q.hotstart_deferred_ = FromI64::from_i64(r.pi("engine.hotstart_deferred_"));']
    for header, cls, pc, pr, label in D.QUALITY_MEMBERS:
        load_struct(header, cls, 'q.' + pc.split('.')[1], label, Q)
    return L, Q

def fn_body(lines, header):
    """The lines of a generated Rust fn after its three opening lines, up to `std::mem::take(o)`."""
    i = lines.index(header) + 3
    j = lines.index('    std::mem::take(o)', i)
    return lines[i:j]

PRELUDE = r'''// generated by tools/ctxload.py from tools/ctxdump.py's field lists: reads what
// examples/common/ctxdump.rs (or gate/diff/ctxdump.hpp) prints back into a fresh
// SimulationContext and SwmmEngine. Dumping them again reproduces the input byte for byte.
// Each line's label is checked; doubles are loaded from their bits (NaN payload and sign
// kept), strings as bytes, vector and map lengths from the dump. Members the dump marks
// "skip" stay at their defaults. The hourly meteorology records and the engine members are
// private in the translation: loading them needs the `gate` feature's hooks.
#![allow(dead_code, unused_parens, unused_mut, unused_variables, unused_assignments, unused_imports, clippy::all)]
use std::collections::{BTreeMap, HashMap};
use swmm6_rust::core::context::{HotstartSaveEntry, PluginSpec, ProcessComponentSpec, SimulationContext};
use swmm6_rust::core::engine::SwmmEngine;
use swmm6_rust::core::options::{Bytes, FilePathPair};
use swmm6_rust::core::plugins::{DefaultInputPlugin, DefaultOutputPlugin, DefaultReportPlugin, DefaultStateIOPlugin};
use swmm6_rust::core::user_flags::{UserFlagDef, UserFlagValue};
use swmm6_rust::data::name_index::NameIndex;
use swmm6_rust::data::table::Table;
use swmm6_rust::hydrology::hourly_meteorology::{HourlyMeteorology, HourlyMetRecord};

/// A cursor over a dump. Every read names the label it expects.
pub struct Rd<'a> {
    b: &'a [u8],
    p: usize,
}

impl<'a> Rd<'a> {
    pub fn new(b: &'a [u8]) -> Self {
        Rd { b, p: 0 }
    }
    pub fn fail(&self, what: &str) -> ! {
        let line = self.b[..self.p].iter().filter(|&&c| c == b'\n').count() + 1;
        let end = self.b[self.p..].iter().position(|&c| c == b'\n').map(|k| self.p + k).unwrap_or(self.b.len());
        let shown = String::from_utf8_lossy(&self.b[self.p..end.min(self.p + 120)]);
        panic!("ctxload: line {}: {} (at \"{}\")", line, what, shown)
    }
    /// The next line is `label ...`.
    pub fn peek(&self, l: &str) -> bool {
        let r = &self.b[self.p..];
        r.len() > l.len() && r.starts_with(l.as_bytes()) && r[l.len()] == b' '
    }
    /// The next line is `label none`.
    pub fn at_none(&self, l: &str) -> bool {
        let r = &self.b[self.p..];
        r.starts_with(l.as_bytes()) && r[l.len()..].starts_with(b" none\n")
    }
    fn label(&mut self, l: &str) {
        if !self.peek(l) {
            self.fail(&format!("expected label {}", l));
        }
        self.p += l.len() + 1;
    }
    fn byte(&mut self, c: u8) {
        if self.p >= self.b.len() || self.b[self.p] != c {
            self.fail(&format!("expected {:?}", c as char));
        }
        self.p += 1;
    }
    fn token(&mut self) -> &'a [u8] {
        let s = self.p;
        while self.p < self.b.len() && !matches!(self.b[self.p], b' ' | b'\n' | b':') {
            self.p += 1;
        }
        &self.b[s..self.p]
    }
    fn dec(&mut self) -> i64 {
        let t = self.token();
        match std::str::from_utf8(t).ok().and_then(|s| s.parse::<i64>().ok()) {
            Some(v) => v,
            None => self.fail("expected a decimal integer"),
        }
    }
    fn hex(&mut self) -> f64 {
        let t = self.token();
        match std::str::from_utf8(t).ok().filter(|s| s.len() == 16).and_then(|s| u64::from_str_radix(s, 16).ok()) {
            Some(u) => f64::from_bits(u),
            None => self.fail("expected 16 hex digits"),
        }
    }
    fn len(&mut self) -> usize {
        let n = self.dec();
        if n < 0 {
            self.fail("negative length");
        }
        n as usize
    }
    fn bytes(&mut self) -> Bytes {
        let n = self.len();
        self.byte(b':');
        if self.p + n > self.b.len() {
            self.fail("string runs past the end");
        }
        let s = self.b[self.p..self.p + n].to_vec();
        self.p += n;
        s
    }
    /// `pf`: label, 16 hex digits of the bits.
    pub fn pf(&mut self, l: &str) -> f64 {
        self.label(l);
        let v = self.hex();
        self.byte(b'\n');
        v
    }
    /// `pi`: label, decimal.
    pub fn pi(&mut self, l: &str) -> i64 {
        self.label(l);
        let v = self.dec();
        self.byte(b'\n');
        v
    }
    /// A count printed by `pi`.
    pub fn count(&mut self, l: &str) -> usize {
        let n = self.pi(l);
        if n < 0 {
            self.fail("negative count");
        }
        n as usize
    }
    /// `ps`: label, length, ':', the bytes.
    pub fn ps(&mut self, l: &str) -> Bytes {
        self.label(l);
        let s = self.bytes();
        self.byte(b'\n');
        s
    }
    fn vhead(&mut self, l: &str) -> usize {
        self.label(l);
        if !self.b[self.p..].starts_with(b"n=") {
            self.fail("expected n=");
        }
        self.p += 2;
        let n = self.len();
        self.byte(b'\n');
        n
    }
    fn vtail(&mut self, n: usize) {
        if n > 0 {
            self.byte(b'\n');
        }
    }
    pub fn vf(&mut self, l: &str) -> Vec<f64> {
        let n = self.vhead(l);
        let mut v = Vec::with_capacity(n);
        for _ in 0..n {
            self.byte(b' ');
            v.push(self.hex());
        }
        self.vtail(n);
        v
    }
    pub fn vi(&mut self, l: &str) -> Vec<i64> {
        let n = self.vhead(l);
        let mut v = Vec::with_capacity(n);
        for _ in 0..n {
            self.byte(b' ');
            v.push(self.dec());
        }
        self.vtail(n);
        v
    }
    pub fn vs(&mut self, l: &str) -> Vec<Bytes> {
        let n = self.vhead(l);
        let mut v = Vec::with_capacity(n);
        for _ in 0..n {
            self.byte(b' ');
            v.push(self.bytes());
        }
        self.vtail(n);
        v
    }
    /// `sk`: a member the dumper does not print. It stays at its default.
    pub fn sk(&mut self, l: &str) {
        self.label("skip");
        let t = self.token();
        if t != l.as_bytes() {
            self.fail(&format!("expected skip {}", l));
        }
        self.byte(b'\n');
    }
    /// `pn`: an empty std::optional.
    pub fn pn(&mut self, l: &str) {
        if !self.at_none(l) {
            self.fail(&format!("expected {} none", l));
        }
        self.p += l.len() + 6;
    }
    /// A plugin of a built-in type (the only kind the translation has).
    pub fn builtin(&mut self, l: &str) {
        if self.pi(l) != 1 {
            self.fail("a plugin that is not built in cannot be loaded");
        }
    }
    /// A line the harness printed before the dump (`open_rc 0`), returned whole.
    pub fn header_line(&mut self, l: &str) -> Option<&'a [u8]> {
        if !self.peek(l) {
            return None;
        }
        let s = self.p;
        while self.b[self.p] != b'\n' {
            self.p += 1;
        }
        self.p += 1;
        Some(&self.b[s..self.p])
    }
    /// Consume one dumped item: its label line and any continuation lines (vector values,
    /// which start with a space).
    pub fn skip_item(&mut self) {
        if self.p >= self.b.len() {
            self.fail("unexpected end of the dump");
        }
        loop {
            while self.p < self.b.len() && self.b[self.p] != b'\n' {
                self.p += 1;
            }
            self.p += 1;
            if self.p >= self.b.len() || self.b[self.p] != b' ' {
                break;
            }
        }
    }
    /// Consume every item whose label is `l` or starts with `l.` / `l[`, and `skip` lines of them.
    pub fn skip_prefix(&mut self, l: &str) {
        loop {
            let mut r = &self.b[self.p..];
            if r.starts_with(b"skip ") {
                r = &r[5..];
            }
            if r.len() > l.len() && r.starts_with(l.as_bytes()) && matches!(r[l.len()], b'.' | b'[' | b' ') {
                self.skip_item();
            } else {
                break;
            }
        }
    }
    pub fn at_end(&self) -> bool {
        self.p == self.b.len()
    }
    pub fn finish(&self) {
        if !self.at_end() {
            self.fail("input continues after the dump");
        }
    }
}

/// Upstream's integer value of an integer, bool or enum member (as the dumpers print it).
pub trait FromI64: Sized {
    fn from_i64(v: i64) -> Self;
}
macro_rules! from_i64_int {
    ($($t:ty),*) => {$(
        impl FromI64 for $t {
            fn from_i64(v: i64) -> Self {
                let x = v as $t;
                if x as i64 != v {
                    panic!("ctxload: {} out of range for {}", v, stringify!($t));
                }
                x
            }
        }
    )*};
}
from_i64_int!(i8, i16, i32, i64, isize, u8, u16, u32, usize);
impl FromI64 for u64 {
    // printed as (long long), so a u64 above i64::MAX arrives negative
    fn from_i64(v: i64) -> Self {
        v as u64
    }
}
impl FromI64 for bool {
    fn from_i64(v: i64) -> Self {
        match v {
            0 => false,
            1 => true,
            _ => panic!("ctxload: {} is not a bool", v),
        }
    }
}
fn bad_enum(name: &str, v: i64) -> ! {
    panic!("ctxload: {} has no variant with value {}", name, v)
}

/// A double member from its bits; a float member from the double it was widened to.
pub trait FromF64 {
    fn from_f64(v: f64) -> Self;
}
impl FromF64 for f64 {
    fn from_f64(v: f64) -> Self {
        v
    }
}
impl FromF64 for f32 {
    fn from_f64(v: f64) -> Self {
        if v.is_nan() {
            // keep sign and the payload's top bits, as (double)f widened them
            let u = v.to_bits();
            f32::from_bits((((u >> 63) as u32) << 31) | 0x7f80_0000 | ((u >> 29) as u32 & 0x007f_ffff))
        } else {
            v as f32
        }
    }
}

/// Replace a sequence's elements (a Vec takes the dump's length; an array must match it).
pub trait SeqSet {
    type Elem;
    fn seq_set(&mut self, v: Vec<Self::Elem>);
}
impl<T> SeqSet for Vec<T> {
    type Elem = T;
    fn seq_set(&mut self, v: Vec<T>) {
        *self = v;
    }
}
impl<T, const N: usize> SeqSet for [T; N] {
    type Elem = T;
    fn seq_set(&mut self, v: Vec<T>) {
        if v.len() != N {
            panic!("ctxload: {} elements for an array of {}", v.len(), N);
        }
        for (d, s) in self.iter_mut().zip(v) {
            *d = s;
        }
    }
}
/// Size a sequence to the dump's count before its elements are filled in place.
pub trait SeqLen {
    fn seq_len(&mut self, n: usize);
}
impl<T: Default> SeqLen for Vec<T> {
    fn seq_len(&mut self, n: usize) {
        self.clear();
        self.resize_with(n, T::default);
    }
}
impl<T, const N: usize> SeqLen for [T; N] {
    fn seq_len(&mut self, n: usize) {
        if n != N {
            panic!("ctxload: {} elements for an array of {}", n, N);
        }
    }
}
fn seq_len_with<T>(v: &mut Vec<T>, n: usize, f: impl FnMut() -> T) {
    v.clear();
    v.resize_with(n, f);
}
fn set_vi<S: SeqSet>(d: &mut S, v: Vec<i64>)
where
    S::Elem: FromI64,
{
    d.seq_set(v.into_iter().map(S::Elem::from_i64).collect());
}
fn set_vf<S: SeqSet>(d: &mut S, v: Vec<f64>)
where
    S::Elem: FromF64,
{
    d.seq_set(v.into_iter().map(S::Elem::from_f64).collect());
}
fn set_va<const N: usize>(d: &mut Vec<[f64; N]>, flat: Vec<f64>) {
    if N == 0 || flat.len() % N != 0 {
        panic!("ctxload: {} values for rows of {}", flat.len(), N);
    }
    *d = flat.chunks_exact(N).map(|c| <[f64; N]>::try_from(c).unwrap()).collect();
}
fn set_fpp(d: &mut Vec<FilePathPair>, o: Vec<Bytes>, a: Vec<Bytes>) {
    if o.len() != a.len() {
        panic!("ctxload: {} paths, {} absolute paths", o.len(), a.len());
    }
    *d = o.into_iter().zip(a).map(|(original, absolute)| FilePathPair { absolute, original }).collect();
}
fn opt_new<T: Default>(_: &Option<T>) -> T {
    T::default()
}
/// std::map members (std::map iterates in key order, as BTreeMap does).
pub trait MapLoad {
    type K: Default;
    type V: Default;
    fn map_clear(&mut self);
    fn map_put(&mut self, k: Self::K, v: Self::V);
}
impl<K: Ord + Default, V: Default> MapLoad for BTreeMap<K, V> {
    type K = K;
    type V = V;
    fn map_clear(&mut self) {
        self.clear();
    }
    fn map_put(&mut self, k: K, v: V) {
        if self.insert(k, v).is_some() {
            panic!("ctxload: duplicate map key");
        }
    }
}
/// std::unordered_map members (dumped sorted by key; the order of a HashMap is not observable).
impl<K: Eq + std::hash::Hash + Default, V: Default> MapLoad for HashMap<K, V> {
    type K = K;
    type V = V;
    fn map_clear(&mut self) {
        self.clear();
    }
    fn map_put(&mut self, k: K, v: V) {
        if self.insert(k, v).is_some() {
            panic!("ctxload: duplicate map key");
        }
    }
}
fn map_k<M: MapLoad>(_: &M) -> M::K {
    M::K::default()
}
fn map_v<M: MapLoad>(_: &M) -> M::V {
    M::V::default()
}
/// A NameIndex from its names in index order (the map is rebuilt by add()).
fn load_names(ix: &mut NameIndex, v: Vec<Bytes>) {
    ix.clear();
    for s in &v {
        if let Err(e) = ix.add(s) {
            panic!("ctxload: {}", e);
        }
    }
}
/// UserFlags::make_key joined object type, name and flag with ':'; any split at two colons
/// rebuilds the same key.
fn split_flag_key(k: &[u8]) -> Option<(&[u8], &[u8], &[u8])> {
    let a = k.iter().position(|&c| c == b':')?;
    let b = a + 1 + k[a + 1..].iter().position(|&c| c == b':')?;
    Some((&k[..a], &k[a + 1..b], &k[b + 1..]))
}
#[cfg(feature = "gate")]
fn set_met_records(_r: &Rd, m: &mut HourlyMeteorology, recs: Vec<HourlyMetRecord>) {
    m.gate_set_records(recs);
}
#[cfg(not(feature = "gate"))]
fn set_met_records(r: &Rd, _m: &mut HourlyMeteorology, _recs: Vec<HourlyMetRecord>) {
    r.fail("hourly meteorology records need the `gate` feature");
}
'''

def main():
    D.collect_enums()
    # 1. what ctxdump.py writes for the Rust dumper now (into a scratch root)
    with tempfile.TemporaryDirectory() as tmp:
        os.makedirs(f'{tmp}/gate/diff')
        saved = D.ROOT
        D.ROOT = tmp
        with contextlib.redirect_stdout(io.StringIO()):
            D.main()
        D.ROOT = saved
        ref = open(f'{tmp}/crates/engine/examples/common/ctxdump.rs').read().split('\n')
    committed = f'{ROOT}/crates/engine/examples/common/ctxdump.rs'
    if os.path.exists(committed) and open(committed).read().split('\n') != ref:
        print(f'ctxload: warning: {committed} is not what tools/ctxdump.py generates now; regenerate it')
    # 2. the walk, with ctxdump's emitters re-deriving the dumper
    D.CNT[0] = 0
    del D.SKIPS[:]
    R, L, RE, LE = walk()
    want = fn_body(ref, 'pub fn dump(c: &SimulationContext) -> Vec<u8> {')
    want_e = fn_body(ref, 'pub fn engdump(e: &SwmmEngine) -> Vec<u8> {')
    for name, got, exp in (('dump', R, want), ('engdump', RE, want_e)):
        if got != exp:
            i = next(i for i in range(max(len(got), len(exp))) if i >= len(got) or i >= len(exp) or got[i] != exp[i])
            raise SystemExit(f'ctxload: ctxdump.py changed ({name}() line {i}):\n  ctxdump: {exp[i] if i < len(exp) else "<end>"}\n  ctxload: {got[i] if i < len(got) else "<end>"}\nupdate tools/ctxload.py')
    # 3. enums: from upstream's integer value back to the translation's variant
    E = []
    for path, names in rust_enums():
        E.append(f'impl FromI64 for {path} {{')
        E.append('    fn from_i64(v: i64) -> Self {')
        E.append(f'        use {path} as E;')
        for n in names:
            E.append(f'        if v == E::{n} as i64 {{ return E::{n}; }}')
        E.append(f'        bad_enum("{path.rsplit("::", 1)[1]}", v)')
        E.append('    }')
        E.append('}')
    out = [re.sub(r'^fn ', 'pub(crate) fn ', PRELUDE, flags=re.M)]
    out += E
    out += ['/// Fill `c` (a fresh context) from the dump at `r`, as `ctxdump::dump` printed it.',
            'pub fn load(r: &mut Rd, c: &mut SimulationContext) {']
    out += L
    out += ['}', '',
            '/// The members `ctxdump::engdump` prints, into a fresh engine (needs the `gate` feature).',
            '#[cfg(feature = "gate")]',
            'pub fn load_engine(r: &mut Rd, e: &mut SwmmEngine) {',
            '    let mut m = e.gate_members();']
    out += LE
    out += ['}', '#[cfg(not(feature = "gate"))]',
            'pub fn load_engine(r: &mut Rd, _e: &mut SwmmEngine) {',
            '    r.fail("the engine members need the `gate` feature");', '}']
    LR, LQ = walk_rt()
    out += ['', '/// The members `rtdump::rtdump` prints, over solvers the init legs built (needs the `gate`',
            '/// feature). See tools/ctxload.py walk_rt.', '#[cfg(feature = "gate")]',
            'pub fn load_rt(r: &mut Rd, e: &mut SwmmEngine) {', '    let g = e.gate_solvers();']
    out += LR
    out += ['}', '', '/// The members `ctxdump::qualdump` prints (needs the `gate` feature).', '#[cfg(feature = "gate")]',
            'pub fn load_qual(r: &mut Rd, e: &mut SwmmEngine) {', '    let q = e.gate_quality();']
    out += LQ
    out += ['}']
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    open(OUT, 'w').write('\n'.join(out) + '\n')
    # 4. the hydrology loader, checked against hydrodump.rs the same way
    with tempfile.TemporaryDirectory() as tmp:
        os.makedirs(f'{tmp}/gate/diff'); os.makedirs(f'{tmp}/crates/engine/examples/common')
        saved = D.ROOT
        D.ROOT = tmp
        try:
            D.CNT[0] = 0
            with contextlib.redirect_stdout(io.StringIO()):
                D.hydro_dumper()
        finally:
            D.ROOT = saved
        href = open(f'{tmp}/crates/engine/examples/common/hydrodump.rs').read().split('\n')
    D.CNT[0] = 0
    HR, HG, HL = hydro_walk()
    i0 = href.index('    let o = &mut o;') + 1
    want_h = href[i0:i0 + len(HR)]
    if HR != want_h:
        i = next(i for i in range(len(HR)) if HR[i] != want_h[i])
        raise SystemExit(f'ctxload: hydrodump changed (line {i}):\n  hydrodump: {want_h[i]}\n  ctxload: {HR[i]}\nupdate tools/ctxload.py')
    # the solvers, then the four hand-written lines (uh_params, uh_name_to_idx_, the
    # rpt_files_ size and keys), then gage_states_
    ig = i0 + len(HR) + 4
    if href[ig:ig + len(HG)] != HG:
        raise SystemExit('ctxload: hydrodump changed (gage_states_ or the hand-written lines before it); update tools/ctxload.py')
    hout = [HPRELUDE, '/// Fill the hydrology solvers of a started engine from the dump at `r`, as',
            '/// `hydrodump::hydrodump` printed it (needs the `gate` feature).', '#[cfg(feature = "gate")]',
            'pub fn load_hydro(r: &mut Rd, e: &mut SwmmEngine) {', '    let mut h = e.gate_hydrology();']
    hout += HL
    hout += ['}', '#[cfg(not(feature = "gate"))]', 'pub fn load_hydro(r: &mut Rd, _e: &mut SwmmEngine) {',
             '    r.fail("the hydrology solvers need the `gate` feature");', '}']
    open(HOUT, 'w').write('\n'.join(hout) + '\n')
    print(f'ctxload: hydroload {len(HL)} loader lines')
    print(f'ctxload: {len(L) + len(LE)} loader lines, {len(E) // 1} enum lines; {len(SKIPPED)} members skipped (not loadable):', '; '.join(SKIPPED))

if __name__ == '__main__':
    main()
