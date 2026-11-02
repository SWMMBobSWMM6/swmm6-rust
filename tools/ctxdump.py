#!/usr/bin/env python3
"""Generate a SimulationContext dumper twice: in C++ for upstream (gate/diff/ctxdump.hpp) and
in Rust for the translation (crates/engine/examples/common/ctxdump.rs). Both print the same
text for the same state, so `cmp` on the two dumps compares every listed value.

Field lists and element types come from upstream's headers, parsed here. Values print as:
  double/float -> hex bits;  integers, bool, enums -> decimal;  std::string -> length:bytes;
  nested std::vector of those -> count then elements.  Anything else prints "skip <name>"
on both sides, and the skip list is reported, so nothing is silently unchecked.
"""
import re, os, sys

UP = '/root/work/swmm6src/upstream/src/engine'
ROOT = os.environ.get('CTXDUMP_ROOT', '/root/work/swmm6-rust')
RUST_KEYWORDS = {'type', 'match', 'loop', 'move', 'ref', 'self', 'fn', 'mod', 'use', 'impl', 'in', 'where', 'box', 'crate', 'super'}

INTS = {'int', 'long', 'bool', 'char', 'uint8_t', 'int8_t', 'int16_t', 'uint16_t', 'int32_t', 'uint32_t', 'int64_t',
        'uint64_t', 'std::size_t', 'size_t', 'long long', 'unsigned', 'unsigned char', 'std::uint8_t', 'std::int8_t', 'std::int32_t',
        'std::int64_t', 'std::uint32_t', 'std::uint64_t', 'std::int16_t', 'short'}
FLOATS = {'double', 'float'}
STRINGS = {'std::string', 'FilePathPair'}

def strip(src):
    src = re.sub(r'/\*.*?\*/', '', src, flags=re.S)
    src = re.sub(r'//[^\n]*', '', src)
    return src

ENUMS = set()
def collect_enums():
    for d, _, fs in os.walk(UP):
        if '/2d' in d:
            continue
        for f in fs:
            if f.endswith('.hpp'):
                for m in re.finditer(r'enum\s+(?:class\s+)?(\w+)\s*(?::\s*[\w:]+)?\s*\{', strip(open(os.path.join(d, f), errors='replace').read())):
                    ENUMS.add(m.group(1))

def struct_body(header, struct):
    src = strip(open(os.path.join(UP, header)).read())
    m = re.search(r'struct\s+' + struct + r'\s*\{', src)
    if not m:
        raise SystemExit(f'no struct {struct} in {header}')
    i = m.end(); d = 1
    while d:
        d += {'{': 1, '}': -1}.get(src[i], 0); i += 1
    return src[m.end():i - 1]

def split_top(s, sep=','):
    """Split on `sep` outside <>, () and []."""
    out = []; d = 0; cur = ''
    for ch in s:
        if ch in '<([': d += 1
        elif ch in '>)]': d -= 1
        if ch == sep and d == 0:
            out.append(cur); cur = ''
        else:
            cur += ch
    out.append(cur)
    return out

def fields(header, struct):
    """(name, type, arraylen) for every data member at depth 0, in order."""
    body = struct_body(header, struct)
    # remove nested braces (methods, nested structs) to depth 0 text, keeping ';' structure
    depth = 0; buf = ''
    for ch in body:
        if ch == '{':
            depth += 1
            if depth == 1:
                buf += ' {} '
            continue
        if ch == '}':
            depth -= 1
            if depth == 0:
                buf += ';'   # an inline method body ends a statement too
            continue
        if depth == 0:
            buf += ch
    res = []
    for stmt in buf.split(';'):
        s = ' '.join(stmt.split())
        if s.startswith(('public:', 'private:', 'protected:')):
            s = s.split(':', 1)[1].strip()
        if not s or s.startswith(('using ', 'static ', 'friend ', 'template', 'enum ', 'struct ', 'class ', 'typedef')):
            continue
        if '(' in s.split('=')[0]:
            continue
        m = re.fullmatch(r'((?:const\s+)?[\w:]+(?:<.+?>+)?(?:\s+long)?)\s+([A-Za-z_].*)', s)
        if not m:
            continue
        ty, names = m.group(1), m.group(2)
        for nm in split_top(names):          # each declarator, initializer dropped
            nm = re.sub(r'\s*=\s*.*$', '', nm.strip())
            nm = re.sub(r'\s*\{\}\s*$', '', nm).strip()
            am = re.fullmatch(r'(\w+)\s*\[\s*(\w+)\s*\]', nm)
            if am:
                res.append((am.group(1), ty, am.group(2)))
            elif re.fullmatch(r'\w+', nm):
                res.append((nm, ty, None))
    return res

SKIPS = []

def rname(n):
    return n + '_' if n in RUST_KEYWORDS else n

# Aggregate member types dumped member by member: type text as written -> (header, struct).
STRUCTS = {
    'ueb::Parameters': ('hydrology/surface/UebSnow.hpp', 'Parameters'),
    'ueb::Site': ('hydrology/surface/UebSnow.hpp', 'Site'),
    'SnowpackConfig': ('hydrology/surface/UebConfig.hpp', 'SnowpackConfig'),
    'RadiationSource': ('hydrology/surface/UebConfig.hpp', 'RadiationSource'),
    'ClimateSite': ('hydrology/surface/UebConfig.hpp', 'ClimateSite'),
    'surface::SnowConfig': ('hydrology/surface/UebConfig.hpp', 'SnowConfig'),
    'MeshSnowRow': ('hydrology/surface/MeshSnowConfig.hpp', 'MeshSnowRow'),
    'surface::MeshSnowConfig': ('hydrology/surface/MeshSnowConfig.hpp', 'MeshSnowConfig'),
    'LidNodeLayer': ('data/LidNodeData.hpp', 'LidNodeLayer'),
    'LidLayerTreatment': ('data/LidNodeData.hpp', 'LidLayerTreatment'),
    'richards::Material': ('hydrology/RichardsColumn.hpp', 'Material'),
    'richards::Options': ('hydrology/RichardsColumn.hpp', 'Options'),
    'Buildup': ('data/MsxSurfaceData.hpp', 'Buildup'),
    'Washoff': ('data/MsxSurfaceData.hpp', 'Washoff'),
    'Loading': ('data/MsxSurfaceData.hpp', 'Loading'),
    'openswmm::treatment::TreatExpr': ('quality/Treatment.hpp', 'TreatExpr'),
    'Token': ('quality/Treatment.hpp', 'Token'),
}

def kind(ty):
    ty = ty.replace('const ', '').strip()
    if ty in FLOATS: return ('f', ty)
    if ty in INTS or ty in ENUMS or ty.split('::')[-1] in ENUMS: return ('i', ty)
    if ty in STRINGS: return ('s', ty)
    if ty in STRUCTS: return ('struct',) + STRUCTS[ty]
    m = re.fullmatch(r'std::vector<std::array<double,\s*(\w+)>>', ty)
    if m:
        return ('va', m.group(1))
    m = re.fullmatch(r'std::(vector|optional|array|map|pair)<(.+)>', ty)
    if not m:
        return None
    args = [a.strip() for a in split_top(m.group(2))]
    sub = [kind(a) for a in args]
    if m.group(1) in ('vector', 'optional'):
        return (m.group(1)[0] if m.group(1) == 'vector' else 'opt', sub[0]) if sub[0] else None
    if m.group(1) == 'array':
        return ('arr', sub[0], args[1]) if sub[0] else None
    if None in sub:
        return None
    return (m.group(1), sub[0], sub[1])

CNT = [0]
def fresh(p):
    CNT[0] += 1
    return f'{p}{CNT[0]}'

def emit(path_c, path_r, label, k, C, R):
    """Append code printing one value of kind k."""
    if k[0] == 'f':
        C.append(f'  pf("{label}", (double)({path_c}));'); R.append(f'    pf(o, "{label}", ({path_r}) as f64);')
    elif k[0] == 'i':
        C.append(f'  pi("{label}", (long long)({path_c}));'); R.append(f'    pi(o, "{label}", ({path_r}) as i64);')
    elif k[0] == 's':
        sc = f'std::string({path_c})' if k[1] == 'FilePathPair' else path_c
        sr = f'&({path_r}).original' if k[1] == 'FilePathPair' else f'&{path_r}'
        C.append(f'  ps("{label}", {sc});'); R.append(f'    ps(o, "{label}", {sr});')
    elif k[0] == 'v' and k[1][0] in 'fis':
        ik = k[1]
        if ik[0] == 'i':
            C.append(f'  {{ std::vector<long long> t; for (auto x : {path_c}) t.push_back((long long)x); vi("{label}", t); }}')
            R.append(f'    vi(o, "{label}", &({path_r}).iter().map(|x| *x as i64).collect::<Vec<i64>>());')
        elif ik[0] == 'f':
            C.append(f'  {{ std::vector<double> t; for (auto x : {path_c}) t.push_back((double)x); vf("{label}", t); }}')
            R.append(f'    vf(o, "{label}", &({path_r}).iter().map(|x| *x as f64).collect::<Vec<f64>>());')
        else:
            if ik[1] == 'FilePathPair':
                C.append(f'  {{ std::vector<std::string> t; for (auto& x : {path_c}) t.push_back(x.original); vs("{label}", t); }}')
                R.append(f'    vs(o, "{label}", &({path_r}).iter().map(|x| x.original.clone()).collect::<Vec<_>>());')
            else:
                C.append(f'  vs("{label}", {path_c});'); R.append(f'    vs(o, "{label}", &{path_r});')
    elif k[0] == 'va':
        C.append(f'  {{ std::vector<double> t; for (auto& a : {path_c}) for (double x : a) t.push_back(x); vf("{label}", t); }}')
        R.append(f'    vf(o, "{label}", &({path_r}).iter().flat_map(|a| a.iter().copied()).collect::<Vec<f64>>());')
    elif k[0] in ('v', 'arr'):
        # count, then each element (nested vectors, arrays, structs, pairs)
        e = fresh('e')
        C.append(f'  pi("{label}.n", (long long)({path_c}).size());'); R.append(f'    pi(o, "{label}.n", ({path_r}).len() as i64);')
        C.append(f'  for (const auto& {e} : {path_c}) {{'); R.append(f'    for {e} in ({path_r}).iter() {{')
        emit(e, f'(*{e})', label + '[]', k[1], C, R)
        C.append('  }'); R.append('    }')
    elif k[0] == 'opt':
        e = fresh('e')
        C.append(f'  if ({path_c}) {{'); R.append(f'    match &({path_r}) {{ Some({e}) => {{')
        emit(f'(*({path_c}))', f'(*{e})', label, k[1], C, R)
        C.append(f'  }} else pn("{label}");'); R.append(f'    }} None => pn(o, "{label}") }}')
    elif k[0] == 'map':
        e = fresh('e')
        C.append(f'  pi("{label}.n", (long long)({path_c}).size());'); R.append(f'    pi(o, "{label}.n", ({path_r}).len() as i64);')
        C.append(f'  for (const auto& {e} : {path_c}) {{'); R.append(f'    for {e} in ({path_r}).iter() {{')
        emit(f'{e}.first', f'(*{e}.0)', label + '.key', k[1], C, R)
        emit(f'{e}.second', f'(*{e}.1)', label + '.val', k[2], C, R)
        C.append('  }'); R.append('    }')
    elif k[0] == 'pair':
        emit(f'{path_c}.first', f'{path_r}.0', label + '.first', k[1], C, R)
        emit(f'{path_c}.second', f'{path_r}.1', label + '.second', k[2], C, R)
    elif k[0] == 'struct':
        dump_struct(k[1], k[2], path_c, path_r, label, C, R)

def dump_struct(header, struct, path_c, path_r, label, C, R):
    for name, ty, arr in fields(header, struct):
        k = kind(ty)
        pc, pr, lb = f'{path_c}.{name}', f'{path_r}.{rname(name)}', f'{label}.{name}'
        if arr:
            if k and k[0] in 'fi':
                C.append(f'  for (int q = 0; q < (int)(sizeof({pc})/sizeof({pc}[0])); ++q) {{')
                R.append(f'    for q in 0..({pr}).len() {{')
                emit(f'{pc}[q]', f'{pr}[q]', lb + '[]', k, C, R)
                C.append('  }'); R.append('    }')
            else:
                SKIPS.append(lb); C.append(f'  sk("{lb}");'); R.append(f'    sk(o, "{lb}");')
            continue
        em = re.fullmatch(r'std::vector<(\w+)>', ty.replace('const ', '').strip())
        if k is None and em and em.group(1) in ELEM_STRUCTS:
            C.append(f'  pi("{lb}.n", (long long)({pc}).size());'); R.append(f'    pi(o, "{lb}.n", ({pr}).len() as i64);')
            C.append(f'  for (auto& el : {pc}) {{'); R.append(f'    for el in &{pr} {{')
            dump_struct(ELEM_STRUCTS[em.group(1)], em.group(1), 'el', 'el', lb + '[]', C, R)
            C.append('  }'); R.append('    }')
            continue
        if k is None:
            SKIPS.append(f'{lb} ({ty})'); C.append(f'  sk("{lb}");'); R.append(f'    sk(o, "{lb}");')
            continue
        emit(pc, pr, lb, k, C, R)

# element structs of std::vector members that are dumped field by field
ELEM_STRUCTS = {'UnitHydEntry': 'data/InflowData.hpp', 'RDIIDecayEntry': 'data/InflowData.hpp'}

# (header, struct, C++ path, Rust path, label)
STORES = [
    ('data/NodeData.hpp', 'NodeData', 'c.nodes', 'c.nodes', 'nodes'),
    ('data/LinkData.hpp', 'LinkData', 'c.links', 'c.links', 'links'),
    ('data/SubcatchData.hpp', 'SubcatchData', 'c.subcatches', 'c.subcatches', 'subcatches'),
    ('data/GageData.hpp', 'GageData', 'c.gages', 'c.gages', 'gages'),
    ('data/PollutantData.hpp', 'PollutantData', 'c.pollutants', 'c.pollutants', 'pollutants'),
    ('data/NodeSubtypes.hpp', 'StorageData', 'c.node_subtypes.storages', 'c.node_subtypes.storages', 'storages'),
    ('data/NodeSubtypes.hpp', 'OutfallData', 'c.node_subtypes.outfalls', 'c.node_subtypes.outfalls', 'outfalls'),
    ('data/NodeSubtypes.hpp', 'DividerData', 'c.node_subtypes.dividers', 'c.node_subtypes.dividers', 'dividers'),
    ('data/LinkSubtypes.hpp', 'ConduitData', 'c.link_subtypes.conduits', 'c.link_subtypes.conduits', 'conduits'),
    ('data/LinkSubtypes.hpp', 'PumpData', 'c.link_subtypes.pumps', 'c.link_subtypes.pumps', 'pumps'),
    ('data/LinkSubtypes.hpp', 'OrificeData', 'c.link_subtypes.orifices', 'c.link_subtypes.orifices', 'orifices'),
    ('data/LinkSubtypes.hpp', 'WeirData', 'c.link_subtypes.weirs', 'c.link_subtypes.weirs', 'weirs'),
    ('data/LinkSubtypes.hpp', 'OutletData', 'c.link_subtypes.outlets', 'c.link_subtypes.outlets', 'outlets'),
    ('data/QualityData.hpp', 'LanduseData', 'c.landuses', 'c.landuses', 'landuses'),
    ('data/QualityData.hpp', 'BuildupData', 'c.buildup', 'c.buildup', 'buildup'),
    ('data/QualityData.hpp', 'WashoffData', 'c.washoff', 'c.washoff', 'washoff'),
    ('data/QualityData.hpp', 'TreatmentData', 'c.treatment', 'c.treatment', 'treatment'),
    ('data/InflowData.hpp', 'ExtInflowData', 'c.ext_inflows', 'c.ext_inflows', 'ext_inflows'),
    ('data/InflowData.hpp', 'DwfData', 'c.dwf_inflows', 'c.dwf_inflows', 'dwf'),
    ('data/InflowData.hpp', 'RDIIAssignData', 'c.rdii_assigns', 'c.rdii_assigns', 'rdii'),
    ('data/InflowData.hpp', 'PatternData', 'c.patterns', 'c.patterns', 'patterns'),
    ('data/InflowData.hpp', 'UnitHydData', 'c.unit_hyds', 'c.unit_hyds', 'unit_hyds'),
    ('data/InflowData.hpp', 'RDIIDecayData', 'c.rdii_decay', 'c.rdii_decay', 'rdii_decay'),
    ('data/InfraData.hpp', 'TransectStore', 'c.transects', 'c.transects', 'transects'),
    ('data/InfraData.hpp', 'StreetStore', 'c.streets', 'c.streets', 'streets'),
    ('data/InfraData.hpp', 'InletStore', 'c.inlets', 'c.inlets', 'inlets'),
    ('data/InfraData.hpp', 'InletUsageStore', 'c.inlet_usages', 'c.inlet_usages', 'inlet_usages'),
    ('data/InfraData.hpp', 'ControlRuleStore', 'c.control_rules', 'c.control_rules', 'control_rules'),
    ('data/HydrologyData.hpp', 'SnowpackStore', 'c.snowpacks', 'c.snowpacks', 'snowpacks'),
    ('data/HydrologyData.hpp', 'AquiferStore', 'c.aquifers', 'c.aquifers', 'aquifers'),
    ('data/HydrologyData.hpp', 'LidControlStore', 'c.lid_controls', 'c.lid_controls', 'lid_controls'),
    ('data/HydrologyData.hpp', 'LidUsageStore', 'c.lid_usage', 'c.lid_usage', 'lid_usage'),
    ('core/SpatialFrame.hpp', 'SpatialFrame', 'c.spatial', 'c.spatial', 'spatial'),
    ('core/SimulationOptions.hpp', 'SimulationOptions', 'c.options', 'c.options', 'options'),
    ('hydraulics/fv/FvOptions.hpp', 'FvOptions', 'c.options.fv', 'c.options.fv', 'options.fv'),
    ('core/SimulationContext.hpp', 'FilesSpec', 'c.files', 'c.files', 'files'),
    ('data/MsxSurfaceData.hpp', 'MsxSurfaceRows', 'c.reactions.surface.rows', 'c.reactions.surface.rows', 'msx_rows'),
    ('data/InitialQualityData.hpp', 'InitialQualityData', 'c.initial_quality', 'c.initial_quality', 'initial_quality'),
]
NAMES = ['node_names', 'link_names', 'subcatch_names', 'gage_names', 'pollutant_names', 'landuse_names',
         'snowpack_names', 'aquifer_names', 'lid_names']
CTX_SCALARS = ['error_code', 'warning_code', 'error_message', 'warnings', 'errors', 'title_notes',
               'adjust_temp', 'adjust_evap', 'adjust_rain', 'adjust_hydcon', 'subcatch_n_perv_pattern',
               'subcatch_d_store_pattern', 'subcatch_infil_pattern', 'base_n_perv', 'base_ds_perv',
               'has_subcatch_adj_patterns', 'xsect_generation', 'gpkg_units_internal',
               'pending_gw_nodes', 'pending_gw_aquifers', 'deferred_section_rows']

def main():
    collect_enums()
    C = ['// generated by tools/ctxdump.py: prints a SimulationContext; the Rust twin prints the same',
         '#pragma once', '#include <cstdio>', '#include <cstring>', '#include <cstdint>', '#include <string>', '#include <vector>',
         '#include <algorithm>', '#include "core/SimulationContext.hpp"', 'namespace ctxdump {',
         'static void pf(const char* l, double v){ uint64_t u; std::memcpy(&u,&v,8); std::printf("%s %016llx\\n", l, (unsigned long long)u); }',
         'static void pi(const char* l, long long v){ std::printf("%s %lld\\n", l, v); }',
         'static void ps(const char* l, const std::string& s){ std::printf("%s %zu:", l, s.size()); std::fwrite(s.data(),1,s.size(),stdout); std::printf("\\n"); }',
         'static void vf(const char* l, const std::vector<double>& v){ std::printf("%s n=%zu\\n", l, v.size()); for (double x : v){ uint64_t u; std::memcpy(&u,&x,8); std::printf(" %016llx", (unsigned long long)u);} if(!v.empty()) std::printf("\\n"); }',
         'static void vi(const char* l, const std::vector<long long>& v){ std::printf("%s n=%zu\\n", l, v.size()); for (long long x : v) std::printf(" %lld", x); if(!v.empty()) std::printf("\\n"); }',
         'static void vs(const char* l, const std::vector<std::string>& v){ std::printf("%s n=%zu\\n", l, v.size()); for (auto& s : v){ std::printf(" %zu:", s.size()); std::fwrite(s.data(),1,s.size(),stdout);} if(!v.empty()) std::printf("\\n"); }',
         'static void sk(const char* l){ std::printf("skip %s\\n", l); }',
         'static void pn(const char* l){ std::printf("%s none\\n", l); }',
         'inline void dump(const openswmm::SimulationContext& c) {']
    R = ['// generated by tools/ctxdump.py: Rust twin of gate/diff/ctxdump.hpp', '#![allow(dead_code, unused_parens, clippy::all)]',
         'use std::io::Write;', 'use swmm6_rust::core::context::SimulationContext;',
         'fn pf(o: &mut Vec<u8>, l: &str, v: f64) { writeln!(o, "{} {:016x}", l, v.to_bits()).unwrap(); }',
         'fn pi(o: &mut Vec<u8>, l: &str, v: i64) { writeln!(o, "{} {}", l, v).unwrap(); }',
         'fn ps(o: &mut Vec<u8>, l: &str, s: &[u8]) { write!(o, "{} {}:", l, s.len()).unwrap(); o.extend_from_slice(s); o.push(b\'\\n\'); }',
         'fn vf(o: &mut Vec<u8>, l: &str, v: &[f64]) { writeln!(o, "{} n={}", l, v.len()).unwrap(); for x in v { write!(o, " {:016x}", x.to_bits()).unwrap(); } if !v.is_empty() { o.push(b\'\\n\'); } }',
         'fn vi(o: &mut Vec<u8>, l: &str, v: &[i64]) { writeln!(o, "{} n={}", l, v.len()).unwrap(); for x in v { write!(o, " {}", x).unwrap(); } if !v.is_empty() { o.push(b\'\\n\'); } }',
         'fn vs(o: &mut Vec<u8>, l: &str, v: &[Vec<u8>]) { writeln!(o, "{} n={}", l, v.len()).unwrap(); for s in v { write!(o, " {}:", s.len()).unwrap(); o.extend_from_slice(s); } if !v.is_empty() { o.push(b\'\\n\'); } }',
         'fn sk(o: &mut Vec<u8>, l: &str) { writeln!(o, "skip {}", l).unwrap(); }',
         'fn pn(o: &mut Vec<u8>, l: &str) { writeln!(o, "{} none", l).unwrap(); }',
         'pub fn dump(c: &SimulationContext) -> Vec<u8> {', '    let mut o: Vec<u8> = Vec::new();', '    let o = &mut o;']
    for n in NAMES:
        C.append(f'  vs("{n}", c.{n}.names());'); R.append(f'    vs(o, "{n}", c.{n}.names());')
    ctx_fields = {n: (t, a) for n, t, a in fields('core/SimulationContext.hpp', 'SimulationContext')}
    for n in CTX_SCALARS:
        t, a = ctx_fields[n]; k = kind(t)
        if a:
            C.append(f'  for (int q = 0; q < 12; ++q) {{'); R.append(f'    for q in 0..12 {{')
            emit(f'c.{n}[q]', f'c.{n}[q]', n + '[]', k, C, R); C.append('  }'); R.append('    }')
        else:
            emit(f'c.{n}', f'c.{n}', n, k, C, R)
    # tables
    C.append('  pi("tables.n", (long long)c.tables.tables.size());'); R.append('    pi(o, "tables.n", c.tables.tables.len() as i64);')
    C.append('  for (auto& t : c.tables.tables) { ps("table.id", t.id); pi("table.type", (long long)t.type); ps("table.comment", t.comment); vf("table.x", t.x); vf("table.y", t.y); pi("table.n_relative", t.n_relative); pf("table.rel_anchor", t.rel_anchor); pi("table.is_file", t.is_file_based); ps("table.file", t.file_path.original); }')
    R.append('    for t in &c.tables.tables { ps(o, "table.id", &t.id); pi(o, "table.type", t.type_ as i64); ps(o, "table.comment", &t.comment); vf(o, "table.x", &t.x); vf(o, "table.y", &t.y); pi(o, "table.n_relative", t.n_relative as i64); pf(o, "table.rel_anchor", t.rel_anchor); pi(o, "table.is_file", t.is_file_based as i64); ps(o, "table.file", &t.file_path.original); }')
    # unknown headers, deferred rows, passthrough
    C.append('  for (auto& h : c.unknown_section_headers) { ps("unknown_header", h.first); pi("unknown_header.line", h.second); }')
    R.append('    for h in &c.unknown_section_headers { ps(o, "unknown_header", &h.0); pi(o, "unknown_header.line", h.1); }')
    C.append('  for (auto& kv : c.passthrough_sections) { ps("passthrough", kv.first); vs("passthrough.lines", kv.second); }')
    R.append('    for (k, v) in &c.passthrough_sections { ps(o, "passthrough", k); vs(o, "passthrough.lines", v); }')
    C.append('  { std::vector<std::pair<std::string,std::string>> e(c.options.ext_options.begin(), c.options.ext_options.end()); std::sort(e.begin(), e.end()); for (auto& kv : e) { ps("ext_option", kv.first); ps("ext_option.value", kv.second); } }')
    R.append('    { let mut e: Vec<_> = c.options.ext_options.iter().collect(); e.sort(); for (k, v) in e { ps(o, "ext_option", k); ps(o, "ext_option.value", v); } }')
    C.append('  { std::vector<long long> t; for (auto x : c.node_subtypes.subtype_row) t.push_back(x); vi("node_subtypes.subtype_row", t); }')
    R.append('    vi(o, "node_subtypes.subtype_row", &c.node_subtypes.subtype_row.iter().map(|x| *x as i64).collect::<Vec<i64>>());')
    C.append('  { std::vector<long long> t; for (auto x : c.link_subtypes.subtype_row) t.push_back(x); vi("link_subtypes.subtype_row", t); }')
    R.append('    vi(o, "link_subtypes.subtype_row", &c.link_subtypes.subtype_row.iter().map(|x| *x as i64).collect::<Vec<i64>>());')
    C.append('  pi("events.n", (long long)c.events.size()); for (auto& e : c.events) { pf("event.start", e.start); pf("event.end", e.end); }')
    R.append('    pi(o, "events.n", c.events.len() as i64); for e in &c.events { pf(o, "event.start", e.start); pf(o, "event.end", e.end); }')
    C.append('  pi("pending_link_nodes.n", (long long)c.pending_link_nodes.size()); for (auto& p : c.pending_link_nodes) { pi("pending_link_node", p.first); ps("pending_link_node.n1", p.second.first); ps("pending_link_node.n2", p.second.second); }')
    R.append('    pi(o, "pending_link_nodes.n", c.pending_link_nodes.len() as i64); for p in &c.pending_link_nodes { pi(o, "pending_link_node", p.0 as i64); ps(o, "pending_link_node.n1", &(p.1).0); ps(o, "pending_link_node.n2", &(p.1).1); }')
    for header, struct, pc, pr, label in STORES:
        dump_struct(header, struct, pc, pr, label, C, R)
    C.append('}'); C.append('} // namespace ctxdump')
    R.append('    std::mem::take(o)'); R.append('}')
    open(f'{ROOT}/gate/diff/ctxdump.hpp', 'w').write('\n'.join(C) + '\n')
    os.makedirs(f'{ROOT}/crates/engine/examples/common', exist_ok=True)
    open(f'{ROOT}/crates/engine/examples/common/ctxdump.rs', 'w').write('\n'.join(R) + '\n')
    print(f'{len(SKIPS)} skipped:', '; '.join(SKIPS))

if __name__ == '__main__':
    main()
