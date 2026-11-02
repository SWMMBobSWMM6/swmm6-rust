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
    m = re.search(r'(?:struct|class)\s+' + struct + r'\s*(?:final\s*)?(?::[^{;]*)?\{', src)
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
        m = re.fullmatch(r'((?:const\s+)?(?:unsigned\s+)?[\w:]+(?:<.+?>+)?(?:\s+long)?)\s+([A-Za-z_].*)', s)
        if not m:
            continue
        ty, names = m.group(1), m.group(2)
        for nm in split_top(names):          # each declarator, initializer dropped
            nm = re.sub(r'\s*=\s*.*$', '', nm.strip())
            nm = re.sub(r'\s*\{\}\s*$', '', nm).strip()
            am = re.fullmatch(r'(\w+)\s*\[([^\]]+)\]', nm)
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
    # LID storage nodes ([LID_NODES], [LID_NODE_OUTLETS]) and their run-time state
    'LidNodeConfig': ('data/LidNodeData.hpp', 'LidNodeConfig'),
    'LidNodeOutlet': ('data/LidNodeData.hpp', 'LidNodeOutlet'),
    'LidNodeState': ('data/LidNodeData.hpp', 'LidNodeState'),
    'LidNodeCell': ('data/LidNodeData.hpp', 'LidNodeCell'),
    'LidWaterTransfer': ('data/LidNodeData.hpp', 'LidWaterTransfer'),
    'LidPortTransfer': ('data/LidNodeData.hpp', 'LidPortTransfer'),
    'GreenAmptState': ('hydrology/Infiltration.hpp', 'GreenAmptState'),
    'richards::Report': ('hydrology/RichardsColumn.hpp', 'Report'),
    # post-parse resolver outputs: gage rain series tables, transect tables
    'TableCursor': ('data/TableData.hpp', 'TableCursor'),
    'TableBlock': ('data/TableData.hpp', 'TableBlock'),
    'xsect::LocateLut': ('hydraulics/XSectLookup.hpp', 'LocateLut'),
    'Transfer': ('hydrology/RichardsColumn.hpp', 'Transfer'),
    # SWMMEngine::open: the reactions component (ReactionData), heat/ARD/water-age config
    'ParameterOverride': ('data/ReactionData.hpp', 'ParameterOverride'),
    'WallArea': ('data/ReactionData.hpp', 'WallArea'),
    'LidLaw': ('data/ReactionData.hpp', 'LidLaw'),
    'LidState': ('data/ReactionData.hpp', 'LidState'),
    'SourcePattern': ('data/ReactionData.hpp', 'SourcePattern'),
    'Source': ('data/ReactionData.hpp', 'Source'),
    'SubsurfaceLaw': ('data/ReactionData.hpp', 'SubsurfaceLaw'),
    'SurfaceLaw': ('data/ReactionData.hpp', 'SurfaceLaw'),
    'MsxRoutingLedger': ('data/ReactionData.hpp', 'MsxRoutingLedger'),
    'RxToken': ('data/ReactionTokens.hpp', 'RxToken'),
    'RxExprSpan': ('data/ReactionTokens.hpp', 'RxExprSpan'),
    'MsxSurfaceData': ('data/MsxSurfaceData.hpp', 'MsxSurfaceData'),
    'MsxSurfaceRows': ('data/MsxSurfaceData.hpp', 'MsxSurfaceRows'),
    'ConductionConfig': ('data/HeatData.hpp', 'ConductionConfig'),
    'SolarConfig': ('data/HeatData.hpp', 'SolarConfig'),
    'CloudConfig': ('data/HeatData.hpp', 'CloudConfig'),
    'RadiativeConfig': ('data/HeatData.hpp', 'RadiativeConfig'),
    'SedimentConfig': ('data/BedZoneData.hpp', 'SedimentConfig'),
    'HeatOverrideData': ('data/HeatData.hpp', 'HeatOverrideData'),
    'HeatOverrideRow': ('data/HeatOverrideData.hpp', 'HeatOverrideRow'),
    'ArdTransportRow': ('data/ArdConfigData.hpp', 'ArdTransportRow'),
}

# `using` aliases of member types, spelled out
ALIASES = {
    'NodeCoupling2DParams': 'std::array<double, kNodeCoupling2DKeys>',
}

def kind(ty):
    ty = ty.replace('const ', '').strip()
    ty = ALIASES.get(ty, ty)
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
        if k[1] == 'FilePathPair':
            C.append(f'  ps("{label}.abs", ({path_c}).absolute);'); R.append(f'    ps(o, "{label}.abs", &({path_r}).absolute);')
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
                C.append(f'  {{ std::vector<std::string> t; for (auto& x : {path_c}) t.push_back(x.absolute); vs("{label}.abs", t); }}')
                R.append(f'    vs(o, "{label}.abs", &({path_r}).iter().map(|x| x.absolute.clone()).collect::<Vec<_>>());')
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
ELEM_STRUCTS = {'UnitHydEntry': 'data/InflowData.hpp', 'RDIIDecayEntry': 'data/InflowData.hpp',
                'Table': 'data/TableData.hpp'}

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
    ('data/ReactionData.hpp', 'ReactionData', 'c.reactions', 'c.reactions', 'reactions'),
    ('data/HeatData.hpp', 'HeatConfigData', 'c.heat_config', 'c.heat_config', 'heat_config'),
    ('data/ArdConfigData.hpp', 'ArdConfigData', 'c.ard_config', 'c.ard_config', 'ard_config'),
    ('data/WaterAgeData.hpp', 'WaterAgeConfigData', 'c.water_age_config', 'c.water_age_config', 'water_age_config'),
    ('data/InitialQualityData.hpp', 'InitialQualityData', 'c.initial_quality', 'c.initial_quality', 'initial_quality'),
    ('data/NodeCoupling2DData.hpp', 'NodeCoupling2DRows', 'c.node_coupling_2d', 'c.node_coupling_2d', 'node_coupling_2d'),
]
NAMES = ['node_names', 'link_names', 'subcatch_names', 'gage_names', 'pollutant_names', 'landuse_names',
         'snowpack_names', 'aquifer_names', 'lid_names']
CTX_SCALARS = ['error_code', 'warning_code', 'error_message', 'warnings', 'errors', 'title_notes',
               'adjust_temp', 'adjust_evap', 'adjust_rain', 'adjust_hydcon', 'subcatch_n_perv_pattern',
               'subcatch_d_store_pattern', 'subcatch_infil_pattern', 'base_n_perv', 'base_ds_perv',
               'has_subcatch_adj_patterns', 'xsect_generation', 'gpkg_units_internal',
               'pending_gw_nodes', 'pending_gw_aquifers', 'deferred_section_rows', 'embedded_component_sections',
               'lid_node_outlets', 'state', 'reported_species_names']

# SWMMEngine members that open() sets, dumped after swmm_engine_open (gate/diff/ctxopen.cpp).
# They are private: the C++ side reads them through the explicit-instantiation access rule
# (a pointer to a private member may be named in an explicit instantiation), which leaves the
# class layout untouched. Types are read from the header.
ENGINE_HEADER = 'core/SWMMEngine.hpp'
ENGINE_MEMBERS = ['lenient_open_', 'metadata_only_', 'rpt_path_', 'out_path_']
# private members of the built-in plugins, the same way: (header, class, member)
PLUGIN_MEMBERS = [('plugins/DefaultOutputPlugin.hpp', 'DefaultOutputPlugin', 'out_path_'),
                  ('plugins/DefaultReportPlugin.hpp', 'DefaultReportPlugin', 'rpt_path_')]

def member_type(header, name):
    src = strip(open(os.path.join(UP, header)).read())
    m = re.search(r'\n\s*((?:const\s+)?[\w:]+(?:<[^;{}]*?>)?)\s+' + name + r'\s*(?:=[^;]*)?;', src)
    if not m:
        raise SystemExit(f'no member {name} in {header}')
    return m.group(1)

def engine_dumper():
    C = ['#ifdef CTXDUMP_ENGINE', '#include "core/SWMMEngine.hpp"', '#include "plugins/DefaultInputPlugin.hpp"',
         '#include "plugins/DefaultOutputPlugin.hpp"', '#include "plugins/DefaultReportPlugin.hpp"',
         '#include "plugins/DefaultStateIOPlugin.hpp"', 'namespace ctxdump_access {',
         'template <typename Tag, typename Tag::type M> struct Rob { friend typename Tag::type get(Tag) { return M; } };']
    R = ['use swmm6_rust::core::engine::SwmmEngine;', 'pub fn engdump(e: &SwmmEngine) -> Vec<u8> {', '    let mut o: Vec<u8> = Vec::new();', '    let o = &mut o;']
    body_c = []
    for n in ENGINE_MEMBERS:
        t = member_type(ENGINE_HEADER, n)
        C.append(f'struct T_{n} {{ typedef {t} openswmm::SWMMEngine::*type; friend type get(T_{n}); }};')
        C.append(f'template struct Rob<T_{n}, &openswmm::SWMMEngine::{n}>;')
        k = kind(t)
        emit(f'(e.*get(T_{n}()))', f'e.{n.rstrip("_")}()', f'engine.{n}', k, body_c, R)
    for h, cls, n in PLUGIN_MEMBERS:
        t = member_type(h, n)
        C.append(f'struct T_{cls}_{n} {{ typedef {t} openswmm::{cls}::*type; friend type get(T_{cls}_{n}); }};')
        C.append(f'template struct Rob<T_{cls}_{n}, &openswmm::{cls}::{n}>;')
    t = member_type(ENGINE_HEADER, 'binary_output_plugin_') if False else None
    C.append('struct T_binary_output_plugin_ { typedef openswmm::DefaultOutputPlugin* openswmm::SWMMEngine::*type; friend type get(T_binary_output_plugin_); };')
    C.append('template struct Rob<T_binary_output_plugin_, &openswmm::SWMMEngine::binary_output_plugin_>;')
    C.append('} // namespace ctxdump_access')
    C.append('namespace ctxdump {')
    C.append('inline void engdump(openswmm::SWMMEngine& e) {')
    C.append('  using namespace ctxdump_access;')
    C += body_c
    # the plugin lists (polymorphic upstream; built-in types only in the translation)
    C.append('  const auto& pf = e.plugin_factory();')
    C.append('  pi("engine.plugins.input.n", (long long)pf.input_plugins().size());')
    C.append('  for (auto* p : pf.input_plugins()) { pi("engine.plugins.input.state", (long long)p->state()); auto* d = dynamic_cast<openswmm::DefaultInputPlugin*>(p); pi("engine.plugins.input.default", d != nullptr); if (d) { vs("engine.plugins.input.skipped", d->skipped_sections()); ps("engine.plugins.input.last_error", d->last_error_message()); } }')
    C.append('  pi("engine.plugins.output.n", (long long)pf.output_plugins().size());')
    C.append('  for (auto* p : pf.output_plugins()) { pi("engine.plugins.output.state", (long long)p->state()); auto* d = dynamic_cast<openswmm::DefaultOutputPlugin*>(p); pi("engine.plugins.output.default", d != nullptr); if (d) ps("engine.plugins.output.path", d->*get(T_DefaultOutputPlugin_out_path_())); }')
    C.append('  pi("engine.plugins.report.n", (long long)pf.report_plugins().size());')
    C.append('  for (auto* p : pf.report_plugins()) { pi("engine.plugins.report.state", (long long)p->state()); auto* d = dynamic_cast<openswmm::DefaultReportPlugin*>(p); pi("engine.plugins.report.default", d != nullptr); if (d) ps("engine.plugins.report.path", d->*get(T_DefaultReportPlugin_rpt_path_())); }')
    C.append('  pi("engine.plugins.state_io.n", (long long)pf.state_io_plugins().size());')
    C.append('  for (auto* p : pf.state_io_plugins()) { pi("engine.plugins.state_io.state", (long long)p->state()); pi("engine.plugins.state_io.default", dynamic_cast<openswmm::DefaultStateIOPlugin*>(p) != nullptr); }')
    C.append('  { long long k = -1; for (std::size_t i = 0; i < pf.output_plugins().size(); ++i) if (pf.output_plugins()[i] == e.*get(T_binary_output_plugin_())) k = (long long)i; pi("engine.binary_output_plugin_", k); }')
    C.append('}')
    C.append('} // namespace ctxdump')
    C.append('#endif // CTXDUMP_ENGINE')
    R.append('    let pf = e.plugin_factory();')
    R.append('    pi(o, "engine.plugins.input.n", pf.input_plugins().len() as i64);')
    R.append('    for p in pf.input_plugins() { pi(o, "engine.plugins.input.state", p.state() as i64); pi(o, "engine.plugins.input.default", 1); vs(o, "engine.plugins.input.skipped", &p.skipped_sections()); ps(o, "engine.plugins.input.last_error", p.last_error_message()); }')
    R.append('    pi(o, "engine.plugins.output.n", pf.output_plugins().len() as i64);')
    R.append('    for p in pf.output_plugins() { pi(o, "engine.plugins.output.state", p.state() as i64); pi(o, "engine.plugins.output.default", 1); ps(o, "engine.plugins.output.path", p.out_path()); }')
    R.append('    pi(o, "engine.plugins.report.n", pf.report_plugins().len() as i64);')
    R.append('    for p in pf.report_plugins() { pi(o, "engine.plugins.report.state", p.state() as i64); pi(o, "engine.plugins.report.default", 1); ps(o, "engine.plugins.report.path", p.rpt_path()); }')
    R.append('    pi(o, "engine.plugins.state_io.n", pf.state_io_plugins().len() as i64);')
    R.append('    for p in pf.state_io_plugins() { pi(o, "engine.plugins.state_io.state", p.state() as i64); pi(o, "engine.plugins.state_io.default", 1); }')
    R.append('    pi(o, "engine.binary_output_plugin_", e.binary_output_plugin().map(|k| k as i64).unwrap_or(-1));')
    R.append('    std::mem::take(o)')
    R.append('}')
    return C, R

def main():
    collect_enums()
    C = ['// generated by tools/ctxdump.py: prints a SimulationContext; the Rust twin prints the same',
         '#pragma once', '#include <cstdio>', '#include <cstring>', '#include <cstdint>', '#include <string>', '#include <vector>',
         '#include <algorithm>', '#include <variant>', '#include "core/SimulationContext.hpp"', 'namespace ctxdump {',
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
         'use std::io::Write;', 'use swmm6_rust::core::context::SimulationContext;', 'use swmm6_rust::core::user_flags::UserFlagValue;',
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
    C.append('  for (auto& t : c.tables.tables) { ps("table.id", t.id); pi("table.type", (long long)t.type); ps("table.comment", t.comment); vf("table.x", t.x); vf("table.y", t.y); pi("table.n_relative", t.n_relative); pf("table.rel_anchor", t.rel_anchor); pi("table.is_file", t.is_file_based); ps("table.file", t.file_path.original); ps("table.file.abs", t.file_path.absolute); }')
    R.append('    for t in &c.tables.tables { ps(o, "table.id", &t.id); pi(o, "table.type", t.type_ as i64); ps(o, "table.comment", &t.comment); vf(o, "table.x", &t.x); vf(o, "table.y", &t.y); pi(o, "table.n_relative", t.n_relative as i64); pf(o, "table.rel_anchor", t.rel_anchor); pi(o, "table.is_file", t.is_file_based as i64); ps(o, "table.file", &t.file_path.original); ps(o, "table.file.abs", &t.file_path.absolute); }')
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
    # specs, user flags and hotstart saves (members dump_struct cannot express)
    C.append('  for (auto& p : c.plugin_specs) { ps("plugin.path", p.path); vs("plugin.init_args", p.init_args); }')
    R.append('    for p in &c.plugin_specs { ps(o, "plugin.path", &p.path); vs(o, "plugin.init_args", &p.init_args); }')
    C.append('  for (auto& p : c.process_component_specs) { ps("pcomp.id", p.id); ps("pcomp.config_path", p.config_path); for (auto& kv : p.args) { ps("pcomp.arg", kv.first); ps("pcomp.arg.value", kv.second); } ps("pcomp.resolved_config_path", p.resolved_config_path); }')
    R.append('    for p in &c.process_component_specs { ps(o, "pcomp.id", &p.id); ps(o, "pcomp.config_path", &p.config_path); for kv in &p.args { ps(o, "pcomp.arg", &kv.0); ps(o, "pcomp.arg.value", &kv.1); } ps(o, "pcomp.resolved_config_path", &p.resolved_config_path); }')
    C.append('  for (auto& d : c.user_flags.all_defs()) { ps("uflag_def.name", d.name); pi("uflag_def.type", (long long)d.type); ps("uflag_def.description", d.description); }')
    R.append('    for d in c.user_flags.all_defs() { ps(o, "uflag_def.name", &d.name); pi(o, "uflag_def.type", d.type_ as i64); ps(o, "uflag_def.description", &d.description); }')
    # the value map is an unordered_map upstream: sorted by stored key on both sides
    C.append('  { std::vector<std::pair<std::string, const openswmm::UserFlagValue*>> e; for (auto& kv : c.user_flags.all_values()) e.emplace_back(kv.first, &kv.second); std::sort(e.begin(), e.end(), [](const auto& a, const auto& b) { return a.first < b.first; }); for (auto& kv : e) { ps("uflag_value", kv.first); pi("uflag_value.kind", (long long)kv.second->index()); if (auto b = std::get_if<bool>(kv.second)) pi("uflag_value.v", *b); else if (auto i = std::get_if<int>(kv.second)) pi("uflag_value.v", *i); else if (auto d = std::get_if<double>(kv.second)) pf("uflag_value.v", *d); else ps("uflag_value.v", std::get<std::string>(*kv.second)); } }')
    R.append('    { let mut e: Vec<_> = c.user_flags.all_values().iter().collect(); e.sort_by(|a, b| a.0.cmp(b.0)); for (k, v) in e { ps(o, "uflag_value", k); match v { UserFlagValue::Bool(b) => { pi(o, "uflag_value.kind", 0); pi(o, "uflag_value.v", *b as i64); } UserFlagValue::Int(i) => { pi(o, "uflag_value.kind", 1); pi(o, "uflag_value.v", *i as i64); } UserFlagValue::Double(d) => { pi(o, "uflag_value.kind", 2); pf(o, "uflag_value.v", *d); } UserFlagValue::String(s) => { pi(o, "uflag_value.kind", 3); ps(o, "uflag_value.v", s); } } } }')
    C.append('  for (auto& h : c.files.hotstart_saves) { ps("hotstart_save.path", h.path.original); ps("hotstart_save.path.abs", h.path.absolute); pf("hotstart_save.datetime", h.datetime); pi("hotstart_save.format", (long long)h.format); }')
    R.append('    for h in &c.files.hotstart_saves { ps(o, "hotstart_save.path", &h.path.original); ps(o, "hotstart_save.path.abs", &h.path.absolute); pf(o, "hotstart_save.datetime", h.datetime); pi(o, "hotstart_save.format", h.format as i64); }')
    for header, struct, pc, pr, label in STORES:
        dump_struct(header, struct, pc, pr, label, C, R)
    # written by input/PostParseResolver.cpp (resolve_cross_references)
    C.append('  pi("transect_tables.n", (long long)c.transect_tables.size());'); R.append('    pi(o, "transect_tables.n", c.transect_tables.len() as i64);')
    C.append('  for (auto& el : c.transect_tables) {'); R.append('    for el in &c.transect_tables {')
    dump_struct('hydraulics/Transect.hpp', 'TransectData', 'el', 'el', 'transect_tables[]', C, R)
    C.append('  }'); R.append('    }')
    C.append('  pi("met.n", (long long)c.snowpacks.meteorology.records().size()); for (auto& r : c.snowpacks.meteorology.records()) { pf("met.date", r.date); for (double v : r.values) pf("met.value", v); }')
    R.append('    pi(o, "met.n", c.snowpacks.meteorology.records().len() as i64); for r in c.snowpacks.meteorology.records().iter() { pf(o, "met.date", r.date); for v in r.values.iter() { pf(o, "met.value", *v); } }')
    # the species registry (private members, read through its accessors)
    C.append('  pi("species.n", c.species_registry.count()); pi("species.pollutants", c.species_registry.pollutant_count()); for (int i = 0; i < c.species_registry.count(); ++i) { ps("species.name", c.species_registry.name(i)); pi("species.kind", (long long)c.species_registry.kind(i)); ps("species.units", c.species_registry.units(i)); }')
    R.append('    pi(o, "species.n", c.species_registry.count() as i64); pi(o, "species.pollutants", c.species_registry.pollutant_count() as i64); for i in 0..c.species_registry.count() { ps(o, "species.name", c.species_registry.name(i)); pi(o, "species.kind", c.species_registry.kind(i) as i64); ps(o, "species.units", c.species_registry.units(i)); }')
    C.append('}'); C.append('} // namespace ctxdump')
    R.append('    std::mem::take(o)'); R.append('}')
    EC, ER = engine_dumper()
    C += EC
    R += ER
    open(f'{ROOT}/gate/diff/ctxdump.hpp', 'w').write('\n'.join(C) + '\n')
    os.makedirs(f'{ROOT}/crates/engine/examples/common', exist_ok=True)
    open(f'{ROOT}/crates/engine/examples/common/ctxdump.rs', 'w').write('\n'.join(R) + '\n')
    print(f'{len(SKIPS)} skipped:', '; '.join(SKIPS))

if __name__ == '__main__':
    main()
