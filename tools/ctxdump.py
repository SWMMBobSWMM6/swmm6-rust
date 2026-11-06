#!/usr/bin/env python3
"""Generate a SimulationContext dumper twice: in C++ for upstream (gate/diff/ctxdump.hpp) and
in Rust for the translation (crates/engine/examples/common/ctxdump.rs). Both print the same
text for the same state, so `cmp` on the two dumps compares every listed value.

Field lists and element types come from upstream's headers, parsed here. Values print as:
  double/float -> hex bits;  integers, bool, enums -> decimal;  std::string -> length:bytes;
  nested std::vector of those -> count then elements; a std::unordered_map -> its entries
  sorted by key (a CiHash map's keys case-folded first, as the translation stores them);
  a C array of structs -> each element.  Anything else prints "skip <name>" on both sides,
  and the skip list is reported, so nothing is silently unchecked.
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
    # remove nested braces (methods, nested structs) to depth 0 text, keeping ';' structure.
    # A nested `struct X { ... } member;` keeps its declarator: it becomes `X member`, so the
    # member is dumped (X in STRUCTS) or reported as a skip rather than dropped unseen.
    depth = 0; buf = ''; nested = None
    for ch in body:
        if ch == '{':
            depth += 1
            if depth == 1:
                m = re.search(r'(?:struct|class)\s+(\w+)\s*$', buf.split(';')[-1])
                nested = m.group(1) if m else None
                buf += ' {} '
            continue
        if ch == '}':
            depth -= 1
            if depth == 0:
                buf += ';'   # an inline method body ends a statement too
                if nested:
                    buf += nested + ' '
                nested = None
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
        if s.startswith('mutable '):
            s = s[len('mutable '):]
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
    # SWMMEngine::initQuality: the LEGACY quality router's FV transport state
    'FvQualityTransport': ('quality/FvQualityTransport.hpp', 'FvQualityTransport'),
    'FvTransportWater': ('quality/FvTransportWater.hpp', 'FvTransportWater'),
    'ConservativeMixing': ('quality/ConservativeMixing.hpp', 'ConservativeMixing'),
    'MixingWater': ('quality/ConservativeMixing.hpp', 'MixingWater'),
    'MixedScalar': ('quality/ConservativeMixing.hpp', 'MixedScalar'),
    # SimulationContext's run-time members (dumped since the initialize wave; CTX_MEMBERS)
    'climate::ClimateState': ('hydrology/Climate.hpp', 'ClimateState'),
    'RadiationState': ('hydrology/Radiation.hpp', 'RadiationState'),
    'MovingAvg7': ('hydrology/Climate.hpp', 'MovingAvg7'),
    'WaterAgeState': ('data/WaterAgeData.hpp', 'WaterAgeState'),
    'LidLayerSpeciesState': ('data/LidLayerSpeciesData.hpp', 'LidLayerSpeciesState'),
    'HeatState': ('data/HeatData.hpp', 'HeatState'),
    'BedZoneState': ('data/BedZoneData.hpp', 'BedZoneState'),
    'ForcingData': ('data/ForcingData.hpp', 'ForcingData'),
    'ElemClimateChannel': ('data/ForcingData.hpp', 'ElemClimateChannel'),
    'QualityOnlyInflow': ('core/SimulationContext.hpp', 'QualityOnlyInflow'),
    'MassBalance2D': ('core/SimulationContext.hpp', 'MassBalance2D'),
    'MaxStats': ('core/SimulationContext.hpp', 'MaxStats'),
    'RoutingStepStats': ('core/SimulationContext.hpp', 'RoutingStepStats'),
    'VJDiag': ('core/SimulationContext.hpp', 'VJDiag'),
    'InletDiag': ('core/SimulationContext.hpp', 'InletDiag'),
    'ControlLogEntry': ('core/SimulationContext.hpp', 'ControlLogEntry'),
    # SnowpackStore's UEB run-time state and daily weather
    'surface::ueb::NativeState': ('hydrology/surface/UebSubcatchments.hpp', 'NativeState'),
    'HeldState': ('hydrology/surface/UebHeld.hpp', 'HeldState'),
    'LumpedState': ('hydrology/surface/UebLumped.hpp', 'LumpedState'),
    'LumpedFluxes': ('hydrology/surface/UebLumped.hpp', 'LumpedFluxes'),
    'State': ('hydrology/surface/UebSnow.hpp', 'State'),
    'Fluxes': ('hydrology/surface/UebSnow.hpp', 'Fluxes'),
    'ThermalSample': ('hydrology/surface/UebSnow.hpp', 'ThermalSample'),
    'climate::SnowDailyWeather': ('hydrology/HourlyMeteorology.hpp', 'SnowDailyWeather'),
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
    m = re.fullmatch(r'std::unordered_map<(.+)>', ty)
    if m:
        args = [a.strip() for a in split_top(m.group(1))]
        kk, kv = kind(args[0]), kind(args[1])
        if kk is None or kv is None or kk[0] not in 'is' or len(args) not in (2, 4):
            return None
        # ('umap', key kind, value kind, fold): a CiHash map's keys are case-folded first
        return ('umap', kk, kv, len(args) == 4 and args[2] == 'CiHash')
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
    elif k[0] == 'umap':
        # iteration order is libstdc++'s bucket order upstream: printed sorted by key instead
        e = fresh('e'); kc = 'std::string' if k[1][0] == 's' else 'long long'
        fc = ' for (auto& ch : k) ch = (ch >= 0x61 && ch <= 0x7a) ? (char)(ch & ~32) : ch;' if k[3] else ''
        fr = '.iter().map(|&b| swmm6_rust::data::name_index::ci_upper(b)).collect::<Vec<u8>>()' if k[3] else '.clone()'
        C.append(f'  {{ std::vector<std::pair<{kc}, const decltype(({path_c}).begin()->second)*>> t; for (const auto& kv : {path_c}) {{ {kc} k = kv.first;{fc} t.emplace_back(k, &kv.second); }}'
                 f' std::sort(t.begin(), t.end(), [](const auto& a, const auto& b) {{ return a.first < b.first; }}); pi("{label}.n", (long long)t.size()); for (const auto& {e} : t) {{')
        R.append(f'    {{ let mut t: Vec<_> = ({path_r}).iter().map(|(k, v)| (k{fr}, v)).collect(); t.sort_by(|a, b| a.0.cmp(&b.0)); pi(o, "{label}.n", t.len() as i64); for {e} in t.iter() {{')
        emit(f'{e}.first', f'{e}.0', label + '.key', k[1], C, R)
        emit(f'(*{e}.second)', f'(*{e}.1)', label + '.val', k[2], C, R)
        C.append('  } }'); R.append('    } }')
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
            elif k and k[0] == 'struct':
                q = fresh('q')
                C.append(f'  for (int {q} = 0; {q} < (int)(sizeof({pc})/sizeof({pc}[0])); ++{q}) {{')
                R.append(f'    for {q} in 0..({pr}).len() {{')
                dump_struct(k[1], k[2], f'{pc}[{q}]', f'{pr}[{q}]', lb + '[]', C, R)
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
                'Table': 'data/TableData.hpp',
                # SWMMEngine::initQuality: land-use parameters, FV transport workspaces
                'BuildupParams': 'quality/Landuse.hpp', 'WashoffParams': 'quality/Landuse.hpp',
                'MixedScalar': 'quality/ConservativeMixing.hpp', 'MixingScalar': 'quality/ConservativeMixing.hpp',
                'WaterTransfer': 'quality/ConservativeMixing.hpp',
                'BoundaryWaterImport': 'quality/FvTransportWater.hpp'}

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
    # set by SWMMEngine::initQuality (resize_quality, qual_init_buildup) and initMassBalance
    ('core/SimulationContext.hpp', 'MassBalance', 'c.mass_balance', 'c.mass_balance', 'mass_balance'),
    ('core/SimulationContext.hpp', 'OutfallExports', 'c.outfall_exports', 'c.outfall_exports', 'outfall_exports'),
    ('core/SimulationContext.hpp', 'NegativeSourceStats', 'c.negsrc', 'c.negsrc', 'negsrc'),
]
NAMES = ['node_names', 'link_names', 'subcatch_names', 'gage_names', 'pollutant_names', 'landuse_names',
         'snowpack_names', 'aquifer_names', 'lid_names']
CTX_SCALARS = ['error_code', 'warning_code', 'error_message', 'warnings', 'errors', 'title_notes',
               'adjust_temp', 'adjust_evap', 'adjust_rain', 'adjust_hydcon', 'subcatch_n_perv_pattern',
               'subcatch_d_store_pattern', 'subcatch_infil_pattern', 'base_n_perv', 'base_ds_perv',
               'has_subcatch_adj_patterns', 'xsect_generation', 'gpkg_units_internal',
               'pending_gw_nodes', 'pending_gw_aquifers', 'deferred_section_rows', 'embedded_component_sections',
               'lid_node_outlets', 'state', 'reported_species_names']
# The context's run-time members, dumped (by ctx_member) after everything else. Left out, as
# not plain data: wall_start (the wall clock at open), state_accessors (closures over the
# engine's solvers), twod_io (pointers into the 2D module, null in a 1D build).
CTX_MEMBERS = ['current_time', 'current_date', 'elapsed_ms', 'old_elapsed_ms', 'next_report_ms', 'runoff_interp_f',
               'dt_controls_remaining', 'climate_state', 'water_age_state', 'lid_layer_state', 'heat_state', 'bed_state',
               'forcing', 'coupled_node', 'coupled_pond_area', 'coupled_pond_band', 'coupling_delivery_remaining',
               'surface_outfall_link_limit', 'surface_outfall_quality', 'link_gain_quality', 'runtime_coupling_used',
               'mass_balance_2d', 'max_courant_crit', 'max_flow_turns', 'max_non_converged', 'max_mass_bal_errs',
               'routing_stats', 'vj_diag', 'inlet_diag', 'control_log', 'control_rule_names', 'control_log_dropped',
               'inp_file_path']

def ctx_member(n, t, a, C, R):
    """One SimulationContext member of CTX_MEMBERS: a value, a struct, or a C array of either."""
    k = kind(t)
    if k is None:
        raise SystemExit(f'ctxdump: no kind for ctx member {n} ({t})')
    if a:
        q = fresh('q')
        C.append(f'  for (int {q} = 0; {q} < (int)(sizeof(c.{n})/sizeof(c.{n}[0])); ++{q}) {{')
        R.append(f'    for {q} in 0..(c.{n}).len() {{')
        if k[0] == 'struct':
            dump_struct(k[1], k[2], f'c.{n}[{q}]', f'c.{n}[{q}]', n + '[]', C, R)
        else:
            emit(f'c.{n}[{q}]', f'c.{n}[{q}]', n + '[]', k, C, R)
        C.append('  }'); R.append('    }')
    else:
        emit(f'c.{n}', f'c.{n}', n, k, C, R)

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
         '#include "plugins/DefaultStateIOPlugin.hpp"', 'namespace ctxdump_access {']
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

# ---------------------------------------------------------------------------
# SWMMEngine::initHydraulics state (gate/diff/inithyd.cpp, examples/inithyd.rs): the
# Router and everything it owns, the engine solvers initHydraulics initializes, and the
# context fields it writes that the context dumper does not cover. Written to its own
# pair of files (gate/diff/rtdump.hpp, crates/engine/examples/common/rtdump.rs) so the
# context dumper's users do not depend on the engine's solver members. The C++ side
# needs `#define private public` before SWMMEngine.hpp (inithyd.cpp does that).
# ---------------------------------------------------------------------------

# member types dumped member by member: type text as written -> (header, struct)
RT_STRUCTS = {
    'kinwave::KWSolver': ('hydraulics/KinematicWave.hpp', 'KWSolver'),
    'dynwave::DWSolver': ('hydraulics/DynamicWave.hpp', 'DWSolver'),
    'DWNodeArrays': ('hydraulics/DynamicWave.hpp', 'DWNodeArrays'),
    'DPSConfig': ('hydraulics/DynamicWave.hpp', 'DPSConfig'),
    'DPSLinkArrays': ('hydraulics/DynamicWave.hpp', 'DPSLinkArrays'),
    'fv::NetworkMeshData': ('hydraulics/fv/NetworkMeshData.hpp', 'NetworkMeshData'),
    'fv::NetworkStateData': ('hydraulics/fv/NetworkMeshData.hpp', 'NetworkStateData'),
    'fv::FvOptions': ('hydraulics/fv/FvOptions.hpp', 'FvOptions'),
    'XSectParams': ('hydraulics/XSectBatch.hpp', 'XSectParams'),
    'FvClosure': ('hydraulics/fv/FvClosureKernels.hpp', 'FvClosure'),
    'ExtInflowSoA': ('hydrology/Inflow.hpp', 'ExtInflowSoA'),
    'DwfInflowSoA': ('hydrology/Inflow.hpp', 'DwfInflowSoA'),
    'RDIIGroupSoA': ('hydrology/RDII.hpp', 'RDIIGroupSoA'),
    'ExfilSoA': ('hydraulics/Exfiltration.hpp', 'ExfilSoA'),
    'InletSoA': ('hydraulics/Inlet.hpp', 'InletSoA'),
    'PumpGroup': ('hydraulics/HydStructures.hpp', 'PumpGroup'),
    'OrificeGroup': ('hydraulics/HydStructures.hpp', 'OrificeGroup'),
    'WeirGroup': ('hydraulics/HydStructures.hpp', 'WeirGroup'),
    'OutletGroup': ('hydraulics/HydStructures.hpp', 'OutletGroup'),
    'InterfaceFileSoA': ('core/InterfaceFile.hpp', 'InterfaceFileSoA'),
    'ParseError': ('controls/Controls.hpp', 'ParseError'),
}
# element structs of std::vector members, dumped field by field: struct -> header
RT_ELEM_STRUCTS = {
    'VJuncPair': 'hydraulics/DynamicWave.hpp', 'NodeTile': 'hydraulics/DynamicWave.hpp',
    'FvGeometry': 'hydraulics/fv/NetworkMeshData.hpp', 'CulvertCurve': 'hydraulics/HydClosureKernels.hpp',
    'TimePattern': 'hydrology/Inflow.hpp', 'UHResponseData': 'hydrology/RDII.hpp',
    'GreenAmptState': 'hydrology/Infiltration.hpp', 'Rule': 'controls/Controls.hpp',
    'Premise': 'controls/Controls.hpp', 'Action': 'controls/Controls.hpp', 'PIDState': 'controls/Controls.hpp',
    'NamedVariable': 'controls/Controls.hpp', 'PremiseSoA': 'controls/Controls.hpp',
    'PendingAction': 'controls/Controls.hpp',
}

def _rt_shape_group(pc, pr, lb, C, R):
    """One XSectGroups ShapeGroup. Its table pointers print the tables they point at."""
    C.append(f'  pi("{lb}.shape", (long long)({pc}).shape); pi("{lb}.count", ({pc}).count);')
    R.append(f'    pi(o, "{lb}.shape", ({pr}).shape as i64); pi(o, "{lb}.count", ({pr}).count as i64);')
    C.append(f'  {{ std::vector<long long> t; for (auto x : ({pc}).link_idx) t.push_back(x); vi("{lb}.link_idx", t); }}')
    R.append(f'    vi(o, "{lb}.link_idx", &({pr}).link_idx.iter().map(|x| *x as i64).collect::<Vec<i64>>());')
    for f in ['y_full', 'a_full', 'r_full', 's_full', 'w_max', 'inv_y_full', 'y_bot', 'a_bot', 's_bot', 'r_bot']:
        C.append(f'  vf("{lb}.{f}", ({pc}).{f});'); R.append(f'    vf(o, "{lb}.{f}", &({pr}).{f});')
    for f in ['area_tables', 'hrad_tables', 'width_tables']:
        C.append(f'  pi("{lb}.{f}.n", (long long)({pc}).{f}.size()); for (const double* t : ({pc}).{f}) {{ if (t) vf("{lb}.{f}[]", std::vector<double>(t, t + ({pc}).transect_tbl_size)); else pn("{lb}.{f}[]"); }}')
        R.append(f'    pi(o, "{lb}.{f}.n", ({pr}).{f}.len() as i64); for t in ({pr}).{f}.iter() {{ match t {{ Some(t) => vf(o, "{lb}.{f}[]", &t[..({pr}).transect_tbl_size as usize]), None => pn(o, "{lb}.{f}[]") }} }}')
    C.append(f'  pi("{lb}.transect_tbl_size", ({pc}).transect_tbl_size);')
    R.append(f'    pi(o, "{lb}.transect_tbl_size", ({pr}).transect_tbl_size as i64);')
    for f in ['buf_d', 'buf_r', 'buf_r2']:
        C.append(f'  vf("{lb}.{f}", ({pc}).{f});'); R.append(f'    vf(o, "{lb}.{f}", &({pr}).{f}.borrow());')

def _rt_xsect_groups(pc, pr, lb, C, R):
    C.append(f'  pi("{lb}.n", (long long)({pc}).groups_.size());')
    R.append(f'    pi(o, "{lb}.n", ({pr}).num_groups() as i64);')
    C.append(f'  for (const auto& g : ({pc}).groups_) {{'); R.append(f'    for gi in 0..({pr}).num_groups() {{ let g = ({pr}).group(gi);')
    _rt_shape_group('g', '(*g)', lb + '[]', C, R)
    C.append('  }'); R.append('    }')
    # packed_groups_ / packed_count_ / mask_active_ are private in the translation (build
    # and attachTransectTables clear them); upstream's are printed as counts only
    C.append(f'  sk("{lb}.packed_groups_"); sk("{lb}.packed_count_"); sk("{lb}.mask_active_");')
    R.append(f'    sk(o, "{lb}.packed_groups_"); sk(o, "{lb}.packed_count_"); sk(o, "{lb}.mask_active_");')

def _rt_uhparams(pc, pr, lb, C, R):
    C.append(f'  pi("{lb}.n", (long long)({pc}).size()); for (const auto& u : ({pc})) {{')
    R.append(f'    pi(o, "{lb}.n", ({pr}).len() as i64); for u in ({pr}).iter() {{')
    for f in ['r', 'tPeak', 'tBase', 'iaMax', 'iaRecov', 'iaInit']:
        C.append(f'    {{ std::vector<double> t; for (int m = 0; m < 12; ++m) for (int k = 0; k < 3; ++k) t.push_back(u.{f}[m][k]); vf("{lb}[].{f}", t); }}')
        R.append(f'        vf(o, "{lb}[].{f}", &u.{f}.iter().flat_map(|a| a.iter().copied()).collect::<Vec<f64>>());')
    C.append('  }'); R.append('    }')

def _rt_decay(pc, pr, lb, C, R):
    C.append(f'  pi("{lb}.n", (long long)({pc}).size()); for (const auto& tr : ({pc})) for (const auto& d : tr) {{')
    R.append(f'    pi(o, "{lb}.n", ({pr}).len() as i64); for tr in ({pr}).iter() {{ for d in tr.iter() {{')
    rt_dump_struct('hydrology/RDII.hpp', 'ExpDecayParams', 'd', 'd', lb + '[]', C, R)
    C.append('  }'); R.append('    } }')

def _rt_uh_map(pc, pr, lb, C, R):
    # iteration order IS the check: upstream walks this map in initHydraulics
    C.append(f'  pi("{lb}.n", (long long)({pc}).size()); for (const auto& kv : ({pc})) {{ ps("{lb}.key", kv.first); pi("{lb}.val", kv.second); }}')
    # the translation keeps insertion order; the helper replays it through libstdc++'s table
    R.append(f'    {{ let v = swmm6_rust::core::std_unordered_map::ci_map_iteration_order(e.rdii_solver().uh_name_index()); pi(o, "{lb}.n", v.len() as i64); for (k, i) in &v {{ ps(o, "{lb}.key", k); pi(o, "{lb}.val", *i as i64); }} }}')

def _rt_expr_index(pc, pr, lb, C, R):
    # looked up only, never iterated: printed sorted by its case-folded key
    C.append(f'  {{ std::vector<std::pair<std::string,int>> t; for (const auto& kv : ({pc})) {{ std::string k = kv.first; for (auto& ch : k) ch = (ch >= 0x61 && ch <= 0x7a) ? (char)(ch & ~32) : ch; t.emplace_back(k, kv.second); }} std::sort(t.begin(), t.end()); pi("{lb}.n", (long long)t.size()); for (auto& kv : t) {{ ps("{lb}.key", kv.first); pi("{lb}.val", kv.second); }} }}')
    R.append(f'    {{ let mut t: Vec<(Vec<u8>, i32)> = ({pr}).iter().map(|(k, v)| (k.clone(), *v)).collect(); t.sort(); pi(o, "{lb}.n", t.len() as i64); for (k, v) in &t {{ ps(o, "{lb}.key", k); pi(o, "{lb}.val", *v as i64); }} }}')

def _rt_expressions(pc, pr, lb, C, R):
    C.append(f'  pi("{lb}.n", (long long)({pc}).size()); for (const auto& x : ({pc})) {{ pi("{lb}[].valid", x.valid); pi("{lb}[].postfix.n", (long long)x.postfix.size()); for (const auto& t : x.postfix) {{ pi("{lb}[].tok.type", (long long)t.type); pf("{lb}[].tok.value", t.value); pi("{lb}[].tok.var_idx", t.var_idx); ps("{lb}[].tok.var_name", t.var_name); }} }}')
    R.append(f'    pi(o, "{lb}.n", ({pr}).len() as i64); for x in ({pr}).iter() {{ pi(o, "{lb}[].valid", x.valid as i64); pi(o, "{lb}[].postfix.n", x.postfix.len() as i64); for t in &x.postfix {{ pi(o, "{lb}[].tok.type", t.type_ as i64); pf(o, "{lb}[].tok.value", t.value); pi(o, "{lb}[].tok.var_idx", t.var_idx as i64); ps(o, "{lb}[].tok.var_name", &t.var_name); }} }}')

def _rt_ptr(rust_flag):
    def f(pc, pr, lb, C, R):
        C.append(f'  pi("{lb}", ({pc}) != nullptr);'); R.append(f'    pi(o, "{lb}", ({rust_flag}) as i64);')
    return f

def _rt_opt(pc, pr, lb, C, R):
    C.append(f'  pi("{lb}", (bool)({pc}));'); R.append(f'    pi(o, "{lb}", ({pr}).is_some() as i64);')

# The CPU FV backend behind Router::fv_solver_ (a unique_ptr<INetworkSolver>): the
# members ExplicitFvSolver::initialize and the Router's init-time setters write, and (from
# 'cell_qlat_' on in each list, 5 Nov 2026) the members only the step writes. Not printed: the
# RK2 snapshots and the species FCT scratch (their paths are refused before the step), the
# PressurizedHeadSolver press_, and fireFaces' / limitPositivity's function-static scratch.
RT_FV_SOLVER = {
    'vf': ['f_mass_', 'f_mom_', 'f_sstar_', 'f_corr_l_', 'f_corr_r_', 'f_phi_l_', 'f_phi_r_', 'f_phi_flux_',
           'cell_slope_', 'cell_eta_slope_', 'cell_u_slope_', 'cell_eta_', 'cell_u_', 'cell_ah_', 'cell_t_',
           'cell_i1_', 'cell_q_int_', 'uf_grad_', 'withdrawal_rates_', 'unmet_withdrawal_', 'save_unmet_',
           'rk_unmet_', 'node_exch_', 'node_in_', 'node_out_', 'flood_vol_', 'node_qdummy_', 'dummy_vol_',
           'relax_beta_', 'acc_a_', 'acc_q_', 'acc_nvol_', 'node_carry_', 'node_vfull_', 'face_a_rim_',
           'face_lat_spill_',
           'cell_qlat_', 'node_qstruct_', 'save_dummy_', 'save_cell_a_', 'save_cell_q_', 'save_cell_phi_',
           'save_node_vol_', 'save_node_head_', 'save_exch_', 'save_in_', 'save_out_', 'save_flood_',
           'save_qint_', 'save_carry_'],
    'vi': ['cell_ho2_', 'tpa_scratch_', 'inlet_control_', 'cell_tier_', 'node_tier_', 'face_tier_', 'cell_active_',
           'active_faces_', 'node_alg_', 'node_pass_static_', 'node_pass_', 'node_pub_static_', 'node_pub_',
           'node_lat_div_', 'save_tpa_', 'halo_prev_', 'live_stamp_', 'node_solve_list_',
           'node_heavy_', 'node_light_'],
    'vvi': ['cells_by_tier_', 'faces_by_tier_', 'nodes_by_tier_', 'due_f_upto_', 'due_c_upto_', 'due_n_upto_'],
    'i': ['hllc_', 'tpa_', 'mixed_wave_', 'mw_rh_', 'deg1_guard_slot_', 'census_pass_far_', 'lts_fit_',
          'water_tracking_', 'conservative_node_balance_', 'lts_tiers_', 'lists_valid_', 'since_rebuild_',
          'census_count_', 'lts_valid_', 'lts_countdown_', 'total_steps_', 'total_flux_', 'active_n_',
          'n_macro_cycles_', 'n_macro_rejected_', 'state_clean_', 'acc_dirty_',
          'live_gen_', 'all_faces_live_', 'div_link_', 'last_nsteps_', 'dt_argmin_pressurized_', 'dt_argmin_band_',
          'dt_argmin_free_', 'dt_argmin_node_', 'press_step_', 'pressure_candidate_'],
    'f': ['dt_cache_', 'dt_census_', 'lts_dt0_', 'min_h_', 'sim_time_', 'active_sum_', 'active_min_', 'active_max_',
          'div_value_', 'last_h_', 'suggested_h_'],
}

def _rt_fv_solver(pc, pr, lb, C, R):
    C.append(f'  if (auto* s = dynamic_cast<const openswmm::fv::ExplicitFvSolver*>(({pc}).get())) {{ pi("{lb}", 1);')
    R.append(f'    match &({pr}) {{ None => pi(o, "{lb}", 0), Some(s) => {{ pi(o, "{lb}", 1);')
    C.append(f'  pi("{lb}.attached", s->mesh_ != nullptr); pi("{lb}.step_sink", s->step_sink_ != nullptr);')
    R.append(f'    pi(o, "{lb}.attached", s.attached_ as i64); pi(o, "{lb}.step_sink", s.step_sink_set_ as i64);')
    rt_dump_struct('hydraulics/fv/FvOptions.hpp', 'FvOptions', 's->opts_', 's.opts_', lb + '.opts_', C, R)
    for f in RT_FV_SOLVER['vf']:
        emit(f's->{f}', f's.{f}', f'{lb}.{f}', ('v', ('f', 'double')), C, R)
    for f in RT_FV_SOLVER['vi']:
        emit(f's->{f}', f's.{f}', f'{lb}.{f}', ('v', ('i', 'int')), C, R)
    for f in RT_FV_SOLVER['vvi']:
        emit(f's->{f}', f's.{f}', f'{lb}.{f}', ('v', ('v', ('i', 'int'))), C, R)
    for f in RT_FV_SOLVER['i']:
        emit(f's->{f}', f's.{f}', f'{lb}.{f}', ('i', 'int'), C, R)
    for f in RT_FV_SOLVER['f']:
        emit(f's->{f}', f's.{f}', f'{lb}.{f}', ('f', 'double'), C, R)
    emit('s->tier_occupancy_', 's.tier_occupancy_', f'{lb}.tier_occupancy_', ('arr', ('i', 'long'), 'kMaxLtsTiers'), C, R)
    for f, st in [('f_state_l_', 'FaceState'), ('f_state_r_', 'FaceState'), ('f_flux_', 'FaceFlux')]:
        C.append(f'  pi("{lb}.{f}.n", (long long)s->{f}.size()); for (const auto& el : s->{f}) {{')
        R.append(f'    pi(o, "{lb}.{f}.n", s.{f}.len() as i64); for el in s.{f}.iter() {{')
        rt_dump_struct('hydraulics/fv/FvKernels.hpp', st, 'el', '(*el)', f'{lb}.{f}[]', C, R)
        C.append('  }'); R.append('    }')
    for f in ['water_', 'save_water_', 'rk_water_']:
        rt_dump_struct('hydraulics/fv/AcceptedWaterVolumes.hpp', 'AcceptedWaterVolumes', f's->{f}', f's.{f}', f'{lb}.{f}', C, R)
    C.append(f'  ps("{lb}.div_what_", s->div_what_ ? std::string(s->div_what_) : std::string());')
    R.append(f'    ps(o, "{lb}.div_what_", &s.div_what_.unwrap_or("").as_bytes().to_vec());')
    C.append(f'  }} else pi("{lb}", 0);')
    R.append('    } }')

# (struct, member) -> emitter; replaces the generic handling of that member
RT_CUSTOM = {
    ('Router', 'groups_'): _rt_xsect_groups,
    ('Router', 'fv_solver_'): _rt_fv_solver,
    ('RDIISolver', 'uh_params'): _rt_uhparams,
    ('RDIISolver', 'decay_params'): _rt_decay,
    ('RDIISolver', 'uh_name_to_idx_'): _rt_uh_map,
    ('ControlEngine', 'expr_index_'): _rt_expr_index,
    ('ControlEngine', 'expressions_'): _rt_expressions,
}
# members the header parser does not list (raw pointers): (struct) -> [(member, emitter)]
RT_EXTRA = {
    'Router': [('checkpoint_context_', 'checkpoint_context_attached_'), ('structures_', 'structures_attached_')],
    'DWSolver': [('groups_', 'groups_attached_')],
    'InterfaceManager': [('infile_', None), ('outfile_', None)],
}

def rt_dump_struct(header, struct, path_c, path_r, label, C, R):
    """dump_struct with the RT_CUSTOM / RT_EXTRA hooks, recursing into itself."""
    for name, ty, arr in fields(header, struct):
        pc, pr, lb = f'{path_c}.{name}', f'{path_r}.{rname(name)}', f'{label}.{name}'
        if (struct, name) in RT_CUSTOM:
            RT_CUSTOM[(struct, name)](pc, pr, lb, C, R)
            continue
        k = kind(ty)
        if arr:
            if k and k[0] in 'fi':
                C.append(f'  for (int q = 0; q < (int)(sizeof({pc})/sizeof({pc}[0])); ++q) {{')
                R.append(f'    for q in 0..({pr}).len() {{')
                emit(f'{pc}[q]', f'{pr}[q]', lb + '[]', k, C, R)
                C.append('  }'); R.append('    }')
            else:
                SKIPS.append(lb); C.append(f'  sk("{lb}");'); R.append(f'    sk(o, "{lb}");')
            continue
        t = ty.replace('const ', '').strip()
        em = re.fullmatch(r'std::vector<(?:\w+::)*(\w+)>', t)
        if k is None and em and em.group(1) in RT_ELEM_STRUCTS:
            C.append(f'  pi("{lb}.n", (long long)({pc}).size());'); R.append(f'    pi(o, "{lb}.n", ({pr}).len() as i64);')
            C.append(f'  for (auto& el : {pc}) {{'); R.append(f'    for el in ({pr}).iter() {{')
            rt_dump_struct(RT_ELEM_STRUCTS[em.group(1)], em.group(1), 'el', '(*el)', lb + '[]', C, R)
            C.append('  }'); R.append('    }')
            continue
        if k is not None and k[0] == 'struct':
            rt_dump_struct(k[1], k[2], pc, pr, lb, C, R)
            continue
        if k is None:
            SKIPS.append(f'{lb} ({ty})'); C.append(f'  sk("{lb}");'); R.append(f'    sk(o, "{lb}");')
            continue
        emit(pc, pr, lb, k, C, R)
    for name, flag in RT_EXTRA.get(struct, []):
        lb = f'{label}.{name}'
        if flag is None:   # FILE*: open or not
            C.append(f'  pi("{lb}", ({path_c}.{name}) != nullptr);'); R.append(f'    pi(o, "{lb}", ({path_r}.{name}).is_some() as i64);')
        else:
            C.append(f'  pi("{lb}", ({path_c}.{name}) != nullptr);'); R.append(f'    pi(o, "{lb}", ({path_r}.{flag}) as i64);')

# (header, struct, C++ path, Rust path, label)
RT_MEMBERS = [
    ('hydraulics/Routing.hpp', 'Router', 'e.router_', 'e.router()', 'router'),
    ('controls/Controls.hpp', 'ControlEngine', 'e.controls_', 'e.control_engine()', 'controls'),
    ('hydrology/Inflow.hpp', 'InflowSolver', 'e.inflow_', 'e.inflow_solver()', 'inflow'),
    ('hydrology/RDII.hpp', 'RDIISolver', 'e.rdii_', 'e.rdii_solver()', 'rdii'),
    ('hydraulics/Exfiltration.hpp', 'ExfilSolver', 'e.exfil_', 'e.exfil_solver()', 'exfil'),
    ('hydraulics/Inlet.hpp', 'InletSolver', 'e.inlet_', 'e.inlet_solver()', 'inlet'),
    ('hydraulics/HydStructures.hpp', 'StructureSolver', 'e.hydstruct_', 'e.structure_solver()', 'hydstruct'),
    ('core/InterfaceFile.hpp', 'InterfaceManager', 'e.iface_', 'e.interface_manager()', 'iface'),
]

def router_dumper():
    saved = dict(STRUCTS)
    STRUCTS.update(RT_STRUCTS)
    nskip = len(SKIPS)
    C = ['// generated by tools/ctxdump.py: prints the state SWMMEngine::initHydraulics builds; the Rust twin prints the same',
         '// include after `#define private public` and core/SWMMEngine.hpp (see gate/diff/inithyd.cpp)',
         '#pragma once', '#include "ctxdump.hpp"', 'namespace ctxdump {',
         'inline void rtdump(const openswmm::SWMMEngine& e) {']
    R = ['// generated by tools/ctxdump.py: Rust twin of gate/diff/rtdump.hpp', '#![allow(dead_code, unused_parens, clippy::all)]',
         'use std::io::Write;', 'use swmm6_rust::core::engine::SwmmEngine;',
         'fn pf(o: &mut Vec<u8>, l: &str, v: f64) { writeln!(o, "{} {:016x}", l, v.to_bits()).unwrap(); }',
         'fn pi(o: &mut Vec<u8>, l: &str, v: i64) { writeln!(o, "{} {}", l, v).unwrap(); }',
         'fn ps(o: &mut Vec<u8>, l: &str, s: &[u8]) { write!(o, "{} {}:", l, s.len()).unwrap(); o.extend_from_slice(s); o.push(b\'\\n\'); }',
         'fn vf(o: &mut Vec<u8>, l: &str, v: &[f64]) { writeln!(o, "{} n={}", l, v.len()).unwrap(); for x in v { write!(o, " {:016x}", x.to_bits()).unwrap(); } if !v.is_empty() { o.push(b\'\\n\'); } }',
         'fn vi(o: &mut Vec<u8>, l: &str, v: &[i64]) { writeln!(o, "{} n={}", l, v.len()).unwrap(); for x in v { write!(o, " {}", x).unwrap(); } if !v.is_empty() { o.push(b\'\\n\'); } }',
         'fn vs(o: &mut Vec<u8>, l: &str, v: &[Vec<u8>]) { writeln!(o, "{} n={}", l, v.len()).unwrap(); for s in v { write!(o, " {}:", s.len()).unwrap(); o.extend_from_slice(s); } if !v.is_empty() { o.push(b\'\\n\'); } }',
         'fn sk(o: &mut Vec<u8>, l: &str) { writeln!(o, "skip {}", l).unwrap(); }',
         'fn pn(o: &mut Vec<u8>, l: &str) { writeln!(o, "{} none", l).unwrap(); }',
         'pub fn rtdump(e: &SwmmEngine) -> Vec<u8> {', '    let mut o: Vec<u8> = Vec::new();', '    let o = &mut o;']
    # the context fields initHydraulics writes (routing_stats, water_age_state) are in
    # ctxdump::dump (CTX_MEMBERS), which the harness prints first
    for header, struct, pc, pr, label in RT_MEMBERS:
        rt_dump_struct(header, struct, pc, pr, label, C, R)
    C.append('  pi("engine.next_event_", e.next_event_); pi("engine.between_events_", e.between_events_);')
    R.append('    pi(o, "engine.next_event_", e.next_event() as i64); pi(o, "engine.between_events_", e.between_events() as i64);')
    C.append('}'); C.append('} // namespace ctxdump')
    R.append('    std::mem::take(o)'); R.append('}')
    STRUCTS.clear(); STRUCTS.update(saved)
    rskips = SKIPS[nskip:]
    del SKIPS[nskip:]
    open(f'{ROOT}/gate/diff/rtdump.hpp', 'w').write('\n'.join(C) + '\n')
    open(f'{ROOT}/crates/engine/examples/common/rtdump.rs', 'w').write('\n'.join(R) + '\n')
    print(f'rtdump: {len(rskips)} skipped:', '; '.join(rskips))

# ---------------------------------------------------------------------------
# SWMMEngine::initHydrology: the hydrology solver members and the context members it
# fills that the main dump does not cover (gate/diff/inithydro.cpp, examples/inithydro.rs).
# The C++ side is compiled with `private`/`protected` defined as `public`, so it reads
# the members directly; the Rust side reads them through SwmmEngine's gate accessors.
# Struct types are looked up in HYDRO_STRUCTS, merged into STRUCTS only while this dumper
# is generated, so the main dump is unchanged ('Token' there is quality/Treatment's).
HYDRO_STRUCTS = {
    'RadiationState': ('hydrology/Radiation.hpp', 'RadiationState'),
    'MovingAvg7': ('hydrology/Climate.hpp', 'MovingAvg7'),
    'HeldState': ('hydrology/surface/UebHeld.hpp', 'HeldState'),
    'LumpedState': ('hydrology/surface/UebLumped.hpp', 'LumpedState'),
    'LumpedFluxes': ('hydrology/surface/UebLumped.hpp', 'LumpedFluxes'),
    'State': ('hydrology/surface/UebSnow.hpp', 'State'),
    'Fluxes': ('hydrology/surface/UebSnow.hpp', 'Fluxes'),
    'ThermalSample': ('hydrology/surface/UebSnow.hpp', 'ThermalSample'),
    'climate::SnowDailyWeather': ('hydrology/HourlyMeteorology.hpp', 'SnowDailyWeather'),
    'RunoffSoA': ('hydrology/Runoff.hpp', 'RunoffSoA'),
    'CompletedSubareaWater': ('hydrology/Runoff.hpp', 'CompletedSubareaWater'),
    'surface::InfilBank': ('hydrology/surface/InfilBank.hpp', 'InfilBank'),
    'HortonState': ('hydrology/Infiltration.hpp', 'HortonState'),
    'CurveNumState': ('hydrology/Infiltration.hpp', 'CurveNumState'),
    'SnowSoA': ('hydrology/surface/DegreeDaySnow.hpp', 'SnowSoA'),
    'GWSoA': ('hydrology/Groundwater.hpp', 'GWSoA'),
    'mathexpr::Expression': ('math/MathExpr.hpp', 'Expression'),
    'Token': ('math/MathExpr.hpp', 'Token'),
    'LIDGroupSoA': ('hydrology/LID.hpp', 'LIDGroupSoA'),
    'gage::GageState': ('hydrology/Gage.hpp', 'GageState'),
    'ExpDecayParams': ('hydrology/RDII.hpp', 'ExpDecayParams'),
    'UHResponseData': ('hydrology/RDII.hpp', 'UHResponseData'),
    'RDIIGroupSoA': ('hydrology/RDII.hpp', 'RDIIGroupSoA'),
}
# The context members initHydrology fills (climate_state, lid_layer_state, snowpacks.runtime,
# snowpacks.daily_weather) are in ctxdump::dump, which the harness prints first.
HYDRO_ENGINE = [
    ('hydrology/Runoff.hpp', 'RunoffSolver', 'e.runoff_', 'e.runoff_solver()', 'runoff_'),
    ('hydrology/surface/DegreeDaySnow.hpp', 'DegreeDaySnow', 'e.snow_', 'e.snow_solver()', 'snow_'),
    ('hydrology/Groundwater.hpp', 'GWSolver', 'e.groundwater_', 'e.groundwater_solver()', 'groundwater_'),
    ('hydrology/LID.hpp', 'LIDSolver', 'e.lid_', 'e.lid_solver()', 'lid_'),
    # rdii_ is initialized by initHydraulics upstream; inithydro's "rdii" mode runs it
    ('hydrology/RDII.hpp', 'RDIISolver', 'e.rdii_', 'e.rdii_solver()', 'rdii_'),
]
HYDRO_VECTORS = ['gw_frac_perv_', 'gw_perv_evap_', 'gw_lid_perv_area_', 'gw_lid_infil_vol_', 'gw_lid_perv_evap_vol_',
                 'gw_2d_recharge_vol_', 'gw_infil_with_2d_', 'gw_infil_rate_', 'snow_rain_', 'snow_snow_',
                 'gw_q_interp_', 'old_gw_rate_']

def hydro_dumper():
    saved = dict(STRUCTS)
    STRUCTS.update(HYDRO_STRUCTS)
    nskip = len(SKIPS)
    C = ['// generated by tools/ctxdump.py: the hydrology members SWMMEngine::initHydrology fills;',
         '// include after core/SWMMEngine.hpp compiled with private/protected defined as public', '#pragma once',
         '#include <cstdio>', '#include "ctxdump.hpp"', 'namespace ctxdump {',
         'inline void hydrodump(openswmm::SWMMEngine& e) {']
    R = ['// generated by tools/ctxdump.py: Rust twin of gate/diff/hydrodump.hpp', '#![allow(dead_code, unused_parens, clippy::all)]',
         'use std::io::Write;', 'use swmm6_rust::core::engine::SwmmEngine;',
         'fn pf(o: &mut Vec<u8>, l: &str, v: f64) { writeln!(o, "{} {:016x}", l, v.to_bits()).unwrap(); }',
         'fn pi(o: &mut Vec<u8>, l: &str, v: i64) { writeln!(o, "{} {}", l, v).unwrap(); }',
         'fn ps(o: &mut Vec<u8>, l: &str, s: &[u8]) { write!(o, "{} {}:", l, s.len()).unwrap(); o.extend_from_slice(s); o.push(b\'\\n\'); }',
         'fn vf(o: &mut Vec<u8>, l: &str, v: &[f64]) { writeln!(o, "{} n={}", l, v.len()).unwrap(); for x in v { write!(o, " {:016x}", x.to_bits()).unwrap(); } if !v.is_empty() { o.push(b\'\\n\'); } }',
         'fn vi(o: &mut Vec<u8>, l: &str, v: &[i64]) { writeln!(o, "{} n={}", l, v.len()).unwrap(); for x in v { write!(o, " {}", x).unwrap(); } if !v.is_empty() { o.push(b\'\\n\'); } }',
         'fn vs(o: &mut Vec<u8>, l: &str, v: &[Vec<u8>]) { writeln!(o, "{} n={}", l, v.len()).unwrap(); for s in v { write!(o, " {}:", s.len()).unwrap(); o.extend_from_slice(s); } if !v.is_empty() { o.push(b\'\\n\'); } }',
         'fn sk(o: &mut Vec<u8>, l: &str) { writeln!(o, "skip {}", l).unwrap(); }',
         'fn pn(o: &mut Vec<u8>, l: &str) { writeln!(o, "{} none", l).unwrap(); }',
         'pub fn hydrodump(e: &mut SwmmEngine) -> Vec<u8> {', '    let mut o: Vec<u8> = Vec::new();', '    let o = &mut o;']
    for header, struct, pc, pr, label in HYDRO_ENGINE:
        dump_struct(header, struct, pc, pr, label, C, R)
    # RDIISolver: uh_params holds 2-D C arrays and uh_name_to_idx_ is an unordered_map
    # (sorted by key here); the generic walk skips both
    C.append('  for (auto& u : e.rdii_.uh_params) for (int m = 0; m < 12; ++m) for (int k = 0; k < 3; ++k) { pf("rdii_.uh_params[].r", u.r[m][k]); pf("rdii_.uh_params[].tPeak", u.tPeak[m][k]); pf("rdii_.uh_params[].tBase", u.tBase[m][k]); pf("rdii_.uh_params[].iaMax", u.iaMax[m][k]); pf("rdii_.uh_params[].iaRecov", u.iaRecov[m][k]); pf("rdii_.uh_params[].iaInit", u.iaInit[m][k]); }')
    R.append('    for u in &e.rdii_solver().uh_params { for m in 0..12 { for k in 0..3 { pf(o, "rdii_.uh_params[].r", u.r[m][k]); pf(o, "rdii_.uh_params[].tPeak", u.tPeak[m][k]); pf(o, "rdii_.uh_params[].tBase", u.tBase[m][k]); pf(o, "rdii_.uh_params[].iaMax", u.iaMax[m][k]); pf(o, "rdii_.uh_params[].iaRecov", u.iaRecov[m][k]); pf(o, "rdii_.uh_params[].iaInit", u.iaInit[m][k]); } } }')
    C.append('  { std::vector<std::pair<std::string, int>> m(e.rdii_.uh_name_to_idx_.begin(), e.rdii_.uh_name_to_idx_.end()); std::sort(m.begin(), m.end()); for (auto& kv : m) { ps("rdii_.uh_name_to_idx_.key", kv.first); pi("rdii_.uh_name_to_idx_.val", kv.second); } }')
    R.append('    { let mut m = e.rdii_solver().uh_name_to_idx_.clone(); m.sort(); for kv in &m { ps(o, "rdii_.uh_name_to_idx_.key", &kv.0); pi(o, "rdii_.uh_name_to_idx_.val", kv.1 as i64); } }')
    # LIDSolver::rpt_files_ holds std::ofstream objects: dumped as its size and its keys
    # (group, unit) in map order
    C.append('  pi("lid_.rpt_files_.size", (long long)e.lid_.rpt_files_.size());')
    C.append('  for (auto& kv : e.lid_.rpt_files_) { pi("lid_.rpt_files_.key.first", kv.first.first); pi("lid_.rpt_files_.key.second", kv.first.second); }')
    R.append('    pi(o, "lid_.rpt_files_.size", e.lid_solver().rpt_files_.len() as i64);')
    R.append('    for (k, v) in e.lid_solver().rpt_files_.iter() { pi(o, "lid_.rpt_files_.key.first", k.0 as i64); pi(o, "lid_.rpt_files_.key.second", k.1 as i64); let _ = v; }')
    emit('e.gage_states_', 'e.gage_states()', 'gage_states_', kind('std::vector<gage::GageState>'), C, R)
    C.append('  pi("last_melt_doy_", e.last_melt_doy_);'); R.append('    pi(o, "last_melt_doy_", e.last_melt_doy() as i64);')
    R.append('    { let v = e.hydrology_vectors();')
    for k, n in enumerate(HYDRO_VECTORS):
        C.append(f'  vf("{n}", e.{n});'); R.append(f'    vf(o, "{n}", v[{k}]);')
    R.append('    }')
    C.append('  { std::vector<long long> t; for (auto x : e.gw_q_node_) t.push_back(x); vi("gw_q_node_", t); }')
    R.append('    vi(o, "gw_q_node_", &e.gw_q_node().iter().map(|x| *x as i64).collect::<Vec<i64>>());')
    # the climate file reader: its FILE* as open/position, then every member
    C.append('  { auto& cf = e.climate_file_; pi("climate_file_.open", cf.file_ != nullptr); pi("climate_file_.pos", cf.file_ ? (long long)std::ftell(cf.file_) : -1LL);'
             ' pi("climate_file_.format_", (long long)cf.format_); pi("climate_file_.unit_system_", cf.unit_system_);'
             ' for (int v = 0; v < 4; ++v) for (int d = 0; d < 32; ++d) pf("climate_file_.file_data_", cf.file_data_[v][d]);'
             ' pi("climate_file_.buf_year_", cf.buf_year_); pi("climate_file_.buf_month_", cf.buf_month_); pf("climate_file_.day_offset_", cf.day_offset_);'
             ' for (int i = 0; i < 4; ++i) pi("climate_file_.field_pos_", cf.field_pos_[i]); pi("climate_file_.date_field_pos_", cf.date_field_pos_);'
             ' pi("climate_file_.wind_type_", cf.wind_type_); pi("climate_file_.temp_units_", cf.temp_units_);'
             ' ps("climate_file_.saved_line_", std::string(cf.saved_line_)); pi("climate_file_.has_saved_line_", cf.has_saved_line_); }')
    R.append('    { let cf = e.climate_file_reader(); pi(o, "climate_file_.open", cf.file_.is_some() as i64); let pos = match cf.file_.as_mut() { Some(f) => f.ftell(), None => -1 }; pi(o, "climate_file_.pos", pos);'
             ' pi(o, "climate_file_.format_", cf.format_ as i64); pi(o, "climate_file_.unit_system_", cf.unit_system_ as i64);'
             ' for v in 0..4 { for d in 0..32 { pf(o, "climate_file_.file_data_", cf.file_data_[v][d]); } }'
             ' pi(o, "climate_file_.buf_year_", cf.buf_year_ as i64); pi(o, "climate_file_.buf_month_", cf.buf_month_ as i64); pf(o, "climate_file_.day_offset_", cf.day_offset_);'
             ' for i in 0..4 { pi(o, "climate_file_.field_pos_", cf.field_pos_[i] as i64); } pi(o, "climate_file_.date_field_pos_", cf.date_field_pos_ as i64);'
             ' pi(o, "climate_file_.wind_type_", cf.wind_type_ as i64); pi(o, "climate_file_.temp_units_", cf.temp_units_ as i64);'
             ' ps(o, "climate_file_.saved_line_", &cf.saved_line_); pi(o, "climate_file_.has_saved_line_", cf.has_saved_line_ as i64); }')
    C.append('}'); C.append('} // namespace ctxdump')
    R.append('    std::mem::take(o)'); R.append('}')
    STRUCTS.clear(); STRUCTS.update(saved)
    hskips = SKIPS[nskip:]
    del SKIPS[nskip:]
    open(f'{ROOT}/gate/diff/hydrodump.hpp', 'w').write('\n'.join(C) + '\n')
    open(f'{ROOT}/crates/engine/examples/common/hydrodump.rs', 'w').write('\n'.join(R) + '\n')
    print(f'hydrodump: {len(hskips)} skipped:', '; '.join(hskips))

# SWMMEngine members that initQuality sets (and hotstart_deferred_, which it reads), dumped by
# gate/diff/initqual.cpp. The C++ side reads them directly: the includer must define
# private/protected as public before core/SWMMEngine.hpp and then define CTXDUMP_QUALITY.
# (header, class, C++ path, Rust path, label)
QUALITY_MEMBERS = [
    ('quality/QualityRouting.hpp', 'QualitySolver', 'e.quality_', 'e.quality()', 'quality'),
    ('quality/Landuse.hpp', 'LanduseSolver', 'e.landuse_solver_', 'e.landuse_solver()', 'landuse_solver'),
    ('quality/Landuse.hpp', 'SurfaceQualitySoA', 'e.surface_quality_', 'e.surface_quality()', 'surface_quality'),
]

def quality_dumper():
    C = ['#ifdef CTXDUMP_QUALITY', 'namespace ctxdump {', 'inline void qualdump(openswmm::SWMMEngine& e) {']
    R = ['pub fn qualdump(e: &SwmmEngine) -> Vec<u8> {', '    let mut o: Vec<u8> = Vec::new();', '    let o = &mut o;']
    C.append('  pi("engine.hotstart_deferred_", (long long)e.hotstart_deferred_);')
    R.append('    pi(o, "engine.hotstart_deferred_", e.hotstart_deferred() as i64);')
    for header, cls, pc, pr, label in QUALITY_MEMBERS:
        dump_struct(header, cls, pc, pr, label, C, R)
    C += ['}', '} // namespace ctxdump', '#endif // CTXDUMP_QUALITY']
    R += ['    std::mem::take(o)', '}']
    return C, R

# ---------------------------------------------------------------------------
# The run lifecycle's SWMMEngine members (initialize, start, step) and the run-time members
# of what start() opens (runoff and RDII interface files, the built-in output and report
# plugins): gate/diff/lifedump.hpp and examples/common/lifedump.rs, used by gate/diff/ctxstep
# and examples/ctxstep. The C++ side is compiled with private/protected defined as public.
# (name, C++ type); the types are given here because several of these members share a
# declaration upstream (`bool step_active_ = false, runtime_pipeline_ = false;`).
LIFE_MEMBERS = [
    ('plugins_prepare_attempted_', 'bool'), ('hotstart_save_done_', 'std::vector<char>'),
    ('hotstart_runoff_', 'std::vector<double>'), ('step_active_', 'bool'), ('runtime_pipeline_', 'bool'),
    ('checkpoint_enabled_', 'bool'), ('checkpoint_ready_', 'bool'), ('checkpoint_output_', 'bool'),
    ('checkpoint_model_', 'std::string'), ('pending_hotstart_checkpoint_', 'std::string'),
    ('hotstart_checkpoint_error_', 'std::string'), ('exchange_target_', 'double'),
    ('old_runoff_time_', 'double'), ('new_runoff_time_', 'double'), ('old_runoff_ms_', 'double'),
    ('prev_runoff_step_sec_', 'double'), ('new_runoff_ms_', 'double'), ('new_rule_time_ms_', 'double'),
    ('prev_lid_infil_vol_', 'double'), ('prev_lid_evap_vol_', 'double'), ('wet_q_interp_', 'std::vector<double>'),
    ('lid_drain_q_interp_', 'std::vector<double>'), ('lid_drain_q_units_', 'std::vector<std::pair<int, double>>'),
    ('lid_vlidin_vol_', 'std::vector<double>'), ('lid_qsurf_cfs_', 'std::vector<double>'),
    ('has_runoff_', 'bool'), ('has_snow_', 'bool'), ('save_results_', 'int'), ('do_routing_', 'bool'),
    ('last_step_steady_', 'bool'), ('routing_prev_rates_', 'std::array<double, 11>'),
]
LIFE_AVG = ['node_depth', 'node_head', 'node_volume', 'node_lat_inflow', 'node_total_inflow', 'node_overflow',
            'link_flow', 'link_depth', 'link_velocity', 'link_volume', 'link_capacity']

def life_dumper():
    C = ['// generated by tools/ctxdump.py: the run lifecycle members of SWMMEngine (initialize, start, step);',
         '// include after core/SWMMEngine.hpp compiled with private/protected defined as public', '#pragma once',
         '#include <cstdio>', '#include "ctxdump.hpp"', '#include "plugins/DefaultOutputPlugin.hpp"',
         '#include "plugins/DefaultReportPlugin.hpp"', 'namespace ctxdump {',
         'inline void lifedump(openswmm::SWMMEngine& e) {']
    R = ['// generated by tools/ctxdump.py: Rust twin of gate/diff/lifedump.hpp', '#![allow(dead_code, unused_parens, clippy::all)]',
         'use std::io::Write;', 'use swmm6_rust::core::engine::SwmmEngine;',
         'fn pf(o: &mut Vec<u8>, l: &str, v: f64) { writeln!(o, "{} {:016x}", l, v.to_bits()).unwrap(); }',
         'fn pi(o: &mut Vec<u8>, l: &str, v: i64) { writeln!(o, "{} {}", l, v).unwrap(); }',
         'fn ps(o: &mut Vec<u8>, l: &str, s: &[u8]) { write!(o, "{} {}:", l, s.len()).unwrap(); o.extend_from_slice(s); o.push(b\'\\n\'); }',
         'fn vf(o: &mut Vec<u8>, l: &str, v: &[f64]) { writeln!(o, "{} n={}", l, v.len()).unwrap(); for x in v { write!(o, " {:016x}", x.to_bits()).unwrap(); } if !v.is_empty() { o.push(b\'\\n\'); } }',
         'fn vi(o: &mut Vec<u8>, l: &str, v: &[i64]) { writeln!(o, "{} n={}", l, v.len()).unwrap(); for x in v { write!(o, " {}", x).unwrap(); } if !v.is_empty() { o.push(b\'\\n\'); } }',
         'fn vs(o: &mut Vec<u8>, l: &str, v: &[Vec<u8>]) { writeln!(o, "{} n={}", l, v.len()).unwrap(); for s in v { write!(o, " {}:", s.len()).unwrap(); o.extend_from_slice(s); } if !v.is_empty() { o.push(b\'\\n\'); } }',
         'fn pn(o: &mut Vec<u8>, l: &str) { writeln!(o, "{} none", l).unwrap(); }',
         'pub fn lifedump(e: &mut SwmmEngine) -> Vec<u8> {', '    let mut o: Vec<u8> = Vec::new();', '    let o = &mut o;']
    for n, t in LIFE_MEMBERS:
        emit(f'e.{n}', f'e.{n}', f'engine.{n}', kind(t), C, R)
    # LegacySpeciesTrailer hotstart_species_
    C.append('  { auto& h = e.hotstart_species_; vs("engine.hotstart_species_.species", h.species); pi("engine.hotstart_species_.n_subcatch", h.n_subcatch); pi("engine.hotstart_species_.n_landuses", h.n_landuses); vf("engine.hotstart_species_.washoff", h.washoff); vf("engine.hotstart_species_.ponded", h.ponded); vf("engine.hotstart_species_.buildup", h.buildup); }')
    R.append('    { let h = &e.hotstart_species_; vs(o, "engine.hotstart_species_.species", &h.species); pi(o, "engine.hotstart_species_.n_subcatch", h.n_subcatch as i64); pi(o, "engine.hotstart_species_.n_landuses", h.n_landuses as i64); vf(o, "engine.hotstart_species_.washoff", &h.washoff); vf(o, "engine.hotstart_species_.ponded", &h.ponded); vf(o, "engine.hotstart_species_.buildup", &h.buildup); }')
    # AvgAccumulator avg_ (float sums)
    for f in LIFE_AVG:
        emit(f'e.avg_.{f}', f'e.avg_.{f}', f'engine.avg_.{f}', kind('std::vector<float>'), C, R)
    C.append('  pi("engine.avg_.n_steps", e.avg_.n_steps);'); R.append('    pi(o, "engine.avg_.n_steps", e.avg_.n_steps as i64);')
    # the XSectParams cache: its size and generation (the entries hold pointers)
    C.append('  pi("engine.xsp_cache_.size", (long long)e.xsp_cache_.size()); pi("engine.xsp_cache_gen_", (long long)e.xsp_cache_gen_);')
    R.append('    pi(o, "engine.xsp_cache_.size", e.xsp_cache_len_ as i64); pi(o, "engine.xsp_cache_gen_", e.xsp_cache_gen_ as i64);')
    # runoff interface file (a unique_ptr upstream; the FILE* as open + position)
    C.append('  if (e.runoff_iface_file_) { auto& f = *e.runoff_iface_file_; pi("engine.runoff_iface_file_.open", f.fp_ != nullptr); pi("engine.runoff_iface_file_.pos", f.fp_ ? (long long)std::ftell(f.fp_) : -1LL);'
             ' pi("engine.runoff_iface_file_.writing_", f.writing_); pi("engine.runoff_iface_file_.n_subcatch_", f.n_subcatch_); pi("engine.runoff_iface_file_.n_pollut_", f.n_pollut_);'
             ' pi("engine.runoff_iface_file_.n_results_", f.n_results_); pi("engine.runoff_iface_file_.flow_units_", f.flow_units_); pi("engine.runoff_iface_file_.unit_system_", f.unit_system_);'
             ' pi("engine.runoff_iface_file_.step_count_", f.step_count_); pi("engine.runoff_iface_file_.max_steps_pos_", f.max_steps_pos_);'
             ' { std::vector<double> t(f.buf_.begin(), f.buf_.end()); vf("engine.runoff_iface_file_.buf_", t); } } else pn("engine.runoff_iface_file_");')
    R.append('    match e.runoff_iface_file_.as_mut() { Some(f) => { pi(o, "engine.runoff_iface_file_.open", f.fp_.is_some() as i64); let pos = f.tell(); pi(o, "engine.runoff_iface_file_.pos", pos);'
             ' pi(o, "engine.runoff_iface_file_.writing_", f.writing_ as i64); pi(o, "engine.runoff_iface_file_.n_subcatch_", f.n_subcatch_ as i64); pi(o, "engine.runoff_iface_file_.n_pollut_", f.n_pollut_ as i64);'
             ' pi(o, "engine.runoff_iface_file_.n_results_", f.n_results_ as i64); pi(o, "engine.runoff_iface_file_.flow_units_", f.flow_units_ as i64); pi(o, "engine.runoff_iface_file_.unit_system_", f.unit_system_ as i64);'
             ' pi(o, "engine.runoff_iface_file_.step_count_", f.step_count_ as i64); pi(o, "engine.runoff_iface_file_.max_steps_pos_", f.max_steps_pos_);'
             ' vf(o, "engine.runoff_iface_file_.buf_", &f.buf_.iter().map(|x| *x as f64).collect::<Vec<f64>>()); } None => pn(o, "engine.runoff_iface_file_") }')
    # RDII interface file
    C.append('  { auto& f = e.rdii_iface_file_; pi("engine.rdii_iface_file_.open", f.fp_ != nullptr); pi("engine.rdii_iface_file_.pos", f.fp_ ? (long long)std::ftell(f.fp_) : -1LL);'
             ' pi("engine.rdii_iface_file_.writing_", f.writing_); pi("engine.rdii_iface_file_.binary_", f.binary_); pi("engine.rdii_iface_file_.rdii_step_", f.rdii_step_); pi("engine.rdii_iface_file_.flow_units_", f.flow_units_);'
             ' { std::vector<long long> t(f.node_idx_.begin(), f.node_idx_.end()); vi("engine.rdii_iface_file_.node_idx_", t); } { std::vector<double> t(f.flows_.begin(), f.flows_.end()); vf("engine.rdii_iface_file_.flows_", t); }'
             ' pf("engine.rdii_iface_file_.start_date_", f.start_date_); pf("engine.rdii_iface_file_.end_date_", f.end_date_); pi("engine.rdii_iface_file_.have_record_", f.have_record_); }')
    R.append('    { let f = &mut e.rdii_iface_file_; pi(o, "engine.rdii_iface_file_.open", f.fp_.is_some() as i64); let pos = f.tell(); pi(o, "engine.rdii_iface_file_.pos", pos);'
             ' pi(o, "engine.rdii_iface_file_.writing_", f.writing_ as i64); pi(o, "engine.rdii_iface_file_.binary_", f.binary_ as i64); pi(o, "engine.rdii_iface_file_.rdii_step_", f.rdii_step_ as i64); pi(o, "engine.rdii_iface_file_.flow_units_", f.flow_units_ as i64);'
             ' vi(o, "engine.rdii_iface_file_.node_idx_", &f.node_idx_.iter().map(|x| *x as i64).collect::<Vec<i64>>()); vf(o, "engine.rdii_iface_file_.flows_", &f.flows_.iter().map(|x| *x as f64).collect::<Vec<f64>>());'
             ' pf(o, "engine.rdii_iface_file_.start_date_", f.start_date_); pf(o, "engine.rdii_iface_file_.end_date_", f.end_date_); pi(o, "engine.rdii_iface_file_.have_record_", f.have_record_ as i64); }')
    # the built-in output and report plugins' run-time members
    C.append('  for (auto* p : e.plugins_.output_plugins()) { auto* d = dynamic_cast<openswmm::DefaultOutputPlugin*>(p); if (!d) continue; pi("engine.out_plugin.out_file_", d->out_file_ != nullptr); pi("engine.out_plugin.flow_units_code_", d->flow_units_code_); pf("engine.out_plugin.ucf_length_", d->ucf_length_); ps("engine.out_plugin.last_error_", d->last_error_); pi("engine.out_plugin.step_count_", (long long)d->step_count_); pi("engine.out_plugin.n_periods_", d->n_periods_); }')
    R.append('    for d in e.plugin_factory().output_plugins() { pi(o, "engine.out_plugin.out_file_", d.out_file_.is_some() as i64); pi(o, "engine.out_plugin.flow_units_code_", d.flow_units_code_ as i64); pf(o, "engine.out_plugin.ucf_length_", d.ucf_length_); ps(o, "engine.out_plugin.last_error_", &d.last_error_); pi(o, "engine.out_plugin.step_count_", d.step_count() as i64); pi(o, "engine.out_plugin.n_periods_", d.n_periods() as i64); }')
    C.append('  for (auto* p : e.plugins_.report_plugins()) { auto* d = dynamic_cast<openswmm::DefaultReportPlugin*>(p); if (!d) continue; pi("engine.rpt_plugin.file_", d->file_ != nullptr); ps("engine.rpt_plugin.last_error_", d->last_error_); }')
    R.append('    for d in e.plugin_factory().report_plugins() { pi(o, "engine.rpt_plugin.file_", d.file_.is_some() as i64); ps(o, "engine.rpt_plugin.last_error_", &d.last_error_); }')
    C.append('}'); C.append('} // namespace ctxdump')
    R.append('    std::mem::take(o)'); R.append('}')
    open(f'{ROOT}/gate/diff/lifedump.hpp', 'w').write('\n'.join(C) + '\n')
    open(f'{ROOT}/crates/engine/examples/common/lifedump.rs', 'w').write('\n'.join(R) + '\n')
    print('lifedump: written')

def main():
    collect_enums()
    C = ['// generated by tools/ctxdump.py: prints a SimulationContext; the Rust twin prints the same',
         '#pragma once', '#include <cstdio>', '#include <cstring>', '#include <cstdint>', '#include <string>', '#include <vector>',
         '#include <algorithm>', '#include <variant>', '#include "core/SimulationContext.hpp"',
         '// private members are read through the explicit-instantiation access rule (a pointer to a',
         '// private member may be named in an explicit instantiation); the class layout is untouched',
         'namespace ctxdump_access {',
         'template <typename Tag, typename Tag::type M> struct Rob { friend typename Tag::type get(Tag) { return M; } };',
         'struct T_UserFlags_def_index_ { typedef std::unordered_map<std::string, std::size_t> openswmm::UserFlags::*type; friend type get(T_UserFlags_def_index_); };',
         'template struct Rob<T_UserFlags_def_index_, &openswmm::UserFlags::def_index_>;',
         '} // namespace ctxdump_access',
         'namespace ctxdump {',
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
    # tables: every Table member (the FILE* handle is not data) and the by-name index
    dump_struct('data/TableData.hpp', 'TableData', 'c.tables', 'c.tables', 'tables', C, R)
    # unknown headers, deferred rows, passthrough
    C.append('  for (auto& h : c.unknown_section_headers) { ps("unknown_header", h.first); pi("unknown_header.line", h.second); }')
    R.append('    for h in &c.unknown_section_headers { ps(o, "unknown_header", &h.0); pi(o, "unknown_header.line", h.1); }')
    C.append('  for (auto& kv : c.passthrough_sections) { ps("passthrough", kv.first); vs("passthrough.lines", kv.second); }')
    R.append('    for (k, v) in &c.passthrough_sections { ps(o, "passthrough", k); vs(o, "passthrough.lines", v); }')
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
    # UserFlags::def_index_ (private upstream: read through ctxdump_access)
    C.append('  { const auto& m = c.user_flags.*get(ctxdump_access::T_UserFlags_def_index_()); std::vector<std::pair<std::string, std::size_t>> t(m.begin(), m.end()); std::sort(t.begin(), t.end()); for (auto& kv : t) { ps("uflag_def_index", kv.first); pi("uflag_def_index.v", (long long)kv.second); } }')
    R.append('    { let mut t: Vec<_> = c.user_flags.def_index().iter().collect(); t.sort(); for (k, v) in t { ps(o, "uflag_def_index", k); pi(o, "uflag_def_index.v", *v as i64); } }')
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
    for n in CTX_MEMBERS:
        t, a = ctx_fields[n]
        ctx_member(n, t, a, C, R)
    C.append('}'); C.append('} // namespace ctxdump')
    R.append('    std::mem::take(o)'); R.append('}')
    EC, ER = engine_dumper()
    C += EC
    R += ER
    QC, QR = quality_dumper()
    C += QC
    R += QR
    open(f'{ROOT}/gate/diff/ctxdump.hpp', 'w').write('\n'.join(C) + '\n')
    os.makedirs(f'{ROOT}/crates/engine/examples/common', exist_ok=True)
    open(f'{ROOT}/crates/engine/examples/common/ctxdump.rs', 'w').write('\n'.join(R) + '\n')
    print(f'{len(SKIPS)} skipped:', '; '.join(SKIPS))
    router_dumper()
    hydro_dumper()
    life_dumper()

if __name__ == '__main__':
    main()
