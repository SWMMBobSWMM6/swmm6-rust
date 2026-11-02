#!/usr/bin/env python3
"""Mechanical first pass for upstream's structure-of-arrays headers (data/*.hpp).

Translates the parts that follow fixed shapes and marks everything else for a person:
  * `enum class X : T { A = 0, ... }`          -> #[repr] Rust enum, CamelCase variants
  * `std::vector<T> name;` / `T name = v;`      -> struct fields (+ Default)
  * method bodies made only of per-field statements:
      F.assign(n, D)  F.resize(n, D)  F.clear()  F.reserve(n)  F.shrink_to_fit()
      g(F, D) [grow lambda]  r(F) [reserve lambda]  e(F) [erase lambda]
      std::fill(F.begin(), F.end(), V)  std::copy(A.begin(), A.end(), B.begin())
Any other line is emitted as a `// MANUAL:` comment with the C++ text and must be
translated by hand; the output does not compile until every MANUAL line is resolved.

This is a tool for the translation, not part of the engine. Usage: soa2rs.py X.hpp > x.rs
"""
import re, sys

TYPES = {'double': 'f64', 'float': 'f32', 'int': 'i32', 'long': 'i64', 'bool': 'bool',
         'uint8_t': 'u8', 'int8_t': 'i8', 'char': 'i8', 'int64_t': 'i64', 'uint64_t': 'u64',
         'int32_t': 'i32', 'uint32_t': 'u32', 'std::size_t': 'usize', 'size_t': 'usize',
         'std::string': 'Bytes', 'std::int8_t': 'i8', 'std::uint8_t': 'u8', 'std::int32_t': 'i32',
         'std::int64_t': 'i64', 'std::uint32_t': 'u32'}

def camel(s):
    return ''.join(p.capitalize() for p in s.lower().split('_')) or s

def rtype(t, enums):
    t = t.strip()
    m = re.fullmatch(r'std::vector<(.+)>', t)
    if m:
        return f'Vec<{rtype(m.group(1), enums)}>'
    if t in TYPES:
        return TYPES[t]
    if t in enums:
        return t
    return f'/*MANUAL type*/ {t}'

def rval(v, ty, enums):
    v = v.strip()
    m = re.fullmatch(r'static_cast<[\w:]+>\((.+)\)', v)
    if m:
        v = m.group(1).strip()
    if v in ('std::string{}', 'std::string()', '""'):
        return 'Vec::new()'
    m = re.fullmatch(r'(\w+)::(\w+)', v)
    if m and m.group(1) in enums:
        return f'{m.group(1)}::{camel(m.group(2))}'
    if v in ('true', 'false'):
        return v
    m = re.fullmatch(r'(-?\d+)[lLuU]*', v)
    if m:
        n = m.group(1)
        return n + '.0' if ty in ('f64', 'f32') else n
    m = re.fullmatch(r'-?(\d+\.\d*|\.\d+|\d+)([eE][-+]?\d+)?', v)
    if m:
        return v if ty in ('f64', 'f32') else f'/*MANUAL value*/ {v}'
    return f'/*MANUAL value*/ {v}'

def strip_comments(src):
    src = re.sub(r'/\*.*?\*/', '', src, flags=re.S)
    src = re.sub(r'//[^\n]*', '', src)
    # clang-format splits long declarations: "std::vector<double>\n    name;"
    return re.sub(r'(std::vector<[\w:<> ]+>)\s*\n\s*(\w+\s*;)', r'\1 \2', src)

def main(path):
    src = strip_comments(open(path).read())
    out = []
    enums = {}
    for m in re.finditer(r'enum\s+class\s+(\w+)\s*:\s*([\w:]+)\s*\{([^}]*)\}', src):
        name, base, body = m.groups()
        enums[name] = TYPES.get(base, 'i32')
    for m in re.finditer(r'enum\s+class\s+(\w+)\s*:\s*([\w:]+)\s*\{([^}]*)\}', src):
        name, base, body = m.groups()
        out.append('#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]')
        out.append(f'#[repr({enums[name]})]')
        out.append(f'pub enum {name} {{')
        nxt = 0
        for item in [x.strip() for x in body.split(',') if x.strip()]:
            if '=' in item:
                k, v = [s.strip() for s in item.split('=')]
                nxt = int(v)
            else:
                k = item
            out.append(f'    {camel(k)} = {nxt},')
            nxt += 1
        out.append('}\n')
    for sm in re.finditer(r'struct\s+(\w+)\s*\{', src):
        name = sm.group(1)
        # find matching brace
        i = sm.end(); depth = 1
        while depth:
            depth += {'{': 1, '}': -1}.get(src[i], 0); i += 1
        body = src[sm.end():i - 1]
        fields = []; methods = []
        j = 0; lines = body.split('\n'); k = 0
        while k < len(lines):
            line = lines[k].strip()
            fm = re.fullmatch(r'(std::vector<[\w:<> ]+>|[\w:]+)\s+(\w+)\s*(=\s*(.+))?;', line)
            mm = re.match(r'(?:\w[\w:<>&*\s]*?)\s+(\w+)\s*\(([^)]*)\)\s*(const)?\s*(noexcept)?\s*\{', line)
            if fm and not mm:
                ty, nm, _, init = fm.groups()
                fields.append((nm, rtype(ty, enums), init))
                k += 1; continue
            if mm and line.count('{') == line.count('}'):
                methods.append([line, '// one-line: ' + line, '}']); k += 1; continue
            if mm:
                d = line.count('{') - line.count('}'); blines = [line]
                while d > 0:
                    k += 1; l2 = lines[k].strip(); blines.append(l2); d += l2.count('{') - l2.count('}')
                methods.append(blines); k += 1; continue
            k += 1
        ftype = {n: t for n, t, _ in fields}
        out.append('#[derive(Debug, Clone, Default)]')
        out.append(f'pub struct {name} {{')
        for n, t, init in fields:
            out.append(f'    pub {n}: {t},')
        out.append('}\n')
        inits = [(n, t, i) for n, t, i in fields if i is not None]
        if inits:
            out.append(f'// MANUAL: member initializers to carry into a Default impl:')
            for n, t, i in inits:
                out.append(f'//   {n} = {rval(i, t, enums)}')
            out.append('')
        out.append(f'impl {name} {{')
        for blines in methods:
            out.extend(translate_method(blines, ftype, enums))
        out.append('}\n')
    print('\n'.join(out))

def elem(t):
    m = re.fullmatch(r'Vec<(.+)>', t)
    return m.group(1) if m else t

def translate_method(bl, ftype, enums):
    head = bl[0]
    mm = re.match(r'(?:\w[\w:<>&*\s]*?)\s+(\w+)\s*\(([^)]*)\)', head)
    name, args = mm.group(1), mm.group(2)
    rargs = []
    for a in [x.strip() for x in args.split(',') if x.strip()]:
        parts = a.split()
        rargs.append(f'{parts[-1]}: {TYPES.get(" ".join(parts[:-1]).replace("const ", ""), "/*MANUAL*/")}')
    o = [f'    pub fn {name}(&mut self{"".join(", " + r for r in rargs)}) {{']
    for line in bl[1:-1]:
        s = line.strip()
        if not s:
            continue
        def F(f): return f'self.{f}'
        m = re.fullmatch(r'(\w+)\.(assign|resize)\((\w+),\s*(.+)\);', s)
        if m and m.group(1) in ftype:
            f, op, n, d = m.groups(); et = elem(ftype[f])
            if op == 'assign':
                o.append(f'        {F(f)}.clear();')
            o.append(f'        {F(f)}.resize({n}, {rval(d, et, enums)});'); continue
        m = re.fullmatch(r'(\w+)\.(clear|shrink_to_fit)\(\);', s)
        if m and m.group(1) in ftype:
            o.append(f'        {F(m.group(1))}.{m.group(2)}();'); continue
        m = re.fullmatch(r'(\w+)\.reserve\((\w+)\);', s)
        if m and m.group(1) in ftype:
            o.append(f'        {F(m.group(1))}.reserve({m.group(2)}.saturating_sub({F(m.group(1))}.len()));'); continue
        m = re.fullmatch(r'g\((\w+),\s*(.+)\);', s)
        if m and m.group(1) in ftype:
            o.append(f'        {F(m.group(1))}.resize(un, {rval(m.group(2), elem(ftype[m.group(1)]), enums)});'); continue
        m = re.fullmatch(r'r\((\w+)\);', s)
        if m and m.group(1) in ftype:
            o.append(f'        {F(m.group(1))}.reserve(un.saturating_sub({F(m.group(1))}.len()));'); continue
        m = re.fullmatch(r'e\((\w+)\);', s)
        if m and m.group(1) in ftype:
            o.append(f'        if ui < {F(m.group(1))}.len() {{ {F(m.group(1))}.remove(ui); }}'); continue
        m = re.fullmatch(r'std::fill\((\w+)\.begin\(\),\s*\1\.end\(\),\s*(.+)\);', s)
        if m and m.group(1) in ftype:
            o.append(f'        {F(m.group(1))}.fill({rval(m.group(2), elem(ftype[m.group(1)]), enums)});'); continue
        m = re.fullmatch(r'std::copy\((\w+)\.begin\(\),\s*\1\.end\(\),\s*(\w+)\.begin\(\)\);', s)
        if m and m.group(1) in ftype:
            a, b = m.groups()
            o.append(f'        {{ let n = {F(a)}.len(); {F(b)}[..n].copy_from_slice(&{F(a)}); }}'); continue
        m = re.fullmatch(r'const auto un = static_cast<std::size_t>\(n\);', s)
        if m:
            o.append('        let un = usize::try_from(n).unwrap_or(0);'); continue
        m = re.fullmatch(r'const auto ui = static_cast<std::size_t>\(idx\);', s)
        if m:
            o.append('        let ui = idx as usize;'); continue
        if re.fullmatch(r'auto [gre] = \[&\].*', s) or s in ('};',):
            continue
        o.append(f'        // MANUAL: {s}')
    o.append('    }\n')
    return o

if __name__ == '__main__':
    main(sys.argv[1])
