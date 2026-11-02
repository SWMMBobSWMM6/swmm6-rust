// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/NodeCoupling2D.cpp and NodeCoupling2D.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// (Upstream's files have no SPDX/license comment block; their doc comments carry the
// copyright and "@license Apache-2.0", reproduced above.)
//
// The [2D_NODE_COUPLING] section reader (nodecoupling2d::readSection) and its helpers:
// the token tables, findToken, choices, number, isOutfallKey, typeToken, parseType,
// keyToken, parseKey, parseValue, keyAllowed and checkType, with defaultType and
// effectiveType beside them. The data types are in crate::data::node_coupling_2d.
//
// Differences in form only: `std::strtod` is core::cstdlib::strtod on the token up to its
// first NUL, as c_str() would give it; `end != c_str` is "consumed > 0" and `*end == '\0'`
// is "consumed it all". `const char*` results are byte-string slices; an out-parameter
// enum is returned as Option.
//
// Not translated here: shortest, formatValue, formatParams and formatRow (the .inp writer's
// row formatting, which needs snprintf's %g).

use crate::core::context::SimulationContext;
use crate::core::cstdlib::{c_str, strtod};
use crate::core::options::Bytes;
use crate::data::node::NodeType;
use crate::data::node_coupling_2d::{unset_node_coupling_2d_params, NodeCoupling2DKey, NodeCoupling2DType, K_NODE_COUPLING_2D_KEYS};
use crate::input::tokenizer;

fn cat(parts: &[&[u8]]) -> Bytes {
    parts.concat()
}

const K_TYPE_TOKENS: [&[u8]; 6] = [b"NONE", b"MANHOLE", b"GRATE", b"CURB", b"OPEN", b"OUTFALL"];

const K_KEY_TOKENS: [&[u8]; K_NODE_COUPLING_2D_KEYS as usize] = [
    b"D",
    b"PERIMETER",
    b"AREA",
    b"CHAMBER_AREA",
    b"LID_OPEN_FRACTION",
    b"CREST",
    b"COEFFS",
    b"CW",
    b"CSW",
    b"COI",
    b"COO",
    b"LAW",
    b"SPREAD",
    b"TAILWATER_RELAX",
];

const K_CREST_TOKENS: [&[u8]; 4] = [b"NODE_RIM", b"MESH_BED", b"MAX", b"MIN"];
const K_COEFFS_TOKENS: [&[u8]; 4] = [b"DEFAULT", b"RUBINATO2017", b"RUBINATO2017_LOWER", b"RUBINATO2017_UPPER"];
const K_LAW_TOKENS: [&[u8]; 2] = [b"ORIFICE", b"WEIR_ORIFICE"];

const TYPES: [NodeCoupling2DType; 6] = [
    NodeCoupling2DType::None,
    NodeCoupling2DType::Manhole,
    NodeCoupling2DType::Grate,
    NodeCoupling2DType::Curb,
    NodeCoupling2DType::Open,
    NodeCoupling2DType::Outfall,
];

const KEYS: [NodeCoupling2DKey; K_NODE_COUPLING_2D_KEYS as usize] = [
    NodeCoupling2DKey::D,
    NodeCoupling2DKey::Perimeter,
    NodeCoupling2DKey::Area,
    NodeCoupling2DKey::ChamberArea,
    NodeCoupling2DKey::LidOpenFraction,
    NodeCoupling2DKey::Crest,
    NodeCoupling2DKey::Coeffs,
    NodeCoupling2DKey::Cw,
    NodeCoupling2DKey::Csw,
    NodeCoupling2DKey::Coi,
    NodeCoupling2DKey::Coo,
    NodeCoupling2DKey::Law,
    NodeCoupling2DKey::Spread,
    NodeCoupling2DKey::TailwaterRelax,
];

/// findToken(): the first index whose literal equals the upper-cased token, or -1.
fn find_token(tokens: &[&[u8]], tok: &[u8]) -> i32 {
    let up = tokenizer::to_upper(tok);
    for (i, t) in tokens.iter().enumerate() {
        if up.as_slice() == *t {
            return i as i32;
        }
    }
    -1
}

/// choices(): the literals joined by '|'.
fn choices(tokens: &[&[u8]]) -> Bytes {
    let mut s = Bytes::new();
    for (i, t) in tokens.iter().enumerate() {
        if i != 0 {
            s.push(b'|');
        }
        s.extend_from_slice(t);
    }
    s
}

/// number(): one complete finite token through strtod. `out` always receives strtod's result.
fn number(tok: &[u8], out: &mut f64) -> bool {
    let cs = c_str(tok);
    let (x, used) = strtod(cs);
    *out = x;
    used != 0 && used == cs.len() && x.is_finite()
}

fn is_outfall_key(k: NodeCoupling2DKey) -> bool {
    k == NodeCoupling2DKey::Spread || k == NodeCoupling2DKey::TailwaterRelax
}

/// typeToken().
pub fn type_token(t: NodeCoupling2DType) -> &'static [u8] {
    let i = t as i32;
    if (0..6).contains(&i) {
        K_TYPE_TOKENS[i as usize]
    } else {
        b"NONE"
    }
}

/// parseType().
pub fn parse_type(token: &[u8]) -> Option<NodeCoupling2DType> {
    let i = find_token(&K_TYPE_TOKENS, token);
    if i < 0 {
        return None;
    }
    Some(TYPES[i as usize])
}

/// keyToken().
pub fn key_token(k: NodeCoupling2DKey) -> &'static [u8] {
    let i = k as i32;
    if i >= 0 && i < K_NODE_COUPLING_2D_KEYS {
        K_KEY_TOKENS[i as usize]
    } else {
        b""
    }
}

/// parseKey().
pub fn parse_key(token: &[u8]) -> Option<NodeCoupling2DKey> {
    let i = find_token(&K_KEY_TOKENS, token);
    if i < 0 {
        return None;
    }
    Some(KEYS[i as usize])
}

/// parseValue(): "" when `token` is a legal value for key `k` (stored in `out`), otherwise
/// the diagnostic. `out` receives strtod's result even when the value is rejected.
pub fn parse_value(k: NodeCoupling2DKey, token: &[u8], out: &mut f64) -> Bytes {
    let key = key_token(k);
    let i;
    match k {
        NodeCoupling2DKey::Crest => {
            i = find_token(&K_CREST_TOKENS, token);
            if i < 0 {
                return cat(&[key, b" must be ", &choices(&K_CREST_TOKENS), b": ", token]);
            }
            *out = i as f64;
            return Bytes::new();
        }
        NodeCoupling2DKey::Coeffs => {
            i = find_token(&K_COEFFS_TOKENS, token);
            if i < 0 {
                return cat(&[key, b" must be ", &choices(&K_COEFFS_TOKENS), b": ", token]);
            }
            *out = i as f64;
            return Bytes::new();
        }
        NodeCoupling2DKey::Law => {
            i = find_token(&K_LAW_TOKENS, token);
            if i < 0 {
                return cat(&[key, b" must be ", &choices(&K_LAW_TOKENS), b": ", token]);
            }
            *out = i as f64;
            return Bytes::new();
        }
        _ => {}
    }
    if !number(token, out) {
        return cat(&[key, b" must be a number: ", token]);
    }
    match k {
        NodeCoupling2DKey::LidOpenFraction => {
            if !(*out > 0.0 && *out <= 1.0) {
                return cat(&[key, b" must be in (0, 1]: ", token]);
            }
        }
        NodeCoupling2DKey::Spread | NodeCoupling2DKey::TailwaterRelax => {
            if *out < 0.0 {
                return cat(&[key, b" must be >= 0: ", token]);
            }
        }
        _ => {
            if !(*out > 0.0) {
                return cat(&[key, b" must be > 0: ", token]);
            }
        }
    }
    Bytes::new()
}

/// keyAllowed().
pub fn key_allowed(t: NodeCoupling2DType, k: NodeCoupling2DKey) -> bool {
    match t {
        NodeCoupling2DType::None => false,
        NodeCoupling2DType::Outfall => is_outfall_key(k),
        _ => !is_outfall_key(k),
    }
}

/// defaultType(): the class default of a node without a row.
pub fn default_type(ctx: &SimulationContext, node: i32) -> NodeCoupling2DType {
    let un = node as usize;
    match ctx.nodes.type_[un] {
        NodeType::Outfall => return NodeCoupling2DType::Outfall,
        NodeType::Storage => return NodeCoupling2DType::Open,
        NodeType::Divider => return NodeCoupling2DType::Manhole,
        NodeType::Junction => {}
    }
    // Inlet junctions already route street flow into the sewer through
    // [INLETS]; coupling them to the mesh too would route it twice (D6).
    // Virtual junctions have no rim.
    let inlet = un < ctx.nodes.is_inlet.len() && ctx.nodes.is_inlet[un] != 0;
    let virtual_ = un < ctx.nodes.is_virtual.len() && ctx.nodes.is_virtual[un] != 0;
    if inlet || virtual_ {
        NodeCoupling2DType::None
    } else {
        NodeCoupling2DType::Manhole
    }
}

/// effectiveType(): the authored row's type, or the class default.
pub fn effective_type(ctx: &SimulationContext, node: i32) -> NodeCoupling2DType {
    let r = ctx.node_coupling_2d.row_of(node);
    if r >= 0 {
        ctx.node_coupling_2d.type_[r as usize]
    } else {
        default_type(ctx, node)
    }
}

/// checkType(): "" when type `t` is legal for `node`, otherwise the reason.
pub fn check_type(ctx: &SimulationContext, node: i32, t: NodeCoupling2DType) -> Bytes {
    let outfall = ctx.nodes.type_[node as usize] == NodeType::Outfall;
    if outfall && t != NodeCoupling2DType::Outfall && t != NodeCoupling2DType::None {
        return cat(&[b"an outfall couples as OUTFALL or NONE, not ", type_token(t)]);
    }
    if !outfall && t == NodeCoupling2DType::Outfall {
        return b"type OUTFALL needs an outfall node".to_vec();
    }
    if t == NodeCoupling2DType::Open && ctx.nodes.type_[node as usize] != NodeType::Storage {
        return b"type OPEN needs a storage node".to_vec();
    }
    Bytes::new()
}

/// nodecoupling2d::readSection: [2D_NODE_COUPLING] `NODE TYPE [KEY=VALUE ...]`.
pub fn read_section(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let t = tokenizer::tokenize(line);
        if t.is_empty() {
            continue;
        }
        let fail = |ctx: &mut SimulationContext, why: &[u8]| {
            ctx.errors.push(cat(&[b"[2D_NODE_COUPLING] ", why, " — line: ".as_bytes(), line]));
        };
        if t.len() < 2 {
            fail(ctx, b"expected NODE TYPE [KEY=VALUE ...]");
            continue;
        }
        let n = ctx.node_names.find(&t[0]);
        if n < 0 {
            ctx.deferred_section_rows.push((b"2D_NODE_COUPLING".to_vec(), line.clone()));
            continue;
        }
        let type_ = match parse_type(&t[1]) {
            Some(ty) => ty,
            None => {
                let why = cat(&[b"unknown type ", &t[1], b" (expected ", &choices(&K_TYPE_TOKENS), b")"]);
                fail(ctx, &why);
                continue;
            }
        };
        let why = check_type(ctx, n, type_);
        if !why.is_empty() {
            let w = cat(&[&t[0], b": ", &why]);
            fail(ctx, &w);
            continue;
        }
        if ctx.node_coupling_2d.row_of(n) >= 0 {
            let w = cat(&[b"duplicate row for node ", &t[0]]);
            fail(ctx, &w);
            continue;
        }

        let mut params = unset_node_coupling_2d_params();
        let mut err = Bytes::new();
        let mut i = 2usize;
        while i < t.len() && err.is_empty() {
            let ti = &t[i];
            let eq = ti.iter().position(|&c| c == b'=');
            let eq = match eq {
                Some(e) if e != 0 && e + 1 != ti.len() => e,
                _ => {
                    err = cat(&[b"expected KEY=VALUE: ", ti]);
                    break;
                }
            };
            let key = match parse_key(&ti[..eq]) {
                Some(k) => k,
                None => {
                    err = cat(&[b"unknown key ", &ti[..eq]]);
                    break;
                }
            };
            if !key_allowed(type_, key) {
                err = cat(&[key_token(key), b" does not apply to type ", type_token(type_)]);
                break;
            }
            let slot = &mut params[key as usize];
            if !slot.is_nan() {
                err = cat(&[b"duplicate key ", key_token(key)]);
                break;
            }
            err = parse_value(key, &ti[eq + 1..], slot);
            i += 1;
        }
        if !err.is_empty() {
            fail(ctx, &err);
            continue;
        }

        let rows = &mut ctx.node_coupling_2d;
        rows.node.push(n);
        rows.type_.push(type_);
        rows.params.push(params);
    }
}
