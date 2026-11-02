// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/PathResolver.cpp and PathResolver.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Lexical path resolution for external-file tokens (namespace openswmm::io). Paths are
// bytes. std::filesystem::path is core::fs_path::FsPath (libstdc++'s POSIX behaviour);
// utf8_path(s) is FsPath::new(s) and path_utf8(p, true) is p.generic_string(), as they
// are on POSIX. std::isalpha/toupper/tolower are the C locale's (ASCII).

use crate::core::fs_path::FsPath;
use crate::core::options::Bytes;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum PathClass {
    Relative = 0,
    AbsoluteSameVolume = 1,
    AbsoluteCrossVolume = 2,
    Invalid = 3,
}

#[derive(Debug, Clone)]
pub struct RelativeResult {
    pub path: Bytes,
    pub classification: PathClass,
    pub up_levels: i32,
    pub warning: Bytes,
}

impl Default for RelativeResult {
    fn default() -> Self {
        RelativeResult { path: Bytes::new(), classification: PathClass::Invalid, up_levels: 0, warning: Bytes::new() }
    }
}

/// `makeRelative`'s default `max_up_levels`.
pub const MAKE_RELATIVE_DEFAULT_MAX_UP_LEVELS: i32 = 16;

// ---- internal helpers (anonymous namespace) ---------------------------------

fn to_forward_slashes(s: &[u8]) -> Bytes {
    s.iter().map(|&c| if c == b'\\' { b'/' } else { c }).collect()
}

fn trim_trailing_slash(mut s: Bytes) -> Bytes {
    if s.len() <= 1 {
        return s;
    }
    if s.len() == 3 && s[1] == b':' && (s[2] == b'/' || s[2] == b'\\') {
        return s;
    }
    if s.len() >= 2 && (s[0] == b'/' && s[1] == b'/') {
        let slashes = s.iter().filter(|&&c| c == b'/').count();
        if slashes < 4 {
            return s;
        }
    }
    while s.len() > 1 && (*s.last().unwrap() == b'/' || *s.last().unwrap() == b'\\') {
        s.pop();
    }
    s
}

fn find(hay: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.len() > hay.len() {
        return None;
    }
    (0..=hay.len() - needle.len()).find(|&i| &hay[i..i + needle.len()] == needle)
}

fn collapse_slashes_and_dots(s: &[u8]) -> Bytes {
    let mut out = Bytes::with_capacity(s.len());
    let unc = s.len() >= 2 && s[0] == b'/' && s[1] == b'/';
    for (i, &c) in s.iter().enumerate() {
        if c == b'/' {
            if unc && i == 1 && out.len() == 1 && *out.last().unwrap() == b'/' {
                out.push(c);
                continue;
            }
            if !out.is_empty() && *out.last().unwrap() == b'/' {
                continue;
            }
            out.push(b'/');
        } else {
            out.push(c);
        }
    }
    loop {
        match find(&out, b"/./") {
            None => break,
            Some(p) => {
                out.drain(p..p + 2);
            }
        }
    }
    if out.len() >= 2 && out[0] == b'.' && out[1] == b'/' {
        out.drain(0..2);
    }
    if out.len() >= 2 && out[out.len() - 2] == b'/' && *out.last().unwrap() == b'.' {
        out.pop();
        out.pop();
    }
    out
}

fn root_name(path: &[u8]) -> Bytes {
    if path.len() >= 2 && path[0].is_ascii_alphabetic() && path[1] == b':' {
        let mut r = path[..2].to_vec();
        r[0] = r[0].to_ascii_uppercase();
        return r;
    }
    if path.len() >= 2 && (path[0] == b'/' || path[0] == b'\\') && (path[1] == b'/' || path[1] == b'\\') {
        let mut i = 2;
        while i < path.len() && path[i] != b'/' && path[i] != b'\\' {
            i += 1;
        }
        if i < path.len() {
            i += 1;
            while i < path.len() && path[i] != b'/' && path[i] != b'\\' {
                i += 1;
            }
        }
        return path[..i].to_vec();
    }
    Bytes::new()
}

fn count_up_levels(normalised: &[u8]) -> i32 {
    let mut n = 0;
    let mut i = 0usize;
    while i + 2 <= normalised.len() && normalised[i] == b'.' && normalised[i + 1] == b'.' {
        if i + 2 == normalised.len() {
            n += 1;
            break;
        }
        if normalised[i + 2] != b'/' {
            break;
        }
        n += 1;
        i += 3;
    }
    n
}

// ---- public API ---------------------------------------------------------------

pub fn is_absolute_path(path: &[u8]) -> bool {
    if path.is_empty() {
        return false;
    }
    if path[0] == b'/' {
        return true;
    }
    if path.len() >= 2 && path[0].is_ascii_alphabetic() && path[1] == b':' {
        return true;
    }
    if path.len() >= 2 && (path[0] == b'\\' || path[0] == b'/') && (path[1] == b'\\' || path[1] == b'/') {
        return true;
    }
    false
}

pub fn normalise_separators(path: &[u8]) -> Bytes {
    if path.is_empty() {
        return path.to_vec();
    }
    let out = to_forward_slashes(path);
    let out = collapse_slashes_and_dots(&out);
    trim_trailing_slash(out)
}

pub fn parent_dir(file_path: &[u8]) -> Bytes {
    if file_path.is_empty() {
        return Bytes::new();
    }
    let last = *file_path.last().unwrap();
    let norm = normalise_separators(file_path);
    if last == b'/' || last == b'\\' {
        return norm;
    }
    let pos = match norm.iter().rposition(|&c| c == b'/') {
        None => return Bytes::new(),
        Some(p) => p,
    };
    if pos < 2 {
        return norm[..pos + 1].to_vec();
    }
    if pos == 2 && norm.len() >= 3 && norm[1] == b':' {
        return norm[..3].to_vec();
    }
    norm[..pos].to_vec()
}

pub fn resolve_relative(stored_token: &[u8], anchor_dir: &[u8]) -> Bytes {
    if stored_token.is_empty() {
        return Bytes::new();
    }
    if is_absolute_path(stored_token) {
        let p = FsPath::new(&to_forward_slashes(stored_token));
        return p.generic_string();
    }
    if anchor_dir.is_empty() {
        let p = FsPath::new(&to_forward_slashes(stored_token));
        return p.lexically_normal().generic_string();
    }
    let joined = FsPath::new(&to_forward_slashes(anchor_dir)).join(&FsPath::new(&to_forward_slashes(stored_token)));
    joined.lexically_normal().generic_string()
}

pub fn make_relative(target_absolute: &[u8], anchor_dir: &[u8], max_up_levels: i32) -> RelativeResult {
    let mut r = RelativeResult::default();

    if target_absolute.is_empty() {
        r.classification = PathClass::Invalid;
        r.warning = b"empty target path".to_vec();
        return r;
    }

    if anchor_dir.is_empty() {
        r.path = normalise_separators(target_absolute);
        r.classification = if is_absolute_path(target_absolute) {
            PathClass::AbsoluteSameVolume
        } else {
            PathClass::Invalid
        };
        if r.classification != PathClass::Relative {
            r.warning = b"no anchor directory provided".to_vec();
        }
        return r;
    }

    let tgt_fs = to_forward_slashes(target_absolute);
    let anc_fs = to_forward_slashes(anchor_dir);

    let tgt_root = root_name(&tgt_fs);
    let anc_root = root_name(&anc_fs);
    let eq_ci = |a: &[u8], b: &[u8]| -> bool {
        if a.len() != b.len() {
            return false;
        }
        a.iter().zip(b.iter()).all(|(x, y)| x.to_ascii_lowercase() == y.to_ascii_lowercase())
    };
    if !eq_ci(&tgt_root, &anc_root) {
        r.path = FsPath::new(&tgt_fs).generic_string();
        r.classification = PathClass::AbsoluteCrossVolume;
        let mut w = b"target on a different volume ('".to_vec();
        w.extend_from_slice(&tgt_root);
        w.extend_from_slice(b"' vs '");
        w.extend_from_slice(&anc_root);
        w.extend_from_slice("') \u{2014} cannot express relatively".as_bytes());
        r.warning = w;
        return r;
    }

    let mut tgt_body = tgt_fs.clone();
    let mut anc_body = anc_fs.clone();
    if !tgt_root.is_empty() {
        tgt_body.drain(0..tgt_root.len());
    }
    if !anc_root.is_empty() {
        anc_body.drain(0..anc_root.len());
    }
    if tgt_body.is_empty() || tgt_body[0] != b'/' {
        tgt_body.insert(0, b'/');
    }
    if anc_body.is_empty() || anc_body[0] != b'/' {
        anc_body.insert(0, b'/');
    }

    let tgt = FsPath::new(&tgt_body);
    let anc = FsPath::new(&anc_body);
    let prox = tgt.lexically_proximate(&anc);

    if prox.empty() || prox.is_absolute() {
        r.path = FsPath::new(&tgt_fs).generic_string();
        r.classification = PathClass::AbsoluteSameVolume;
        r.warning = b"lexical relative form unavailable".to_vec();
        return r;
    }

    let out = prox.generic_string();
    let mut out = collapse_slashes_and_dots(&out);
    if out.is_empty() {
        out = b".".to_vec();
    }

    let up = count_up_levels(&out);
    if up > max_up_levels {
        r.path = FsPath::new(&tgt_fs).generic_string();
        r.classification = PathClass::AbsoluteSameVolume;
        r.warning = format!("relative form exceeds {} '..' levels (would be {})", max_up_levels, up).into_bytes();
        return r;
    }

    r.path = out;
    r.classification = PathClass::Relative;
    r.up_levels = up;
    r
}
