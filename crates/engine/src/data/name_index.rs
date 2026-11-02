// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/NameIndex.hpp and
// src/engine/core/StringCase.hpp (HydroCouple/openswmm.engine, branch swmm6_rel,
// commit 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the Apache License,
// Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Upstream keys an unordered_map with a case-folding hash and equality (ASCII a-z only).
// Here the key is the name folded the same way; lookups fold the query. Nothing iterates
// the map, so its order cannot reach a result. Index order is insertion order, held in
// `names`, which keeps the spelling as first written.

use std::collections::HashMap;

/// `ci_upper`: ASCII a-z to A-Z, every other byte unchanged.
#[inline]
pub fn ci_upper(c: u8) -> u8 {
    if c.is_ascii_lowercase() {
        c & !32
    } else {
        c
    }
}

/// Case-insensitive equality (`ieq`).
pub fn ieq(a: &[u8], b: &[u8]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(&x, &y)| ci_upper(x) == ci_upper(y))
}

fn fold(s: &[u8]) -> Vec<u8> {
    s.iter().map(|&c| ci_upper(c)).collect()
}

#[derive(Debug, Clone, Default)]
pub struct NameIndex {
    map: HashMap<Vec<u8>, i32>,
    names: Vec<Vec<u8>>,
    bulk_depth: i32,
}

impl NameIndex {
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a name; a duplicate is an error (upstream throws std::invalid_argument).
    pub fn add(&mut self, name: &[u8]) -> Result<i32, String> {
        let idx = self.try_add(name);
        if idx < 0 {
            return Err(format!("NameIndex: duplicate name '{}'", String::from_utf8_lossy(name)));
        }
        Ok(idx)
    }

    /// Add a name and return its index, or -1 if a name equal to it ignoring case exists.
    pub fn try_add(&mut self, name: &[u8]) -> i32 {
        let key = fold(name);
        if self.map.contains_key(&key) {
            return -1;
        }
        let idx = self.names.len() as i32;
        self.map.insert(key, idx);
        self.names.push(name.to_vec());
        idx
    }

    pub fn find(&self, name: &[u8]) -> i32 {
        *self.map.get(&fold(name)).unwrap_or(&-1)
    }

    pub fn try_find(&self, name: &[u8]) -> Option<i32> {
        self.map.get(&fold(name)).copied()
    }

    /// The spelling the name was first added with.
    pub fn canonical(&self, name: &[u8]) -> Option<&[u8]> {
        self.try_find(name).map(|i| self.names[i as usize].as_slice())
    }

    /// Upstream uses vector::at, which throws on a bad index; here it panics.
    pub fn name_of(&self, idx: i32) -> &[u8] {
        &self.names[usize::try_from(idx).expect("NameIndex::name_of: negative index")]
    }

    pub fn size(&self) -> i32 {
        self.names.len() as i32
    }

    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    pub fn reserve(&mut self, n: usize) {
        self.map.reserve(n);
        self.names.reserve(n);
    }

    pub fn clear(&mut self) {
        self.map.clear();
        self.names.clear();
    }

    pub fn pop_back(&mut self) {
        if let Some(tail) = self.names.pop() {
            self.map.remove(&fold(&tail));
        }
    }

    pub fn rename(&mut self, idx: i32, new_name: &[u8]) -> bool {
        if idx < 0 || idx >= self.names.len() as i32 {
            return false;
        }
        let key = fold(new_name);
        if let Some(&other) = self.map.get(&key) {
            if other != idx {
                return false;
            }
        }
        let old = fold(&self.names[idx as usize]);
        self.map.remove(&old);
        self.names[idx as usize] = new_name.to_vec();
        self.map.insert(key, idx);
        true
    }

    pub fn remove_at(&mut self, idx: i32) {
        if idx < 0 || idx >= self.names.len() as i32 {
            return;
        }
        self.names.remove(idx as usize);
        if self.bulk_depth > 0 {
            return;
        }
        self.rebuild_map();
    }

    pub fn begin_bulk_remove(&mut self) {
        self.bulk_depth += 1;
    }

    pub fn end_bulk_remove(&mut self) {
        if self.bulk_depth == 0 {
            return;
        }
        self.bulk_depth -= 1;
        if self.bulk_depth > 0 {
            return;
        }
        self.rebuild_map();
    }

    pub fn names(&self) -> &[Vec<u8>] {
        &self.names
    }

    fn rebuild_map(&mut self) {
        // map_[name] = i overwrites, so with duplicates (which add() prevents) the last wins
        self.map.clear();
        self.map.reserve(self.names.len());
        for (i, n) in self.names.iter().enumerate() {
            self.map.insert(fold(n), i as i32);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn case_insensitive_first_spelling_kept() {
        let mut ix = NameIndex::new();
        assert_eq!(ix.try_add(b"Node1"), 0);
        assert_eq!(ix.try_add(b"NODE1"), -1);
        assert_eq!(ix.find(b"node1"), 0);
        assert_eq!(ix.canonical(b"nOdE1"), Some(&b"Node1"[..]));
        // only ASCII folds: Latin-1 0xE9 and 0xC9 stay distinct, as in upstream
        assert_eq!(ix.try_add(b"\xe9"), 1);
        assert_eq!(ix.try_add(b"\xc9"), 2);
        assert!(ix.rename(0, b"J1"));
        assert_eq!(ix.find(b"node1"), -1);
        assert_eq!(ix.find(b"j1"), 0);
        ix.remove_at(0);
        assert_eq!(ix.find(b"\xc9"), 1);
    }
}
