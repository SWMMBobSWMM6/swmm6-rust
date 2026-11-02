// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/UserFlags.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// User-defined model flags: a schema from [USER_FLAGS] (exact-case names) and values from
// [USER_FLAG_VALUES] keyed by the composite "OBJECTTYPE:OBJECTNAME:FLAGNAME", looked up
// ignoring ASCII case and stored with the first spelling seen. Strings are byte strings;
// std::variant<bool, int, double, std::string> is the enum UserFlagValue. Where upstream
// throws std::out_of_range (get_def, get_value) this returns Err with the same message
// bytes.
//
// The schema index (`def_index_`, an unordered_map) is a HashMap: it is private, and the
// only loop over it (undefine's renumbering) gives the same result in any order.
//
// The value map is NOT a faithful copy of upstream's order. Upstream's is an
// std::unordered_map<std::string, UserFlagValue, CiHash, CiEqual>, and InpWriter.cpp
// iterates all_values() to write [USER_FLAG_VALUES], so the written line order is
// libstdc++'s bucket order (FNV-1a over the upper-cased key, prime bucket counts). Here
// ValueMap iterates in insertion order, which is deterministic but differs from
// upstream's whenever a model has two or more flag values. Lookups, the stored spelling,
// counts and erasures agree with upstream; only iteration order differs.

use crate::core::options::Bytes;
use crate::data::name_index::ci_upper;
use std::collections::HashMap;

/// Value type for a user flag.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum UserFlagType {
    /// YES/NO/TRUE/FALSE/1/0
    Boolean = 0,
    /// Signed integer.
    Integer = 1,
    /// Double-precision floating-point.
    Real = 2,
    /// Arbitrary string.
    String = 3,
}

/// Schema entry for a single user-defined flag (from [USER_FLAGS]).
#[derive(Debug, Clone, PartialEq)]
pub struct UserFlagDef {
    /// Flag identifier (uppercase).
    pub name: Bytes,
    /// Value type.
    pub type_: UserFlagType,
    /// Human-readable description (may be empty).
    pub description: Bytes,
}

// `type` has no member initializer upstream; every upstream construction assigns it.
// Default gives BOOLEAN, what value-initialization (`UserFlagDef{}`) gives.
impl Default for UserFlagDef {
    fn default() -> Self {
        UserFlagDef { name: Vec::new(), type_: UserFlagType::Boolean, description: Vec::new() }
    }
}

/// A concrete flag value (std::variant<bool, int, double, std::string>).
#[derive(Debug, Clone, PartialEq)]
pub enum UserFlagValue {
    Bool(bool),
    Int(i32),
    Double(f64),
    String(Bytes),
}

// A default-constructed std::variant holds its first alternative, value-initialized.
impl Default for UserFlagValue {
    fn default() -> Self {
        UserFlagValue::Bool(false)
    }
}

/// A (object_type, object_name, flag_name) -> value assignment (from [USER_FLAG_VALUES]).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct UserFlagAssignment {
    /// e.g. "NODE", "LINK", "SUBCATCHMENT".
    pub object_type: Bytes,
    /// Object identifier as it appears in the model.
    pub object_name: Bytes,
    /// Must match a name in the schema.
    pub flag_name: Bytes,
    /// Typed runtime value.
    pub value: UserFlagValue,
}

fn fold(s: &[u8]) -> Bytes {
    s.iter().map(|&c| ci_upper(c)).collect()
}

/// `std::string::find(c, from)`: npos (None) when `from` is at or past the end.
fn find_byte(s: &[u8], c: u8, from: usize) -> Option<usize> {
    if from >= s.len() {
        return None;
    }
    s[from..].iter().position(|&b| b == c).map(|p| p + from)
}

/// The composite-key value map (upstream UserFlags::ValueMap). Keys compare ignoring ASCII
/// case and keep the first spelling stored; iteration is insertion order (see the file
/// header: upstream's is libstdc++ bucket order).
#[derive(Debug, Clone, Default)]
pub struct ValueMap {
    entries: Vec<(Bytes, UserFlagValue)>,
    index: HashMap<Bytes, usize>,
}

impl ValueMap {
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn get(&self, key: &[u8]) -> Option<&UserFlagValue> {
        self.index.get(&fold(key)).map(|&i| &self.entries[i].1)
    }

    pub fn contains_key(&self, key: &[u8]) -> bool {
        self.index.contains_key(&fold(key))
    }

    /// `map[key] = value`: an existing equal-ignoring-case key keeps its stored spelling.
    pub fn assign(&mut self, key: Bytes, value: UserFlagValue) {
        let k = fold(&key);
        if let Some(&i) = self.index.get(&k) {
            self.entries[i].1 = value;
        } else {
            self.index.insert(k, self.entries.len());
            self.entries.push((key, value));
        }
    }

    /// `map.erase(key) != 0`.
    pub fn erase(&mut self, key: &[u8]) -> bool {
        match self.index.remove(&fold(key)) {
            None => false,
            Some(i) => {
                self.entries.remove(i);
                for v in self.index.values_mut() {
                    if *v > i {
                        *v -= 1;
                    }
                }
                true
            }
        }
    }

    /// Keep the entries whose stored key satisfies `keep`.
    pub fn retain_keys<F: FnMut(&[u8]) -> bool>(&mut self, mut keep: F) {
        self.entries.retain(|(k, _)| keep(k));
        self.index.clear();
        for (i, (k, _)) in self.entries.iter().enumerate() {
            self.index.insert(fold(k), i);
        }
    }

    pub fn clear(&mut self) {
        self.entries.clear();
        self.index.clear();
    }

    /// (stored key, value) pairs in insertion order.
    pub fn iter(&self) -> impl Iterator<Item = (&Bytes, &UserFlagValue)> {
        self.entries.iter().map(|(k, v)| (k, v))
    }
}

/// Schema definitions plus per-object values.
#[derive(Debug, Clone, Default)]
pub struct UserFlags {
    defs: Vec<UserFlagDef>,
    /// name -> defs index (exact case).
    def_index: HashMap<Bytes, usize>,
    /// composite key -> value (case-insensitive).
    values: ValueMap,
}

impl UserFlags {
    // ---- schema (from [USER_FLAGS])

    /// Insert or replace an exact-case schema entry; replacing keeps its slot and values.
    pub fn define(&mut self, def: UserFlagDef) {
        if let Some(&i) = self.def_index.get(&def.name) {
            self.defs[i] = def;
        } else {
            self.def_index.insert(def.name.clone(), self.defs.len());
            self.defs.push(def);
        }
    }

    /// Remove a schema entry, renumber the ones after it, and remove every value whose key
    /// suffix after the second colon equals `name` exactly (case-sensitive).
    pub fn undefine(&mut self, name: &[u8]) -> bool {
        let idx = match self.def_index.get(name) {
            None => return false,
            Some(&i) => i,
        };
        self.defs.remove(idx);
        self.def_index.remove(name);
        // Reindex definitions that followed the erased slot.
        for v in self.def_index.values_mut() {
            if *v > idx {
                *v -= 1;
            }
        }
        // Remove orphaned per-object values (key = "TYPE:NAME:FLAG"; same left-to-right
        // decomposition as InpWriter).
        self.values.retain_keys(|key| {
            let p1 = find_byte(key, b':', 0);
            let p2 = match p1 {
                None => None,
                Some(p) => find_byte(key, b':', p + 1),
            };
            !matches!(p2, Some(p) if &key[p + 1..] == name)
        });
        true
    }

    pub fn is_defined(&self, name: &[u8]) -> bool {
        self.def_index.contains_key(name)
    }

    /// The exact-case schema entry, or Err (upstream throws std::out_of_range).
    pub fn get_def(&self, name: &[u8]) -> Result<&UserFlagDef, Bytes> {
        match self.def_index.get(name) {
            Some(&i) => Ok(&self.defs[i]),
            None => {
                let mut msg = b"UserFlags: no definition for flag '".to_vec();
                msg.extend_from_slice(name);
                msg.push(b'\'');
                Err(msg)
            }
        }
    }

    /// Schema definitions in insertion order.
    pub fn all_defs(&self) -> &[UserFlagDef] {
        &self.defs
    }

    pub fn def_count(&self) -> usize {
        self.defs.len()
    }

    // ---- per-object values (from [USER_FLAG_VALUES])

    pub fn set(
        &mut self,
        object_type: &[u8],
        object_name: &[u8],
        flag_name: &[u8],
        value: UserFlagValue,
    ) {
        self.values.assign(Self::make_key(object_type, object_name, flag_name), value);
    }

    pub fn unset(&mut self, object_type: &[u8], object_name: &[u8], flag_name: &[u8]) -> bool {
        self.values.erase(&Self::make_key(object_type, object_name, flag_name))
    }

    pub fn has_value(&self, object_type: &[u8], object_name: &[u8], flag_name: &[u8]) -> bool {
        self.values.contains_key(&Self::make_key(object_type, object_name, flag_name))
    }

    /// The stored value, or Err (upstream throws std::out_of_range).
    pub fn get_value(
        &self,
        object_type: &[u8],
        object_name: &[u8],
        flag_name: &[u8],
    ) -> Result<&UserFlagValue, Bytes> {
        match self.values.get(&Self::make_key(object_type, object_name, flag_name)) {
            Some(v) => Ok(v),
            None => {
                let mut msg = b"UserFlags: no value for ".to_vec();
                msg.extend_from_slice(object_type);
                msg.push(b':');
                msg.extend_from_slice(object_name);
                msg.push(b':');
                msg.extend_from_slice(flag_name);
                Err(msg)
            }
        }
    }

    /// A copy of the stored value, or None.
    pub fn try_get_value(
        &self,
        object_type: &[u8],
        object_name: &[u8],
        flag_name: &[u8],
    ) -> Option<UserFlagValue> {
        self.values.get(&Self::make_key(object_type, object_name, flag_name)).cloned()
    }

    pub fn all_values(&self) -> &ValueMap {
        &self.values
    }

    pub fn value_count(&self) -> usize {
        self.values.len()
    }

    // ---- reset

    pub fn clear(&mut self) {
        self.defs.clear();
        self.def_index.clear();
        self.values.clear();
    }

    // ---- composite-key helper (public for tests)

    /// object_type + ':' + object_name + ':' + flag_name, unescaped and unnormalized.
    pub fn make_key(object_type: &[u8], object_name: &[u8], flag_name: &[u8]) -> Bytes {
        let mut key = Vec::with_capacity(object_type.len() + 1 + object_name.len() + 1 + flag_name.len());
        key.extend_from_slice(object_type);
        key.push(b':');
        key.extend_from_slice(object_name);
        key.push(b':');
        key.extend_from_slice(flag_name);
        key
    }
}
