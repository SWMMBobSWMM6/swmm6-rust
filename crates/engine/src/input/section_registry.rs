// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/SectionRegistry.hpp and SectionRegistry.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Section tags are upper-cased (ASCII) on the way in and on lookup. A custom handler
// shadows a built-in one of the same tag.
//
// `registered_tags` iterates two hash maps in upstream, so its order is libstdc++'s
// bucket order. Here it is built-in tags in registration order, then custom tags not
// already listed, and the order is not promised to match upstream's. Nothing on the run
// path uses it (see PLAN.md, translation hazards).

use std::collections::HashMap;

use crate::core::context::SimulationContext;
use crate::core::options::Bytes;
use crate::input::tokenizer;

/// upstream: std::function<void(SimulationContext&, const std::vector<std::string>&)>
pub type SectionHandler = Box<dyn Fn(&mut SimulationContext, &[Bytes])>;

#[derive(Default)]
pub struct SectionRegistry {
    builtin: HashMap<Bytes, SectionHandler>,
    custom: HashMap<Bytes, SectionHandler>,
    builtin_order: Vec<Bytes>,
    custom_order: Vec<Bytes>,
}

impl SectionRegistry {
    pub fn new() -> Self {
        Self::default()
    }

    fn normalize_tag(tag: &[u8]) -> Bytes {
        tokenizer::to_upper(tag)
    }

    pub fn register_builtin(&mut self, tag: &[u8], handler: impl Fn(&mut SimulationContext, &[Bytes]) + 'static) {
        let k = Self::normalize_tag(tag);
        if self.builtin.insert(k.clone(), Box::new(handler)).is_none() {
            self.builtin_order.push(k);
        }
    }

    pub fn register_custom(&mut self, tag: &[u8], handler: impl Fn(&mut SimulationContext, &[Bytes]) + 'static) {
        let k = Self::normalize_tag(tag);
        if self.custom.insert(k.clone(), Box::new(handler)).is_none() {
            self.custom_order.push(k);
        }
    }

    pub fn has(&self, tag: &[u8]) -> bool {
        let key = Self::normalize_tag(tag);
        self.custom.contains_key(&key) || self.builtin.contains_key(&key)
    }

    pub fn is_custom(&self, tag: &[u8]) -> bool {
        self.custom.contains_key(&Self::normalize_tag(tag))
    }

    pub fn dispatch(&self, tag: &[u8], ctx: &mut SimulationContext, lines: &[Bytes]) {
        let key = Self::normalize_tag(tag);
        if let Some(h) = self.custom.get(&key) {
            h(ctx, lines);
            return;
        }
        if let Some(h) = self.builtin.get(&key) {
            h(ctx, lines);
        }
    }

    pub fn registered_tags(&self) -> Vec<Bytes> {
        let mut tags = self.builtin_order.clone();
        for t in &self.custom_order {
            if !tags.contains(t) {
                tags.push(t.clone());
            }
        }
        tags
    }
}
