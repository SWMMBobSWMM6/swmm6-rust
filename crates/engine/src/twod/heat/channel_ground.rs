// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/heat/ChannelGround.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// ChannelGroundState: the MESH_CONTACTS channel/ground thermal contacts. Initializing them is
// phase 2b (SurfaceRouter2D::initializeChannelGround refuses it before it starts), so the state
// is only ever cleared here: active() is false. The private members the inactive state carries
// are not modelled; only their emptiness is observable on the 2a path.

#![allow(clippy::all)]

#[derive(Debug, Clone, Default)]
pub struct ChannelGroundState {
    /// The number of channel/ground contacts (contacts_.size() upstream): always 0 in 2a.
    pub n_contacts: usize,
}

impl ChannelGroundState {
    pub fn clear(&mut self) {
        *self = ChannelGroundState::default();
    }
    pub fn active(&self) -> bool {
        self.n_contacts != 0
    }
}
