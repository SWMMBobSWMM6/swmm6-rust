// SPDX-License-Identifier: Apache-2.0
//
// swmm6-rust: the OpenSWMM Engine 1D sources (HydroCouple/openswmm.engine, branch
// swmm6_rel, commit 88329d6df2f8dd63fe7e7c601954ec2bde3573e6) translated into Rust.
// Each module names the upstream files it derives from and carries their notices.
// This is not an OpenSWMM or HydroCouple release. See NOTICE.

#![forbid(unsafe_code)]

pub mod core;
pub mod data;
pub mod edit;
pub mod hydraulics;
pub mod hydrology;
pub mod input;
pub mod math;
pub mod transport;
