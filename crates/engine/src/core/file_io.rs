// SPDX-License-Identifier: Apache-2.0
//
// Paths. Upstream keeps paths as std::string (bytes) and, on Windows, converts UTF-8 to
// wide characters (core/FileIO.hpp, utf8_path). On Unix a path is its bytes.
// Written for swmm6-rust, 2 November 2026, by Claude (Anthropic) for Robert Dickinson.

use std::path::PathBuf;

#[cfg(unix)]
pub fn path_from_bytes(b: &[u8]) -> PathBuf {
    use std::os::unix::ffi::OsStrExt;
    PathBuf::from(std::ffi::OsStr::from_bytes(b))
}

#[cfg(not(unix))]
pub fn path_from_bytes(b: &[u8]) -> PathBuf {
    PathBuf::from(String::from_utf8_lossy(b).into_owned())
}
