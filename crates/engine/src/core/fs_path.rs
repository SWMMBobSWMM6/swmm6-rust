// SPDX-License-Identifier: Apache-2.0
//
// The lexical operations of std::filesystem::path that core/PathResolver.cpp relies on,
// as libstdc++ (GCC 13, POSIX) performs them: construction (component split),
// operator/=, has_filename, filename, remove_filename, parent_path, has_relative_path,
// has_root_directory, is_absolute, lexically_normal, lexically_relative,
// lexically_proximate and generic_string. Written for swmm6-rust, 2 November 2026, by
// Claude (Anthropic) for Robert Dickinson, from the documented behaviour of those
// members, and checked against the reference build's libstdc++ by a differential test.
//
// A path is its byte string plus the component list libstdc++ keeps: when the path has
// one component (or none) the list is empty and `t` is that component's kind; otherwise
// `t` is Multi and each component carries its bytes, kind and offset. The offsets matter:
// lexically_normal truncates the string at a component's offset. libstdc++ can leave a
// Multi path with a single component (lexically_normal's ".." step on "/a/"); that state
// is kept as it is, since later calls see it.
//
// POSIX has no root-name, so RootName never occurs; it is kept for the member tests that
// name it.

use crate::core::options::Bytes;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PathType {
    Filename,
    RootName,
    RootDir,
    Multi,
}

#[derive(Debug, Clone)]
struct Cmpt {
    s: Bytes,
    t: PathType,
    pos: usize,
}

#[derive(Debug, Clone)]
pub struct FsPath {
    pathname: Bytes,
    cmpts: Vec<Cmpt>,
    t: PathType,
}

impl Default for FsPath {
    fn default() -> Self {
        FsPath { pathname: Bytes::new(), cmpts: Vec::new(), t: PathType::Filename }
    }
}

fn is_dot_s(s: &[u8]) -> bool {
    s.len() == 1 && s[0] == b'.'
}

fn is_dotdot_s(s: &[u8]) -> bool {
    s.len() == 2 && s[0] == b'.' && s[1] == b'.'
}

impl FsPath {
    /// `path(std::string)`: the bytes, split into components.
    pub fn new(s: &[u8]) -> FsPath {
        let mut p = FsPath { pathname: s.to_vec(), cmpts: Vec::new(), t: PathType::Filename };
        p.split_cmpts();
        p
    }

    /// A component as a path of its own kind (libstdc++'s private `path(str, type)`).
    fn single(s: &[u8], t: PathType) -> FsPath {
        FsPath { pathname: s.to_vec(), cmpts: Vec::new(), t }
    }

    /// `_M_split_cmpts` with the POSIX parser (one separator, '/').
    fn split_cmpts(&mut self) {
        self.cmpts.clear();
        let input = &self.pathname;
        if input.is_empty() {
            self.t = PathType::Filename;
            return;
        }
        if input.len() == 1 && input[0] == b'/' {
            self.t = PathType::RootDir;
            return;
        }
        let len = input.len();
        let find_not_sep = |from: usize| -> Option<usize> { (from..len).find(|&k| input[k] != b'/') };
        let mut buf: Vec<Cmpt> = Vec::new();
        // root_path()
        let mut pos: Option<usize> = Some(0);
        let mut last_type = PathType::Multi;
        if input[0] == b'/' {
            buf.push(Cmpt { s: vec![b'/'], t: PathType::RootDir, pos: 0 });
            pos = find_not_sep(1);
            last_type = PathType::RootDir;
        }
        // next() until invalid
        loop {
            let last_pos = pos;
            let mut f: Option<Cmpt> = None;
            if let Some(p0) = pos {
                match find_not_sep(p0) {
                    Some(start) => {
                        let end = (start..len).find(|&k| input[k] == b'/');
                        let stop = end.unwrap_or(len);
                        f = Some(Cmpt { s: input[start..stop].to_vec(), t: PathType::Filename, pos: start });
                        pos = end;
                    }
                    None => {
                        pos = None;
                        if last_type == PathType::Filename || (last_pos == Some(0) && !input.is_empty()) {
                            f = Some(Cmpt { s: Bytes::new(), t: PathType::Filename, pos: len });
                        }
                    }
                }
            }
            match f {
                Some(c) => {
                    last_type = c.t;
                    buf.push(c);
                }
                None => break,
            }
        }
        if buf.len() == 1 {
            self.t = buf[0].t;
            return;
        }
        self.t = PathType::Multi;
        self.cmpts = buf;
    }

    pub fn native(&self) -> &[u8] {
        &self.pathname
    }

    pub fn empty(&self) -> bool {
        self.pathname.is_empty()
    }

    pub fn clear(&mut self) {
        self.pathname.clear();
        self.split_cmpts();
    }

    /// The elements iteration yields: the components of a Multi path, else the path itself
    /// (nothing for an empty path).
    fn elements(&self) -> Vec<FsPath> {
        if self.t == PathType::Multi {
            self.cmpts.iter().map(|c| FsPath::single(&c.s, c.t)).collect()
        } else if self.empty() {
            Vec::new()
        } else {
            vec![FsPath::single(&self.pathname, self.t)]
        }
    }

    pub fn has_filename(&self) -> bool {
        if self.empty() {
            return false;
        }
        if self.t == PathType::Filename {
            return !self.pathname.is_empty();
        }
        if self.t == PathType::Multi {
            if *self.pathname.last().unwrap() == b'/' {
                return false;
            }
            let b = self.cmpts.last().unwrap();
            return b.t == PathType::Filename && !b.s.is_empty();
        }
        false
    }

    pub fn filename(&self) -> FsPath {
        if self.empty() {
            return FsPath::default();
        }
        if self.t == PathType::Filename {
            return self.clone();
        }
        if self.t == PathType::Multi {
            if *self.pathname.last().unwrap() == b'/' {
                return FsPath::default();
            }
            let last = self.cmpts.last().unwrap();
            if last.t == PathType::Filename {
                return FsPath::single(&last.s, last.t);
            }
        }
        FsPath::default()
    }

    pub fn has_relative_path(&self) -> bool {
        if self.t == PathType::Filename && !self.pathname.is_empty() {
            return true;
        }
        if !self.cmpts.is_empty() {
            let mut it = 0;
            if self.cmpts[it].t == PathType::RootName {
                it += 1;
            }
            if it < self.cmpts.len() && self.cmpts[it].t == PathType::RootDir {
                it += 1;
            }
            if it < self.cmpts.len() && !self.cmpts[it].s.is_empty() {
                return true;
            }
        }
        false
    }

    pub fn has_root_directory(&self) -> bool {
        if self.t == PathType::RootDir {
            return true;
        }
        if !self.cmpts.is_empty() {
            let mut it = 0;
            if self.cmpts[it].t == PathType::RootName {
                it += 1;
            }
            if it < self.cmpts.len() && self.cmpts[it].t == PathType::RootDir {
                return true;
            }
        }
        false
    }

    /// POSIX: a path is absolute when it has a root directory.
    pub fn is_absolute(&self) -> bool {
        self.has_root_directory()
    }

    pub fn remove_filename(&mut self) {
        if self.t == PathType::Multi {
            if !self.cmpts.is_empty() {
                let last = self.cmpts.len() - 1;
                if self.cmpts[last].t == PathType::Filename && !self.cmpts[last].s.is_empty() {
                    self.pathname.truncate(self.cmpts[last].pos);
                    let prev_t = self.cmpts[last - 1].t;
                    if prev_t == PathType::RootDir || prev_t == PathType::RootName {
                        self.cmpts.pop();
                        if self.cmpts.len() == 1 {
                            self.t = self.cmpts[0].t;
                            self.cmpts.clear();
                        }
                    } else {
                        self.cmpts[last].s.clear();
                        self.cmpts[last].t = PathType::Filename;
                    }
                }
            }
        } else if self.t == PathType::Filename {
            self.clear();
        }
    }

    pub fn parent_path(&self) -> FsPath {
        if !self.has_relative_path() {
            self.clone()
        } else if self.cmpts.len() >= 2 {
            let parent = &self.cmpts[self.cmpts.len() - 2];
            let l = parent.pos + parent.s.len();
            FsPath::new(&self.pathname[..l])
        } else {
            FsPath::default()
        }
    }

    /// `*this /= p` (the POSIX operator/=).
    pub fn div_assign(&mut self, p: &FsPath) {
        if p.is_absolute() || self.empty() {
            *self = p.clone();
            return;
        }
        let sep: &[u8] = if self.has_filename() {
            b"/"
        } else if p.empty() {
            return;
        } else {
            b""
        };
        let orig_pathlen = self.pathname.len();
        let orig_type = self.t;
        self.pathname.extend_from_slice(sep);
        let basepos = self.pathname.len();
        self.pathname.extend_from_slice(&p.pathname);

        let mut cm: Vec<Cmpt> = if orig_type == PathType::Multi {
            let mut c = std::mem::take(&mut self.cmpts);
            if c.last().map_or(false, |b| b.s.is_empty()) {
                c.pop();
            }
            c
        } else if orig_pathlen != 0 {
            vec![Cmpt { s: self.pathname[..orig_pathlen].to_vec(), t: orig_type, pos: 0 }]
        } else {
            Vec::new()
        };
        if p.t == PathType::Multi {
            for c in &p.cmpts {
                cm.push(Cmpt { s: c.s.clone(), t: PathType::Filename, pos: c.pos + basepos });
            }
        } else if !p.empty() || !sep.is_empty() {
            cm.push(Cmpt { s: p.pathname.clone(), t: p.t, pos: basepos });
        }
        self.cmpts = cm;
        self.t = PathType::Multi;
    }

    /// `lhs / rhs`.
    pub fn join(&self, rhs: &FsPath) -> FsPath {
        let mut r = self.clone();
        r.div_assign(rhs);
        r
    }

    pub fn lexically_normal(&self) -> FsPath {
        let mut ret = FsPath::default();
        if self.empty() {
            return ret;
        }
        for p in self.elements() {
            if is_dotdot_s(&p.pathname) {
                if ret.has_filename() {
                    if !is_dotdot_s(&ret.filename().pathname) {
                        ret.remove_filename();
                    } else {
                        ret.div_assign(&p);
                    }
                } else if !ret.has_relative_path() {
                    if !ret.has_root_directory() {
                        ret.div_assign(&p);
                    }
                } else {
                    let ei = ret.cmpts.len() - 2;
                    let elem_has_filename = ret.cmpts[ei].t == PathType::Filename && !ret.cmpts[ei].s.is_empty();
                    if elem_has_filename && !is_dotdot_s(&ret.cmpts[ei].s) {
                        if ei == 0 {
                            ret.clear();
                        } else {
                            let epos = ret.cmpts[ei].pos;
                            ret.pathname.truncate(epos);
                            ret.cmpts.pop();
                            if ret.cmpts[ei - 1].t == PathType::Filename {
                                ret.cmpts[ei].s.clear();
                                ret.cmpts[ei].t = PathType::Filename;
                            } else {
                                ret.cmpts.pop();
                            }
                        }
                    } else {
                        ret.div_assign(&p);
                    }
                }
            } else if is_dot_s(&p.pathname) {
                ret.div_assign(&FsPath::default());
            } else {
                ret.div_assign(&p);
            }
        }
        if ret.cmpts.len() >= 2 {
            let n = ret.cmpts.len();
            if ret.cmpts[n - 1].s.is_empty() && is_dotdot_s(&ret.cmpts[n - 2].s) {
                ret = ret.parent_path();
            }
        } else if ret.empty() {
            ret = FsPath::new(b".");
        }
        ret
    }

    /// Element equality as path::compare() == 0 decides it for two components.
    fn element_eq(a: &FsPath, b: &FsPath) -> bool {
        a.pathname == b.pathname || (a.t == PathType::RootDir && b.t == PathType::RootDir)
    }

    pub fn lexically_relative(&self, base: &FsPath) -> FsPath {
        let mut ret = FsPath::default();
        // root_name() is empty on both sides on POSIX.
        if self.is_absolute() != base.is_absolute() {
            return ret;
        }
        if !self.has_root_directory() && base.has_root_directory() {
            return ret;
        }
        let a_el = self.elements();
        let b_el = base.elements();
        let mut i = 0;
        while i < a_el.len() && i < b_el.len() && Self::element_eq(&a_el[i], &b_el[i]) {
            i += 1;
        }
        if i == a_el.len() && i == b_el.len() {
            ret = FsPath::new(b".");
        } else {
            let mut n: i32 = 0;
            for p in &b_el[i..] {
                if is_dotdot_s(&p.pathname) {
                    n -= 1;
                } else if !p.empty() && !is_dot_s(&p.pathname) {
                    n += 1;
                }
            }
            if n == 0 && (i == a_el.len() || a_el[i].empty()) {
                ret = FsPath::new(b".");
            } else if n >= 0 {
                let dotdot = FsPath::new(b"..");
                while n > 0 {
                    n -= 1;
                    ret.div_assign(&dotdot);
                }
                for a in &a_el[i..] {
                    ret.div_assign(a);
                }
            }
        }
        ret
    }

    pub fn lexically_proximate(&self, base: &FsPath) -> FsPath {
        let rel = self.lexically_relative(base);
        if rel.empty() {
            self.clone()
        } else {
            rel
        }
    }

    /// generic_string(): the elements joined by single '/' (a root directory is "/").
    pub fn generic_string(&self) -> Bytes {
        if self.t == PathType::RootDir {
            return vec![b'/'];
        }
        let mut out = Bytes::with_capacity(self.pathname.len());
        let mut add_slash = false;
        for e in self.elements() {
            if add_slash {
                out.push(b'/');
            }
            out.extend_from_slice(&e.pathname);
            add_slash = e.t == PathType::Filename;
        }
        out
    }
}
