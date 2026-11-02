// SPDX-License-Identifier: Apache-2.0
//
// The iteration order of libstdc++'s std::unordered_map (GCC 13 _Hashtable with unique
// keys and _Prime_rehash_policy) under insertion only. Written for swmm6-rust, 2 November
// 2026, by Claude (Anthropic) for Robert Dickinson, from the documented algorithm
// (_M_insert_unique_node, _M_insert_bucket_begin, _M_rehash_aux for unique keys,
// _Prime_rehash_policy::_M_need_rehash and _M_next_bkt). The prime table is the head of
// libstdc++'s std::__detail::__prime_list, read from the gate machine's libstdc++.so.
//
// Why: where upstream loops over an unordered_map and the order reaches a result (PLAN.md,
// "Iterating an unordered_map"), the order is the hash table's internal node list. That
// order depends only on each key's hash code, the insertion sequence and the bucket-count
// growth, so it is reproduced here without storing the keys: the caller keeps its keys
// in insertion order, does its own (unique-key) lookup, and calls `insert` with the hash
// code of each NEW key. `order()` gives the insertion indices in iteration order.
//
// Not modelled: erase, rehash/reserve calls by the caller, a non-default max_load_factor,
// and bucket counts past the last listed prime (beyond 2,019,773,507 buckets).

/// libstdc++ __prime_list, the entries this table can reach.
const PRIMES: &[u64] = &[
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 103, 109, 113, 127, 137, 139, 149, 157, 167, 179,
    193, 199, 211, 227, 241, 257, 277, 293, 313, 337, 359, 383, 409, 439, 467, 503, 541, 577, 619, 661, 709, 761, 823, 887, 953, 1031, 1109, 1193,
    1289, 1381, 1493, 1613, 1741, 1879, 2029, 2179, 2357, 2549, 2753, 2971, 3209, 3469, 3739, 4027, 4349, 4703, 5087, 5503, 5953, 6427, 6949, 7517,
    8123, 8783, 9497, 10273, 11113, 12011, 12983, 14033, 15173, 16411, 17749, 19183, 20753, 22447, 24281, 26267, 28411, 30727, 33223, 35933, 38873,
    42043, 45481, 49201, 53201, 57557, 62233, 67307, 72817, 78779, 85229, 92203, 99733, 107897, 116731, 126271, 136607, 147793, 159871, 172933,
    187091, 202409, 218971, 236897, 256279, 277261, 299951, 324503, 351061, 379787, 410857, 444487, 480881, 520241, 562841, 608903, 658753, 712697,
    771049, 834181, 902483, 976369, 1056323, 1142821, 1236397, 1337629, 1447153, 1565659, 1693859, 1832561, 1982627, 2144977, 2320627, 2510653,
    2716249, 2938679, 3179303, 3439651, 3721303, 4026031, 4355707, 4712381, 5098259, 5515729, 5967347, 6456007, 6984629, 7556579, 8175383, 8844859,
    9569143, 10352717, 11200489, 12117689, 13109983, 14183539, 15345007, 16601593, 17961079, 19431899, 21023161, 22744717, 24607243, 26622317,
    28802401, 31160981, 33712729, 36473443, 39460231, 42691603, 46187573, 49969847, 54061849, 58488943, 63278561, 68460391, 74066549, 80131819,
    86693767, 93793069, 101473717, 109783337, 118773397, 128499677, 139022417, 150406843, 162723577, 176048909, 190465427, 206062531, 222936881,
    241193053, 260944219, 282312799, 305431229, 330442829, 357502601, 386778277, 418451333, 452718089, 489790921, 529899637, 573292817, 620239453,
    671030513, 725980837, 785430967, 849749479, 919334987, 994618837, 1076067617, 1164186217, 1259520799, 1362662261, 1474249943, 1594975441,
    1725587117, 1866894511, 2019773507,
];

/// `_M_buckets[b]`: the node before bucket b's first node: the list head or a node.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Before {
    Head,
    Node(usize),
}

const NIL: usize = usize::MAX;

/// Node order of a libstdc++ unordered_map with unique keys, max_load_factor 1.
pub struct UnorderedMapOrder {
    hash: Vec<u64>,
    next: Vec<usize>,
    head_next: usize,
    buckets: Vec<Option<Before>>,
    element_count: usize,
    next_resize: usize,
}

impl Default for UnorderedMapOrder {
    fn default() -> Self {
        Self::new()
    }
}

impl UnorderedMapOrder {
    /// A default-constructed map: one (single) bucket, next_resize 0.
    pub fn new() -> Self {
        UnorderedMapOrder { hash: Vec::new(), next: Vec::new(), head_next: NIL, buckets: vec![None], element_count: 0, next_resize: 0 }
    }

    pub fn len(&self) -> usize {
        self.element_count
    }

    pub fn is_empty(&self) -> bool {
        self.element_count == 0
    }

    /// _Prime_rehash_policy::_M_next_bkt
    fn next_bkt(&mut self, n: usize) -> usize {
        const FAST_BKT: [u8; 14] = [2, 2, 2, 3, 5, 5, 7, 7, 11, 11, 11, 11, 13, 13];
        if n < FAST_BKT.len() {
            if n == 0 {
                return 1;
            }
            self.next_resize = (FAST_BKT[n] as f64 * 1.0).floor() as usize;
            return FAST_BKT[n] as usize;
        }
        // lower_bound over [prime_list + 6, last_prime)
        let last = PRIMES.len() - 1;
        let mut k = 6;
        while k < last && (PRIMES[k] as usize) < n {
            k += 1;
        }
        if k == last {
            self.next_resize = usize::MAX;
        } else {
            self.next_resize = (PRIMES[k] as f64 * 1.0).floor() as usize;
        }
        PRIMES[k] as usize
    }

    /// _Prime_rehash_policy::_M_need_rehash
    fn need_rehash(&mut self, n_bkt: usize, n_elt: usize, n_ins: usize) -> Option<usize> {
        if n_elt + n_ins > self.next_resize {
            let floor_hint = if self.next_resize != 0 { 0 } else { 11 };
            let min_bkts = std::cmp::max(n_elt + n_ins, floor_hint) as f64 / 1.0;
            if min_bkts >= n_bkt as f64 {
                let want = std::cmp::max(min_bkts.floor() as usize + 1, n_bkt * 2);
                return Some(self.next_bkt(want));
            }
            self.next_resize = (n_bkt as f64 * 1.0).floor() as usize;
            return None;
        }
        None
    }

    /// _M_rehash_aux(n, true_type)
    fn rehash(&mut self, n: usize) {
        let mut nb: Vec<Option<Before>> = vec![None; n];
        let mut p = self.head_next;
        self.head_next = NIL;
        let mut bbegin_bkt = 0usize;
        while p != NIL {
            let nxt = self.next[p];
            let bkt = (self.hash[p] % n as u64) as usize;
            match nb[bkt] {
                None => {
                    self.next[p] = self.head_next;
                    self.head_next = p;
                    nb[bkt] = Some(Before::Head);
                    if self.next[p] != NIL {
                        nb[bbegin_bkt] = Some(Before::Node(p));
                    }
                    bbegin_bkt = bkt;
                }
                Some(b) => {
                    let after = self.after(b);
                    self.next[p] = after;
                    self.set_after(b, p);
                }
            }
            p = nxt;
        }
        self.buckets = nb;
    }

    fn after(&self, b: Before) -> usize {
        match b {
            Before::Head => self.head_next,
            Before::Node(i) => self.next[i],
        }
    }

    fn set_after(&mut self, b: Before, v: usize) {
        match b {
            Before::Head => self.head_next = v,
            Before::Node(i) => self.next[i] = v,
        }
    }

    /// Insert a key known not to be present, with hash code `code`; returns its
    /// insertion index (the caller's key index).
    pub fn insert(&mut self, code: u64) -> usize {
        let node = self.hash.len();
        self.hash.push(code);
        self.next.push(NIL);
        if let Some(n) = self.need_rehash(self.buckets.len(), self.element_count, 1) {
            self.rehash(n);
        }
        let bkt = (code % self.buckets.len() as u64) as usize;
        // _M_insert_bucket_begin
        match self.buckets[bkt] {
            Some(b) => {
                let after = self.after(b);
                self.next[node] = after;
                self.set_after(b, node);
            }
            None => {
                self.next[node] = self.head_next;
                self.head_next = node;
                if self.next[node] != NIL {
                    let nb = (self.hash[self.next[node]] % self.buckets.len() as u64) as usize;
                    self.buckets[nb] = Some(Before::Node(node));
                }
                self.buckets[bkt] = Some(Before::Head);
            }
        }
        self.element_count += 1;
        node
    }

    /// Insertion indices in iteration (begin() to end()) order.
    pub fn order(&self) -> Vec<usize> {
        let mut v = Vec::with_capacity(self.element_count);
        let mut p = self.head_next;
        while p != NIL {
            v.push(p);
            p = self.next[p];
        }
        v
    }
}

/// std::hash<int>: the value widened to size_t (sign-extended).
pub fn hash_int(k: i32) -> u64 {
    k as i64 as u64
}

/// openswmm::CiHash: 64-bit FNV-1a over the ASCII-uppercased bytes.
pub fn ci_hash(s: &[u8]) -> u64 {
    let mut h: u64 = 14695981039346656037;
    for &c in s {
        h ^= c.to_ascii_uppercase() as u64;
        h = h.wrapping_mul(1099511628211);
    }
    h
}

/// Iteration order of an `unordered_map<std::string, V, CiHash, CiEqual>` filled only by
/// inserting `entries` in this order, each key new (a later spelling of a key would have
/// been a lookup, not an insert): the entries in the order a range-for visits them.
pub fn ci_map_iteration_order<V: Clone>(entries: &[(Vec<u8>, V)]) -> Vec<(Vec<u8>, V)> {
    let mut m = UnorderedMapOrder::new();
    for (k, _) in entries {
        m.insert(ci_hash(k));
    }
    m.order().into_iter().map(|i| entries[i].clone()).collect()
}
