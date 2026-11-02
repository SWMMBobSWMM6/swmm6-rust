// SPDX-License-Identifier: Apache-2.0
//
// libstdc++'s std::map (GCC 13 _Rb_tree) for the operations find and emplace, kept as
// the red-black tree itself. Written for swmm6-rust, 2 November 2026, by Claude
// (Anthropic) for Robert Dickinson, from the documented algorithm (_M_lower_bound,
// _M_get_insert_unique_pos, _Rb_tree_insert_and_rebalance, _Rb_tree_decrement).
//
// Why a tree and not a search: when the key's operator< is not a strict weak order (a
// double key holding NaN compares neither less nor greater than anything), what find()
// returns depends on the tree's shape. Reproducing the insert and rebalance steps
// reproduces that shape, so the answers match upstream's for every key, NaN included.
// For a strict weak order this is an ordinary ordered map.

const NIL: usize = usize::MAX;
const HEADER: usize = 0;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Color {
    Red,
    Black,
}

struct Node<K, V> {
    color: Color,
    parent: usize,
    left: usize,
    right: usize,
    kv: Option<(K, V)>,
}

/// An ordered map with libstdc++'s tree, compared by `less`.
pub struct StdMap<K, V> {
    nodes: Vec<Node<K, V>>,
    less: fn(&K, &K) -> bool,
}

impl<K, V> StdMap<K, V> {
    pub fn new(less: fn(&K, &K) -> bool) -> Self {
        // The header: parent = root, left = leftmost, right = rightmost; red.
        let header = Node { color: Color::Red, parent: NIL, left: HEADER, right: HEADER, kv: None };
        StdMap { nodes: vec![header], less }
    }

    fn key(&self, n: usize) -> &K {
        &self.nodes[n].kv.as_ref().unwrap().0
    }

    fn root(&self) -> usize {
        self.nodes[HEADER].parent
    }

    pub fn len(&self) -> usize {
        self.nodes.len() - 1
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn lower_bound(&self, k: &K) -> usize {
        let mut x = self.root();
        let mut y = HEADER;
        while x != NIL {
            if !(self.less)(self.key(x), k) {
                y = x;
                x = self.nodes[x].left;
            } else {
                x = self.nodes[x].right;
            }
        }
        y
    }

    /// `find(k)`: the mapped value, or None for end().
    pub fn find(&self, k: &K) -> Option<&V> {
        let j = self.lower_bound(k);
        if j == HEADER || (self.less)(k, self.key(j)) {
            None
        } else {
            Some(&self.nodes[j].kv.as_ref().unwrap().1)
        }
    }

    fn decrement(&self, x: usize) -> usize {
        let n = &self.nodes[x];
        if n.color == Color::Red && n.parent != NIL && self.nodes[n.parent].parent == x {
            return n.right; // the header
        }
        if n.left != NIL {
            let mut y = n.left;
            while self.nodes[y].right != NIL {
                y = self.nodes[y].right;
            }
            return y;
        }
        let mut x = x;
        let mut y = n.parent;
        while x == self.nodes[y].left {
            x = y;
            y = self.nodes[y].parent;
        }
        y
    }

    /// `emplace(k, v)`: true when inserted, false when an equivalent key is present.
    pub fn emplace(&mut self, k: K, v: V) -> bool {
        // _M_get_insert_unique_pos
        let mut x = self.root();
        let mut y = HEADER;
        let mut comp = true;
        while x != NIL {
            y = x;
            comp = (self.less)(&k, self.key(x));
            x = if comp { self.nodes[x].left } else { self.nodes[x].right };
        }
        let mut j = y;
        let leftmost = self.nodes[HEADER].left;
        let mut insert_at: Option<usize> = None;
        if comp {
            if j == leftmost {
                insert_at = Some(y);
            } else {
                j = self.decrement(j);
            }
        }
        if insert_at.is_none() {
            if (self.less)(self.key(j), &k) {
                insert_at = Some(y);
            } else {
                return false;
            }
        }
        let p = insert_at.unwrap();
        // _M_insert_node with x == 0
        let insert_left = p == HEADER || (self.less)(&k, self.key(p));
        let z = self.nodes.len();
        self.nodes.push(Node { color: Color::Red, parent: p, left: NIL, right: NIL, kv: Some((k, v)) });
        self.insert_and_rebalance(insert_left, z, p);
        true
    }

    fn rotate_left(&mut self, x: usize) {
        let y = self.nodes[x].right;
        self.nodes[x].right = self.nodes[y].left;
        if self.nodes[y].left != NIL {
            let yl = self.nodes[y].left;
            self.nodes[yl].parent = x;
        }
        self.nodes[y].parent = self.nodes[x].parent;
        let xp = self.nodes[x].parent;
        if x == self.root() {
            self.nodes[HEADER].parent = y;
        } else if x == self.nodes[xp].left {
            self.nodes[xp].left = y;
        } else {
            self.nodes[xp].right = y;
        }
        self.nodes[y].left = x;
        self.nodes[x].parent = y;
    }

    fn rotate_right(&mut self, x: usize) {
        let y = self.nodes[x].left;
        self.nodes[x].left = self.nodes[y].right;
        if self.nodes[y].right != NIL {
            let yr = self.nodes[y].right;
            self.nodes[yr].parent = x;
        }
        self.nodes[y].parent = self.nodes[x].parent;
        let xp = self.nodes[x].parent;
        if x == self.root() {
            self.nodes[HEADER].parent = y;
        } else if x == self.nodes[xp].right {
            self.nodes[xp].right = y;
        } else {
            self.nodes[xp].left = y;
        }
        self.nodes[y].right = x;
        self.nodes[x].parent = y;
    }

    fn insert_and_rebalance(&mut self, insert_left: bool, x0: usize, p: usize) {
        let mut x = x0;
        if insert_left {
            self.nodes[p].left = x;
            if p == HEADER {
                self.nodes[HEADER].parent = x;
                self.nodes[HEADER].right = x;
            } else if p == self.nodes[HEADER].left {
                self.nodes[HEADER].left = x;
            }
        } else {
            self.nodes[p].right = x;
            if p == self.nodes[HEADER].right {
                self.nodes[HEADER].right = x;
            }
        }
        while x != self.root() && self.nodes[self.nodes[x].parent].color == Color::Red {
            let xp = self.nodes[x].parent;
            let xpp = self.nodes[xp].parent;
            if xp == self.nodes[xpp].left {
                let y = self.nodes[xpp].right;
                if y != NIL && self.nodes[y].color == Color::Red {
                    self.nodes[xp].color = Color::Black;
                    self.nodes[y].color = Color::Black;
                    self.nodes[xpp].color = Color::Red;
                    x = xpp;
                } else {
                    if x == self.nodes[xp].right {
                        x = xp;
                        self.rotate_left(x);
                    }
                    let xp = self.nodes[x].parent;
                    self.nodes[xp].color = Color::Black;
                    self.nodes[xpp].color = Color::Red;
                    self.rotate_right(xpp);
                }
            } else {
                let y = self.nodes[xpp].left;
                if y != NIL && self.nodes[y].color == Color::Red {
                    self.nodes[xp].color = Color::Black;
                    self.nodes[y].color = Color::Black;
                    self.nodes[xpp].color = Color::Red;
                    x = xpp;
                } else {
                    if x == self.nodes[xp].left {
                        x = xp;
                        self.rotate_right(x);
                    }
                    let xp = self.nodes[x].parent;
                    self.nodes[xp].color = Color::Black;
                    self.nodes[xpp].color = Color::Red;
                    self.rotate_left(xpp);
                }
            }
        }
        let r = self.root();
        self.nodes[r].color = Color::Black;
    }
}
