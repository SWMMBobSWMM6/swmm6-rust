// Rust side of gate/diff/umaporder.cpp: the same insertion sequences through
// core::std_unordered_map::UnorderedMapOrder with std::hash<int> and CiHash codes.
use swmm6_rust::core::std_unordered_map::{ci_hash, hash_int, UnorderedMapOrder};
fn main() {
    let mut s: u32 = 12345;
    let mut rnd = || {
        s = s.wrapping_mul(1103515245).wrapping_add(12345);
        (s >> 8) & 0xffffff
    };
    for trial in 0..200 {
        let n = if trial % 7 == 0 { rnd() % 300 } else { rnd() % 40 } as i32;
        let mut mi = UnorderedMapOrder::new();
        let mut keys: Vec<i32> = vec![];
        for _ in 0..n {
            let k = (rnd() % if trial % 3 != 0 { 60 } else { 100000 }) as i32;
            if !keys.contains(&k) {
                mi.insert(hash_int(k));
                keys.push(k);
            }
        }
        print!("I {}:", keys.len());
        for i in mi.order() {
            print!(" {}", i);
        }
        println!();
        let mut ms = UnorderedMapOrder::new();
        let mut sk: Vec<Vec<u8>> = vec![];
        for i in 0..n {
            let mut k = format!("uh{}", rnd() % 500).into_bytes();
            if i % 5 == 0 {
                k[0] = b'U';
            }
            if !sk.iter().any(|x| x.eq_ignore_ascii_case(&k)) {
                ms.insert(ci_hash(&k));
                sk.push(k);
            }
        }
        print!("S {}:", sk.len());
        for i in ms.order() {
            print!(" {}", i);
        }
        println!();
    }
}
