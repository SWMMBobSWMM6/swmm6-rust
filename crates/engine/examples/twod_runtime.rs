// Rust twin of gate/diff/twod_runtime.cpp: the same pseudo-random sequence of operations on the
// translation's RuntimeSources and RuntimeForcings (twod/coupling/runtime_coupling.rs), the same
// printed results and store dumps.
// Build: cargo build --release --features twod --example twod_runtime
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
use std::io::Write;

use swmm6_rust::twod::coupling::runtime_coupling::{RuntimeForcing, RuntimeForcings, RuntimeReceipt, RuntimeSource, RuntimeSources};

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545F4914F6CDD1D)
    }
    fn ri(&mut self, n: i32) -> i32 {
        (self.next() % n as u64) as i32
    }
    fn rd(&mut self) -> f64 {
        (self.ri(4001) - 2000) as f64 / 16.0
    }
    fn runtil(&mut self) -> f64 {
        if self.ri(4) == 0 {
            f64::INFINITY
        } else {
            self.ri(400) as f64 / 4.0
        }
    }
    fn receipt(&mut self) -> RuntimeReceipt {
        let mut r = RuntimeReceipt::default();
        r.requested = self.rd();
        r.applied = self.rd();
        r.water_in = self.rd();
        r.water_out = self.rd();
        r.heat_requested = self.rd();
        r.heat_applied = self.rd();
        r.last_flow = self.rd();
        r.last_dt = self.rd();
        let mut k = self.ri(3);
        while k > 0 {
            let x = self.rd();
            r.requested_mass.push(x);
            k -= 1;
        }
        let mut k = self.ri(3);
        while k > 0 {
            let x = self.rd();
            r.applied_mass.push(x);
            k -= 1;
        }
        r
    }
}

const IDS: [&str; 4] = ["a", "b", "c", "d"];

fn hx(o: &mut Vec<u8>, v: f64) {
    write!(o, " {:016x}", v.to_bits()).unwrap();
}
fn vec(o: &mut Vec<u8>, n: &str, v: &[f64]) {
    write!(o, " {}{}", n, v.len()).unwrap();
    for &x in v {
        hx(o, x);
    }
}
fn receipt(o: &mut Vec<u8>, r: &RuntimeReceipt) {
    write!(o, " rc").unwrap();
    for x in [r.requested, r.applied, r.water_in, r.water_out, r.heat_requested, r.heat_applied, r.last_flow, r.last_dt] {
        hx(o, x);
    }
    vec(o, "rm", &r.requested_mass);
    vec(o, "am", &r.applied_mass);
}

fn dump_s(o: &mut Vec<u8>, s: &RuntimeSources) {
    for r in &s.rows {
        o.extend_from_slice(b"row ");
        o.extend_from_slice(&r.id);
        o.push(b' ');
        o.extend_from_slice(&r.owner);
        write!(o, " {}", r.cell).unwrap();
        hx(o, r.flow);
        hx(o, r.heat);
        write!(o, " {}", r.boundary_edge).unwrap();
        hx(o, r.boundary_head);
        hx(o, r.until);
        write!(o, " {}", if r.enabled { 1 } else { 0 }).unwrap();
        vec(o, "c", &r.concentrations);
        vec(o, "r", &r.rates);
        receipt(o, &r.receipt);
        o.push(b'\n');
    }
    write!(o, "cells {}:", s.cells.len()).unwrap();
    for c in &s.cells {
        write!(o, " [").unwrap();
        for k in c {
            write!(o, " {}", k).unwrap();
        }
        write!(o, " ]").unwrap();
    }
    write!(o, "\nfind").unwrap();
    for id in IDS {
        for cell in -1..9 {
            match s.find(id.as_bytes(), cell) {
                Some(k) => write!(o, " {}", k).unwrap(),
                None => write!(o, " -1").unwrap(),
            }
        }
    }
    o.push(b'\n');
}

fn index_of(f: &RuntimeForcings, r: Option<&RuntimeForcing>) -> isize {
    match r {
        None => -1,
        Some(r) => f.rows.iter().position(|x| std::ptr::eq(x, r)).map_or(-1, |k| k as isize),
    }
}

fn dump_f(o: &mut Vec<u8>, f: &RuntimeForcings) {
    for r in &f.rows {
        write!(o, "frow {} {} {}", r.channel, r.cell, r.mode).unwrap();
        hx(o, r.rate);
        hx(o, r.until);
        write!(o, " {} ", if r.enabled { 1 } else { 0 }).unwrap();
        o.extend_from_slice(&r.owner);
        vec(o, "c", &r.concentrations);
        receipt(o, &r.receipt);
        o.push(b'\n');
    }
    write!(o, "fcells {}:", f.cells.len()).unwrap();
    for c in &f.cells {
        write!(o, " {},{},{},{}", c[0], c[1], c[2], c[3]).unwrap();
    }
    write!(o, " uniform {},{},{},{}\nffind", f.uniform[0], f.uniform[1], f.uniform[2], f.uniform[3]).unwrap();
    for ch in 1..=4 {
        for cell in -1..8 {
            match f.find(ch, cell) {
                Some(k) => write!(o, " {}", k).unwrap(),
                None => write!(o, " -1").unwrap(),
            }
        }
    }
    o.push(b'\n');
}

fn main() {
    let mut g = Rng(0x9E3779B97F4A7C15);
    let mut s = RuntimeSources::default();
    let mut snap = RuntimeSources::default();
    let mut f = RuntimeForcings::default();
    let mut o: Vec<u8> = Vec::new();
    for step in 0..6000 {
        let op = g.ri(16);
        write!(o, "op {} {}:", step, op).unwrap();
        match op {
            0..=2 => {
                let mut v = RuntimeSource::default();
                v.id = IDS[g.ri(4) as usize].as_bytes().to_vec();
                v.owner = if g.ri(2) != 0 { b"own".to_vec() } else { Vec::new() };
                v.cell = g.ri(9);
                v.flow = g.rd();
                v.heat = g.rd();
                v.boundary_edge = g.ri(3) - 1;
                if g.ri(2) != 0 {
                    v.boundary_head = g.rd();
                }
                v.until = g.runtil();
                v.enabled = g.ri(5) != 0;
                let mut k = g.ri(3);
                while k > 0 {
                    let c = g.rd();
                    v.concentrations.push(c);
                    let r = g.rd();
                    v.rates.push(r);
                    k -= 1;
                }
                if g.ri(2) != 0 {
                    v.receipt = g.receipt();
                }
                o.extend_from_slice(b" replace ");
                o.extend_from_slice(&v.id);
                writeln!(o, " {}", v.cell).unwrap();
                s.replace(v);
            }
            3 => {
                let id: &str = if g.ri(5) == 4 { "zz" } else { IDS[g.ri(4) as usize] };
                let cell = g.ri(10) - 1;
                let r = s.clear(id.as_bytes(), cell);
                writeln!(o, " clear {} {} -> {}", id, cell, if r { 1 } else { 0 }).unwrap();
            }
            4 => {
                let count = 3 + g.ri(6);
                write!(o, " rebuild {} ->", count).unwrap();
                match s.rebuild(count) {
                    Ok(()) => writeln!(o, " ok").unwrap(),
                    Err(w) => {
                        o.extend_from_slice(b" out_of_range ");
                        o.extend_from_slice(&w);
                        o.push(b'\n');
                    }
                }
            }
            5 => {
                writeln!(o, " prune").unwrap();
                s.prune();
            }
            6 => {
                let now = g.ri(400) as f64 / 4.0;
                write!(o, " expire").unwrap();
                hx(&mut o, now);
                let r = s.expire(now);
                writeln!(o, " -> {}", if r { 1 } else { 0 }).unwrap();
            }
            7 => {
                let now = g.ri(400) as f64 / 4.0;
                write!(o, " nextTime").unwrap();
                hx(&mut o, now);
                write!(o, " ->").unwrap();
                hx(&mut o, s.next_time(now));
                writeln!(o, " active {}", if s.active() { 1 } else { 0 }).unwrap();
            }
            8 => {
                writeln!(o, " snapshot").unwrap();
                snap = s.clone();
            }
            9 => {
                writeln!(o, " syncReceipts").unwrap();
                s.sync_receipts(&snap);
            }
            10 | 11 => {
                let mut r = RuntimeForcing::default();
                r.channel = 1 + g.ri(4);
                r.cell = g.ri(8) - 1;
                r.mode = g.ri(3);
                r.rate = g.rd();
                r.until = g.runtil();
                r.enabled = r.mode != 0;
                r.owner = if g.ri(2) != 0 { b"own".to_vec() } else { Vec::new() };
                let mut k = g.ri(3);
                while k > 0 {
                    let c = g.rd();
                    r.concentrations.push(c);
                    k -= 1;
                }
                if g.ri(2) != 0 {
                    r.receipt = g.receipt();
                }
                writeln!(o, " freplace {} {}", r.channel, r.cell).unwrap();
                f.replace(r);
            }
            12 => {
                let count = 3 + g.ri(6);
                write!(o, " frebuild {} ->", count).unwrap();
                match f.rebuild(count) {
                    Ok(()) => writeln!(o, " ok").unwrap(),
                    Err(w) => {
                        o.extend_from_slice(b" out_of_range ");
                        o.extend_from_slice(&w);
                        o.push(b'\n');
                    }
                }
            }
            13 => {
                let now = g.ri(400) as f64 / 4.0;
                write!(o, " fexpire").unwrap();
                hx(&mut o, now);
                f.expire(now);
                o.push(b'\n');
            }
            14 => {
                let now = g.ri(400) as f64 / 4.0;
                write!(o, " fnextTime").unwrap();
                hx(&mut o, now);
                write!(o, " ->").unwrap();
                hx(&mut o, f.next_time(now));
                writeln!(o, " active {}", if f.active() { 1 } else { 0 }).unwrap();
            }
            _ => {
                write!(o, " fget").unwrap();
                for ch in 1..=4 {
                    for cell in -1..9 {
                        let k = index_of(&f, f.get(ch, cell));
                        write!(o, " {}", k).unwrap();
                        hx(&mut o, f.rate(ch, cell, 0.5 * cell as f64));
                    }
                }
                o.push(b'\n');
            }
        }
        dump_s(&mut o, &s);
        dump_f(&mut o, &f);
    }
    std::io::stdout().lock().write_all(&o).unwrap();
}
