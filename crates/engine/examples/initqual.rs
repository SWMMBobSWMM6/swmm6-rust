// Rust side of gate/diff/initqual.cpp: open the deck in the translated engine as the harness
// does (report and output "/dev/null", no input plugin library), print open_rc, run one init
// leg, then dump the context, the members open() sets and the quality members.
// Usage: initqual deck.inp quality|geometry|massbal [side-file]
// In massbal mode the side file (written by `gate/diff/initqual deck massref side`, run first)
// carries upstream's hydrology solver state, which is installed into the translation's
// hydrology solvers (runoff_, lid_, snow_, groundwater_, through the `gate` feature's
// SwmmEngine::gate_hydrology) before init_quality, init_geometry and init_mass_balance (see
// the notes in gate/diff/initqual.cpp). As upstream's harness does, the LID groups are
// rebuilt from the members msxlid::initializeUnits and storedVolume() read, and
// storedVolume() is computed from them, not copied.
// Build: cargo build --release --features gate --example initqual
#[path = "common/ctxdump.rs"]
mod ctxdump;
use std::io::Write;
use swmm6_rust::core::engine::SwmmEngine;
use swmm6_rust::hydrology::lid::{LIDGroupSoA, LIDType};

/// The side file read in the C++ harness's order; None for a missing or incomplete file.
struct Side<'a> {
    tok: std::str::SplitAsciiWhitespace<'a>,
}

impl<'a> Side<'a> {
    fn tag(&mut self, t: &str) -> Option<()> {
        if self.tok.next()? == t {
            Some(())
        } else {
            None
        }
    }
    fn int(&mut self, t: &str) -> Option<i32> {
        self.tag(t)?;
        self.tok.next()?.parse().ok()
    }
    fn bits(&mut self) -> Option<f64> {
        Some(f64::from_bits(u64::from_str_radix(self.tok.next()?, 16).ok()?))
    }
    fn vec(&mut self, t: &str) -> Option<Vec<f64>> {
        self.tag(t)?;
        let n: usize = self.tok.next()?.parse().ok()?;
        let mut v = Vec::with_capacity(n);
        for _ in 0..n {
            v.push(self.bits()?);
        }
        Some(v)
    }
}

fn lid_type(t: i32) -> LIDType {
    match t {
        0 => LIDType::BioCell,
        1 => LIDType::RainGarden,
        2 => LIDType::GreenRoof,
        3 => LIDType::InfilTrench,
        4 => LIDType::PermPavement,
        5 => LIDType::RainBarrel,
        6 => LIDType::VegSwale,
        7 => LIDType::RoofDiscon,
        _ => panic!("LIDType {}", t),
    }
}

/// Install the side file's state; Err(rc) when it is missing, incomplete (-2) or not ok.
fn install(e: &mut SwmmEngine, path: &str) -> Result<i32, i32> {
    let text = std::fs::read_to_string(path).map_err(|_| -2)?;
    let mut s = Side { tok: text.split_ascii_whitespace() };
    let ok = s.int("ok").ok_or(-2)?;
    let rc = s.int("rc").ok_or(-2)?;
    if ok == 0 {
        return Err(rc);
    }
    let r = (|| {
        let area = s.vec("runoff_area")?;
        s.tag("lid_stored_volume")?;
        let lid_v = s.bits()?;
        let en = s.vec("snow_enabled")?;
        let (wsnow, fw, farea) = (s.vec("snow_wsnow")?, s.vec("snow_fw")?, s.vec("snow_fArea")?);
        let (td, ld, th, po) = (s.vec("gw_total_depth")?, s.vec("gw_lower_depth")?, s.vec("gw_theta")?, s.vec("gw_porosity")?);
        let ng = s.int("lid_groups")?;
        let mut groups = Vec::new();
        for _ in 0..ng {
            let mut g = LIDGroupSoA::default();
            g.type_ = lid_type(s.int("lidg_type")?);
            g.count = s.int("lidg_count")?;
            g.subcatch_idx = s.vec("lidg_subcatch_idx")?.iter().map(|&x| x as i32).collect();
            g.control_idx = s.vec("lidg_control_idx")?.iter().map(|&x| x as i32).collect();
            g.stor_covered = s.vec("lidg_stor_covered")?.iter().map(|&x| x as i32).collect();
            g.area = s.vec("lidg_area")?;
            g.surf_store = s.vec("lidg_surf_store")?;
            g.soil_thick = s.vec("lidg_soil_thick")?;
            g.stor_void = s.vec("lidg_stor_void")?;
            g.pave_void = s.vec("lidg_pave_void")?;
            g.pave_imperv_frac = s.vec("lidg_pave_imperv_frac")?;
            g.surf_void_frac = s.vec("lidg_surf_void_frac")?;
            g.surf_side_slope = s.vec("lidg_surf_side_slope")?;
            g.unit_width = s.vec("lidg_unit_width")?;
            g.surf_depth = s.vec("lidg_surf_depth")?;
            g.soil_moist = s.vec("lidg_soil_moist")?;
            g.stor_depth = s.vec("lidg_stor_depth")?;
            g.pave_depth = s.vec("lidg_pave_depth")?;
            groups.push(g);
        }
        s.tag("end")?;
        Some((area, lid_v, [en, wsnow, fw, farea], [td, ld, th, po], groups))
    })();
    let (area, _lid_v, snow, gw, groups) = r.ok_or(-2)?;
    let [en, wsnow, fw, farea] = snow;
    let [td, ld, th, po] = gw;
    let h = e.gate_hydrology();
    h.runoff_.soa_.area = area;
    h.lid_.groups_ = groups; // initializeUnits and storedVolume read groups_ only
    {
        let sn = &mut h.snow_.soa_;
        sn.enabled = en.iter().map(|&x| x as u8).collect();
        sn.wsnow = wsnow;
        sn.fw = fw;
        sn.fArea = farea;
    }
    {
        let g = &mut h.groundwater_.soa_;
        g.total_depth = td;
        g.lower_depth = ld;
        g.theta = th;
        g.porosity = po;
    }
    Ok(rc)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mode = if args.len() > 2 { args[2].as_str() } else { "quality" };
    let mut e = SwmmEngine::new();
    let rc = e.open(Some(args[1].as_bytes()), Some(b"/dev/null"), Some(b"/dev/null"), None, false);
    // Buffered like the C++ twin's piped stdout: a terminate loses everything printed.
    let mut out: Vec<u8> = Vec::new();
    writeln!(out, "open_rc {}", rc).unwrap();
    if rc == 0 {
        let r = match mode {
            "quality" => e.gate_init_quality(),
            "geometry" => {
                e.gate_init_geometry();
                Ok(())
            }
            "massbal" => {
                let r = install(&mut e, &args[3]);
                writeln!(out, "massbal side rc {}", match r { Ok(c) | Err(c) => c }).unwrap();
                if r.is_ok() {
                    e.gate_init_quality().map(|()| {
                        e.gate_init_geometry();
                        e.gate_init_mass_balance();
                    })
                } else {
                    writeln!(out, "massbal skipped").unwrap();
                    Ok(())
                }
            }
            _ => panic!("mode {}", mode),
        };
        if let Err(t) = r {
            eprintln!("terminate called after throwing an instance of '{}'\n  what():  {}", t.exception, t.what);
            std::process::abort();
        }
    }
    out.extend_from_slice(&ctxdump::dump(e.context()));
    out.extend_from_slice(&ctxdump::engdump(&e));
    out.extend_from_slice(&ctxdump::qualdump(&e));
    std::io::stdout().lock().write_all(&out).unwrap();
}
