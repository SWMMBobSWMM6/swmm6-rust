// Rust side of gate/diff/ucfdump.cpp: the same lines, byte for byte.
#![allow(clippy::needless_range_loop, clippy::field_reassign_with_default)]
use std::io::Write;
use swmm6_rust::core::options::{FlowUnits, SimulationOptions};
use swmm6_rust::core::unit_conversion::*;

fn hx(out: &mut Vec<u8>, v: f64) {
    out.extend_from_slice(format!(" {:016x}", v.to_bits()).as_bytes());
}

fn du(out: &mut Vec<u8>, tag: &str, fu: i32, d: &DisplayUnits) {
    out.extend_from_slice(format!("{} {} {} {}", tag, fu, d.unit_system, d.flow_units).as_bytes());
    for v in [
        d.flow, d.length, d.volume, d.rainfall, d.raindepth, d.evaprate, d.landarea, d.mvol,
        d.landvol,
    ] {
        hx(out, v);
    }
    out.extend_from_slice(b" |");
    for w in [
        d.flow_word, d.length_word, d.depth_word, d.vel_word, d.mvol_word, d.landvol_word,
        d.storage_vol_word,
    ] {
        out.extend_from_slice(w);
        out.extend_from_slice(b"|");
    }
    out.extend_from_slice(b"\n");
    let ts = [
        32.0, 212.0, -40.0, 0.0, -0.0, 98.6, 1e308, -1e308, 5e-324, 33.8, 50.0, f64::INFINITY,
        f64::NEG_INFINITY, f64::NAN,
    ];
    out.extend_from_slice(format!("TEMP {}", fu).as_bytes());
    for t in ts {
        hx(out, d.temperature(t));
    }
    out.extend_from_slice(b"\n");
}

fn main() {
    let mut out: Vec<u8> = Vec::new();
    for q in 0..10 {
        out.extend_from_slice(format!("Ucf {}", q).as_bytes());
        hx(&mut out, UCF[q][0]);
        hx(&mut out, UCF[q][1]);
        out.extend_from_slice(b"\n");
    }
    for q in 0..10 {
        out.extend_from_slice(format!("Ucf_inv {}", q).as_bytes());
        hx(&mut out, UCF_INV[q][0]);
        hx(&mut out, UCF_INV[q][1]);
        out.extend_from_slice(b"\n");
    }
    for f in 0..6 {
        out.extend_from_slice(format!("Qcf {}", f).as_bytes());
        hx(&mut out, QCF[f]);
        hx(&mut out, QCF_INV[f]);
        out.extend_from_slice(b" ");
        out.extend_from_slice(FLOW_UNIT_WORDS[f]);
        out.extend_from_slice(b"\n");
    }
    out.extend_from_slice(b"CONST");
    hx(&mut out, ACRES_TO_FT2);
    hx(&mut out, FT3_TO_MGAL);
    hx(&mut out, SEC_PER_DAY);
    out.extend_from_slice(b"\n");
    out.extend_from_slice(
        format!(
            "Q {} {} {} {} {} {} {} {} {} {} {}\n",
            RAINFALL, RAINDEPTH, EVAPRATE, LENGTH, LANDAREA, VOLUME, WINDSPEED, TEMPERATURE, MASS,
            GWFLOW, FLOW
        )
        .as_bytes(),
    );
    for fu in [-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, i32::MIN, i32::MAX] {
        out.extend_from_slice(format!("US {} {}\n", fu, get_unit_system(fu)).as_bytes());
    }
    let units = [
        FlowUnits::Cfs,
        FlowUnits::Gpm,
        FlowUnits::Mgd,
        FlowUnits::Cms,
        FlowUnits::Lps,
        FlowUnits::Mld,
    ];
    for (fu, &u) in units.iter().enumerate() {
        let fu = fu as i32;
        let mut opts = SimulationOptions::default();
        opts.flow_units = u;
        for q in [-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 100, i32::MIN, i32::MAX] {
            out.extend_from_slice(format!("UCF {} {}", fu, q).as_bytes());
            hx(&mut out, ucf(q, &opts));
            hx(&mut out, ucf_inv(q, &opts));
            out.extend_from_slice(b"\n");
        }
        du(&mut out, "DU", fu, &DisplayUnits::from(&opts));
    }
    du(&mut out, "DU0", -1, &DisplayUnits::default());
    {
        let opts = SimulationOptions::default();
        du(&mut out, "DUdef", -1, &DisplayUnits::from(&opts));
    }
    std::io::stdout().write_all(&out).unwrap();
}
