// Rust side of gate/diff/tabfile.cpp: TableData.cpp's file reader and validate_table, and
// MultiColumnSeriesFile, over every file in the list (argv[1]) that `tabfile gen` wrote.
use std::io::Write;
use swmm6_rust::core::options::Bytes;
use swmm6_rust::data::table::*;
use swmm6_rust::input::multi_column_series_file::*;

type O = Vec<u8>;

fn hx(o: &mut O, v: f64) {
    write!(o, "{:016x} ", v.to_bits()).unwrap();
}
fn bs(o: &mut O, s: &[u8]) {
    write!(o, "{}:", s.len()).unwrap();
    o.extend_from_slice(s);
    o.push(b' ');
}
fn fnv(v: &[f64], mut h: u64) -> u64 {
    for d in v {
        let u = d.to_bits();
        for i in 0..8 {
            h ^= (u >> (8 * i)) & 0xff;
            h = h.wrapping_mul(1099511628211);
        }
    }
    h
}
const FNV0: u64 = 1469598103934665603;

fn blk(o: &mut O, name: &str, b: &TableBlock) {
    writeln!(o, "{} nc={} frs={} n={} ny={}", name, b.num_cols, b.file_row_start, b.x.len(), b.y.len()).unwrap();
    o.extend_from_slice(b" x ");
    for &d in &b.x {
        hx(o, d);
    }
    o.push(b'\n');
    o.extend_from_slice(b" y ");
    for &d in &b.y {
        hx(o, d);
    }
    o.push(b'\n');
}
fn cache_digest(o: &mut O, t: &Table) {
    let c = &t.cache;
    write!(o, "cache nc={} frs={} n={} ny={} h={:016x}", c.num_cols, c.file_row_start, c.x.len(), c.y.len(), fnv(&c.y, fnv(&c.x, FNV0))).unwrap();
    if !c.x.is_empty() {
        o.extend_from_slice(b" first=");
        hx(o, c.x[0]);
        o.extend_from_slice(b"last=");
        hx(o, c.x[c.x.len() - 1]);
    }
    o.push(b'\n');
}
fn fields(o: &mut O, t: &Table) {
    writeln!(o, "fb={} fh={} nc={} tr={} ncr={} dso={}", t.is_file_based as i32, t.file_handle.is_some() as i32, t.num_cols, t.total_rows, t.num_cache_rows, t.data_start_offset).unwrap();
    o.extend_from_slice(b"ro");
    for off in &t.row_offsets {
        write!(o, " {}", off).unwrap();
    }
    o.push(b'\n');
    o.extend_from_slice(b"cid ");
    for s in &t.column_ids {
        bs(o, s);
    }
    o.push(b'\n');
    let mut m: Vec<(&Bytes, &i32)> = t.column_map.iter().collect();
    m.sort();
    o.extend_from_slice(b"cmap ");
    for (k, v) in m {
        bs(o, k);
        write!(o, "={} ", v).unwrap();
    }
    o.push(b'\n');
    o.extend_from_slice(b"xmin ");
    hx(o, t.x_min);
    o.extend_from_slice(b"xmax ");
    hx(o, t.x_max);
    o.extend_from_slice(b"dxmin ");
    hx(o, t.dx_min);
    o.push(b'\n');
    blk(o, "first", &t.first_boundary);
    blk(o, "last", &t.last_boundary);
    cache_digest(o, t);
    writeln!(o, "cursor {} {}", t.cursor.index, t.cursor.direction).unwrap();
}
fn val(o: &mut O, v: &TableValidation, t: &Table) {
    write!(o, "valid={} ne={} nw={} dxmin ", v.valid as i32, v.errors.len(), v.warnings.len()).unwrap();
    hx(o, t.dx_min);
    o.push(b'\n');
    for e in &v.errors {
        o.extend_from_slice(b" E ");
        bs(o, e);
        o.push(b'\n');
    }
    for w in &v.warnings {
        o.extend_from_slice(b" W ");
        bs(o, w);
        o.push(b'\n');
    }
}

fn queries(t: &Table, full: bool) -> Vec<f64> {
    let mut qs = Vec::new();
    let (lo, hi) = (t.x_min, t.x_max);
    if hi > lo {
        let span = hi - lo;
        for k in 0..=40 {
            qs.push((lo - 0.05 * span) + (k as f64) * (1.1 * span / 40.0));
        }
        for k in (0..=40).rev() {
            qs.push((lo + (k as f64) * (span / 40.0)) + span / 97.0);
        }
    } else {
        qs.extend_from_slice(&[-1.0e6, -1.0, 0.0, lo, lo + 0.5, lo + 1.0, 1.0e9]);
    }
    qs.push(lo);
    qs.push(hi);
    qs.push(lo - 1.0);
    qs.push(hi + 1.0);
    let fx = &t.first_boundary.x;
    for i in 0..fx.len().min(6) {
        qs.push(fx[i]);
    }
    let mut i = 0;
    while i + 1 < fx.len() && i < 5 {
        qs.push(0.5 * (fx[i] + fx[i + 1]));
        i += 1;
    }
    let lx = &t.last_boundary.x;
    let nlb = lx.len();
    for i in (if nlb > 6 { nlb - 6 } else { 0 })..nlb {
        qs.push(lx[i]);
    }
    qs.push(f64::NAN);
    qs.push(f64::INFINITY);
    qs.push(f64::NEG_INFINITY);
    if full {
        return qs;
    }
    qs.iter().step_by(3).copied().collect()
}
fn cols_for(nc: i32) -> Vec<i32> {
    if nc <= 8 {
        (-1..=nc).collect()
    } else {
        vec![-1, 0, 1, nc / 2, nc - 1, nc]
    }
}
fn lookups(o: &mut O, t: &mut Table, full: bool) {
    let qs = queries(t, full);
    for col in cols_for(t.num_cols) {
        writeln!(o, "col {}", col).unwrap();
        for &q in &qs {
            let a = table_lookup_column(t, col, q);
            hx(o, q);
            hx(o, a);
            write!(o, "{} {} {} {} ", t.cursor.index, t.cursor.direction, t.cache.file_row_start, t.cache.x.len()).unwrap();
            let b = table_step_column(t, col, q);
            hx(o, b);
            writeln!(o, "{} {} {} {}", t.cursor.index, t.cursor.direction, t.cache.file_row_start, t.cache.x.len()).unwrap();
        }
    }
    cache_digest(o, t);
}
fn loads(o: &mut O, t: &mut Table) {
    let tr = t.total_rows;
    let starts = [0, 1, 4095, 4096, 4097, tr / 2, if tr > 0 { tr - 1 } else { 0 }, tr, tr + 1, tr + 9000];
    for s in starts {
        let r = table_load_cache(t, s);
        write!(o, "load {} -> {} ", s, r).unwrap();
        cache_digest(o, t);
    }
}
fn run_table(o: &mut O, path: &[u8], br: usize, full: bool) {
    writeln!(o, "-- table br={}", br).unwrap();
    let mut t = Table::new(b"TS1", TableType::Timeseries);
    t.file_path.set(path);
    let ok = table_open_file(&mut t, br);
    writeln!(o, "open={}", ok as i32).unwrap();
    fields(o, &t);
    let v = validate_table(&mut t);
    val(o, &v, &t);
    t.x = vec![0.0]; // so empty() is false and the file-backed checks run
    let v = validate_table(&mut t);
    val(o, &v, &t);
    t.type_ = TableType::CurveStorage;
    let v = validate_table(&mut t);
    val(o, &v, &t);
    t.type_ = TableType::Timeseries;
    t.x.clear();
    lookups(o, &mut t, full);
    loads(o, &mut t);
    fields(o, &t);
    t.file_handle = None;
}

fn print_parsed(o: &mut O, pf: &ParsedSeriesFile) {
    write!(o, "headers {} ", pf.headers.len()).unwrap();
    for h in &pf.headers {
        bs(o, h);
    }
    o.push(b'\n');
    o.extend_from_slice(b"first ");
    hx(o, pf.first_date);
    o.extend_from_slice(b"last ");
    hx(o, pf.last_date);
    writeln!(o, "unparsed {} ndates {} ncols {}", pf.unparsed_rows, pf.dates.len(), pf.columns.len()).unwrap();
    for i in 0..pf.col_first_date.len() {
        write!(o, "cs {} ", i).unwrap();
        hx(o, pf.col_first_date[i]);
        hx(o, pf.col_last_date[i]);
        writeln!(o, "{} {}", pf.col_periods_precip[i], pf.col_unparsed_cells[i]).unwrap();
    }
    o.extend_from_slice(b"dates ");
    for &d in &pf.dates {
        hx(o, d);
    }
    o.push(b'\n');
    for (i, c) in pf.columns.iter().enumerate() {
        write!(o, "c{} {} ", i, c.len()).unwrap();
        for &d in c {
            hx(o, d);
        }
        o.push(b'\n');
    }
}
fn run_mcsf(o: &mut O, path: &[u8]) {
    o.extend_from_slice(b"-- mcsf\n");
    writeln!(o, "looks={}", looks_like_multicolumn_series_file(path) as i32).unwrap();
    let mut pf = ParsedSeriesFile::default();
    pf.headers.push(b"stale".to_vec());
    pf.unparsed_rows = 77;
    let mut errs: Vec<Bytes> = vec![b"pre".to_vec()];
    let mut st = SeriesFileStatus::FormatFailed;
    let ok = parse_multicolumn_series_file(path, &mut pf, &mut errs, Some(&mut st));
    writeln!(o, "ok={} st={} ne={}", ok as i32, st as i32, errs.len()).unwrap();
    for e in &errs {
        o.extend_from_slice(b" E ");
        bs(o, e);
        o.push(b'\n');
    }
    print_parsed(o, &pf);
    o.extend_from_slice(b"find");
    for h in &pf.headers {
        write!(o, " {} {} {}", pf.find_column(h), pf.find_column(&h.to_ascii_uppercase()), pf.find_column(&h.to_ascii_lowercase())).unwrap();
    }
    writeln!(o, " | {} {} fdc {}", pf.find_column(b"nosuch"), pf.find_column(b""), pf.first_data_column()).unwrap();
    let mut pf2 = ParsedSeriesFile::default();
    let mut e2: Vec<Bytes> = Vec::new();
    let ok2 = parse_multicolumn_series_file(path, &mut pf2, &mut e2, None);
    writeln!(o, "ok2={} ne2={} same={}", ok2 as i32, e2.len(), (fnv(&pf2.dates, FNV0) == fnv(&pf.dates, FNV0)) as i32).unwrap();
    let mut cache = MultiColumnFileCache::default();
    let mut ce: Vec<Bytes> = Vec::new();
    let mut cs = SeriesFileStatus::FormatFailed;
    let p1 = cache.get_or_parse(path, &mut ce, Some(&mut cs)).map(|p| (p as *const ParsedSeriesFile as usize, p.dates.len(), fnv(&p.dates, FNV0)));
    write!(o, "cache1 {} st={} ne={} pc={}", p1.is_some() as i32, cs as i32, ce.len(), cache.parse_count()).unwrap();
    if let Some((_, n, h)) = p1 {
        write!(o, " n={} h={:016x}", n, h).unwrap();
    }
    o.push(b'\n');
    let p2 = cache.get_or_parse(path, &mut ce, None).map(|p| p as *const ParsedSeriesFile as usize);
    writeln!(o, "cache2 {} same={} ne={} pc={}", p2.is_some() as i32, (p1.map(|t| t.0) == p2) as i32, ce.len(), cache.parse_count()).unwrap();
    let mut px = path.to_vec();
    px.extend_from_slice(b"/x");
    let p3 = cache.get_or_parse(&px, &mut ce, Some(&mut cs)).is_some();
    writeln!(o, "cache3 {} st={} ne={} pc={}", p3 as i32, cs as i32, ce.len(), cache.parse_count()).unwrap();
}

fn inmem(o: &mut O) {
    o.extend_from_slice(b"== INMEM\n");
    let (nan, inf) = (f64::NAN, f64::INFINITY);
    let cases: Vec<(&str, TableType, Vec<f64>, Vec<f64>, i32)> = vec![
        ("E", TableType::Timeseries, vec![], vec![], 1),
        ("DUP", TableType::Timeseries, vec![1.0, 2.0, 2.0, 3.0], vec![1.0, 2.0, 3.0, 4.0], 1),
        ("SEQ", TableType::Timeseries, vec![1.0, 3.0, 2.0], vec![1.0, 2.0, 3.0], 1),
        ("SEQ2", TableType::Timeseries, vec![45000.75, 45000.25], vec![1.0, 2.0], 1),
        ("CSEQ", TableType::CurveRating, vec![0.0, 1.0, 0.5], vec![1.0, 2.0, 3.0], 1),
        ("NAN", TableType::Timeseries, vec![1.0, nan, 3.0], vec![1.0, 2.0, 3.0], 1),
        ("INF", TableType::Timeseries, vec![1.0, 2.0, inf], vec![1.0, 2.0, 3.0], 1),
        ("MIS", TableType::Timeseries, vec![1.0, 2.0], vec![1.0, 2.0, 3.0], 2),
        ("ONE", TableType::Timeseries, vec![40000.0], vec![7.0, 8.0], 2),
        ("OK2", TableType::Timeseries, vec![40000.0, 40000.5, 40000.5, 40001.25], vec![1.0, 10.0, 2.0, 20.0, 3.0, 30.0, 4.0, 40.0], 2),
    ];
    for (id, ty, x, y, nc) in cases {
        let mut t = Table::new(id.as_bytes(), ty);
        t.x = x;
        t.y = y;
        t.num_cols = nc;
        writeln!(o, "case {}", id).unwrap();
        let v = validate_table(&mut t);
        val(o, &v, &t);
        if nc as usize * t.x.len() != t.y.len() {
            continue;
        }
        let qs = [39999.0, 40000.0, 40000.25, 40000.5, 40000.75, 40001.25, 40002.0, nan, -inf, inf, 2.5, 0.5];
        for col in -1..=nc {
            for &q in &qs {
                let a = table_lookup_column(&mut t, col, q);
                hx(o, a);
                write!(o, "{} {} ", t.cursor.index, t.cursor.direction).unwrap();
                let b = table_step_column(&mut t, col, q);
                hx(o, b);
                writeln!(o, "{} {}", t.cursor.index, t.cursor.direction).unwrap();
            }
        }
    }
}
fn probes(o: &mut O) {
    o.extend_from_slice(b"== PROBES\n");
    let dts: [&str; 30] = [
        "2020-01-02", "2020-01-02 03:04", "2020-01-02T03:04:05", "2020-01-02  T 3:4", "2020-01-02x03:04",
        "2020-1-2 25:61:61", "01/02/2020", "01/02/2020 3", "01/02/2020 03:04", "01/02/2020 03:04:05 PM",
        "01/02/2020 12:00 am", "01/02/2020 12:00 PM", "01/02/2020 13:00 PM", "01/02/2020 00:30 AM",
        "1/2/2020 7 pm", "13/45/2020", "2/29/2021 1:00", "-5-3-2", "", "   ", "x", "01/02", "2020-01",
        "99999999999/1/2020", "1/2/2020 99999:00", " \t01/02/2020 04:05", "01/02/2020 04:05AM", "01/02/2020 AM",
        "2020-01-02 03:04 PM", "01/02/2020 3:4:5 XM",
    ];
    for d in dts {
        let mut out = -7.0;
        let r = parse_series_datetime(d.as_bytes(), &mut out);
        bs(o, d.as_bytes());
        write!(o, "{} ", r as i32).unwrap();
        hx(o, out);
        o.push(b'\n');
    }
    let toks = [
        "rain.csv", "rain.csv:G1", "C:\\dir\\f.csv", "C:\\dir\\f.csv:G2", "/data/Rich:2024/f.csv",
        "/data/Rich:2024/f.csv:col", "f.csv:", ":col", "1:x", "a:b:c", "dir:sub/f", "x:y\\z", "",
    ];
    for tk in toks {
        let (mut p, mut c) = (b"old".to_vec(), b"old".to_vec());
        let r = split_series_file_token(tk.as_bytes(), &mut p, &mut c);
        bs(o, tk.as_bytes());
        write!(o, "{} ", r as i32).unwrap();
        bs(o, &p);
        bs(o, &c);
        o.push(b'\n');
    }
}

fn main() {
    let list = std::fs::read(std::env::args().nth(1).expect("usage: tabfile LIST")).expect("read list");
    let mut paths: Vec<&[u8]> = Vec::new();
    let mut p = 0usize;
    while p < list.len() {
        let e = list[p..].iter().position(|&c| c == b'\n').map(|k| p + k).unwrap_or(list.len());
        paths.push(&list[p..e]);
        p = e + 1;
    }
    let mut o: O = Vec::new();
    for path in paths {
        o.extend_from_slice(b"== FILE ");
        bs(&mut o, path);
        o.push(b'\n');
        run_table(&mut o, path, 128, true);
        run_table(&mut o, path, 3, false);
        run_table(&mut o, path, 0, false);
        run_mcsf(&mut o, path);
    }
    inmem(&mut o);
    probes(&mut o);
    writeln!(o, "total parses {}", multicolumn_parse_count_total()).unwrap();
    std::io::stdout().write_all(&o).unwrap();
}
