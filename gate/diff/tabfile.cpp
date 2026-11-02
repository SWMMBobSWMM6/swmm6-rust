// Differential for data/TableData.cpp (validate_table and the paged external-file reader)
// and input/MultiColumnSeriesFile.{hpp,cpp}, against upstream's compiled engine.
// Rust side: crates/engine/examples/tabfile.rs.
//
//   U=/root/work/swmm6src/upstream/src/engine B=/root/work/swmm6src/build
//   g++ -std=c++20 -O3 -ffp-contract=off -mavx2 -mfma -DNDEBUG -DSWMM_OS_LINUX -I$U -I$B/include \
//       tabfile.cpp $B/src/engine/libopenswmm.engine.so -Wl,-rpath,$B/src/engine -o tabfile
//   D=$(mktemp -d); ./tabfile gen $D > $D/files.txt      # adversarial files + full list
//   ./tabfile run $D/files.txt > tabfile_cpp.txt
//   cargo run --release --example tabfile -- $D/files.txt > tabfile_rs.txt
//   cmp tabfile_cpp.txt tabfile_rs.txt
//
// `gen` snapshots every .dat/.rain/.csv/.tsf/.tsv/.txt (not CMakeLists.txt) under the
// corpus and upstream's tests and examples into DIR/snap (DIR/snap/MANIFEST maps them back),
// writes the adversarial files, and lists all of them, then a missing file, a directory
// and an empty path.
#include "data/TableData.hpp"
#include "input/MultiColumnSeriesFile.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <numeric>
#include <string>
#include <vector>
using namespace openswmm;
namespace fs = std::filesystem;

// ---------------------------------------------------------------- output helpers
static void hx(double v) { uint64_t u; std::memcpy(&u, &v, 8); std::printf("%016llx ", (unsigned long long)u); }
static void bs(const std::string& s) { std::printf("%zu:", s.size()); std::fwrite(s.data(), 1, s.size(), stdout); std::printf(" "); }
static void nl() { std::printf("\n"); }
static uint64_t fnv(const std::vector<double>& v, uint64_t h = 1469598103934665603ULL) {
    for (double d : v) { uint64_t u; std::memcpy(&u, &d, 8); for (int i = 0; i < 8; ++i) { h ^= (u >> (8 * i)) & 0xff; h *= 1099511628211ULL; } }
    return h;
}
static void blk(const char* name, const TableBlock& b) {
    std::printf("%s nc=%d frs=%zu n=%zu ny=%zu\n", name, b.num_cols, b.file_row_start, b.x.size(), b.y.size());
    std::printf(" x "); for (double d : b.x) hx(d); nl();
    std::printf(" y "); for (double d : b.y) hx(d); nl();
}
static void cache_digest(const Table& t) {
    const TableBlock& c = t.cache;
    std::printf("cache nc=%d frs=%zu n=%zu ny=%zu h=%016llx", c.num_cols, c.file_row_start, c.x.size(), c.y.size(),
                (unsigned long long)fnv(c.y, fnv(c.x)));
    if (!c.x.empty()) { std::printf(" first="); hx(c.x.front()); std::printf("last="); hx(c.x.back()); }
    nl();
}
static void fields(const Table& t) {
    std::printf("fb=%d fh=%d nc=%d tr=%zu ncr=%zu dso=%ld\n", (int)t.is_file_based, (int)(t.file_handle != nullptr),
                t.num_cols, t.total_rows, t.num_cache_rows, t.data_start_offset);
    std::printf("ro"); for (long o : t.row_offsets) std::printf(" %ld", o); nl();
    std::printf("cid "); for (auto& s : t.column_ids) bs(s); nl();
    std::vector<std::pair<std::string, int>> m(t.column_map.begin(), t.column_map.end());
    std::sort(m.begin(), m.end());
    std::printf("cmap "); for (auto& kv : m) { bs(kv.first); std::printf("=%d ", kv.second); } nl();
    std::printf("xmin "); hx(t.x_min); std::printf("xmax "); hx(t.x_max); std::printf("dxmin "); hx(t.dx_min); nl();
    blk("first", t.first_boundary); blk("last", t.last_boundary); cache_digest(t);
    std::printf("cursor %d %d\n", t.cursor.index, t.cursor.direction);
}
static void val(const TableValidation& v, const Table& t) {
    std::printf("valid=%d ne=%zu nw=%zu dxmin ", (int)v.valid, v.errors.size(), v.warnings.size()); hx(t.dx_min); nl();
    for (auto& e : v.errors) { std::printf(" E "); bs(e); nl(); }
    for (auto& w : v.warnings) { std::printf(" W "); bs(w); nl(); }
}

// ---------------------------------------------------------------- table runs
static std::vector<double> queries(const Table& t, bool full) {
    std::vector<double> qs;
    double lo = t.x_min, hi = t.x_max;
    if (hi > lo) {
        double span = hi - lo;
        for (int k = 0; k <= 40; ++k) qs.push_back((lo - 0.05 * span) + (double)k * (1.1 * span / 40.0));
        for (int k = 40; k >= 0; --k) qs.push_back((lo + (double)k * (span / 40.0)) + span / 97.0);
    } else {
        double s[] = {-1.0e6, -1.0, 0.0, lo, lo + 0.5, lo + 1.0, 1.0e9};
        qs.insert(qs.end(), s, s + 7);
    }
    qs.push_back(lo); qs.push_back(hi); qs.push_back(lo - 1.0); qs.push_back(hi + 1.0);
    for (size_t i = 0; i < t.first_boundary.x.size() && i < 6; ++i) qs.push_back(t.first_boundary.x[i]);
    for (size_t i = 0; i + 1 < t.first_boundary.x.size() && i < 5; ++i)
        qs.push_back(0.5 * (t.first_boundary.x[i] + t.first_boundary.x[i + 1]));
    size_t nlb = t.last_boundary.x.size();
    for (size_t i = (nlb > 6 ? nlb - 6 : 0); i < nlb; ++i) qs.push_back(t.last_boundary.x[i]);
    qs.push_back(std::nan("")); qs.push_back(INFINITY); qs.push_back(-INFINITY);
    if (full) return qs;
    std::vector<double> r;
    for (size_t i = 0; i < qs.size(); i += 3) r.push_back(qs[i]);
    return r;
}
static std::vector<int> cols_for(int nc) {
    std::vector<int> c;
    if (nc <= 8) { for (int i = -1; i <= nc; ++i) c.push_back(i); }
    else { int s[] = {-1, 0, 1, nc / 2, nc - 1, nc}; c.assign(s, s + 6); }
    return c;
}
static void lookups(Table& t, bool full) {
    std::vector<double> qs = queries(t, full);
    for (int col : cols_for(t.num_cols)) {
        std::printf("col %d\n", col);
        for (double q : qs) {
            double a = table_lookup_column(t, col, q);
            hx(q); hx(a); std::printf("%d %d %zu %zu ", t.cursor.index, t.cursor.direction, t.cache.file_row_start, t.cache.x.size());
            double b = table_step_column(t, col, q);
            hx(b); std::printf("%d %d %zu %zu\n", t.cursor.index, t.cursor.direction, t.cache.file_row_start, t.cache.x.size());
        }
    }
    cache_digest(t);
}
static void loads(Table& t) {
    size_t tr = t.total_rows;
    size_t starts[] = {0, 1, 4095, 4096, 4097, tr / 2, tr ? tr - 1 : 0, tr, tr + 1, tr + 9000};
    for (size_t s : starts) { size_t r = table_load_cache(t, s); std::printf("load %zu -> %zu ", s, r); cache_digest(t); }
}
static void run_table(const std::string& path, size_t br, bool full) {
    std::printf("-- table br=%zu\n", br);
    Table t; t.id = "TS1"; t.type = TableType::TIMESERIES; t.file_path = path;
    bool ok = table_open_file(t, br);
    std::printf("open=%d\n", (int)ok);
    fields(t);
    TableValidation v = validate_table(t); val(v, t);
    t.x.assign(1, 0.0); // so empty() is false and the file-backed checks run
    v = validate_table(t); val(v, t);
    t.type = TableType::CURVE_STORAGE; v = validate_table(t); val(v, t);
    t.type = TableType::TIMESERIES;
    t.x.clear();
    lookups(t, full);
    loads(t);
    fields(t);
    if (t.file_handle) std::fclose(t.file_handle);
    t.file_handle = nullptr;
}

// ---------------------------------------------------------------- MultiColumnSeriesFile runs
static std::string upper(std::string s) { for (auto& c : s) c = (char)std::toupper((unsigned char)c); return s; }
static std::string lower(std::string s) { for (auto& c : s) c = (char)std::tolower((unsigned char)c); return s; }
static void print_parsed(const input::ParsedSeriesFile& pf) {
    std::printf("headers %zu ", pf.headers.size()); for (auto& h : pf.headers) bs(h); nl();
    std::printf("first "); hx(pf.first_date); std::printf("last "); hx(pf.last_date); std::printf("unparsed %ld ndates %zu ncols %zu\n", pf.unparsed_rows, pf.dates.size(), pf.columns.size());
    for (size_t i = 0; i < pf.col_first_date.size(); ++i) {
        std::printf("cs %zu ", i); hx(pf.col_first_date[i]); hx(pf.col_last_date[i]);
        std::printf("%ld %ld\n", pf.col_periods_precip[i], pf.col_unparsed_cells[i]);
    }
    std::printf("dates "); for (double d : pf.dates) hx(d); nl();
    for (size_t i = 0; i < pf.columns.size(); ++i) { std::printf("c%zu %zu ", i, pf.columns[i].size()); for (double d : pf.columns[i]) hx(d); nl(); }
}
static void run_mcsf(const std::string& path) {
    std::printf("-- mcsf\n");
    std::printf("looks=%d\n", (int)input::looks_like_multicolumn_series_file(path));
    input::ParsedSeriesFile pf; pf.headers.push_back("stale"); pf.unparsed_rows = 77;
    std::vector<std::string> errs{"pre"};
    input::SeriesFileStatus st = input::SeriesFileStatus::FORMAT_FAILED;
    bool ok = input::parse_multicolumn_series_file(path, pf, errs, &st);
    std::printf("ok=%d st=%d ne=%zu\n", (int)ok, (int)st, errs.size());
    for (auto& e : errs) { std::printf(" E "); bs(e); nl(); }
    print_parsed(pf);
    std::printf("find");
    for (auto& h : pf.headers) std::printf(" %d %d %d", pf.find_column(h), pf.find_column(upper(h)), pf.find_column(lower(h)));
    std::printf(" | %d %d fdc %d\n", pf.find_column("nosuch"), pf.find_column(""), pf.first_data_column());
    input::ParsedSeriesFile pf2; std::vector<std::string> e2;
    bool ok2 = input::parse_multicolumn_series_file(path, pf2, e2, nullptr);
    std::printf("ok2=%d ne2=%zu same=%d\n", (int)ok2, e2.size(), (int)(fnv(pf2.dates) == fnv(pf.dates)));
    input::MultiColumnFileCache cache; std::vector<std::string> ce;
    input::SeriesFileStatus cs = input::SeriesFileStatus::FORMAT_FAILED;
    const input::ParsedSeriesFile* p1 = cache.get_or_parse(path, ce, &cs);
    std::printf("cache1 %d st=%d ne=%zu pc=%d", (int)(p1 != nullptr), (int)cs, ce.size(), cache.parse_count());
    if (p1) std::printf(" n=%zu h=%016llx", p1->dates.size(), (unsigned long long)fnv(p1->dates));
    nl();
    const input::ParsedSeriesFile* p2 = cache.get_or_parse(path, ce, nullptr);
    std::printf("cache2 %d same=%d ne=%zu pc=%d\n", (int)(p2 != nullptr), (int)(p1 == p2), ce.size(), cache.parse_count());
    const input::ParsedSeriesFile* p3 = cache.get_or_parse(path + "/x", ce, &cs);
    std::printf("cache3 %d st=%d ne=%zu pc=%d\n", (int)(p3 != nullptr), (int)cs, ce.size(), cache.parse_count());
}

// ---------------------------------------------------------------- in-memory tables and probes
static void inmem() {
    std::printf("== INMEM\n");
    struct Case { const char* id; TableType ty; std::vector<double> x, y; int nc; };
    double nan = std::nan(""), inf = INFINITY;
    std::vector<Case> cases = {
        {"E", TableType::TIMESERIES, {}, {}, 1},
        {"DUP", TableType::TIMESERIES, {1, 2, 2, 3}, {1, 2, 3, 4}, 1},
        {"SEQ", TableType::TIMESERIES, {1, 3, 2}, {1, 2, 3}, 1},
        {"SEQ2", TableType::TIMESERIES, {45000.75, 45000.25}, {1, 2}, 1},
        {"CSEQ", TableType::CURVE_RATING, {0, 1, 0.5}, {1, 2, 3}, 1},
        {"NAN", TableType::TIMESERIES, {1, nan, 3}, {1, 2, 3}, 1},
        {"INF", TableType::TIMESERIES, {1, 2, inf}, {1, 2, 3}, 1},
        {"MIS", TableType::TIMESERIES, {1, 2}, {1, 2, 3}, 2},
        {"ONE", TableType::TIMESERIES, {40000.0}, {7, 8}, 2},
        {"OK2", TableType::TIMESERIES, {40000.0, 40000.5, 40000.5, 40001.25}, {1, 10, 2, 20, 3, 30, 4, 40}, 2},
    };
    for (auto& c : cases) {
        Table t; t.id = c.id; t.type = c.ty; t.x = c.x; t.y = c.y; t.num_cols = c.nc;
        std::printf("case %s\n", c.id);
        TableValidation v = validate_table(t); val(v, t);
        if ((size_t)c.nc * c.x.size() != c.y.size()) continue;
        double qs[] = {39999, 40000, 40000.25, 40000.5, 40000.75, 40001.25, 40002, nan, -inf, inf, 2.5, 0.5};
        for (int col = -1; col <= c.nc; ++col)
            for (double q : qs) {
                double a = table_lookup_column(t, col, q); hx(a); std::printf("%d %d ", t.cursor.index, t.cursor.direction);
                double b = table_step_column(t, col, q); hx(b); std::printf("%d %d\n", t.cursor.index, t.cursor.direction);
            }
    }
}
static void probes() {
    std::printf("== PROBES\n");
    const char* dts[] = {"2020-01-02", "2020-01-02 03:04", "2020-01-02T03:04:05", "2020-01-02  T 3:4", "2020-01-02x03:04",
                         "2020-1-2 25:61:61", "01/02/2020", "01/02/2020 3", "01/02/2020 03:04", "01/02/2020 03:04:05 PM",
                         "01/02/2020 12:00 am", "01/02/2020 12:00 PM", "01/02/2020 13:00 PM", "01/02/2020 00:30 AM",
                         "1/2/2020 7 pm", "13/45/2020", "2/29/2021 1:00", "-5-3-2", "", "   ", "x", "01/02", "2020-01",
                         "99999999999/1/2020", "1/2/2020 99999:00", " \t01/02/2020 04:05", "01/02/2020 04:05AM", "01/02/2020 AM",
                         "2020-01-02 03:04 PM", "01/02/2020 3:4:5 XM"};
    for (const char* d : dts) {
        double out = -7.0; bool r = input::parse_series_datetime(d, out);
        bs(d); std::printf("%d ", (int)r); hx(out); nl();
    }
    const char* toks[] = {"rain.csv", "rain.csv:G1", "C:\\dir\\f.csv", "C:\\dir\\f.csv:G2", "/data/Rich:2024/f.csv",
                          "/data/Rich:2024/f.csv:col", "f.csv:", ":col", "1:x", "a:b:c", "dir:sub/f", "x:y\\z", ""};
    for (const char* tk : toks) {
        std::string p = "old", c = "old"; bool r = input::split_series_file_token(tk, p, c);
        bs(tk); std::printf("%d ", (int)r); bs(p); bs(c); nl();
    }
}

// ---------------------------------------------------------------- run
static int run(const char* list) {
    std::FILE* f = std::fopen(list, "rb");
    if (!f) return 1;
    std::string all; char buf[65536]; size_t k;
    while ((k = std::fread(buf, 1, sizeof buf, f)) > 0) all.append(buf, k);
    std::fclose(f);
    std::vector<std::string> paths;
    size_t p = 0;
    while (p < all.size()) {
        size_t e = all.find('\n', p);
        if (e == std::string::npos) e = all.size();
        paths.push_back(all.substr(p, e - p));
        p = e + 1;
    }
    for (auto& path : paths) {
        std::printf("== FILE "); bs(path); nl();
        run_table(path, 128, true);
        run_table(path, 3, false);
        run_table(path, 0, false);
        run_mcsf(path);
    }
    inmem();
    probes();
    std::printf("total parses %ld\n", input::multicolumn_parse_count_total());
    return 0;
}

// ---------------------------------------------------------------- gen
static std::string D;
static std::vector<std::string> made;
static void wr(const std::string& name, const std::string& content) {
    std::string path = D + "/" + name;
    std::FILE* f = std::fopen(path.c_str(), "wb");
    std::fwrite(content.data(), 1, content.size(), f);
    std::fclose(f);
    made.push_back(path);
}
static uint64_t rs = 88172645463325252ULL;
static uint64_t rnd() { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
// minutes after 2020-01-01 00:00 -> calendar fields
static void cal(long minutes, int& y, int& mo, int& d, int& h, int& mi) {
    static const int dm[12] = {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}; // 2020
    long day = minutes / 1440; long rem = minutes % 1440;
    h = (int)(rem / 60); mi = (int)(rem % 60); y = 2020; mo = 1;
    while (day >= dm[mo - 1]) { day -= dm[mo - 1]; ++mo; if (mo > 12) { mo = 1; ++y; } }
    d = (int)day + 1;
}
static std::string us(long minutes) { int y, mo, d, h, mi; cal(minutes, y, mo, d, h, mi); char b[64]; std::snprintf(b, sizeof b, "%02d/%02d/%04d %02d:%02d", mo, d, y, h, mi); return b; }
static std::string iso(long minutes) { int y, mo, d, h, mi; cal(minutes, y, mo, d, h, mi); char b[64]; std::snprintf(b, sizeof b, "%04d-%02d-%02d %02d:%02d", y, mo, d, h, mi); return b; }
static std::string num(double v) { char b[64]; std::snprintf(b, sizeof b, "%.6g", v); return b; }

static std::vector<int> aval; static int agas, ansolid, acand;
static std::vector<std::pair<int, int>> acalls;
static bool acmp(int x, int y) { // McIlroy, "A Killer Adversary for Quicksort" (1999)
    acalls.emplace_back(x, y);
    if (aval[x] == agas && aval[y] == agas) { if (x == acand) aval[x] = ansolid++; else aval[y] = ansolid++; }
    if (aval[x] == agas) acand = x; else if (aval[y] == agas) acand = y;
    return aval[x] < aval[y];
}
static std::vector<int> tval;
static std::vector<std::pair<int, int>> rcalls;
static bool rcmp(int x, int y) { rcalls.emplace_back(x, y); return tval[x] < tval[y]; }

static int gen(const char* dir) {
    D = dir;
    std::vector<std::string> found;
    const char* roots[] = {"/root/work/corpus", "/root/work/swmm6src/upstream/tests", "/root/work/swmm6src/upstream/examples"};
    for (const char* r : roots) {
        std::error_code ec;
        for (auto it = fs::recursive_directory_iterator(r, ec); it != fs::recursive_directory_iterator(); it.increment(ec)) {
            if (ec || !it->is_regular_file()) continue;
            std::string name = it->path().filename().string(), ext = lower(it->path().extension().string());
            if (ext == ".dat" || ext == ".rain" || ext == ".csv" || ext == ".tsf" || ext == ".tsv" || (ext == ".txt" && name != "CMakeLists.txt"))
                found.push_back(it->path().string());
        }
    }
    std::sort(found.begin(), found.end());
    // Snapshot them: the shared corpus can change while the two sides run.
    {
        fs::create_directories(D + "/snap");
        std::string manifest;
        std::vector<std::string> snap;
        for (size_t i = 0; i < found.size(); ++i) {
            char pre[32]; std::snprintf(pre, sizeof pre, "%03zu_", i);
            std::string dst = D + "/snap/" + pre + fs::path(found[i]).filename().string();
            std::error_code ec;
            if (!fs::copy_file(found[i], dst, fs::copy_options::overwrite_existing, ec)) continue;
            snap.push_back(dst);
            manifest += dst + "\t" + found[i] + "\n";
        }
        std::FILE* mf = std::fopen((D + "/snap/MANIFEST").c_str(), "wb");
        std::fwrite(manifest.data(), 1, manifest.size(), mf);
        std::fclose(mf);
        found = snap;
    }

    std::string s;
    wr("adv_crlf.dat", "; SWMM .dat with CRLF\r\n\r\n01/01/2020 00:00 0.10 1.5\r\n01/01/2020 00:15 0.20 2.5\r\n# hash comment\r\n"
                       "01/01/2020 00:30:30 0.30 3.5\r\n01/01/2020 01:00 0.40\r\n01/01/2020 01:15 0.5 6.5 7.5\r\n");
    wr("adv_blank_comments.dat", "\n\n; c1\n   ; indented comment\n\t\t\n  01/02/2020 00:00 1.0\n# x\n01/02/2020 00:30 2.0\n\n"
                                 "01/02/2020 01:00 3.0 ; trailing\n   \n01/02/2020 01:30 4.0\n");
    wr("adv_out_of_order.dat", "01/03/2020 02:00 2\n01/03/2020 01:00 1\n01/03/2020 03:00 3\n01/03/2020 03:00 33\n01/03/2020 00:30 0.5\n"
                               "01/03/2020 04:00 4\n01/03/2020 03:30 3.5\n");
    wr("adv_time_only.dat", "0:00 1.0\n0:15 2.0\n1.5 3.0\n01/01/2020 02:00 4.0\n2:30 5.0\n01/01/2020 03:00 6.0\n");
    wr("adv_iso.csv", "DateTime,Rain,Flow\n2020-01-01 00:00,0.1,1.0\n2020-01-01T01:00,0.2,2.0\n2020-01-01   02:00:30,0.3,3\n"
                      "2020-01-02,0.4,4\n2020-1-3 5:7,0.5,5\n2020-01-03 06:00:00,0.6\n");
    wr("adv_us_missing.csv", "Time, A , B ,C\n01/02/2020 13:30:15,1.5,,3\n01/02/2020 14:00,,2\n01/02/2020 15:00\n"
                             "1/2/2020 16:00:00,4,5,6,7\nbad,1,2\n13/45/2020 00:00,1,1\n01/02/2020 17:00, 8 ,9 , 10\n");
    wr("adv_multi.tsv", "Date\tG1\tG2\tG3\n01/05/2020 00:00\t1\t2\t3\n01/05/2020 01:00\t\t5\t6\n01/05/2020 02:00\t7\n"
                        "2020-01-05T03:00:00\t8\t9\t10\n");
    wr("adv_ampm.tsf", "IDs:\tG1\tG2\nDate/Time\tRainfall\tRainfall\n\tin.\tin.\n01/01/2020 12:00:00 AM\t0.1\t0.2\n"
                       "01/01/2020 12:30:00 PM\t0.3\t0.4\n01/01/2020 01:00 PM\t0.5\t0.6\n01/01/2020 13:00 PM\t0.7\t0.8\n"
                       "01/01/2020 00:15 AM\t0.9\t1.0\n01/01/2020 3 PM\t1\t1\n01/01/2020 11:59:59 pm\t2\t2\n");
    wr("adv_ids_lower.tsf", "ids:\tX\n\n\n01/01/2020 01:00:00 AM\t5\n");
    wr("adv_quoted.csv", "\"Date\",\"Rain, total\",\"Flow\"\n\"2020-01-01 00:00\",\"1.5\",2\n\"2020-01-01 01:00\",\"2.5\",\"3\"\n"
                         "2020-01-01 02:00,3.5,4\n\"2020-01-01 03:00,4.5,5\n");
    s.clear();
    for (int i = 0; i < 10000; ++i) {
        if (i % 1000 == 500) s += "; comment row\n\n";
        s += us(5L * i) + " " + num((i % 97) * 0.01) + " " + num(i * 0.5) + "\n";
    }
    wr("adv_big.dat", s);
    s = "Date,Q1,Q2,Q3\n";
    for (int i = 0; i < 9000; ++i) s += iso(10L * i) + "," + num(i * 0.25) + "," + num((i % 13) * 1.5) + "," + num(9000 - i) + "\n";
    wr("adv_big.csv", s);
    for (int nrow : {4096, 8192, 8193}) {
        s.clear();
        for (int i = 0; i < nrow; ++i) s += us(60L * i) + " " + num(i * 0.001) + "\n";
        wr("adv_" + std::to_string(nrow) + ".dat", s);
    }
    s = "01/01/2020 00:00 1.0 2.0\n";
    int lens[] = {1022, 1023, 1024, 1025, 2046, 2047, 2048, 2049, 2700};
    for (int j = 0; j < 9; ++j) {
        std::string line = us(15L * (j + 1));
        while ((int)line.size() < lens[j]) line += " 1.234567";
        line.resize(lens[j]);
        s += line + "\n";
    }
    s += "01/01/2020 03:00 9.0 9.5\n";
    wr("adv_longline.dat", s);
    s = "Date";
    for (int j = 0; j < 300; ++j) { char b[16]; std::snprintf(b, sizeof b, ",Col%04d", j); s += b; }
    s += "\n2020-01-01 00:00,1,2,3\n2020-01-01 01:00,4,5,6\n";
    wr("adv_longheader.csv", s);
    wr("adv_bom.csv", "\xEF\xBB\xBF" "Date,R\n2020-01-01 00:00,1\n2020-01-01 01:00,2\n");
    wr("adv_bom.dat", "\xEF\xBB\xBF" "01/01/2020 00:00 1\n01/01/2020 01:00 2\n");
    wr("adv_bom_short", "\xEF\xBB");
    wr("adv_empty.dat", "");
    wr("adv_comments_only.dat", "; only comments\n# here\n\n");
    wr("adv_comment_csvlike.dat", "; a,b\n");
    {
        const char* cells[] = {"1", " 2.5 ", "+1e5", "nan(123)", "-nan", "0x1.8p3", "1e400", "-1e400", "1e-400", "4.9e-324",
                               "2.5e-324", "2e-324", "inf", "-Infinity", "infinit", "1.5abc", ".5", "5.", "1e", "1e+", "abc", "",
                               " ", "0x", "-0x.8p1", "nan(0x7ffffffffffff)", "nan(abc)", "nan(08)", "nan()", "nan(", "\"3.5\"",
                               "1d5", "0e999999", "00012", "-0", "+-1", "\v2", "\f3", "0x1p-1074", "0x1.fffffffffffffp1023",
                               "0x1.fffffffffffff8p1023", "1e-320", "123456789012345678901234567890", "0.1e-5",
                               "9007199254740993", "2.2250738585072011e-308", "-.5e-2x", "0X1P4", "0x1.000000000000081p0",
                               "0x10000000000000000001p-76", "-0x0.00000000000000000000001p0", "1.7976931348623158e308",
                               "1.7976931348623159e308", "NAN(7)", "-INF", "infinityx", "0.000000000000000000000000000001e30",
                               // hexadecimal rounding: ties to even, sticky bits, subnormal ties
                               "0x1.00000000000008p0", "0x1.00000000000018p0", "0x1.000000000000080000000001p0",
                               "0x1.00000000000007ffp0", "0x1.0000000000000800000000000000000001p0", "0x0.0000000000001p-1022",
                               "0x0.00000000000008p-1022", "0x0.00000000000018p-1022", "0x0.00000000000008000001p-1022",
                               "0x1.fffffffffffff7p1023", "-0x1.00000000000008p-1074", "0x3.fffffffffffffcp-1076"};
        s = "Date,V,W\n";
        int i = 0;
        for (const char* c : cells) { s += iso(i) + "," + c + "," + c + "\n"; ++i; }
        s += iso(i++) + ",0x0." + std::string(500, '0') + "1p2004,0x" + std::string(400, '0') + "1" + std::string(300, '0') + "p-1200\n";
        wr("adv_strtod.csv", s);
    }
    // NUL inside the date/time tokens; \f and \v between multi-column values; blank "\r\r" lines
#define LIT(x) std::string(x, sizeof(x) - 1)
    wr("adv_nul_date.dat", LIT("01/01/20\00020 00:00 1\n01/01/2020 00:\00030 2\n01/01/2020 01:00 3\n\0\n"));
    wr("adv_nul_date.csv", LIT("Date,A\n2020-01-01 0\0000:00,1\n2020-01-0\0001 01:00,2\n2020-01-01 02:00,3\n"));
    // %d is strtol (LONG_MIN/LONG_MAX on overflow) truncated to int: 4294969316 -> 2020, -1e20 -> 0
    wr("adv_bigint.dat", "1/1/4294969316 00:00 1\n1/1/2020 -99999999999999999999:30 2\n1/1/2020 3:-99999999999999999999 3\n"
                         "4294967297/2/2020 4:00 4\n1/1/2020 99999999999999999999:00 5\n");
    wr("adv_bigint.csv", "Date,A\n2020-01-01 -99999999999999999999:30,1\n4294969316-01-02 01:00,2\n"
                         "1/3/4294969316 -99999999999999999999:5,3\n2020-01-04T05:-99999999999999999999,4\n");
    wr("adv_ff_multi.dat", "01/01/2020 00:00 1 2 3\n01/01/2020 01:00 4\f5 6\n01/01/2020 02:00 7\v8 9\n01/01/2020 03:00 10 11 12\n");
    wr("adv_crcr.csv", "Date,A\r\r\n2020-01-01 00:00,1\r\n\r\r\n2020-01-01 01:00,2\r\n \r\n\r\n");
    s.clear(); // > 4096 rows, some longer than 1023 bytes: the count index and the row skips disagree
    for (int i = 0; i < 5000; ++i) {
        std::string line = us(30L * i) + " " + num(i * 0.1) + " " + num(i * 2.0);
        if (i % 100 == 7) line.resize(1500, ' ');
        s += line + "\n";
    }
    wr("adv_long_big.dat", s);
    wr("adv_noeol.dat", "01/01/2020 00:00 1\n01/01/2020 01:00 2");
    wr("adv_noeol.csv", "Date,A\n2020-01-01 00:00,1\n2020-01-01 01:00,2");
    wr("adv_dates.dat", "01/01/2020 00:00:30 1.0\n2/29/2021 00:00 1\n13/45/2020 00:00 1\n99999999999/1/2020 00:00 1\n"
                        "-1/5/2020 00:00 1\n01/01/2020 -1:30 1\n01/01/2020 99999:00 1\n01/01/2020 25:61:61 1\n1/1/20 0:0 1\n"
                        "01/01/2020 00 1\n01/01/2020x 00:00 1\n01-01-2020 00:00 1\n01/01/2020\t00:00\t2.5\n"
                        "01/01/2020 00:00 1.0\f2.0\n01/01/2020 00:00 1e400 2\n01/01/2020 00:00 inf nan -nan(12)\n"
                        "01/01/2020000000000000000000000000000000 00:00 1\n01/01/2020 00:00:00:00 3\n+1/+2/+2020 +3:+4 5\n"
                        "01/01/2020 00:00 1 x 2\n01/01/2020 00:00 -.5 .5 5. 1e 1e+2\n");
    {
        s = "Date,V,W\n";
        std::vector<int> idx;
        for (int i = 0; i < 300; ++i) { int k = (int)(rnd() % 40); s += iso(15L * k) + "," + std::to_string(i) + "," + num(i * 0.5) + "\n"; }
        wr("adv_shuffled.csv", s);
        s = "Date\tV\n";
        for (int i = 0; i < 5000; ++i) { int k = (int)(rnd() % 600); s += us(15L * k) + "\t" + std::to_string(i) + "\n"; }
        wr("adv_shuffled_big.tsv", s);
        int n = 2000;
        aval.assign(n, n); agas = n; ansolid = 0; acand = 0;
        std::vector<int> ptr(n); std::iota(ptr.begin(), ptr.end(), 0);
        std::sort(ptr.begin(), ptr.end(), acmp);
        for (int i = 0; i < n; ++i) if (aval[i] == agas) aval[i] = ansolid++;
        s = "Date,V\n";
        for (int i = 0; i < n; ++i) s += iso(aval[i]) + "," + std::to_string(i) + "\n";
        wr("adv_killer.csv", s);
        // Same comparison path with ties: value-adjacent elements merge when no call ever
        // asked "earlier < later" of them (those calls are the only ones a tie would flip).
        std::vector<char> asked((size_t)n * n, 0);
        for (auto& c : acalls) asked[(size_t)c.first * n + c.second] = 1;
        std::vector<int> by_val(n);
        for (int i = 0; i < n; ++i) by_val[aval[i]] = i;
        tval.assign(n, 0);
        std::vector<int> group;
        int g = 0;
        for (int v = 0; v < n; ++v) {
            int e = by_val[v];
            bool join = !group.empty();
            for (int m : group) if (asked[(size_t)m * n + e]) { join = false; break; }
            if (!join && !group.empty()) { ++g; group.clear(); }
            group.push_back(e); tval[e] = g;
        }
        std::iota(ptr.begin(), ptr.end(), 0);
        std::sort(ptr.begin(), ptr.end(), rcmp);
        if (rcalls != acalls) std::fprintf(stderr, "adv_killer_ties: replay path differs\n");
        std::fprintf(stderr, "adv_killer_ties: %d distinct of %d, %zu comparisons\n", g + 1, n, rcalls.size());
        s = "Date,V\n";
        for (int i = 0; i < n; ++i) s += iso(tval[i]) + "," + std::to_string(i) + "\n";
        wr("adv_killer_ties.csv", s);
    }
    wr("adv_numeric_header.csv", "Date,1,2\n2020-01-01 00:00,1,2\n2020-01-01 01:00,3,4\n");
    wr("adv_numeric_header_us.csv", "Time,1,2\n01/01/2020 00:00, 1, 2\n01/01/2020 01:00,3,4\n01/01/2020 02:00 5 6\n");
    wr("adv_tabs.dat", "01/01/2020\t00:00\t1.5\t2.5\n01/01/2020\t01:00\t3.5\t4.5\n");
    wr("adv_nul.dat", LIT("01/01/2020 00:00 1.0\0 junk 7\n01/01/2020 01:00 2.0 3.0\n\0\n01/01/2020 02:00 4\n"));
    wr("adv_nul.csv", LIT("Date,A,B\n2020-01-01 00:00,1\0,5,6\n2020-01-01 01:00,2,3\n"));
    wr("adv_leading_ws_comment.csv", "Date,A\n  ; comment\n2020-01-01 00:00,1\n\t# c\n   \n2020-01-01 01:00,2\n");
    wr("adv_header_only.csv", "Date,A\n");
    wr("adv_header_onecol.csv", "Date\n2020-01-01 00:00\n2020-01-01 01:00\n");
    wr("adv_mixed_eol.csv", "Date,A\r\n2020-01-01 00:00,1\n2020-01-01 01:00,2\r2020-01-01 02:00,3\r\n2020-01-01 03:00,4\r\r\n");
    wr("adv_dup_headers.csv", "Date,A,a,A\n2020-01-01 00:00,1,2,3\n");
    wr("adv_trailing_delim.csv", "Date,A,B,\n2020-01-01 00:00,1,2,\n2020-01-01 01:00,3,4,\n");
    wr("adv_ws_iso.tsv", "When\tX\n2020-01-01T00:00\t1\n2020-01-01   01:00\t2\n2020-01-01 T 02:00\t3\n");
    wr("adv_tab_dat_sniff.txt", "01/01/2020 00:00\t1\n01/01/2020 01:00\t2\n");
    wr("adv_space_header.txt", "Date Rain\n01/01/2020 00:00 1\n");
    wr("adv_semicolon_first.csv", ";Date,A\nDate,A\n2020-01-01 00:00,1\n");

    for (auto& f : found) std::printf("%s\n", f.c_str());
    for (auto& f : made) std::printf("%s\n", f.c_str());
    std::printf("%s/no_such_file.csv\n", dir);
    std::printf("%s\n", dir);
    std::printf("\n");
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 3 && std::strcmp(argv[1], "gen") == 0) return gen(argv[2]);
    if (argc == 3 && std::strcmp(argv[1], "run") == 0) return run(argv[2]);
    std::fprintf(stderr, "usage: tabfile gen DIR | run LIST\n");
    return 2;
}
