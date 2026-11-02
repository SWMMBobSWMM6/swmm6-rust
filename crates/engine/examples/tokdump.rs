// The Rust side of gate/diff/tokdump.cpp: the same output from the translated tokenizer.
use std::io::Write;
use swmm6_rust::input::tokenizer as tk;
fn main() {
    let mut out: Vec<u8> = Vec::new();
    for path in std::env::args().skip(1) {
        let data = std::fs::read(&path).expect("read");
        // std::getline: split on '\n'; a final line without '\n' still counts, an empty tail does not
        let mut lines: Vec<&[u8]> = data.split(|&c| c == b'\n').collect();
        if data.last() == Some(&b'\n') || data.is_empty() { lines.pop(); }
        for mut line in lines {
            if line.last() == Some(&b'\r') { line = &line[..line.len() - 1]; }
            for t in tk::tokenize(line) {
                out.extend_from_slice(&t); out.push(0x1f);
                out.push(if tk::is_numeric(&t) { b'N' } else { b'-' });
                out.push(if tk::is_boolean(&t) { if tk::parse_boolean(&t) { b'T' } else { b'F' } } else { b'-' });
                out.push(0x1e);
            }
            out.push(b'|');
            for v in tk::tokenize_views(line) { out.extend_from_slice(v); out.push(0x1e); }
            out.push(b'|');
            out.extend_from_slice(tk::trim(tk::strip_comment(line)));
            out.push(b'\n');
        }
    }
    std::io::stdout().write_all(&out).unwrap();
}
