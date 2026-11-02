// Differential harness: OpenSWMM 6's own Tokenizer (compiled from upstream source) over
// every line of the given files. One output line per input line: tokenize() tokens, then
// tokenize_views(), then is_numeric/is_boolean/parse_boolean per token.
#include "Tokenizer.hpp"
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
using openswmm::input::Tokenizer;
int main(int argc, char** argv) {
    std::string out;
    for (int a = 1; a < argc; ++a) {
        std::ifstream in(argv[a], std::ios::binary);
        std::string line;
        while (std::getline(in, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            for (auto& t : Tokenizer::tokenize(line)) {
                out += t; out += '\x1f';
                out += Tokenizer::is_numeric(t) ? 'N' : '-';
                out += Tokenizer::is_boolean(t) ? (Tokenizer::parse_boolean(t) ? 'T' : 'F') : '-';
                out += '\x1e';
            }
            out += '|';
            for (auto v : Tokenizer::tokenize_views(line)) { out.append(v.data(), v.size()); out += '\x1e'; }
            out += '|';
            out.append(Tokenizer::trim(Tokenizer::strip_comment(line)));
            out += '\n';
        }
    }
    fwrite(out.data(), 1, out.size(), stdout);
}
