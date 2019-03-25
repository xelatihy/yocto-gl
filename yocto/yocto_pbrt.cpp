//
// Implementation for Yocto/Pbrt loader.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

#include "yocto_pbrt.h"

#include <memory>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

using std::unique_ptr;

}

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF FAST PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// normalize obj line for simpler parsing
inline void normalize_pbrt_line(char* str, char comment_char = '#') {
    auto has_content = false;
    auto start       = str;
    while (*str) {
        if (*str == comment_char) {
            *str = 0;
            break;
        } else if (*str == ' ' || *str == '\t' || *str == '\r' ||
                   *str == '\n') {
            *str = ' ';
        } else {
            has_content = true;
        }
        str++;
    }
    if (!has_content) *start = 0;
}

// Parse values from a string
inline void parse_value(char*& str, int& value) {
    char* end = nullptr;
    value     = (int)strtol(str, &end, 10);
    if (str == end) throw pbrtio_error("cannot parse value");
    str = end;
}
inline void parse_value(char*& str, bool& value) {
    auto valuei = 0;
    parse_value(str, valuei);
    value = (bool)valuei;
}
inline void parse_value(char*& str, float& value) {
    char* end = nullptr;
    value     = strtof(str, &end);
    if (str == end) throw pbrtio_error("cannot parse value");
    str = end;
}
inline void parse_value(char*& str, string& value, bool ok_if_empty = false) {
    value = "";
    while (*str == ' ') str++;
    if (!*str && !ok_if_empty) {
        throw pbrtio_error("cannot parse value");
    }
    while (*str && *str != ' ') {
        value += *str;
        str++;
    }
}
template <typename T, int N>
inline void parse_value(char*& str, vec<T, N>& value) {
    for (auto i = 0; i < N; i++) parse_value(str, value[i]);
}
template <typename T, int N>
inline void parse_value(char*& str, frame<T, N>& value) {
    for (auto i = 0; i < N + 1; i++) parse_value(str, value[i]);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PBRT CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

// Load pbrt scene
void load_pbrt(const string& filename, const pbrt_callbacks& cb,
    const load_pbrt_options& options) {
    // open file
    auto fs = fopen(filename.c_str(), "rt");
    if (!fs) throw pbrtio_error("cannot load obj " + filename);
    auto fs_guard = unique_ptr<FILE, void (*)(FILE*)>{
        fs, [](FILE* f) { fclose(f); }};

    // read the file line by line
    char buffer[4096];
    while (fgets(buffer, sizeof(buffer), fs)) {
        // line
        auto line = buffer;
        normalize_pbrt_line(line);
        if (!*line) continue;
    }
}

}  // namespace yocto
