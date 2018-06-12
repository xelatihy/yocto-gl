//
// # Yocto/Utils: Tiny C++ Library of String, File and CLI utilities
//
// Yocto/Utils is a collection of generic utility functions mostly for building
// command line applications.
//
// 1. Path-like path operations: `path_dirname()`, `path_extension()`,
//    `path_basename()`, `path_filename()`, `replace_path_extension()`,
//    `prepend_path_extension()`, `split_path()`
// 2. Python-like format strings (only support for position arguments and no
//    formatting commands): `format()`, `print()`
// 3. simple logger with support for console and file streams:
//     1. create a `logger`
//     2. add more streams with `addconsole_stream()` or `add_file_stream()`
//     3. write log messages with `log_msg()` and its variants
//     4. you can also use a global default logger with the free functions
//        `log_XXX()`
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
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
//

#ifndef _YGL_UTILS_H_
#define _YGL_UTILS_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_math.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <ctime>
#include <map>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// STRING, PATH AND FILE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Get directory name (including '/').
inline std::string path_dirname(const std::string& filename) {
    auto pos = filename.rfind('/');
    if (pos == std::string::npos) pos = filename.rfind('\\');
    if (pos == std::string::npos) return "";
    return filename.substr(0, pos + 1);
}

// Get extension (including '.').
inline std::string path_extension(const std::string& filename) {
    auto pos = filename.rfind('.');
    if (pos == std::string::npos) return "";
    return filename.substr(pos);
}

// Get file basename.
inline std::string path_basename(const std::string& filename) {
    auto dirname = path_dirname(filename);
    auto extension = path_extension(filename);
    return filename.substr(
        dirname.size(), filename.size() - dirname.size() - extension.size());
}

// Get filename without directory (equiv to get_basename() +
// get_extension()).
inline std::string path_filename(const std::string& filename) {
    return path_basename(filename) + path_extension(filename);
}

// Replace extension.
inline std::string replace_path_extension(
    const std::string& filename, const std::string& ext) {
    return path_dirname(filename) + path_basename(filename) + ext;
}

// Prepend a string to the extension.
inline std::string prepend_path_extension(
    const std::string& filename, const std::string& prep) {
    return path_dirname(filename) + path_basename(filename) + prep +
           path_extension(filename);
}

// Really-minimal Python like string format. The implementation is not fast
// nor memory efficient. But it is good enough for some needs.
inline std::string format(
    const std::string& fmt, const std::vector<std::string>& args) {
    auto open = false;
    auto cur = 0;
    auto str = std::string();
    for (auto c : fmt) {
        if (c == '{') {
            str += args[cur++];
            open = true;
        } else if (c == '}') {
            if (!open) throw std::runtime_error("bad format");
            open = false;
        } else {
            str += c;
        }
    }
    return str;
}

// format value
inline std::string format_value(const std::string& val) { return val; }
inline std::string format_value(const char* val) { return val; }
inline std::string format_value(char* val) { return val; }
inline std::string format_value(const int& val) { return std::to_string(val); }
inline std::string format_value(const uint64_t& val) {
    return std::to_string(val);
}
inline std::string format_value(vec2i val) {
    return std::to_string(val.x) + " " + std::to_string(val.y);
}
inline std::string format_value(vec3i val) {
    return std::to_string(val.x) + " " + std::to_string(val.y) + " " +
           std::to_string(val.z);
}
inline std::string format_value(vec4i val) {
    return std::to_string(val.x) + " " + std::to_string(val.y) + " " +
           std::to_string(val.z) + " " + std::to_string(val.w);
}
inline std::string format_value(const float& val) {
    return std::to_string(val);
}
inline std::string format_value(vec2f val) {
    return std::to_string(val.x) + " " + std::to_string(val.y);
}
inline std::string format_value(vec3f val) {
    return std::to_string(val.x) + " " + std::to_string(val.y) + " " +
           std::to_string(val.z);
}
inline std::string format_value(vec4f val) {
    return std::to_string(val.x) + " " + std::to_string(val.y) + " " +
           std::to_string(val.z) + " " + std::to_string(val.w);
}
inline std::string format_value(const bool& val) {
    return (val) ? "true" : "false";
}

// Implementation of the function below.
inline void _format_one(std::vector<std::string>& vals) {}
template <typename Arg, typename... Args>
inline void _format_one(
    std::vector<std::string>& vals, const Arg& arg, const Args&... args) {
    vals.push_back(format_value(arg));
    _format_one(vals, args...);
}

// Really-minimal Python like string format. Internally uses streams for
// generality and supports for now only the '{}' operator. The implementation
// is not fast nor memory efficient. But it is good enough for some needs.
template <typename... Args>
inline std::string format(const std::string& fmt, const Args&... args) {
    auto vals = std::vector<std::string>();
    _format_one(vals, args...);
    return format(fmt, vals);
}

// Wrapper for `format()` that prints to stdout with/without ending newline.
template <typename... Args>
inline void print(const std::string& fmt, const Args&... args) {
    printf("%s", format(fmt, args...).c_str());
}
template <typename... Args>
inline void println(const std::string& fmt, const Args&... args) {
    printf("%s\n", format(fmt, args...).c_str());
}

// Log verbosity and output stream
inline bool& log_verbose() {
    static bool verbose = true;
    return verbose;
}

// Log callback function
using log_func = void (*)(const char* time, const char* tag, const char* msg);
inline log_func& log_callback() {
    static log_func func = [](const char* time, const char* tag,
                               const char* msg) {
        printf("%s %s %s\n", time, tag, msg);
    };
    return func;
}

// Implementation for logging functions.
inline void _log_msg(const std::string& msg, const char* tag) {
    using namespace std::chrono;
    auto time = system_clock::to_time_t(system_clock::now());
    auto tm = std::localtime(&time);
    char tmstr[64];
    strftime(tmstr, 64, "%T", tm);
    log_callback()(tmstr, tag, msg.c_str());
}

// Wrapper for `format()` that logs to a stream (default to stdout).
template <typename... Args>
inline void log_info(const std::string& msg, const Args&... args) {
    if (!log_verbose()) return;
    _log_msg(format(msg, args...), "INFO ");
}
template <typename... Args>
inline void log_warning(const std::string& msg, const Args&... args) {
    if (!log_verbose()) return;
    _log_msg(format(msg, args...), "WARN ");
}
template <typename... Args>
inline void log_error(const std::string& msg, const Args&... args) {
    _log_msg(format(msg, args...), "ERROR");
}
template <typename... Args>
inline void log_fatal(const std::string& msg, const Args&... args) {
    _log_msg(format(msg, args...), "FATAL");
    exit(1);
}

// Time-based logging
inline std::vector<std::pair<std::string, int64_t>>& _log_timed_stack() {
    static auto stack = std::vector<std::pair<std::string, int64_t>>();
    return stack;
}
template <typename... Args>
inline void log_info_begin(const std::string& msg, const Args&... args) {
    if (!log_verbose()) return;
    _log_msg(format("begin " + msg, args...), "INFO ");
    auto start =
        std::chrono::high_resolution_clock::now().time_since_epoch().count();
    _log_timed_stack().push_back({format(msg, args...), start});
}
inline void log_info_end() {
    if (!log_verbose()) return;
    auto start = _log_timed_stack().back().second;
    auto stop =
        std::chrono::high_resolution_clock::now().time_since_epoch().count();
    auto msg = _log_timed_stack().back().first;
    _log_timed_stack().pop_back();
    auto elapsed = (stop - start) / 1000000;
    auto hours = (int)(elapsed / 3600000);
    elapsed %= 3600000;
    auto mins = (int)(elapsed / 60000);
    elapsed %= 60000;
    auto secs = (int)(elapsed / 1000);
    auto msecs = (int)(elapsed % 1000);
    char buf[256];
    sprintf(buf, "%02d:%02d:%02d.%03d", hours, mins, secs, msecs);
    _log_msg("done " + msg + " [" + std::string(buf) + "]", "INFO ");
}

}  // namespace ygl

#endif
