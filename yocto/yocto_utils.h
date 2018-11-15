//
// # Yocto/Utils: Tiny collection of utilities to support Yocto/GL
//
//
// Yocto/Utils is a collection of utilities used in writing other Yocto/GL
// libraries and example applications. We support printing and parsing builting
// and Yocto/Math values, parsing command line arguments, simple path
// manipulation, file lading/saving and basic concurrency utilities.
//
//
// ## Printing and parsing values
//
// Use `format()` to format a string using `{}` as placeholder and `print()`
// to print it. Use `parse()` to parse a value from a string.
//
//
// ## Command-Line Parsing
//
// We provide a simple, immediate-mode, command-line parser. The parser
// works in an immediate-mode manner since it reads each value as you call each
// function, rather than building a data structure and parsing offline. We
// support option and position arguments, automatic help generation, and
// error checking.
//
// 1. initialize the parser with `make_cmdline_parser(argc, argv, help)`
// 2. read a value with `value = parse_argument(parser, name, default, help)`
//    - is name starts with '--' or '-' then it is an option
//    - otherwise it is a positional arguments
//    - options and arguments may be intermixed
//    - the type of each option is determined by the default value `default`
//    - the value is parsed on the stop
// 3. finished parsing with `check_cmdline(parser)`
//    - if an error occurred, the parser will exit and print a usage message
//
//
// ## Path manipulation
//
// We define a few path manipulation utilities to split and join path components.
//
//
// ## File IO
//
// 1. load and save text files with `load_text()` and `save_text()`
// 2. load and save binary files with `load_binary()` and `save_binary()`
//
//
// ## Concurrency utilities
//
// C++ has very basic supprt for concurrency and most of it is still platform
// dependent. We provide here very basic support for concurrency utlities
// built on top of C++ low-level threading and synchronization.
//
// 1. use `concurrent_queue()` for communicationing values between threads
// 2. use `parallel_for()` for basic parallel for loops
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

#ifndef _YOCTO_UTILS_H_
#define _YOCTO_UTILS_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_math.h"

#include <cctype>
#include <chrono>
#include <cstdio>
#include <deque>
#include <mutex>
#include <string>
#include <string_view>
#include <thread>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

using std::deque;
using std::lock_guard;
using std::mutex;
using std::string_view;
using std::thread;
using namespace std::chrono_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// PRINT/PARSE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Formats a string `fmt` with values taken from `args`. Uses `{}` as
// placeholder.
template <typename... Args>
inline string format(const string& fmt, const Args&... args);

// Converts to string.
template <typename T>
inline string to_string(const T& value);

// Prints a formatted string to stdout or file.
template <typename... Args>
inline bool print(FILE* fs, const string& fmt, const Args&... args);
template <typename... Args>
inline bool print(const string& fmt, const Args&... args);

// Format duration string from nanoseconds
inline string format_duration(int64_t duration);
// Format a large integer number in human readable form
inline string format_num(uint64_t num);

// Parse a list of space separated values.
template <typename... Args>
inline bool parse(const string& str, Args&... args);

// get time in nanoseconds - useful only to compute difference of times
inline int64_t get_time();

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOGGING UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Log info/error/fatal/trace message
template <typename... Args>
inline void log_info(const string& fmt, const Args&... args);
template <typename... Args>
inline void log_error(const string& fmt, const Args&... args);
template <typename... Args>
inline void log_warning(const string& fmt, const Args&... args);
template <typename... Args>
inline void log_fatal(const string& fmt, const Args&... args);

// log levels
enum struct log_level {
    fatal   = 0,
    error   = 1,
    warning = 2,
    info    = 3,
    trace   = 4
};

// Setup logging
inline void set_log_level(log_level level);
inline void set_log_console(bool enabled);
inline void set_log_file(const string& filename, bool append = false);

// Log traces for timing and program debugging
struct log_scope;
template <typename... Args>
inline void log_trace(const string& fmt, const Args&... args);
template <typename... Args>
inline log_scope log_trace_begin(const string& fmt, const Args&... args);
template <typename... Args>
inline void log_trace_end(log_scope& scope);
template <typename... Args>
inline log_scope log_trace_scoped(const string& fmt, const Args&... args);

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMMEDIATE-MODE COMMAND LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// Command line parser data. All data should be considered private.
struct cmdline_parser {
    vector<string> args              = {};    // command line arguments
    string         help_command      = "";    // program name
    string         help_usage        = "";    // program help
    string         help_options      = "";    // options help
    string         help_arguments    = "";    // arguments help
    string         error             = "";    // current parse error
    bool           add_help_flag     = true;  // adding help flag
    bool           add_logging_flags = true;  // adding logging flags
};

// Initialize a command line parser.
inline cmdline_parser make_cmdline_parser(int argc, char** argv,
    const string& usage, const string& cmd = "", bool add_help_flag = true,
    bool add_logging_flags = true);
// check if any error occurred and exit appropriately
inline void check_cmdline(cmdline_parser& parser);

// Parse an int, float, string, vecXX and bool option or positional argument.
// Options's names starts with "--" or "-", otherwise they are arguments.
// vecXX options use space-separated values but all in one argument
// (use " or ' from the common line).
// Boolean flags are indicated with a pair of names "--name/--no-name", so
// that we have both options available. You can also use the parse flag function
// in which case only one name is used and the flag will flip the value passed.
template <typename T>
inline T    parse_argument(cmdline_parser& parser, const string& name, T def,
       const string& usage, bool req = false);
inline bool parse_argument_flag(
    cmdline_parser& parser, const string& name, bool def, const string& usage);
// Parse all arguments left on the command line.
template <typename T>
inline vector<T> parse_arguments(cmdline_parser& parser, const string& name,
    const vector<T>& def, const string& usage, bool req = false);
// Parse a labeled enum, with enum values that are successive integers.
template <typename T>
inline T parse_argument(cmdline_parser& parser, const string& name, T def,
    const string& usage, const vector<string>& labels, bool req = false);

// Parse an int, float, string, vecXX and bool option or positional argument.
// Options's names starts with "--" or "-", otherwise they are arguments.
// vecXX options use space-separated values but all in one argument
// (use " or ' from the common line). Booleans are flags.
// Boolean flags are indicated with a pair of names "--name/--no-name", so
// that we have both options available. You can also use the parse flag function
// in which case only one name is used and the flag will flip the value passed.
template <typename T>
inline bool parse_argument_ref(cmdline_parser& parser, const string& name,
    T& val, const string& usage, bool req = false);
inline bool parse_argument_flag(
    cmdline_parser& parser, const string& name, bool& val, const string& usage);
// Parse all arguments left on the command line.
template <typename T>
inline bool parse_arguments_ref(cmdline_parser& parser, const string& name,
    vector<T>& val, const string& usage, bool req = false);
// Parse a labeled enum, with enum values that are successive integers.
template <typename T>
inline bool parse_argument_ref(cmdline_parser& parser, const string& name,
    T& val, const string& usage, const vector<string>& labels, bool req = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// PATH UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Normalize path delimiters.
inline string normalize_path(const string& filename);
// Get directory name (not including '/').
inline string get_dirname(const string& filename);
// Get extension (not including '.').
inline string get_extension(const string& filename);
// Get filename without directory.
inline string get_filename(const string& filename);
// Replace extension.
inline string replace_extension(const string& filename, const string& ext);

// Check if a file can be opened for reading.
inline bool exists_file(const string& filename);

}  // namespace yocto

// -----------------------------------------------------------------------------
// FILE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a text file
inline bool load_text(const string& filename, string& str);
inline bool save_text(const string& filename, const string& str);

// Load/save a binary file
inline bool load_binary(const string& filename, vector<byte>& data);
inline bool save_binary(const string& filename, const vector<byte>& data);

}  // namespace yocto

// -----------------------------------------------------------------------------
// CONCURRENCY UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// a simple concurrent queue that locks at every call
template <typename T>
struct concurrent_queue {
    concurrent_queue();
    concurrent_queue(const concurrent_queue& other);
    concurrent_queue& operator=(const concurrent_queue& other);

    bool empty();
    void clear();
    void push(const T& value);
    bool try_pop(T& value);

   private:
    mutex    _mutex;
    deque<T> _queue;
};

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename Func>
inline void parallel_for(int begin, int end, const Func& func,
    atomic<bool>* cancel = nullptr, bool serial = false);
template <typename Func>
inline void parallel_for(int num, const Func& func,
    atomic<bool>* cancel = nullptr, bool serial = false) {
    parallel_for(0, num, func, cancel, serial);
}

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes a reference to a `T`.
template <typename T, typename Func>
inline void parallel_foreach(vector<T>& values, const Func& func,
    atomic<bool>* cancel = nullptr, bool serial = false) {
    parallel_for(0, (int)values.size(),
        [&func, &values](int idx) { func(values[idx]); }, cancel, serial);
}
template <typename T, typename Func>
inline void parallel_foreach(const vector<T>& values, const Func& func,
    atomic<bool>* cancel = nullptr, bool serial = false) {
    parallel_for(0, (int)values.size(),
        [&func, &values](int idx) { func(values[idx]); }, cancel, serial);
}

}  // namespace yocto

// ---------------------------------------------------------------------------//
//                                                                            //
//                             IMPLEMENTATION                                 //
//                                                                            //
// ---------------------------------------------------------------------------//

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF STRING/TIME UTILITIES FOR CLI APPLICATIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Prints basic types
inline bool print_value(string& str, const string& value) {
    str += value;
    return true;
}
inline bool print_value(string& str, const char* value) {
    str += value;
    return true;
}
inline bool print_value(string& str, int value) {
    str += std::to_string(value);
    return true;
}
inline bool print_value(string& str, float value) {
    str += std::to_string(value);
    return true;
}
inline bool print_value(string& str, double value) {
    str += std::to_string(value);
    return true;
}
template <typename T>
inline bool print_value(string& str, const T* value) {
    char buffer[512];
    sprintf(buffer, "%p", value);
    str += buffer;
    return true;
}

// Print compound types.
inline bool print_value(string& str, const vec2f& v) {
    print_value(str, v.x);
    print_value(str, " ");
    print_value(str, v.y);
    return true;
}
inline bool print_value(string& str, const vec3f& v) {
    print_value(str, v.x);
    print_value(str, " ");
    print_value(str, v.y);
    print_value(str, " ");
    print_value(str, v.z);
    return true;
}
inline bool print_value(string& str, const vec4f& v) {
    print_value(str, v.x);
    print_value(str, " ");
    print_value(str, v.y);
    print_value(str, " ");
    print_value(str, v.z);
    print_value(str, " ");
    print_value(str, v.w);
    return true;
}
inline bool print_value(string& str, const vec2i& v) {
    print_value(str, v.x);
    print_value(str, " ");
    print_value(str, v.y);
    return true;
}
inline bool print_value(string& str, const vec3i& v) {
    print_value(str, v.x);
    print_value(str, " ");
    print_value(str, v.y);
    print_value(str, " ");
    print_value(str, v.z);
    return true;
}
inline bool print_value(string& str, const vec4i& v) {
    print_value(str, v.x);
    print_value(str, " ");
    print_value(str, v.y);
    print_value(str, " ");
    print_value(str, v.z);
    print_value(str, " ");
    print_value(str, v.w);
    return true;
}
inline bool print_value(string& str, const mat2f& v) {
    print_value(str, v.x);
    print_value(str, " ");
    print_value(str, v.y);
    return true;
}
inline bool print_value(string& str, const mat3f& v) {
    print_value(str, v.x);
    print_value(str, " ");
    print_value(str, v.y);
    print_value(str, " ");
    print_value(str, v.z);
    return true;
}
inline bool print_value(string& str, const mat4f& v) {
    print_value(str, v.x);
    print_value(str, " ");
    print_value(str, v.y);
    print_value(str, " ");
    print_value(str, v.z);
    print_value(str, " ");
    print_value(str, v.w);
    return true;
}
inline bool print_value(string& str, const frame2f& v) {
    print_value(str, v.x);
    print_value(str, " ");
    print_value(str, v.y);
    print_value(str, " ");
    print_value(str, v.o);
    return true;
}
inline bool print_value(string& str, const frame3f& v) {
    print_value(str, v.x);
    print_value(str, " ");
    print_value(str, v.y);
    print_value(str, " ");
    print_value(str, v.z);
    print_value(str, " ");
    print_value(str, v.o);
    return true;
}
inline bool print_value(string& str, const bbox1f& v) {
    print_value(str, v.min);
    print_value(str, " ");
    print_value(str, v.max);
    return true;
}
inline bool print_value(string& str, const bbox2f& v) {
    print_value(str, v.min);
    print_value(str, " ");
    print_value(str, v.max);
    return true;
}
inline bool print_value(string& str, const bbox3f& v) {
    print_value(str, v.min);
    print_value(str, " ");
    print_value(str, v.max);
    return true;
}
inline bool print_value(string& str, const bbox4f& v) {
    print_value(str, v.min);
    print_value(str, " ");
    print_value(str, v.max);
    return true;
}
inline bool print_value(string& str, const bbox1i& v) {
    print_value(str, v.min);
    print_value(str, " ");
    print_value(str, v.max);
    return true;
}
inline bool print_value(string& str, const bbox2i& v) {
    print_value(str, v.min);
    print_value(str, " ");
    print_value(str, v.max);
    return true;
}
inline bool print_value(string& str, const bbox3i& v) {
    print_value(str, v.min);
    print_value(str, " ");
    print_value(str, v.max);
    return true;
}
inline bool print_value(string& str, const bbox4i& v) {
    print_value(str, v.min);
    print_value(str, " ");
    print_value(str, v.max);
    return true;
}
inline bool print_value(string& str, const ray2f& v) {
    print_value(str, v.o);
    print_value(str, " ");
    print_value(str, v.d);
    print_value(str, " ");
    print_value(str, v.tmin);
    print_value(str, " ");
    print_value(str, v.tmax);
    return true;
}

// Prints a string.
inline bool print_next(string& str, const string& fmt) {
    return print_value(str, fmt);
}
template <typename Arg, typename... Args>
inline bool print_next(
    string& str, const string& fmt, const Arg& arg, const Args&... args) {
    auto pos = fmt.find("{}");
    if (pos == string::npos) return print_value(str, fmt);
    if (!print_value(str, fmt.substr(0, pos))) return false;
    if (!print_value(str, arg)) return false;
    return print_next(str, fmt.substr(pos + 2), args...);
}

// Formats a string `fmt` with values taken from `args`. Uses `{}` as
// placeholder.
template <typename... Args>
inline string format(const string& fmt, const Args&... args) {
    auto str = string();
    print_next(str, fmt, args...);
    return str;
}

// Prints a string.
template <typename... Args>
inline bool print(FILE* fs, const string& fmt, const Args&... args) {
    auto str = format(fmt, args...);
    return fprintf(fs, "%s", str.c_str()) >= 0;
}
template <typename... Args>
inline bool print(const string& fmt, const Args&... args) {
    return print(stdout, fmt, args...);
}

// Converts to string.
template <typename T>
inline string to_string(const T& value) {
    auto str = string();
    print_value(str, value);
    return str;
}

// Prints basic types to string
inline bool parse_value(string_view& str, string& value) {
    auto pos = str.find_first_not_of(" \t\r\n");
    if(pos == string_view::npos) return false;
    str.remove_prefix(pos);
    pos = str.find_first_of(" \t\r\n");
    if(pos == string_view::npos) {
        value = str;
        str.remove_prefix(str.length());
    } else {
        value = str.substr(0, pos);
        str.remove_prefix(pos);
    }
    return true;
}
inline bool parse_value(string_view& str, int& value) {
    char* end = nullptr;
    value     = (int)strtol(data(str), &end, 10);
    if (data(str) == end) return false;
    str.remove_prefix(end - data(str));
    // auto n = 0;
    // if (sscanf(str.str, "%d%n", &value, &n) != 1) return false;
    // str.str += n;
    return true;
}
inline bool parse_value(string_view& str, float& value) {
    char* end = nullptr;
    value     = strtof(data(str), &end);
    if (data(str) == end) return false;
    str.remove_prefix(end - data(str));
    // auto n = 0;
    // if (sscanf(str.str, "%f%n", &value, &n) != 1) return false;
    // str.str += n;
    return true;
}
inline bool parse_value(string_view& str, double& value) {
    char* end = nullptr;
    value     = strtod(data(str), &end);
    if (data(str) == end) return false;
    str.remove_prefix(end - data(str));
    // auto n = 0;
    // if (sscanf(str.str, "%lf%n", &value, &n) != 1) return false;
    // str.str += n;
    return true;
}
inline bool parse_value(string_view& str, bool& value) {
    auto ivalue = 0;
    if (!parse_value(str, ivalue)) return false;
    value = (bool)ivalue;
    return true;
}

// Print compound types
template <typename T, size_t N>
inline bool parse_value(string_view& str, array<T, N>& value) {
    for (auto i = 0; i < N; i++) {
        if (!parse_value(str, value[i])) return false;
    }
    return true;
}
template <typename T>
inline bool parse_values(string_view& str, T* values, int N) {
    for (auto i = 0; i < N; i++) {
        if (!parse_value(str, values[i])) return false;
    }
    return true;
}

// Data acess
inline bool parse_value(string_view& str, vec2f& v) {
    if(!parse_value(str, v.x)) return false;
    if(!parse_value(str, v.y)) return false;
    return true;
}
inline bool parse_value(string_view& str, vec3f& v) {
    if(!parse_value(str, v.x)) return false;
    if(!parse_value(str, v.y)) return false;
    if(!parse_value(str, v.z)) return false;
    return true;
}
inline bool parse_value(string_view& str, vec4f& v) {
    if(!parse_value(str, v.x)) return false;
    if(!parse_value(str, v.y)) return false;
    if(!parse_value(str, v.z)) return false;
    if(!parse_value(str, v.w)) return false;
    return true;
}
inline bool parse_value(string_view& str, vec2i& v) {
    if(!parse_value(str, v.x)) return false;
    if(!parse_value(str, v.y)) return false;
    return true;
}
inline bool parse_value(string_view& str, vec3i& v) {
    if(!parse_value(str, v.x)) return false;
    if(!parse_value(str, v.y)) return false;
    if(!parse_value(str, v.z)) return false;
    return true;
}
inline bool parse_value(string_view& str, vec4i& v) {
    if(!parse_value(str, v.x)) return false;
    if(!parse_value(str, v.y)) return false;
    if(!parse_value(str, v.z)) return false;
    if(!parse_value(str, v.w)) return false;
    return true;
}
inline bool parse_value(string_view& str, mat2f& v) {
    if(!parse_value(str, v.x)) return false;
    if(!parse_value(str, v.y)) return false;
    return true;
}
inline bool parse_value(string_view& str, mat3f& v) {
    if(!parse_value(str, v.x)) return false;
    if(!parse_value(str, v.y)) return false;
    if(!parse_value(str, v.z)) return false;
    return true;
}
inline bool parse_value(string_view& str, mat4f& v) {
    if(!parse_value(str, v.x)) return false;
    if(!parse_value(str, v.y)) return false;
    if(!parse_value(str, v.z)) return false;
    if(!parse_value(str, v.w)) return false;
    return true;
}
inline bool parse_value(string_view& str, frame2f& v) {
    if(!parse_value(str, v.x)) return false;
    if(!parse_value(str, v.y)) return false;
    if(!parse_value(str, v.o)) return false;
    return true;
}
inline bool parse_value(string_view& str, frame3f& v) {
    if(!parse_value(str, v.x)) return false;
    if(!parse_value(str, v.y)) return false;
    if(!parse_value(str, v.z)) return false;
    if(!parse_value(str, v.o)) return false;
    return true;
}
inline bool parse_value(string_view& str, bbox1f& v) {
    if(!parse_value(str, v.min)) return false;
    if(!parse_value(str, v.max)) return false;
    return true;
}
inline bool parse_value(string_view& str, bbox2f& v) {
    if(!parse_value(str, v.min)) return false;
    if(!parse_value(str, v.max)) return false;
    return true;
}
inline bool parse_value(string_view& str, bbox3f& v) {
    if(!parse_value(str, v.min)) return false;
    if(!parse_value(str, v.max)) return false;
    return true;
}
inline bool parse_value(string_view& str, bbox4f& v) {
    if(!parse_value(str, v.min)) return false;
    if(!parse_value(str, v.max)) return false;
    return true;
}
inline bool parse_value(string_view& str, bbox1i& v) {
    if(!parse_value(str, v.min)) return false;
    if(!parse_value(str, v.max)) return false;
    return true;
}
inline bool parse_value(string_view& str, bbox2i& v) {
    if(!parse_value(str, v.min)) return false;
    if(!parse_value(str, v.max)) return false;
    return true;
}
inline bool parse_value(string_view& str, bbox3i& v) {
    if(!parse_value(str, v.min)) return false;
    if(!parse_value(str, v.max)) return false;
    return true;
}
inline bool parse_value(string_view& str, bbox4i& v) {
    if(!parse_value(str, v.min)) return false;
    if(!parse_value(str, v.max)) return false;
    return true;
}
inline bool parse_value(string_view& str, ray2f& v) {
    if(!parse_value(str, v.o)) return false;
    if(!parse_value(str, v.d)) return false;
    if(!parse_value(str, v.tmin)) return false;
    if(!parse_value(str, v.tmax)) return false;
    return true;
}
inline bool parse_value(string_view& str, ray3f& v) {
    if(!parse_value(str, v.o)) return false;
    if(!parse_value(str, v.d)) return false;
    if(!parse_value(str, v.tmin)) return false;
    if(!parse_value(str, v.tmax)) return false;
    return true;
}

// Prints a string.
inline bool parse_next(string_view& str) { return true; }
template <typename Arg, typename... Args>
inline bool parse_next(string_view& str, Arg& arg, Args&... args) {
    if (!parse_value(str, arg)) return false;
    return parse_next(str, args...);
}

// Returns trus if this is white space
inline bool is_whitespace(string_view str) {
    return str.find_first_not_of(" \t\r\n") == string_view::npos;
}

// Parse a list of space separated values.
template <typename... Args>
inline bool parse(const string& str, Args&... args) {
    auto view = string_view{str.c_str()};
    if (!parse_next(view, args...)) return false;
    return is_whitespace(view);
}

// get time in nanoseconds - useful only to compute difference of times
inline int64_t get_time() {
    return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF LOGGING UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Logging configutation
inline bool& _log_console() {
    static auto _log_console = true;
    return _log_console;
}
inline FILE*& _log_filestream() {
    static auto _log_filestream = (FILE*)nullptr;
    return _log_filestream;
}
inline log_level& _log_level() {
    static auto _log_level = log_level::info;
    return _log_level;
}
inline bool is_log_level_skipped(log_level level) {
    return level > _log_level();
}

// Logs a message
inline void log_message(log_level level, const char* msg) {
    static const char* labels[] = {"FATAL", "ERROR", "WARN ", "INFO ", "TRACE"};
    if (_log_console()) {
        printf("%s\n", msg);
        fflush(stdout);
    }
    if (_log_filestream()) {
        fprintf(_log_filestream(), "%s %s\n", labels[(int)level], msg);
        fflush(_log_filestream());
    }
}

// Log info/error/fatal/trace message
template <typename... Args>
inline void log_info(const string& fmt, const Args&... args) {
    if (is_log_level_skipped(log_level::info)) return;
    log_message(log_level::info, format(fmt, args...).c_str());
}
template <typename... Args>
inline void log_error(const string& fmt, const Args&... args) {
    if (is_log_level_skipped(log_level::error)) return;
    log_message(log_level::error, format(fmt, args...).c_str());
}
template <typename... Args>
inline void log_warning(const string& fmt, const Args&... args) {
    if (is_log_level_skipped(log_level::warning)) return;
    log_message(log_level::warning, format(fmt, args...).c_str());
}
template <typename... Args>
inline void log_fatal(const string& fmt, const Args&... args) {
    if (is_log_level_skipped(log_level::fatal)) return;
    log_message(log_level::fatal, format(fmt, args...).c_str());
    exit(1);
}

// Log traces for timing and program debugging
struct log_scope {
    string  message    = "";
    int64_t start_time = -1;
    bool    scoped     = false;
    ~log_scope();
};
template <typename... Args>
inline void log_trace(const string& fmt, const Args&... args) {
    if (is_log_level_skipped(log_level::trace)) return;
    log_message(log_level::trace, format(fmt, args...).c_str());
}
template <typename... Args>
inline log_scope log_trace_begin(const string& fmt, const Args&... args) {
    if (is_log_level_skipped(log_level::trace)) return {"", -1, false};
    auto message = format(fmt, args...);
    log_trace(message + " [started]");
    return {message, get_time(), false};
}
template <typename... Args>
inline void log_trace_end(log_scope& scope) {
    if (is_log_level_skipped(log_level::trace)) return;
    if (scope.start_time >= 0) {
        log_trace(scope.message + " [ended: " +
                  format_duration(get_time() - scope.start_time) + "]");
    } else {
        log_trace(scope.message + " [ended]");
    }
}
template <typename... Args>
inline log_scope log_trace_scoped(const string& fmt, const Args&... args) {
    if (is_log_level_skipped(log_level::trace)) return {"", -1, false};
    auto message = format(fmt, args...);
    log_trace(message + " [started]");
    return {message, get_time(), true};
}
inline log_scope::~log_scope() {
    if (scoped) log_trace_end(*this);
}

// Configure the logging
inline void set_log_level(log_level level) { _log_level() = level; }
inline void set_log_console(bool enabled) { _log_console() = enabled; }
inline void set_log_file(const string& filename, bool append) {
    if (_log_filestream()) {
        fclose(_log_filestream());
        _log_filestream() = nullptr;
    }
    if (empty(filename)) return;
    _log_filestream() = fopen(filename.c_str(), append ? "at" : "wt");
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF STRING FORMAT UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Format duration string from nanoseconds
inline string format_duration(int64_t duration) {
    auto elapsed = duration / 1000000;  // milliseconds
    auto hours   = (int)(elapsed / 3600000);
    elapsed %= 3600000;
    auto mins = (int)(elapsed / 60000);
    elapsed %= 60000;
    auto secs  = (int)(elapsed / 1000);
    auto msecs = (int)(elapsed % 1000);
    char buffer[256];
    sprintf(buffer, "%02d:%02d:%02d.%03d", hours, mins, secs, msecs);
    return buffer;
}
// Format a large integer number in human readable form
inline string format_num(uint64_t num) {
    auto rem = num % 1000;
    auto div = num / 1000;
    if (div > 0) return format_num(div) + "," + std::to_string(rem);
    return std::to_string(rem);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF COMMAND-LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// initialize a command line parser
inline cmdline_parser make_cmdline_parser(int argc, char** argv,
    const string& usage, const string& cmd, bool add_help_flag,
    bool add_logging_flags) {
    auto parser              = cmdline_parser{};
    parser.args              = {argv + 1, argv + argc};
    parser.help_command      = (empty(cmd)) ? argv[0] : cmd;
    parser.help_usage        = usage;
    parser.add_help_flag     = add_help_flag;
    parser.add_logging_flags = add_logging_flags;
    return parser;
}

// check if option or argument
inline bool is_optional_argument(const string& name) {
    return name.size() > 1 && name.front() == '-';
}

// check if flag
inline bool is_optional_flag(const string& name) {
    return name.size() > 1 && name.front() == '-' && name.find('/') != name.npos;
}

// get names from string
inline vector<string> get_option_names(const string& name_) {
    auto names = vector<string>();
    auto name  = name_;
    while (name.find(',') != name.npos) {
        names.push_back(name.substr(0, name.find(',')));
        name = name.substr(name.find(',') + 1);
    }
    names.push_back(name);
    return names;
}

// get names from string
inline vector<pair<string, string>> get_flag_names(const string& name) {
    auto names = vector<pair<string, string>>();
    for (auto& name : get_option_names(name)) {
        if (name.find('/')) {
            names.push_back({name.substr(0, name.find('/')),
                name.substr(name.find('/') + 1)});
        } else {
            names.push_back({name, ""});
        }
    }
    return names;
}

// get default string
template <typename T>
inline string get_option_default_string(const T& value) {
    return to_string(value);
}
inline string get_option_default_string(const bool& value) {
    return (value) ? "true"s : "false"s;
}
template <typename T>
inline string get_option_default_string(const vector<T>& values) {
    auto defs = string();
    for (auto& d : values) defs += " " + d;
    return defs;
}

// get option typename
template <typename T>
inline string get_option_typename() {
    return "value";
}
template <>
inline string get_option_typename<int>() {
    return "int";
}
template <>
inline string get_option_typename<bool>() {
    return "";
}
template <>
inline string get_option_typename<float>() {
    return "float";
}
template <>
inline string get_option_typename<string>() {
    return "string";
}

// add help
template <typename T>
inline string get_option_usage(const string& name, const string& usage,
    const T& def_, bool req, const vector<string>& choices) {
    auto def = ""s;
    if (!req) {
        def = get_option_default_string(def_);
        if (def != "") def = "[" + def + "]";
    }
    auto nametype = name;
    if (get_option_typename<T>() != "") {
        nametype += " <" + get_option_typename<T>() + ">";
    }
    char buffer[4096];
    sprintf(buffer, "  %-24s %s %s\n", nametype.c_str(), usage.c_str(),
        def.c_str());
    auto usagelines = string(buffer);
    if (!empty(choices)) {
        usagelines += "        accepted values:";
        for (auto& c : choices) usagelines += " " + c;
        usagelines += "\n";
    }
    return usagelines;
}

// print cmdline help
inline void print_cmdline_usage(const cmdline_parser& parser) {
    printf("%s: %s\n", parser.help_command.c_str(), parser.help_usage.c_str());
    printf("usage: %s %s %s\n\n", parser.help_command.c_str(),
        (empty(parser.help_options)) ? "" : "[options]",
        (empty(parser.help_arguments)) ? "" : "arguments");
    if (!empty(parser.help_options)) {
        printf("options:\n");
        printf("%s\n", parser.help_options.c_str());
    }
    if (!empty(parser.help_arguments)) {
        printf("arguments:\n");
        printf("%s\n", parser.help_arguments.c_str());
    }
}

// Parse a flag. Name should start with either "--" or "-".
inline bool parse_flag_argument(cmdline_parser& parser, const string& name,
    bool& value, const string& usage, bool req);

// check if any error occurred and exit appropriately
inline void check_cmdline(cmdline_parser& parser) {
    if (parser.add_help_flag) {
        auto help = false;
        if (parse_flag_argument(parser, "--help,-?", help, "print help", false)) {
            print_cmdline_usage(parser);
            exit(0);
        }
    }
    if (parser.add_logging_flags) {
        auto verbose = false, quiet = false;
        if (parse_flag_argument(
                parser, "--verbose,-v", verbose, "use verbose output", false)) {
            set_log_level(log_level::trace);
        }
        if (parse_flag_argument(
                parser, "--quiet,-q", quiet, "use quiet output", false)) {
            set_log_level(log_level::error);
        }
    }
    if (!empty(parser.args)) {
        auto found = false;
        for (auto& name : parser.args) {
            if (is_optional_argument(name)) {
                parser.error += "unknow option " + name + "\n";
                found = true;
                break;
            }
        }
        if (!found) parser.error += "unmatched arguments remaining\n";
    }
    if (!empty(parser.error)) {
        printf("error: %s\n", parser.error.c_str());
        print_cmdline_usage(parser);
        exit(1);
    }
}

// Parse an option string. Name should start with "--" or "-".
template <typename T>
inline bool parse_option_argument(cmdline_parser& parser, const string& name,
    T& value, const string& usage, bool req, const vector<string>& choices) {
    parser.help_options += get_option_usage(name, usage, value, req, choices);
    if (parser.error != "") return false;
    auto names = get_option_names(name);
    auto pos   = parser.args.end();
    for (auto& name : names) {
        pos = std::min(
            pos, std::find(parser.args.begin(), parser.args.end(), name));
    }
    if (pos == parser.args.end()) {
        if (req) parser.error += "missing value for " + name + "\n";
        return false;
    }
    if (pos == parser.args.end() - 1) {
        parser.error += "missing value for " + name + "\n";
        return false;
    }
    auto vals = *(pos + 1);
    parser.args.erase(pos, pos + 2);
    if (!empty(choices) &&
        std::find(choices.begin(), choices.end(), vals) == choices.end()) {
        parser.error += "bad value for " + name + "\n";
        return false;
    }
    auto new_value = value;
    if (!parse(vals, new_value)) {
        parser.error += "bad value for " + name + "\n";
        return false;
    }
    value = new_value;
    return true;
}

// Parse an argument string. Name should not start with "--" or "-".
template <typename T>
inline bool parse_positional_argument(cmdline_parser& parser, const string& name,
    T& value, const string& usage, bool req, const vector<string>& choices) {
    parser.help_arguments += get_option_usage(name, usage, value, req, choices);
    if (parser.error != "") return false;
    auto pos = std::find_if(parser.args.begin(), parser.args.end(),
        [](auto& v) { return v[0] != '-'; });
    if (pos == parser.args.end()) {
        if (req) parser.error += "missing value for " + name + "\n";
        return false;
    }
    auto vals = *pos;
    parser.args.erase(pos);
    if (!empty(choices) &&
        std::find(choices.begin(), choices.end(), vals) == choices.end()) {
        parser.error += "bad value for " + name + "\n";
        return false;
    }
    auto new_value = value;
    if (!parse(vals, new_value)) {
        parser.error += "bad value for " + name + "\n";
        return false;
    }
    value = new_value;
    return true;
}

// Parse all left argument strings. Name should not start with "--" or "-".
template <typename T>
inline bool parse_positional_arguments(cmdline_parser& parser,
    const string& name, vector<T>& values, const string& usage, bool req) {
    parser.help_arguments += get_option_usage(name, usage, values, req, {});
    if (parser.error != "") return false;
    auto pos = std::find_if(parser.args.begin(), parser.args.end(),
        [](auto& v) { return v[0] != '-'; });
    if (pos == parser.args.end()) {
        if (req) parser.error += "missing value for " + name + "\n";
        return false;
    }
    auto vals = vector<string>{pos, parser.args.end()};
    parser.args.erase(pos, parser.args.end());
    auto new_values = values;
    new_values.resize(vals.size());
    for (auto i = 0; i < vals.size(); i++) {
        if (!parse(vals[i], new_values[i])) {
            parser.error += "bad value for " + name + "\n";
            return false;
        }
    }
    values = new_values;
    return true;
}

// Parse a flag. Name should start with either "--" or "-".
inline bool parse_flag_argument(cmdline_parser& parser, const string& name,
    bool& value, const string& usage, bool req) {
    parser.help_options += get_option_usage(name, usage, false, false, {});
    if (parser.error != "") return false;
    auto names     = get_flag_names(name);
    auto pos       = parser.args.end();
    auto new_value = value;
    for (auto& [name_on, name_off] : names) {
        pos = std::min(
            pos, std::find(parser.args.begin(), parser.args.end(), name_on));
        if (pos != parser.args.end()) {
            new_value = true;
            break;
        }
        pos = std::min(pos,
            std::find(parser.args.begin(), parser.args.end(), name_off));
        if (pos != parser.args.end()) {
            new_value = false;
            break;
        }
    }
    if (pos == parser.args.end()) {
        if (req) parser.error += "missing value for " + name + "\n";
        return false;
    }
    parser.args.erase(pos);
    value = new_value;
    return true;
}

// Parse an integer, float, string. If name starts with "--" or "-", then it is
// an option, otherwise it is a position argument.
template <typename T>
inline bool parse_argument_ref(cmdline_parser& parser, const string& name,
    T& value, const string& usage, bool req) {
    if (is_optional_argument(name)) {
        return parse_option_argument(parser, name, value, usage, req, {});
    } else {
        return parse_positional_argument(parser, name, value, usage, req, {});
    }
}
template <>
inline bool parse_argument_ref<bool>(cmdline_parser& parser, const string& name,
    bool& value, const string& usage, bool req) {
    if (is_optional_flag(name)) {
        return parse_flag_argument(parser, name, value, usage, req);
    } else if (is_optional_argument(name)) {
        return parse_option_argument(parser, name, value, usage, req, {});
    } else {
        return parse_positional_argument(parser, name, value, usage, req, {});
    }
}

// Parse a boolean flag.
inline bool parse_argument_flag_ref(cmdline_parser& parser, const string& name,
    bool& value, const string& usage) {
    auto new_name = (!value) ? name : "/" + name;
    return parse_flag_argument(parser, new_name, value, usage, false);
}

template <typename T>
inline bool parse_argument_ref(cmdline_parser& parser, const string& name,
    T& value, const string& usage, const vector<string>& labels, bool req) {
    auto values = labels.at((int)value);
    auto parsed = false;
    if (is_optional_argument(name)) {
        parsed = parse_option_argument(parser, name, values, usage, req, labels);
    } else {
        parsed = parse_positional_argument(
            parser, name, values, usage, req, labels);
    }
    if (!parsed) return false;
    auto pos = std::find(labels.begin(), labels.end(), values);
    if (pos == labels.end()) return false;
    value = (T)(pos - labels.begin());
    return true;
}

// Parser an argument
template <typename T>
inline bool parse_arguments_ref(cmdline_parser& parser, const string& name,
    vector<T>& values, const string& usage, bool req) {
    return parse_positional_arguments(parser, name, values, usage, req);
}

// Parse an integer, float, string. If name starts with "--" or "-", then it is
// an option, otherwise it is a position argument.
template <typename T>
inline T parse_argument(cmdline_parser& parser, const string& name, T def,
    const string& usage, bool req) {
    auto value = def;
    if (!parse_argument_ref(parser, name, value, usage, req)) return def;
    return value;
}

// Parse a boolean flag.
inline bool parse_argument_flag(
    cmdline_parser& parser, const string& name, bool def, const string& usage) {
    auto value = def;
    if (!parse_argument_flag_ref(parser, name, value, usage)) return def;
    return value;
}

template <typename T>
inline T parse_argument(cmdline_parser& parser, const string& name, T def,
    const string& usage, const vector<string>& labels, bool req) {
    auto value = def;
    if (!parse_argument_ref(parser, name, value, usage, labels, req))
        return def;
    return value;
}

// Parser an argument
template <typename T>
inline vector<T> parse_arguments(cmdline_parser& parser, const string& name,
    const vector<T>& def, const string& usage, bool req) {
    auto values = vector<T>{};
    if (!parse_arguments_ref(parser, name, values, usage, req)) return def;
    return values;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF PATH UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

string normalize_path(const string& filename_) {
    auto filename = filename_;
    for (auto& c : filename)
        if (c == '\\') c = '/';
    if (filename.size() > 1 && filename[0] == '/' && filename[1] == '/') {
        log_error("absolute paths are not supported");
        return filename_;
    }
    if (filename.size() > 3 && filename[1] == ':' && filename[2] == '/' &&
        filename[3] == '/') {
        log_error("absolute paths are not supported");
        return filename_;
    }
    auto pos = (size_t)0;
    while ((pos = filename.find("//")) != filename.npos)
        filename = filename.substr(0, pos) + filename.substr(pos + 1);
    return filename;
}

// Get directory name (not including '/').
string get_dirname(const string& filename_) {
    auto filename = normalize_path(filename_);
    auto pos      = filename.rfind('/');
    if (pos == string::npos) return "";
    return filename.substr(0, pos);
}

// Get extension (not including '.').
string get_extension(const string& filename_) {
    auto filename = normalize_path(filename_);
    auto pos      = filename.rfind('.');
    if (pos == string::npos) return "";
    return filename.substr(pos + 1);
}

// Get filename without directory.
string get_filename(const string& filename_) {
    auto filename = normalize_path(filename_);
    auto pos      = filename.rfind('/');
    if (pos == string::npos) return "";
    return filename.substr(pos + 1);
}

// Replace extension.
string replace_extension(const string& filename_, const string& ext_) {
    auto filename = normalize_path(filename_);
    auto ext      = normalize_path(ext_);
    if (ext.at(0) == '.') ext = ext.substr(1);
    auto pos = filename.rfind('.');
    if (pos == string::npos) return filename;
    return filename.substr(0, pos) + "." + ext;
}

// Check if a file can be opened for reading.
bool exists_file(const string& filename) {
    auto f = fopen(filename.c_str(), "r");
    if (!f) return false;
    fclose(f);
    return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF FILE READING
// -----------------------------------------------------------------------------
namespace yocto {

// log io error
template <typename... Args>
inline void log_io_error(const string& fmt, const Args&... args) {
    log_error(fmt, args...);
}

// File stream wrapper
struct file_stream {
    string filename = "";
    string mode     = "";
    FILE*  fs       = nullptr;

    file_stream()                   = default;
    file_stream(const file_stream&) = delete;
    file_stream& operator=(const file_stream&) = delete;
    file_stream(file_stream&&)                 = default;
    file_stream& operator=(file_stream&&) = default;

    ~file_stream() {
        if (fs) {
            fclose(fs);
            fs = nullptr;
        }
    }

    operator bool() const { return fs; }
};

// Opens a file
inline file_stream open(const string& filename, const string& mode) {
    auto fs = fopen(filename.c_str(), mode.c_str());
    if (!fs) {
        log_io_error("cannot open {}", filename);
        return {};
    }
    return {filename, mode, fs};
}

// Close a file
inline bool close(file_stream& fs) {
    if (!fs) {
        log_io_error("cannot close {}", fs.filename);
        return false;
    }
    fclose(fs.fs);
    fs.fs = nullptr;
    return true;
}

// Gets the length of a file
inline size_t get_length(file_stream& fs) {
    if (!fs) return 0;
    fseek(fs.fs, 0, SEEK_END);
    auto fsize = ftell(fs.fs);
    fseek(fs.fs, 0, SEEK_SET);
    return fsize;
}

// Print to file
inline bool write_text(file_stream& fs, const string& str) {
    if (!fs) return false;
    if (fprintf(fs.fs, "%s", str.c_str()) < 0) {
        log_io_error("cannot write to {}", fs.filename);
        return false;
    }
    return true;
}

// Write to file
template <typename T>
inline bool write_value(file_stream& fs, const T& value) {
    if (!fs) return false;
    if (fwrite(&value, sizeof(T), 1, fs.fs) != 1) {
        log_io_error("cannot write to {}", fs.filename);
        return false;
    }
    return true;
}

// Write to file
template <typename T>
inline bool write_values(file_stream& fs, const vector<T>& vals) {
    if (!fs) return false;
    if (fwrite(data(vals), sizeof(T), vals.size(), fs.fs) != vals.size()) {
        log_io_error("cannot write to {}", fs.filename);
        return false;
    }
    return true;
}

// Write to file
template <typename T>
inline bool write_values(file_stream& fs, size_t num, const T* vals) {
    if (!fs) return false;
    if (fwrite(vals, sizeof(T), num, fs.fs) != num) {
        log_io_error("cannot write to {}", fs.filename);
        return false;
    }
    return true;
}

// Print shortcut
template <typename... Args>
inline bool print(file_stream& fs, const string& fmt, const Args&... args) {
    if (!fs) return false;
    return write_text(fs, format(fmt, args...));
}

// Read binary data to fill the whole buffer
inline bool read_line(file_stream& fs, string& value) {
    if (!fs) return false;
    // TODO: make lkne as large as possible
    value = "";
    char buffer[4096];
    if (!fgets(buffer, 4096, fs.fs)) return false;
    value = string(buffer);
    return true;
}

// Read binary data to fill the whole buffer
template <typename T>
inline bool read_value(file_stream& fs, T& value) {
    if (!fs) return false;
    if (fread(&value, sizeof(T), 1, fs.fs) != 1) {
        log_io_error("cannot read from {}", fs.filename);
        return false;
    }
    return true;
}

// Read binary data to fill the whole buffer
template <typename T>
inline bool read_values(file_stream& fs, vector<T>& vals) {
    if (!fs) return false;
    if (fread(data(vals), sizeof(T), vals.size(), fs.fs) != vals.size()) {
        log_io_error("cannot read from {}", fs.filename);
        return false;
    }
    return true;
}

// Read binary data to fill the whole buffer
template <typename T>
inline bool read_values(file_stream& fs, size_t num, T* vals) {
    if (!fs) return false;
    if (fread(vals, sizeof(T), num, fs.fs) != num) {
        log_io_error("cannot read from {}", fs.filename);
        return false;
    }
    return true;
}

// Load a text file
inline bool load_text(const string& filename, string& str) {
    // https://stackoverflow.com/questions/2602013/read-whole-ascii-file-into-c-stdstring
    auto fs = open(filename, "rb");
    if (!fs) return false;
    auto buffer = vector<char>(get_length(fs));
    if (!read_values(fs, buffer)) return false;
    str = string{buffer.begin(), buffer.end()};
    return true;
}

// Save a text file
inline bool save_text(const string& filename, const string& str) {
    auto fs = open(filename, "wt");
    if (!fs) return false;
    if (!write_text(fs, str)) return false;
    return true;
}

// Load a binary file
inline bool load_binary(const string& filename, vector<byte>& data) {
    // https://stackoverflow.com/questions/2602013/read-whole-ascii-file-into-c-stdstring
    auto fs = open(filename, "rb");
    if (!fs) return false;
    data = vector<byte>(get_length(fs));
    if (!read_values(fs, data)) return false;
    return true;
}

// Save a binary file
inline bool save_binary(const string& filename, const vector<byte>& data) {
    auto fs = open(filename.c_str(), "wb");
    if (!fs) return false;
    if (!write_values(fs, data)) return false;
    return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR CONCURRENCY UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// a simple concurrent queue that locks at every call
template <typename T>
inline concurrent_queue<T>::concurrent_queue() {}
template <typename T>
inline concurrent_queue<T>::concurrent_queue(const concurrent_queue<T>& other) {
    if (!empty(other._queue)) log_error("cannot copy full queue");
    clear();
}
template <typename T>
inline concurrent_queue<T>& concurrent_queue<T>::operator=(
    const concurrent_queue<T>& other) {
    if (!empty(other._queue)) log_error("cannot copy full queue");
    clear();
}

template <typename T>
inline bool concurrent_queue<T>::empty() {
    lock_guard<mutex> lock(_mutex);
    return _queue.empty();
}
template <typename T>
inline void concurrent_queue<T>::clear() {
    lock_guard<mutex> lock(_mutex);
    _queue.clear();
}
template <typename T>
inline void concurrent_queue<T>::push(const T& value) {
    lock_guard<mutex> lock(_mutex);
    _queue.push_back(value);
}
template <typename T>
inline bool concurrent_queue<T>::try_pop(T& value) {
    lock_guard<mutex> lock(_mutex);
    if (_queue.empty()) return false;
    value = _queue.front();
    _queue.pop_front();
    return true;
}

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms.
template <typename Func>
inline void parallel_for(
    int begin, int end, const Func& func, atomic<bool>* cancel, bool serial) {
    if (serial) {
        for (auto idx = begin; idx < end; idx++) {
            if (cancel && *cancel) break;
            func(idx);
        }
    } else {
        auto        threads  = vector<thread>{};
        auto        nthreads = thread::hardware_concurrency();
        atomic<int> next_idx(begin);
        for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
            threads.emplace_back([&func, &next_idx, cancel, end]() {
                while (true) {
                    if (cancel && *cancel) break;
                    auto idx = next_idx.fetch_add(1);
                    if (idx >= end) break;
                    func(idx);
                }
            });
        }
        for (auto& t : threads) t.join();
    }
}

}  // namespace yocto

#endif
