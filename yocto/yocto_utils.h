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
// ## Python-like iterators and collection helpers
//
// This library includes a set of functions to help use C++ collections with
// more ease, inspired by Python. All functions and operators are defined in
// the yocto namespace so they will not affect the code outside. But within
// the Yocto/GL collection they are the best way to do this.
//
// 1. use `range()` to iterato over an integer sequence
// 2. use `enumerate()` to iteratare over a vector and number its elements
// 3. use opeartors + to either concatenate two vectors or a vector and an 
//    element
// 4. use operators += to append an element or a vector to a given vector
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
// 1. Get paths components with `get_dirname()`, `get_filename()` and 
//   `get_extension()`
// 2. Replace the extension with `replace_path_extension()`
// 3. check if a file exists with `exists_file()`
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
#include <thread>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ostream>
#include <istream>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

using std::deque;
using std::lock_guard;
using std::mutex;
using std::thread;
using std::ostream;
using std::istream;
using std::stringstream;
using std::ifstream;
using std::ofstream;
using std::fstream;
using std::cout;
using namespace std::chrono_literals;

using std::getline;

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
inline bool print(const string& fmt, const Args&... args);
template <typename... Args>
inline bool print(FILE* fs, const string& fmt, const Args&... args);
template <typename... Args>
inline bool print(ostream& stream, const string& fmt, const Args&... args);

// Format duration string from nanoseconds
inline string format_duration(int64_t duration);
// Format a large integer number in human readable form
inline string format_num(uint64_t num);

// Parse a list of space separated values.
template <typename... Args>
inline bool parse(const string& str, Args&... args);
template <typename... Args>
inline bool parse(const istream& stream, Args&... args);

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
// IOSTREAM UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Iostream utilities for basic types
inline ostream& operator<<(ostream& os, const vec2f& value);
inline ostream& operator<<(ostream& os, const vec3f& value);
inline ostream& operator<<(ostream& os, const vec4f& value);
inline ostream& operator<<(ostream& os, const vec2i& value);
inline ostream& operator<<(ostream& os, const vec3i& value);
inline ostream& operator<<(ostream& os, const vec4i& value);
inline ostream& operator<<(ostream& os, const mat2f& value);
inline ostream& operator<<(ostream& os, const mat3f& value);
inline ostream& operator<<(ostream& os, const mat4f& value);
inline ostream& operator<<(ostream& os, const frame2f& value);
inline ostream& operator<<(ostream& os, const frame3f& value);
inline ostream& operator<<(ostream& os, const ray2f& value);
inline ostream& operator<<(ostream& os, const ray3f& value);
inline ostream& operator<<(ostream& os, const bbox1f& value);
inline ostream& operator<<(ostream& os, const bbox2f& value);
inline ostream& operator<<(ostream& os, const bbox3f& value);
inline ostream& operator<<(ostream& os, const bbox4f& value);

// Iostream utilities for basic types
inline istream& operator>>(istream& is, vec2f& value);
inline istream& operator>>(istream& is, vec3f& value);
inline istream& operator>>(istream& is, vec4f& value);
inline istream& operator>>(istream& is, vec2i& value);
inline istream& operator>>(istream& is, vec3i& value);
inline istream& operator>>(istream& is, vec4i& value);
inline istream& operator>>(istream& is, mat2f& value);
inline istream& operator>>(istream& is, mat3f& value);
inline istream& operator>>(istream& is, mat4f& value);
inline istream& operator>>(istream& is, frame2f& value);
inline istream& operator>>(istream& is, frame3f& value);
inline istream& operator>>(istream& is, ray2f& value);
inline istream& operator>>(istream& is, ray3f& value);
inline istream& operator>>(istream& is, bbox1f& value);
inline istream& operator>>(istream& is, bbox2f& value);
inline istream& operator>>(istream& is, bbox3f& value);
inline istream& operator>>(istream& is, bbox4f& value);

}

// -----------------------------------------------------------------------------
// PYTHON-LIKE ITERATORS
// -----------------------------------------------------------------------------
namespace yocto {

// Range helpper (this should not be used directly)
struct _range_helper {
    struct _iterator {
        int _pos = 0;
        _iterator& operator++() { _pos ++; return *this;  }
        bool operator!=(const _iterator& other) const { return _pos != other._pos; }
        int operator*() const { return _pos; }
    };
    int _start = 0, _end = 0;
    _iterator begin() const { return {_start}; }
    _iterator end() const { return {_end}; }
};

// Python `range()` equivalent. Construct an object to iterate over a sequence.
inline auto range(int max) { return _range_helper{0, max}; }
inline auto range(int min, int max) { return _range_helper{min, max}; }

// Enumerate helper (this should not be used directly)
template<typename T>
struct _enumerate_helper {
    struct _iterator {
        T* _data = nullptr;
        int _pos = 0;
        _iterator& operator++() { _pos ++; return *this;  }
        bool operator!=(const _iterator& other) const { return _pos != other._pos; }
        pair<int&, T&> operator*() const { return {_pos, *(_data + _pos)}; }
    };
    T* _data = nullptr;
    int _size = 0;
    _iterator begin() const { return {_data, 0}; }
    _iterator end() const { return {_data, _size}; }
};

// Python `enumerate()` equivalent. Construct an object that iteraterates over a
// sequence of elements and numbers them.
template<typename T>
inline auto enumerate(const vector<T>& vals) { return _enumerate_helper<const T>{vals.data(), vals.size()}; };
template<typename T>
inline auto enumerate(vector<T>& vals) { return _enumerate_helper<T>{vals.data(), vals.size()}; };

// Vector append and concatenation
template<typename T>
inline vector<T>& operator+=(vector<T>& a, const vector<T>& b) {
    a.insert(a.end(), b.begin(), b.end());
    return a;
}
template<typename T>
inline vector<T>& operator+=(vector<T>& a, const T& b) {
    a.push_back(b);
    return a;
}
template<typename T>
inline vector<T> operator+(const vector<T>& a, const vector<T>& b) {
    auto c = a;
    return c += b;
}
template<typename T>
inline vector<T> operator+(const vector<T>& a, const T& b) {
    auto c = a;
    return c += b;
}

}

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
// Get directory name (including '/').
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

// Print a value
template<typename T>
inline bool print_value(stringstream& stream, const T& value) {
    stream << value;
    return (bool)stream;
}

// Prints a string.
inline bool print_next(stringstream& stream, const string& fmt) {
    return print_value(stream, fmt);
}
template <typename Arg, typename... Args>
inline bool print_next(
    stringstream& stream, const string& fmt, const Arg& arg, const Args&... args) {
    auto pos = fmt.find("{}");
    if (pos == string::npos) return print_value(stream, fmt);
    if (!print_value(stream, fmt.substr(0, pos))) return false;
    if (!print_value(stream, arg)) return false;
    return print_next(stream, fmt.substr(pos + 2), args...);
}

// Formats a string `fmt` with values taken from `args`. Uses `{}` as
// placeholder.
template <typename... Args>
inline string format(const string& fmt, const Args&... args) {
    auto stream = stringstream();
    print_next(stream, fmt, args...);
    return stream.str();
}

// Prints a string.
template <typename... Args>
inline bool print(const string& fmt, const Args&... args) {
    return print(stdout, fmt, args...);
}
template <typename... Args>
inline bool print(FILE* fs, const string& fmt, const Args&... args) {
    auto str = format(fmt, args...);
    return fprintf(fs, "%s", str.c_str()) >= 0;
}
template <typename... Args>
inline bool print(ostream& stream, const string& fmt, const Args&... args) {
    auto str = format(fmt, args...);
    stream << str;
    return (bool)stream;
}

// Converts to string.
template <typename T>
inline string to_string(const T& value) {
    auto stream = stringstream();
    stream << value;
    return stream.str();
}

// Parse a list of space separated values.
template <typename... Args>
inline bool parse(const string& str, Args&... args) {
    auto stream = stringstream{str};
    if (!parse_next(stream, args...)) return false;
    stream >> std::ws;
    return stream.get() == EOF;
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
inline ofstream& _log_filestream() {
    static auto _log_filestream = ofstream();
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
        cout << msg << "\n";
        fflush(stdout);
    }
    if (_log_filestream()) {
        _log_filestream() << labels[(int)level] << " " << msg << "\n";
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
        _log_filestream().close();
        _log_filestream() = {};
    }
    if (empty(filename)) return;
    _log_filestream().open(filename, append ? std::ios::app : std::ios::out);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IOSTREAM UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Iostream utilities for basic types
inline ostream& operator<<(ostream& os, const vec2f& value) {
    return os << value.x << " " << value.y;
}
inline ostream& operator<<(ostream& os, const vec3f& value) {
    return os << value.x << " " << value.y << " " << value.z;
}
inline ostream& operator<<(ostream& os, const vec4f& value) {
    return os << value.x << " " << value.y << " " << value.z << " " << value.w;
}
inline ostream& operator<<(ostream& os, const vec2i& value) {
    return os << value.x << " " << value.y;
}
inline ostream& operator<<(ostream& os, const vec3i& value) {
    return os << value.x << " " << value.y << " " << value.z;
}
inline ostream& operator<<(ostream& os, const vec4i& value) {
    return os << value.x << " " << value.y << " " << value.z << " " << value.w;
}
inline ostream& operator<<(ostream& os, const mat2f& value) {
    return os << value.x << " " << value.y;
}
inline ostream& operator<<(ostream& os, const mat3f& value) {
    return os << value.x << " " << value.y << " " << value.z;
}
inline ostream& operator<<(ostream& os, const mat4f& value) {
    return os << value.x << " " << value.y << " " << value.z << " " << value.w;
}
inline ostream& operator<<(ostream& os, const frame2f& value) {
    return os << value.x << " " << value.y << " " << value.o;
}
inline ostream& operator<<(ostream& os, const frame3f& value) {
    return os << value.x << " " << value.y << " " << value.z << " " << value.o;
}
inline ostream& operator<<(ostream& os, const ray2f& value) {
    return os << value.o << " " << value.d << " " << value.tmin << " " << value.tmax;
}
inline ostream& operator<<(ostream& os, const ray3f& value) {
    return os << value.o << " " << value.d << " " << value.tmin << " " << value.tmax;
}
inline ostream& operator<<(ostream& os, const bbox1f& value) {
    return os << value.min << " " << value.max;
}
inline ostream& operator<<(ostream& os, const bbox2f& value) {
    return os << value.min << " " << value.max;
}
inline ostream& operator<<(ostream& os, const bbox3f& value) {
    return os << value.min << " " << value.max;
}
inline ostream& operator<<(ostream& os, const bbox4f& value) {
    return os << value.min << " " << value.max;
}

// Iostream utilities for basic types
inline istream& operator>>(istream& is, vec2f& value) {
    return is >> value.x >> value.y;
}
inline istream& operator>>(istream& is, vec3f& value) {
    return is >> value.x >> value.y >> value.z;
}
inline istream& operator>>(istream& is, vec4f& value) {
    return is >> value.x >> value.y >> value.z >> value.w;
}
inline istream& operator>>(istream& is, vec2i& value) {
    return is >> value.x >> value.y;
}
inline istream& operator>>(istream& is, vec3i& value) {
    return is >> value.x >> value.y >> value.z;
}
inline istream& operator>>(istream& is, vec4i& value) {
    return is >> value.x >> value.y >> value.z >> value.w;
}
inline istream& operator>>(istream& is, mat2f& value) {
    return is >> value.x >> value.y;
}
inline istream& operator>>(istream& is, mat3f& value) {
    return is >> value.x >> value.y >> value.z;
}
inline istream& operator>>(istream& is, mat4f& value) {
    return is >> value.x >> value.y >> value.z >> value.w;
}
inline istream& operator>>(istream& is, frame2f& value) {
    return is >> value.x >> value.y >> value.o;
}
inline istream& operator>>(istream& is, frame3f& value) {
    return is >> value.x >> value.y >> value.z >> value.o;
}
inline istream& operator>>(istream& is, ray2f& value) {
    return is >> value.o >> value.d >> value.tmin >> value.tmax;
}
inline istream& operator>>(istream& is, ray3f& value) {
    return is >> value.o >> value.d >> value.tmin >> value.tmax;
}
inline istream& operator>>(istream& is, bbox1f& value) {
    return is >> value.min >> value.max;
}
inline istream& operator>>(istream& is, bbox2f& value) {
    return is >> value.min >> value.max;
}
inline istream& operator>>(istream& is, bbox3f& value) {
    return is >> value.min >> value.max;
}
inline istream& operator>>(istream& is, bbox4f& value) {
    return is >> value.min >> value.max;
}

}

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
    auto stream = stringstream{};
    stream << value;
    return stream.str();
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
    auto usage = ""s;
    usage += parser.help_command + ": " + parser.help_usage + "\n";
    usage += "usage: " + parser.help_command;
    if(!(empty(parser.help_options))) usage += "[options] ";
    if(!(empty(parser.help_arguments))) usage += "arguments";
    usage += "\n\n";
    if (!empty(parser.help_options))
        usage += "options:\n" + parser.help_options + "\n";
    if (!empty(parser.help_arguments))
        usage += "arguments:\n" + parser.help_arguments + "\n";
    cout << usage;
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
        cout << "error: " + parser.error + "\n";
        print_cmdline_usage(parser);
        exit(1);
    }
}

// Parse option value
inline bool parse_option_value(const string& vals, string& val) {
    val = vals;
    return true;
}
template<typename T>
inline bool parse_option_value(const string& vals, T& val) {
    auto stream = stringstream{vals};
    stream >> val;
    if(stream.fail()) return false;
    return true;
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
    if (!parse_option_value(vals, new_value)) {
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
    if (!parse_option_value(vals, new_value)) {
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
        if (!parse_option_value(vals[i], new_values[i])) {
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

// Get directory name (including '/').
string get_dirname(const string& filename_) {
    auto filename = normalize_path(filename_);
    auto pos      = filename.rfind('/');
    if (pos == string::npos) return "";
    return filename.substr(0, pos+1);
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

// write value to a stream
template<typename T>
inline ostream& write_value(ostream& stream, const T& value) {
    return stream.write((char*)&value, sizeof(T));
}

// write values to a stream
template<typename T>
inline ostream& write_values(ostream& stream, const vector<T>& values) {
    if(values.empty()) return stream;
    return stream.write((char*)values.data(), values.size()*sizeof(T));
}

// Read binary data to fill the whole buffer
template <typename T>
inline istream& read_value(istream& stream, T& value) {
    return stream.read((char*)&value, sizeof(T));
}

// Read binary data to fill the whole buffer
template <typename T>
inline istream& read_values(istream& stream, vector<T>& values) {
    if(values.empty()) return stream;
    return stream.read((char*)values.data(), values.size()*sizeof(T));
}

// Load a text file
inline bool load_text(const string& filename, string& str) {
    // https://stackoverflow.com/questions/2602013/read-whole-ascii-file-into-c-stdstring
    auto stream = ifstream(filename);
    if(!stream) {
        log_io_error("cannot open file {}", filename);
        return false;
    }
    stringstream buffer;
    buffer << stream.rdbuf();
    if(stream.fail()) {
        log_io_error("cannot read file {}", filename);
        return false;
    }
    str = buffer.str();
    return true;
}

// Save a text file
inline bool save_text(const string& filename, const string& str) {
    auto stream = ofstream(filename);
    if(!stream) {
        log_io_error("cannot open file {}", filename);
        return false;
    }
    stream << str;
    if(!stream) {
        log_io_error("cannot write file {}", filename);
        return false;
    }
    return true;
}

// Load a binary file
inline bool load_binary(const string& filename, vector<byte>& data) {
    // https://stackoverflow.com/questions/2602013/read-whole-ascii-file-into-c-stdstring
    auto stream = ifstream(filename, std::ios::binary);
    if(!stream) {
        log_io_error("cannot open file {}", filename);
        return false;
    }
    stringstream buffer;
    buffer << stream.rdbuf();
    if(stream.fail()) {
        log_io_error("cannot read file {}", filename);
        return false;
    }
    auto str = buffer.str();
    data = vector<byte>((byte*)str.data(), (byte*)str.data() + str.size());
    return true;
}

// Save a binary file
inline bool save_binary(const string& filename, const vector<byte>& data) {
    auto stream = ofstream(filename, std::ios::binary);
    if(!stream) {
        log_io_error("cannot open file {}", filename);
        return false;
    }
    stream.write((char*)data.data(), data.size());
    if(!stream) {
        log_io_error("cannot write file {}", filename);
        return false;
    }
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
