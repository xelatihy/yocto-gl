//
// # Yocto/Utils: Tiny collection of utilities to support Yocto/GL
//
//
// Yocto/Utils is a collection of utilities used in writing other Yocto/GL
// libraries and example applications. We support printing builtin and
// Yocto/Math values, parsing command line arguments, simple path
// manipulation, file lading/saving and basic concurrency utilities.
//
//
// ## Printing and parsing values
//
// Use `print_value()` to write a string in a stream or `println_values()`
// to print a line of values. Use `format_duraction()` and `format_num()`
// for pretty printing times and numbers. These will change once lib `fmt`
// is accepted in the standard.
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
// 1. initialize the parser with `init_cmdline_parser(parser, argc, argv, help)`
// 2. read a value with `value = parse_cmdline_argument(parser, name, default,
// help)`
//    - is name starts with '--' or '-' then it is an option
//    - otherwise it is a positional arguments
//    - options and arguments may be intermixed
//    - the type of each option is determined by the default value `default`
//    - the value is parsed on the stop
// 3. finished parsing with `check_cmdline_parser(parser)`
//    - if an error occurred, the parser will exit and print a usage message
//
//
// ## Path manipulation
//
// We define a few path manipulation utilities to split and join path
// components.
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
// 3. use `input_file` and `output_file` as RIIA FILE* wrappers
// 4. use `read_XXX()` and `write_XXX()` to read/write to this files
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
#include <future>
#include <mutex>
#include <string>
#include <thread>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

using std::deque;
using std::lock_guard;
using std::mutex;
using std::thread;
using std::future;
using std::async;
using namespace std::chrono_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// APPLICATION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Log info/error/fatal/trace message
inline void exit_error(const char* msg) {
    printf("%s\n", msg);
    exit(1);
}
inline void exit_error(const string& msg) {
    printf("%s\n", msg.c_str());
    exit(1);
}

// Format duration string from nanoseconds
inline string format_duration(int64_t duration);
// Format a large integer number in human readable form
inline string format_num(uint64_t num);

// get time in nanoseconds - useful only to compute difference of times
inline int64_t get_time() {
    return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// FILE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// io error
struct io_error : runtime_error {
    explicit io_error(const char* msg) : runtime_error{msg} {}
    explicit io_error(const std::string& msg) : runtime_error{msg} {}
};

// file inout stream
struct input_file {
    input_file(const string& filename, bool binary = false) {
        this->filename = filename;
        file           = fopen(filename.c_str(), binary ? "rb" : "rt");
        if (!file) throw io_error("could not open " + filename);
    }
    input_file(FILE* fs) {
        file  = fs;
        owned = false;
    }

    input_file(const input_file&) = delete;
    input_file& operator=(const input_file&) = delete;

    ~input_file() {
        if (file && owned) fclose(file);
    }

    string filename = "";
    FILE*  file     = nullptr;
    bool   owned    = true;
};

// file writer
struct output_file {
    output_file(const string& filename, bool binary = false) {
        this->filename = filename;
        file           = fopen(filename.c_str(), binary ? "wb" : "wt");
        if (!file) throw io_error("could not open " + filename);
    }
    output_file(FILE* fs) {
        file  = fs;
        owned = false;
    }

    output_file(const output_file&) = delete;
    output_file& operator=(const output_file&) = delete;

    ~output_file() {
        if (file && owned) fclose(file);
    }

    string filename = "";
    FILE*  file     = nullptr;
    bool   owned    = true;
};

// write a value to a file
template <typename T>
inline void write_value(const output_file& fs, const T& value) {
    if (fwrite(&value, sizeof(value), 1, fs.file) != 1) {
        throw io_error("cannot write to " + fs.filename);
    }
}

// write values to a file
template <typename T>
inline void write_values(const output_file& fs, const vector<T>& values) {
    if (values.empty()) return;
    if (fwrite(values.data(), sizeof(values[0]), values.size(), fs.file) !=
        values.size()) {
        throw io_error("cannot write to " + fs.filename);
    }
}
template <typename T>
inline void write_values(const output_file& fs, const T* values, size_t count) {
    if (!count) return;
    if (fwrite(values, sizeof(values[0]), count, fs.file) != count) {
        throw io_error("cannot write to " + fs.filename);
    }
}

// write text to a file
inline void write_text(const output_file& fs, const std::string& str) {
    if (fprintf(fs.file, "%s", str.c_str()) < 0) {
        throw io_error("cannot write to " + fs.filename);
    }
}

// read a value from a file
template <typename T>
inline void read_value(const input_file& fs, T& value) {
    if (fread(&value, sizeof(value), 1, fs.file) != 1) {
        throw io_error("cannot read from " + fs.filename);
    }
}

// read values from a file
template <typename T>
inline void read_values(const input_file& fs, T* values, size_t count) {
    if (!count) return;
    if (fread(values, sizeof(values[0]), count, fs.file) != count) {
        throw io_error("cannot read from " + fs.filename);
    }
}
template <typename T>
inline void read_values(const input_file& fs, vector<T>& values) {
    if (values.empty()) return;
    if (fread(values.data(), sizeof(values[0]), values.size(), fs.file) !=
        values.size()) {
        throw io_error("cannot read from " + fs.filename);
    }
}
// read characters from a file
inline void read_values(const input_file& fs, string& values) {
    if (values.empty()) return;
    if (fread(values.data(), sizeof(values[0]), values.size(), fs.file) !=
        values.size()) {
        throw io_error("cannot read from " + fs.filename);
    }
}

// read a line of text
inline bool read_line(const input_file& fs, string& str) {
    char buffer[4096];
    if (fgets(buffer, sizeof(buffer), fs.file) == nullptr) return false;
    str = buffer;
    return true;
}
inline bool read_line(const input_file& fs, char* buffer, size_t size) {
    if (fgets(buffer, size, fs.file) == nullptr) return false;
    return true;
}

// Printing values
inline void print_value(const output_file& fs, int value) {
    if (fprintf(fs.file, "%d", value) < 0)
        throw io_error("cannot write to file " + fs.filename);
}
inline void print_value(const output_file& fs, bool value) {
    if (fprintf(fs.file, "%d", (int)value) < 0)
        throw io_error("cannot write to file " + fs.filename);
}
inline void print_value(const output_file& fs, float value) {
    if (fprintf(fs.file, "%g", value) < 0)
        throw io_error("cannot write to file " + fs.filename);
}
inline void print_value(const output_file& fs, char value) {
    if (fprintf(fs.file, "%c", value) < 0)
        throw io_error("cannot write to file " + fs.filename);
}
inline void print_value(const output_file& fs, const char* value) {
    if (fprintf(fs.file, "%s", value) < 0)
        throw io_error("cannot write to file " + fs.filename);
}
inline void print_value(const output_file& fs, const string& value) {
    if (fprintf(fs.file, "%s", value.c_str()) < 0)
        throw io_error("cannot write to file " + fs.filename);
}
template <typename T, int N>
inline void print_value(const output_file& fs, const vec<T, N>& value) {
    for (auto i = 0; i < N; i++) {
        if (i) print_value(fs, ' ');
        print_value(fs, value[i]);
    }
}
template <typename T, int N, int M>
inline void print_value(const output_file& fs, const mat<T, N, M>& value) {
    for (auto i = 0; i < M; i++) {
        if (i) print_value(fs, ' ');
        print_value(fs, value[i]);
    }
}
template <typename T, int N>
inline void print_value(const output_file& fs, const frame<T, N>& value) {
    for (auto i = 0; i < N + 1; i++) {
        if (i) print_value(fs, ' ');
        print_value(fs, value[i]);
    }
}

// print values to file
template <typename Arg, typename... Args>
inline void println_values(
    const output_file& fs, const Arg& value, const Args&... values) {
    print_value(fs, value);
    if constexpr (sizeof...(values) > 0) {
        print_value(fs, ' ');
        println_values(fs, values...);
    } else {
        print_value(fs, '\n');
    }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PYTHON-LIKE ITERATORS
// -----------------------------------------------------------------------------
namespace yocto {

// Range helpper (this should not be used directly)
struct _range_helper {
    struct _iterator {
        int        _pos = 0;
        _iterator& operator++() {
            _pos++;
            return *this;
        }
        bool operator!=(const _iterator& other) const {
            return _pos != other._pos;
        }
        int operator*() const { return _pos; }
    };
    int       _start = 0, _end = 0;
    _iterator begin() const { return {_start}; }
    _iterator end() const { return {_end}; }
};

// Python `range()` equivalent. Construct an object to iterate over a sequence.
inline auto range(int max) { return _range_helper{0, max}; }
inline auto range(int min, int max) { return _range_helper{min, max}; }

// Enumerate helper (this should not be used directly)
template <typename T>
struct _enumerate_helper {
    struct _iterator {
        T*         _data = nullptr;
        int        _pos  = 0;
        _iterator& operator++() {
            _pos++;
            return *this;
        }
        bool operator!=(const _iterator& other) const {
            return _pos != other._pos;
        }
        pair<int&, T&> operator*() const { return {_pos, *(_data + _pos)}; }
    };
    T*        _data = nullptr;
    int       _size = 0;
    _iterator begin() const { return {_data, 0}; }
    _iterator end() const { return {_data, _size}; }
};

// Python `enumerate()` equivalent. Construct an object that iteraterates over a
// sequence of elements and numbers them.
template <typename T>
inline auto enumerate(const vector<T>& vals) {
    return _enumerate_helper<const T>{vals.data(), vals.size()};
};
template <typename T>
inline auto enumerate(vector<T>& vals) {
    return _enumerate_helper<T>{vals.data(), vals.size()};
};

// Vector append and concatenation
template <typename T>
inline vector<T>& operator+=(vector<T>& a, const vector<T>& b) {
    a.insert(a.end(), b.begin(), b.end());
    return a;
}
template <typename T>
inline vector<T>& operator+=(vector<T>& a, const T& b) {
    a.push_back(b);
    return a;
}
template <typename T>
inline vector<T> operator+(const vector<T>& a, const vector<T>& b) {
    auto c = a;
    return c += b;
}
template <typename T>
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
inline void init_cmdline_parser(cmdline_parser& parser, int argc, char** argv,
    const string& usage, const string& cmd = "", bool add_help_flag = true,
    bool add_logging_flags = true);
// check if any error occurred and exit appropriately
inline void check_cmdline_parser(cmdline_parser& parser);

// Parse an int, float, string, vecXX and bool option or positional argument.
// Options's names starts with "--" or "-", otherwise they are arguments.
// vecXX options use space-separated values but all in one argument
// (use " or ' from the common line).
// Boolean flags are indicated with a pair of names "--name/--no-name", so
// that we have both options available. You can also use the parse flag function
// in which case only one name is used and the flag will flip the value passed.
template <typename T>
inline T    parse_cmdline_argument(cmdline_parser& parser, const string& name,
       T def, const string& usage, bool req = false);
inline bool parse_cmdline_argument_flag(
    cmdline_parser& parser, const string& name, bool def, const string& usage);
// Parse all arguments left on the command line.
template <typename T>
inline vector<T> parse_cmdline_arguments(cmdline_parser& parser,
    const string& name, const vector<T>& def, const string& usage,
    bool req = false);
// Parse a labeled enum, with enum values that are successive integers.
template <typename T>
inline T parse_cmdline_argument(cmdline_parser& parser, const string& name,
    T def, const string& usage, const vector<string>& labels, bool req = false);

// Parse an int, float, string, vecXX and bool option or positional argument.
// Options's names starts with "--" or "-", otherwise they are arguments.
// vecXX options use space-separated values but all in one argument
// (use " or ' from the common line). Booleans are flags.
// Boolean flags are indicated with a pair of names "--name/--no-name", so
// that we have both options available. You can also use the parse flag function
// in which case only one name is used and the flag will flip the value passed.
template <typename T>
inline bool parse_cmdline_argument_ref(cmdline_parser& parser,
    const string& name, T& val, const string& usage, bool req = false);
inline bool parse_cmdline_argument_flag(
    cmdline_parser& parser, const string& name, bool& val, const string& usage);
// Parse all arguments left on the command line.
template <typename T>
inline bool parse_cmdline_arguments_ref(cmdline_parser& parser,
    const string& name, vector<T>& val, const string& usage, bool req = false);
// Parse a labeled enum, with enum values that are successive integers.
template <typename T>
inline bool parse_cmdline_argument_ref(cmdline_parser& parser,
    const string& name, T& val, const string& usage,
    const vector<string>& labels, bool req = false);

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
inline void load_text(const string& filename, string& str);
inline void save_text(const string& filename, const string& str);

// Load/save a binary file
inline void load_binary(const string& filename, vector<byte>& data);
inline void save_binary(const string& filename, const vector<byte>& data);

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
    parallel_for(
        0, (int)values.size(), [&func, &values](int idx) { func(values[idx]); },
        cancel, serial);
}
template <typename T, typename Func>
inline void parallel_foreach(const vector<T>& values, const Func& func,
    atomic<bool>* cancel = nullptr, bool serial = false) {
    parallel_for(
        0, (int)values.size(), [&func, &values](int idx) { func(values[idx]); },
        cancel, serial);
}

}  // namespace yocto

// ---------------------------------------------------------------------------//
//                                                                            //
//                             IMPLEMENTATION                                 //
//                                                                            //
// ---------------------------------------------------------------------------//

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
inline void init_cmdline_parser(cmdline_parser& parser, int argc, char** argv,
    const string& usage, const string& cmd, bool add_help_flag,
    bool add_logging_flags) {
    parser               = {};
    parser.args          = {argv + 1, argv + argc};
    parser.help_command  = (cmd.empty()) ? argv[0] : cmd;
    parser.help_usage    = usage;
    parser.add_help_flag = add_help_flag;
}

// check if option or argument
inline bool is_optional_argument(const string& name) {
    return name.size() > 1 && name.front() == '-';
}

// check if flag
inline bool is_optional_flag(const string& name) {
    return name.size() > 1 && name.front() == '-' &&
           name.find('/') != name.npos;
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
    return std::to_string(value);
}
inline string get_option_default_string(const bool& value) {
    return (value) ? "true"s : "false"s;
}
inline string get_option_default_string(const string& value) { return value; }
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
    if (!choices.empty()) {
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
    if (!(parser.help_options.empty())) usage += "[options] ";
    if (!(parser.help_arguments.empty())) usage += "arguments";
    usage += "\n\n";
    if (!parser.help_options.empty())
        usage += "options:\n" + parser.help_options + "\n";
    if (!parser.help_arguments.empty())
        usage += "arguments:\n" + parser.help_arguments + "\n";
    printf("%s\n", usage.c_str());
}

// Parse a flag. Name should start with either "--" or "-".
inline bool parse_flag_argument(cmdline_parser& parser, const string& name,
    bool& value, const string& usage, bool req);

// check if any error occurred and exit appropriately
inline void check_cmdline_parser(cmdline_parser& parser) {
    if (parser.add_help_flag) {
        auto help = false;
        if (parse_flag_argument(
                parser, "--help,-?", help, "print help", false)) {
            print_cmdline_usage(parser);
            exit(0);
        }
    }
    if (!parser.args.empty()) {
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
    if (!parser.error.empty()) {
        printf("error: %s\n", parser.error.c_str());
        print_cmdline_usage(parser);
        exit(1);
    }
}

// Parse option value
inline bool parse_option_value(const string& vals, string& val) {
    val = vals;
    return true;
}
inline bool parse_option_value(const string& vals, int& val) {
    return sscanf(vals.c_str(), "%d", &val) > 0;
}
inline bool parse_option_value(const string& vals, bool& val) {
    auto vali = 0;
    if (sscanf(vals.c_str(), "%d", &vali) < 0) return false;
    val = (bool)vali;
    return true;
}
inline bool parse_option_value(const string& vals, float& val) {
    return sscanf(vals.c_str(), "%f", &val) > 0;
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
    if (!choices.empty() &&
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
inline bool parse_positional_argument(cmdline_parser& parser,
    const string& name, T& value, const string& usage, bool req,
    const vector<string>& choices) {
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
    if (!choices.empty() &&
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
        pos = std::min(
            pos, std::find(parser.args.begin(), parser.args.end(), name_off));
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
inline bool parse_cmdline_argument_ref(cmdline_parser& parser,
    const string& name, T& value, const string& usage, bool req) {
    if (is_optional_argument(name)) {
        return parse_option_argument(parser, name, value, usage, req, {});
    } else {
        return parse_positional_argument(parser, name, value, usage, req, {});
    }
}
template <>
inline bool parse_cmdline_argument_ref<bool>(cmdline_parser& parser,
    const string& name, bool& value, const string& usage, bool req) {
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
inline bool parse_cmdline_argument_ref(cmdline_parser& parser,
    const string& name, T& value, const string& usage,
    const vector<string>& labels, bool req) {
    auto values = labels.at((int)value);
    auto parsed = false;
    if (is_optional_argument(name)) {
        parsed = parse_option_argument(
            parser, name, values, usage, req, labels);
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
inline bool parse_cmdline_arguments_ref(cmdline_parser& parser,
    const string& name, vector<T>& values, const string& usage, bool req) {
    return parse_positional_arguments(parser, name, values, usage, req);
}

// Parse an integer, float, string. If name starts with "--" or "-", then it is
// an option, otherwise it is a position argument.
template <typename T>
inline T parse_cmdline_argument(cmdline_parser& parser, const string& name,
    T def, const string& usage, bool req) {
    auto value = def;
    if (!parse_cmdline_argument_ref(parser, name, value, usage, req))
        return def;
    return value;
}

// Parse a boolean flag.
inline bool parse_cmdline_argument_flag(
    cmdline_parser& parser, const string& name, bool def, const string& usage) {
    auto value = def;
    if (!parse_argument_flag_ref(parser, name, value, usage)) return def;
    return value;
}

template <typename T>
inline T parse_cmdline_argument(cmdline_parser& parser, const string& name,
    T def, const string& usage, const vector<string>& labels, bool req) {
    auto value = def;
    if (!parse_cmdline_argument_ref(parser, name, value, usage, labels, req))
        return def;
    return value;
}

// Parser an argument
template <typename T>
inline vector<T> parse_cmdline_arguments(cmdline_parser& parser,
    const string& name, const vector<T>& def, const string& usage, bool req) {
    auto values = vector<T>{};
    if (!parse_cmdline_arguments_ref(parser, name, values, usage, req))
        return def;
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
        throw std::invalid_argument("absolute paths are not supported");
        return filename_;
    }
    if (filename.size() > 3 && filename[1] == ':' && filename[2] == '/' &&
        filename[3] == '/') {
        throw std::invalid_argument("absolute paths are not supported");
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
    return filename.substr(0, pos + 1);
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

// Load a text file
inline void load_text(const string& filename, string& str) {
    // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
    auto fs = input_file(filename);
    fseek(fs.file, 0, SEEK_END);
    auto length = ftell(fs.file);
    fseek(fs.file, 0, SEEK_SET);
    str.resize(length);
    read_values(fs, str);
}

// Save a text file
inline void save_text(const string& filename, const string& str) {
    auto fs = output_file(filename);
    write_text(fs, str);
}

// Load a binary file
inline void load_binary(const string& filename, vector<byte>& data) {
    // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
    auto fs = input_file(filename, true);
    fseek(fs.file, 0, SEEK_END);
    auto length = ftell(fs.file);
    fseek(fs.file, 0, SEEK_SET);
    data.resize(length);
    read_values(fs, data);
}

// Save a binary file
inline void save_binary(const string& filename, const vector<byte>& data) {
    auto fs = output_file(filename, true);
    write_values(fs, data);
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
    if (!other._queue.empty())
        throw std::invalid_argument("cannot copy full queue");
    clear();
}
template <typename T>
inline concurrent_queue<T>& concurrent_queue<T>::operator=(
    const concurrent_queue<T>& other) {
    if (!other._queue.empty())
        throw std::invalid_argument("cannot copy full queue");
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
        auto        futures  = vector<future<void>>{};
        auto        nthreads = thread::hardware_concurrency();
        atomic<int> next_idx(begin);
        for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
            futures.emplace_back(
                async(std::launch::async, [&func, &next_idx, cancel, end]() {
                    while (true) {
                        if (cancel && *cancel) break;
                        auto idx = next_idx.fetch_add(1);
                        if (idx >= end) break;
                        func(idx);
                    }
                }));
        }
        for (auto& f : futures) f.get();
    }
}

}  // namespace yocto

#endif
