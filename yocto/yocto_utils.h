///
/// # Yocto/Utils
///
/// Utilities for writing command line applications, mostly
/// a simple to use command line parser, a logger, a thread-pool,
/// and string, path and file functions, and a timer. All functions are defined
/// in separate namespaces.
///
/// ## Usage for Command Line Parsing
///
/// 1. namespace logging
/// 2. create a parser object
///     - an option for printing help is automatically added
///     `parser = make_parser(argc, argv, program description)`
/// 3. for each option, parse it calling the functions `parse_opt()`
///     - options are parsed on the fly and a comprehensive help is
///       automatically generated
///     - supports bool (flags), int, float, double, std::string, enums
///     - options names are "--longname" for longname and "-s" for short
///     - command line format is "--longname value", "-s v" for all but flags
///     - values are parsed with `iostream <<` operators
///     - for general use `opt = parse_opt<type>()`
///     - for boolean flags is `parse_flag()`
///     - for enums use `parse_opte()`
/// 4. for each unnamed argument, parse it calling the functions parse_arg()
///     - names are only used for help
///     - supports types as above
///     - for general use `arg = parse_arg<type>()`
///     - to parse all remaining values use `args = parse_arga<type>(...)`
/// 5. end cmdline parsing with `check_parser()` to check for unsued values,
///    missing arguments and print help if needed
/// 6. since arguments are parsed immediately, one can easily implement
///    subcommands by just branching the command line code based on a read
///    argument without any need for complex syntax
///
/// Notes: the end of this file contains a test function that also
/// illustreates the library usage.
///
/// ## Usage for Logging
///
/// 1. namespace logging
/// 2. create loggers with `make_file_logger()`, `make_stderr_logger()`,
///    `make_stdout_logger()`
/// 3. you can set default loggers with `get_default_loggers()`; note that none
///    are set by default
/// 4. write log messages with `log_msg()` and its variants.
///
/// ## Usage for Concurrent Execution
///
/// 1. namespace concurrent
/// 2. either create a thread pool `make_thread_pool` or use the global one
/// 3. run tasks in parallel `thread_pool_for()`
/// 4. run tasks asynchronously `thread_pool_async()`
///
/// ## Utilities
///
/// 1. filename splitting functions in namespace `path`
/// 2. loading and save entire files in namespace `file`
/// 3. Python-line string manipulation in namespace `string`
/// 4. Python-like operators for standard containers in namespace `operators`
/// 5. simple timer in namespace `timer`
/// 6. hashing functions in the `hashing` namespace (SHA1 and xxHash)
///
///
/// ## History
///
/// - v 0.22: simpler logging
/// - v 0.21: move to header-only mode
/// - v 0.20: simpler logging
/// - v 0.19: some containers ops
/// - v 0.18: timer
/// - v 0.17: renamed to yocto utils
/// - v 0.16: split into namespaces
/// - v 0.15: remove inline compilation
/// - v 0.14: Python-like operator for std::vector
/// - v 0.13: more file and string utilities
/// - v 0.12: better thread pool implementation
/// - v 0.11: added a few more path utilities
/// - v 0.10: changed default name for help option; better help printing
/// - v 0.9: C-like string formatting
/// - v 0.8: switch to .h/.cpp pair (templated functions are specialized)
/// - v 0.7: logging
/// - v 0.6: doxygen comments
/// - v 0.5: added a few python-like std::string manipulation functions.
/// - v 0.4: [major API change] move to modern C++ interface
/// - v 0.3: adding a few functions for path splitting.
/// - v 0.2: [API change] C++ API
/// - v 0.1: C++ implementation
/// - v 0.0: initial release in C99
///
namespace yu {}

//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
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
// LICENSE OF INCLUDED SOFTWARE for ThreadPool code from LLVM code base
//
// Copyright (c) 2003-2016 University of Illinois at Urbana-Champaign.
// All rights reserved.
//
// Developed by:
//
//     LLVM Team
//
//     University of Illinois at Urbana-Champaign
//
//     http://llvm.org
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// with the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
//     * Redistributions of source code must retain the above copyright notice,
//       this list of conditions and the following disclaimers.
//
//     * Redistributions in binary form must reproduce the above copyright
//     notice,
//       this list of conditions and the following disclaimers in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the names of the LLVM Team, University of Illinois at
//       Urbana-Champaign, nor the names of its contributors may be used to
//       endorse or promote products derived from this Software without specific
//       prior written permission.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
// CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH
// THE SOFTWARE.

#ifndef _YU_H_
#define _YU_H_

#include <cassert>
#include <chrono>
#include <deque>
#include <functional>
#include <future>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-security"
#endif

///
/// General utilities to write command line applications
///
namespace yu {

///
/// Command line parsing
///
namespace cmdline {

///
/// Command line parser.
///
struct parser;

///
/// Inits a command line parser.
///
inline parser make_parser(const std::vector<std::string>& args,
    const std::string& name = "", const std::string& help = "");

///
/// Inits a command line parser.
///
inline parser make_parser(
    int argc, char* argv[], const std::string& name, const std::string& help);

///
/// Ends parsing checking for error for unused options or arguments.
/// Exit if needed.
///
inline bool check_parser(parser& par);

///
/// Parses an optional flag as described in the intro.
///
inline bool parse_flag(parser& par, const std::string& longname,
    const std::string& shortname, const std::string& help, bool def = false);

///
/// Parses an option as described in the intro.
///
template <typename T>
inline T parse_opt(parser& par, const std::string& longname,
    const std::string& shortname, const std::string& help, const T& def,
    bool required = false, const std::vector<T>& choices = {});

///
/// Specialization of parse_opt()
///
inline std::string parse_opts(parser& par, const std::string& longname,
    const std::string& shortname, const std::string& help,
    const std::string& def, bool required = false,
    const std::vector<std::string>& choices = {}) {
    return parse_opt<std::string>(
        par, longname, shortname, help, def, required, choices);
}

///
/// Specialization of parse_opt()
///
inline int parse_opti(parser& par, const std::string& longname,
    const std::string& shortname, const std::string& help, int def,
    bool required = false, const std::vector<int>& choices = {}) {
    return parse_opt<int>(
        par, longname, shortname, help, def, required, choices);
}

///
/// Specialization of parse_opt()
///
inline float parse_optf(parser& par, const std::string& longname,
    const std::string& shortname, const std::string& help, float def,
    bool required = false, const std::vector<float>& choices = {}) {
    return parse_opt<float>(
        par, longname, shortname, help, def, required, choices);
}

///
/// Specialization of parse_opt()
///
inline double parse_optd(parser& par, const std::string& longname,
    const std::string& shortname, const std::string& help, double def,
    bool required = false, const std::vector<double>& choices = {}) {
    return parse_opt<double>(
        par, longname, shortname, help, def, required, choices);
}

///
/// Parses an option enum as described in the intro.
///
inline int parse_opte(parser& par, const std::string& longname,
    const std::string& shortname, const std::string& help, int def,
    const std::vector<std::pair<std::string, int>>& vals,
    bool required = false);

///
/// Parses an option enum as described in the intro.
///
template <typename T>
inline T parse_opte(parser& par, const std::string& longname,
    const std::string& shortname, const std::string& help, T def,
    const std::vector<std::pair<std::string, T>>& vals, bool required = false) {
    return (T)parse_opte(par, longname, shortname, help, (int)def,
        (const std::vector<std::pair<std::string, int>>&)vals, required);
}

///
/// Parses an option array as described in the intro.
///
template <typename T>
inline T parse_opta(parser& par, const std::string& longname,
    const std::string& shortname, const std::string& help, const T& def,
    int nargs, bool required = false, const std::vector<T>& choices = {});

///
/// Parses an argument as described in the intro.
///
template <typename T>
inline T parse_arg(parser& par, const std::string& longname,
    const std::string& help, const T& def, bool required = false,
    const std::vector<T>& choices = {});

///
/// Specialization of parse_arg()
///
inline std::string parse_args(parser& par, const std::string& longname,
    const std::string& help, const std::string& def, bool required = false,
    const std::vector<std::string>& choices = {}) {
    return parse_arg<std::string>(par, longname, help, def, required, choices);
}

///
/// Specialization of parse_arg()
///
inline int parse_argi(parser& par, const std::string& longname,
    const std::string& help, int def, bool required = false,
    const std::vector<int>& choices = {}) {
    return parse_arg<int>(par, longname, help, def, required, choices);
}

///
/// Specialization of parse_arg()
///
inline float parse_argf(parser& par, const std::string& longname,
    const std::string& help, float def, bool required = false,
    const std::vector<float>& choices = {}) {
    return parse_arg<float>(par, longname, help, def, required, choices);
}

///
/// Specialization of parse_arg()
///
inline double parse_argd(parser& par, const std::string& longname,
    const std::string& help, double def, bool required = false,
    const std::vector<double>& choices = {}) {
    return parse_arg<double>(par, longname, help, def, required, choices);
}

///
/// Parses an argument array as described in the intro.
///
template <typename T>
inline std::vector<T> parse_arga(parser& par, const std::string& longname,
    const std::string& help, const std::vector<T>& def, int nargs = -1,
    bool required = false, const std::vector<T>& choices = {});

///
/// Specialization of parse_arga()
///
inline std::vector<std::string> parse_argas(parser& par,
    const std::string& longname, const std::string& help,
    const std::vector<std::string>& def, int nargs = -1, bool required = false,
    const std::vector<std::string>& choices = {}) {
    return parse_arga<std::string>(
        par, longname, help, def, nargs, required, choices);
}

}  // namespace cmdline

///
/// File loading and saving
///
namespace file {

///
/// Loads the contents of a binary file in an in-memory array.
///
inline std::vector<unsigned char> load_binfile(const std::string& filename) {
    auto file = fopen(filename.c_str(), "rb");
    if (!file) return {};
    auto ret = std::vector<unsigned char>();
    char buf[4096];
    int bufn;
    while ((bufn = (int)fread(buf, 1, sizeof(buf), file))) {
        ret.insert(ret.end(), buf, buf + bufn);
    }
    fclose(file);
    return ret;
}

///
/// Loads the contents of a text file into a std::string.
///
inline std::string load_txtfile(const std::string& filename) {
    auto file = fopen(filename.c_str(), "rt");
    if (!file) return "";
    auto ret = std::string();
    char buf[4096];
    int bufn;
    while ((bufn = (int)fread(buf, 1, sizeof(buf) - 1, file))) {
        buf[bufn] = 0;
        ret += buf;
    }
    fclose(file);
    return ret;
}

///
/// Saves binary data to a file.
///
inline bool save_binfile(
    const std::string& filename, const std::vector<unsigned char>& data) {
    auto f = fopen(filename.c_str(), "wt");
    if (!f) return false;
    fwrite(data.data(), 1, data.size(), f);
    fclose(f);
    return true;
}

///
/// Saves a string to a text file.
///
inline bool save_txtfile(const std::string& filename, const std::string& str) {
    auto f = fopen(filename.c_str(), "wt");
    if (!f) return false;
    fwrite(str.c_str(), 1, str.length(), f);
    fclose(f);
    return true;
}

}  // namespace file

///
/// Path manipulation
///
namespace path {

///
/// Get directory name (including '/').
///
inline std::string get_dirname(const std::string& filename) {
    auto pos = filename.rfind('/');
    if (pos == std::string::npos) pos = filename.rfind('\\');
    if (pos == std::string::npos) return "";
    return filename.substr(0, pos + 1);
}

///
/// Get extension (including '.').
///
inline std::string get_extension(const std::string& filename) {
    auto pos = filename.rfind('.');
    if (pos == std::string::npos) return "";
    return filename.substr(pos);
}

///
/// Get file basename.
///
inline std::string get_basename(const std::string& filename) {
    auto dirname = get_dirname(filename);
    auto extension = get_extension(filename);
    return filename.substr(
        dirname.size(), filename.size() - dirname.size() - extension.size());
}

///
/// Get filename without directory (equiv to get_basename() +
/// get_extension()).
///
inline std::string get_filename(const std::string& filename) {
    return get_basename(filename) + get_extension(filename);
}

///
/// Replace extension.
///
inline std::string replace_extension(
    const std::string& filename, const std::string& ext) {
    return get_dirname(filename) + get_basename(filename) + ext;
}

///
/// Prepend a string to the extension.
///
inline std::string prepend_extension(
    const std::string& filename, const std::string& prep) {
    return get_dirname(filename) + get_basename(filename) + prep +
           get_extension(filename);
}

///
/// Splits a path calling the above functions.
///
inline void split_path(const std::string& filename, std::string& dirname,
    std::string& basename, std::string& ext) {
    dirname = get_dirname(filename);
    basename = get_basename(filename);
    ext = get_extension(filename);
}

}  // namespace path

///
/// String manipulation
///
namespace string {

///
/// Checks if a std::string starts with a prefix.
///
inline bool startswith(const std::string& str, const std::string& substr) {
    if (str.length() < substr.length()) return false;
    for (auto i = 0; i < substr.length(); i++)
        if (str[i] != substr[i]) return false;
    return true;
}

///
/// Checks if a std::string ends with a prefix.
///
inline bool endswith(const std::string& str, const std::string& substr) {
    if (str.length() < substr.length()) return false;
    auto offset = str.length() - substr.length();
    for (auto i = 0; i < substr.length(); i++)
        if (str[i + offset] != substr[i]) return false;
    return true;
}

///
/// Check is a string contains a substring.
///
inline bool contains(const std::string& str, const std::string& substr) {
    return str.find(substr) != str.npos;
}

///
/// Splits a std::string into lines at the '\n' character. The line
/// terminator is kept if keep_newline. This function does not work on
/// Window if keep_newline is true.
///
inline std::vector<std::string> splitlines(
    const std::string& str, bool keep_newline = false) {
    if (str.empty()) return {};
    auto lines = std::vector<std::string>();
    auto line = std::vector<char>();
    for (auto c : str) {
        if (c == '\n') {
            if (keep_newline) line.push_back(c);
            lines.push_back(std::string(line.begin(), line.end()));
            line.clear();
        } else {
            line.push_back(c);
        }
    }
    if (!line.empty()) lines.push_back(std::string(line.begin(), line.end()));
    return lines;
}

///
/// Partition the string.
///
inline std::vector<std::string> partition(
    const std::string& str, const std::string& split) {
    auto pos = str.find(split);
    if (pos == str.npos) return {str, "", ""};
    return {str.substr(0, pos), split, str.substr(pos + split.length())};
}

///
/// Splits the string.
///
inline std::vector<std::string> split(const std::string& str) {
    auto ret = std::vector<std::string>();
    ret.push_back("");
    for (auto c : str) {
        if (c == ' ' || c == '\t' || c == '\r' || c == '\n') {
            if (!ret.back().empty()) ret.push_back("");
            continue;
        }
        ret.back() += c;
    }
    if (ret.back().empty()) ret.pop_back();
    return ret;
}

///
/// Strip the string.
///
inline std::string rstrip(const std::string& str) {
    auto pos = str.find_last_not_of(" \t\r\n");
    if (pos == str.npos) return "";
    return str.substr(0, pos + 1);
}

///
/// Strip the string.
///
inline std::string lstrip(const std::string& str) {
    auto pos = str.find_first_not_of(" \t\r\n");
    if (pos == str.npos) return "";
    return str.substr(pos);
}

///
/// Strip the string.
///
inline std::string strip(const std::string& str) { return rstrip(lstrip(str)); }

///
/// Joins a list of std::string with a std::string as separator.
///
inline std::string join(
    const std::vector<std::string>& strs, const std::string& sep) {
    auto ret = std::string();
    auto first = true;
    for (auto& str : strs) {
        if (!first) ret += sep;
        ret += str;
        first = false;
    }
    return ret;
}

///
/// Converts an ASCII string to lowercase.
///
inline std::string lower(const std::string& str) {
    auto s = str;
    for (auto& c : s) c = std::tolower(c);
    return s;
}

///
/// Converts an ASCII string to uppercase.
///
inline std::string upper(const std::string& str) {
    auto s = str;
    for (auto& c : s) c = std::toupper(c);
    return s;
}

///
/// Check if a string is space.
///
inline bool isspace(const std::string& str) {
    for (auto c : str) {
        if (c != ' ' && c != '\n' && c != '\t' && c != '\r') return false;
    }
    return true;
}

///
/// Replace s1 with s2 in str.
///
inline std::string replace(
    const std::string& str, const std::string& s1, const std::string& s2) {
    auto s = std::string();
    auto last = 0;
    auto pos = (int)str.find(s1);
    while (pos != str.npos) {
        s += str.substr(last, pos - last);
        s += s2;
        last = pos + (int)s1.length();
        pos = (int)str.find(s1, last);
    }
    s += str.substr(last);
    return s;
}

//
// Argument conversion for format
//
inline int format_arg(int v) { return v; }
inline float format_arg(float v) { return v; }
inline double format_arg(double v) { return v; }
inline const char* format_arg(bool v) { return (v) ? "true" : "false"; }
inline const char* format_arg(const char* v) { return v; }
inline const char* format_arg(const std::string& v) { return v.c_str(); }

///
/// C-like string formatting. This is only meant for short strings with max
/// length 10000 chars. Memory corruption will happen for longer strings.
///
template <typename... Args>
inline std::string formatf(const std::string& fmt, const Args&... args) {
    char buffer[1024 * 16];
    sprintf(buffer, fmt.c_str(), format_arg(args)...);
    return buffer;
}

}  // namespace string

///
/// Python-like STL functions
///
namespace containers {

///
/// Checks if a containers contains a value
///
template <typename T>
inline bool contains(const std::vector<T>& v, const T& vv) {
    return std::find(v.begin(), v.end(), vv) != v.end();
}

///
/// Checks if a containers contains a value
///
template <typename K, typename V>
inline bool contains(const std::map<K, V>& v, const K& vv) {
    return v.find(vv) != v.end();
}

///
/// Checks if a containers contains a value
///
template <typename K, typename V>
inline bool contains(const std::unordered_map<K, V>& v, const K& vv) {
    return v.find(vv) != v.end();
}

}  // namespace containers

///
/// Python-like STL operators
///
namespace operators {

///
/// Append an element to a vector
///
template <typename T>
inline std::vector<T> operator+(const std::vector<T>& v, const T& vv) {
    auto vc = std::vector<T>();
    vc.reserve(v.size() + 1);
    vc.insert(vc.end(), v.begin(), v.end());
    vc.push_back(vv);
    return vc;
}

///
/// Append an element to a vector
///
template <typename T>
inline std::vector<T>& operator+=(std::vector<T>& v, const T& vv) {
    v.push_back(vv);
    return v;
}

///
/// Append an element to a vector
///
template <typename T, typename ET>
inline std::vector<T> operator+(const std::vector<T>& v, const ET& vv) {
    auto vc = std::vector<T>();
    vc.reserve(v.size() + 1);
    vc.insert(vc.end(), v.begin(), v.end());
    vc.push_back(vv);
    return vc;
}

///
/// Append an element to a vector
///
template <typename T, typename ET>
inline std::vector<T>& operator+=(std::vector<T>& v, const ET& vv) {
    v.push_back(vv);
    return v;
}

///
/// Append a vector to a vector
///
template <typename T>
inline std::vector<T> operator+(
    const std::vector<T>& v, const std::vector<T>& vv) {
    auto vc = std::vector<T>();
    vc.reserve(v.size() + vv.size());
    vc.insert(vc.end(), v.begin(), v.end());
    vc.insert(vc.end(), vv.begin(), vv.end());
    return vc;
}

///
/// Append a vector to a vector
///
template <typename T>
inline std::vector<T>& operator+=(std::vector<T>& v, const std::vector<T>& vv) {
    v.insert(v.end(), vv.begin(), vv.end());
    return v;
}

}  // namespace operators

///
/// Simple Logging
///
namespace logging {

///
/// Logging level
///
enum struct log_level : int {
    /// verbose
    verbose = -3,
    /// trace
    trace = -2,
    /// debug
    debug = -1,
    /// info
    info = 0,
    /// warning
    warning = 1,
    /// error
    error = 2,
    /// fatal
    fatal = 3,
};

///
/// Logger object. A logger can output messages to multiple streams.
/// Use add streams commands for it.
///
struct logger;

///
/// Make a logger with an optional console stream.
///
inline logger make_logger(
    const std::string& name, bool add_console_stream = false);

///
/// Set logger default name
///
inline void set_logger_name(logger& lgr, const std::string& name);

///
/// Add a file stream to a logger.
///
/// - Parameters:
///     - lgr: logger
///     - filename: filename
///     - append: append or write open mode for file logger
///     - short_message: whether to use a short message version
///     - output_level: output level
///     - flush_level: output level
/// - Returns:
///     - true if ok
///
inline bool add_file_stream(logger& lgr, const std::string& filename,
    bool append, bool short_message = false,
    log_level output_level = log_level::info,
    log_level flush_level = log_level::info);

///
/// Add a console stream to a logger.
///
/// - Parameters:
///     - lgr: logger
///     - filename: logger filename or stderr if empty
///     - use_std_error: use standard error instead of standard out
///     - short_message: whether to use a short message version
///     - output_level: output level
///     - flush_level: output level
/// - Returns:
///     - true if ok
///
inline bool add_console_stream(logger& lgr, bool use_std_error = false,
    bool short_message = true, log_level output_level = log_level::info,
    log_level flush_level = log_level::info);

///
/// Get default logger.
/// By default a non-verbose stdout logger is creater.
///
inline logger& get_default_logger();

///
/// Set default logger name
///
inline void set_logger_name(const std::string& name) {
    set_logger_name(get_default_logger(), name);
}

///
/// Add a file logger to the default loggers.
///
inline void add_file_stream(const std::string& filename, bool append,
    bool short_message = false, log_level output_level = log_level::info,
    log_level flush_level = log_level::info) {
    add_file_stream(get_default_logger(), filename, append, short_message,
        output_level, flush_level);
}

///
/// Log a message
///
/// - Parameters:
///     - lgr: logger
///     - level: message level
///     - code: message code (5 chars)
///     - msg: message
///
inline void log_msg(logger& lgr, log_level level, const std::string& name,
    const std::string& msg);

///
/// Log a message formatted ala printf.
///
/// - Parameters:
///     - lgr: logger
///     - level: message level
///     - code: message code (5 chars)
///     - msg: message
///
template <typename... Args>
inline void log_msg(logger& lgr, log_level level, const std::string& name,
    const std::string& msg, const Args&... args) {
    log_msg(lgr, level, name, string::formatf(msg, args...));
}

///
/// Logs a message to the default loggers
///
template <typename... Args>
inline void log_msg(log_level level, const char* name, const std::string& msg,
    const Args&... args) {
    log_msg(get_default_logger(), level, name, msg, args...);
}

///
/// Logs a message to the default loggers
///
template <typename... Args>
inline void log_info(const std::string& msg, const Args&... args) {
    log_msg(log_level::info, "", msg, args...);
}

///
/// Logs a message to the default loggers
///
template <typename... Args>
inline void log_error(const std::string& msg, const Args&... args) {
    log_msg(log_level::fatal, "", msg, args...);
}

///
/// Logs a message to the default loggers
///
template <typename... Args>
inline void log_fatal(const std::string& msg, const Args&... args) {
    log_msg(log_level::fatal, "", msg, args...);
}

///
/// Timer for logging
///
struct log_timer {};

///
/// Log a message and start a timer
///
template <typename... Args>
inline log_timer log_timed(const std::string& msg, const Args&... args) {}

}  // namespace logging

///
/// Concurrent  and async execution based on a thread pool
///
namespace concurrent {

///
/// Forward declaration of thread pool.
///
struct thread_pool;

///
/// Initialize a thread pool with a certain number of threads (0 for
/// defatul).
///
inline thread_pool* make_pool(int nthread = 0);

///
/// Free the thread pool
///
inline void free_pool(thread_pool*& pool);

///
/// Wait for all jobs to finish
///
inline void wait_pool(thread_pool* pool);

///
/// Clear all jobs
///
inline void clear_pool(thread_pool* pool);

///
/// Parallel for implementation
///
inline void parallel_for(
    thread_pool* pool, int count, const std::function<void(int idx)>& task);

///
/// Runs a task asynchronously onto a thread pool
///
inline std::shared_future<void> run_async(
    thread_pool* pool, const std::function<void()>& task);

///
/// Wait for all jobs to finish on a global thread pool
///
inline void wait_pool();

///
/// Clear all jobs on a global thread pool
///
inline void clear_pool();

///
/// Runs a task asynchronously onto a global thread pool
///
inline std::shared_future<void> run_async(const std::function<void()>& task);

///
/// Parallel for implementation on a global thread pool
///
inline void parallel_for(int count, const std::function<void(int idx)>& task);

}  // namespace concurrent

// -----------------------------------------------------------------------------
// TIMER
// -----------------------------------------------------------------------------

///
/// Simple timer for performance measumrents.
///
namespace timer {

///
/// A simple wrapper for std::chrono.
///
struct timer {
    /// initialize a timer and start it if necessary
    timer(bool autostart = true) {
        if (autostart) start();
    }

    /// start a timer
    void start() {
        _start = std::chrono::steady_clock::now();
        _started = true;
    }

    /// stops a timer
    void stop() {
        _end = std::chrono::steady_clock::now();
        _started = false;
    }

    /// elapsed time
    double elapsed() {
        if (_started) stop();
        std::chrono::duration<double> diff = (_end - _start);
        return diff.count();
    }

   private:
    bool _started = false;
    std::chrono::time_point<std::chrono::steady_clock> _start, _end;
};

}  // namespace timer

}  // namespace yu

// -----------------------------------------------------------------------------
// PRIVATE IMPLEMENTATION
// -----------------------------------------------------------------------------

//
// This code and data structure are meant to be porivate and will change at any
// moment
//
namespace yu {

namespace cmdline {

//
// Command-line parsing state
//
struct parser {
    // saved arguments ----------------------------------------
    std::vector<std::string> args;

    // help std::strings -------------------------------------------
    std::string help_prog;   // program description
    std::string help_usage;  // usage lines
    std::string help_opts;   // options lines
    std::string help_args;   // unnamed arguments lines

    // parse data ---------------------------------------------
    bool error;             // whether a parsing error occurred
    std::string error_msg;  // error message
    bool exit_on_error;     // exit on error
    bool print_help;        // print help
};

//
// Inits the parser.
//
inline parser make_parser(const std::vector<std::string>& args,
    const std::string& name, const std::string& help) {
    // clears parser and copy argument data
    auto par = parser();
    par.args = std::vector<std::string>(args.begin() + 1, args.end());
    par.error = false;
    par.exit_on_error = true;
    par.help_prog = (name == "") ? args[0] : name;
    par.help_usage += help;
    par.print_help = parse_flag(par, "--help", "-?", "print help", false);
    return par;
}

//
// Inits the parser.
//
inline parser make_parser(
    int argc, char* argv[], const std::string& name, const std::string& help) {
    return make_parser(std::vector<std::string>(argv, argv + argc), name, help);
}

//
// Print help based on the help lines collected during parsing.
//
static inline void _print_help(parser& par) {
    auto help = std::string();
    help += "usage: " + par.help_prog;
    if (!par.help_opts.empty()) help += "[options] ";
    if (!par.help_args.empty()) help += "<arguments> ";
    help += "\n    " + par.help_usage + "\n\n";
    if (!par.help_opts.empty()) {
        help += "options:\n";
        help += par.help_opts;
    }
    if (!par.help_args.empty()) {
        help += "arguments:\n";
        help += par.help_args;
    }
    printf("%s\n", help.c_str());
}

//
// Ends parsing checking for error for unused options or arguments.
// Exit if needed.
//
inline bool check_parser(parser& par) {
    // check for error
    if (!par.error && par.args.size() > 0) {
        par.error = true;
        if (par.args[0][0] == '-')
            par.error_msg = "unknown option " + par.args[0];
        else
            par.error_msg = "unsued values";
    };
    // check whether we need to print help and exit
    if (par.error) printf("error: %s\n", par.error_msg.c_str());
    if (par.print_help || par.error) {
        _print_help(par);
        if (par.exit_on_error) exit(EXIT_FAILURE);
    }
    auto ok = !par.error;
    return ok;
}

//
// Check if a std::string starts with another.
//
static inline bool _startswith(
    const std::string& str, const std::string& start) {
    if (str.length() < start.length()) return false;
    for (auto i = 0; i < start.length(); i++) {
        if (str[i] != start[i]) return false;
    }
    return true;
}

//
// Check if an option name is valid
//
static inline void _check_name(parser& par, const std::string& longname,
    const std::string& shortname, bool opt, int nargs) {
    // check name
    assert(!longname.empty());
    if (opt) {
        if (!shortname.empty())
            assert(_startswith(shortname, "-") && shortname.length() > 1);
        assert(_startswith(longname, "--") && longname.length() > 2);
    } else {
        assert(shortname.empty());
        assert(longname[0] != '-');
    }
    assert((opt && nargs >= 0) || (!opt && (nargs == -1 || nargs > 0)));
}

//
// Convert a type to std::string
//
static inline std::string _typename(bool) { return "bool"; }
static inline std::string _typename(int) { return "int"; }
static inline std::string _typename(float) { return "float"; }
static inline std::string _typename(double) { return "double"; }
static inline std::string _typename(const std::string&) { return "string"; }
// template <typename T>
// static inline std::string _typename(const std::vector<T>&) {
//     return "";
// }

//
// Converts a value to a std::string
//
template <typename T>
static inline std::string _tostring(const T& val) {
    std::ostringstream stream;
    stream << val;
    return stream.str();
}
template <>
inline std::string _tostring<bool>(const bool& val) {
    return (val) ? "true" : "false";
}

//
// Add a formatted help line
//
template <typename T>
static inline void _add_help(parser& par, const std::string& longname,
    const std::string& shortname, const std::string& help, bool opt, bool req,
    int nargs, const std::vector<T>& def, const std::vector<T>& choices = {}) {
    // dummy variable for function overload
    T _dummy = {};

    // full name
    auto help_fullname = longname;
    if (!shortname.empty()) help_fullname += "/" + shortname;
    if (nargs != 0) help_fullname += " " + _typename(_dummy);

    // default
    auto help_def = std::string();
    if (!req) {
        if (nargs == 0 || nargs == 1) {
            if (!def.empty())
                help_def = "[" + _tostring((const T&)def[0]) + "]";
        } else if (nargs == -1 || nargs > 1) {
            help_def = "[";
            for (auto i = 0; i < def.size(); i++) {
                if (i) help_def += ",";
                help_def += _tostring((const T&)def[i]);
            }
            help_def += "]";
        } else {
            // left empty
        }
    }

    // choices
    auto help_choice = std::string();
    if (!choices.empty()) {
        help_choice += "(";
        for (auto i = 0; i < choices.size(); i++) {
            if (i) help_choice += ",";
            help_choice += _tostring((const T&)choices[i]);
        }
        help_choice += ")";
    }

    // print help line
    char buf[10000] = {0};
    sprintf(buf, "  %-24s  %s %s\n", help_fullname.c_str(), help.c_str(),
        help_def.c_str());
    auto help_line = std::string(buf);
    if (!help_choice.empty()) {
        sprintf(buf, "  %-24s  %s\n", "", help_choice.c_str());
        help_line += buf;
    }

    // add line to proper help
    if (opt)
        par.help_opts += help_line;
    else
        par.help_args += help_line;
}

//
// Parsing routine for arrays of values
//
template <typename T>
static inline std::vector<T> _parse_vals(parser& par,
    const std::string& longname, const std::string& shortname,
    const std::string& help, bool opt, bool req, int nargs,
    const std::vector<T>& def, const std::vector<T>& choices = {}) {
    // prepare default empty vec
    auto vals = std::vector<T>();

    // check whether the name is good
    _check_name(par, longname, shortname, opt, nargs);

    // add help
    _add_help(par, longname, shortname, help, opt, req, nargs, def, choices);

    // skip if alreasy in error
    if (par.error) return vals;

    // find the value position
    auto val_pos = -1;
    if (opt) {
        // find option name
        for (auto i = 0; i < par.args.size() && val_pos < 0; i++) {
            if (shortname == par.args[i]) val_pos = i;
            if (longname == par.args[i]) val_pos = i;
        }

        // remove the option name
        if (val_pos >= 0) { par.args.erase(par.args.begin() + val_pos); }
    } else {
        // check if arg is present
        if (!par.args.empty()) {
            if (par.args[0][0] != '-') { val_pos = 0; }
        }
    }

    // handle not found
    if (val_pos < 0) {
        if (req) {
            par.error = true;
            par.error_msg = "missing value for " + longname;
            return vals;
        } else
            return def;
    }

    // check if value is present
    if (nargs == -1) nargs = std::max(1, (int)par.args.size());
    if (val_pos + nargs > par.args.size()) {
        par.error = true;
        par.error_msg = "missing value for " + longname;
        return vals;
    }

    // loop over values
    for (auto i = 0; i < nargs; i++) {
        // grab value
        auto val_str = par.args[val_pos];
        par.args.erase(par.args.begin() + val_pos);

        // parse value
        auto stream = std::istringstream(val_str);
        auto val = T();
        if (!(stream >> val)) {
            par.error = true;
            par.error_msg = "incorrect value for " + longname;
            return vals;
        }
        // check choices
        if (!choices.empty()) {
            auto in_choices = false;
            for (auto&& c : choices) {
                if (val == c) in_choices = true;
            }
            if (!in_choices) {
                par.error = true;
                par.error_msg = "incorrect value for " + longname;
                return vals;
            }
        }

        // add
        vals.push_back(val);
    }

    // done
    return vals;
}

//
// Parsing routine for values
//
template <typename T>
static inline T _parse_val(parser& par, const std::string& longname,
    const std::string& shortname, const std::string& help, bool opt,
    const T& def, bool req, const std::vector<T>& choices = {}) {
    // parse values
    auto vals = _parse_vals(
        par, longname, shortname, help, opt, req, 1, {def}, choices);

    // return value if present
    if (vals.size())
        return vals[0];
    else
        return {};
}

//
// Parses an optional flag as described in the intro.
//
inline bool parse_flag(parser& par, const std::string& longname,
    const std::string& shortname, const std::string& help, bool def) {
    // parse values
    auto vals = _parse_vals<bool>(
        par, longname, shortname, help, true, false, 0, {def});

    // return value if present
    if (vals.size())
        return def;
    else
        return !def;
}

//
// Parses an option as described in the intro.
//
template <typename T>
inline T parse_opt(parser& par, const std::string& longname,
    const std::string& shortname, const std::string& help, const T& def,
    bool req, const std::vector<T>& choices) {
    return _parse_val<T>(
        par, longname, shortname, help, true, def, req, choices);
}

//
// Parses an argument as described in the intro.
//
template <typename T>
inline T parse_arg(parser& par, const std::string& longname,
    const std::string& help, const T& def, bool req,
    const std::vector<T>& choices) {
    return _parse_val<T>(par, longname, "", help, false, def, req, choices);
}

//
// Parses an option array as described in the intro.
//
template <typename T>
inline std::vector<T> parse_opta(parser& par, const std::string& longname,
    const std::string& help, const std::vector<T>& def, int nargs,
    bool required, const std::vector<T>& choices) {
    return _parse_vals(
        par, longname, "", help, true, required, nargs, def, choices);
}

//
// Parses an argument array as described in the intro.
//
template <typename T>
inline std::vector<T> parse_arga(parser& par, const std::string& longname,
    const std::string& help, const std::vector<T>& def, int nargs,
    bool required, const std::vector<T>& choices) {
    return _parse_vals(
        par, longname, "", help, false, required, nargs, def, choices);
}

//
// Parses an option enum as described in the intro.
//
inline int parse_opte(parser& par, const std::string& longname,
    const std::string& shortname, const std::string& help, int def,
    const std::vector<std::pair<std::string, int>>& vals, bool required) {
    auto choices = std::vector<std::string>();
    auto def_s = std::string();
    for (auto&& kv : vals) {
        choices.push_back(kv.first);
        if (kv.second == def) def_s = kv.first;
    }
    auto val_s =
        parse_opt(par, longname, shortname, help, def_s, required, choices);
    for (auto&& kv : vals) {
        if (kv.first == val_s) return kv.second;
    }
    assert(false);
    return -1;
}

}  // namespace cmdline

namespace logging {

//
// Logger Stream
//
struct logger_stream {
    // stream for the logger
    FILE* file = nullptr;
    // whether the logger uses a short message
    bool short_message = false;
    // log level for output filtering
    log_level output_level = log_level::info;
    // log level for flush filtering
    log_level flush_level = log_level::error;
    // log id generated randomly
    unsigned int guid = std::random_device()();

    // cleanup
    ~logger_stream() {
        if (file == stderr) return;
        if (file == stdout) return;
        if (file) {
            fclose(file);
            file = nullptr;
        }
    }
};

//
// Logger implemented as a collection of filtered stream
//
struct logger {
    /// name
    std::string name;
    /// streams
    std::vector<logger_stream> streams;
};

//
// Make a logger with an optional console stream.
//
inline logger make_logger(const std::string& name, bool add_stream) {
    auto lgr = logger();
    lgr.name = name;
    if (add_stream) add_console_stream(lgr);
    return lgr;
}

//
// set logger name
//
inline void set_logger_name(logger& lgr, const std::string& name) {
    lgr.name = name;
}

//
// default logger
//
static logger default_logger = make_logger("<app>", true);

//
// default logger
//
logger& get_default_logger() { return default_logger; }

//
// Create a file logger
//
inline bool add_file_stream(logger& lgr, const std::string& filename,
    bool append, bool short_message, log_level output_level,
    log_level flush_level) {
    auto file = fopen(filename.c_str(), (append) ? "at" : "wt");
    if (!file) return false;
    lgr.streams.push_back(logger_stream());
    lgr.streams.back().file = file;
    lgr.streams.back().short_message = short_message;
    lgr.streams.back().output_level = output_level;
    lgr.streams.back().flush_level = flush_level;
    return true;
}

//
// Create a stderr logger
//
inline bool add_console_stream(logger& lgr, bool use_std_error,
    bool short_message, log_level output_level, log_level flush_level) {
    lgr.streams.push_back(logger_stream());
    lgr.streams.back().file = (use_std_error) ? stderr : stdout;
    lgr.streams.back().short_message = short_message;
    lgr.streams.back().output_level = output_level;
    lgr.streams.back().flush_level = flush_level;
    return true;
}

//
// Log a message
//
inline void log_msg(logger& lgr, log_level level, const std::string& name,
    const std::string& msg) {
    // type string
    static const char* types[] = {"VERB", "INFO", "WARN", "ERRN"};
    const char* type = types[std::max(0, std::min(3, (int)level + 1))];

    for (auto&& stream : lgr.streams) {
        if (level < stream.output_level) continue;
        if (stream.short_message) {
            // time string
            char time_buf[1024];
            auto tm = time(nullptr);
            auto ttm = localtime(&tm);  // TODO: use thread safe version
            strftime(time_buf, 1024, "%H:%M:%S", ttm);

            // output message
            fprintf(stream.file, "%s %s %s\n", time_buf, type, msg.c_str());
        } else {
            // time string
            char time_buf[1024];
            auto tm = time(nullptr);
            auto ttm = localtime(&tm);  // TODO: use thread safe version
            strftime(time_buf, 1024, "%Y-%m-%d %H:%M:%S", ttm);

            // name
            auto rname = (name == "") ? lgr.name : name;

            // output message
            fprintf(stream.file, "%s %s %4x %-16s %s\n", time_buf, type,
                stream.guid, rname.c_str(), msg.c_str());
        }

        // flush if needed
        if (level < stream.flush_level) return;
        fflush(stream.file);
    }

    // exit if needed
    if (level >= log_level::fatal) exit(EXIT_FAILURE);
}

}  // namespace logging

namespace concurrent {

//
// Thread pool code derived from LLVM codebase
//

// thread pool
struct ThreadPool {
    // thread pool initialized with nthreads
    ThreadPool(int nthreads = std::thread::hardware_concurrency())
        : working_threads(0), stop_flag(false) {
        threads.reserve(nthreads);
        for (auto tid = 0; tid < nthreads; tid++) {
            threads.emplace_back([this] { thread_proc(); });
        }
    }

    // destructor
    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock_guard(queue_lock);
            stop_flag = true;
        }
        queue_condition.notify_all();
        for (auto& Worker : threads) Worker.join();
    }

    // empty the queue
    void clear() {
        {
            std::unique_lock<std::mutex> lock_guard(queue_lock);
            tasks.clear();
        }
        queue_condition.notify_all();
    }

    // schedule an asynchronous taks
    std::shared_future<void> async(std::function<void()> task) {
        // Wrap the Task in a packaged_task to return a future object.
        std::packaged_task<void()> packaged_task(std::move(task));
        auto future = packaged_task.get_future();
        {
            std::unique_lock<std::mutex> lock_guard(queue_lock);
            assert(
                !stop_flag && "Queuing a thread during ThreadPool destruction");
            tasks.push_back(std::move(packaged_task));
        }
        queue_condition.notify_one();
        return future.share();
    }

    // wait for all tasks to finish
    void wait() {
        std::unique_lock<std::mutex> lock_guard(completion_lock);
        completion_condition.wait(
            lock_guard, [&] { return tasks.empty() && !working_threads; });
    }

    // implementation -------------------------------------------------
   private:
    void thread_proc() {
        while (true) {
            std::packaged_task<void()> task;
            {
                std::unique_lock<std::mutex> lock_guard(queue_lock);
                queue_condition.wait(
                    lock_guard, [&] { return stop_flag || !tasks.empty(); });

                if (stop_flag && tasks.empty()) return;

                {
                    working_threads++;
                    std::unique_lock<std::mutex> lock_guard(completion_lock);
                }
                task = std::move(tasks.front());
                tasks.pop_front();
            }

            task();

            {
                std::unique_lock<std::mutex> lock_guard(completion_lock);
                working_threads--;
            }

            completion_condition.notify_all();
        }
    }

    std::vector<std::thread> threads;
    std::deque<std::packaged_task<void()>> tasks;
    std::mutex queue_lock;
    std::condition_variable queue_condition;
    std::mutex completion_lock;
    std::condition_variable completion_condition;
    std::atomic<unsigned> working_threads;
    bool stop_flag = false;
};

//
// End of thread pool code derived from LLVM codebase
//

//
// Forward declaration of thread pool.
//
struct thread_pool {
    ThreadPool* tp = nullptr;
    ~thread_pool() {
        if (tp) delete tp;
    }
};

//
// Initialize a thread pool with a certain number of threads (0 for defatul).
//
inline thread_pool* make_pool(int nthread) {
    if (!nthread) nthread = std::thread::hardware_concurrency();
    auto pool = new thread_pool();
    pool->tp = new ThreadPool(nthread);
    return pool;
}

//
// Clear thread pool
//
inline void free_pool(thread_pool*& pool) {
    if (pool) {
        clear_pool(pool);
        delete pool;
    }
    pool = nullptr;
}

//
// Enqueue a job
//
inline std::shared_future<void> run_async(
    thread_pool* pool, const std::function<void()>& task) {
    return pool->tp->async(task);
}

//
// Wait for jobs to finish
//
inline void wait_pool(thread_pool* pool) { pool->tp->wait(); }

//
// Wait for jobs to finish
//
inline void clear_pool(thread_pool* pool) { pool->tp->clear(); }

//
// Parallel for implementation
//
inline void parallel_for(
    thread_pool* pool, int count, const std::function<void(int idx)>& task) {
    for (auto idx = 0; idx < count; idx++) {
        run_async(pool, [&task, idx]() { task(idx); });
    }
    wait_pool(pool);
}

//
// Global pool
//
static auto global_pool = (thread_pool*)nullptr;

//
// Make the global thread pool
//
static inline void make_global_thread_pool() {
    if (!global_pool) global_pool = make_pool();
}

//
// Enqueue a job
//
inline std::shared_future<void> run_async(const std::function<void()>& task) {
    if (!global_pool) make_global_thread_pool();
    return run_async(global_pool, task);
}

//
// Wait for jobs to finish
//
inline void wait_pool() {
    if (global_pool) global_pool->tp->wait();
}

//
// Wait for jobs to finish
//
inline void clear_pool() {
    if (global_pool) clear_pool(global_pool);
}

//
// Parallel for implementation
//
inline void parallel_for(int count, const std::function<void(int idx)>& task) {
    if (!global_pool) make_global_thread_pool();
    parallel_for(global_pool, count, task);
}

}  // namespace concurrent

}  // namespace yu

#ifndef _WIN32
#pragma GCC diagnostic pop
#endif

#ifdef YU_TEST
//
// Test
//
void run_test() {
    using namespace yu::cmdline;

    // test empty
    auto test0_argc = 1;
    const char* test0_argv[] = {"test0"};
    auto par0 = make_parser(test0_argc, (char**)test0_argv, "test0");
    par0->exit_on_error = false;
    assert(check_parser(par0) == true);

    // test exit on help
    auto test1_argc = 2;
    const char* test1_argv[] = {"test1", "--help"};
    auto par1 = make_parser(test1_argc, (char**)test1_argv, "test1");
    par1->exit_on_error = false;
    assert(check_parser(par1) == true);

    // test opts
    auto test2_argc = 10;
    const char* test2_argv[] = {"test2", "--int", "10", "--float", "3.14",
        "--double", "6.28", "--str", "bob", "--flag"};
    auto par2 = make_parser(test2_argc, (char**)test2_argv, "test2");
    par2->exit_on_error = false;
    assert(parse_flag(par2, "--flag", "", "", false) == true);
    assert(parse_opt<int>(par2, "--int", "", "", 0) == 10);
    assert(fabsf(parse_opt<float>(par2, "--float", "", "", 0) - 3.14f) < 0.01);
    assert(fabs(parse_opt<double>(par2, "--double", "", "", 0) - 6.28) < 0.01);
    assert(parse_opt<std::string>(par2, "--str", "", "", "mike") == "bob");
    assert(parse_flag(par2, "--flag_def", "", "", false) == false);
    assert(parse_opt<int>(par2, "--int_def", "", "", 5) == 5);
    assert(parse_opt<float>(par2, "--float_def", "", "", 2.67f) == 2.67f);
    assert(parse_opt<double>(par2, "--double_def", "", "", 9.54) == 9.54);
    assert(parse_opt<std::string>(par2, "--str_def", "", "", "alex") == "alex");
    assert(check_parser(par2) == true);

    // test args
    auto test3_argc = 3;
    const char* test3_argv[] = {"test3", "10", "bob"};
    auto par3 = make_parser(test3_argc, (char**)test3_argv, "test3");
    par3->exit_on_error = false;
    assert(parse_arg<int>(par3, "int", "", 0, true) == 10);
    assert(parse_arg<std::string>(par3, "str", "", "mike", true) == "bob");
    assert(check_parser(par3) == true);

    // test bad opts
    auto test4_argc = 3;
    const char* test4_argv[] = {"test4", "--int", "bob"};
    auto par4 = make_parser(test4_argc, (char**)test4_argv, "test4");
    par4->exit_on_error = false;
    assert(parse_opt<int>(par4, "--int", "", "", 0) == 0);
    assert(check_parser(par4) == false);
}
#endif

#endif
