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
///
///
/// ## History
///
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

#include <chrono>
#include <functional>
#include <future>
#include <string>
#include <unordered_map>
#include <vector>

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
parser* make_parser(
    const std::vector<std::string>& args, const std::string& help = "");

///
/// Inits a command line parser.
///
parser* make_parser(int argc, char* argv[], const char* help);

///
/// Ends parsing checking for error for unused options or arguments.
/// Exit if needed.
///
bool check_parser(parser* prs);

///
/// Parses an optional flag as described in the intro.
///
bool parse_flag(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, bool def = false);

///
/// Parses an option as described in the intro.
///
template <typename T>
T parse_opt(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, const T& def,
    bool required = false, const std::vector<T>& choices = {});

///
/// Specialization of parse_opt()
///
std::string parse_opts(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help,
    const std::string& def, bool required = false,
    const std::vector<std::string>& choices = {});

///
/// Specialization of parse_opt()
///
int parse_opti(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, int def,
    bool required = false, const std::vector<int>& choices = {});

///
/// Specialization of parse_opt()
///
float parse_optf(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, float def,
    bool required = false, const std::vector<float>& choices = {});

///
/// Specialization of parse_opt()
///
double parse_optd(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, double def,
    bool required = false, const std::vector<double>& choices = {});

///
/// Specialization of parse_opt()
///
bool parse_optb(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, bool def,
    bool required = false);

///
/// Parses an option enum as described in the intro.
///
int parse_opte(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, int def,
    const std::vector<std::pair<std::string, int>>& vals,
    bool required = false);

///
/// Parses an option array as described in the intro.
///
template <typename T>
T parse_opta(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, const T& def,
    int nargs, bool required = false, const std::vector<T>& choices = {});

///
/// Parses an argument as described in the intro.
///
template <typename T>
T parse_arg(parser* par, const std::string& longname, const std::string& help,
    const T& def, bool required = false, const std::vector<T>& choices = {});

///
/// Specialization of parse_arg()
///
std::string parse_args(parser* par, const std::string& longname,
    const std::string& help, const std::string& def, bool required = false,
    const std::vector<std::string>& choices = {});

///
/// Specialization of parse_arg()
///
int parse_argi(parser* par, const std::string& longname,
    const std::string& help, int def, bool required = false,
    const std::vector<int>& choices = {});

///
/// Specialization of parse_arg()
///
float parse_argf(parser* par, const std::string& longname,
    const std::string& help, float def, bool required = false,
    const std::vector<float>& choices = {});

///
/// Specialization of parse_arg()
///
double parse_argd(parser* par, const std::string& longname,
    const std::string& help, double def, bool required = false,
    const std::vector<double>& choices = {});

///
/// Parses an argument array as described in the intro.
///
template <typename T>
std::vector<T> parse_arga(parser* par, const std::string& longname,
    const std::string& help, const std::vector<T>& def, int nargs = -1,
    bool required = false, const std::vector<T>& choices = {});

///
/// Specialization of parse_arga()
///
std::vector<std::string> parse_argas(parser* par, const std::string& longname,
    const std::string& help, const std::vector<std::string>& def,
    int nargs = -1, bool required = false,
    const std::vector<std::string>& choices = {});

///
/// Specialization of parse_arga()
///
std::vector<int> parse_argai(parser* par, const std::string& longname,
    const std::string& help, const std::vector<int>& def, int nargs = -1,
    bool required = false, const std::vector<int>& choices = {});

///
/// Specialization of parse_arga()
///
std::vector<float> parse_argaf(parser* par, const std::string& longname,
    const std::string& help, const std::vector<float>& def, int nargs = -1,
    bool required = false, const std::vector<float>& choices = {});

///
/// Specialization of parse_arga()
///
std::vector<double> parse_argad(parser* par, const std::string& longname,
    const std::string& help, const std::vector<double>& def, int nargs = -1,
    bool required = false, const std::vector<double>& choices = {});

}  // namespace cmdline

///
/// File loading and saving
///
namespace file {

///
/// Loads the contents of a binary file in an in-memory array.
///
std::vector<unsigned char> load_binfile(const std::string& filename);

///
/// Loads the contents of a text file into a std::string.
///
std::string load_txtfile(const std::string& filename);

///
/// Saves binary data to a file.
///
void save_binfile(
    const std::string& filename, const std::vector<unsigned char>& data);

///
/// Saves a string to a text file.
///
void save_txtfile(const std::string& filename, const std::string& str);

}  // namespace file

///
/// Path manipulation
///
namespace path {

///
/// Get directory name (including '/').
///
std::string get_dirname(const std::string& filename);

///
/// Get file basename.
///
std::string get_basename(const std::string& filename);

///
/// Get extension (including '.').
///
std::string get_extension(const std::string& filename);

///
/// Get filename without directory (equiv to get_basename() +
/// get_extension()).
///
std::string get_filename(const std::string& filename);

///
/// Replace extension.
///
std::string replace_extension(
    const std::string& filename, const std::string& ext);

///
/// Prepend a string to the extension.
///
std::string prepend_extension(
    const std::string& filename, const std::string& prep);

///
/// Splits a path calling the above functions.
///
void split_path(const std::string& filename, std::string& dirname,
    std::string& basename, std::string& ext);

}  // namespace path

///
/// String manipulation
///
namespace string {

///
/// Checks if a std::string starts with a prefix.
///
bool startswith(const std::string& str, const std::string& substr);

///
/// Checks if a std::string ends with a prefix.
///
bool endswith(const std::string& str, const std::string& substr);

///
/// Check is a string contains a substring.
///
bool contains(const std::string& str, const std::string& substr);

///
/// Splits a std::string into lines at the '\n' character. The line
/// terminator is kept if keep_newline. This function does not work on
/// Window if keep_newline is true.
///
std::vector<std::string> splitlines(
    const std::string& str, bool keep_newline = false);

///
/// Partition the string.
///
std::vector<std::string> partition(
    const std::string& str, const std::string& split);

///
/// Splits the string.
///
std::vector<std::string> split(const std::string& str);

///
/// Strip the string.
///
std::string strip(const std::string& str);

///
/// Strip the string.
///
std::string rstrip(const std::string& str);

///
/// Strip the string.
///
std::string lstrip(const std::string& str);

///
/// Joins a list of std::string with a std::string as separator.
///
std::string join(const std::vector<std::string>& strs, const std::string& sep);

///
/// Converts an ASCII string to lowercase.
///
std::string lower(const std::string& str);

///
/// Converts an ASCII string to uppercase.
///
std::string upper(const std::string& str);

///
/// Strung is space.
///
bool isspace(const std::string& str);

///
/// Replace s1 with s2 in str.
///
std::string replace(
    const std::string& str, const std::string& s1, const std::string& s2);

///
/// C-like string formatting. This is only meant for short strings with max
/// length 10000 chars. Memory corruption will happen for longer strings.
///
std::string format(const char* fmt, va_list args);

///
/// C-like string formatting. See format_str(fmt,args);
///
std::string format(const char* fmt, ...);

}  // namespace string

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
/// Logger object
///
struct logger;

///
/// Make a file logger. See set_logger() for Parameters.
///
/// - Parameters:
///     - filename: logger filename or stderr if empty
///     - apend: append or write open mode for file logger
///
///
logger* make_file_logger(const std::string& filename, bool append,
    log_level output_level = log_level::info,
    log_level flush_level = log_level::info);

///
/// Make a stderr logger. See set_logger() for Parameters.
///
logger* make_stderr_logger(log_level output_level = log_level::info,
    log_level flush_level = log_level::info);

///
/// Make a stderr logger. See set_logger() for Parameters.
///
logger* make_stdout_logger(log_level output_level = log_level::info,
    log_level flush_level = log_level::info);

///
/// Clear a logger
///
void clear_logger(logger* lgr);

///
/// Get default loggers. This is a modifiable reference.
///
std::vector<logger*>* get_default_loggers();

///
/// Set logger level
///
/// - Parameters:
///     - lgr: logger
///     - name : logger name
///     - output_level: output level
///     - flush_level: output level
///
void set_logger(logger* lgr, log_level output_level,
    log_level flush_level = log_level::error);

///
/// Log a message
///
/// - Parameters:
///     - lgr: logger
///     - level: message level
///     - code: message code (5 chars)
///     - msg: message
///
void log_msg(logger* lgr, log_level level, const std::string& name,
    const std::string& msg);

///
/// Log a message
///
/// - Parameters:
///     - lgr: logger
///     - level: message level
///     - code: message code (5 chars)
///     - msg: message
///
void log_msg(logger* lgr, log_level level, const char* name, const char* msg);

///
/// Log a message formatted ala printf.
///
/// - Parameters:
///     - lgr: logger
///     - level: message level
///     - code: message code (5 chars)
///     - msg: message
///
void log_msgf(
    logger* lgr, log_level level, const char* name, const char* msg, ...);

///
/// Logs a message to the default loggers
///
void log_msgfv(logger* lgr, log_level level, const char* name, const char* msg,
    va_list args);

///
/// Logs a message to the default loggers
///
void log_msg(log_level level, const std::string& name, const std::string& msg);

///
/// Logs a message to the default loggers
///
void log_msgf(log_level level, const char* name, const char* msg, ...);

///
/// Logs a message to the default loggers
///
void log_msgfv(
    log_level level, const char* name, const char* msg, va_list args);

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
thread_pool* make_pool(int nthread = 0);

///
/// Clear thread pool
///
void clear_pool(thread_pool* pool);

///
/// Wait for all jobs to finish
///
void wait_pool(thread_pool* pool);

///
/// Parallel for implementation
///
void parallel_for(
    thread_pool* pool, int count, const std::function<void(int idx)>& task);

///
/// Runs a task asynchronously onto a thread pool
///
std::shared_future<void> run_async(
    thread_pool* pool, const std::function<void()>& task);

///
/// Wait for all jobs to finish on a global thread pool
///
void wait_pool();

///
/// Runs a task asynchronously onto a global thread pool
///
std::shared_future<void> run_async(const std::function<void()>& task);

///
/// Parallel for implementation on a global thread pool
///
void parallel_for(int count, const std::function<void(int idx)>& task);

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

#endif
