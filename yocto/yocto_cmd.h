///
/// YOCTO_CMD: utilities for writing command line applications, mostly
/// a simple to use command line parser and a few string, path and file
/// functions.
///
///
/// USAGE FOR COMMAND LINE PARSING:
///
/// 1. create a parser object
///     - an option for printing help is automatically added
///     parser = make_parser(argc, argv, program description)
/// 2. for each option, parse it calling the functions parse_opt()
///     - options are parsed on the fly and a comprehensive help is
///       automatically generated
///     - supports bool (flags), int, float, double, std::string, enums
///     - options names are "--longname" for longname and "-s" for short
///     - command line format is "--longname value", "-s v" for all but flags
///     - values are parsed with iostream << operators
///     - for general use opt = parse_opt<type>(...)
///     - for boolean flags is parse_flag()
///     - for enums use parse_opte(...)
/// 3. for each unnamed argument, parse it calling the functions parse_arg()
///     - names are only used for help
///     - supports types as above
///     - for general use arg = parse_arg<type>(...)
///     - to parse all remaining values use args = parse_arga<type>(...)
/// 4. end cmdline parsing with check_parser() to check for unsued values,
///    missing arguments and print help if needed
/// 5. since arguments are parsed immediately, one can easily implement
///    subcommands by just branching the command line code based on a read
///    argument without any need for complex syntax
///
/// Notes: the end of this file contains a test function that also
/// illustreates the library usage.
///
/// USAGE FOR COMMAND LINE PARSING:
///
/// 1. create loggers with make_file_logger(), make_stderr_logger() ,
///    make_stdout_logger()
/// 2. you can set default loggers with get_default_loggers(); note that none
///    are set by default
/// 3. write log messages with log_msg() and its variants.
///
/// UTILITIES
///
/// 1. filename splitting with get_dirname(), get_basename(), get_extension(),
///    split_path(), replace_extension(), prepend_extension()
/// 2. loading entire files with load_txtfile() and load_binfile()
/// 3. string manipulation with split_lines()
///
///
/// The interface for each function is described in details in the interface
/// section of this file.
///
///
/// COMPILATION:
///
/// To use the library include the .h and compile the .cpp. To use this library
/// as a header-only library, define YCMD_INLINE before including this file.
///
///
/// HISTORY:
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
namespace ycmd {}

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

#ifndef _YCMD_H_
#define _YCMD_H_

// compilation options
#ifdef YCMD_INLINE
#define YCMD_API inline
#else
#define YCMD_API
#endif

#include <functional>
#include <future>
#include <string>
#include <unordered_map>
#include <vector>

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

namespace ycmd {

// COMMAND LINE PARSING
// --------------------------------------------------------

///
/// Command line parser.
///
struct parser;

///
/// Inits a command line parser.
///
YCMD_API parser* make_parser(
    const std::vector<std::string>& args, const std::string& help = "");

///
/// Inits a command line parser.
///
YCMD_API parser* make_parser(int argc, char* argv[], const char* help);

///
/// Ends parsing checking for error for unused options or arguments.
/// Exit if needed.
///
YCMD_API bool check_parser(parser* prs);

///
/// Parses an optional flag as described in the intro.
///
YCMD_API bool parse_flag(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, bool def = false);

///
/// Parses an option as described in the intro.
///
template <typename T>
YCMD_API T parse_opt(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, const T& def,
    bool required = false, const std::vector<T>& choices = {});

///
/// Specialization of parse_opt()
///
YCMD_API std::string parse_opts(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help,
    const std::string& def, bool required = false,
    const std::vector<std::string>& choices = {});

///
/// Specialization of parse_opt()
///
YCMD_API int parse_opti(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, int def,
    bool required = false, const std::vector<int>& choices = {});

///
/// Specialization of parse_opt()
///
YCMD_API float parse_optf(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, float def,
    bool required = false, const std::vector<float>& choices = {});

///
/// Specialization of parse_opt()
///
YCMD_API double parse_optd(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, double def,
    bool required = false, const std::vector<double>& choices = {});

///
/// Specialization of parse_opt()
///
YCMD_API bool parse_optb(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, bool def,
    bool required = false);

///
/// Parses an option enum as described in the intro.
///
YCMD_API int parse_opte(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, int def,
    const std::vector<std::pair<std::string, int>>& vals,
    bool required = false);

///
/// Parses an option array as described in the intro.
///
template <typename T>
YCMD_API T parse_opta(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, const T& def,
    int nargs, bool required = false, const std::vector<T>& choices = {});

///
/// Parses an argument as described in the intro.
///
template <typename T>
YCMD_API T parse_arg(parser* par, const std::string& longname,
    const std::string& help, const T& def, bool required = false,
    const std::vector<T>& choices = {});

///
/// Specialization of parse_arg()
///
YCMD_API std::string parse_args(parser* par, const std::string& longname,
    const std::string& help, const std::string& def, bool required = false,
    const std::vector<std::string>& choices = {});

///
/// Specialization of parse_arg()
///
YCMD_API int parse_argi(parser* par, const std::string& longname,
    const std::string& help, int def, bool required = false,
    const std::vector<int>& choices = {});

///
/// Specialization of parse_arg()
///
YCMD_API float parse_argf(parser* par, const std::string& longname,
    const std::string& help, float def, bool required = false,
    const std::vector<float>& choices = {});

///
/// Specialization of parse_arg()
///
YCMD_API double parse_argd(parser* par, const std::string& longname,
    const std::string& help, double def, bool required = false,
    const std::vector<double>& choices = {});

///
/// Parses an argument array as described in the intro.
///
template <typename T>
YCMD_API std::vector<T> parse_arga(parser* par, const std::string& longname,
    const std::string& help, const std::vector<T>& def, int nargs = -1,
    bool required = false, const std::vector<T>& choices = {});

///
/// Specialization of parse_arga()
///
YCMD_API std::vector<std::string> parse_argas(parser* par,
    const std::string& longname, const std::string& help,
    const std::vector<std::string>& def, int nargs = -1, bool required = false,
    const std::vector<std::string>& choices = {});

///
/// Specialization of parse_arga()
///
YCMD_API std::vector<int> parse_argai(parser* par, const std::string& longname,
    const std::string& help, const std::vector<int>& def, int nargs = -1,
    bool required = false, const std::vector<int>& choices = {});

///
/// Specialization of parse_arga()
///
YCMD_API std::vector<float> parse_argaf(parser* par,
    const std::string& longname, const std::string& help,
    const std::vector<float>& def, int nargs = -1, bool required = false,
    const std::vector<float>& choices = {});

///
/// Specialization of parse_arga()
///
YCMD_API std::vector<double> parse_argad(parser* par,
    const std::string& longname, const std::string& help,
    const std::vector<double>& def, int nargs = -1, bool required = false,
    const std::vector<double>& choices = {});

// FILE LOADING
// ----------------------------------------------------------------

///
/// Loads the contents of a binary file in an in-memory array.
///
YCMD_API std::vector<unsigned char> load_binfile(const std::string& filename);

///
/// Loads the contents of a text file into a std::string.
///
YCMD_API std::string load_txtfile(const std::string& filename);

// PATH MANIPULATION
// -----------------------------------------------------------

///
/// Get directory name (including '/').
///
YCMD_API std::string get_dirname(const std::string& filename);

///
/// Get file basename.
///
YCMD_API std::string get_basename(const std::string& filename);

///
/// Get extension (including '.').
///
YCMD_API std::string get_extension(const std::string& filename);

///
/// Get filename without directory (equiv to get_basename() +
/// get_extension()).
///
YCMD_API std::string get_filename(const std::string& filename);

///
/// Replace extension.
///
YCMD_API std::string replace_extension(
    const std::string& filename, const std::string& ext);

///
/// Prepend a string to the extension.
///
YCMD_API std::string prepend_extension(
    const std::string& filename, const std::string& prep);

///
/// Splits a path calling the above functions.
///
YCMD_API void split_path(const std::string& filename, std::string& dirname,
    std::string& basename, std::string& ext);

// STRING UTILITIES
// ------------------------------------------------------------

///
/// Checks if a std::string starts with a prefix.
///
YCMD_API bool starts_with(const std::string& str, const std::string& substr);

///
/// Checks if a std::string ends with a prefix.
///
YCMD_API bool ends_with(const std::string& str, const std::string& substr);

///
/// Splits a std::string into lines at the '\n' character. The line
/// terminator is kept if keep_newline. This function does not work on
/// Window if keep_newline is true.
///
YCMD_API std::vector<std::string> split_lines(
    const std::string& str, bool keep_newline = false);

///
/// Joins a list of std::string with a std::string as separator.
///
YCMD_API std::string join_strings(
    const std::vector<std::string>& strs, const std::string& sep);

///
/// C-like string formatting. This is only meant for short strings with max
/// length 10000 chars. Memory corruption will happen for longer strings.
///
YCMD_API std::string format_str(const char* fmt, va_list args);

///
/// C-like string formatting. See format_str(fmt,args);
///
YCMD_API std::string format_str(const char* fmt, ...);

// SIMPLE LOGGING
// --------------------------------------------------------------

///
/// Logging level
///
enum log_level : int {
    /// verbose
    log_level_verbose = -3,
    /// trace
    log_level_trace = -2,
    /// debug
    log_level_debug = -1,
    /// info
    log_level_info = 0,
    /// warning
    log_level_warning = 1,
    /// error
    log_level_error = 2,
    /// fatal
    log_level_fatal = 3,
};

///
/// Logger object
///
struct logger;

///
/// Make a file logger. See set_logger() for parameters.
///
/// Parameters:
/// - filename: logger filename or stderr if empty
/// - apend: append or write open mode for file logger
///
///
YCMD_API logger* make_file_logger(const std::string& filename, bool append,
    int output_level = log_level_info, int flush_level = log_level_info);

///
/// Make a stderr logger. See set_logger() for parameters.
///
YCMD_API logger* make_stderr_logger(
    int output_level = log_level_info, int flush_level = log_level_info);

///
/// Make a stderr logger. See set_logger() for parameters.
///
YCMD_API logger* make_stdout_logger(
    int output_level = log_level_info, int flush_level = log_level_info);
///
/// Clear a logger
///
YCMD_API void clear_logger(logger* lgr);

///
/// Get default loggers. This is a modifiable reference.
///
YCMD_API std::vector<logger*>* get_default_loggers();

///
/// Set logger level
///
/// Parameters:
/// - lgr: logger
/// - name : logger name
/// - output_level: output level
/// - flush_level: output level
///
YCMD_API void set_logger(
    logger* lgr, int output_level, int flush_level = log_level_error);

///
/// Log a message
///
/// Parameters:
/// - lgr: logger
/// - level: message level
/// - code: message code (5 chars)
/// - msg: message
///
YCMD_API void log_msg(
    logger* lgr, int level, const std::string& name, const std::string& msg);

///
/// Log a message
///
/// Parameters:
/// - lgr: logger
/// - level: message level
/// - code: message code (5 chars)
/// - msg: message
///
YCMD_API void log_msg(
    logger* lgr, int level, const char* name, const char* msg);

///
/// Log a message formatted ala printf.
///
/// Parameters:
/// - lgr: logger
/// - level: message level
/// - code: message code (5 chars)
/// - msg: message
///
YCMD_API void log_msgf(
    logger* lgr, int level, const char* name, const char* msg, ...);

///
/// Logs a message to the default loggers
///
YCMD_API void log_msgfv(
    logger* lgr, int level, const char* name, const char* msg, va_list args);

///
/// Logs a message to the default loggers
///
YCMD_API void log_msg(
    int level, const std::string& name, const std::string& msg);

///
/// Logs a message to the default loggers
///
YCMD_API void log_msgf(int level, const char* name, const char* msg, ...);

///
/// Logs a message to the default loggers
///
YCMD_API void log_msgfv(
    int level, const char* name, const char* msg, va_list args);

// THREAD POOL
// -----------------------------------------------------------------

///
/// Forward declaration of thread pool.
///
struct thread_pool;

///
/// Initialize a thread pool with a certain number of threads (0 for
/// defatul).
///
YCMD_API thread_pool* make_thread_pool(int nthread = 0);

///
/// Clear thread pool
///
YCMD_API void clear_thread_pool(thread_pool* pool);

///
/// Wait for all jobs to finish
///
YCMD_API void thread_pool_wait(thread_pool* pool);

///
/// Parallel for implementation
///
YCMD_API void thread_pool_for(
    thread_pool* pool, int count, const std::function<void(int idx)>& task);

///
/// Runs a task asynchronously onto a thread pool
///
YCMD_API std::shared_future<void> thread_pool_async(
    thread_pool* pool, const std::function<void()>& task);

///
/// Wait for all jobs to finish on a global thread pool
///
YCMD_API void thread_pool_wait();

///
/// Runs a task asynchronously onto a global thread pool
///
YCMD_API std::shared_future<void> thread_pool_async(
    const std::function<void()>& task);

///
/// Parallel for implementation on a global thread pool
///
YCMD_API void thread_pool_for(
    int count, const std::function<void(int idx)>& task);

}  // namespace

// -----------------------------------------------------------------------------
// INCLUDE FOR HEADER-ONLY MODE
// -----------------------------------------------------------------------------

#ifdef YCMD_INLINE
#include "yocto_cmd.cpp"
#endif

#endif
