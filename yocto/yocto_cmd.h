///
/// YOCTO_CMD: utilities for writing command line applications, mostly
/// a simple to use command line parser and a few string, path and file
/// functions.
///
///
/// USAGE FOR COMMAND LINE PARSING:
///
/// 1. include this file (more compilation options below)
/// 2. create a parser object
///     - an option for printing help is automatically added
///     parser = make_parser(argc, argv, program description)
/// 3. for each option, parse it calling the functions parse_opt<type>
///     - options are parsed on the fly and a comprehensive help is
///       automatically generated
///     - supports bool (flags), int, float, double, std::string, enums
///     - options names are "--longname" for longname and "-s" for short
///     - command line format is "--longname value", "-s v" for all but flags
///     - values are parsed with iostream << operators
///     - for general use opt = parse_opt<type>(...)
///     - for boolean flags is parse_flag(...)
///     - for enums use parse_opte(...)
/// 4. for each unnamed argument, parse it calling the functions parse_arg<type>
///     - names are only used for help
///     - supports types as above
///     - for general use arg = parse_arg<type>(...)
///     - to parse all remaining values use args = parse_args<type>(...)
/// 5. end cmdline parsing with done_parser
///    - check for unsued values, missing arguments and print help if needed
///    check_parser(parser)
/// 6. since arguments are parsed immediately, one can easily implement
///    subcommands by just branching the command line code based on a read
///    argument without any need for complex syntax
///
/// The interface for each function is described in details in the interface
/// section of this file.
///
/// Notes: the end of this file contains a test function that also
/// illustreates the library usage.
///
///
/// UTILITIES
///
/// 1. filename splitting with get_dirname, get_basename, get_extension
/// 2. loading entire files with load_txtfile and load_binfile
///
/// The interface for each function is described in details in the interface
/// section of this file.
///
/// COMPILATION:
///
/// All functions in this library are inlined, so just inlucde the header.
///
///
/// HISTORY:
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
// Copyright (c) 2016 Fabio Pellacini
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

#ifndef _YCMD_H_
#define _YCMD_H_

#include <string>
#include <unordered_map>
#include <vector>

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

namespace ycmd {

// COMMAND LINE PARSING --------------------------------------------------------

///
/// Command line parser.
///
struct parser;

///
/// Inits a command line parser.
///
inline parser* make_parser(const std::vector<std::string>& args,
                           const std::string& help = "");

///
/// Inits a command line parser.
///
inline parser* make_parser(int argc, char* argv[], const char* help);

///
/// Ends parsing checking for error for unused options or arguments.
/// Exit if needed.
///
inline bool check_parser(parser* prs);

///
/// Parses an optional flag as described in the intro.
///
inline bool parse_flag(parser* par, const std::string& longname,
                       const std::string& shortname, const std::string& help,
                       bool def = false);

///
/// Parses an option as described in the intro.
///
template <typename T>
inline T parse_opt(parser* par, const std::string& longname,
                   const std::string& shortname, const std::string& help,
                   const T& def, bool required = false,
                   const std::vector<T>& choices = {});

///
/// Parses an option enum as described in the intro.
///
template <typename T>
inline T parse_opte(parser* par, const std::string& longname,
                    const std::string& shortname, const std::string& help,
                    const T& def,
                    const std::unordered_map<std::string, T>& vals,
                    bool required = false);

///
/// Parses an option array as described in the intro.
///
template <typename T>
inline T parse_opta(parser* par, const std::string& longname,
                    const std::string& shortname, const std::string& help,
                    const T& def, int nargs, bool required = false,
                    const std::vector<T>& choices = {});

///
/// Parses an argument as described in the intro.
///
template <typename T>
inline T parse_arg(parser* par, const std::string& longname,
                   const std::string& help, const T& def, bool required = false,
                   const std::vector<T>& choices = {});

///
/// Parses an argument array as described in the intro.
///
template <typename T>
inline std::vector<T> parse_arga(parser* par, const std::string& longname,
                                 const std::string& help,
                                 const std::vector<T>& def, int nargs = -1,
                                 bool required = false,
                                 const std::vector<T>& choices = {});

// FILE LOADING ----------------------------------------------------------------

///
/// Loads the contents of a binary file in an in-memory array.
///
inline std::vector<unsigned char> load_binfile(const std::string& filename);

///
/// Loads the contents of a text file into a std::string.
///
inline std::string load_txtfile(const std::string& filename);

// PATH MANIPULATION -----------------------------------------------------------

///
/// Get directory name (including '/').
///
inline std::string get_dirname(const std::string& filename);

///
/// Get file basename.
///
inline std::string get_basename(const std::string& filename);

///
/// Get extension (including '.').
///
inline std::string get_extension(const std::string& filename);

///
/// Get filename without directory (equiv to get_basename() + get_dirname().
///
inline std::string get_filename(const std::string& filename);

///
/// Splits a path calling the above functions.
///
inline void split_path(const std::string& filename, std::string& dirname,
                       std::string& basename, std::string& ext);

// STRING UTILITIES ------------------------------------------------------------

///
/// Checks if a std::string starts with a prefix.
///
inline bool starts_with(const std::string& str, const std::string& substr);

///
/// Checks if a std::string ends with a prefix.
///
inline bool ends_with(const std::string& str, const std::string& substr);

///
/// Splits a std::string into lines at the '\n' character. The line terminator
/// is kept if keep_newline. This function does not work on Window if
/// keep_newline is true.
///
inline std::vector<std::string> split_lines(const std::string& str,
                                            bool keep_newline = false);

///
/// Joins a list of std::string with a std::string as separator.
///
inline std::string join_strings(const std::vector<std::string>& strs,
                                const std::string& sep);

}  // namespace

// -----------------------------------------------------------------------------
// IMPLEMENTATION
// -----------------------------------------------------------------------------

#include <cassert>
#include <cstdio>
#include <cstring>
#include <functional>
#include <sstream>

namespace ycmd {

//
// Command-line parsing state
//
struct parser {
    // saved arguments ----------------------------------------
    std::vector<std::string> args;

    // help std::strings -------------------------------------------
    std::string help_prog;   // program description
    std::string help_usage;  // sage lines
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
inline parser* make_parser(const std::vector<std::string>& args,
                           const std::string& help) {
    // clears parser and copy argument data
    auto par = new parser();
    par->args = std::vector<std::string>(args.begin() + 1, args.end());
    par->error = false;
    par->exit_on_error = true;
    par->help_prog = args[0];
    par->help_usage += help;
    par->print_help = parse_flag(par, "--help", "-h", "print help", false);
    return par;
}

//
// Inits the parser.
//
inline parser* make_parser(int argc, char* argv[], const char* help) {
    return make_parser(std::vector<std::string>(argv, argv + argc), help);
}

//
// Print help based on the help lines collected during parsing.
//
static inline void _print_help(parser* par) {
    auto help = std::string();
    help += "usage: " + par->help_prog;
    if (not par->help_opts.empty()) help += "[options] ";
    if (not par->help_args.empty()) help += "<arguments> ";
    help += "\n    " + par->help_usage + "\n\n";
    if (!par->help_opts.empty()) {
        help += "options:\n";
        help += par->help_opts;
    }
    if (!par->help_args.empty()) {
        help += "arguments:\n";
        help += par->help_args;
    }
    printf("%s\n", help.c_str());
}

//
// Ends parsing checking for error for unused options or arguments.
// Exit if needed.
//
inline bool check_parser(parser* par) {
    // check for error
    if (!par->error && par->args.size() > 0) {
        par->error = true;
        if (par->args[0][0] == '-')
            par->error_msg = "unknown option " + par->args[0];
        else
            par->error_msg = "unsued values";
    };
    // check whether we need to print help and exit
    if (par->error) printf("error: %s\n", par->error_msg.c_str());
    if (par->print_help || par->error) {
        _print_help(par);
        if (par->exit_on_error) exit(EXIT_FAILURE);
    }
    auto ok = !par->error;
    delete par;
    return ok;
}

//
// Check if a std::string starts with another.
//
static inline bool _startswith(const std::string& str,
                               const std::string& start) {
    if (str.length() < start.length()) return false;
    for (auto i = 0; i < start.length(); i++) {
        if (str[i] != start[i]) return false;
    }
    return true;
}

//
// Check if an option name is valid
//
static inline void _check_name(parser* par, const std::string& longname,
                               const std::string& shortname, bool opt,
                               int nargs) {
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
template <typename T>
static inline std::string _typename() {
    return "";
}
template <>
inline std::string _typename<bool>() {
    return "bool";
}
template <>
inline std::string _typename<int>() {
    return "int";
}
template <>
inline std::string _typename<float>() {
    return "float";
}
template <>
inline std::string _typename<double>() {
    return "double";
}
template <>
inline std::string _typename<std::string>() {
    return "std::string";
}

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
static inline void _add_help(parser* par, const std::string& longname,
                             const std::string& shortname,
                             const std::string& help, bool opt, bool req,
                             int nargs, const std::vector<T>& def,
                             const std::vector<T>& choices = {}) {
    // full name
    auto help_fullname = longname;
    if (!shortname.empty()) help_fullname += "/" + shortname;
    if (nargs != 0) help_fullname += " " + _typename<T>();

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
        par->help_opts += help_line;
    else
        par->help_args += help_line;
}

//
// Parsing routine for arrays of values
//
template <typename T>
static inline std::vector<T> _parse_vals(
    parser* par, const std::string& longname, const std::string& shortname,
    const std::string& help, bool opt, bool req, int nargs,
    const std::vector<T>& def, const std::vector<T>& choices = {}) {
    // prepare default empty vec
    auto vals = std::vector<T>();

    // check whether the name is good
    _check_name(par, longname, shortname, opt, nargs);

    // add help
    _add_help(par, longname, shortname, help, opt, req, nargs, def, choices);

    // skip if alreasy in error
    if (par->error) return vals;

    // find the value position
    auto val_pos = -1;
    if (opt) {
        // find option name
        for (auto i = 0; i < par->args.size() && val_pos < 0; i++) {
            if (shortname == par->args[i]) val_pos = i;
            if (longname == par->args[i]) val_pos = i;
        }

        // remove the option name
        if (val_pos >= 0) {
            par->args.erase(par->args.begin() + val_pos);
        }
    } else {
        // check if arg is present
        if (!par->args.empty()) {
            if (par->args[0][0] != '-') {
                val_pos = 0;
            }
        }
    }

    // handle not found
    if (val_pos < 0) {
        if (req) {
            par->error = true;
            par->error_msg = "missing value for " + longname;
            return vals;
        } else
            return def;
    }

    // check if value is present
    if (nargs == -1) nargs = std::max(1, (int)par->args.size());
    if (val_pos + nargs > par->args.size()) {
        par->error = true;
        par->error_msg = "missing value for " + longname;
        return vals;
    }

    // loop over values
    for (auto i = 0; i < nargs; i++) {
        // grab value
        auto val_str = par->args[val_pos];
        par->args.erase(par->args.begin() + val_pos);

        // parse value
        auto stream = std::istringstream(val_str);
        auto val = T();
        if (!(stream >> val)) {
            par->error = true;
            par->error_msg = "incorrect value for " + longname;
            return vals;
        }
        // check choices
        if (!choices.empty()) {
            auto in_choices = false;
            for (auto&& c : choices) {
                if (val == c) in_choices = true;
            }
            if (!in_choices) {
                par->error = true;
                par->error_msg = "incorrect value for " + longname;
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
static inline T _parse_val(parser* par, const std::string& longname,
                           const std::string& shortname,
                           const std::string& help, bool opt, const T& def,
                           bool req, const std::vector<T>& choices = {}) {
    // parse values
    auto vals = _parse_vals(par, longname, shortname, help, opt, req, 1, {def},
                            choices);

    // return value if present
    if (vals.size())
        return vals[0];
    else
        return {};
}

//
// Parses an optional flag as described in the intro.
//
inline bool parse_flag(parser* par, const std::string& longname,
                       const std::string& shortname, const std::string& help,
                       bool def) {
    // parse values
    auto vals = _parse_vals<bool>(par, longname, shortname, help, true, false,
                                  0, {def});

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
inline T parse_opt(parser* par, const std::string& longname,
                   const std::string& shortname, const std::string& help,
                   const T& def, bool req, const std::vector<T>& choices) {
    return _parse_val<T>(par, longname, shortname, help, true, def, req,
                         choices);
}

//
// Parses an argument as described in the intro.
//
template <typename T>
inline T parse_arg(parser* par, const std::string& longname,
                   const std::string& help, const T& def, bool req,
                   const std::vector<T>& choices) {
    return _parse_val<T>(par, longname, "", help, false, def, req, choices);
}

//
// Parses an option array as described in the intro.
//
template <typename T>
inline std::vector<T> parse_opta(parser* par, const std::string& longname,
                                 const std::string& help,
                                 const std::vector<T>& def, int nargs,
                                 bool required, const std::vector<T>& choices) {
    return _parse_vals(par, longname, "", help, true, required, nargs, def,
                       choices);
}

//
// Parses an argument array as described in the intro.
//
template <typename T>
inline std::vector<T> parse_arga(parser* par, const std::string& longname,
                                 const std::string& help,
                                 const std::vector<T>& def, int nargs,
                                 bool required, const std::vector<T>& choices) {
    return _parse_vals(par, longname, "", help, false, required, nargs, def,
                       choices);
}

//
// Parses an option enum as described in the intro.
//
template <typename T>
inline T parse_opte(parser* par, const std::string& longname,
                    const std::string& shortname, const std::string& help,
                    const T& def,
                    const std::unordered_map<std::string, T>& vals,
                    bool required) {
    auto choices = std::vector<std::string>();
    auto def_s = std::string();
    for (auto kv : vals) {
        choices.push_back(kv.first);
        if (kv.second == def) def_s = kv.first;
    }
    auto val_s =
        parse_opt(par, longname, shortname, help, def_s, required, choices);
    return vals.at(val_s);
}
//
// Parse a binary file a chuck at a time and loads its content full in memory.
// Does not attempt to determine file size upfront to handle all cases.
//
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

//
// Parse a text file a chuck at a time and loads its content full in memory.
// Does not attempt to determine file size upfront to handle all cases.
//
inline std::string load_txtfile(const std::string& filename) {
    auto file = fopen(filename.c_str(), "rt");
    if (!file) return 0;
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

//
// Get directory name (including '/').
//
inline std::string get_dirname(const std::string& filename) {
    auto pos = filename.rfind('/');
    if (pos == std::string::npos) pos = filename.rfind('\\');
    if (pos == std::string::npos) return "";
    return filename.substr(0, pos + 1);
}

//
// Get file basename.
//
inline std::string get_basename(const std::string& filename) {
    auto dirname = get_dirname(filename);
    auto extension = get_extension(filename);
    return filename.substr(dirname.size(),
                           filename.size() - dirname.size() - extension.size());
}

//
// Get extension (including '.').
//
inline std::string get_extension(const std::string& filename) {
    auto pos = filename.rfind('.');
    if (pos == std::string::npos) return "";
    return filename.substr(pos);
}

//
// Get filename without directory (equiv to get_basename() + get_dirname().
//
inline std::string get_filename(const std::string& filename) {
    return get_basename(filename) + get_extension(filename);
}

//
// Splits a path. Public interface.
//
inline void split_path(const std::string& filename, std::string& dirname,
                       std::string& basename, std::string& ext) {
    dirname = get_dirname(filename);
    basename = get_basename(filename);
    ext = get_extension(filename);
}

//
// Checks if a std::string starts with a prefix.
//
inline bool starts_with(const std::string& str, const std::string& substr) {
    if (str.length() < substr.length()) return false;
    for (auto i = 0; i < substr.length(); i++)
        if (str[i] != substr[i]) return false;
    return true;
}

//
// Checks if a std::string ends with a prefix.
//
inline bool ends_with(const std::string& str, const std::string& substr) {
    if (str.length() < substr.length()) return false;
    auto offset = str.length() - substr.length();
    for (auto i = 0; i < substr.length(); i++)
        if (str[i + offset] != substr[i]) return false;
    return true;
}

//
// Splits a std::string into lines at the '\n' character.
//
inline std::vector<std::string> split_lines(const std::string& str,
                                            bool keep_newline) {
    if (str.empty()) return {};
    auto lines = std::vector<std::string>();
    auto last = 0;
    auto pos = str.find('\n');
    while (pos != str.npos) {
        lines.push_back(str.substr(last, pos + ((keep_newline) ? 1 : 0)));
        last = (int)pos + 1;
        pos = str.find('\n', last + 1);
    }
    return lines;
}

//
// Joins a list of std::string with a std::string as separator.
//
inline std::string join_strings(const std::vector<std::string>& strs,
                                const std::string& sep) {
    auto ret = std::string();
    auto first = true;
    for (auto& str : strs) {
        if (!first) ret += sep;
        ret += str;
        first = false;
    }
    return ret;
}

}  // namespace

// -----------------------------------------------------------------------------

#ifdef YGL_TEST_CMD
#include <math.h>

int main(int argc, char** argv) {
    // test empty
    auto test0_argc = 1;
    const char* test0_argv[] = {"test0"};
    auto par0 = ycmd::make_parser(test0_argc, (char**)test0_argv, "test0");
    par0->exit_on_error = false;
    assert(check_parser(par0) == true);

    // test exit on help
    auto test1_argc = 2;
    const char* test1_argv[] = {"test1", "-h"};
    auto par1 = ycmd::make_parser(test1_argc, (char**)test1_argv, "test1");
    par1->exit_on_error = false;
    assert(check_parser(par1) == true);

    // test opts
    auto test2_argc = 10;
    const char* test2_argv[] = {"test2", "--int",    "10",   "--float",
                                "3.14",  "--double", "6.28", "--str",
                                "bob",   "--flag"};
    auto par2 = ycmd::make_parser(test2_argc, (char**)test2_argv, "test2");
    par2->exit_on_error = false;
    assert(parse_flag(par2, "--flag", "", "", false) == true);
    assert(ycmd::parse_opt<int>(par2, "--int", "", "", 0) == 10);
    assert(fabsf(ycmd::parse_opt<float>(par2, "--float", "", "", 0) - 3.14f) <
           0.01);
    assert(fabs(ycmd::parse_opt<double>(par2, "--double", "", "", 0) - 6.28) <
           0.01);
    assert(ycmd::parse_opt<std::string>(par2, "--str", "", "", "mike") ==
           "bob");
    assert(parse_flag(par2, "--flag_def", "", "", false) == false);
    assert(ycmd::parse_opt<int>(par2, "--int_def", "", "", 5) == 5);
    assert(ycmd::parse_opt<float>(par2, "--float_def", "", "", 2.67f) == 2.67f);
    assert(ycmd::parse_opt<double>(par2, "--double_def", "", "", 9.54) == 9.54);
    assert(ycmd::parse_opt<std::string>(par2, "--str_def", "", "", "alex") ==
           "alex");
    assert(check_parser(par2) == true);

    // test args
    auto test3_argc = 3;
    const char* test3_argv[] = {"test3", "10", "bob"};
    auto par3 = ycmd::make_parser(test3_argc, (char**)test3_argv, "test3");
    par3->exit_on_error = false;
    assert(ycmd::parse_arg<int>(par3, "int", "", 0, true) == 10);
    assert(ycmd::parse_arg<std::string>(par3, "str", "", "mike", true) ==
           "bob");
    assert(check_parser(par3) == true);

    // test bad opts
    auto test4_argc = 3;
    const char* test4_argv[] = {"test4", "--int", "bob"};
    auto par4 = ycmd::make_parser(test4_argc, (char**)test4_argv, "test4");
    par4->exit_on_error = false;
    assert(ycmd::parse_opt<int>(par4, "--int", "", "", 0) == 0);
    assert(check_parser(par4) == false);
}

#endif

#endif
