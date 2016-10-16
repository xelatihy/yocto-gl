//
// DEPRECATED: Use YOCTO_UTILS as a substitute.
//
// YOCTO_CMDLINE: utilities for writing command line applications, mostly
// a simple to use command line parser.
//

//
// USAGE:
//
// 0. include this file (more compilation options below)
// 1. create a parser object
//     - an option for printing help is automatically added
//     parser = yc_init_parser(argc, argv, program description)
// 2. for each option, parse it calling the functions yc_parse_optXXX
//     - options are parsed on the flh and a comprehensive help is
//       automatically generated
//     - supports bool (flags), int, float, double, string, enums
//     - options names are "--longname" for longname and "-s" for short
//     - command line format is "--longname value", "-s v" for all but flags
//     opt = yc_parse_opti(parser, longname, flag, description, default value)
// 3. for each unnamed argument, parse it calling the functions yc_parse_argXXX
//     - names are only used for help
//     - supports above type and all a special one that comsumes all args
//     saving them in a string array
//     arg = yc_parse_args(parser, longname, description, deafult value,
//               required)
// 4. end cmdline parsing with yc_done_parser
//    - check for unsued values, missing arguments and print help if needed
//    yc_done_parser(parser)
// 5. since arguments are parsed immediately, one can easily implement
//    subcommands by just branching the command line code based on a read
//    argument without any need for complex syntax
// 6. contains also utilities to load whole files in memory and split paths
//
// The interface for each function is described in details in the interface
// section of this file.
//
// Notes: the end of this file contains a test function that also
// illustreates the library usage.
//

//
// COMPILATION:
//
// All functions in this library are inlined by default for ease of use in
// C++. To use the library as a .h/.cpp pair do the following:
// - to use as a .h, just #define YGL_DECLARATION before including this file
// - to build as a .cpp, just #define YGL_IMPLEMENTATION before including this
// file into only one file that you can either link directly or pack as a lib.
//
// This file depends on yocto_math.h.
//

//
// HISTORY:
// - v 0.2: [API change] C++ API
// - v 0.1: C++ implementation
// - v 0.0: initial release in C99
//

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

#ifndef _YC_H_
#define _YC_H_

#include <string>
#include <unordered_map>
#include <vector>

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

// COMMAND LINE PARSING --------------------------------------------------------

//
// Forward declaration of an opaque parser object.
//
struct yc_parser;

//
// Inits a command line parser.
//
inline yc_parser* yc_init_parser(const std::vector<std::string>& args,
                                 const std::string& help = "");

//
// Inits a command line parser.
//
inline yc_parser* yc_init_parser(int argc, char* argv[], const char* help);

//
// Ends parsing checking for error for unused options or arguments.
// Exit if needed.
//
inline bool yc_done_parser(yc_parser* parser);

//
// Parses an optional flag as described in the intro.
//
inline bool yc_parse_flag(yc_parser* parser, const std::string& longname,
                          const std::string& shortname, const std::string& help,
                          bool def = false);

//
// Parses an option as described in the intro.
//
template <typename T>
inline T yc_parse_opt(yc_parser* parser, const std::string& longname,
                      const std::string& shortname, const std::string& help,
                      const T& def, bool required = false,
                      const std::vector<T>& choices = {});

//
// Parses an option enum as described in the intro.
//
template <typename T>
inline T yc_parse_opte(yc_parser* parser, const std::string& longname,
                       const std::string& shortname, const std::string& help,
                       const T& def,
                       const std::unordered_map<std::string, T>& vals,
                       bool required = false);

//
// Parses an option array as described in the intro.
//
template <typename T>
inline T yc_parse_opta(yc_parser* parser, const std::string& longname,
                       const std::string& shortname, const std::string& help,
                       const T& def, int nargs, bool required = false,
                       const std::vector<T>& choices = {});

//
// Parses an argument as described in the intro.
//
template <typename T>
inline T yc_parse_arg(yc_parser* parser, const std::string& longname,
                      const std::string& help, const T& def,
                      bool required = false,
                      const std::vector<T>& choices = {});

//
// Parses an argument array as described in the intro.
//
template <typename T>
inline std::vector<T> yc_parse_arga(yc_parser* parser,
                                    const std::string& longname,
                                    const std::string& help,
                                    const std::vector<T>& def, int nargs = -1,
                                    bool required = false,
                                    const std::vector<T>& choices = {});

// FILE LOADING ----------------------------------------------------------------

//
// Loads the contents of a binary file in an in-memory array.
//
inline std::vector<unsigned char> yc_load_binfile(const std::string& filename);

//
// Loads the contents of a text file into a string.
//
inline std::string yc_load_txtfile(const std::string& filename);

// PATH MANIPULATION -----------------------------------------------------------

//
// Splits paths into directory name (including '/'), basename, and extension
// (including '.').
// If any the pointers is missing, that component will be skipped.
//
inline void yc_split_path(const std::string& filename, std::string* dirname,
                          std::string* basename, std::string* ext);

// -----------------------------------------------------------------------------
// DEPRECATED INTERFACE (WRAPPERS TO ABOVE FUNCTIONS)
// -----------------------------------------------------------------------------

//
// Wrapper. REquries preallocated strings.
//
inline void yc_split_path(const char* filename, char* dirname, char* basename,
                          char* ext);

//
// Wrapper.
//
inline bool yc_parse_optb(yc_parser* parser, const char* longname,
                          const char* shortname, const char* help, bool def);

//
// Wrapper.
//
inline int yc_parse_opti(yc_parser* parser, const char* longname,
                         const char* shortname, const char* help, int def);

//
// Wrapper.
//
inline float yc_parse_optf(yc_parser* parser, const char* longname,
                           const char* shortname, const char* help, float def);

//
// Wrapper.
//
inline double yc_parse_optd(yc_parser* parser, const char* longname,
                            const char* shortname, const char* help,
                            double def);

//
// Wrapper. Need to cleanup returned memory.
//
inline const char* yc_parse_opts(yc_parser* parser, const char* longname,
                                 const char* shortname, const char* help,
                                 const char* def);

//
// Wrapper.
//
inline int yc_parse_opte(yc_parser* parser, const char* longname,
                         const char* shortname, const char* help, int def,
                         const char* vals[]);

//
// Wrapper. Need to cleanup returned memory.
//
inline char* yc_parse_args(yc_parser* parser, const char* longname,
                           const char* help, const char* def, bool required);

//
// Wrapper. Need to cleanup returned memory.
//
inline char** yc_parse_argsa(yc_parser* parser, const char* longname,
                             const char* help, const char* def, bool required,
                             char** strs, int max_strs);

// -----------------------------------------------------------------------------
// IMPLEMENTATION
// -----------------------------------------------------------------------------

#include <cassert>
#include <cstdio>
#include <cstring>
#include <functional>
#include <sstream>

//
// Command-line parsing state
//
struct yc_parser {
    // saved arguments ----------------------------------------
    std::vector<std::string> args;

    // help strings -------------------------------------------
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
inline yc_parser* yc_init_parser(const std::vector<std::string>& args,
                                 const std::string& help) {
    // Allocates parser and copy argument data
    auto parser = new yc_parser();
    parser->args = std::vector<std::string>(args.begin() + 1, args.end());
    parser->error = false;
    parser->exit_on_error = true;
    parser->help_prog = args[0];
    parser->help_usage += help;
    parser->print_help =
        yc_parse_flag(parser, "--help", "-h", "print help", false);
    return parser;
}

//
// Inits the parser.
//
inline yc_parser* yc_init_parser(int argc, char* argv[], const char* help) {
    return yc_init_parser(std::vector<std::string>(argv, argv + argc), help);
}

//
// Print help based on the help lines collected during parsing.
//
static inline void yc__print_help(yc_parser* parser) {
    auto help = std::string();
    help += "usage: " + parser->help_prog;
    if (not parser->help_opts.empty()) help += "[options] ";
    if (not parser->help_args.empty()) help += "<arguments> ";
    help += "\n    " + parser->help_usage + "\n\n";
    if (!parser->help_opts.empty()) {
        help += "options:\n";
        help += parser->help_opts;
    }
    if (!parser->help_args.empty()) {
        help += "arguments:\n";
        help += parser->help_args;
    }
    printf("%s\n", help.c_str());
}

//
// Ends parsing checking for error for unused options or arguments.
// Exit if needed.
//
inline bool yc_done_parser(yc_parser* parser) {
    // check for error
    if (!parser->error && parser->args.size() > 0) {
        parser->error = true;
        if (parser->args[0][0] == '-')
            parser->error_msg = "unknown option " + parser->args[0];
        else
            parser->error_msg = "unsued values";
    };
    // check whether we need to print help and exit
    if (parser->error) printf("error: %s\n", parser->error_msg.c_str());
    if (parser->print_help || parser->error) {
        yc__print_help(parser);
        if (parser->exit_on_error) exit(EXIT_FAILURE);
    }
    return !parser->error;
}

//
// Check if a string starts with another.
//
static inline bool yc__startswith(const std::string& str,
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
static inline void yc__check_name(yc_parser* parser,
                                  const std::string& longname,
                                  const std::string& shortname, bool opt,
                                  int nargs) {
    // check name
    assert(!longname.empty());
    if (opt) {
        if (!shortname.empty())
            assert(yc__startswith(shortname, "-") && shortname.length() > 1);
        assert(yc__startswith(longname, "--") && longname.length() > 2);
    } else {
        assert(shortname.empty());
        assert(longname[0] != '-');
    }
    assert((opt && nargs >= 0) || (!opt && (nargs == -1 || nargs > 0)));
}

//
// Convert a type to string
//
template <typename T>
static inline std::string yc__typename() {
    return "";
}
template <>
inline std::string yc__typename<bool>() {
    return "bool";
}
template <>
inline std::string yc__typename<int>() {
    return "int";
}
template <>
inline std::string yc__typename<float>() {
    return "float";
}
template <>
inline std::string yc__typename<double>() {
    return "double";
}
template <>
inline std::string yc__typename<std::string>() {
    return "string";
}

//
// Converts a value to a string
//
template <typename T>
static inline std::string yc__tostring(const T& val) {
    std::ostringstream stream;
    stream << val;
    return stream.str();
}
template <>
inline std::string yc__tostring<bool>(const bool& val) {
    return (val) ? "true" : "false";
}

//
// Add a formatted help line
//
template <typename T>
static inline void yc__add_help(yc_parser* parser, const std::string& longname,
                                const std::string& shortname,
                                const std::string& help, bool opt, bool req,
                                int nargs, const std::vector<T>& def,
                                const std::vector<T>& choices = {}) {
    // full name
    auto help_fullname = longname;
    if (!shortname.empty()) help_fullname += "/" + shortname;
    if (nargs != 0) help_fullname += " " + yc__typename<T>();

    // default
    auto help_def = std::string();
    if (!req) {
        if (nargs == 0 || nargs == 1) {
            if (!def.empty())
                help_def = "[" + yc__tostring((const T&)def[0]) + "]";
        } else if (nargs == -1 || nargs > 1) {
            help_def = "[";
            for (auto i = 0; i < def.size(); i++) {
                if (i) help_def += ",";
                help_def += yc__tostring((const T&)def[i]);
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
            help_choice += yc__tostring((const T&)choices[i]);
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
        parser->help_opts += help_line;
    else
        parser->help_args += help_line;
}

//
// Parsing routine for arrays of values
//
template <typename T>
static inline std::vector<T> yc__parse_vals(
    yc_parser* parser, const std::string& longname,
    const std::string& shortname, const std::string& help, bool opt, bool req,
    int nargs, const std::vector<T>& def, const std::vector<T>& choices = {}) {
    // prepare default empty vec
    auto vals = std::vector<T>();

    // check whether the name is good
    yc__check_name(parser, longname, shortname, opt, nargs);

    // add help
    yc__add_help(parser, longname, shortname, help, opt, req, nargs, def,
                 choices);

    // skip if alreasy in error
    if (parser->error) return vals;

    // find the value position
    auto val_pos = -1;
    if (opt) {
        // find option name
        for (auto i = 0; i < parser->args.size() && val_pos < 0; i++) {
            if (shortname == parser->args[i]) val_pos = i;
            if (longname == parser->args[i]) val_pos = i;
        }

        // remove the option name
        if (val_pos >= 0) {
            parser->args.erase(parser->args.begin() + val_pos);
        }
    } else {
        // check if arg is present
        if (!parser->args.empty()) {
            if (parser->args[0][0] != '-') {
                val_pos = 0;
            }
        }
    }

    // handle not found
    if (val_pos < 0) {
        if (req) {
            parser->error = true;
            parser->error_msg = "missing value for " + longname;
            return vals;
        } else
            return def;
    }

    // check if value is present
    if (nargs == -1) nargs = std::max(1, (int)parser->args.size());
    if (val_pos + nargs > parser->args.size()) {
        parser->error = true;
        parser->error_msg = "missing value for " + longname;
        return vals;
    }

    // loop over values
    for (auto i = 0; i < nargs; i++) {
        // grab value
        auto val_str = parser->args[val_pos];
        parser->args.erase(parser->args.begin() + val_pos);

        // parse value
        auto stream = std::istringstream(val_str);
        auto val = T();
        if (!(stream >> val)) {
            parser->error = true;
            parser->error_msg = "incorrect value for " + longname;
            return vals;
        }
        // check choices
        if (!choices.empty()) {
            auto in_choices = false;
            for (auto&& c : choices) {
                if (val == c) in_choices = true;
            }
            if (!in_choices) {
                parser->error = true;
                parser->error_msg = "incorrect value for " + longname;
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
static inline T yc__parse_val(yc_parser* parser, const std::string& longname,
                              const std::string& shortname,
                              const std::string& help, bool opt, const T& def,
                              bool req, const std::vector<T>& choices = {}) {
    // parse values
    auto vals = yc__parse_vals(parser, longname, shortname, help, opt, req, 1,
                               {def}, choices);

    // return value if present
    if (vals.size())
        return vals[0];
    else
        return {};
}

//
// Parses an optional flag as described in the intro.
//
inline bool yc_parse_flag(yc_parser* parser, const std::string& longname,
                          const std::string& shortname, const std::string& help,
                          bool def) {
    // parse values
    auto vals = yc__parse_vals<bool>(parser, longname, shortname, help, true,
                                     false, 0, {def});

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
inline T yc_parse_opt(yc_parser* parser, const std::string& longname,
                      const std::string& shortname, const std::string& help,
                      const T& def, bool req, const std::vector<T>& choices) {
    return yc__parse_val<T>(parser, longname, shortname, help, true, def, req,
                            choices);
}

//
// Parses an argument as described in the intro.
//
template <typename T>
inline T yc_parse_arg(yc_parser* parser, const std::string& longname,
                      const std::string& help, const T& def, bool req,
                      const std::vector<T>& choices) {
    return yc__parse_val<T>(parser, longname, "", help, false, def, req,
                            choices);
}

//
// Parses an option array as described in the intro.
//
template <typename T>
inline std::vector<T> yc_parse_opta(yc_parser* parser,
                                    const std::string& longname,
                                    const std::string& help,
                                    const std::vector<T>& def, int nargs,
                                    bool required,
                                    const std::vector<T>& choices) {
    return yc__parse_vals(parser, longname, "", help, true, required, nargs,
                          def, choices);
}

//
// Parses an argument array as described in the intro.
//
template <typename T>
inline std::vector<T> yc_parse_arga(yc_parser* parser,
                                    const std::string& longname,
                                    const std::string& help,
                                    const std::vector<T>& def, int nargs,
                                    bool required,
                                    const std::vector<T>& choices) {
    return yc__parse_vals(parser, longname, "", help, false, required, nargs,
                          def, choices);
}

//
// Parses an option enum as described in the intro.
//
template <typename T>
inline T yc_parse_opte(yc_parser* parser, const std::string& longname,
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
    auto val_s = yc_parse_opt(parser, longname, shortname, help, def_s,
                              required, choices);
    return vals.at(val_s);
}
//
// Parse a binary file a chuck at a time and loads its content full in memory.
// Does not attempt to determine file size upfront to handle all cases.
//
inline std::vector<unsigned char> yc_load_binfile(const std::string& filename) {
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
inline std::string yc_load_txtfile(const std::string& filename) {
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
// Splits a path. Public interface.
//
inline void yc_split_path(const std::string& filename, std::string* dirname,
                          std::string* basename, std::string* ext) {
    // walk till end keeping the position of '/', '\\' and '.'
    auto path_sep = -1, ext_sep = -1;
    for (auto i = 0; i < filename.length(); i++) {
        if (filename[i] == '/' || filename[i] == '\\') path_sep = i;
        if (filename[i] == '.') ext_sep = i;
    }

    // copy strings
    if (dirname) {
        if (path_sep >= 0) {
            *dirname = filename.substr(0, path_sep + 1);
        } else {
            *dirname = "";
        }
    }
    if (basename) {
        auto start = (path_sep >= 0) ? path_sep + 1 : 0;
        if (ext_sep >= 0) {
            *basename = filename.substr(start, ext_sep);
        } else {
            *basename = filename.substr(start);
        }
    }
    if (ext) {
        if (ext_sep >= 0) {
            *ext = filename.substr(ext_sep);
        } else {
            *ext = "";
        }
    }
}

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF DEPRECATED INTERFACE
// -----------------------------------------------------------------------------

//
// Wrapper. REquries preallocated strings.
//
inline void yc_split_path(const char* filename, char* dirname, char* basename,
                          char* ext) {
    std::string dirname_s, basename_s, ext_s;
    yc_split_path(std::string(filename), (dirname) ? &dirname_s : nullptr,
                  (basename) ? &basename_s : nullptr, (ext) ? &ext_s : nullptr);
    if (dirname) strcpy(dirname, dirname_s.c_str());
    if (basename) strcpy(basename, basename_s.c_str());
    if (ext) strcpy(ext, ext_s.c_str());
}

//
// Wrapper.
//
inline bool yc_parse_optb(yc_parser* parser, const char* longname,
                          const char* shortname, const char* help, bool def) {
    return yc_parse_flag(parser, (longname) ? longname : "",
                         (shortname) ? shortname : "", (help) ? help : "", def);
}

//
// Wrapper.
//
inline int yc_parse_opti(yc_parser* parser, const char* longname,
                         const char* shortname, const char* help, int def) {
    return yc_parse_opt<int>(parser, (longname) ? longname : "",
                             (shortname) ? shortname : "", (help) ? help : "",
                             def);
}

//
// Wrapper.
//
inline float yc_parse_optf(yc_parser* parser, const char* longname,
                           const char* shortname, const char* help, float def) {
    return yc_parse_opt<float>(parser, (longname) ? longname : "",
                               (shortname) ? shortname : "", (help) ? help : "",
                               def);
}

//
// Wrapper.
//
inline double yc_parse_optd(yc_parser* parser, const char* longname,
                            const char* shortname, const char* help,
                            double def) {
    return yc_parse_opt<double>(parser, (longname) ? longname : "",
                                (shortname) ? shortname : "",
                                (help) ? help : "", def);
}

//
// Wrapper. Need to cleanup returned memory.
//
inline const char* yc_parse_opts(yc_parser* parser, const char* longname,
                                 const char* shortname, const char* help,
                                 const char* def) {
    auto s = yc_parse_opt<std::string>(parser, (longname) ? longname : "",
                                       (shortname) ? shortname : "",
                                       (help) ? help : "", def);
    auto str = new char[s.length() + 1];
    strcpy(str, s.c_str());
    return str;
}

//
// Wrapper. Need to cleanup returned memory.
//
inline char* yc_parse_args(yc_parser* parser, const char* longname,
                           const char* help, const char* def, bool required) {
    auto s = yc_parse_arg<std::string>(parser, (longname) ? longname : "",
                                       (help) ? help : "", (def) ? def : "",
                                       required);
    auto str = new char[s.length() + 1];
    strcpy(str, s.c_str());
    return str;
}

//
// Wrapper. Need to cleanup returned memory.
//
inline char** yc_parse_argsa(yc_parser* parser, const char* longname,
                             const char* help, const char* def, bool required,
                             char** strs, int max_strs) {
    auto s = yc_parse_arga<std::string>(parser, (longname) ? longname : "",
                                        (help) ? help : "", {}, -1, required);
    assert((int)s.size() < max_strs);
    for (auto i = 0; i < s.size(); i++) {
        strs[i] = new char[s[i].length() + 1];
        strcpy(strs[i], s[i].c_str());
    }
    strs[s.size()] = 0;
    return strs;
}

//
// Wrapper.
//
inline int yc_parse_opte(yc_parser* parser, const char* longname,
                         const char* shortname, const char* help, int def,
                         const char* vals[]) {
    auto vals_m = std::unordered_map<std::string, int>();
    for (auto i = 0;; i++) {
        if (!vals[i]) break;
        vals_m[vals[i]] = i;
    }
    return yc_parse_opte(parser, (longname) ? longname : "",
                         (shortname) ? shortname : "", (help) ? help : "", def,
                         vals_m);
}

// -----------------------------------------------------------------------------

#ifdef YGL_TEST_CMDLINE
#include <math.h>

int main(int argc, char** argv) {
    // test empty
    auto test0_argc = 1;
    const char* test0_argv[] = {"test0"};
    auto parser0 = yc_init_parser(test0_argc, (char**)test0_argv, "test0");
    parser0->exit_on_error = false;
    assert(yc_done_parser(parser0) == true);

    // test exit on help
    auto test1_argc = 2;
    const char* test1_argv[] = {"test1", "-h"};
    auto parser1 = yc_init_parser(test1_argc, (char**)test1_argv, "test1");
    parser1->exit_on_error = false;
    assert(yc_done_parser(parser1) == true);

    // test opts
    auto test2_argc = 10;
    const char* test2_argv[] = {"test2", "--int",    "10",   "--float",
                                "3.14",  "--double", "6.28", "--str",
                                "bob",   "--flag"};
    auto parser2 = yc_init_parser(test2_argc, (char**)test2_argv, "test2");
    parser2->exit_on_error = false;
    assert(yc_parse_flag(parser2, "--flag", "", "", false) == true);
    assert(yc_parse_opt<int>(parser2, "--int", "", "", 0) == 10);
    assert(fabsf(yc_parse_opt<float>(parser2, "--float", "", "", 0) - 3.14f) <
           0.01);
    assert(fabs(yc_parse_opt<double>(parser2, "--double", "", "", 0) - 6.28) <
           0.01);
    assert(yc_parse_opt<std::string>(parser2, "--str", "", "", "mike") ==
           "bob");
    assert(yc_parse_flag(parser2, "--flag_def", "", "", false) == false);
    assert(yc_parse_opt<int>(parser2, "--int_def", "", "", 5) == 5);
    assert(yc_parse_opt<float>(parser2, "--float_def", "", "", 2.67f) == 2.67f);
    assert(yc_parse_opt<double>(parser2, "--double_def", "", "", 9.54) == 9.54);
    assert(yc_parse_opt<std::string>(parser2, "--str_def", "", "", "alex") ==
           "alex");
    assert(yc_done_parser(parser2) == true);

    // test args
    auto test3_argc = 3;
    const char* test3_argv[] = {"test3", "10", "bob"};
    auto parser3 = yc_init_parser(test3_argc, (char**)test3_argv, "test3");
    parser3->exit_on_error = false;
    assert(yc_parse_arg<int>(parser3, "int", "", 0, true) == 10);
    assert(yc_parse_arg<std::string>(parser3, "str", "", "mike", true) ==
           "bob");
    assert(yc_done_parser(parser3) == true);

    // test bad opts
    auto test4_argc = 3;
    const char* test4_argv[] = {"test4", "--int", "bob"};
    auto parser4 = yc_init_parser(test4_argc, (char**)test4_argv, "test4");
    parser4->exit_on_error = false;
    assert(yc_parse_opt<int>(parser4, "--int", "", "", 0) == 0);
    assert(yc_done_parser(parser4) == false);
}

#endif

#endif
