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
// 4. immediate mode command line parsing
//     1. create a parser with `make_cmdline()`
//     2. for each option, parse it calling the functions `parse_opt()`
//         - options are parsed on the fly and a comprehensive help is
//           automatically generated
//         - supports bool (flags), int, float, double, string, enums, vecXX
//         - options names are "--longname" for longname and "-s" for short
//         - command line format is "--longname value", "-s v" for all but flags
//         - for general use `opt = parse_opt<type>()`
//         - for boolean flags is `parse_flag()`
//     3. for each unnamed argument, parse it calling the functions parse_arg()
//         - names are only used for help
//         - supports types as above
//         - for general use `arg = parse_arg<type>()`
//         - to parse all remaining values use `args = parse_args<type>(...)`
//     4. end cmdline parsing with `check_parsing()` to check for unused values,
//        missing arguments
//     5. to check for error use `should_exit()` and to print the message use
//        `get_message()`
//     6. since arguments are parsed immediately, one can easily implement
//        subcommands by just branching the command line code based on a read
//        argument without any need for complex syntax
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
inline std::string format_value(const vec2i& val) {
    return std::to_string(val.x) + " " + std::to_string(val.y);
}
inline std::string format_value(const vec3i& val) {
    return std::to_string(val.x) + " " + std::to_string(val.y) + " " +
           std::to_string(val.z);
}
inline std::string format_value(const vec4i& val) {
    return std::to_string(val.x) + " " + std::to_string(val.y) + " " +
           std::to_string(val.z) + " " + std::to_string(val.w);
}
inline std::string format_value(const float& val) {
    return std::to_string(val);
}
inline std::string format_value(const vec2f& val) {
    return std::to_string(val.x) + " " + std::to_string(val.y);
}
inline std::string format_value(const vec3f& val) {
    return std::to_string(val.x) + " " + std::to_string(val.y) + " " +
           std::to_string(val.z);
}
inline std::string format_value(const vec4f& val) {
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

// Gets/sets log verbosity and output stream
inline bool& log_verbose() {
    static bool verbose = true;
    return verbose;
}

// Implementation for logging functions.
inline void _log_msg(const std::string& msg, const char* tag) {
    using namespace std::chrono;
    auto time = system_clock::to_time_t(system_clock::now());
    auto tm = std::localtime(&time);
    char tmstr[64];
    strftime(tmstr, 64, "%T", tm);
    printf("%s %s %s\n", tmstr, tag, msg.c_str());
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

// -----------------------------------------------------------------------------
// IMMEDIATE MODE COMMAND LINE PARSER
// -----------------------------------------------------------------------------
namespace ygl {

// Immediate mode command line parser. Members are not part of the public API.
struct cmdline_parser {
    std::vector<std::string> _to_parse;    // args left to parse
    std::vector<std::string> _used_names;  // used names for check
    std::string _usage_prog;               // usage prog line
    std::string _usage_help;               // usage help line
    std::string _usage_opts;               // usage option lines
    std::string _usage_args;               // usage argument lines
    bool _usage = false;                   // help option triggered
    std::string _error;                    // parse error
};

// Initialize a command line parser.
inline cmdline_parser make_parser(
    int argc, char** argv, const std::string& prog, const std::string& help);

// Check unused arguments.
inline bool should_exit(cmdline_parser& parser);

// Returns the usage string.
inline std::string get_usage(const cmdline_parser& parser);

// Pase a flag from the command line.
inline bool parse_flag(cmdline_parser& parser, const std::string& name,
    const std::string& flag, const std::string& help, bool def = false,
    bool req = false);

// Pase an option from the command line.
template <typename T>
inline T parse_opt(cmdline_parser& parser, const std::string& name,
    const std::string& flag, const std::string& help, const T& def = {},
    bool req = false, const std::vector<T>& choices = {});

// Parse an enum option from the command line.
template <typename T>
inline T parse_opt(cmdline_parser& parser, const std::string& name,
    const std::string& flag, const std::string& help,
    const std::map<T, std::string>& key_values, const T& def, bool req = false,
    const std::vector<T>& choices = {});

// Parse positional argument from the command line.
template <typename T>
inline T parse_arg(cmdline_parser& parser, const std::string& name,
    const std::string& help, const T& def = {}, bool req = true,
    const std::vector<T>& choices = {});

// Parse all remaining positional argument from the command line.
template <typename T>
inline std::vector<T> parse_args(cmdline_parser& parser,
    const std::string& name, const std::string& help,
    const std::vector<T>& def = {}, bool req = true,
    const std::vector<T>& choices = {});

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMMEDIATE MODE COMMAND LINE PARSER
// -----------------------------------------------------------------------------
namespace ygl {

// cmdline implementation
inline void _check_name(cmdline_parser& parser, const std::string& name,
    const std::string& flag, bool opt) {
    if (opt) {
        if (name.size() < 3 || name[0] != '-' || name[1] != '-' ||
            name[2] == '-')
            throw std::runtime_error("bad name " + name);
    } else {
        if (name.size() < 1 || name[0] == '-')
            throw std::runtime_error("bad name " + name);
    }
    if (find(parser._used_names.begin(), parser._used_names.end(), name) !=
        parser._used_names.end())
        throw std::runtime_error("already used " + name);
    parser._used_names.push_back(name);
    if (flag.empty()) return;
    if (flag.size() < 2 || flag[0] != '-' || flag[1] == '-')
        throw std::runtime_error("bad name " + flag);
    if (find(parser._used_names.begin(), parser._used_names.end(), flag) !=
        parser._used_names.end())
        throw std::runtime_error("already used " + flag);
    parser._used_names.push_back(flag);
}

// cmdline implementation
template <typename T>
inline void _add_usage_str(cmdline_parser& parser, const std::string& name,
    const std::string& flag, bool opt, const std::string& metavar,
    const std::string& help, const std::string& def, bool req,
    const std::vector<T>& choices) {
    auto str = std::string();
    str += "  " + name;
    if (!flag.empty()) str += "/" + flag;
    if (!metavar.empty()) str += " " + metavar;
    while (str.length() < 32) str += " ";
    str += help + " ";
    if (!req && !def.empty()) str += "[" + def + "]";
    if (req) str += "(required)";
    str += "\n";
    if (!choices.empty()) {
        for (auto i = 0; i < 32; i++) str += " ";
        str += "(";
        auto first = true;
        for (auto&& c : choices) {
            if (!first) str += ",";
            str += format_value(c);
            first = false;
        }
        str += ")";
        str += "\n";
    }
    if (opt)
        parser._usage_opts += str;
    else
        parser._usage_args += str;
}

// cmdline implementation
inline void _set_error(cmdline_parser& parser, const std::string& err) {
    if (parser._error.empty()) parser._error = err;
}

// Initialize a command line parser.
inline cmdline_parser make_parser(
    int argc, char** argv, const std::string& prog, const std::string& help) {
    auto parser = cmdline_parser();
    parser._to_parse = std::vector<std::string>(argv + 1, argv + argc);
    parser._usage_prog = (prog.empty()) ? std::string(argv[0]) : prog;
    parser._usage_help = help;
    parser._usage =
        parse_flag(parser, "--help", "-?", "prints and help message");
    return parser;
}

// Check unused arguments.
inline bool should_exit(cmdline_parser& parser) {
    for (auto&& v : parser._to_parse) {
        if (v[0] == '-')
            _set_error(parser, "unknown option " + v);
        else
            _set_error(parser, "unknown argument " + v);
    }
    return parser._usage || !parser._error.empty();
}

// Returns the usage string.
inline std::string get_usage(const cmdline_parser& parser) {
    auto str = std::string();
    if (!parser._error.empty()) str += "error: " + parser._error + "\n\n";
    str += parser._usage_prog;
    if (!parser._usage_opts.empty()) str += " [options]";
    if (!parser._usage_args.empty()) str += " <arguments>";
    str += "\n";
    // while (str.size() < 32) str += " ";
    str += parser._usage_help + "\n\n";
    if (!parser._usage_opts.empty())
        str += "options:\n" + parser._usage_opts + "\n";
    if (!parser._usage_args.empty())
        str += "arguments:\n" + parser._usage_args + "\n";
    return str;
}

// Pase a flag from the command line.
inline bool parse_flag(cmdline_parser& parser, const std::string& name,
    const std::string& flag, const std::string& help, bool def, bool req) {
    // check names
    _check_name(parser, name, flag, true);
    // update usage
    _add_usage_str(
        parser, name, flag, true, "", help, "", req, std::vector<bool>{});
    // skip if error
    if (!parser._error.empty()) return def;
    // find location of option
    auto pos = find(parser._to_parse.begin(), parser._to_parse.end(), name);
    if (pos == parser._to_parse.end())
        pos = find(parser._to_parse.begin(), parser._to_parse.end(), flag);
    if (pos == parser._to_parse.end()) {
        if (req) _set_error(parser, "missing required flag " + name);
        return def;
    }
    // remove parsed arg
    parser._to_parse.erase(pos, pos + 1);
    // done
    return !def;
}

// Parse a value
inline bool parse_value(const std::string& str, int& val) {
    return sscanf(str.c_str(), "%d", &val) == 1;
}
inline bool parse_value(const std::string& str, float& val) {
    return sscanf(str.c_str(), "%f", &val) == 1;
}
inline bool parse_value(const std::string& str, vec2f& val) {
    return sscanf(str.c_str(), "%f%f", &val.x, &val.y) == 2;
}
inline bool parse_value(const std::string& str, vec3f& val) {
    return sscanf(str.c_str(), "%f%f%f", &val.x, &val.y, &val.z) == 3;
}
inline bool parse_value(const std::string& str, vec4f& val) {
    return sscanf(str.c_str(), "%f%f%f%f", &val.x, &val.y, &val.z, &val.w) == 4;
}
inline bool parse_value(const std::string& str, std::string& val) {
    val = str;
    return true;
}
inline bool parse_value(const std::string& str, bool& val) {
    if (str == "true") {
        val = true;
        return true;
    }
    if (str == "false") {
        val = false;
        return true;
    }
    return false;
}

// Pase an option from the command line.
template <typename T>
inline T parse_opt(cmdline_parser& parser, const std::string& name,
    const std::string& flag, const std::string& help, const T& def, bool req,
    const std::vector<T>& choices) {
    // check names
    _check_name(parser, name, flag, true);
    // update usage
    _add_usage_str(parser, name, flag, true, "<val>", help, format_value(def),
        req, choices);
    // skip if error
    if (!parser._error.empty()) return def;
    // find location of option
    auto pos = find(parser._to_parse.begin(), parser._to_parse.end(), name);
    if (pos == parser._to_parse.end())
        pos = find(parser._to_parse.begin(), parser._to_parse.end(), flag);
    if (pos == parser._to_parse.end()) {
        if (req) _set_error(parser, "missing option " + name);
        return def;
    }
    // check if value exists
    if (pos == parser._to_parse.end() - 1) {
        _set_error(parser, "no value for parameter " + name);
        return def;
    }
    // get value
    auto val = def;
    const auto& arg = *(pos + 1);
    // parse
    if (!parse_value(arg, val)) {
        _set_error(
            parser, "incorrect value \"" + arg + "\" for option " + name);
        return def;
    }
    // validate if necessary
    if (!choices.empty()) {
        if (find(choices.begin(), choices.end(), val) == choices.end())
            _set_error(
                parser, "incorrect value \"" + arg + "\" for option " + name);
    }
    // remove parsed arg
    parser._to_parse.erase(pos, pos + 2);
    // done
    return val;
}

// Parse an enum option from the command line.
template <typename T>
inline T parse_opt(cmdline_parser& parser, const std::string& name,
    const std::string& flag, const std::string& help,
    const std::map<T, std::string>& key_values, const T& def, bool req,
    const std::vector<T>& choices) {
    auto keys = std::vector<std::string>{};
    auto key_def = key_values.at(def);
    for (auto&& kv : key_values) keys.push_back(kv.second);
    auto key =
        parse_opt<std::string>(parser, name, flag, help, key_def, req, keys);
    if (!parser._error.empty()) return def;
    auto val = def;
    for (auto&& kv : key_values) {
        if (kv.second == key) val = kv.first;
    }
    return val;
}

// Parse positional argument from the command line.
template <typename T>
inline T parse_arg(cmdline_parser& parser, const std::string& name,
    const std::string& help, const T& def, bool req,
    const std::vector<T>& choices) {
    // check names
    _check_name(parser, name, "", false);
    // update usage
    _add_usage_str(
        parser, name, "", false, "", help, format_value(def), req, choices);
    // skip if error
    if (!parser._error.empty()) return def;
    // find location of argument
    auto pos = std::find_if(parser._to_parse.begin(), parser._to_parse.end(),
        [](const auto& s) { return s.size() > 0 && s[0] != '-'; });
    if (pos == parser._to_parse.end()) {
        if (req) _set_error(parser, "missing argument " + name);
        return def;
    }
    // get value
    auto val = def;
    const auto& arg = *(pos);
    // parse
    if (!parse_value(arg, val)) {
        _set_error(
            parser, "incorrect value \"" + arg + "\" for argument " + name);
    }
    // validate if necessary
    if (!choices.empty()) {
        if (find(choices.begin(), choices.end(), val) == choices.end())
            _set_error(
                parser, "incorrect value \"" + arg + "\" for argument " + name);
    }
    // remove parsed arg
    parser._to_parse.erase(pos, pos + 1);
    // done
    return val;
}

// Parse all remaining positional argument from the command line.
template <typename T>
inline std::vector<T> parse_args(cmdline_parser& parser,
    const std::string& name, const std::string& help, const std::vector<T>& def,
    bool req, const std::vector<T>& choices) {
    // check names
    _check_name(parser, name, "", false);
    // update usage
    _add_usage_str(parser, name, "", false, "", help, "", req, choices);
    // skip if error
    if (!parser._error.empty()) return def;
    // search for all params
    auto vals = std::vector<T>();
    while (true) {
        // find location of argument
        auto pos =
            std::find_if(parser._to_parse.begin(), parser._to_parse.end(),
                [](const auto& s) { return s.size() > 0 && s[0] != '-'; });
        if (pos == parser._to_parse.end()) break;
        // get value
        auto val = T{};
        const auto& arg = *(pos);
        // parse
        if (!parse_value(arg, val)) {
            _set_error(
                parser, "incorrect value \"" + arg + "\" for argument " + name);
        }
        // validate if necessary
        if (!choices.empty()) {
            if (find(choices.begin(), choices.end(), val) == choices.end())
                _set_error(parser,
                    "incorrect value \"" + arg + "\" for argument " + name);
        }
        // remove parsed arg
        parser._to_parse.erase(pos, pos + 1);
        // append value
        vals.push_back(val);
    }
    // check missing
    if (vals.empty()) {
        if (req) _set_error(parser, "missing argument " + name);
        return def;
    }
    // done
    return vals;
}

}  // namespace ygl

#endif
