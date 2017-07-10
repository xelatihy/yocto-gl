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

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF YOCTO_UTILS
// -----------------------------------------------------------------------------

#include "yocto_utils.h"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <condition_variable>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <functional>
#include <mutex>
#include <queue>
#include <random>
#include <sstream>
#include <thread>

//
// Thread pool code derived from LLVM codebase
//

// thread pool
struct ThreadPool {
    ThreadPool(int nthreads = std::thread::hardware_concurrency());
    ~ThreadPool();

    std::shared_future<void> async(std::function<void()> task);
    void wait();

    void _thread_proc();

    std::vector<std::thread> threads;
    std::queue<std::packaged_task<void()>> tasks;
    std::mutex queue_lock;
    std::condition_variable queue_condition;
    std::mutex completion_lock;
    std::condition_variable completion_condition;
    std::atomic<unsigned> working_threads;
    bool stop_flag = false;
};

inline ThreadPool::ThreadPool(int nthreads)
    : working_threads(0), stop_flag(false) {
    threads.reserve(nthreads);
    for (auto tid = 0; tid < nthreads; tid++) {
        threads.emplace_back([this] { _thread_proc(); });
    }
}

inline void ThreadPool::_thread_proc() {
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
            tasks.pop();
        }

        task();

        {
            std::unique_lock<std::mutex> lock_guard(completion_lock);
            working_threads--;
        }

        completion_condition.notify_all();
    }
}

inline ThreadPool::~ThreadPool() {
    {
        std::unique_lock<std::mutex> lock_guard(queue_lock);
        stop_flag = true;
    }
    queue_condition.notify_all();
    for (auto& Worker : threads) Worker.join();
}

inline std::shared_future<void> ThreadPool::async(std::function<void()> task) {
    // Wrap the Task in a packaged_task to return a future object.
    std::packaged_task<void()> packaged_task(std::move(task));
    auto future = packaged_task.get_future();
    {
        std::unique_lock<std::mutex> lock_guard(queue_lock);
        assert(!stop_flag && "Queuing a thread during ThreadPool destruction");
        tasks.push(std::move(packaged_task));
    }
    queue_condition.notify_one();
    return future.share();
}

inline void ThreadPool::wait() {
    std::unique_lock<std::mutex> lock_guard(completion_lock);
    completion_condition.wait(
        lock_guard, [&] { return tasks.empty() && !working_threads; });
}

//
// End of thread pool code derived from LLVM codebase
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
parser* make_parser(
    const std::vector<std::string>& args, const std::string& help) {
    // clears parser and copy argument data
    auto par = new parser();
    par->args = std::vector<std::string>(args.begin() + 1, args.end());
    par->error = false;
    par->exit_on_error = true;
    par->help_prog = args[0];
    par->help_usage += help;
    par->print_help = parse_flag(par, "--help", "-?", "print help", false);
    return par;
}

//
// Inits the parser.
//
parser* make_parser(int argc, char* argv[], const char* help) {
    return make_parser(std::vector<std::string>(argv, argv + argc), help);
}

//
// Print help based on the help lines collected during parsing.
//
static inline void _print_help(parser* par) {
    auto help = std::string();
    help += "usage: " + par->help_prog;
    if (!par->help_opts.empty()) help += "[options] ";
    if (!par->help_args.empty()) help += "<arguments> ";
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
bool check_parser(parser* par) {
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
static inline void _check_name(parser* par, const std::string& longname,
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
static inline void _add_help(parser* par, const std::string& longname,
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
        par->help_opts += help_line;
    else
        par->help_args += help_line;
}

//
// Parsing routine for arrays of values
//
template <typename T>
static inline std::vector<T> _parse_vals(parser* par,
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
        if (val_pos >= 0) { par->args.erase(par->args.begin() + val_pos); }
    } else {
        // check if arg is present
        if (!par->args.empty()) {
            if (par->args[0][0] != '-') { val_pos = 0; }
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
bool parse_flag(parser* par, const std::string& longname,
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
T parse_opt(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, const T& def,
    bool req, const std::vector<T>& choices) {
    return _parse_val<T>(
        par, longname, shortname, help, true, def, req, choices);
}

//
// Explicit instantiations
//
#ifndef YCMD_INLINE
template int parse_opt(parser*, const std::string&, const std::string&,
    const std::string&, const int&, bool, const std::vector<int>&);
template bool parse_opt(parser*, const std::string&, const std::string&,
    const std::string&, const bool&, bool, const std::vector<bool>&);
template float parse_opt(parser*, const std::string&, const std::string&,
    const std::string&, const float&, bool, const std::vector<float>&);
template std::string parse_opt(parser*, const std::string&, const std::string&,
    const std::string&, const std::string&, bool,
    const std::vector<std::string>&);
#endif

//
// Specialization of parse_opt()
//
std::string parse_opts(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help,
    const std::string& def, bool required,
    const std::vector<std::string>& choices) {
    return parse_opt<std::string>(
        par, longname, shortname, help, def, required, choices);
}

//
// Specialization of parse_opt()
//
int parse_opti(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, int def,
    bool required, const std::vector<int>& choices) {
    return parse_opt<int>(
        par, longname, shortname, help, def, required, choices);
}

//
// Specialization of parse_opt()
//
float parse_optf(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, float def,
    bool required, const std::vector<float>& choices) {
    return parse_opt<float>(
        par, longname, shortname, help, def, required, choices);
}

//
// Specialization of parse_opt()
//
double parse_optd(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, double def,
    bool required, const std::vector<double>& choices) {
    return parse_opt<double>(
        par, longname, shortname, help, def, required, choices);
}

//
// Specialization of parse_opt()
//
bool parse_optb(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, bool def,
    bool required) {
    return parse_opt<bool>(par, longname, shortname, help, def, required);
}

//
// Parses an argument as described in the intro.
//
template <typename T>
T parse_arg(parser* par, const std::string& longname, const std::string& help,
    const T& def, bool req, const std::vector<T>& choices) {
    return _parse_val<T>(par, longname, "", help, false, def, req, choices);
}

//
// Explicit instantiations
//
#ifndef YCMD_INLINE
template int parse_arg(parser* par, const std::string& longname,
    const std::string& help, const int& def, bool req,
    const std::vector<int>& choices);
template bool parse_arg(parser* par, const std::string& longname,
    const std::string& help, const bool& def, bool req,
    const std::vector<bool>& choices);
template float parse_arg(parser* par, const std::string& longname,
    const std::string& help, const float& def, bool req,
    const std::vector<float>& choices);
template std::string parse_arg(parser* par, const std::string& longname,
    const std::string& help, const std::string& def, bool req,
    const std::vector<std::string>& choices);
#endif

//
// Parses an argument as described in the intro.
//
std::string parse_args(parser* par, const std::string& longname,
    const std::string& help, const std::string& def, bool required,
    const std::vector<std::string>& choices) {
    return parse_arg<std::string>(par, longname, help, def, required, choices);
}

//
// Parses an argument as described in the intro.
//
int parse_argi(parser* par, const std::string& longname,
    const std::string& help, int def, bool required,
    const std::vector<int>& choices) {
    return parse_arg<int>(par, longname, help, def, required, choices);
}

//
// Parses an argument as described in the intro.
//
float parse_argf(parser* par, const std::string& longname,
    const std::string& help, float def, bool required,
    const std::vector<float>& choices) {
    return parse_arg<float>(par, longname, help, def, required, choices);
}

//
// Parses an argument as described in the intro.
//
double parse_argd(parser* par, const std::string& longname,
    const std::string& help, double def, bool required,
    const std::vector<double>& choices) {
    return parse_arg<double>(par, longname, help, def, required, choices);
}

//
// Parses an option array as described in the intro.
//
template <typename T>
std::vector<T> parse_opta(parser* par, const std::string& longname,
    const std::string& help, const std::vector<T>& def, int nargs,
    bool required, const std::vector<T>& choices) {
    return _parse_vals(
        par, longname, "", help, true, required, nargs, def, choices);
}

//
// Explicit instantiations
//
#ifndef YCMD_INLINE
template std::vector<int> parse_opta(parser* par, const std::string& longname,
    const std::string& help, const std::vector<int>& def, int nargs,
    bool required, const std::vector<int>& choices);
template std::vector<bool> parse_opta(parser* par, const std::string& longname,
    const std::string& help, const std::vector<bool>& def, int nargs,
    bool required, const std::vector<bool>& choices);
template std::vector<float> parse_opta(parser* par, const std::string& longname,
    const std::string& help, const std::vector<float>& def, int nargs,
    bool required, const std::vector<float>& choices);
template std::vector<std::string> parse_opta(parser* par,
    const std::string& longname, const std::string& help,
    const std::vector<std::string>& def, int nargs, bool required,
    const std::vector<std::string>& choices);
#endif

//
// Parses an argument array as described in the intro.
//
template <typename T>
std::vector<T> parse_arga(parser* par, const std::string& longname,
    const std::string& help, const std::vector<T>& def, int nargs,
    bool required, const std::vector<T>& choices) {
    return _parse_vals(
        par, longname, "", help, false, required, nargs, def, choices);
}

//
// Specialization of parse_arga()
//
std::vector<std::string> parse_argas(parser* par, const std::string& longname,
    const std::string& help, const std::vector<std::string>& def, int nargs,
    bool required, const std::vector<std::string>& choices) {
    return parse_arga<std::string>(
        par, longname, help, def, nargs, required, choices);
}

//
// Specialization of parse_arga()
//
std::vector<int> parse_argai(parser* par, const std::string& longname,
    const std::string& help, const std::vector<int>& def, int nargs,
    bool required, const std::vector<int>& choices) {
    return parse_arga<int>(par, longname, help, def, nargs, required, choices);
}

//
// Specialization of parse_arga()
//
std::vector<float> parse_argaf(parser* par, const std::string& longname,
    const std::string& help, const std::vector<float>& def, int nargs,
    bool required, const std::vector<float>& choices) {
    return parse_arga<float>(
        par, longname, help, def, nargs, required, choices);
}

//
// Specialization of parse_arga()
//
std::vector<double> parse_argad(parser* par, const std::string& longname,
    const std::string& help, const std::vector<double>& def, int nargs,
    bool required, const std::vector<double>& choices) {
    return parse_arga<double>(
        par, longname, help, def, nargs, required, choices);
}

//
// Parses an option enum as described in the intro.
//
int parse_opte(parser* par, const std::string& longname,
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

namespace file {
//
// Parse a binary file a chuck at a time and loads its content full in memory.
// Does not attempt to determine file size upfront to handle all cases.
//
std::vector<unsigned char> load_binfile(const std::string& filename) {
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
std::string load_txtfile(const std::string& filename) {
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

//
// Saves binary data to a file.
//
void save_binfile(
    const std::string& filename, const std::vector<unsigned char>& data) {
    auto f = fopen(filename.c_str(), "wt");
    if (!f) return;
    fwrite(data.data(), 1, data.size(), f);
    fclose(f);
}

//
// Saves a string to a text file.
//
void save_txtfile(const std::string& filename, const std::string& str) {
    auto f = fopen(filename.c_str(), "wt");
    if (!f) return;
    fwrite(str.c_str(), 1, str.length(), f);
    fclose(f);
}

}  // namespace file

namespace path {

//
// Get directory name (including '/').
//
std::string get_dirname(const std::string& filename) {
    auto pos = filename.rfind('/');
    if (pos == std::string::npos) pos = filename.rfind('\\');
    if (pos == std::string::npos) return "";
    return filename.substr(0, pos + 1);
}

//
// Get file basename.
//
std::string get_basename(const std::string& filename) {
    auto dirname = get_dirname(filename);
    auto extension = get_extension(filename);
    return filename.substr(
        dirname.size(), filename.size() - dirname.size() - extension.size());
}

//
// Get extension (including '.').
//
std::string get_extension(const std::string& filename) {
    auto pos = filename.rfind('.');
    if (pos == std::string::npos) return "";
    return filename.substr(pos);
}

//
// Get filename without directory (equiv to get_basename() + get_dirname().
//
std::string get_filename(const std::string& filename) {
    return get_basename(filename) + get_extension(filename);
}

//
// Replace extension.
//
std::string replace_extension(
    const std::string& filename, const std::string& ext) {
    return get_dirname(filename) + get_basename(filename) + ext;
}

///
/// Prepend a string to the extension.
///
std::string prepend_extension(
    const std::string& filename, const std::string& prep) {
    return get_dirname(filename) + get_basename(filename) + prep +
           get_extension(filename);
}

//
// Splits a path. Public interface.
//
void split_path(const std::string& filename, std::string& dirname,
    std::string& basename, std::string& ext) {
    dirname = get_dirname(filename);
    basename = get_basename(filename);
    ext = get_extension(filename);
}

}  // namespace path

namespace string {

//
// Checks if a std::string starts with a prefix.
//
bool startswith(const std::string& str, const std::string& substr) {
    if (str.length() < substr.length()) return false;
    for (auto i = 0; i < substr.length(); i++)
        if (str[i] != substr[i]) return false;
    return true;
}

//
// Checks if a std::string ends with a prefix.
//
bool endswith(const std::string& str, const std::string& substr) {
    if (str.length() < substr.length()) return false;
    auto offset = str.length() - substr.length();
    for (auto i = 0; i < substr.length(); i++)
        if (str[i + offset] != substr[i]) return false;
    return true;
}

//
// Check is a string contains a substring.
//
bool contains(const std::string& str, const std::string& substr) {
    return str.find(substr) != str.npos;
}

//
// Splits a std::string into lines at the '\n' character.
//
std::vector<std::string> splitlines(const std::string& str, bool keep_newline) {
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

//
// Partition the string.
//
std::vector<std::string> partition(
    const std::string& str, const std::string& split) {
    auto pos = str.find(split);
    if (pos == str.npos) return {str, "", ""};
    return {str.substr(0, pos), split, str.substr(pos + split.length())};
}

//
// Splits the string.
//
std::vector<std::string> split(const std::string& str) {
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

//
// Strip the string.
//
std::string rstrip(const std::string& str) {
    auto pos = str.find_last_not_of(" \t\r\n");
    if (pos == str.npos) return "";
    return str.substr(0, pos + 1);
}

//
// Strip the string.
//
std::string lstrip(const std::string& str) {
    auto pos = str.find_first_not_of(" \t\r\n");
    if (pos == str.npos) return "";
    return str.substr(pos);
}

//
// Strip the string.
//
std::string strip(const std::string& str) { return rstrip(lstrip(str)); }

//
// Joins a list of std::string with a std::string as separator.
//
std::string join(const std::vector<std::string>& strs, const std::string& sep) {
    auto ret = std::string();
    auto first = true;
    for (auto& str : strs) {
        if (!first) ret += sep;
        ret += str;
        first = false;
    }
    return ret;
}

//
// Converts an ASCII string to lowercase.
//
std::string lower(const std::string& str) {
    auto s = str;
    for (auto& c : s) c = std::tolower(c);
    return s;
}

//
// Converts an ASCII string to uppercase.
//
std::string upper(const std::string& str) {
    auto s = str;
    for (auto& c : s) c = std::toupper(c);
    return s;
}

//
// Strung is space.
//
bool isspace(const std::string& str) {
    for (auto c : str) {
        if (c != ' ' && c != '\n' && c != '\t' && c != '\r') return false;
    }
    return true;
}

//
// Replace s1 with s2 in str.
//
std::string replace(
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
// C-like string formatting
//
std::string format(const char* fmt, va_list args) {
    char buffer[1024 * 16];
    vsprintf(buffer, fmt, args);
    return buffer;
}

//
// C-like string formatting
//
std::string format(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    auto s = format(fmt, args);
    va_end(args);
    return s;
}

}  // namespace string

namespace logging {

//
// Logger
//
struct logger {
    FILE* file = nullptr;
    log_level output_level = log_level::info;
    log_level flush_level = log_level::error;
    unsigned int guid = std::random_device()();

    ~logger() {
        if (file == stderr) return;
        if (file == stdout) return;
        if (file) {
            fclose(file);
            file = nullptr;
        }
    }
};

//
// Skip a message
//
static inline bool _log_skip(logger* lgr, log_level level) {
    return level < lgr->output_level;
}

//
// Get default logger
//
std::vector<logger*>* get_default_loggers() {
    static std::vector<logger*> lgrs = {};
    return &lgrs;
}

//
// Create a stream logger
//
static inline logger* _make_stream_logger(
    FILE* file, log_level output_level, log_level flush_level) {
    auto lgr = new logger();
    lgr->file = file;
    set_logger(lgr, output_level, flush_level);
    return lgr;
}

//
// Create a file logger
//
logger* make_file_logger(const std::string& filename, bool append,
    log_level output_level, log_level flush_level) {
    auto file = fopen(filename.c_str(), (append) ? "at" : "wt");
    if (!file) { throw std::runtime_error("cannot open log file " + filename); }
    return _make_stream_logger(file, output_level, flush_level);
}

//
// Create a stderr logger
//
logger* make_stderr_logger(log_level output_level, log_level flush_level) {
    return _make_stream_logger(stderr, output_level, flush_level);
}

//
// Create a stderr logger
//
logger* make_stdout_logger(log_level output_level, log_level flush_level) {
    return _make_stream_logger(stdout, output_level, flush_level);
}

//
// Clear a logger
//
void clear_logger(logger* lgr) {
    if (lgr) delete lgr;
}

//
// Set logger level
//
void set_logger(logger* lgr, log_level output_level, log_level flush_level) {
    lgr->output_level = output_level;
    lgr->flush_level = flush_level;
}

//
// Log a message
//
void log_msg(logger* lgr, log_level level, const std::string& name,
    const std::string& msg) {
    log_msg(lgr, level, name.c_str(), msg.c_str());
}

//
// Log a message
//
void log_msg(logger* lgr, log_level level, const char* name, const char* msg) {
    // skip if not needed
    if (_log_skip(lgr, level)) return;

    // type string
    const char* types[] = {"VERB", "INFO", "WARN", "ERRN"};
    const char* type = types[std::max(0, std::min(3, (int)level + 1))];

    // time string
    char time_buf[1024];
    auto tm = time(nullptr);
    auto ttm = localtime(&tm);  // TODO: use thread safe version
    strftime(time_buf, 1024, "%Y-%m-%d %H:%M:%S", ttm);

    // output message
    fprintf(lgr->file, "%s %s %4x %-16s %s\n", time_buf, type, lgr->guid, name,
        msg);

    // flush if needed
    if (level < lgr->flush_level) return;
    fflush(lgr->file);
}

//
// Log a message formatted
//
void log_msgfv(logger* lgr, log_level level, const char* name, const char* msg,
    va_list args) {
    // skip if not needed
    if (_log_skip(lgr, level)) return;

    // make message
    char msg_buf[1024 * 16];
    vsprintf(msg_buf, msg, args);

    // log
    log_msg(lgr, level, name, msg_buf);
}

//
// Log a message formatted
//
void log_msgf(
    logger* lgr, log_level level, const char* name, const char* msg, ...) {
    // skip if not needed
    if (_log_skip(lgr, level)) return;

    // make message
    va_list args;
    va_start(args, msg);
    log_msgfv(lgr, level, name, msg, args);
    va_end(args);
}

//
// Log a message to the default logger
//
void log_msg(log_level level, const char* name, const char* msg) {
    for (auto lgr : *get_default_loggers()) { log_msg(lgr, level, name, msg); }
}

//
// Log a message to the default logger
//
void log_msgfv(
    log_level level, const char* name, const char* msg, va_list args) {
    for (auto lgr : *get_default_loggers()) {
        log_msgfv(lgr, level, name, msg, args);
    }
}

//
// Log a message to the default logger
//
void log_msgf(log_level level, const char* name, const char* msg, ...) {
    for (auto lgr : *get_default_loggers()) {
        // skip if not needed
        if (_log_skip(lgr, level)) return;
        va_list args;
        va_start(args, msg);
        log_msgfv(lgr, level, name, msg, args);
        va_end(args);
    }
}

//
// Log a message to the default logger
//
void log_msg(log_level level, const std::string& name, const std::string& msg) {
    for (auto lgr : *get_default_loggers()) {
        log_msg(lgr, level, name.c_str(), msg.c_str());
    }
}

}  // namespace logging

namespace concurrent {

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
thread_pool* make_pool(int nthread) {
    if (!nthread) nthread = std::thread::hardware_concurrency();
    auto pool = new thread_pool();
    pool->tp = new ThreadPool(nthread);
    return pool;
}

//
// Clear thread pool
//
void clear_pool(thread_pool* pool) {
    if (pool) delete pool;
}

//
// Enqueue a job
//
std::shared_future<void> run_async(
    thread_pool* pool, const std::function<void()>& task) {
    return pool->tp->async(task);
}

//
// Wait for jobs to finish
//
void wait_pool(thread_pool* pool) { pool->tp->wait(); }

//
// Parallel for implementation
//
void parallel_for(
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
std::shared_future<void> run_async(const std::function<void()>& task) {
    if (!global_pool) make_global_thread_pool();
    return run_async(global_pool, task);
}

//
// Wait for jobs to finish
//
void wait_pool() {
    if (global_pool) global_pool->tp->wait();
}

//
// Parallel for implementation
//
void parallel_for(int count, const std::function<void(int idx)>& task) {
    if (!global_pool) make_global_thread_pool();
    parallel_for(global_pool, count, task);
}

}  // namespace concurrent

}  // namespace yu

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
