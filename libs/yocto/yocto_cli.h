//
// # Yocto/CLI: Utilities for writing command-line apps
//
// Yocto/CLI is a collection of utilities used in writing command-line
// applications, including parsing command line arguments, printing values,
// timers, progress bars, handling errors via exceptions.
// Yocto/CLI is implemented in `yocto_cli.h` and `yocto_cli.cpp`.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
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

#ifndef _YOCTO_CLI_H_
#define _YOCTO_CLI_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <array>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "yocto_math.h"

#define JSON_USE_IMPLICIT_CONVERSIONS 0
#include "ext/json.hpp"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::string;
using std::vector;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// PRINT/LOG/FORMATTING UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Print a message to the console
void print_info(const string& message);
void print_info(const char* message);
// Prints a message to the console and exit with an error. Returns error code.
void print_fatal(const string& message);
void print_fatal(const char* message);
template <typename T>
inline void print_fatal(const T& error) {
  print_fatal(error.error);
}

// Timer that prints as scope end. Create with `print_timed` and print with
// `print_elapsed`.
struct print_timer {
  int64_t start_time = -1;
  ~print_timer();  // print time if scope ends
};
// Print traces for timing and program debugging
print_timer print_timed(const string& message);
int64_t     print_elapsed(print_timer& timer);

// Print progress
void print_progress(const string& message, int current, int total);
void print_progress_begin(const string& message, int total = 1);
void print_progress_end();
void print_progress_next();

// Format duration string from nanoseconds
string format_duration(int64_t duration);
// Format a large integer number in human readable form
string format_num(uint64_t num);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIMPLE TIMER
// -----------------------------------------------------------------------------
namespace yocto {

// Simple timer
struct simple_timer {
  int64_t start = -1, stop = -1;
  simple_timer();
};

// Timer opreations
void    start_timer(simple_timer& timer);
void    stop_timer(simple_timer& timer);
int64_t elapsed_nanoseconds(simple_timer& timer);
double  elapsed_seconds(simple_timer& timer);
string  elapsed_formatted(simple_timer& timer);

}  // namespace yocto

// -----------------------------------------------------------------------------
// ERROR HANDLING VIA EXCEPTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// handle errors
template <typename Func>
inline void handle_errors(Func&& run);
inline void handle_errors(
    void (*run)(const vector<string>&), const vector<string>& args);

}  // namespace yocto

// -----------------------------------------------------------------------------
// COMMAND LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// Initialize a command line parser.
struct cli_state;
cli_state make_cli(const string& cmd, const string& usage);
// parse arguments, checks for errors
struct cli_status {
  string   error = "";
  explicit operator bool() const { return error.empty(); }
};
cli_status parse_cli(cli_state& cli, const vector<string>& args);
// parse arguments, checks for errors
bool parse_cli(cli_state& cli, const vector<string>& args, string& error);
// a convenience function that packs args to strings
vector<string> make_cli_args(int argc, const char** argv);

// get usage
string get_usage(const cli_state& cli);
// get help
bool get_help(const cli_state& cli);

// Add a subcommand
struct cli_command;
cli_command add_command(
    const cli_command& cli, const string& name, const string& usage);
void                set_command_var(const cli_command& cli, string& value);
void                set_help_var(const cli_command& cli, bool& value);
[[deprecated]] void add_command_name(const cli_command& cli, const string& name,
    string& value, const string& usage);

// Add a command with variables. Calls the function add_options for the type.
template <typename T>
inline cli_state make_cli(const string& name, T& values, const string& usage);
template <typename T>
inline cli_command add_command(
    const cli_command& cli, const string& name, T& values, const string& usage);

// Add an optional argument. Supports strings, numbers, and boolean flags.
// Optional arguments will be parsed with name `--<name>` and `-<alt>`.
// Optional booleans will support both `--<name>` and `--no-<name>` to enabled
// and disable the flag.
void add_option(const cli_command& cli, const string& name, int& value,
    const string& usage, const vector<int>& minmax = {}, const string& alt = "",
    bool req = false);
void add_option(const cli_command& cli, const string& name, float& value,
    const string& usage, const vector<float>& minmax = {},
    const string& alt = "", bool req = false);
void add_option(const cli_command& cli, const string& name, bool& value,
    const string& usage, const vector<string>& choices = {},
    const string& alt = "", bool req = false);
void add_option(const cli_command& cli, const string& name, string& value,
    const string& usage, const vector<string>& choices = {},
    const string& alt = "", bool req = false);
void add_option_with_config(const cli_command& cli, const string& name,
    string& value, const string& usage, const string& config,
    const string& alt = "", bool req = false);
// Add a positional argument. Supports strings, numbers, and boolean flags.
void add_argument(const cli_command& cli, const string& name, int& value,
    const string& usage, const vector<int>& minmax = {}, bool req = true);
void add_argument(const cli_command& cli, const string& name, float& value,
    const string& usage, const vector<float>& minmax = {}, bool req = true);
void add_argument(const cli_command& cli, const string& name, bool& value,
    const string& usage, const vector<string>& choices = {}, bool req = true);
void add_argument(const cli_command& cli, const string& name, string& value,
    const string& usage, const vector<string>& choices = {}, bool req = true);
void add_argument_with_config(const cli_command& cli, const string& name,
    string& value, const string& usage, const string& config, bool req = true);
// Add an optional argument with values as labels. Supports integers, enums
// and strings.
void add_option(const cli_command& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, const string& alt = "",
    bool req = false);
template <typename T, typename = std::enable_if_t<std::is_enum_v<T>>>
inline void add_option(const cli_command& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, const string& alt = "",
    bool req = false);
// Add a positional argument with values as labels. Supports string, integers
// and enums.
void add_argument(const cli_command& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, bool req = true);
template <typename T, typename = std::enable_if_t<std::is_enum_v<T>>>
inline void add_argument(const cli_command& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, bool req = true);
// Add a positional argument that consumes all arguments left.
// Supports strings and enums.
void add_argument(const cli_command& cli, const string& name,
    vector<int>& value, const string& usage, const vector<int>& minmax,
    bool req = true);
void add_argument(const cli_command& cli, const string& name,
    vector<float>& value, const string& usage, const vector<float>& minmax,
    bool req = true);
void add_argument(const cli_command& cli, const string& name,
    vector<int>& value, const string& usage, const vector<string>& choices = {},
    bool req = true);
void add_argument(const cli_command& cli, const string& name,
    vector<string>& value, const string& usage,
    const vector<string>& choices = {}, bool req = true);

// Add an optional argument. Supports basic math types.
void add_option(const cli_command& cli, const string& name, vec2i& value,
    const string& usage, const vector<int>& minmax = {}, const string& alt = "",
    bool req = false);
void add_option(const cli_command& cli, const string& name, vec3i& value,
    const string& usage, const vector<int>& minmax = {}, const string& alt = "",
    bool req = false);
void add_option(const cli_command& cli, const string& name, vec4i& value,
    const string& usage, const vector<int>& minmax = {}, const string& alt = "",
    bool req = false);
void add_option(const cli_command& cli, const string& name, vec2f& value,
    const string& usage, const vector<float>& minmax = {},
    const string& alt = "", bool req = false);
void add_option(const cli_command& cli, const string& name, vec3f& value,
    const string& usage, const vector<float>& minmax = {},
    const string& alt = "", bool req = false);
void add_option(const cli_command& cli, const string& name, vec4f& value,
    const string& usage, const vector<float>& minmax = {},
    const string& alt = "", bool req = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF COMMAND-LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// Json value
using json_value_ = nlohmann::ordered_json;

// Command line setter.
using cli_setter = bool (*)(
    const json_value_&, void*, const vector<string>& choices);
// Command line variable.
struct cli_variable {
  void*                                       value     = nullptr;
  cli_setter                                  setter    = nullptr;
  vector<string>                              choices   = {};
  nlohmann::ordered_map<string, cli_variable> variables = {};
};
// Command line state.
struct cli_state {
  json_value_  defaults  = {};
  json_value_  schema    = {};
  json_value_  value     = {};
  cli_variable variables = {};
};
// Command line command.
struct cli_command {
  cli_state* state = nullptr;
  string     path  = "";
  cli_command() {}
  cli_command(cli_state* state, const string& path)
      : state{state}, path{path} {}
  cli_command(const cli_state& state) : state{(cli_state*)&state}, path{""} {}
};

// Command line error.
struct cli_error : std::runtime_error {
  cli_error(const string& message) : std::runtime_error(message), _usage() {}
  cli_error(const string& message, const string& usage)
      : std::runtime_error(message), _usage(usage) {}

  string usage() const { return _usage; }
  void   set_usage(const string& usage) { _usage = usage; }

 private:
  string _usage = "";
};
struct cli_help : std::runtime_error {
  cli_help(const string& message) : std::runtime_error(message), _usage() {}
  cli_help(const string& message, const string& usage)
      : std::runtime_error(message), _usage(usage) {}

  string usage() const { return _usage; }
  void   set_usage(const string& usage) { _usage = usage; }

 private:
  string _usage = "";
};

// Add a command with variables. Calls the function add_options for the type.
template <typename T>
inline cli_state make_cli(const string& name, T& values, const string& usage) {
  auto cli = make_cli(name, usage);
  add_options(cli, values);
  return cli;
}
template <typename T>
inline cli_command add_command(const cli_command& cli, const string& name,
    T& values, const string& usage) {
  auto cmd = add_command(cli, name, usage);
  add_options(cmd, values);
  return cli;
}

template <typename T, typename>
inline void add_option(const cli_command& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, const string& alt,
    bool req) {
  return add_option(
      cli, name, (std::underlying_type_t<T>&)value, usage, choices, alt, req);
}
template <typename T, typename>
inline void add_argument(const cli_command& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, bool req) {
  return add_argument(
      cli, name, (std::underlying_type_t<T>&)value, usage, choices, req);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// ERROR HANDLING VIA EXCEPTIONS
// -----------------------------------------------------------------------------
namespace yocto {

template <typename Func>
inline void handle_errors(Func&& run) {
  try {
    run();
  } catch (const cli_error& error) {
    print_info("");
    print_info("error: " + string{error.what()});
    print_info("");
    print_info(error.usage());
    exit(1);
  } catch (const cli_help& error) {
    print_info("");
    print_info(error.usage());
    exit(0);
  } catch (const std::exception& error) {
    print_info("");
    print_info(error.what());
    exit(1);
  }
}

// handle errors
inline void handle_errors(
    void (*run)(const vector<string>&), const vector<string>& args) {
  handle_errors([&]() { run(args); });
}

}  // namespace yocto

#endif
