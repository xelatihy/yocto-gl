//
// # Yocto/CommonIO: Utilities for writing command-line apps
//
// Yocto/CommonIO is a collection of utilities used in writing command-line
// applications, including parsing command line arguments, simple path
// manipulation, file lading and saving, and printing values, timers and
// progress bars.
// Yocto/CommonIO is implemented in `yocto_commonio.h` and `yocto_commonio.cpp`.
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

#ifndef _YOCTO_COMMONIO_H_
#define _YOCTO_COMMONIO_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <chrono>
#include <cstdio>
#include <functional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <vector>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::function;
using std::string;
using std::unordered_set;
using std::vector;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// PRINT/FORMATTING UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Print a message to the console
inline void print_info(const string& msg);
// Prints a message to the console and exit with an error. Returns error code.
inline int print_fatal(const string& msg);

// Timer that prints as scope end. Create with `print_timed` and print with
// `print_elapsed`.
struct print_timer {
  int64_t start_time = -1;
  ~print_timer();  // print time if scope ends
};
// Print traces for timing and program debugging
inline print_timer print_timed(const string& msg);
inline int64_t     print_elapsed(print_timer& timer);

// Print progress
inline void print_progress(const string& message, int current, int total);

// Format duration string from nanoseconds
inline string format_duration(int64_t duration);
// Format a large integer number in human readable form
inline string format_num(uint64_t num);

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
inline void    start_timer(simple_timer& timer);
inline void    stop_timer(simple_timer& timer);
inline int64_t elapsed_nanoseconds(simple_timer& timer);
inline double  elapsed_seconds(simple_timer& timer);
inline string  elapsed_formatted(simple_timer& timer);

}  // namespace yocto

// -----------------------------------------------------------------------------
// COMMAND LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// Initialize a command line parser.
struct cli_command;
inline cli_command make_cli(const string& cmd, const string& usage);
// parse arguments, checks for errors, and exits on error or help
inline void parse_cli(cli_command& cli, int argc, const char** argv);
// parse arguments and checks for errors
inline bool parse_cli(
    cli_command& cli, int argc, const char** argv, string& error);
// gets usage message
inline string get_usage(const cli_command& cli);
// gets whether help was invoked
inline bool get_help(const cli_command& cli);
// gets the set command
inline string get_command(const cli_command& cli);

// Add an optional argument. Supports strings, numbers, and boolean flags.
// Optional arguments will be parsed with name `--<name>` and `-<alt>`.
// Optional booleans will support both `--<name>` and `--no-<name>` to enabled
// and disable the flag.
inline void add_optional(cli_command& cli, const string& name, int& value,
    const string& usage, const vector<int>& minmax = {}, const string& alt = "",
    bool req = false);
inline void add_optional(cli_command& cli, const string& name, float& value,
    const string& usage, const vector<float>& minmax = {},
    const string& alt = "", bool req = false);
inline void add_optional(cli_command& cli, const string& name, bool& value,
    const string& usage, const vector<string>& choices = {},
    const string& alt = "", bool req = false);
inline void add_optional(cli_command& cli, const string& name, string& value,
    const string& usage, const vector<string>& choices = {},
    const string& alt = "", bool req = false);
// Add a positional argument. Supports strings, numbers, and boolean flags.
inline void add_positional(cli_command& cli, const string& name, int& value,
    const string& usage, const vector<int>& minmax = {}, bool req = true);
inline void add_positional(cli_command& cli, const string& name, float& value,
    const string& usage, const vector<float>& minmax = {}, bool req = true);
inline void add_positional(cli_command& cli, const string& name, bool& value,
    const string& usage, const vector<string>& choices = {}, bool req = true);
inline void add_positional(cli_command& cli, const string& name, string& value,
    const string& usage, const vector<string>& choices = {}, bool req = true);
// Add an optional argument with values as labels. Supports integers, enums and
// strings.
inline void add_optional(cli_command& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, const string& alt = "",
    bool req = false);
template <typename T, typename = std::enable_if_t<std::is_enum_v<T>>>
inline void add_optional(cli_command& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, const string& alt = "",
    bool req = false);
// Add a positional argument with values as labels. Supports string, integers
// and enums.
inline void add_positional(cli_command& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, bool req = true);
template <typename T, typename = std::enable_if_t<std::is_enum_v<T>>>
inline void add_positional(cli_command& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, bool req = true);
// Add a positional argument that consumes all arguments left.
// Supports strings and enums.
inline void add_positional(cli_command& cli, const string& name,
    vector<int>& value, const string& usage, const vector<int>& minmax,
    bool req = true);
inline void add_positional(cli_command& cli, const string& name,
    vector<float>& value, const string& usage, const vector<float>& minmax,
    bool req = true);
inline void add_positional(cli_command& cli, const string& name,
    vector<int>& value, const string& usage, const vector<string>& choices = {},
    bool req = true);
inline void add_positional(cli_command& cli, const string& name,
    vector<string>& value, const string& usage,
    const vector<string>& choices = {}, bool req = true);

// Add a subcommand
inline cli_command& add_command(
    cli_command& cli, const string& name, const string& usage);
inline void add_command_name(
    cli_command& cli, const string& name, string& value, const string& usage);

// Parses an optional or positional argument. Optional arguments' names start
// with "--" or "-", otherwise they are arguments. Supports strings, numbers,
// boolean flags and enums.
// Many names, separated by commas, can be used for each argument.
// Boolean flags are indicated with a pair of names "--name/--no-name", so that
// both options are explicitly specified.
inline void add_option(cli_command& cli, const string& name, int& value,
    const string& usage, bool req = false);
inline void add_option(cli_command& cli, const string& name, float& value,
    const string& usage, bool req = false);
inline void add_option(cli_command& cli, const string& name, bool& value,
    const string& usage, bool req = false);
inline void add_option(cli_command& cli, const string& name, string& value,
    const string& usage, bool req = false);
// Parses an optional or positional argument where values can only be within a
// set of choices. Supports strings, integers and enums.
inline void add_option(cli_command& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, bool req = false);
inline void add_option(cli_command& cli, const string& name, string& value,
    const string& usage, const vector<string>& choices, bool req = false);
// Parse all arguments left on the command line. Can only be used as argument.
inline void add_option(cli_command& cli, const string& name,
    vector<string>& value, const string& usage, bool req = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// FORMATTING
// -----------------------------------------------------------------------------
namespace yocto {

// This is a very crude replacement for `std::format()` that will be used when
// available on all platforms.
template <typename... Args>
inline string format(const string& fmt, Args&&... args);

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF COMMAND-LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// Command line value type
enum struct cli_type { integer, uinteger, number, boolean, string };
// Command line value
struct cli_value {
  int64_t  integer  = 0;
  uint64_t uinteger = 0;
  double   number   = 0;
  string   text     = "";
};
// Command line option. All data should be considered private.
struct cli_option {
  string                            name       = "";
  string                            alt        = "";
  bool                              positional = false;
  cli_type                          type       = cli_type::string;
  bool                              req        = false;
  int                               nargs      = 0;
  string                            usage      = "";
  vector<cli_value>                 minmax     = {};
  vector<string>                    choices    = {};
  vector<cli_value>                 value      = {};
  vector<cli_value>                 def        = {};
  bool                              set        = false;
  function<bool(const cli_option&)> set_value  = {};
};
// Command line command. All data should be considered private.
struct cli_command {
  string                        name        = "";
  string                        usage       = "";
  vector<cli_command>           commands    = {};
  vector<cli_option>            options     = {};
  bool                          help        = false;
  string                        command     = "";
  function<void(const string&)> set_command = {};
};

template <typename T, typename>
inline void add_optional(cli_command& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, const string& alt,
    bool req) {
  return add_optional(
      cli, name, (std::underlying_type_t<T>&)value, usage, choices, alt, req);
}
template <typename T, typename>
inline void add_positional(cli_command& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, bool req) {
  return add_positional(
      cli, name, (std::underlying_type_t<T>&)value, usage, choices, req);
}

// Backward compatibility
using cli_state [[deprecated]] = cli_command;

}  // namespace yocto

// -----------------------------------------------------------------------------
// PRINT/FORMATTING UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Print a message to the console
inline void print_info(const string& msg) { printf("%s\n", msg.c_str()); }
// Prints a messgae to the console and exit with an error.
inline int print_fatal(const string& msg) {
  printf("\n%s\n", msg.c_str());
  exit(1);
  return 1;
}

// get time in nanoseconds - useful only to compute difference of times
inline int64_t get_time_() {
  return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}

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
  snprintf(
      buffer, sizeof(buffer), "%02d:%02d:%02d.%03d", hours, mins, secs, msecs);
  return buffer;
}

// Format a large integer number in human readable form
inline string format_num(uint64_t num) {
  auto rem = num % 1000;
  auto div = num / 1000;
  if (div > 0) return format_num(div) + "," + std::to_string(rem);
  return std::to_string(rem);
}

// Print traces for timing and program debugging
inline print_timer print_timed(const string& msg) {
  printf("%s", msg.c_str());
  fflush(stdout);
  // print_info(fmt + " [started]", args...);
  return print_timer{get_time_()};
}
inline int64_t print_elapsed(print_timer& timer) {
  if (timer.start_time < 0) return -1;
  auto elapsed = get_time_() - timer.start_time;
  printf(" in %s\n", format_duration(elapsed).c_str());
  timer.start_time = -1;
  return elapsed;
}
inline print_timer::~print_timer() { print_elapsed(*this); }

// Print progress
inline void print_progress(const string& message, int current, int total) {
  static auto pad = [](const string& str, int n) -> string {
    return string(std::max(0, n - (int)str.size()), '0') + str;
  };
  static auto pade = [](const string& str, int n) -> string {
    return str + string(std::max(0, n - (int)str.size()), ' ');
  };
  static auto pads = [](const string& str, int n) -> string {
    return string(std::max(0, n - (int)str.size()), ' ') + str;
  };
  using clock               = std::chrono::high_resolution_clock;
  static int64_t start_time = 0;
  if (current == 0) start_time = clock::now().time_since_epoch().count();
  auto elapsed = clock::now().time_since_epoch().count() - start_time;
  elapsed /= 1000000;  // millisecs
  auto mins  = pad(std::to_string(elapsed / 60000), 2);
  auto secs  = pad(std::to_string((elapsed % 60000) / 1000), 2);
  auto msecs = pad(std::to_string((elapsed % 60000) % 1000), 3);
  auto cur   = pads(std::to_string(current), 4);
  auto tot   = pads(std::to_string(total), 4);
  auto n     = (int)(20 * (float)current / (float)total);
  auto bar   = "[" + pade(string(n, '='), 20) + "]";
  auto line  = bar + " " + cur + "/" + tot + " " + mins + ":" + secs + "." +
              msecs + " " + pade(message, 30);
  printf("\r%s\r", line.c_str());
  if (current == total) printf("\n");
  fflush(stdout);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIMPLE TIMER
// -----------------------------------------------------------------------------
namespace yocto {

// Simple timer
inline simple_timer::simple_timer() {
  start = get_time_();
  stop  = -1;
}

// Timer opreations
inline void start_timer(simple_timer& timer) {
  timer.start = get_time_();
  timer.stop  = -1;
}
inline void    stop_timer(simple_timer& timer) { timer.stop = get_time_(); }
inline int64_t elapsed_nanoseconds(simple_timer& timer) {
  return get_time_() - timer.start;
}
inline double elapsed_seconds(simple_timer& timer) {
  return (double)(get_time_() - timer.start) * 1e-9;
}
inline string elapsed_formatted(simple_timer& timer) {
  return format_duration(get_time_() - timer.start);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF COMMAND-LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto {

template <typename T>
inline cli_type get_cli_type() {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  if constexpr (std::is_same_v<T, string>) {
    return cli_type::string;
  } else if constexpr (std::is_same_v<T, bool>) {
    return cli_type::boolean;
  } else if constexpr (std::is_enum_v<T>) {
    return cli_type::integer;
  } else if constexpr (std::is_integral_v<T> && !std::is_unsigned_v<T>) {
    return cli_type::integer;
  } else if constexpr (std::is_integral_v<T> && std::is_unsigned_v<T>) {
    return cli_type::uinteger;
  } else if constexpr (std::is_floating_point_v<T>) {
    return cli_type::number;
  } else {
    // probably should be an error
    return cli_type::string;
  }
}

template <typename T>
inline vector<cli_value> make_cli_values(const vector<T>& values) {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  auto cvalues = vector<cli_value>{};
  for (auto value : values) {  // needs to copy for booleans
    auto& cvalue = cvalues.emplace_back();
    if constexpr (std::is_same_v<T, string>) {
      cvalue.text = value;
    } else if constexpr (std::is_same_v<T, bool>) {
      cvalue.integer = value ? 1 : 0;
    } else if constexpr (std::is_enum_v<T>) {
      cvalue.integer = (int64_t)value;
    } else if constexpr (std::is_integral_v<T> && !std::is_unsigned_v<T>) {
      cvalue.integer = value;
    } else if constexpr (std::is_integral_v<T> && std::is_unsigned_v<T>) {
      cvalue.uinteger = value;
    } else if constexpr (std::is_floating_point_v<T>) {
      cvalue.number = value;
    } else {
      // probably should be an error
      // pass
    }
  }
  return cvalues;
}
template <typename T>
inline vector<cli_value> make_cli_values(const T& value) {
  return make_cli_values(vector<T>{value});
}

template <typename T>
inline bool get_value(const cli_option& option, T& value) {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  if (option.value.size() != 1) throw std::out_of_range{"bad option size"};
  auto& cvalue = option.value[0];
  if constexpr (std::is_same_v<T, string>) {
    if (option.type != cli_type::string) return false;
    value = cvalue.text;
  } else if constexpr (std::is_same_v<T, bool>) {
    if (option.type != cli_type::boolean) return false;
    value = cvalue.integer != 0;
  } else if constexpr (std::is_enum_v<T>) {
    if (option.type != cli_type::integer) return false;
    value = (T)cvalue.integer;
  } else if constexpr (std::is_integral_v<T> && !std::is_unsigned_v<T>) {
    if (option.type != cli_type::integer) return false;
    value = (T)cvalue.integer;
  } else if constexpr (std::is_integral_v<T> && std::is_unsigned_v<T>) {
    if (option.type != cli_type::uinteger) return false;
    value = (T)cvalue.uinteger;
  } else if constexpr (std::is_floating_point_v<T>) {
    if (option.type != cli_type::number) return false;
    value = (T)cvalue.number;
  } else {
    return false;
  }
  return true;
}

template <typename T>
inline bool get_value(const cli_option& option, vector<T>& values) {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  values.clear();
  for (auto& cvalue : option.value) {
    auto& value = values.emplace_back();
    if constexpr (std::is_same_v<T, string>) {
      if (option.type != cli_type::string) return false;
      value = cvalue.text;
    } else if constexpr (std::is_same_v<T, bool>) {
      if (option.type != cli_type::boolean) return false;
      value = cvalue.integer != 0;
    } else if constexpr (std::is_enum_v<T>) {
      if (option.type != cli_type::integer) return false;
      value = (T)cvalue.integer;
    } else if constexpr (std::is_integral_v<T> && !std::is_unsigned_v<T>) {
      if (option.type != cli_type::integer) return false;
      value = (T)cvalue.integer;
    } else if constexpr (std::is_integral_v<T> && std::is_unsigned_v<T>) {
      if (option.type != cli_type::uinteger) return false;
      value = (T)cvalue.uinteger;
    } else if constexpr (std::is_floating_point_v<T>) {
      if (option.type != cli_type::number) return false;
      value = (T)cvalue.number;
    } else {
      return false;
    }
  }
  return true;
}

template <typename T>
inline void add_optional_impl(cli_command& cli, const string& name, T& value,
    const string& usage, const vector<T>& minmax, const vector<string>& choices,
    const string& alt, bool req) {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  auto& option      = cli.options.emplace_back();
  option.name       = name;
  option.alt        = alt;
  option.positional = false;
  option.type       = get_cli_type<T>();
  option.req        = req;
  option.nargs      = std::is_same_v<T, bool> ? 0 : 1;
  option.usage      = usage;
  option.minmax     = make_cli_values(minmax);
  option.choices    = choices;
  option.value      = make_cli_values(value);
  option.def        = make_cli_values(value);
  option.set_value  = [&value](const cli_option& option) -> bool {
    return get_value(option, value);
  };
}

template <typename T>
inline void add_positional_impl(cli_command& cli, const string& name, T& value,
    const string& usage, const vector<T>& minmax, const vector<string>& choices,
    bool req) {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  auto& option      = cli.options.emplace_back();
  option.name       = name;
  option.alt        = "";
  option.positional = true;
  option.type       = get_cli_type<T>();
  option.req        = req;
  option.nargs      = 1;
  option.usage      = usage;
  option.minmax     = make_cli_values(value);
  option.choices    = choices;
  option.value      = make_cli_values(value);
  option.def        = make_cli_values(value);
  option.set_value  = [&value](const cli_option& option) -> bool {
    return get_value(option, value);
  };
}

template <typename T>
inline void add_positionalv_impl(cli_command& cli, const string& name,
    vector<T>& value, const string& usage, const vector<T>& minmax,
    const vector<string>& choices, bool req) {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  auto& option      = cli.options.emplace_back();
  option.name       = name;
  option.alt        = "";
  option.positional = true;
  option.type       = get_cli_type<T>();
  option.req        = req;
  option.nargs      = -1;
  option.usage      = usage;
  option.minmax     = make_cli_values(value);
  option.choices    = choices;
  option.value      = make_cli_values(value);
  option.def        = make_cli_values(value);
  option.set_value  = [&value](const cli_option& option) -> bool {
    return get_value(option, value);
  };
}

// Add an optional argument. Supports strings, numbers, and boolean flags.
// Optional arguments will be parsed with name `--<name>` and `-<alt>`.
// Optional booleans will support both `--<name>` and `--no-<name>` to enabled
// and disable the flag.
inline void add_optional(cli_command& cli, const string& name, int& value,
    const string& usage, const vector<int>& minmax, const string& alt,
    bool req) {
  return add_optional_impl(
      cli, name, value, usage, {minmax[0], minmax[1]}, {}, alt, req);
}
inline void add_optional(cli_command& cli, const string& name, float& value,
    const string& usage, const vector<float>& minmax, const string& alt,
    bool req) {
  return add_optional_impl(cli, name, value, usage, minmax, {}, alt, req);
}
inline void add_optional(cli_command& cli, const string& name, bool& value,
    const string& usage, const vector<string>& choices, const string& alt,
    bool req) {
  return add_optional_impl(cli, name, value, usage, {}, choices, alt, req);
}
inline void add_optional(cli_command& cli, const string& name, string& value,
    const string& usage, const vector<string>& choices, const string& alt,
    bool req) {
  return add_optional_impl(cli, name, value, usage, {}, choices, alt, req);
}
inline void add_optional(cli_command& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, const string& alt,
    bool req) {
  return add_optional_impl(cli, name, value, usage, {}, choices, alt, req);
}
// Add a positional argument. Supports strings, numbers, and boolean flags.
inline void add_positional(cli_command& cli, const string& name, int& value,
    const string& usage, const vector<int>& minmax, bool req) {
  return add_positional_impl(cli, name, value, usage, minmax, {}, req);
}
inline void add_positional(cli_command& cli, const string& name, float& value,
    const string& usage, const vector<float>& minmax, bool req) {
  return add_positional_impl(cli, name, value, usage, minmax, {}, req);
}
inline void add_positional(cli_command& cli, const string& name, bool& value,
    const string& usage, const vector<string>& choices, bool req) {
  return add_positional_impl(cli, name, value, usage, {}, choices, req);
}
inline void add_positional(cli_command& cli, const string& name, string& value,
    const string& usage, const vector<string>& choices, bool req) {
  return add_positional_impl(cli, name, value, usage, {}, choices, req);
}
inline void add_positional(cli_command& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, bool req) {
  return add_positional_impl(cli, name, value, usage, {}, choices, req);
}
// Add a positional argument that consumes all arguments left.
// Supports strings and enums.
inline void add_positional(cli_command& cli, const string& name,
    vector<int>& value, const string& usage, const vector<int>& minmax,
    bool req) {
  return add_positionalv_impl(cli, name, value, usage, minmax, {}, req);
}
inline void add_positional(cli_command& cli, const string& name,
    vector<float>& value, const string& usage, const vector<float>& minmax,
    bool req) {
  return add_positionalv_impl(cli, name, value, usage, minmax, {}, req);
}
inline void add_positional(cli_command& cli, const string& name,
    vector<int>& value, const string& usage, const vector<string>& choices,
    bool req) {
  return add_positionalv_impl(cli, name, value, usage, {}, choices, req);
}
inline void add_positional(cli_command& cli, const string& name,
    vector<string>& value, const string& usage, const vector<string>& choices,
    bool req) {
  return add_positionalv_impl(cli, name, value, usage, {}, choices, req);
}

inline vector<string> split_cli_names(const string& name_) {
  auto name  = name_;
  auto split = vector<string>{};
  if (name.empty()) throw std::invalid_argument("option name cannot be empty");
  if (name.find_first_of(" \t\r\n") != string::npos)
    throw std::invalid_argument("option name cannot contain whitespaces");
  while (name.find_first_of(",/") != string::npos) {
    auto pos = name.find_first_of(",/");
    if (pos > 0) split.push_back(name.substr(0, pos));
    name = name.substr(pos + 1);
  }
  if (!name.empty()) split.push_back(name);
  if (split.empty()) throw std::invalid_argument("option name cannot be empty");
  for (auto& name : split)
    if ((split[0][0] == '-') != (name[0] == '-'))
      throw std::invalid_argument("inconsistent option names for " + name);
  return split;
}

template <typename T>
inline void add_option_impl(cli_command& cli, const string& name, T& value,
    const string& usage, bool req) {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  auto& option = cli.options.emplace_back();
  auto  names  = split_cli_names(name);
  if (names[0].find('-') == 0) {
    option.name = names[0].substr(1);
    if (option.name.find('-') == 0) option.name = option.name.substr(1);
    if (names.size() > 1 && names.at(1).find("--no-") != 0) {
      option.alt = names[1].substr(1);
      if (option.alt.find('-') == 0) option.alt = option.alt.substr(1);
    }
  } else {
    option.name       = name;
    option.alt        = "";
    option.positional = true;
  }
  option.type      = get_cli_type<T>();
  option.req       = req;
  option.nargs     = !std::is_same_v<T, bool> ? 1 : 0;
  option.usage     = usage;
  option.value     = make_cli_values(value);
  option.def       = make_cli_values(value);
  option.choices   = {};
  option.set_value = [&value](const cli_option& option) -> bool {
    return get_value(option, value);
  };
}

template <typename T>
inline void add_option_impl(cli_command& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, bool req) {
  static_assert(
      std::is_same_v<T, string> || std::is_integral_v<T> || std::is_enum_v<T>,
      "unsupported type");
  auto& option = cli.options.emplace_back();
  auto  names  = split_cli_names(name);
  if (names[0].find('-') == 0) {
    option.name = names[0].substr(1);
    if (option.name.find('-') == 0) option.name = option.name.substr(1);
    if (names.size() > 1 && names.at(1).find("--no-") != 0) {
      option.alt = names[1].substr(1);
      if (option.alt.find('-') == 0) option.alt = option.alt.substr(1);
    }
  } else {
    option.name       = name;
    option.alt        = "";
    option.positional = true;
  }
  option.type      = get_cli_type<T>();
  option.req       = req;
  option.nargs     = 1;
  option.usage     = usage;
  option.value     = make_cli_values(value);
  option.def       = make_cli_values(value);
  option.choices   = choices;
  option.set_value = [&value](const cli_option& option) -> bool {
    return get_value(option, value);
  };
}

template <typename T>
inline void add_option_impl(cli_command& cli, const string& name,
    vector<T>& value, const string& usage, bool req) {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  auto& option = cli.options.emplace_back();
  auto  names  = split_cli_names(name);
  if (names[0].find('-') == 0) {
    option.name = names[0].substr(1);
    if (option.name.find('-') == 0) option.name = option.name.substr(1);
    if (names.size() > 1 && names.at(1).find("--no-") != 0) {
      option.alt = names[1].substr(1);
      if (option.alt.find('-') == 0) option.alt = option.alt.substr(1);
    }
  } else {
    option.name       = name;
    option.alt        = "";
    option.positional = true;
  }
  option.type      = get_cli_type<T>();
  option.req       = req;
  option.nargs     = -1;
  option.usage     = usage;
  option.value     = make_cli_values(value);
  option.def       = make_cli_values(value);
  option.choices   = {};
  option.set_value = [&value](const cli_option& option) -> bool {
    return get_value(option, value);
  };
}

// Parses an optional or positional argument. Optional arguments' names start
// with "--" or "-", otherwise they are arguments. Supports strings, numbers,
// boolean flags and enums.
// Many names, separated by commas, can be used for each argument.
// Boolean flags are indicated with a pair of names "--name/--no-name", so that
// both options are explicitly specified.
inline void add_option(cli_command& cli, const string& name, int& value,
    const string& usage, bool req) {
  return add_option_impl(cli, name, value, usage, req);
}
inline void add_option(cli_command& cli, const string& name, float& value,
    const string& usage, bool req) {
  return add_option_impl(cli, name, value, usage, req);
}
inline void add_option(cli_command& cli, const string& name, bool& value,
    const string& usage, bool req) {
  return add_option_impl(cli, name, value, usage, req);
}
inline void add_option(cli_command& cli, const string& name, string& value,
    const string& usage, bool req) {
  return add_option_impl(cli, name, value, usage, req);
}
// Parses an optional or positional argument where values can only be within a
// set of choices. Supports strings, integers and enums.
inline void add_option(cli_command& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, bool req) {
  return add_option_impl(cli, name, value, usage, choices, req);
}
inline void add_option(cli_command& cli, const string& name, string& value,
    const string& usage, const vector<string>& choices, bool req) {
  return add_option_impl(cli, name, value, usage, choices, req);
}
// Parse all arguments left on the command line. Can only be used as argument.
inline void add_option(cli_command& cli, const string& name,
    vector<string>& value, const string& usage, bool req) {
  return add_option_impl(cli, name, value, usage, req);
}

// initialize a command line parser
inline cli_command make_cli(const string& name, const string& usage) {
  auto cli  = cli_command{};
  cli.name  = name;
  cli.usage = usage;
  cli.commands.reserve(256);
  return cli;
}

// add command
inline cli_command& add_command(
    cli_command& cli, const string& name, const string& usage) {
  for (auto& cmd : cli.commands) {
    if (cmd.name == name) {
      throw std::invalid_argument{"cannot add two commands with the same name"};
    }
  }
  auto& cmd = cli.commands.emplace_back();
  cmd.name  = name;
  cmd.usage = usage;
  cmd.commands.reserve(256);
  return cmd;
}

inline void add_command_name(
    cli_command& cli, const string& name, string& value, const string& usage) {
  cli.set_command = [&value](const string& cvalue) { value = cvalue; };
}

inline void validate_names(const cli_command& cmd) {
  // check for errors
  auto used = unordered_set<string>{};
  for (auto& option : cmd.options) {
    if (option.name.empty())
      throw std::invalid_argument("name cannot be empty");
    auto names = split_cli_names(option.name);
    if (names.empty()) throw std::invalid_argument("name cannot be empty");
    for (auto& name : names) {
      if (used.find(name) != used.end())
        throw std::invalid_argument("option name " + name + " already in use");
      used.insert(name);
      if ((name[0] == '-') != (option.name[0] == '-'))
        throw std::invalid_argument("inconsistent option type for " + name);
    }
  }
  for (auto& scmd : cmd.commands) validate_names(scmd);
}

inline bool get_help(const cli_command& cli) {
  if (cli.help) return true;
  for (auto& cmd : cli.commands) return get_help(cmd);
  return false;
}

inline string get_usage(const cli_command& root, const cli_command& cli) {
  auto type_name = [](const cli_option& option) -> string {
    auto str = string{};
    str += "<";
    if (option.nargs < 0) str += "[";
    if (!option.choices.empty()) str += "string";
    switch (option.type) {
      case cli_type::integer: str += "integer"; break;
      case cli_type::uinteger: str += "uinteger"; break;
      case cli_type::number: str += "number"; break;
      case cli_type::string: str += "string"; break;
      case cli_type::boolean: str += "boolean"; break;
    }
    if (option.nargs < 0) str += "]";
    str += ">";
    return str;
  };
  auto def_string = [](const cli_option& option) -> string {
    if (option.req) return string{"[required]"};
    auto str = string{};
    str += "[";
    for (auto& value : option.def) {
      switch (option.type) {
        case cli_type::integer:
          str += option.choices.empty() ? std::to_string(value.integer)
                                        : option.choices[value.integer];
          break;
        case cli_type::uinteger:
          str += option.choices.empty() ? std::to_string(value.uinteger)
                                        : option.choices[value.uinteger];
          break;
        case cli_type::number: str += std::to_string(value.number); break;
        case cli_type::string: str += '\"' + value.text + '\"'; break;
        case cli_type::boolean: str += value.integer ? "true" : "false"; break;
      }
    }
    str += "]";
    return str;
  };

  if (!cli.command.empty()) {
    for (auto& subcommand : cli.commands)
      if (cli.command == subcommand.name) return get_usage(root, subcommand);
  }

  auto message      = string{};
  auto has_optional = false, has_positional = false, has_commands = false;
  auto usage_optional = string{}, usage_positional = string{},
       usage_command = string{};
  for (auto& option : cli.options) {
    auto line = string{};
    if (option.positional) {
      line += "  " + option.name;
    } else {
      line += "  --" + option.name;
      if (!option.alt.empty()) line += ", -" + option.alt;
    }
    if (option.nargs > 0) line += " " + type_name(option);
    while (line.size() < 32) line += " ";
    line += option.usage;
    line += " " + def_string(option) + "\n";
    if (!option.choices.empty()) {
      line += "    with choices: ";
      auto len = 16;
      for (auto& choice : option.choices) {
        if (len + choice.size() + 2 > 78) {
          line += "\n                  ";
          len = 16;
        }
        line += choice + ", ";
        len += choice.size() + 2;
      }
      line = line.substr(0, line.size() - 2);
      line += "\n";
    }
    if (option.positional) {
      has_positional = true;
      usage_positional += line;
    } else {
      has_optional = true;
      usage_optional += line;
    }
  }
  for (auto& scmd : cli.commands) {
    has_commands = true;
    auto line    = "  " + scmd.name;
    while (line.size() < 32) line += " ";
    line += scmd.usage + "\n";
    usage_command += line;
  }
  auto is_command = &cli != &root;
  message += "usage: " + root.name + (is_command ? " " + cli.name : "") +
             (has_commands ? " command" : "") +
             (has_optional ? " [options]" : "") +
             (has_positional ? " <arguments>" : "") + "\n";
  message += cli.usage + "\n\n";
  if (has_commands) {
    message += "commands:\n" + usage_command + "\n";
  }
  if (has_optional) {
    message += "options:\n" + usage_optional + "\n";
  }
  if (has_positional) {
    message += "arguments:\n" + usage_positional + "\n";
  }
  return message;
}

inline string get_usage(const cli_command& cli) { return get_usage(cli, cli); }

inline string get_command(const cli_command& cli) { return cli.command; }

inline static bool parse_value(
    cli_option& option, const vector<string>& args, size_t start) {
  option.value.resize(option.nargs > 0 ? option.nargs : (args.size() - start));
  for (auto idx = (size_t)0; idx < option.value.size(); idx++) {
    auto& value   = option.value.at(idx);
    auto& arg     = args.at(start + idx);
    auto& choices = option.choices;
    if (!choices.empty()) {
      if (std::find(choices.begin(), choices.end(), arg) == choices.end())
        return false;
    }
    switch (option.type) {
      case cli_type::string: {
        value.text = arg;
      } break;
      case cli_type::boolean: {
        if (arg == "true" || arg == "1") {
          value.integer = 1;
        } else if (arg == "false" || arg == "0") {
          value.integer = 0;
        } else {
          return false;
        }
      } break;
      case cli_type::integer: {
        if (choices.empty()) {
          auto end      = (char*)nullptr;
          value.integer = (int)strtol(arg.c_str(), &end, 10);
          if (end == nullptr) return false;
        } else {
          value.integer = (int64_t)(
              std::find(choices.begin(), choices.end(), arg) - choices.begin());
        }
      } break;
      case cli_type::uinteger: {
        if (choices.empty()) {
          auto end       = (char*)nullptr;
          value.uinteger = (int)strtoul(arg.c_str(), &end, 10);
          if (end == nullptr) return false;
        } else {
          value.uinteger = (uint64_t)(
              std::find(choices.begin(), choices.end(), arg) - choices.begin());
        }
      } break;
      case cli_type::number: {
        auto end     = (char*)nullptr;
        value.number = strtod(arg.c_str(), &end);
        if (end == nullptr) return false;
      } break;
    }
  }
  return true;
}

inline bool parse_cli(cli_command& cli, vector<string>& args, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  // current parsing state
  auto commands    = vector<cli_command*>{&cli};
  auto positionals = vector<int>{0};

  // parse arguments
  for (auto idx = (size_t)0; idx < args.size(); idx++) {
    auto& cmd = *commands.back();
    auto& arg = args[idx];
    if (arg == "--help") {
      cmd.help = true;
      break;
    }
    auto is_positional = args[idx].find('-') != 0;
    if (!cmd.commands.empty() && is_positional) {
      auto pos = std::find_if(cmd.commands.begin(), cmd.commands.end(),
          [&arg](auto& command) { return command.name == arg; });
      if (pos == cmd.commands.end()) return cli_error("unknown command " + arg);
      cmd.command = arg;
      commands.push_back(&(*pos));
      positionals.push_back(0);
      continue;
    } else if (is_positional) {
      auto pos   = cmd.options.end();
      auto count = 0;
      for (auto it = cmd.options.begin(); it != cmd.options.end(); ++it) {
        auto& option = *it;
        if (!option.positional) continue;
        if (count == positionals.back()) {
          pos = it;
          positionals.back()++;
          break;
        }
        count++;
      }
      if (pos == cmd.options.end())
        return cli_error("too many positional arguments");
      auto& option = *pos;
      option.set   = true;
      if (option.nargs > 0) {
        if (idx + (size_t)option.nargs > args.size())
          return cli_error("missing value for " + option.name);
        if (!parse_value(option, args, idx))
          return cli_error("bad value for " + option.name);
        idx += option.nargs - 1;
      } else if (option.nargs < 0) {
        if (!parse_value(option, args, idx))
          return cli_error("bad value for " + option.name);
        idx += args.size();
      } else {
        throw std::invalid_argument{"unsupported number of arguments"};
      }
    } else {
      auto pos = std::find_if(
          cmd.options.begin(), cmd.options.end(), [&arg](auto& option) {
            if (option.positional) return false;
            return arg == "--" + option.name || arg == "-" + option.alt;
          });
      if (pos == cmd.options.end()) return cli_error("unknown option " + arg);
      auto& option = *pos;
      option.set   = true;
      if (option.nargs == 0) {
        if (option.type != cli_type::boolean)
          throw std::invalid_argument{"unsupported flag type"};
        option.value.resize(1);
        option.value[0].integer = 1;
      } else if (option.nargs > 0) {
        if (idx + (size_t)option.nargs >= args.size())
          return cli_error("missing value for " + option.name);
        if (!parse_value(option, args, idx + 1))
          return cli_error("bad value for " + option.name);
        idx += option.nargs;
      } else {
        throw std::invalid_argument{"unsupported number of arguments"};
      }
    }
  }

  // check for help
  for (auto command_ptr : commands) {
    auto& command = *command_ptr;
    if (command.help) return true;
  }

  // check for required, set defaults and set references
  for (auto command_ptr : commands) {
    auto& command = *command_ptr;
    if (!command.commands.empty() && command.command.empty())
      return cli_error("command not set for " + command.name);
    if (command.set_command) command.set_command(command.command);
    for (auto& option : command.options) {
      if (option.req && !option.set)
        return cli_error("missing value for " + option.name);
      if (!option.set) option.value = option.def;
      if (option.set_value) {
        if (!option.set_value(option)) {
          return cli_error("bad value for " + option.name);
        }
      }
    }
  }

  // done
  return true;
}

inline bool parse_cli(
    cli_command& cli, int argc, const char** argv, string& error) {
  // validate names
  validate_names(cli);
  // prepare args
  auto args = vector<string>{argv + 1, argv + argc};
  // parse
  return parse_cli(cli, args, error);
}

inline void parse_cli(cli_command& cli, int argc, const char** argv) {
  auto error = string{};
  if (!parse_cli(cli, argc, argv, error)) {
    print_info("error: " + error);
    print_info("");
    print_info(get_usage(cli));
    exit(1);
  } else if (get_help(cli)) {
    print_info(get_usage(cli));
    exit(0);
  }
}

}  // namespace yocto

#endif
