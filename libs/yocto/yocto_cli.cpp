//
// Implementation for Yocto/Bvh
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

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_cli.h"

#include <array>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "ext/json.hpp"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::pair;
using std::string;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// PRINT/FORMATTING UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Print a message to the console
void print_info(const string& message) { printf("%s\n", message.c_str()); }
// Prints a message to the console and exit with an error.
int print_fatal(const string& message) {
  printf("\n%s\n", message.c_str());
  exit(1);
  return 1;
}

// get time in nanoseconds - useful only to compute difference of times
int64_t get_time_() {
  return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}

// Format duration string from nanoseconds
string format_duration(int64_t duration) {
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
string format_num(uint64_t num) {
  auto rem = num % 1000;
  auto div = num / 1000;
  if (div > 0) return format_num(div) + "," + std::to_string(rem);
  return std::to_string(rem);
}

// Print traces for timing and program debugging
print_timer print_timed(const string& message) {
  printf("%s", message.c_str());
  fflush(stdout);
  // print_info(fmt + " [started]", args...);
  return print_timer{get_time_()};
}
int64_t print_elapsed(print_timer& timer) {
  if (timer.start_time < 0) return -1;
  auto elapsed = get_time_() - timer.start_time;
  printf(" in %s\n", format_duration(elapsed).c_str());
  timer.start_time = -1;
  return elapsed;
}
print_timer::~print_timer() { print_elapsed(*this); }

// Print progress
void print_progress(const string& message, int current, int total) {
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

static int    print_progress_current = 0;
static int    print_progress_total   = 0;
static string print_progress_message = "";
void          print_progress_begin(const string& message, int total) {
  print_progress_current = 0;
  print_progress_total   = total;
  print_progress_message = message;
  print_progress(
      print_progress_message, print_progress_current, print_progress_total);
}
void print_progress_end() {
  print_progress_current = print_progress_total;
  print_progress(
      print_progress_message, print_progress_current, print_progress_total);
}
void print_progress_next() {
  print_progress_current += 1;
  print_progress(
      print_progress_message, print_progress_current, print_progress_total);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIMPLE TIMER
// -----------------------------------------------------------------------------
namespace yocto {

// Simple timer
simple_timer::simple_timer() {
  start = get_time_();
  stop  = -1;
}

// Timer opreations
void start_timer(simple_timer& timer) {
  timer.start = get_time_();
  timer.stop  = -1;
}
void    stop_timer(simple_timer& timer) { timer.stop = get_time_(); }
int64_t elapsed_nanoseconds(simple_timer& timer) {
  return get_time_() - timer.start;
}
double elapsed_seconds(simple_timer& timer) {
  return (double)(get_time_() - timer.start) * 1e-9;
}
string elapsed_formatted(simple_timer& timer) {
  return format_duration(get_time_() - timer.start);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF COMMAND-LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto {

using cli_json = nlohmann::ordered_json;

template <typename T>
struct cli_is_array {
  static const bool value = false;
};
template <typename T, size_t N>
struct cli_is_array<array<T, N>> {
  static const bool value = true;
};

template <typename T>
struct cli_is_vector {
  static const bool value = false;
};
template <typename T, typename A>
struct cli_is_vector<vector<T, A>> {
  static const bool value = true;
};

template <typename T>
constexpr bool cli_is_integer_v = std::is_same_v<T, int32_t> ||
                                  std::is_same_v<T, int64_t>;
template <typename T>
constexpr bool cli_is_unsigned_v = std::is_same_v<T, uint32_t> ||
                                   std::is_same_v<T, uint64_t>;
template <typename T>
constexpr bool cli_is_floating_v = std::is_same_v<T, float> ||
                                   std::is_same_v<T, double>;
template <typename T>
constexpr bool cli_is_boolean_v = std::is_same_v<T, bool>;
template <typename T>
constexpr bool cli_is_string_v = std::is_same_v<T, string>;
template <typename T>
constexpr bool cli_is_array_v = cli_is_array<T>::value;
template <typename T>
constexpr bool cli_is_vector_v = cli_is_vector<T>::value;

template <typename T>
static void cli_to_value(
    cli_value& cvalue, const T& value, const vector<string>& choices) {
  static_assert(cli_is_integer_v<T> || cli_is_unsigned_v<T> ||
                    cli_is_floating_v<T> || cli_is_boolean_v<T> ||
                    cli_is_string_v<T> || cli_is_array_v<T> ||
                    cli_is_vector_v<T>,
      "unsupported type");
  if constexpr (cli_is_integer_v<T>) {
    if (!choices.empty()) {
      cvalue.type    = cli_type::string;
      cvalue.string_ = choices.at(value);
    } else {
      cvalue.type    = cli_type::integer;
      cvalue.integer = value;
    }
  } else if constexpr (cli_is_unsigned_v<T>) {
    if (!choices.empty()) {
      cvalue.type    = cli_type::string;
      cvalue.string_ = choices.at(value);
    } else {
      cvalue.type      = cli_type::unsigned_;
      cvalue.unsigned_ = value;
    }
  } else if constexpr (cli_is_floating_v<T>) {
    if (!choices.empty()) {
      throw std::invalid_argument{"invalid argument"};
    } else {
      cvalue.type   = cli_type::number;
      cvalue.number = value;
    }
  } else if constexpr (cli_is_boolean_v<T>) {
    if (!choices.empty()) {
      cvalue.type    = cli_type::string;
      cvalue.string_ = choices.at(value ? 1 : 0);
    } else {
    }
  } else if constexpr (cli_is_string_v<T>) {
    if (!choices.empty()) {
      if (std::find(choices.begin(), choices.end(), value) == choices.end())
        throw std::out_of_range{"bad value"};
      cvalue.type    = cli_type::string;
      cvalue.string_ = value;
    } else {
      cvalue.type    = cli_type::string;
      cvalue.string_ = value;
    }
  } else {
    throw std::runtime_error{"type not supported"};
  }
}

template <typename T>
static bool cli_from_value(
    const cli_value& cvalue, T& value, const vector<string>& choices) {
  static_assert(std::is_same_v<T, int32_t> || std::is_same_v<T, int64_t> ||
                    std::is_same_v<T, int32_t> || std::is_same_v<T, int64_t> ||
                    std::is_same_v<T, float> || std::is_same_v<T, double> ||
                    std::is_same_v<T, bool> || std::is_same_v<T, string> ||
                    cli_is_array_v<T> || cli_is_vector_v<T>,
      "unsupported type");
  if constexpr (cli_is_integer_v<T>) {
    if (cvalue.type == cli_type::integer) {
      value = (T)cvalue.integer;
      return true;
    } else if (cvalue.type == cli_type::unsigned_) {
      value = (T)cvalue.unsigned_;
      return true;
    } else if (cvalue.type == cli_type::number) {
      value = (T)cvalue.number;
      return true;
    } else if (cvalue.type == cli_type::string) {
      if (std::find(choices.begin(), choices.end(), cvalue.string_) ==
          choices.end()) {
        return false;
      }
      value = (T)(std::find(choices.begin(), choices.end(), cvalue.string_) -
                  choices.begin());
      return true;
    } else {
      return false;
    }
  } else if constexpr (cli_is_unsigned_v<T>) {
    if (cvalue.type == cli_type::integer) {
      value = (T)cvalue.integer;
      return true;
    } else if (cvalue.type == cli_type::unsigned_) {
      value = (T)cvalue.unsigned_;
      return true;
    } else if (cvalue.type == cli_type::number) {
      value = (T)cvalue.number;
      return true;
    } else if (cvalue.type == cli_type::string) {
      if (std::find(choices.begin(), choices.end(), cvalue.string_) ==
          choices.end()) {
        return false;
      }
      value = (T)(std::find(choices.begin(), choices.end(), cvalue.string_) -
                  choices.begin());
      return true;
    } else {
      return false;
    }
  } else if constexpr (cli_is_floating_v<T>) {
    if (cvalue.type == cli_type::integer) {
      value = (T)cvalue.integer;
      return true;
    } else if (cvalue.type == cli_type::unsigned_) {
      value = (T)cvalue.unsigned_;
      return true;
    } else if (cvalue.type == cli_type::number) {
      value = (T)cvalue.number;
      return true;
    } else {
      return false;
    }
  } else if constexpr (cli_is_boolean_v<T>) {
    if (cvalue.type == cli_type::boolean) {
      value = (T)cvalue.boolean;
      return true;
    } else if (cvalue.type == cli_type::string) {
      if (std::find(choices.begin(), choices.end(), cvalue.string_) ==
          choices.end()) {
        return false;
      }
      value = choices[0] == cvalue.string_ ? false : true;
      return true;
    } else {
      return false;
    }
  } else if constexpr (cli_is_string_v<T>) {
    if (cvalue.type == cli_type::string) {
      if (!choices.empty() && std::find(choices.begin(), choices.end(),
                                  cvalue.string_) == choices.end()) {
        return false;
      }
      value = cvalue.string_;
      return true;
    } else {
      return false;
    }
  } else if constexpr (cli_is_array_v<T>) {
    if (cvalue.type == cli_type::array) {
      if (cvalue.array.size() != value.size()) return false;
      for (auto idx = (size_t)0; idx < cvalue.array.size(); idx++) {
        if (!cli_from_value(cvalue.array[idx], value[idx], choices))
          return false;
      }
      return true;
    } else {
      return false;
    }
  } else if constexpr (cli_is_vector_v<T>) {
    if (cvalue.type == cli_type::array) {
      value.resize(cvalue.array.size());
      for (auto idx = (size_t)0; idx < cvalue.array.size(); idx++) {
        if (!cli_from_value(cvalue.array[idx], value[idx], choices))
          return false;
      }
      return true;
    } else {
      return false;
    }
  } else {
    throw std::runtime_error{"type not supported"};
    return false;
  }
}

template <typename T>
static void cli_to_schema(cli_schema& schema, const T& value,
    const vector<string>& choices, const string& name, const string& usage,
    bool req, bool positional) {
  static_assert(std::is_same_v<T, int32_t> || std::is_same_v<T, int64_t> ||
                    std::is_same_v<T, int32_t> || std::is_same_v<T, int64_t> ||
                    std::is_same_v<T, float> || std::is_same_v<T, double> ||
                    std::is_same_v<T, bool> || std::is_same_v<T, string> ||
                    cli_is_array_v<T> || cli_is_vector_v<T>,
      "unsupported type");
  schema.cli_name       = name;
  schema.description    = usage;
  schema.cli_required   = req;
  schema.cli_positional = positional;
  cli_to_value(schema.default_, value, choices);
  if (!choices.empty()) {
    schema.type  = cli_type::string;
    schema.enum_ = choices;
  } else {
    if constexpr (std::is_same_v<T, int32_t> || std::is_same_v<T, int64_t>) {
      schema.type = cli_type::integer;
    } else if constexpr (std::is_same_v<T, uint32_t> ||
                         std::is_same_v<T, uint64_t>) {
      schema.type = cli_type::unsigned_;
    } else if constexpr (std::is_same_v<T, float> ||
                         std::is_same_v<T, double>) {
      schema.type = cli_type::number;
    } else if constexpr (std::is_same_v<T, bool>) {
      schema.type = cli_type::boolean;
    } else if constexpr (std::is_same_v<T, string>) {
      schema.type = cli_type::string;
    } else if constexpr (cli_is_array_v<T>) {
      schema.type      = cli_type::array;
      schema.min_items = value.size();
      schema.max_items = value.size();
    } else if constexpr (cli_is_vector_v<T>) {
      schema.type      = cli_type::array;
      schema.min_items = 0;
      schema.max_items = std::numeric_limits<size_t>::max();
    }
  }
}

template <typename T>
static void cli_to_schema(cli_json& js, const T& value,
    const vector<string>& choices, const string& name, const string& usage) {
  static_assert(std::is_same_v<T, int32_t> || std::is_same_v<T, int64_t> ||
                    std::is_same_v<T, int32_t> || std::is_same_v<T, int64_t> ||
                    std::is_same_v<T, float> || std::is_same_v<T, double> ||
                    std::is_same_v<T, bool> || std::is_same_v<T, string> ||
                    cli_is_array_v<T> || cli_is_vector_v<T>,
      "unsupported type");
  js["cli_name"]    = name;
  js["description"] = usage;
  cli_to_json(js["default"], value, choices);
  if (!choices.empty()) {
    js["type"] = "string";
    js["enum"] = choices;
  } else {
    if constexpr (std::is_same_v<T, int32_t> || std::is_same_v<T, int64_t>) {
      js["type"] = "integer";
    } else if constexpr (std::is_same_v<T, int32_t> ||
                         std::is_same_v<T, int64_t>) {
      js["type"] = "integer";
    } else if constexpr (std::is_same_v<T, float> ||
                         std::is_same_v<T, double>) {
      js["type"] = "number";
    } else if constexpr (std::is_same_v<T, bool>) {
      js["type"] = "boolean";
    } else if constexpr (std::is_same_v<T, string>) {
      js["type"] = "string";
    } else if constexpr (cli_is_array_v<T>) {
      js["type"]     = "array";
      js["minItems"] = value.size();
      js["maxItems"] = value.size();
    } else if constexpr (cli_is_vector_v<T>) {
      js["type"] = "array";
    }
  }
}

template <typename T>
static bool cli_from_value_(
    const cli_value& cvalue, void* value, const vector<string>& choices) {
  return cli_from_value(cvalue, *(T*)value, choices);
}

static cli_value& get_defaults(const cli_command& cli) {
  if (cli.path.empty())
    return cli.state->defaults;
  else
    return cli.state->defaults.object.at(cli.path);
}
static cli_schema& get_schema(const cli_command& cli) {
  if (cli.path.empty())
    return cli.state->schema;
  else
    return cli.state->schema.properties.at(cli.path);
}

template <typename T>
static void add_option_impl(const cli_command& cli, const string& name,
    T& value, const string& usage, const vector<T>& minmax,
    const vector<string>& choices, const string& alt, bool req) {
  auto& defaults = get_defaults(cli);
  auto& schema   = get_schema(cli);
  cli_to_value(defaults.object[name], value, choices);
  cli_to_schema(
      schema.properties[name], value, choices, name, usage, req, false);
  if (req) schema.required.push_back(name);
  cli.state->variables.push_back(
      cli_variable{cli.path.empty() ? name : (cli.path + "/" + name), &value,
          cli_from_value_<T>, choices});
}

template <typename T, size_t N>
static void add_option_impl(const cli_command& cli, const string& name,
    array<T, N>& value, const string& usage, const vector<T>& minmax,
    const vector<string>& choices, const string& alt, bool req) {
  auto& defaults = get_defaults(cli);
  auto& schema   = get_schema(cli);
  cli_to_value(defaults.object[name], value, choices);
  cli_to_schema(
      schema.properties[name], value, choices, name, usage, req, false);
  if (req) schema.required.push_back(name);
  cli.state->variables.push_back(
      cli_variable{cli.path.empty() ? name : (cli.path + "/" + name), &value,
          cli_from_value_<T>, choices});
}
template <typename T>
static void add_argument_impl(const cli_command& cli, const string& name,
    T& value, const string& usage, const vector<T>& minmax,
    const vector<string>& choices, bool req) {
  auto& defaults = get_defaults(cli);
  auto& schema   = get_schema(cli);
  cli_to_value(defaults.object[name], value, choices);
  cli_to_schema(
      schema.properties[name], value, choices, name, usage, req, true);
  if (req) schema.required.push_back(name);
  schema.cli_positionals.push_back(name);
  cli.state->variables.push_back(
      cli_variable{cli.path.empty() ? name : (cli.path + "/" + name), &value,
          cli_from_value_<T>, choices});
}

template <typename T, size_t N>
static void add_argument_impl(const cli_command& cli, const string& name,
    array<T, N>& value, const string& usage, const vector<T>& minmax,
    const vector<string>& choices, bool req) {
  auto& defaults = get_defaults(cli);
  auto& schema   = get_schema(cli);
  cli_to_value(defaults.object[name], value, choices);
  cli_to_schema(
      schema.properties[name], value, choices, name, usage, req, true);
  if (req) schema.required.push_back(name);
  schema.cli_positionals.push_back(name);
  cli.state->variables.push_back(
      cli_variable{cli.path.empty() ? name : (cli.path + "/" + name), &value,
          cli_from_value_<T>, choices});
}

template <typename T>
static void add_argumentv_impl(const cli_command& cli, const string& name,
    vector<T>& value, const string& usage, const vector<T>& minmax,
    const vector<string>& choices, bool req) {
  auto& defaults = get_defaults(cli);
  auto& schema   = get_schema(cli);
  cli_to_value(defaults.object[name], value, choices);
  cli_to_schema(
      schema.properties[name], value, choices, name, usage, req, true);
  if (req) schema.required.push_back(name);
  schema.cli_positionals.push_back(name);
  cli.state->variables.push_back(
      cli_variable{cli.path.empty() ? name : (cli.path + "/" + name), &value,
          cli_from_value_<T>, choices});
}

// Add an optional argument. Supports strings, numbers, and boolean flags.
// Optional arguments will be parsed with name `--<name>` and `-<alt>`.
// Optional booleans will support both `--<name>` and `--no-<name>` to enabled
// and disable the flag.
void add_option(const cli_command& cli, const string& name, int& value,
    const string& usage, const vector<int>& minmax, const string& alt,
    bool req) {
  add_option_impl(cli, name, value, usage, minmax, {}, alt, req);
}
void add_option(const cli_command& cli, const string& name, float& value,
    const string& usage, const vector<float>& minmax, const string& alt,
    bool req) {
  add_option_impl(cli, name, value, usage, minmax, {}, alt, req);
}
void add_option(const cli_command& cli, const string& name, bool& value,
    const string& usage, const vector<string>& choices, const string& alt,
    bool req) {
  add_option_impl(cli, name, value, usage, {}, choices, alt, req);
}
void add_option(const cli_command& cli, const string& name, string& value,
    const string& usage, const vector<string>& choices, const string& alt,
    bool req) {
  add_option_impl(cli, name, value, usage, {}, choices, alt, req);
}
void add_option_with_config(const cli_command& cli, const string& name,
    string& value, const string& usage, const string& config, const string& alt,
    bool req) {
  add_option_impl(cli, name, value, usage, {}, {}, alt, req);
  if (!config.empty()) {
    get_schema(cli).properties.at(name).cli_config = config;
  }
}
void add_option(const cli_command& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, const string& alt,
    bool req) {
  add_option_impl(cli, name, value, usage, {}, choices, alt, req);
}
// Add a positional argument. Supports strings, numbers, and boolean flags.
void add_argument(const cli_command& cli, const string& name, int& value,
    const string& usage, const vector<int>& minmax, bool req) {
  add_argument_impl(cli, name, value, usage, minmax, {}, req);
}
void add_argument(const cli_command& cli, const string& name, float& value,
    const string& usage, const vector<float>& minmax, bool req) {
  add_argument_impl(cli, name, value, usage, minmax, {}, req);
}
void add_argument(const cli_command& cli, const string& name, bool& value,
    const string& usage, const vector<string>& choices, bool req) {
  add_argument_impl(cli, name, value, usage, {}, choices, req);
}
void add_argument(const cli_command& cli, const string& name, string& value,
    const string& usage, const vector<string>& choices, bool req) {
  add_argument_impl(cli, name, value, usage, {}, choices, req);
}
void add_argument_with_config(const cli_command& cli, const string& name,
    string& value, const string& usage, const string& config, bool req) {
  add_argument_impl(cli, name, value, usage, {}, {}, req);
  if (!config.empty()) {
    get_schema(cli).properties.at(name).cli_config = config;
  }
}
void add_argument(const cli_command& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, bool req) {
  add_argument_impl(cli, name, value, usage, {}, choices, req);
}
// Add a positional argument that consumes all arguments left.
// Supports strings and enums.
void add_argument(const cli_command& cli, const string& name,
    vector<int>& value, const string& usage, const vector<int>& minmax,
    bool req) {
  add_argumentv_impl(cli, name, value, usage, minmax, {}, req);
}
void add_argument(const cli_command& cli, const string& name,
    vector<float>& value, const string& usage, const vector<float>& minmax,
    bool req) {
  add_argumentv_impl(cli, name, value, usage, minmax, {}, req);
}
void add_argument(const cli_command& cli, const string& name,
    vector<int>& value, const string& usage, const vector<string>& choices,
    bool req) {
  add_argumentv_impl(cli, name, value, usage, {}, choices, req);
}
void add_argument(const cli_command& cli, const string& name,
    vector<string>& value, const string& usage, const vector<string>& choices,
    bool req) {
  add_argumentv_impl(cli, name, value, usage, {}, choices, req);
}

// Add an optional argument. Supports basic math types.
void add_option(const cli_command& cli, const string& name, vec2i& value,
    const string& usage, const vector<int>& minmax, const string& alt,
    bool req) {
  add_option_impl(
      cli, name, (array<int, 2>&)value, usage, minmax, {}, alt, req);
}
void add_option(const cli_command& cli, const string& name, vec3i& value,
    const string& usage, const vector<int>& minmax, const string& alt,
    bool req) {
  add_option_impl(
      cli, name, (array<int, 3>&)value, usage, minmax, {}, alt, req);
}
void add_option(const cli_command& cli, const string& name, vec4i& value,
    const string& usage, const vector<int>& minmax, const string& alt,
    bool req) {
  add_option_impl(
      cli, name, (array<int, 4>&)value, usage, minmax, {}, alt, req);
}
void add_option(const cli_command& cli, const string& name, vec2f& value,
    const string& usage, const vector<float>& minmax, const string& alt,
    bool req) {
  add_option_impl(
      cli, name, (array<float, 2>&)value, usage, minmax, {}, alt, req);
}
void add_option(const cli_command& cli, const string& name, vec3f& value,
    const string& usage, const vector<float>& minmax, const string& alt,
    bool req) {
  add_option_impl(
      cli, name, (array<float, 3>&)value, usage, minmax, {}, alt, req);
}
void add_option(const cli_command& cli, const string& name, vec4f& value,
    const string& usage, const vector<float>& minmax, const string& alt,
    bool req) {
  add_option_impl(
      cli, name, (array<float, 4>&)value, usage, minmax, {}, alt, req);
}

static string schema_to_usage(
    const cli_schema& schema_, const string& command, const string& program) {
  auto type_to_string = [](const cli_schema& schema) -> string {
    auto type = schema.type;
    switch (type) {
      case cli_type::none: return "";
      case cli_type::integer: return "<integer>";
      case cli_type::unsigned_: return "<integer>";
      case cli_type::number: return "<number>";
      case cli_type::boolean: return "<bool>";
      case cli_type::string: return "<string>";
      case cli_type::array: return "<array>";
      case cli_type::object: return "<object>";
    }
  };
  auto default_to_string = [](const cli_schema& schema) -> string {
    auto& value = schema.default_;
    switch (value.type) {
      case cli_type::none: return "";
      case cli_type::integer: return "[" + std::to_string(value.integer) + "]";
      case cli_type::unsigned_:
        return "[" + std::to_string(value.unsigned_) + "]";
      case cli_type::number: return "[" + std::to_string(value.number) + "]";
      case cli_type::boolean:
        return "[" + (value.boolean ? string{"true"} : string{"false"}) + "]";
      case cli_type::string: return "[" + value.string_ + "]";
      case cli_type::array: return "[]";
      case cli_type::object: return "";
    }
  };
  auto& schema   = command.empty() ? schema_ : schema_.properties.at(command);
  auto  progname = program;
  if (progname.rfind('/') != string::npos) {
    progname = progname.substr(progname.rfind('/') + 1);
  }
  if (progname.rfind('\\') != string::npos) {
    progname = progname.substr(progname.rfind('\\') + 1);
  }
  auto message       = string{};
  auto usage_options = string{}, usage_arguments = string{},
       usage_commands = string{};
  for (auto& [key, property] : schema.properties) {
    if (property.type == cli_type::object) {
      auto line = "  " + key;
      while (line.size() < 32) line += " ";
      line += property.description + "\n";
      usage_commands += line;
    } else {
      auto is_positional = std::find(schema.cli_positionals.begin(),
                               schema.cli_positionals.end(),
                               key) != schema.cli_positionals.end();
      auto line          = (is_positional ? "  " : "  --") + key;
      if (property.type != cli_type::boolean || is_positional) {
        line += " " + type_to_string(property);
      }
      while (line.size() < 32) line += " ";
      line += property.description;
      if (property.default_.type != cli_type::none) {
        line += " " + default_to_string(property);
      }
      line += "\n";
      if (!property.enum_.empty()) {
        line += "    with choices: ";
        auto len = 16;
        for (auto& choice : property.enum_) {
          if (len + choice.size() + 2 > 78) {
            line += "\n                  ";
            len = 16;
          }
          line += (string)choice + ", ";
          len += (int)choice.size() + 2;
        }
        line = line.substr(0, line.size() - 2);
        line += "\n";
      }
      if (is_positional)
        usage_arguments += line;
      else
        usage_options += line;
    }
  }
  usage_options += "  --help                        Prints help. [false]\n";
  usage_options +=
      "  --config <string>             Load configuration. ["
      "]\n";
  message += "usage: " + progname + (command.empty() ? "" : (" " + command)) +
             (!usage_commands.empty() ? " command" : "") +
             (!usage_options.empty() ? " [options]" : "") +
             (!usage_arguments.empty() ? " <arguments>" : "") + "\n";
  message += schema.description + "\n\n";
  if (!usage_commands.empty()) {
    message += "commands:\n" + usage_commands + "\n";
  }
  if (!usage_options.empty()) {
    message += "options:\n" + usage_options + "\n";
  }
  if (!usage_arguments.empty()) {
    message += "arguments:\n" + usage_arguments + "\n";
  }
  return message;
}

static bool arg_to_value(cli_value& value, const cli_schema& schema,
    const string& name, bool positional, const vector<string>& args,
    size_t& idx, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  if (!schema.enum_.empty()) {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    value.type    = cli_type::string;
    value.string_ = args[idx++];
  } else if (schema.type == cli_type::integer) {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    auto end      = (char*)nullptr;
    value.type    = cli_type::integer;
    value.integer = strtol(args[idx++].c_str(), &end, 10);
    if (end == nullptr) return cli_error("bad value for " + name);
  } else if (schema.type == cli_type::unsigned_) {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    auto end      = (char*)nullptr;
    value.type    = cli_type::unsigned_;
    value.integer = strtoul(args[idx++].c_str(), &end, 10);
    if (end == nullptr) return cli_error("bad value for " + name);
  } else if (schema.type == cli_type::number) {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    auto end     = (char*)nullptr;
    value.type   = cli_type::number;
    value.number = strtod(args[idx++].c_str(), &end);
    if (end == nullptr) return cli_error("bad value for " + name);
  } else if (schema.type == cli_type::boolean) {
    if (positional) {
      if (idx >= args.size()) return cli_error("missing value for " + name);
      value.type    = cli_type::string;
      value.string_ = args[idx++] == "true" ? true : false;
    } else {
      value.type    = cli_type::boolean;
      value.boolean = true;
    }
  } else if (schema.type == cli_type::string) {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    value.type    = cli_type::string;
    value.string_ = args[idx++];
  } else if (schema.type == cli_type::array) {
    value.type = cli_type::array;
    if (idx + schema.min_items >= args.size())
      return cli_error("missing value for " + name);
    auto end = std::min(idx + schema.max_items, args.size());
    while (idx < end) {
      if (!arg_to_value(value.array.emplace_back(), schema.items.at(0), name,
              positional, args, idx, error))
        return false;
    }
  } else {
    throw std::runtime_error("unsupported type");
  }
  return true;
}

static bool args_to_value(cli_value& value, const cli_schema& schema,
    const vector<string>& args, size_t idx, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  auto get_try_config = [](const string& base_, const string& name) {
    auto base = base_;
    if (base.rfind('/') != base.npos) {
      base = base.substr(0, base.rfind('/'));
    } else if (base.rfind('\\') != base.npos) {
      base = base.substr(0, base.rfind('\\'));
    } else {
      base = ".";
    }
    return base + "/" + name;
  };

  auto contains = [](const auto& object, const string& name) {
    return object.find(name) != object.end();
  };

  // init
  value.type   = cli_type::object;
  value.object = {};

  // add things to schema
  auto commands = vector<string>{}, positionals = vector<string>{};
  for (auto& [key, value] : schema.properties) {
    if (value.type == cli_type::object) commands.push_back(key);
  }
  if (!schema.cli_positionals.empty()) {
    for (auto& key : schema.cli_positionals) {
      positionals.push_back((string)key);
    }
  }

  // current parsing state
  auto positional = 0;

  // parse arguments
  while (idx < args.size()) {
    auto& arg = args[idx++];
    if (arg == "--help") {
      value.object["help"].type    = cli_type::boolean;
      value.object["help"].boolean = true;
      continue;
    }
    if (arg == "--config") {
      if (idx >= args.size()) return cli_error("missing value for config");
      value.object["config"].type    = cli_type::string;
      value.object["config"].string_ = args[idx++];
      continue;
    }
    auto is_positional = arg.find('-') != 0;
    if (!commands.empty() && is_positional) {
      if (std::find(commands.begin(), commands.end(), arg) != commands.end()) {
        value.object["command"].type    = cli_type::string;
        value.object["command"].string_ = arg;
        if (!args_to_value(
                value.object[arg], schema.properties.at(arg), args, idx, error))
          return false;
        break;
      } else {
        return cli_error("unknown command " + arg);
      }
    } else if (is_positional) {
      if (positional >= positionals.size())
        return cli_error("too many positional arguments");
      auto  name    = positionals[positional++];
      auto& oschema = schema.properties.at(name);
      idx--;
      if (!arg_to_value(value.object[name], oschema, arg, is_positional, args,
              idx, error))
        return false;
      if (!oschema.cli_config.empty() && !contains(value.object, "config")) {
        value.object["config"].type = cli_type::string;
        value.object["config"].string_ =
            "try:" +
            get_try_config(value.object.at(name).string_, oschema.cli_config);
      }
    } else {
      auto name = string{};
      for (auto& [key, value] : schema.properties) {
        if ("--" + key == arg && value.type != cli_type::object) {
          name = key;
          break;
        }
      }
      if (name == "") return cli_error("unknown option " + arg);
      auto& oschema = schema.properties.at(name);
      if (!arg_to_value(value.object[name], oschema, name, is_positional, args,
              idx, error))
        return false;
      if (!oschema.cli_config.empty() && !contains(value.object, "config")) {
        value.object["config"].type = cli_type::string;
        value.object["config"].string_ =
            "try:" +
            get_try_config(value.object.at(name).string_, oschema.cli_config);
      }
    }
  }

  // done
  return true;
}

static bool validate_value(const cli_value& value, const cli_schema& schema,
    const string& name, bool check_required, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  auto contains = [](const auto& object, const string& name) {
    return object.find(name) != object.end();
  };

  switch (value.type) {
    case cli_type::none: {
      return cli_error("bad value for " + name);
    } break;
    case cli_type::integer: {
      if (schema.type != cli_type::integer &&
          schema.type != cli_type::unsigned_ && schema.type != cli_type::number)
        return cli_error("bad value for " + name);
      return true;
    } break;
    case cli_type::unsigned_: {
      if (schema.type != cli_type::integer &&
          schema.type != cli_type::unsigned_ && schema.type != cli_type::number)
        return cli_error("bad value for " + name);
      return true;
    } break;
    case cli_type::number: {
      if (schema.type != cli_type::number)
        return cli_error("bad value for " + name);
      return true;
    } break;
    case cli_type::boolean: {
      if (schema.type != cli_type::boolean)
        return cli_error("bad value for " + name);
      return true;
    } break;
    case cli_type::string: {
      if (schema.type != cli_type::string)
        return cli_error("bad value for " + name);
      return true;
    } break;
    case cli_type::array: {
      if (schema.type != cli_type::array)
        return cli_error("bad value for " + name);
      if (schema.min_items > value.array.size())
        return cli_error("bad value for " + name);
      if (schema.max_items < value.array.size())
        return cli_error("bad value for " + name);
      for (auto& item : value.array)
        if (!validate_value(item, schema.items.at(0), name, false, error))
          return false;
      return true;
    } break;
    case cli_type::object: {
      if (schema.type != cli_type::object)
        return cli_error("bad value for " + name);
      for (auto& [key, property] : value.object) {
        if (key == "help") {
          if (property.type != cli_type::boolean)
            return cli_error("bad value for " + key);
        } else if (key == "config") {
          if (property.type != cli_type::string)
            return cli_error("bad value for " + key);
        } else if (key == "command") {
          if (property.type != cli_type::string)
            return cli_error("bad value for " + key);
        } else {
          if (!contains(schema.properties, key))
            return cli_error("unknown option " + key);
          auto selected_command = contains(value.object, "command") &&
                                  value.object.at("command").type ==
                                      cli_type::string &&
                                  value.object.at("command").string_ == key;
          if (!validate_value(property, schema.properties.at(key), key,
                  check_required && selected_command, error))
            return false;
        }
      }
      if (check_required) {
        for (auto& req : schema.required) {
          if (!contains(value.object, req))
            return cli_error("missing value for " + (string)req);
        }
      }
      return true;
    } break;
    default: {
      return cli_error("bad value for " + name);
    } break;
  }

  return true;
}

// update json objects
void update_value_objects(cli_value& value, const cli_value& update) {
  auto contains = [](const auto& object, const string& name) {
    return object.find(name) != object.end();
  };

  if (value.type != cli_type::object) return;
  for (auto& [key, property] : update.object) {
    if (property.type == cli_type::object && contains(value.object, key) &&
        value.object.at(key).type == cli_type::object) {
      update_value_objects(value.object.at(key), value);
    } else {
      value.object[key] = property;
    }
  }
}

// set variables
static bool value_to_variables(const cli_value& value,
    const vector<cli_variable>& variables, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  auto contains = [](const auto& object, const string& name) {
    return object.find(name) != object.end();
  };

  for (auto& variable : variables) {
    if (variable.path.find('/') != string::npos) {
      auto name = variable.path;
      if (!contains(value.object, name)) continue;
      if (!variable.setter(
              value.object.at(name), variable.value, variable.choices))
        return cli_error("bad value for " + variable.path);
    } else {
      auto name    = variable.path.substr(variable.path.find('/') + 1);
      auto command = variable.path.substr(0, variable.path.find('/'));
      if (!contains(value.object, command)) continue;
      if (value.object.at(command).type != cli_type::object) continue;
      if (!contains(value.object.at(command).object, name)) continue;
      if (!variable.setter(value.object.at(command).object.at(name),
              variable.value, variable.choices))
        return cli_error("bad value for " + variable.path);
    }
  }
  return true;
}

// update json objects
void update_json_objects(cli_json& js, const cli_json& update) {
  if (!js.is_object()) return;
  for (auto& [key, value] : update.items()) {
    if (value.is_object() && js.contains(key) && js.at(key).is_object()) {
      update_json_objects(js.at(key), value);
    } else {
      js[key] = value;
    }
  }
}

void from_json(const cli_json& js, cli_value& value) {
  switch (js.type()) {
    case cli_json::value_t::null: {
      value.type = cli_type::none;
    } break;
    case cli_json::value_t::number_integer: {
      value.type    = cli_type::integer;
      value.integer = (int64_t)js;
    } break;
    case cli_json::value_t::number_unsigned: {
      value.type      = cli_type::unsigned_;
      value.unsigned_ = (uint64_t)js;
    } break;
    case cli_json::value_t::number_float: {
      value.type   = cli_type::number;
      value.number = (double)js;
    } break;
    case cli_json::value_t::boolean: {
      value.type    = cli_type::boolean;
      value.boolean = (bool)js;
    } break;
    case cli_json::value_t::string: {
      value.type    = cli_type::string;
      value.string_ = (string)js;
    } break;
    case cli_json::value_t::array: {
      value.type = cli_type::array;
      for (auto& jitem : js) from_json(jitem, value.array.emplace_back());
    } break;
    case cli_json::value_t::object: {
      value.type = cli_type::object;
      for (auto& [key, jitem] : js.items()) from_json(jitem, value.object[key]);
    } break;
    case cli_json::value_t::binary: {
      value.type = cli_type::none;
    } break;
    case cli_json::value_t::discarded: {
      value.type = cli_type::none;
    } break;
  }
}

// grabs a configuration and update json arguments
static bool config_to_value(cli_value& value, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  auto get_config = [](const cli_value& value) -> string {
    auto contains = [](const auto& object, const string& name) {
      return object.find(name) != object.end();
    };

    auto current = &value;
    while (true) {
      if (current->type != cli_type::object) break;
      if (contains(current->object, "config") &&
          current->object.at("config").type == cli_type::string)
        return current->object.at("config").string_;
      if (contains(current->object, "command") &&
          current->object.at("command").type == cli_type::string) {
        current = &current->object.at(current->object.at("command").string_);
      } else {
        break;
      }
    }
    return "";
  };

  auto try_config = false;
  auto config     = get_config(value);
  if (config.empty()) return true;
  if (config.find("try:") == 0) {
    config     = config.substr(4);
    try_config = true;
  }

  auto js = cli_json{};
  auto fs = std::ifstream(config);
  if (!fs) {
    if (try_config) return true;
    return cli_error("missing configuration file " + config);
  }
  try {
    js = cli_json::parse(fs);
  } catch (...) {
    return cli_error("error loading configuration " + config);
  }

  auto cvalue = cli_value{};
  try {
    from_json(js, cvalue);
  } catch (...) {
    return cli_error("error converting configuration " + config);
  }

  update_value_objects(cvalue, value);
  value = cvalue;
  return true;
}

// initialize a command line parser
cli_state make_cli(const string& name, const string& usage) {
  auto cli               = cli_state{};
  cli.defaults.type      = cli_type::object;
  cli.defaults.object    = {};
  cli.schema.cli_name    = name;
  cli.schema.description = usage;
  cli.schema.type        = cli_type::object;
  cli.schema.properties  = {};
  return cli;
}

// add command
cli_command add_command(
    const cli_command& cli, const string& name, const string& usage) {
  auto& defaults                      = get_defaults(cli);
  defaults.object[name].type          = cli_type::object;
  auto& schema                        = get_schema(cli);
  schema.properties[name].cli_name    = name;
  schema.properties[name].description = usage;
  schema.properties[name].type        = cli_type::object;
  return {cli.state, cli.path.empty() ? name : (cli.path + "/" + name)};
}

void set_command_var(const cli_command& cli, string& value) {
  cli.state->variables.push_back(
      cli_variable{cli.path.empty() ? "command" : (cli.path + "/command"),
          &value, cli_from_value_<string>, {}});
}

void set_help_var(const cli_command& cli, bool& value) {
  cli.state->variables.push_back(
      cli_variable{cli.path.empty() ? "help" : (cli.path + "/help"), &value,
          cli_from_value_<string>, {}});
}

void add_command_name(const cli_command& cli, const string& name, string& value,
    const string& usage) {
  cli.state->variables.push_back(
      cli_variable{cli.path.empty() ? "command" : (cli.path + "/command"),
          &value, cli_from_value_<string>, {}});
}

static bool parse_cli(cli_state& cli, const vector<string>& args,
    cli_value& values, string& error) {
  auto idx = (size_t)1;
  if (!args_to_value(values, get_schema(cli), args, idx, error)) return false;
  if (!config_to_value(values, error)) return false;
  if (!validate_value(values, get_schema(cli), "", true, error)) return false;
  if (!value_to_variables(values, cli.variables, error)) return false;
  return true;
}

bool parse_cli(cli_state& cli, const vector<string>& args, string& error) {
  auto values = cli_value{};
  return parse_cli(cli, args, values, error);
}

void parse_cli(cli_state& cli, const vector<string>& args) {
  auto get_command = [](const cli_value& value) -> string {
    auto contains = [](const auto& object, const string& name) {
      return object.find(name) != object.end();
    };

    auto command = string{};
    auto current = &value;
    while (true) {
      if (!contains(current->object, "command")) break;
      command += (command.empty() ? "" : "/") +
                 current->object.at("command").string_;
      current = &current->object.at(current->object.at("command").string_);
    }
    return command;
  };
  auto get_help = [&](const cli_value& value) -> bool {
    auto contains = [](const auto& object, const string& name) {
      return object.find(name) != object.end();
    };

    auto current = &value;
    while (true) {
      if (contains(current->object, "help"))
        return current->object.at("help").boolean;
      if (contains(current->object, "command")) {
        current = &current->object.at(current->object.at("command").string_);
      } else {
        break;
      }
    }
    return false;
  };
  auto values = cli_value{};
  auto error  = string{};
  if (!parse_cli(cli, args, values, error)) {
    print_info("error: " + error);
    print_info("");
    print_info(schema_to_usage(get_schema(cli), get_command(values), args[0]));
    exit(1);
  } else if (get_help(values)) {
    print_info(schema_to_usage(get_schema(cli), get_command(values), args[0]));
    exit(0);
  }
}

void parse_cli(cli_state& cli, int argc, const char** argv) {
  parse_cli(cli, vector<string>{argv, argv + argc});
}

}  // namespace yocto
