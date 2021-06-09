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
static void cli_to_json(
    cli_json& js, const T& value, const vector<string>& choices) {
  if (!choices.empty()) {
    if constexpr (std::is_same_v<T, int32_t> || std::is_same_v<T, int64_t> ||
                  std::is_same_v<T, uint32_t> || std::is_same_v<T, uint64_t>) {
      js = choices.at((int)value);
    } else if constexpr (std::is_same_v<T, string>) {
      js = value;
    } else if constexpr (std::is_same_v<T, bool>) {
      js = choices.at(value ? 1 : 0);
    } else {
      throw std::runtime_error{"type not supported"};
    }
  } else {
    js = value;
  }
}

template <typename T>
struct cli_is_array {
  static const bool value = false;
};
template <typename T, size_t N>
struct cli_is_array<array<T, N>> {
  static const bool value = true;
};
template <typename T>
constexpr bool cli_is_array_v = cli_is_array<T>::value;

template <typename T>
struct cli_is_vector {
  static const bool value = false;
};
template <typename T, typename A>
struct cli_is_vector<vector<T, A>> {
  static const bool value = true;
};
template <typename T>
constexpr bool cli_is_vector_v = cli_is_vector<T>::value;

template <typename T>
static bool cli_from_json(
    const cli_json& js, T& value, const vector<string>& choices) {
  static_assert(std::is_same_v<T, int32_t> || std::is_same_v<T, int64_t> ||
                    std::is_same_v<T, int32_t> || std::is_same_v<T, int64_t> ||
                    std::is_same_v<T, float> || std::is_same_v<T, double> ||
                    std::is_same_v<T, bool> || std::is_same_v<T, string> ||
                    cli_is_array_v<T> || cli_is_vector_v<T>,
      "unsupported type");
  if (!choices.empty()) {
    auto values = (string)js;
    if (std::find(choices.begin(), choices.end(), values) == choices.end()) {
      return false;
    }
    if constexpr (std::is_same_v<T, int32_t> || std::is_same_v<T, int64_t> ||
                  std::is_same_v<T, uint32_t> || std::is_same_v<T, uint64_t>) {
      value = (T)(std::find(choices.begin(), choices.end(), values) -
                  choices.begin());
    } else if constexpr (std::is_same_v<T, string>) {
      value = values;
    } else if constexpr (std::is_same_v<T, bool>) {
      value = values == choices[0] ? false : true;
    } else {
      throw std::runtime_error{"type not supported"};
    }
  } else {
    value = (T)js;
  }
  return true;
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
static bool cli_from_json_(
    const void* js, void* value, const vector<string>& choices) {
  return cli_from_json(*(const cli_json*)js, *(T*)value, choices);
}

static cli_json& get_defaults(const cli_command& cli) {
  auto& root = *(cli_json*)cli.state->defaults.get();
  if (cli.path.empty())
    return root;
  else
    return root.at(cli.path);
}
static cli_json& get_schema(const cli_command& cli) {
  auto& root = *(cli_json*)cli.state->schema.get();
  if (cli.path.empty())
    return root;
  else
    return root.at("properties").at(cli.path);
}

template <typename T>
static void add_option_impl(const cli_command& cli, const string& name,
    T& value, const string& usage, const vector<T>& minmax,
    const vector<string>& choices, const string& alt, bool req) {
  auto& defaults = get_defaults(cli);
  auto& schema   = get_schema(cli);
  cli_to_json(defaults[name], value, choices);
  cli_to_schema(schema["properties"][name], value, choices, name, usage);
  if (req) schema["required"].push_back(name);
  cli.state->variables.push_back(
      cli_variable{cli.path.empty() ? name : (cli.path + "/" + name), &value,
          cli_from_json_<T>, choices});
}

template <typename T, size_t N>
static void add_option_impl(const cli_command& cli, const string& name,
    array<T, N>& value, const string& usage, const vector<T>& minmax,
    const vector<string>& choices, const string& alt, bool req) {
  auto& defaults = get_defaults(cli);
  auto& schema   = get_schema(cli);
  cli_to_json(defaults[name], value, choices);
  cli_to_schema(schema["properties"][name], value, choices, name, usage);
  if (req) schema["required"].push_back(name);
  cli.state->variables.push_back(
      cli_variable{cli.path.empty() ? name : (cli.path + "/" + name), &value,
          cli_from_json_<T>, choices});
}
template <typename T>
static void add_argument_impl(const cli_command& cli, const string& name,
    T& value, const string& usage, const vector<T>& minmax,
    const vector<string>& choices, bool req) {
  auto& defaults = get_defaults(cli);
  auto& schema   = get_schema(cli);
  cli_to_json(defaults[name], value, choices);
  cli_to_schema(schema["properties"][name], value, choices, name, usage);
  if (req) schema["required"].push_back(name);
  if (!schema.contains("cli_positionals"))
    schema["cli_positionals"] = cli_json::array();
  schema["cli_positionals"].push_back(name);
  cli.state->variables.push_back(
      cli_variable{cli.path.empty() ? name : (cli.path + "/" + name), &value,
          cli_from_json_<T>, choices});
}

template <typename T, size_t N>
static void add_argument_impl(const cli_command& cli, const string& name,
    array<T, N>& value, const string& usage, const vector<T>& minmax,
    const vector<string>& choices, bool req) {
  auto& defaults = get_defaults(cli);
  auto& schema   = get_schema(cli);
  cli_to_json(defaults[name], value, choices);
  cli_to_schema(schema["properties"][name], value, choices, name, usage);
  if (req) schema["required"].push_back(name);
  if (!schema.contains("cli_positionals"))
    schema["cli_positionals"] = cli_json::array();
  schema["cli_positionals"].push_back(name);
  cli.state->variables.push_back(
      cli_variable{cli.path.empty() ? name : (cli.path + "/" + name), &value,
          cli_from_json_<T>, choices});
}

template <typename T>
static void add_argumentv_impl(const cli_command& cli, const string& name,
    vector<T>& value, const string& usage, const vector<T>& minmax,
    const vector<string>& choices, bool req) {
  auto& defaults = get_defaults(cli);
  auto& schema   = get_schema(cli);
  cli_to_json(defaults[name], value, choices);
  cli_to_schema(schema["properties"][name], value, choices, name, usage);
  if (req) schema["required"].push_back(name);
  if (!schema.contains("cli_positionals"))
    schema["positionals"] = cli_json::array();
  schema["positionals"].push_back(name);
  cli.state->variables.push_back(
      cli_variable{cli.path.empty() ? name : (cli.path + "/" + name), &value,
          cli_from_json_<T>, choices});
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
    get_schema(cli).at("properties").at(name)["cli_config"] = config;
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
    get_schema(cli).at("properties").at(name)["cli_config"] = config;
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
    const cli_json& schema_, const string& command, const string& program) {
  auto& schema   = command.empty() ? schema_
                                   : schema_.at("properties").at(command);
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
  for (auto& [key, value] : schema.at("properties").items()) {
    if (value.value("type", "") == "object") {
      auto line = "  " + key;
      while (line.size() < 32) line += " ";
      line += value.value("description", "") + "\n";
      usage_commands += line;
    } else {
      auto is_positional = false;
      auto line          = (is_positional ? "  " : "  --") + key;
      if (value.value("type", "") != "boolean" || is_positional) {
        line += " <" + value.value("type", "") + ">";
      }
      while (line.size() < 32) line += " ";
      line += value.value("description", "");
      if (value.contains("default")) {
        auto def = value.at("default").dump();
        if (value.value("type", "") == "string")
          def = def.substr(1, def.length() - 2);
        line += " [" + def + "]\n";
        // TODO: required
      } else {
        line += "\n";
      }
      if (value.contains("enum")) {
        line += "    with choices: ";
        auto len = 16;
        for (auto& choice : value.at("enum")) {
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
  message += schema.value("description", "") + "\n\n";
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

static bool arg_to_json(cli_json& js, const cli_json& schema,
    const string& name, bool positional, const vector<string>& args,
    size_t& idx, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  if (schema.contains("enum")) {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    js = args[idx++];
  } else if (schema.value("type", "") == "integer") {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    auto end = (char*)nullptr;
    js       = strtol(args[idx++].c_str(), &end, 10);
    if (end == nullptr) return cli_error("bad value for " + name);
  } else if (schema.value("type", "") == "number") {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    auto end = (char*)nullptr;
    js       = strtod(args[idx++].c_str(), &end);
    if (end == nullptr) return cli_error("bad value for " + name);
  } else if (schema.value("type", "") == "boolean") {
    if (positional) {
      if (idx >= args.size()) return cli_error("missing value for " + name);
      js = args[idx++] == "true" ? true : false;
    } else {
      js = true;
    }
  } else if (schema.value("type", "") == "string") {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    js = args[idx++];
  } else if (schema.value("type", "") == "array") {
    js = cli_json::array();
    if (schema.contains("minItems")) {
      if (idx + schema.value("minItems", (size_t)0) >= args.size())
        return cli_error("missing value for " + name);
    }
    auto end = schema.contains("maxItems")
                   ? (idx + (size_t)schema.at("maxItems"))
                   : args.size();
    while (idx < end) {
      if (!arg_to_json(js.emplace_back(), schema.at("items"), name, positional,
              args, idx, error))
        return false;
    }
  } else {
    throw std::runtime_error("unsupported type");
  }
  return true;
}

static bool args_to_json(cli_json& js, const cli_json& schema,
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

  // init
  js = cli_json::object();

  // add things to schema
  auto commands = vector<string>{}, positionals = vector<string>{};
  for (auto& [key, value] : schema.at("properties").items()) {
    if (value.value("type", ""s) == "object") commands.push_back(key);
  }
  if (schema.contains("cli_positionals")) {
    for (auto& key : schema.at("cli_positionals")) {
      positionals.push_back((string)key);
    }
  }

  // current parsing state
  auto positional = 0;

  // parse arguments
  while (idx < args.size()) {
    auto& arg = args[idx++];
    if (arg == "--help") {
      js["help"] = true;
      continue;
    }
    if (arg == "--config") {
      if (idx >= args.size()) return cli_error("missing value for config");
      js["config"] = args[idx++];
      continue;
    }
    auto is_positional = arg.find('-') != 0;
    if (!commands.empty() && is_positional) {
      if (std::find(commands.begin(), commands.end(), arg) != commands.end()) {
        js["command"] = arg;
        if (!args_to_json(
                js[arg], schema.at("properties").at(arg), args, idx, error))
          return false;
        break;
      } else {
        return cli_error("unknown command " + arg);
      }
    } else if (is_positional) {
      if (positional >= positionals.size())
        return cli_error("too many positional arguments");
      auto  name    = positionals[positional++];
      auto& oschema = schema.at("properties").at(name);
      idx--;
      if (!arg_to_json(js[name], oschema, arg, is_positional, args, idx, error))
        return false;
      if (oschema.contains("cli_config") && !js.contains("config")) {
        js["config"] = "try:" +
                       get_try_config(js.at(name), oschema.at("cli_config"));
      }
    } else {
      auto name = string{};
      for (auto& [key, value] : schema.at("properties").items()) {
        if ("--" + key == arg && value.value("type", "") != "object") {
          name = key;
          break;
        }
      }
      if (name == "") return cli_error("unknown option " + arg);
      auto& oschema = schema.at("properties").at(name);
      if (!arg_to_json(
              js[name], oschema, name, is_positional, args, idx, error))
        return false;
      if (oschema.contains("cli_config") && !js.contains("config")) {
        js["config"] = "try:" +
                       get_try_config(js.at(name), oschema.at("cli_config"));
      }
    }
  }

  // done
  return true;
}

static bool validate_json(const cli_json& js, const cli_json& schema,
    const string& name, bool check_required, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  if (js.is_null()) {
    return cli_error("bad value for " + name);
  } else if (js.is_number_integer()) {
    if (schema.at("type") != "integer" && schema.at("type") != "number")
      return cli_error("bad value for " + name);
  } else if (js.is_number()) {
    if (schema.at("type") != "number")
      return cli_error("bad value for " + name);
  } else if (js.is_boolean()) {
    if (schema.at("type") != "boolean")
      return cli_error("bad value for " + name);
  } else if (js.is_string()) {
    if (schema.at("type") != "string")
      return cli_error("bad value for " + name);
  } else if (js.is_array()) {
    if (schema.at("type") != "array") return cli_error("bad value for " + name);
    if (schema.contains("minItems") &&
        (size_t)schema.at("minItems") > js.size())
      return cli_error("bad value for " + name);
    if (schema.contains("maxItems") &&
        (size_t)schema.at("maxItems") < js.size())
      return cli_error("bad value for " + name);
    for (auto& jsv : js)
      if (!validate_json(jsv, schema.at("items"), name, false, error))
        return false;
  } else if (js.is_object()) {
    if (schema.at("type") != "object")
      return cli_error("bad value for " + name);
    for (auto& [key, jsv] : js.items()) {
      if (key == "help") {
        if (!jsv.is_boolean()) return cli_error("bad value for " + key);
      } else if (key == "config") {
        if (!jsv.is_string()) return cli_error("bad value for " + key);
      } else if (key == "command") {
        if (!jsv.is_string()) return cli_error("bad value for " + key);
      } else {
        if (!schema.at("properties").contains(key))
          return cli_error("unknown option " + key);
        if (!validate_json(jsv, schema.at("properties").at(key), key,
                check_required && js.value("command", "") == key, error))
          return false;
      }
    }
    if (check_required && schema.contains("required")) {
      for (auto& req : schema.at("required")) {
        if (!js.contains((string)req))
          return cli_error("missing value for " + (string)req);
      }
    }
  } else {
    return cli_error("bad value for " + name);
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

// set variables
static bool json_to_variables(
    const cli_json& js, const vector<cli_variable>& variables, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };
  for (auto& variable : variables) {
    if (!js.contains(cli_json::json_pointer{"/" + variable.path})) continue;
    if (!variable.setter(&js.at(cli_json::json_pointer{"/" + variable.path}),
            variable.value, variable.choices))
      return cli_error("bad value for " + variable.path);
  }
  return true;
}

// grabs a configution and update json arguments
static bool config_to_json(cli_json& js, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };
  auto get_config = [](const cli_json& jargs) -> string {
    auto current = &jargs;
    while (true) {
      if (current->contains("config")) return (string)current->at("config");
      if (current->contains("command")) {
        current = &current->at((string)current->at("command"));
      } else {
        break;
      }
    }
    return "";
  };

  auto try_config = false;
  auto config     = get_config(js);
  if (config.empty()) return true;
  if (config.find("try:") == 0) {
    config     = config.substr(4);
    try_config = true;
  }

  auto cjs = cli_json{};
  auto fs  = std::ifstream(config);
  if (!fs) {
    if (try_config) return true;
    return cli_error("missing configuration file " + config);
  }
  try {
    cjs = cli_json::parse(fs);
  } catch (...) {
    return cli_error("error loading configuration " + config);
  }

  update_json_objects(cjs, js);
  js = cjs;
  return true;
}

// initialize a command line parser
cli_state make_cli(const string& name, const string& usage) {
  auto defaults         = cli_json::object();
  auto schema           = cli_json::object();
  schema["cli_name"]    = name;
  schema["description"] = usage;
  schema["type"]        = "object";
  schema["properties"]  = cli_json::object();
  auto cli              = cli_state{};
  // we force the call to copy constructors since there is a bug in the JSON
  // library causing GCC and clang to behave differently with copy constructors
  // initialized with {} --- original code commented here
  // cli.defaults = {
  //     new cli_json{defaults}, [](void* json) { delete (cli_json*)json; }};
  // cli.schema = {
  //     new cli_json{schema}, [](void* json) { delete (cli_json*)json; }};
  cli.defaults = {
      new cli_json(defaults), [](void* json) { delete (cli_json*)json; }};
  cli.schema = {
      new cli_json(schema), [](void* json) { delete (cli_json*)json; }};
  return cli;
}

// add command
cli_command add_command(
    const cli_command& cli, const string& name, const string& usage) {
  auto& defaults                            = get_defaults(cli);
  defaults[name]                            = cli_json::object();
  auto& schema                              = get_schema(cli);
  schema["properties"][name]                = cli_json::object();
  schema["properties"][name]["cli_name"]    = name;
  schema["properties"][name]["description"] = usage;
  schema["properties"][name]["type"]        = "object";
  schema["properties"][name]["properties"]  = cli_json::object();
  return {cli.state, cli.path.empty() ? name : (cli.path + "/" + name)};
}

void set_command_var(const cli_command& cli, string& value) {
  cli.state->variables.push_back(
      cli_variable{cli.path.empty() ? "command" : (cli.path + "/command"),
          &value, cli_from_json_<string>, {}});
}

void set_help_var(const cli_command& cli, bool& value) {
  cli.state->variables.push_back(
      cli_variable{cli.path.empty() ? "help" : (cli.path + "/help"), &value,
          cli_from_json_<string>, {}});
}

void add_command_name(const cli_command& cli, const string& name, string& value,
    const string& usage) {
  cli.state->variables.push_back(
      cli_variable{cli.path.empty() ? "command" : (cli.path + "/command"),
          &value, cli_from_json_<string>, {}});
}

static bool parse_cli(cli_state& cli, const vector<string>& args,
    cli_json& jargs, string& error) {
  auto idx = (size_t)1;
  if (!args_to_json(jargs, get_schema(cli), args, idx, error)) return false;
  if (!config_to_json(jargs, error)) return false;
  if (!validate_json(jargs, get_schema(cli), "", true, error)) return false;
  if (!json_to_variables(jargs, cli.variables, error)) return false;
  return true;
}

bool parse_cli(cli_state& cli, const vector<string>& args, string& error) {
  auto jargs = cli_json{};
  return parse_cli(cli, args, jargs, error);
}

void parse_cli(cli_state& cli, const vector<string>& args) {
  auto get_command = [](const cli_json& jargs) -> string {
    auto command = string{};
    auto current = &jargs;
    while (true) {
      if (!current->contains("command")) break;
      command += (command.empty() ? "" : "/") + current->value("command", "");
      current = &current->at(current->value("command", ""));
    }
    return command;
  };
  auto get_help = [&](const cli_json& jargs) -> bool {
    auto current = &jargs;
    while (true) {
      if (current->contains("help")) return (bool)current->at("help");
      if (current->contains("command")) {
        current = &current->at((string)current->at("command"));
      } else {
        break;
      }
    }
    return false;
  };
  auto jargs = cli_json{};
  auto error = string{};
  if (!parse_cli(cli, args, jargs, error)) {
    print_info("error: " + error);
    print_info("");
    print_info(schema_to_usage(get_schema(cli), get_command(jargs), args[0]));
    exit(1);
  } else if (get_help(jargs)) {
    print_info(schema_to_usage(get_schema(cli), get_command(jargs), args[0]));
    exit(0);
  }
}

void parse_cli(cli_state& cli, int argc, const char** argv) {
  parse_cli(cli, vector<string>{argv, argv + argc});
}

}  // namespace yocto
