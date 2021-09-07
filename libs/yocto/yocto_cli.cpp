//
// Implementation for Yocto/Cli
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

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>

// -----------------------------------------------------------------------------
// PRINT/FORMATTING UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Print a message to the console
void print_info(const string& message) { printf("%s\n", message.c_str()); }
void print_info(const char* message) { printf("%s\n", message); }
// Prints a message to the console and exit with an error.
void print_fatal(const string& message) {
  printf("\n%s\n", message.c_str());
  exit(1);
}
void print_fatal(const char* message) {
  printf("\n%s\n", message);
  exit(1);
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
static void cli_to_json(
    json_value_& json, const T& value, const vector<string>& choices) {
  if (choices.empty()) {
    json = value;
  } else {
    if constexpr (cli_is_array_v<T>) {
      json = json_value_::array();
      for (auto& item : value) {
        cli_to_json(json.emplace_back(), item, choices);
      }
    } else if constexpr (cli_is_vector_v<T>) {
      json = json_value_::array();
      for (auto& item : value) {
        cli_to_json(json.emplace_back(), item, choices);
      }
    } else {
      if constexpr (std::is_integral_v<T> || std::is_floating_point_v<T>) {
        json = choices.at((uint64_t)value);
      } else if constexpr (std::is_same_v<T, string>) {
        if (std::find(choices.begin(), choices.end(), value) == choices.end())
          throw std::invalid_argument{"bad label"};
        json = value;
      } else {
        throw std::invalid_argument{"not supported"};
      }
    }
  }
}

template <typename T>
static void cli_to_schema(json_value_& schema, const T& value,
    const vector<string>& choices, const string& name, const string& usage) {
  schema                = json_value_::object();
  schema["title"]       = name;
  schema["description"] = usage;
  cli_to_json(schema["default"], value, choices);
  if constexpr (std::is_same_v<T, int32_t> || std::is_same_v<T, int64_t>) {
    if (!choices.empty()) {
      schema["type"] = "string";
      schema["enum"] = choices;
    } else {
      schema["type"] = "integer";
    }
  } else if constexpr (std::is_same_v<T, uint32_t> ||
                       std::is_same_v<T, uint64_t>) {
    if (!choices.empty()) {
      schema["type"] = "string";
      schema["enum"] = choices;
    } else {
      schema["type"] = "integer";
    }
  } else if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>) {
    if (!choices.empty()) {
      schema["type"] = "string";
      schema["enum"] = choices;
    } else {
      schema["type"] = "number";
    }
  } else if constexpr (std::is_same_v<T, bool>) {
    if (!choices.empty()) {
      schema["type"] = "string";
      schema["enum"] = choices;
    } else {
      schema["type"] = "boolean";
    }
  } else if constexpr (std::is_same_v<T, string>) {
    if (!choices.empty()) {
      schema["type"] = "string";
      schema["enum"] = choices;
    } else {
      schema["type"] = "string";
    }
  } else if constexpr (cli_is_array_v<T>) {
    using V            = typename T::value_type;
    schema["type"]     = "array";
    schema["minItems"] = value.size();
    schema["maxItems"] = value.size();
    cli_to_schema(schema["items"], V{}, choices, "item", "");
  } else if constexpr (cli_is_vector_v<T>) {
    using V            = typename T::value_type;
    schema["type"]     = "array";
    schema["minItems"] = (size_t)0;
    schema["maxItems"] = (size_t)(1024 * 1024);
    cli_to_schema(schema["items"], V{}, choices, "item", "");
  }
}

template <typename T>
[[nodiscard]] static bool cli_from_json(
    const json_value_& json, T& value, const vector<string>& choices) {
  if (choices.empty()) {
    try {
      value = (T)json;
      return true;
    } catch (...) {
      return false;
    }
  } else {
    if constexpr (cli_is_array_v<T>) {
      if (json.size() != value.size()) return false;
      for (auto idx = (size_t)0; idx < value.size(); idx++) {
        if (!cli_from_json(json[idx], value[idx], choices)) return false;
      }
      return true;
    } else if constexpr (cli_is_vector_v<T>) {
      value.resize(json.size());
      for (auto idx = (size_t)0; idx < value.size(); idx++) {
        if (!cli_from_json(json[idx], value[idx], choices)) return false;
      }
      return true;
    } else {
      auto values = (string)json;
      if (std::find(choices.begin(), choices.end(), values) == choices.end())
        return false;
      if constexpr (std::is_integral_v<T> || std::is_floating_point_v<T>) {
        value = (T)(std::find(choices.begin(), choices.end(), values) -
                    choices.begin());
        return true;
      } else if constexpr (std::is_same_v<T, string>) {
        value = values;
        return true;
      } else {
        throw std::invalid_argument{"not supported"};
      }
    }
  }
}

template <typename T>
[[nodiscard]] static bool cli_from_json_(
    const json_value_& cvalue, void* value, const vector<string>& choices) {
  return cli_from_json(cvalue, *(T*)value, choices);
}

static json_value_& get_defaults(const cli_command& cli) {
  if (cli.path.empty())
    return cli.state->defaults;
  else
    return cli.state->defaults[cli.path];
}
static json_value_& get_value(const cli_command& cli) {
  if (cli.path.empty())
    return cli.state->value;
  else
    return cli.state->value[cli.path];
}
static json_value_& get_schema(const cli_command& cli) {
  if (cli.path.empty())
    return cli.state->schema;
  else
    return cli.state->schema["properties"].at(cli.path);
}
static cli_variable& get_variables(const cli_command& cli) {
  if (cli.path.empty())
    return cli.state->variables;
  else
    return cli.state->variables.variables.at(cli.path);
}

template <typename T>
static void add_option_impl(const cli_command& cli, const string& name,
    T& value, const string& usage, const vector<T>& minmax,
    const vector<string>& choices, const string& alt, bool req,
    bool positional) {
  auto& defaults  = get_defaults(cli);
  auto& schema    = get_schema(cli);
  auto& variables = get_variables(cli);
  cli_to_json(defaults[name], value, choices);
  cli_to_schema(schema["properties"][name], value, choices, name, usage);
  if (req) schema["required"].push_back(name);
  if (positional) schema["clipositional"].push_back(name);
  variables.variables[name] = {&value, cli_from_json_<T>, choices, {}};
}

template <typename T, size_t N>
static void add_option_impl(const cli_command& cli, const string& name,
    array<T, N>& value, const string& usage, const vector<T>& minmax,
    const vector<string>& choices, const string& alt, bool req,
    bool positional) {
  auto& defaults  = get_defaults(cli);
  auto& schema    = get_schema(cli);
  auto& variables = get_variables(cli);
  cli_to_json(defaults[name], value, choices);
  cli_to_schema(schema["properties"][name], value, choices, name, usage);
  if (req) schema["required"].push_back(name);
  if (positional) schema["clipositional"].push_back(name);
  variables.variables[name] = {
      &value, cli_from_json_<array<T, N>>, choices, {}};
}

template <typename T>
static void add_option_impl(const cli_command& cli, const string& name,
    vector<T>& value, const string& usage, const vector<T>& minmax,
    const vector<string>& choices, const string& alt, bool req,
    bool positional) {
  auto& defaults  = get_defaults(cli);
  auto& schema    = get_schema(cli);
  auto& variables = get_variables(cli);
  cli_to_json(defaults[name], value, choices);
  cli_to_schema(schema["properties"][name], value, choices, name, usage);
  if (req) schema["required"].push_back(name);
  if (positional) schema["clipositional"].push_back(name);
  variables.variables[name] = {&value, cli_from_json_<vector<T>>, choices, {}};
}

// Add an optional argument. Supports strings, numbers, and boolean flags.
// Optional arguments will be parsed with name `--<name>` and `-<alt>`.
// Optional booleans will support both `--<name>` and `--no-<name>` to enabled
// and disable the flag.
void add_option(const cli_command& cli, const string& name, int& value,
    const string& usage, const vector<int>& minmax, const string& alt,
    bool req) {
  add_option_impl(cli, name, value, usage, minmax, {}, alt, req, false);
}
void add_option(const cli_command& cli, const string& name, float& value,
    const string& usage, const vector<float>& minmax, const string& alt,
    bool req) {
  add_option_impl(cli, name, value, usage, minmax, {}, alt, req, false);
}
void add_option(const cli_command& cli, const string& name, bool& value,
    const string& usage, const vector<string>& choices, const string& alt,
    bool req) {
  add_option_impl(
      cli, name, value, usage, vector<bool>{}, choices, alt, req, false);
}
void add_option(const cli_command& cli, const string& name, string& value,
    const string& usage, const vector<string>& choices, const string& alt,
    bool req) {
  add_option_impl(cli, name, value, usage, {}, choices, alt, req, false);
}
void add_option_with_config(const cli_command& cli, const string& name,
    string& value, const string& usage, const string& config, const string& alt,
    bool req) {
  add_option_impl(cli, name, value, usage, {}, {}, alt, req, false);
  if (!config.empty()) {
    get_schema(cli)["properties"][name]["cliconfig"] = config;
  }
}
void add_option(const cli_command& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, const string& alt,
    bool req) {
  add_option_impl(cli, name, value, usage, {}, choices, alt, req, false);
}
// Add a positional argument. Supports strings, numbers, and boolean flags.
void add_argument(const cli_command& cli, const string& name, int& value,
    const string& usage, const vector<int>& minmax, bool req) {
  add_option_impl(cli, name, value, usage, minmax, {}, "", req, false);
}
void add_argument(const cli_command& cli, const string& name, float& value,
    const string& usage, const vector<float>& minmax, bool req) {
  add_option_impl(cli, name, value, usage, minmax, {}, "", req, true);
}
void add_argument(const cli_command& cli, const string& name, bool& value,
    const string& usage, const vector<string>& choices, bool req) {
  add_option_impl(cli, name, value, usage, {}, choices, "", req, true);
}
void add_argument(const cli_command& cli, const string& name, string& value,
    const string& usage, const vector<string>& choices, bool req) {
  add_option_impl(cli, name, value, usage, {}, choices, "", req, true);
}
void add_argument_with_config(const cli_command& cli, const string& name,
    string& value, const string& usage, const string& config, bool req) {
  add_option_impl(cli, name, value, usage, {}, {}, "", req, true);
  if (!config.empty()) {
    get_schema(cli)["properties"][name]["cliconfig"] = config;
  }
}
void add_argument(const cli_command& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, bool req) {
  add_option_impl(cli, name, value, usage, {}, choices, "", req, true);
}
// Add a positional argument that consumes all arguments left.
// Supports strings and enums.
void add_argument(const cli_command& cli, const string& name,
    vector<int>& value, const string& usage, const vector<int>& minmax,
    bool req) {
  add_option_impl(cli, name, value, usage, minmax, {}, "", req, true);
}
void add_argument(const cli_command& cli, const string& name,
    vector<float>& value, const string& usage, const vector<float>& minmax,
    bool req) {
  add_option_impl(cli, name, value, usage, minmax, {}, "", req, true);
}
void add_argument(const cli_command& cli, const string& name,
    vector<int>& value, const string& usage, const vector<string>& choices,
    bool req) {
  add_option_impl(
      cli, name, value, usage, vector<int>{}, choices, "", req, true);
}
void add_argument(const cli_command& cli, const string& name,
    vector<string>& value, const string& usage, const vector<string>& choices,
    bool req) {
  add_option_impl(
      cli, name, value, usage, vector<string>{}, choices, "", req, true);
}

// Add an optional argument. Supports basic math types.
void add_option(const cli_command& cli, const string& name, vec2i& value,
    const string& usage, const vector<int>& minmax, const string& alt,
    bool req) {
  add_option_impl(
      cli, name, (array<int, 2>&)value, usage, minmax, {}, alt, req, false);
}
void add_option(const cli_command& cli, const string& name, vec3i& value,
    const string& usage, const vector<int>& minmax, const string& alt,
    bool req) {
  add_option_impl(
      cli, name, (array<int, 3>&)value, usage, minmax, {}, alt, req, false);
}
void add_option(const cli_command& cli, const string& name, vec4i& value,
    const string& usage, const vector<int>& minmax, const string& alt,
    bool req) {
  add_option_impl(
      cli, name, (array<int, 4>&)value, usage, minmax, {}, alt, req, false);
}
void add_option(const cli_command& cli, const string& name, vec2f& value,
    const string& usage, const vector<float>& minmax, const string& alt,
    bool req) {
  add_option_impl(
      cli, name, (array<float, 2>&)value, usage, minmax, {}, alt, req, false);
}
void add_option(const cli_command& cli, const string& name, vec3f& value,
    const string& usage, const vector<float>& minmax, const string& alt,
    bool req) {
  add_option_impl(
      cli, name, (array<float, 3>&)value, usage, minmax, {}, alt, req, false);
}
void add_option(const cli_command& cli, const string& name, vec4f& value,
    const string& usage, const vector<float>& minmax, const string& alt,
    bool req) {
  add_option_impl(
      cli, name, (array<float, 4>&)value, usage, minmax, {}, alt, req, false);
}

static string schema_to_usage_default(const json_value_& value) {
  if (value.is_object()) return "";
  return value.dump();
}

static string schema_to_usage_typename(const json_value_& schema) {
  auto type = schema.at("type").get<string>();
  if (type == "integer") {
    return "integer";
  } else if (type == "number") {
    return "number";
  } else if (type == "boolean") {
    return "boolean";
  } else if (type == "string") {
    return "string";
  } else if (type == "array") {
    if ((uint64_t)schema.at("minItems") == (uint64_t)schema.at("maxItems"))
      return schema_to_usage_typename(schema.at("items")) + "[" +
             std::to_string((uint64_t)schema.at("minItems")) + "]";
    else
      return schema_to_usage_typename(schema.at("items")) + "[]";
  } else {
    throw std::invalid_argument{"invalid type"};
  }
}

static string schema_to_usage(
    const json_value_& schema, const string& command, const string& program) {
  auto progname = program;
  if (progname.rfind('/') != string::npos) {
    progname = progname.substr(progname.rfind('/') + 1);
  }
  if (progname.rfind('\\') != string::npos) {
    progname = progname.substr(progname.rfind('\\') + 1);
  }
  auto commandname = command;
  for (auto& c : commandname)
    if (c == '/') c = ' ';
  auto message       = string{};
  auto usage_options = string{}, usage_arguments = string{},
       usage_commands = string{};
  for (auto& [key, property] : schema["properties"].items()) {
    if (property["type"] == "object") {
      auto line = "  " + key;
      while (line.size() < 32) line += " ";
      line += (string)property["description"] + "\n";
      usage_commands += line;
    } else {
      auto is_positional = false;
      for (auto& positional : schema["clipositional"]) {
        if (positional == key) is_positional = true;
      }
      auto line = (is_positional ? "  " : "  --") + key;
      if (property["type"] != "boolean" || is_positional) {
        line += " <" + schema_to_usage_typename(property) + ">";
      }
      while (line.size() < 32) line += " ";
      line += (string)property["description"];
      if (!property["default"].is_null()) {
        line += " (" + schema_to_usage_default(property["default"]) + ")";
      }
      line += "\n";
      if (property.contains("enum")) {
        line += "    with choices: ";
        auto len = 16;
        for (auto& choice : property["enum"]) {
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
  usage_options += "  --help                        Prints help. (false)\n";
  usage_options +=
      "  --config <string>             Load configuration. (\"\")\n";
  message += "usage: " + progname +
             (commandname.empty() ? "" : (" " + commandname)) +
             (!usage_commands.empty() ? " command" : "") +
             (!usage_options.empty() ? " [options]" : "") +
             (!usage_arguments.empty() ? " <arguments>" : "") + "\n";
  message += (string)schema["description"] + "\n\n";
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

[[nodiscard]] static bool arg_to_json(json_value_& value,
    const json_value_& schema, const string& name, bool positional,
    const vector<string>& args, size_t& idx, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };
  if (schema.contains("enum")) {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    value = args[idx++];
    return true;
  } else if (schema["type"] == "integer") {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    auto end = (char*)nullptr;
    value    = (int64_t)strtol(args[idx++].c_str(), &end, 10);
    if (end == nullptr) return cli_error("bad value for " + name);
    return true;
  } else if (schema["type"] == "uinteger") {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    auto end = (char*)nullptr;
    value    = (uint64_t)strtoul(args[idx++].c_str(), &end, 10);
    if (end == nullptr) return cli_error("bad value for " + name);
    return true;
  } else if (schema["type"] == "number") {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    auto end = (char*)nullptr;
    value    = strtod(args[idx++].c_str(), &end);
    if (end == nullptr) return cli_error("bad value for " + name);
    return true;
  } else if (schema["type"] == "boolean") {
    if (positional) {
      if (idx >= args.size()) return cli_error("missing value for " + name);
      value = args[idx++] == "true" ? true : false;
    } else {
      value = true;
    }
    return true;
  } else if (schema["type"] == "string") {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    value = args[idx++];
    return true;
  } else if (schema["type"] == "array") {
    value = json_value_::array();
    if (idx + (size_t)schema["minItems"] > args.size())
      return cli_error("missing value for " + name);
    auto end = std::min(idx + (size_t)schema["maxItems"], args.size());
    while (idx < end) {
      if (!arg_to_json(value.emplace_back(), schema["items"], name, positional,
              args, idx, error))
        return false;
    }
    return true;
  } else {
    throw std::runtime_error("unsupported type");
    return false;
  }
}

[[nodiscard]] static bool args_to_json(json_value_& value,
    const json_value_& schema, const vector<string>& args, size_t idx,
    string& error) {
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
  value = json_value_::object();

  // add things to schema
  auto commands = vector<string>{}, positionals = vector<string>{};
  for (auto& [key, property] : schema["properties"].items()) {
    if (property["type"] == "object") commands.push_back(key);
  }
  if (!schema["clipositional"].empty()) {
    for (auto& key : schema["clipositional"]) {
      positionals.push_back((string)key);
    }
  }

  // current parsing state
  auto positional = 0;

  // parse arguments
  while (idx < args.size()) {
    auto& arg = args[idx++];
    if (arg == "--help") {
      value["help"] = true;
      continue;
    }
    if (arg == "--config") {
      if (idx >= args.size()) return cli_error("missing value for config");
      value["config"] = args[idx++];
      continue;
    }
    auto is_positional = arg.find('-') != 0;
    if (!commands.empty() && is_positional) {
      if (std::find(commands.begin(), commands.end(), arg) != commands.end()) {
        value["command"] = arg;
        return args_to_json(
            value[arg], schema["properties"][arg], args, idx, error);
      } else {
        return cli_error("unknown command " + arg);
      }
    } else if (is_positional) {
      if (positional >= positionals.size())
        return cli_error("too many positional arguments");
      auto  name    = positionals[positional++];
      auto& oschema = schema["properties"][name];
      idx--;
      if (!arg_to_json(
              value[name], oschema, arg, is_positional, args, idx, error))
        return false;
      if (oschema.contains("cliconfig") && !value.contains("config")) {
        value["config"] = "try:" + get_try_config((string)value[name],
                                       (string)oschema["cliconfig"]);
      }
    } else {
      auto name = string{};
      for (auto& [key, value] : schema["properties"].items()) {
        if ("--" + key == arg && value.at("type") != "object") {
          name = key;
          break;
        }
      }
      if (name == "") return cli_error("unknown option " + arg);
      auto& oschema = schema["properties"][name];
      if (!arg_to_json(
              value[name], oschema, name, is_positional, args, idx, error))
        return false;
      if (oschema.contains("cliconfig") && !value.contains("config")) {
        value["config"] = "try:" + get_try_config((string)value[name],
                                       (string)oschema["cliconfig"]);
      }
    }
  }

  // done
  return true;
}

static bool validate_json(const json_value_& value, const json_value_& schema,
    const string& name, bool check_required, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  switch (value.type()) {
    case json_value_::value_t::null: {
      return cli_error("bad value for " + name);
    } break;
    case json_value_::value_t::number_integer: {
      if (schema["type"] != "number" && schema["type"] != "integer")
        return cli_error("bad value for " + name);
      return true;
    } break;
    case json_value_::value_t::number_unsigned: {
      if (schema["type"] != "number" && schema["type"] != "integer")
        return cli_error("bad value for " + name);
      return true;
    } break;
    case json_value_::value_t::number_float: {
      if (schema["type"] != "number" && schema["type"] != "integer")
        return cli_error("bad value for " + name);
      return true;
    } break;
    case json_value_::value_t::boolean: {
      if (schema["type"] != "boolean")
        return cli_error("bad value for " + name);
      return true;
    } break;
    case json_value_::value_t::string: {
      if (schema["type"] != "string") return cli_error("bad value for " + name);
      return true;
    } break;
    case json_value_::value_t::array: {
      if (schema["type"] != "array") return cli_error("bad value for " + name);
      if ((size_t)schema["minItems"] > value.size())
        return cli_error("bad value for " + name);
      if ((size_t)schema["maxItems"] < value.size())
        return cli_error("bad value for " + name);
      for (auto& item : value)
        if (!validate_json(item, schema["items"], name, false, error))
          return false;
      return true;
    } break;
    case json_value_::value_t::object: {
      if (schema["type"] != "object") return cli_error("bad value for " + name);
      for (auto& [key, property] : value.items()) {
        if (key == "help") {
          if (!property.is_boolean()) return cli_error("bad value for " + key);
        } else if (key == "config") {
          if (!property.is_string()) return cli_error("bad value for " + key);
        } else if (key == "command") {
          if (!property.is_string()) return cli_error("bad value for " + key);
        } else {
          if (!schema["properties"].contains(key))
            return cli_error("unknown option " + key);
          auto selected_command = value.contains("command") &&
                                  value["command"].is_string() &&
                                  value["command"] == key;
          if (!validate_json(property, schema["properties"][key], key,
                  check_required && selected_command, error))
            return false;
        }
      }
      if (check_required) {
        for (auto& req : schema["required"]) {
          if (!value.contains((string)req))
            return cli_error("missing value for " + (string)req);
        }
      }
      return true;
    } break;
    default: {
      return cli_error("bad value for " + name);
    } break;
  }
}

// update json objects
void update_value_objects(json_value_& value, const json_value_& update) {
  if (!value.is_object()) return;
  for (auto& [key, property] : update.items()) {
    if (property.is_object() && value.contains(key) && value[key].is_object()) {
      update_value_objects(value[key], value);
    } else {
      value[key] = property;
    }
  }
}

// set variables
[[nodiscard]] static bool json_to_variable(const json_value_& json,
    cli_variable& variable, const string& name, string& error) {
  if (variable.setter) {
    if (!variable.setter(json, variable.value, variable.choices)) {
      error = "bad value for " + name;
      return false;
    }
  }
  for (auto& [key, property] : variable.variables) {
    if (json.contains(key)) {
      if (!json_to_variable(json.at(key), property, key, error)) return false;
    }
  }
  return true;
}

// grabs a configuration and update json arguments
static string get_config(const json_value_& json) {
  if (!json.is_object()) return "";
  if (json.contains("config") && json.at("config").is_string())
    return json.at("config").get<string>();
  if (json.contains("command") && json.at("command").is_string()) {
    return get_config(json.at(json.at("command").get<string>()));
  }
  return "";
}

[[nodiscard]] static bool config_to_json(json_value_& json, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  auto filename = get_config(json);
  if (filename.empty()) return true;
  if (filename.find("try:") == 0) {
    filename = filename.substr(4);
    if (!exists(std::filesystem::u8path(filename))) return true;
  }

  auto stream = std::ifstream(filename);
  if (!stream) return cli_error("missing configuration file " + filename);
  auto config = json_value_::parse(stream, nullptr, false);
  if (config.is_discarded())
    return cli_error("error converting configuration " + filename);

  update_value_objects(config, json);
  json = config;
  return true;
}

// initialize a command line parser
cli_state make_cli(const string& name, const string& usage) {
  auto cli                    = cli_state{};
  cli.defaults                = json_value_::object();
  cli.value                   = json_value_::object();
  cli.schema                  = json_value_::object();
  cli.schema["title"]         = name;
  cli.schema["description"]   = usage;
  cli.schema["type"]          = "object";
  cli.schema["properties"]    = json_value_::object();
  cli.schema["required"]      = json_value_::array();
  cli.schema["clipositional"] = json_value_::array();
  return cli;
}

// add command
cli_command add_command(
    const cli_command& cli, const string& name, const string& usage) {
  auto& defaults            = get_defaults(cli);
  defaults[name]            = json_value_::object();
  auto& value               = get_value(cli);
  value[name]               = json_value_::object();
  auto& schema              = get_schema(cli);
  auto& property            = schema["properties"][name];
  property                  = json_value_::object();
  property["title"]         = name;
  property["description"]   = usage;
  property["type"]          = "object";
  property["properties"]    = json_value_::object();
  property["required"]      = json_value_::array();
  property["clipositional"] = json_value_::array();
  auto& variables           = get_variables(cli);
  variables.variables[name] = {};
  return {cli.state, cli.path.empty() ? name : (cli.path + "/" + name)};
}

void set_command_var(const cli_command& cli, string& value) {
  auto& variables                = get_variables(cli);
  variables.variables["command"] = {&value, cli_from_json_<string>, {}};
}

void set_help_var(const cli_command& cli, bool& value) {
  auto& variables             = get_variables(cli);
  variables.variables["help"] = {&value, cli_from_json_<bool>, {}};
}

void add_command_name(const cli_command& cli, const string& name, string& value,
    const string& usage) {
  auto& variables                = get_variables(cli);
  variables.variables["command"] = {&value, cli_from_json_<string>, {}};
}

// get usage
static string get_usage(const json_value_& json, const json_value_& schema,
    const string& commandname, const string& progname) {
  if (json.contains("command") && !json.at("command").empty()) {
    auto command = json.at("command").get<string>();
    return get_usage(json.at(command), schema.at("properties").at(command),
        commandname.empty() ? command : (commandname + " " + command),
        progname);
  } else {
    return schema_to_usage(schema, commandname, progname);
  }
}
string get_usage(const cli_state& cli) {
  return get_usage(
      cli.value, cli.schema, "", cli.schema.at("title").get<string>());
}

// get help
static bool get_help(const json_value_& json, const json_value_& schema) {
  if (json.contains("help")) return json.at("help").get<bool>();
  if (json.contains("command") && !json.at("command").empty()) {
    return get_help(json.at(json.at("command").get<string>()),
        schema.at("properties").at(json.at("command").get<string>()));
  }
  return false;
}
bool get_help(const cli_state& cli) { return get_help(cli.value, cli.schema); }

// get command
static string get_command(const json_value_& json, const json_value_& schema) {
  if (json.contains("command")) {
    auto command   = json.at("command").get<string>();
    auto subcommnd = get_command(json.at(command), schema.at(command));
    return command + (subcommnd.empty() ? "" : ("/" + subcommnd));
  } else {
    return "";
  }
}
string get_command(const cli_state& cli) {
  return get_command(cli.value, cli.schema);
}

bool parse_cli(cli_state& cli, const vector<string>& args, string& error) {
  auto cli_error = [&error, &cli]() {
    error += "\n" + get_usage(cli);
    return false;
  };
  auto idx = (size_t)1;
  if (!args_to_json(cli.value, cli.schema, args, idx, error))
    return cli_error();
  if (!config_to_json(cli.value, error)) return false;
  if (!validate_json(cli.value, cli.schema, "", true, error))
    return cli_error();
  if (!json_to_variable(cli.value, cli.variables, "", error))
    return cli_error();
  if (get_help(cli)) return cli_error();
  return true;
}

// a convenience function that packs args to strings
vector<string> make_cli_args(int argc, const char** argv) {
  return vector<string>(argv, argv + argc);
}

cli_status parse_cli(cli_state& cli, const vector<string>& args) {
  auto error = string{};
  if (!parse_cli(cli, args, error)) return cli_status{error};
  return {};
}

}  // namespace yocto
