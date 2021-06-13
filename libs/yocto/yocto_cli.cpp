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
// FILE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Opens a file with a utf8 file name
static FILE* fopen_utf8(const char* filename, const char* mode) {
#ifdef _WIN32
  auto path8    = std::filesystem::u8path(filename);
  auto str_mode = string{mode};
  auto wmode    = std::wstring(str_mode.begin(), str_mode.end());
  return _wfopen(path8.c_str(), wmode.c_str());
#else
  return fopen(filename, mode);
#endif
}

// Load a text file
static bool load_text(const string& filename, string& str, string& error) {
  // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
  auto fs = fopen_utf8(filename.c_str(), "rb");
  if (!fs) {
    error = filename + ": file not found";
    return false;
  }
  fseek(fs, 0, SEEK_END);
  auto length = ftell(fs);
  fseek(fs, 0, SEEK_SET);
  str.resize(length);
  if (fread(str.data(), 1, length, fs) != length) {
    fclose(fs);
    error = filename + ": read error";
    return false;
  }
  fclose(fs);
  return true;
}

// Save a text file
static bool save_text(
    const string& filename, const string& str, string& error) {
  auto fs = fopen_utf8(filename.c_str(), "wt");
  if (!fs) {
    error = filename + ": file not found";
    return false;
  }
  if (fprintf(fs, "%s", str.c_str()) < 0) {
    fclose(fs);
    error = filename + ": write error";
    return false;
  }
  fclose(fs);
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON I/O
// -----------------------------------------------------------------------------
namespace yocto {

// Load json
bool load_json(const string& filename, json_value& json, string& error) {
  // load file
  auto text = string{};
  if (!load_text(filename, text, error)) return false;
  // parse
  return parse_json(text, json, error);
}

static void from_json(const nlohmann::json& js, json_value& value) {
  switch (js.type()) {
    case nlohmann::json::value_t::null: {
      value.set_null();
    } break;
    case nlohmann::json::value_t::number_integer: {
      value.set_integer((int64_t)js);
    } break;
    case nlohmann::json::value_t::number_unsigned: {
      value.set_uinteger((uint64_t)js);
    } break;
    case nlohmann::json::value_t::number_float: {
      value.set_number((double)js);
    } break;
    case nlohmann::json::value_t::boolean: {
      value.set_boolean((bool)js);
    } break;
    case nlohmann::json::value_t::string: {
      value.set_string((string)js);
    } break;
    case nlohmann::json::value_t::array: {
      value.set_array();
      for (auto& jitem : js) from_json(jitem, value.get_array().emplace_back());
    } break;
    case nlohmann::json::value_t::object: {
      value.set_object();
      for (auto& [key, jitem] : js.items())
        from_json(jitem, value.get_object()[key]);
    } break;
    case nlohmann::json::value_t::binary: {
      value.set_null();
    } break;
    case nlohmann::json::value_t::discarded: {
      value.set_null();
    } break;
  }
}

// Pars json
bool parse_json(const string& text, json_value& json, string& error) {
  try {
    auto njson = nlohmann::json::parse(text);
    from_json(njson, json);
    return true;
  } catch (...) {
    error = "error parsing json";
    return false;
  }
}

// Save json
bool save_json(const string& filename, json_value& json, string& error) {
  // convert to string
  auto text = string{};
  if (!dump_json(text, json, error)) return false;
  // save file
  return save_text(filename, text, error);
}

static void to_json(nlohmann::json& js, const json_value& value) {
  switch (value.type()) {
    case json_type::null: {
      js = {};
    } break;
    case json_type::integer: {
      js = value.get_integer();
    } break;
    case json_type::uinteger: {
      js = value.get_uinteger();
    } break;
    case json_type::number: {
      js = value.get_number();
    } break;
    case json_type::boolean: {
      js = value.get_boolean();
    } break;
    case json_type::string: {
      js = value.get_string();
    } break;
    case json_type::array: {
      js = nlohmann::json::array();
      for (auto& item : value.get_array()) to_json(js.emplace_back(), item);
    } break;
    case json_type::object: {
      js = nlohmann::json::object();
      for (auto& [key, item] : value.get_object()) to_json(js[key], item);
    } break;
  }
}

// Dump json
bool dump_json(string& text, const json_value& json, string& error) {
  try {
    auto njson = nlohmann::json{};
    to_json(njson, json);
    text = njson.dump(2);
    return true;
  } catch (...) {
    error = "error dumping json";
    return false;
  }
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
    json_value& json, const T& value, const vector<string>& choices) {
  if (choices.empty()) {
    json.set(value);
  } else {
    if constexpr (cli_is_array_v<T>) {
      json.set_array(value.size());
      for (auto idx = (size_t)0; idx < value.size(); idx++) {
        cli_to_json(json.at(idx), value.at(idx), choices);
      }
    } else if constexpr (cli_is_vector_v<T>) {
      json.set_array(value.size());
      for (auto idx = (size_t)0; idx < value.size(); idx++) {
        cli_to_json(json.at(idx), value.at(idx), choices);
      }
    } else {
      if constexpr (std::is_integral_v<T> || std::is_floating_point_v<T>) {
        json.set(choices.at((uint64_t)value));
      } else if constexpr (std::is_same_v<T, string>) {
        if (std::find(choices.begin(), choices.end(), value) == choices.end())
          throw std::out_of_range{"bad value"};
        json.set(value);
      } else {
        throw std::invalid_argument{"not supported"};
      }
    }
  }
}

template <typename T>
static void cli_to_schema(json_value& schema, const T& value,
    const vector<string>& choices, const string& name, const string& usage) {
  schema.set_object();
  schema["title"].set(name);
  schema["description"].set(usage);
  cli_to_json(schema["default"], value, choices);
  if constexpr (std::is_same_v<T, int32_t> || std::is_same_v<T, int64_t>) {
    if (!choices.empty()) {
      schema["type"].set("string");
      schema["enum"].set(choices);
    } else {
      schema["type"].set("integer");
    }
  } else if constexpr (std::is_same_v<T, uint32_t> ||
                       std::is_same_v<T, uint64_t>) {
    if (!choices.empty()) {
      schema["type"].set("string");
      schema["enum"].set(choices);
    } else {
      schema["type"].set("integer");
    }
  } else if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>) {
    if (!choices.empty()) {
      schema["type"].set("string");
      schema["enum"].set(choices);
    } else {
      schema["type"].set("number");
    }
  } else if constexpr (std::is_same_v<T, bool>) {
    if (!choices.empty()) {
      schema["type"].set("string");
      schema["enum"].set(choices);
    } else {
      schema["type"].set("boolean");
    }
  } else if constexpr (std::is_same_v<T, string>) {
    if (!choices.empty()) {
      schema["type"].set("string");
      schema["enum"].set(choices);
    } else {
      schema["type"].set("string");
    }
  } else if constexpr (cli_is_array_v<T>) {
    schema["type"].set("array");
    schema["min_items"].set(value.size());
    schema["max_items"].set(value.size());
    cli_to_schema(schema["item"], T{}, choices, "item", "");
  } else if constexpr (cli_is_vector_v<T>) {
    schema["type"].set("array");
    schema["min_items"].set((size_t)0);
    schema["max_items"].set(std::numeric_limits<size_t>::max());
    cli_to_schema(schema["item"], T{}, choices, "item", "");
  }
}

template <typename T>
static void cli_from_json(
    const json_value& json, T& value, const vector<string>& choices) {
  if (choices.empty()) {
    value = json.get<T>();
  } else {
    if constexpr (cli_is_array_v<T>) {
      if (json.size() != value.size()) throw json_error{"bad array size"};
      for (auto idx = (size_t)0; idx < value.size(); idx++) {
        cli_from_json(json.at(idx), value.at(idx), choices);
      }
    } else if constexpr (cli_is_vector_v<T>) {
      value.resize(json.size());
      for (auto idx = (size_t)0; idx < value.size(); idx++) {
        cli_from_json(json.at(idx), value.at(idx), choices);
      }
    } else {
      auto values = json.get<string>();
      if (std::find(choices.begin(), choices.end(), values) == choices.end())
        throw json_error{"invalid label"};
      if constexpr (std::is_integral_v<T> || std::is_floating_point_v<T>) {
        value = (T)(std::find(choices.begin(), choices.end(), values) -
                    choices.begin());
      } else if constexpr (std::is_same_v<T, string>) {
        value = values;
      } else {
        throw std::invalid_argument{"not supported"};
      }
    }
  }
}

template <typename T>
static bool cli_from_json_(
    const json_value& cvalue, void* value, const vector<string>& choices) {
  try {
    cli_from_json(cvalue, *(T*)value, choices);
    return true;
  } catch (...) {
    return false;
  }
}

static json_value& get_defaults(const cli_command& cli) {
  if (cli.path.empty())
    return cli.state->defaults;
  else
    return cli.state->defaults.get_object().at(cli.path);
}
static json_value& get_schema(const cli_command& cli) {
  if (cli.path.empty())
    return cli.state->schema;
  else
    return cli.state->schema.at("properties").at(cli.path);
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
  cli_to_json(defaults.insert(name), value, choices);
  cli_to_schema(
      schema.at("properties").insert(name), value, choices, name, usage);
  if (req) schema.at("required").append(name);
  if (positional) schema.at("clipositional").append(name);
  variables.variables[name] = {&value, cli_from_json_<T>, choices};
}

template <typename T, size_t N>
static void add_option_impl(const cli_command& cli, const string& name,
    array<T, N>& value, const string& usage, const vector<T>& minmax,
    const vector<string>& choices, const string& alt, bool req,
    bool positional) {
  auto& defaults  = get_defaults(cli);
  auto& schema    = get_schema(cli);
  auto& variables = get_variables(cli);
  cli_to_json(defaults.insert(name), value, choices);
  cli_to_schema(
      schema.at("properties").insert(name), value, choices, name, usage);
  if (req) schema.at("required").append(name);
  if (positional) schema.at("clipositional").append(name);
  variables.variables[name] = {&value, cli_from_json_<array<T, N>>, choices};
}

template <typename T>
static void add_option_impl(const cli_command& cli, const string& name,
    vector<T>& value, const string& usage, const vector<T>& minmax,
    const vector<string>& choices, const string& alt, bool req,
    bool positional) {
  auto& defaults  = get_defaults(cli);
  auto& schema    = get_schema(cli);
  auto& variables = get_variables(cli);
  cli_to_json(defaults.insert(name), value, choices);
  cli_to_schema(
      schema.at("properties").insert(name), value, choices, name, usage);
  if (req) schema.at("required").append(name);
  if (positional) schema.at("clipositional").append(name);
  variables.variables[name] = {&value, cli_from_json_<vector<T>>, choices};
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
    get_schema(cli).at("properties").at(name)["cliconfig"].set(config);
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
    get_schema(cli).at("properties").at(name)["cliconfig"].set(config);
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

static string schema_to_usage(
    const json_value& schema_, const string& command, const string& program) {
  auto default_to_string = [](const json_value& value) -> string {
    switch (value.type()) {
      case json_type::null: return "";
      case json_type::integer:
        return "[" + std::to_string(value.get_integer()) + "]";
      case json_type::uinteger:
        return "[" + std::to_string(value.get_uinteger()) + "]";
      case json_type::number:
        return "[" + std::to_string(value.get_number()) + "]";
      case json_type::boolean:
        return "[" + (value.get_boolean() ? string{"true"} : string{"false"}) +
               "]";
      case json_type::string: return "[" + value.get_string() + "]";
      case json_type::array: return "[]";
      case json_type::object: return "";
    }
  };
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
  for (auto& [key, property] : schema.at("properties").get_object()) {
    if (property.at("type").get<string>() == "object") {
      auto line = "  " + key;
      while (line.size() < 32) line += " ";
      line += property.at("description").get<string>() + "\n";
      usage_commands += line;
    } else {
      auto is_positional = false;
      for (auto& positional : schema.at("clipositional").get_array()) {
        if (positional.get<string>() == key) is_positional = true;
      }
      auto line = (is_positional ? "  " : "  --") + key;
      if (property.at("type").get<string>() != "boolean" || is_positional) {
        line += " <" + property.at("type").get<string>() + ">";
      }
      while (line.size() < 32) line += " ";
      line += property.at("description").get<string>();
      if (!property.at("default").is_null()) {
        line += " " + default_to_string(property.at("default"));
      }
      line += "\n";
      if (property.contains("enum")) {
        line += "    with choices: ";
        auto len = 16;
        for (auto& choice : property.at("enum").get_array()) {
          if (len + choice.get<string>().size() + 2 > 78) {
            line += "\n                  ";
            len = 16;
          }
          line += choice.get<string>() + ", ";
          len += (int)choice.get<string>().size() + 2;
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
  usage_options += "  --config <string>             Load configuration. []\n";
  message += "usage: " + progname + (command.empty() ? "" : (" " + command)) +
             (!usage_commands.empty() ? " command" : "") +
             (!usage_options.empty() ? " [options]" : "") +
             (!usage_arguments.empty() ? " <arguments>" : "") + "\n";
  message += schema.at("description").get<string>() + "\n\n";
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

static bool arg_to_json(json_value& value, const json_value& schema,
    const string& name, bool positional, const vector<string>& args,
    size_t& idx, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  if (schema.contains("enum")) {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    value.set_string(args[idx++]);
  } else if (schema.at("type").get<string>() == "integer") {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    auto end = (char*)nullptr;
    value.set_integer(strtol(args[idx++].c_str(), &end, 10));
    if (end == nullptr) return cli_error("bad value for " + name);
  } else if (schema.at("type").get<string>() == "uinteger") {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    auto end = (char*)nullptr;
    value.set_uinteger(strtoul(args[idx++].c_str(), &end, 10));
    if (end == nullptr) return cli_error("bad value for " + name);
  } else if (schema.at("type").get<string>() == "number") {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    auto end = (char*)nullptr;
    value.set_number(strtod(args[idx++].c_str(), &end));
    if (end == nullptr) return cli_error("bad value for " + name);
  } else if (schema.at("type").get<string>() == "boolean") {
    if (positional) {
      if (idx >= args.size()) return cli_error("missing value for " + name);
      value.set_boolean(args[idx++] == "true" ? true : false);
    } else {
      value.set_boolean(true);
    }
  } else if (schema.at("type").get<string>() == "string") {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    value.set_string(args[idx++]);
  } else if (schema.at("type").get<string>() == "array") {
    value.set_array();
    if (idx + schema.at("min_items").get<size_t>() >= args.size())
      return cli_error("missing value for " + name);
    auto end = std::min(
        idx + schema.at("max_items").get<size_t>(), args.size());
    while (idx < end) {
      if (!arg_to_json(value.get_array().emplace_back(), schema.at("items"),
              name, positional, args, idx, error))
        return false;
    }
  } else {
    throw std::runtime_error("unsupported type");
  }
  return true;
}

static bool args_to_json(json_value& value, const json_value& schema,
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
  value.set_object();

  // add things to schema
  auto commands = vector<string>{}, positionals = vector<string>{};
  for (auto& [key, property] : schema.at("properties").get_object()) {
    if (property.at("type").get<string>() == "object") commands.push_back(key);
  }
  if (!schema.at("clipositional").empty()) {
    for (auto& key : schema.at("clipositional").get_array()) {
      positionals.push_back(key.get<string>());
    }
  }

  // current parsing state
  auto positional = 0;

  // parse arguments
  while (idx < args.size()) {
    auto& arg = args[idx++];
    if (arg == "--help") {
      value.get_object()["help"].set_boolean(true);
      continue;
    }
    if (arg == "--config") {
      if (idx >= args.size()) return cli_error("missing value for config");
      value["config"].set_string(args[idx++]);
      continue;
    }
    auto is_positional = arg.find('-') != 0;
    if (!commands.empty() && is_positional) {
      if (std::find(commands.begin(), commands.end(), arg) != commands.end()) {
        value["command"].set(arg);
        if (!args_to_json(value.insert(arg), schema.at("properties").at(arg),
                args, idx, error))
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
      if (!arg_to_json(value.insert(name), oschema, arg, is_positional, args,
              idx, error))
        return false;
      if (oschema.contains("cliconfig") &&
          !contains(value.get_object(), "config")) {
        value.get_object()["config"].set(
            "try:" + get_try_config(value.at(name).get<string>(),
                         oschema.at("cliconfig").get<string>()));
      }
    } else {
      auto name = string{};
      for (auto& [key, value] : schema.at("properties").get_object()) {
        if ("--" + key == arg && !value.is_object()) {
          name = key;
          break;
        }
      }
      if (name == "") return cli_error("unknown option " + arg);
      auto& oschema = schema.at("properties").at(name);
      if (!arg_to_json(value.get_object()[name], oschema, name, is_positional,
              args, idx, error))
        return false;
      if (!oschema.at("cliconfig").empty() &&
          !contains(value.get_object(), "config")) {
        value.get_object()["config"].set(
            "try:" + get_try_config(value.at(name).get<string>(),
                         oschema.at("cliconfig").get<string>()));
      }
    }
  }

  // done
  return true;
}

static bool validate_json(const json_value& value, const json_value& schema,
    const string& name, bool check_required, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  auto contains = [](const auto& object, const string& name) {
    return object.find(name) != object.end();
  };

  switch (value.type()) {
    case json_type::null: {
      return cli_error("bad value for " + name);
    } break;
    case json_type::integer: {
      if (schema.at("type").get<string>() != "number" &&
          schema.at("type").get<string>() != "integer")
        return cli_error("bad value for " + name);
      return true;
    } break;
    case json_type::uinteger: {
      if (schema.at("type").get<string>() != "number" &&
          schema.at("type").get<string>() != "integer")
        return cli_error("bad value for " + name);
      return true;
    } break;
    case json_type::number: {
      if (schema.at("type").get<string>() != "number" &&
          schema.at("type").get<string>() != "integer")
        return cli_error("bad value for " + name);
      return true;
    } break;
    case json_type::boolean: {
      if (schema.at("type").get<string>() != "boolean")
        return cli_error("bad value for " + name);
      return true;
    } break;
    case json_type::string: {
      if (schema.at("type").get<string>() != "string")
        return cli_error("bad value for " + name);
      return true;
    } break;
    case json_type::array: {
      if (schema.at("type").get<string>() != "array")
        return cli_error("bad value for " + name);
      if (schema.at("min_items").get<size_t>() > value.size())
        return cli_error("bad value for " + name);
      if (schema.at("max_items").get<size_t>() < value.size())
        return cli_error("bad value for " + name);
      for (auto& item : value.get_array())
        if (!validate_json(item, schema.at("items"), name, false, error))
          return false;
      return true;
    } break;
    case json_type::object: {
      if (schema.at("type").get<string>() != "object")
        return cli_error("bad value for " + name);
      for (auto& [key, property] : value.get_object()) {
        if (key == "help") {
          if (property.type() != json_type::boolean)
            return cli_error("bad value for " + key);
        } else if (key == "config") {
          if (property.type() != json_type::string)
            return cli_error("bad value for " + key);
        } else if (key == "command") {
          if (property.type() != json_type::string)
            return cli_error("bad value for " + key);
        } else {
          if (!schema.at("properties").contains(key))
            return cli_error("unknown option " + key);
          auto selected_command = contains(value.get_object(), "command") &&
                                  value.at("command").type() ==
                                      json_type::string &&
                                  value.at("command").get_string() == key;
          if (!validate_json(property, schema.at("properties").at(key), key,
                  check_required && selected_command, error))
            return false;
        }
      }
      if (check_required) {
        for (auto& req : schema.at("required").get_array()) {
          if (!contains(value.get_object(), req.get<string>()))
            return cli_error("missing value for " + req.get<string>());
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
void update_value_objects(json_value& value, const json_value& update) {
  auto contains = [](const auto& object, const string& name) {
    return object.find(name) != object.end();
  };

  if (value.type() != json_type::object) return;
  for (auto& [key, property] : update.get_object()) {
    if (property.type() == json_type::object &&
        contains(value.get_object(), key) &&
        value.get_object().at(key).type() == json_type::object) {
      update_value_objects(value.get_object().at(key), value);
    } else {
      value.get_object()[key] = property;
    }
  }
}

// set variables
static bool value_to_variable(
    const json_value& value, cli_variable& variable, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  auto contains = [](const auto& object, const string& name) {
    return object.find(name) != object.end();
  };

  if (variable.setter) {
    if (!variable.setter(value, variable.value, variable.choices)) return false;
  }

  for (auto& [key, property] : variable.variables) {
    if (value.type() != json_type::object)
      throw std::runtime_error{"something went wrong"};
    if (contains(value.get_object(), key)) {
      if (!value_to_variable(value.get_object().at(key), property, error)) {
        return cli_error("bad value for " + key);
      }
    }
  }

  return true;
}

// grabs a configuration and update json arguments
static bool config_to_value(json_value& value, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  auto get_config = [](const json_value& value) -> string {
    auto contains = [](const auto& object, const string& name) {
      return object.find(name) != object.end();
    };

    auto current = &value;
    while (true) {
      if (current->type() != json_type::object) break;
      if (contains(current->get_object(), "config") &&
          current->get_object().at("config").type() == json_type::string)
        return current->get_object().at("config").get_string();
      if (contains(current->get_object(), "command") &&
          current->get_object().at("command").type() == json_type::string) {
        current = &current->get_object().at(
            current->get_object().at("command").get_string());
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

  auto cvalue = json_value{};
  if (try_config) {
    if (!load_json(config, cvalue, error)) return true;
  } else {
    if (!load_json(config, cvalue, error))
      return cli_error("error converting configuration " + config);
  }

  update_value_objects(cvalue, value);
  value = cvalue;
  return true;
}

// initialize a command line parser
cli_state make_cli(const string& name, const string& usage) {
  auto cli = cli_state{};
  cli.defaults.set_object();
  cli.schema.set_object();
  cli.schema["title"].set(name);
  cli.schema["description"].set(usage);
  cli.schema["type"].set("object");
  cli.schema["properties"].set_object();
  cli.schema["required"].set_array();
  cli.schema["clipositional"].set_array();
  return cli;
}

// add command
cli_command add_command(
    const cli_command& cli, const string& name, const string& usage) {
  auto& defaults = get_defaults(cli);
  defaults[name].set_object();
  auto& schema   = get_schema(cli);
  auto& property = schema.at("properties")[name];
  property.set_object();
  property["title"].set(name);
  property["description"].set(usage);
  property["type"].set("object");
  property["properties"].set_object();
  property["required"].set_array();
  property["clipositional"].set_array();
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

static bool parse_cli(cli_state& cli, const vector<string>& args,
    json_value& value, string& error) {
  auto idx = (size_t)1;
  if (!args_to_json(value, get_schema(cli), args, idx, error)) return false;
  if (!config_to_value(value, error)) return false;
  if (!validate_json(value, get_schema(cli), "", true, error)) return false;
  if (!value_to_variable(value, cli.variables, error)) return false;
  return true;
}

bool parse_cli(cli_state& cli, const vector<string>& args, string& error) {
  auto values = json_value{};
  return parse_cli(cli, args, values, error);
}

void parse_cli(cli_state& cli, const vector<string>& args) {
  auto get_command = [](const json_value& value) -> string {
    auto contains = [](const auto& object, const string& name) {
      return object.find(name) != object.end();
    };

    auto command = string{};
    auto current = &value;
    while (true) {
      if (!contains(current->get_object(), "command")) break;
      command += (command.empty() ? "" : "/") +
                 current->get_object().at("command").get_string();
      current = &current->get_object().at(
          current->get_object().at("command").get_string());
    }
    return command;
  };
  auto get_help = [&](const json_value& value) -> bool {
    auto contains = [](const auto& object, const string& name) {
      return object.find(name) != object.end();
    };

    auto current = &value;
    while (true) {
      if (contains(current->get_object(), "help"))
        return current->get_object().at("help").get_boolean();
      if (contains(current->get_object(), "command")) {
        current = &current->get_object().at(
            current->get_object().at("command").get_string());
      } else {
        break;
      }
    }
    return false;
  };
  auto values = json_value{};
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
