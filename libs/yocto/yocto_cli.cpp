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
json_value load_json(const string& filename) {
  auto json = json_value{};
  load_json(filename, json);
  return json;
}
void load_json(const string& filename, json_value& json) {
  auto error = std::string{};
  auto text  = string{};
  if (!load_text(filename, text, error))
    throw json_error{"cannot load", nullptr};
  parse_json(text, json);
}
void save_json(const string& filename, const json_value& json) {
  auto text  = dump_json(json);
  auto error = std::string{};
  if (!save_text(filename, text, error))
    throw json_error{"cannot save", nullptr};
}

[[maybe_unused]] static void from_json(
    const nlohmann::json& js, json_value& value) {
  switch (js.type()) {
    case nlohmann::json::value_t::null: {
      value = nullptr;
    } break;
    case nlohmann::json::value_t::number_integer: {
      value = (int64_t)js;
    } break;
    case nlohmann::json::value_t::number_unsigned: {
      value = (uint64_t)js;
    } break;
    case nlohmann::json::value_t::number_float: {
      value = (double)js;
    } break;
    case nlohmann::json::value_t::boolean: {
      value = (bool)js;
    } break;
    case nlohmann::json::value_t::string: {
      value = (string)js;
    } break;
    case nlohmann::json::value_t::array: {
      value = json_array{};
      for (auto& jitem : js) from_json(jitem, value.emplace_back());
    } break;
    case nlohmann::json::value_t::object: {
      value = json_object{};
      for (auto& [key, jitem] : js.items()) from_json(jitem, value[key]);
    } break;
    case nlohmann::json::value_t::binary: {
      value = nullptr;
    } break;
    case nlohmann::json::value_t::discarded: {
      value = nullptr;
    } break;
  }
}

struct json_sax {
  explicit json_sax(json_value* root, bool allow_exceptions = true)
      : _root(root), _allow_exceptions(allow_exceptions) {}

  bool null() {
    handle_value(nullptr);
    return true;
  }

  bool boolean(bool val) {
    handle_value(val);
    return true;
  }

  bool number_integer(int64_t val) {
    handle_value(val);
    return true;
  }

  bool number_unsigned(uint64_t val) {
    handle_value(val);
    return true;
  }

  bool number_float(double val, const std::string& /*unused*/) {
    handle_value(val);
    return true;
  }

  bool string(std::string& val) {
    handle_value(val);
    return true;
  }

  bool binary(vector<byte>& val) {
    handle_value(std::move(val));
    return true;
  }

  bool start_object(std::size_t len) {
    _stack.push_back(handle_value(json_object{}));
    return true;
  }

  bool key(std::string& val) {
    // add null at given key and store the reference for later
    // _object_item = &(_stack.back()->operator[](val));
    _object_item = &(_stack.back()->insert_back(val));
    return true;
  }

  bool end_object() {
    _stack.pop_back();
    return true;
  }

  bool start_array(std::size_t len) {
    _stack.push_back(handle_value(json_array{}));
    return true;
  }

  bool end_array() {
    _stack.pop_back();
    return true;
  }

  template <class Exception>
  bool parse_error(std::size_t /*unused*/, const std::string& /*unused*/,
      const Exception& ex) {
    _errored = true;
    static_cast<void>(ex);
    if (_allow_exceptions) {
      if (_stack.empty()) throw json_error{"parse error", _root};
      if (_stack.back()->is_array())
        throw json_error{"parse error", _stack.back()};
      if (_stack.back()->is_object())
        throw json_error{"parse error", _stack.back()};
    }
    return false;
  }

  bool is_errored() const { return _errored; }

 private:
  template <typename Value>
  json_value* handle_value(Value&& v) {
    if (_stack.empty()) {
      *_root = json_value(std::forward<Value>(v));
      return _root;
    }

    assert(_stack.back()->is_array() || _stack.back()->is_object());

    if (_stack.back()->is_array()) {
      _stack.back()->emplace_back(std::forward<Value>(v));
      return &(_stack.back()->back());
    }

    assert(_stack.back()->is_object());
    assert(_object_item);
    *_object_item = json_value(std::forward<Value>(v));
    return _object_item;
  }

  json_value*         _root;
  vector<json_value*> _stack            = {};
  json_value*         _object_item      = nullptr;
  bool                _errored          = false;
  bool                _allow_exceptions = true;
};

// Parse json
json_value parse_json(const string& text) {
  auto json = json_value{};
  parse_json(text, json);
  return json;
}
void parse_json(const string& text, json_value& json) {
  json     = json_value{};
  auto sax = json_sax{&json};
  nlohmann::json::sax_parse(text, &sax);
}

// Dump json
static void to_json(nlohmann::json& js, const json_value& value) {
  switch (value.type()) {
    case json_type::null: js = {}; break;
    case json_type::integer: js = (int64_t)value; break;
    case json_type::uinteger: js = (uint64_t)value; break;
    case json_type::number: js = (double)value; break;
    case json_type::boolean: js = (bool)value; break;
    case json_type::string: js = (string)value; break;
    case json_type::array: {
      js = nlohmann::json::array();
      for (auto& item : value) to_json(js.emplace_back(), item);
    } break;
    case json_type::object: {
      js = nlohmann::json::object();
      for (auto& [key, item] : value.items()) to_json(js[key], item);
    } break;
  }
}

// Dump json
string dump_json(const json_value& json) {
  auto text = string{};
  dump_json(text, json);
  return text;
}
void dump_json(string& text, const json_value& json) {
  auto njson = nlohmann::json{};
  to_json(njson, json);
  text = njson.dump(2);
}

// Parse json
bool load_json(const string& filename, json_value& json, string& error) {
  try {
    load_json(filename, json);
    return true;
  } catch (std::exception& exception) {
    error = exception.what();
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

// Parse json
bool parse_json(const string& text, json_value& json, string& error) {
  try {
    parse_json(text, json);
    return true;
  } catch (std::exception& exception) {
    error = exception.what();
    return false;
  }
}

// Dump json
bool dump_json(string& text, const json_value& json, string& error) {
  try {
    dump_json(text, json);
    return true;
  } catch (std::exception& exception) {
    error = exception.what();
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
    json = value;
  } else {
    if constexpr (cli_is_array_v<T>) {
      json = json_array(value.size());
      for (auto idx = (size_t)0; idx < value.size(); idx++) {
        cli_to_json(json[idx], value[idx], choices);
      }
    } else if constexpr (cli_is_vector_v<T>) {
      json = json_array(value.size());
      for (auto idx = (size_t)0; idx < value.size(); idx++) {
        cli_to_json(json[idx], value[idx], choices);
      }
    } else {
      if constexpr (std::is_integral_v<T> || std::is_floating_point_v<T>) {
        json = choices.at((uint64_t)value);
      } else if constexpr (std::is_same_v<T, string>) {
        if (std::find(choices.begin(), choices.end(), value) == choices.end())
          throw std::out_of_range{"bad value"};
        json = value;
      } else {
        throw std::invalid_argument{"not supported"};
      }
    }
  }
}

template <typename T>
static void cli_to_schema(json_value& schema, const T& value,
    const vector<string>& choices, const string& name, const string& usage) {
  schema                = json_object{};
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
    schema["type"]      = "array";
    schema["min_items"] = value.size();
    schema["max_items"] = value.size();
    cli_to_schema(schema["item"], T{}, choices, "item", "");
  } else if constexpr (cli_is_vector_v<T>) {
    schema["type"]      = "array";
    schema["min_items"] = (size_t)0;
    schema["max_items"] = std::numeric_limits<size_t>::max();
    cli_to_schema(schema["item"], T{}, choices, "item", "");
  }
}

template <typename T>
static void cli_from_json(
    const json_value& json, T& value, const vector<string>& choices) {
  if (choices.empty()) {
    value = (T)json;
  } else {
    if constexpr (cli_is_array_v<T>) {
      if (json.size() != value.size())
        throw json_error{"bad array size", &json};
      for (auto idx = (size_t)0; idx < value.size(); idx++) {
        cli_from_json(json[idx], value[idx], choices);
      }
    } else if constexpr (cli_is_vector_v<T>) {
      value.resize(json.size());
      for (auto idx = (size_t)0; idx < value.size(); idx++) {
        cli_from_json(json[idx], value[idx], choices);
      }
    } else {
      auto values = (string)json;
      if (std::find(choices.begin(), choices.end(), values) == choices.end())
        throw json_error{"invalid label", &json};
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
    return cli.state->defaults[cli.path];
}
static json_value& get_schema(const cli_command& cli) {
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
  cli_to_json(defaults[name], value, choices);
  cli_to_schema(schema["properties"][name], value, choices, name, usage);
  if (req) schema["required"].push_back(name);
  if (positional) schema["clipositional"].push_back(name);
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
  cli_to_json(defaults[name], value, choices);
  cli_to_schema(schema["properties"][name], value, choices, name, usage);
  if (req) schema["required"].push_back(name);
  if (positional) schema["clipositional"].push_back(name);
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

static string schema_to_usage(
    const json_value& schema_, const string& command, const string& program) {
  auto default_to_string = [](const json_value& value) -> string {
    switch (value.type()) {
      case json_type::null: return "";
      case json_type::integer:
        return "[" + std::to_string((int64_t)value) + "]";
      case json_type::uinteger:
        return "[" + std::to_string((uint64_t)value) + "]";
      case json_type::number: return "[" + std::to_string((double)value) + "]";
      case json_type::boolean:
        return "[" + ((bool)value ? string{"true"} : string{"false"}) + "]";
      case json_type::string: return "[" + (string)value + "]";
      case json_type::array: return "[]";
      case json_type::object: return "";
    }
  };
  auto& schema   = command.empty() ? schema_ : schema_["properties"][command];
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
        line += " <" + (string)property["type"] + ">";
      }
      while (line.size() < 32) line += " ";
      line += (string)property["description"];
      if (!property["default"].is_null()) {
        line += " " + default_to_string(property["default"]);
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
  usage_options += "  --help                        Prints help. [false]\n";
  usage_options += "  --config <string>             Load configuration. []\n";
  message += "usage: " + progname + (command.empty() ? "" : (" " + command)) +
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

static bool arg_to_json(json_value& value, const json_value& schema,
    const string& name, bool positional, const vector<string>& args,
    size_t& idx, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  if (schema.contains("enum")) {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    value = args[idx++];
  } else if (schema["type"] == "integer") {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    auto end = (char*)nullptr;
    value    = (int64_t)strtol(args[idx++].c_str(), &end, 10);
    if (end == nullptr) return cli_error("bad value for " + name);
  } else if (schema["type"] == "uinteger") {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    auto end = (char*)nullptr;
    value    = (uint64_t)strtoul(args[idx++].c_str(), &end, 10);
    if (end == nullptr) return cli_error("bad value for " + name);
  } else if (schema["type"] == "number") {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    auto end = (char*)nullptr;
    value    = strtod(args[idx++].c_str(), &end);
    if (end == nullptr) return cli_error("bad value for " + name);
  } else if (schema["type"] == "boolean") {
    if (positional) {
      if (idx >= args.size()) return cli_error("missing value for " + name);
      value = args[idx++] == "true" ? true : false;
    } else {
      value = true;
    }
  } else if (schema["type"] == "string") {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    value = args[idx++];
  } else if (schema["type"] == "array") {
    value = json_array{};
    if (idx + (size_t)schema["min_items"] >= args.size())
      return cli_error("missing value for " + name);
    auto end = std::min(idx + (size_t)schema["max_items"], args.size());
    while (idx < end) {
      if (!arg_to_json(value.emplace_back(), schema["items"], name, positional,
              args, idx, error))
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

  // init
  value = json_object{};

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
        if (!args_to_json(
                value[arg], schema["properties"][arg], args, idx, error))
          return false;
        break;
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

static bool validate_json(const json_value& value, const json_value& schema,
    const string& name, bool check_required, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  switch (value.type()) {
    case json_type::null: {
      return cli_error("bad value for " + name);
    } break;
    case json_type::integer: {
      if (schema["type"] != "number" && schema["type"] != "integer")
        return cli_error("bad value for " + name);
      return true;
    } break;
    case json_type::uinteger: {
      if (schema["type"] != "number" && schema["type"] != "integer")
        return cli_error("bad value for " + name);
      return true;
    } break;
    case json_type::number: {
      if (schema["type"] != "number" && schema["type"] != "integer")
        return cli_error("bad value for " + name);
      return true;
    } break;
    case json_type::boolean: {
      if (schema["type"] != "boolean")
        return cli_error("bad value for " + name);
      return true;
    } break;
    case json_type::string: {
      if (schema["type"] != "string") return cli_error("bad value for " + name);
      return true;
    } break;
    case json_type::array: {
      if (schema["type"] != "array") return cli_error("bad value for " + name);
      if ((size_t)schema["min_items"] > value.size())
        return cli_error("bad value for " + name);
      if ((size_t)schema["max_items"] < value.size())
        return cli_error("bad value for " + name);
      for (auto& item : value)
        if (!validate_json(item, schema["items"], name, false, error))
          return false;
      return true;
    } break;
    case json_type::object: {
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

  return true;
}

// update json objects
void update_value_objects(json_value& value, const json_value& update) {
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
static bool value_to_variable(
    const json_value& value, cli_variable& variable, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  if (variable.setter) {
    if (!variable.setter(value, variable.value, variable.choices)) return false;
  }

  for (auto& [key, property] : variable.variables) {
    if (!value.is_object()) throw std::runtime_error{"something went wrong"};
    if (value.contains(key)) {
      if (!value_to_variable(value[key], property, error)) {
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
    auto current_ = &value;
    while (true) {
      auto& current = *current_;
      if (!current.is_object()) break;
      if (current.contains("config") && current["config"].is_string())
        return (string)current["config"];
      if (current.contains("command") && current["command"].is_string()) {
        current_ = &current[(string)current["command"]];
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
  auto cli                    = cli_state{};
  cli.defaults                = json_object{};
  cli.schema                  = json_object{};
  cli.schema["title"]         = name;
  cli.schema["description"]   = usage;
  cli.schema["type"]          = "object";
  cli.schema["properties"]    = json_object{};
  cli.schema["required"]      = json_array{};
  cli.schema["clipositional"] = json_array{};
  return cli;
}

// add command
cli_command add_command(
    const cli_command& cli, const string& name, const string& usage) {
  auto& defaults            = get_defaults(cli);
  defaults[name]            = json_object{};
  auto& schema              = get_schema(cli);
  auto& property            = schema["properties"][name];
  property                  = json_object{};
  property["title"]         = name;
  property["description"]   = usage;
  property["type"]          = "object";
  property["properties"]    = json_object{};
  property["required"]      = json_array{};
  property["clipositional"] = json_array{};
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
    auto command  = string{};
    auto current_ = &value;
    while (true) {
      auto& current = *current_;
      if (!current.contains("command")) break;
      command += (command.empty() ? "" : "/") + (string)current["command"];
      current_ = &current[(string)current["command"]];
    }
    return command;
  };
  auto get_help = [&](const json_value& value) -> bool {
    auto current_ = &value;
    while (true) {
      auto& current = *current_;
      if (current.contains("help")) return (bool)current["help"];
      if (current.contains("command")) {
        current_ = &current[(string)current["command"]];
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
