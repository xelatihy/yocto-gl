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
// IMPLEMENTATION OF COMMAND-LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto {

inline size_t json_empty(const json_value& json) {
  if (json.get_type() == json_type::array) {
    return json.get_array().empty();
  } else if (json.get_type() == json_type::object) {
    return json.get_object().empty();
  } else {
    throw json_error{"array or object expected"};
  }
}
inline size_t json_size(const json_value& json) {
  if (json.get_type() == json_type::array) {
    return json.get_array().size();
  } else if (json.get_type() == json_type::object) {
    return json.get_object().size();
  } else {
    throw json_error{"array or object expected"};
  }
}
inline const json_value& json_at(const json_value& json, size_t idx) {
  if (json.get_type() == json_type::array) {
    if (idx >= json.get_array().size()) throw json_error{"out of bounds"};
    return json.get_array().at(idx);
  } else {
    throw json_error{"array expected"};
  }
}
inline json_value& json_at(json_value& json, size_t idx) {
  if (json.get_type() == json_type::array) {
    if (idx >= json.get_array().size()) throw json_error{"out of bounds"};
    return json.get_array().at(idx);
  } else {
    throw json_error{"array expected"};
  }
}
inline const json_value& json_at(const json_value& json, const string& key) {
  if (json.get_type() == json_type::array) {
    auto it = json.get_object().find(key);
    if (it == json.get_object().end()) throw json_error{"missing key " + key};
    return it->second;
  } else {
    throw json_error{"object expected"};
  }
}
inline json_value& json_at(json_value& json, const string& key) {
  if (json.get_type() == json_type::array) {
    auto it = json.get_object().find(key);
    if (it == json.get_object().end()) throw json_error{"missing key " + key};
    return it->second;
  } else {
    throw json_error{"object expected"};
  }
}
inline json_value& json_append(json_value& json) {
  if (json.get_type() == json_type::array) {
    return json.get_array().emplace_back();
  } else {
    throw json_error{"array expected"};
  }
}
inline json_value& json_insert(json_value& json, const string& key) {
  if (json.get_type() == json_type::object) {
    return json.get_object()[key];
  } else {
    throw json_error{"object expected"};
  }
}

template <typename T>
inline json_value to_json(const T& value) {
  auto json = json_value{};
  to_json(json, value);
  return json;
}

inline void to_json(json_value& json, int64_t value) {
  json.set_integer(value);
}
inline void to_json(json_value& json, int32_t value) {
  json.set_integer(value);
}
inline void to_json(json_value& json, uint64_t value) {
  json.set_unsigned(value);
}
inline void to_json(json_value& json, uint32_t value) {
  json.set_unsigned(value);
}
inline void to_json(json_value& json, double value) { json.set_number(value); }
inline void to_json(json_value& json, float value) { json.set_number(value); }
inline void to_json(json_value& json, bool value) { json.set_boolean(value); }
inline void to_json(json_value& json, const char* value) {
  json.set_string(value);
}
inline void to_json(json_value& json, const string& value) {
  json.set_string(value);
}
inline void to_json(json_value& json, const json_array& value) {
  json.set_array(value);
}
inline void to_json(json_value& json, const json_object& value) {
  json.set_object(value);
}
template <typename T, size_t N>
inline void to_json(json_value& json, const array<T, N>& value) {
  json.set_array(value.size());
  for (auto idx = 0; idx < value.size(); idx++)
    to_json(json.get_array().at(idx), value.at(idx));
}
template <typename T>
inline void to_json(json_value& json, const vector<T>& value) {
  json.set_array(value.size());
  for (auto idx = 0; idx < value.size(); idx++)
    to_json(json.get_array().at(idx), value.at(idx));
}

inline void to_json_array(json_value& json) { json.set_array(); }
inline void to_json_array(json_value& json, size_t size) {
  json.set_array(size);
}
inline json_value& to_json_append(json_value& json) {
  if (json.type != json_type::array) json.set_array();
  return json.get_array().emplace_back();
}
inline void        to_json_object(json_value& json) { json.set_object(); }
inline json_value& to_json_insert(json_value& json, const string& key) {
  if (json.type != json_type::object) json.set_object();
  return json.get_object()[key];
}

template <typename T>
inline T from_json(const json_value& json) {
  auto value = T{};
  from_json(json, value);
  return value;
}

inline void from_json(const json_value& json, int64_t& value) {
  if (json.get_type() == json_type::integer) {
    value = (int64_t)json.get_integer();
  } else if (json.get_type() == json_type::unsigned_) {
    value = (int64_t)json.get_unsigned();
  } else {
    throw json_error{"integer expected"};
  }
}
inline void from_json(const json_value& json, int32_t& value) {
  value = (int32_t)from_json<int64_t>(json);
}
inline void from_json(const json_value& json, uint64_t& value) {
  if (json.get_type() == json_type::integer) {
    value = (uint64_t)json.get_integer();
  } else if (json.get_type() == json_type::unsigned_) {
    value = (uint64_t)json.get_unsigned();
  } else {
    throw json_error{"integer expected"};
  }
}
inline void from_json(const json_value& json, uint32_t& value) {
  value = (uint32_t)from_json<uint64_t>(json);
}
inline void from_json(const json_value& json, double& value) {
  if (json.get_type() == json_type::integer) {
    value = (double)json.get_integer();
  } else if (json.get_type() == json_type::unsigned_) {
    value = (double)json.get_unsigned();
  } else if (json.get_type() == json_type::number) {
    value = (double)json.get_number();
  } else {
    throw json_error{"number expected"};
  }
}
inline void from_json(const json_value& json, float& value) {
  value = (double)from_json<double>(json);
}
inline void from_json(const json_value& json, bool& value) {
  if (json.get_type() == json_type::boolean) {
    value = json.get_boolean();
  } else {
    throw json_error{"boolean expected"};
  }
}
inline void from_json(const json_value& json, string& value) {
  if (json.get_type() == json_type::string) {
    value = json.get_string();
  } else {
    throw json_error{"string expected"};
  }
}
template <typename T, size_t N>
inline void from_json(const json_value& json, array<T, N>& value) {
  if (json.get_type() == json_type::array) {
    if (json.get_array().size() != N)
      throw json_error{"array of size " + std::to_string(N) + " expected"};
    for (auto idx = (size_t)0; idx < value.size(); idx++)
      from_json(json.get_array().at(idx), value.at(idx));
  } else {
    throw json_error{"array expected"};
  }
}
template <typename T>
inline void from_json(const json_value& json, vector<T>& value) {
  if (json.get_type() == json_type::array) {
    value.resize(json.get_array().size());
    for (auto idx = (size_t)0; idx < value.size(); idx++)
      from_json(json.get_array().at(idx), value.at(idx));
  } else {
    throw json_error{"array expected"};
  }
}

inline void to_schema(json_schema& schema, int64_t value, const string& title,
    const string& description) {
  schema.type        = to_json(value).type;
  schema.title       = title;
  schema.description = description;
  schema.default_    = to_json(value);
}
inline void to_schema(json_schema& schema, int64_t value, const string& title,
    const string& description, int64_t min, int64_t max) {
  schema.type        = to_json(value).type;
  schema.title       = title;
  schema.description = description;
  schema.default_    = to_json(value);
  schema.min         = to_json(min);
  schema.max         = to_json(max);
}
inline void to_schema(json_schema& schema, int32_t value, const string& title,
    const string& description) {
  to_schema(schema, (int64_t)value, title, description);
}
inline void to_schema(json_schema& schema, int32_t value, const string& title,
    const string& description, int32_t min, int32_t max) {
  to_schema(
      schema, (int64_t)value, title, description, (int64_t)min, (int64_t)max);
}
inline void to_schema(json_schema& schema, uint64_t value, const string& title,
    const string& description) {
  schema.type        = to_json(value).type;
  schema.title       = title;
  schema.description = description;
  schema.default_    = to_json(value);
}
inline void to_schema(json_schema& schema, uint64_t value, const string& title,
    const string& description, uint64_t min, uint64_t max) {
  schema.type        = to_json(value).type;
  schema.title       = title;
  schema.description = description;
  schema.default_    = to_json(value);
  schema.min         = to_json(min);
  schema.max         = to_json(max);
}
inline void to_schema(json_schema& schema, uint32_t value, const string& title,
    const string& description) {
  to_schema(schema, (uint64_t)value, title, description);
}
inline void to_schema(json_schema& schema, uint32_t value, const string& title,
    const string& description, uint32_t min, uint32_t max) {
  to_schema(schema, (uint64_t)value, title, description, (uint64_t)min,
      (uint64_t)max);
}
inline void to_schema(json_schema& schema, double value, const string& title,
    const string& description) {
  schema.type        = to_json(value).type;
  schema.title       = title;
  schema.description = description;
  schema.default_    = to_json(value);
}
inline void to_schema(json_schema& schema, double value, const string& title,
    const string& description, double min, double max) {
  schema.type        = to_json(value).type;
  schema.title       = title;
  schema.description = description;
  schema.default_    = to_json(value);
  schema.min         = to_json(min);
  schema.max         = to_json(max);
}
inline void to_schema(json_schema& schema, float value, const string& title,
    const string& description) {
  to_schema(schema, (double)value, title, description);
}
inline void to_schema(json_schema& schema, float value, const string& title,
    const string& description, float min, float max) {
  to_schema(
      schema, (double)value, title, description, (double)min, (double)max);
}
inline void to_schema(json_schema& schema, bool value, const string& title,
    const string& description) {
  schema.type        = to_json(value).type;
  schema.title       = title;
  schema.description = description;
  schema.default_    = to_json(value);
}
inline void to_schema(json_schema& schema, const string& value,
    const string& title, const string& description) {
  schema.type        = to_json(value).type;
  schema.title       = title;
  schema.description = description;
  schema.default_    = to_json(value);
}
inline void to_schema(json_schema& schema, const string& value,
    const string& title, const string& description,
    const vector<string>& enum_) {
  schema.type        = to_json(value).type;
  schema.title       = title;
  schema.description = description;
  schema.default_    = to_json(value);
  schema.enum_       = enum_;
}
template <typename T, size_t N>
inline void to_schema(json_schema& schema, const array<T, N>& value,
    const string& title, const string& description) {
  schema.type        = to_json(value).type;
  schema.title       = title;
  schema.description = description;
  schema.default_    = to_json(value);
  schema.min_items   = N;
  schema.max_items   = N;
  to_schema(schema.items.emplace_back(), T{}, "item", "");
}
template <typename T>
inline void to_schema(json_schema& schema, const vector<T>& value,
    const string& title, const string& description) {
  schema.type        = to_json(value).type;
  schema.title       = title;
  schema.description = description;
  schema.default_    = to_json(value);
  schema.min_items   = 0;
  schema.max_items   = std::numeric_limits<size_t>::max();
  to_schema(schema.items.emplace_back(), T{}, "item", "");
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
static void cli_to_json(
    json_value& json, const T& value, const vector<string>& choices) {
  if (choices.empty()) {
    to_json(json, value);
  } else {
    if constexpr (cli_is_array_v<T>) {
      to_json_array(json, value.size());
      for (auto idx = (size_t)0; idx < value.size(); idx++) {
        cli_to_json(json_at(json, idx), value.at(idx), choices);
      }
    } else if constexpr (cli_is_vector_v<T>) {
      to_json_array(json, value.size());
      for (auto idx = (size_t)0; idx < value.size(); idx++) {
        cli_to_json(json_at(json, idx), value.at(idx), choices);
      }
    } else {
      if constexpr (std::is_integral_v<T> || std::is_floating_point_v<T>) {
        to_json(json, choices.at((uint64_t)value));
      } else if constexpr (std::is_same_v<T, string>) {
        if (std::find(choices.begin(), choices.end(), value) == choices.end())
          throw std::out_of_range{"bad value"};
        to_json(json, value);
      } else {
        throw std::invalid_argument{"not supported"};
      }
    }
  }
}

template <typename T>
static void cli_to_schema(json_schema& schema, const T& value,
    const vector<string>& choices, const string& name, const string& usage) {
  schema.title       = name;
  schema.description = usage;
  cli_to_json(schema.default_, value, choices);
  if constexpr (std::is_same_v<T, int32_t> || std::is_same_v<T, int64_t>) {
    if (!choices.empty()) {
      schema.type  = json_type::string;
      schema.enum_ = choices;
    } else {
      schema.type = json_type::integer;
    }
  } else if constexpr (std::is_same_v<T, uint32_t> ||
                       std::is_same_v<T, uint64_t>) {
    if (!choices.empty()) {
      schema.type  = json_type::string;
      schema.enum_ = choices;
    } else {
      schema.type = json_type::unsigned_;
    }
  } else if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>) {
    if (!choices.empty()) {
      schema.type  = json_type::string;
      schema.enum_ = choices;
    } else {
      schema.type = json_type::number;
    }
  } else if constexpr (std::is_same_v<T, bool>) {
    if (!choices.empty()) {
      schema.type  = json_type::string;
      schema.enum_ = choices;
    } else {
      schema.type = json_type::boolean;
    }
  } else if constexpr (std::is_same_v<T, string>) {
    if (!choices.empty()) {
      schema.type  = json_type::string;
      schema.enum_ = choices;
    } else {
      schema.type = json_type::string;
    }
  } else if constexpr (cli_is_array_v<T>) {
    schema.type      = json_type::array;
    schema.min_items = value.size();
    schema.max_items = value.size();
    cli_to_schema(schema.items.emplace_back(), T{}, choices, "item", "");
  } else if constexpr (cli_is_vector_v<T>) {
    schema.type      = json_type::array;
    schema.min_items = 0;
    schema.max_items = std::numeric_limits<size_t>::max();
    cli_to_schema(schema.items.emplace_back(), T{}, choices, "item", "");
  }
}

template <typename T>
static void cli_from_json(
    const json_value& json, T& value, const vector<string>& choices) {
  if (choices.empty()) {
    from_json(json, value);
  } else {
    if constexpr (cli_is_array_v<T>) {
      auto size = json_size(json);
      if (size != value.size()) throw json_error{"bad array size"};
      for (auto idx = (size_t)0; idx < value.size(); idx++) {
        cli_from_json(json_at(json, idx), value.at(idx), choices);
      }
    } else if constexpr (cli_is_vector_v<T>) {
      auto size = json_size(json);
      value.resize(size);
      for (auto idx = (size_t)0; idx < value.size(); idx++) {
        cli_from_json(json_at(json, idx), value.at(idx), choices);
      }
    } else {
      auto values = from_json<string>(json);
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
    return cli.state->defaults.object.at(cli.path);
}
static json_schema& get_schema(const cli_command& cli) {
  if (cli.path.empty())
    return cli.state->schema;
  else
    return cli.state->schema.properties.at(cli.path);
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
  cli_to_json(defaults.object[name], value, choices);
  cli_to_schema(schema.properties[name], value, choices, name, usage);
  if (req) schema.required.push_back(name);
  if (positional) schema.cli_positionals.push_back(name);
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
  cli_to_json(defaults.object[name], value, choices);
  cli_to_schema(schema.properties[name], value, choices, name, usage);
  if (req) schema.required.push_back(name);
  if (positional) schema.cli_positionals.push_back(name);
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
  cli_to_json(defaults.object[name], value, choices);
  cli_to_schema(schema.properties[name], value, choices, name, usage);
  if (req) schema.required.push_back(name);
  if (positional) schema.cli_positionals.push_back(name);
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
    get_schema(cli).properties.at(name).cli_config = config;
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
    get_schema(cli).properties.at(name).cli_config = config;
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
    const json_schema& schema_, const string& command, const string& program) {
  auto type_to_string = [](const json_schema& schema) -> string {
    auto type = schema.type;
    switch (type) {
      case json_type::none: return "";
      case json_type::integer: return "<integer>";
      case json_type::unsigned_: return "<integer>";
      case json_type::number: return "<number>";
      case json_type::boolean: return "<bool>";
      case json_type::string: return "<string>";
      case json_type::array: return "<array>";
      case json_type::object: return "<object>";
    }
  };
  auto default_to_string = [](const json_schema& schema) -> string {
    auto& value = schema.default_;
    switch (value.type) {
      case json_type::none: return "";
      case json_type::integer: return "[" + std::to_string(value.integer) + "]";
      case json_type::unsigned_:
        return "[" + std::to_string(value.unsigned_) + "]";
      case json_type::number: return "[" + std::to_string(value.number) + "]";
      case json_type::boolean:
        return "[" + (value.boolean ? string{"true"} : string{"false"}) + "]";
      case json_type::string: return "[" + value.string_ + "]";
      case json_type::array: return "[]";
      case json_type::object: return "";
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
    if (property.type == json_type::object) {
      auto line = "  " + key;
      while (line.size() < 32) line += " ";
      line += property.description + "\n";
      usage_commands += line;
    } else {
      auto is_positional = std::find(schema.cli_positionals.begin(),
                               schema.cli_positionals.end(),
                               key) != schema.cli_positionals.end();
      auto line          = (is_positional ? "  " : "  --") + key;
      if (property.type != json_type::boolean || is_positional) {
        line += " " + type_to_string(property);
      }
      while (line.size() < 32) line += " ";
      line += property.description;
      if (property.default_.type != json_type::none) {
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

static bool arg_to_value(json_value& value, const json_schema& schema,
    const string& name, bool positional, const vector<string>& args,
    size_t& idx, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  if (!schema.enum_.empty()) {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    value.type    = json_type::string;
    value.string_ = args[idx++];
  } else if (schema.type == json_type::integer) {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    auto end      = (char*)nullptr;
    value.type    = json_type::integer;
    value.integer = strtol(args[idx++].c_str(), &end, 10);
    if (end == nullptr) return cli_error("bad value for " + name);
  } else if (schema.type == json_type::unsigned_) {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    auto end      = (char*)nullptr;
    value.type    = json_type::unsigned_;
    value.integer = strtoul(args[idx++].c_str(), &end, 10);
    if (end == nullptr) return cli_error("bad value for " + name);
  } else if (schema.type == json_type::number) {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    auto end     = (char*)nullptr;
    value.type   = json_type::number;
    value.number = strtod(args[idx++].c_str(), &end);
    if (end == nullptr) return cli_error("bad value for " + name);
  } else if (schema.type == json_type::boolean) {
    if (positional) {
      if (idx >= args.size()) return cli_error("missing value for " + name);
      value.type    = json_type::string;
      value.string_ = args[idx++] == "true" ? true : false;
    } else {
      value.type    = json_type::boolean;
      value.boolean = true;
    }
  } else if (schema.type == json_type::string) {
    if (idx >= args.size()) return cli_error("missing value for " + name);
    value.type    = json_type::string;
    value.string_ = args[idx++];
  } else if (schema.type == json_type::array) {
    value.type = json_type::array;
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

static bool args_to_value(json_value& value, const json_schema& schema,
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
  value.type   = json_type::object;
  value.object = {};

  // add things to schema
  auto commands = vector<string>{}, positionals = vector<string>{};
  for (auto& [key, value] : schema.properties) {
    if (value.type == json_type::object) commands.push_back(key);
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
      value.object["help"].type    = json_type::boolean;
      value.object["help"].boolean = true;
      continue;
    }
    if (arg == "--config") {
      if (idx >= args.size()) return cli_error("missing value for config");
      value.object["config"].type    = json_type::string;
      value.object["config"].string_ = args[idx++];
      continue;
    }
    auto is_positional = arg.find('-') != 0;
    if (!commands.empty() && is_positional) {
      if (std::find(commands.begin(), commands.end(), arg) != commands.end()) {
        value.object["command"].type    = json_type::string;
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
        value.object["config"].type = json_type::string;
        value.object["config"].string_ =
            "try:" +
            get_try_config(value.object.at(name).string_, oschema.cli_config);
      }
    } else {
      auto name = string{};
      for (auto& [key, value] : schema.properties) {
        if ("--" + key == arg && value.type != json_type::object) {
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
        value.object["config"].type = json_type::string;
        value.object["config"].string_ =
            "try:" +
            get_try_config(value.object.at(name).string_, oschema.cli_config);
      }
    }
  }

  // done
  return true;
}

static bool validate_value(const json_value& value, const json_schema& schema,
    const string& name, bool check_required, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  auto contains = [](const auto& object, const string& name) {
    return object.find(name) != object.end();
  };

  switch (value.type) {
    case json_type::none: {
      return cli_error("bad value for " + name);
    } break;
    case json_type::integer: {
      if (schema.type != json_type::integer &&
          schema.type != json_type::unsigned_ &&
          schema.type != json_type::number)
        return cli_error("bad value for " + name);
      return true;
    } break;
    case json_type::unsigned_: {
      if (schema.type != json_type::integer &&
          schema.type != json_type::unsigned_ &&
          schema.type != json_type::number)
        return cli_error("bad value for " + name);
      return true;
    } break;
    case json_type::number: {
      if (schema.type != json_type::number)
        return cli_error("bad value for " + name);
      return true;
    } break;
    case json_type::boolean: {
      if (schema.type != json_type::boolean)
        return cli_error("bad value for " + name);
      return true;
    } break;
    case json_type::string: {
      if (schema.type != json_type::string)
        return cli_error("bad value for " + name);
      return true;
    } break;
    case json_type::array: {
      if (schema.type != json_type::array)
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
    case json_type::object: {
      if (schema.type != json_type::object)
        return cli_error("bad value for " + name);
      for (auto& [key, property] : value.object) {
        if (key == "help") {
          if (property.type != json_type::boolean)
            return cli_error("bad value for " + key);
        } else if (key == "config") {
          if (property.type != json_type::string)
            return cli_error("bad value for " + key);
        } else if (key == "command") {
          if (property.type != json_type::string)
            return cli_error("bad value for " + key);
        } else {
          if (!contains(schema.properties, key))
            return cli_error("unknown option " + key);
          auto selected_command = contains(value.object, "command") &&
                                  value.object.at("command").type ==
                                      json_type::string &&
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
void update_value_objects(json_value& value, const json_value& update) {
  auto contains = [](const auto& object, const string& name) {
    return object.find(name) != object.end();
  };

  if (value.type != json_type::object) return;
  for (auto& [key, property] : update.object) {
    if (property.type == json_type::object && contains(value.object, key) &&
        value.object.at(key).type == json_type::object) {
      update_value_objects(value.object.at(key), value);
    } else {
      value.object[key] = property;
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
    if (value.type != json_type::object)
      throw std::runtime_error{"something went wrong"};
    if (contains(value.object, key)) {
      if (!value_to_variable(value.object.at(key), property, error)) {
        return cli_error("bad value for " + key);
      }
    }
  }

  return true;
}

using ordered_json = nlohmann::ordered_json;
void from_json(const ordered_json& js, json_value& value) {
  switch (js.type()) {
    case ordered_json::value_t::null: {
      value.type = json_type::none;
    } break;
    case ordered_json::value_t::number_integer: {
      value.type    = json_type::integer;
      value.integer = (int64_t)js;
    } break;
    case ordered_json::value_t::number_unsigned: {
      value.type      = json_type::unsigned_;
      value.unsigned_ = (uint64_t)js;
    } break;
    case ordered_json::value_t::number_float: {
      value.type   = json_type::number;
      value.number = (double)js;
    } break;
    case ordered_json::value_t::boolean: {
      value.type    = json_type::boolean;
      value.boolean = (bool)js;
    } break;
    case ordered_json::value_t::string: {
      value.type    = json_type::string;
      value.string_ = (string)js;
    } break;
    case ordered_json::value_t::array: {
      value.type = json_type::array;
      for (auto& jitem : js) from_json(jitem, value.array.emplace_back());
    } break;
    case ordered_json::value_t::object: {
      value.type = json_type::object;
      for (auto& [key, jitem] : js.items()) from_json(jitem, value.object[key]);
    } break;
    case ordered_json::value_t::binary: {
      value.type = json_type::none;
    } break;
    case ordered_json::value_t::discarded: {
      value.type = json_type::none;
    } break;
  }
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
      if (current->type != json_type::object) break;
      if (contains(current->object, "config") &&
          current->object.at("config").type == json_type::string)
        return current->object.at("config").string_;
      if (contains(current->object, "command") &&
          current->object.at("command").type == json_type::string) {
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

  auto js = nlohmann::ordered_json{};
  auto fs = std::ifstream(config);
  if (!fs) {
    if (try_config) return true;
    return cli_error("missing configuration file " + config);
  }
  try {
    js = nlohmann::ordered_json::parse(fs);
  } catch (...) {
    return cli_error("error loading configuration " + config);
  }

  auto cvalue = json_value{};
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
  cli.defaults.type      = json_type::object;
  cli.defaults.object    = {};
  cli.schema.title       = name;
  cli.schema.description = usage;
  cli.schema.type        = json_type::object;
  cli.schema.properties  = {};
  return cli;
}

// add command
cli_command add_command(
    const cli_command& cli, const string& name, const string& usage) {
  auto& defaults                      = get_defaults(cli);
  defaults.object[name].type          = json_type::object;
  auto& schema                        = get_schema(cli);
  schema.properties[name].title       = name;
  schema.properties[name].description = usage;
  schema.properties[name].type        = json_type::object;
  auto& variables                     = get_variables(cli);
  variables.variables[name]           = {};
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
  if (!args_to_value(value, get_schema(cli), args, idx, error)) return false;
  if (!config_to_value(value, error)) return false;
  if (!validate_value(value, get_schema(cli), "", true, error)) return false;
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
      if (!contains(current->object, "command")) break;
      command += (command.empty() ? "" : "/") +
                 current->object.at("command").string_;
      current = &current->object.at(current->object.at("command").string_);
    }
    return command;
  };
  auto get_help = [&](const json_value& value) -> bool {
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
