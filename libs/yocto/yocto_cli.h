//
// # Yocto/CLI: Utilities for writing command-line apps
//
// Yocto/CLI is a collection of utilities used in writing command-line
// applications, including parsing command line arguments, printing values,
// timers and progress bars.
// Yocto/CLI is implemented in `yocto_cli.h`.
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
#include <string>
#include <type_traits>
#include <vector>

#include "yocto_math.h"

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
// Prints a message to the console and exit with an error. Returns error code.
int print_fatal(const string& message);

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
// COMMAND LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// Initialize a command line parser.
struct cli_state;
cli_state make_cli(const string& cmd, const string& usage);
// parse arguments, checks for errors, and exits on error or help
void parse_cli(cli_state& cli, int argc, const char** argv);
void parse_cli(cli_state& cli, const vector<string>& args);
// parse arguments, checks for errors
bool parse_cli(cli_state& cli, const vector<string>& args, string& error);

// Add a subcommand
struct cli_command;
cli_command add_command(
    const cli_command& cli, const string& name, const string& usage);
void                set_command_var(const cli_command& cli, string& value);
void                set_help_var(const cli_command& cli, bool& value);
[[deprecated]] void add_command_name(const cli_command& cli, const string& name,
    string& value, const string& usage);

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

// Simple ordered map
template <typename Key, typename Value>
struct ordered_map {
  size_t size() const { return _data.size(); }
  bool   empty() const { return _data.empty(); }

  Value& operator[](const Key& key) {
    if (auto it = _search(key); it) return it->second;
    _data.emplace_back(key, Value{});
    return _data.back().second;
  }
  Value& at(const Key& key) {
    if (auto it = _search(key); it) return it->second;
    throw std::out_of_range{"missing key for " + key};
  }
  const Value& at(const Key& key) const {
    if (auto it = _search(key); it) return it->second;
    throw std::out_of_range{"missing key for " + key};
  }

  pair<Key, Value>* find(const string& key) {
    if (auto it = _search(key); it) return it;
    return end();
  }
  const pair<Key, Value>* find(const string& key) const {
    if (auto it = _search(key); it) return it;
    return end();
  }

  pair<Key, Value>*       begin() { return _data.data(); }
  pair<Key, Value>*       end() { return _data.data() + _data.size(); }
  const pair<Key, Value>* begin() const { return _data.data(); }
  const pair<Key, Value>* end() const { return _data.data() + _data.size(); }

 private:
  vector<pair<Key, Value>> _data;
  pair<Key, Value>*        _search(const string& key) {
    for (auto& item : _data)
      if (key == item.first) return &item;
    return nullptr;
  }
  const pair<Key, Value>* _search(const string& key) const {
    for (auto& item : _data)
      if (key == item.first) return &item;
    return nullptr;
  }
};

// Command line type.
enum struct cli_type {
  none,
  integer,
  unsigned_,
  number,
  boolean,
  string,
  array,
  object
};
// Command line value.
struct cli_value;
using cli_array  = vector<cli_value>;
using cli_object = ordered_map<string, cli_value>;
struct cli_value {
  cli_type   type      = cli_type::none;
  int64_t    integer   = 0;
  uint64_t   unsigned_ = 0;
  double     number    = 0;
  bool       boolean   = false;
  string     string_   = {};
  cli_array  array     = {};
  cli_object object    = {};

  cli_value() : type{cli_type::none} {}
  cli_value(cli_type type) : type{type} {}
  cli_value(const cli_value& other) {
    type = other.type;
    switch (type) {
      case cli_type::none: break;
      case cli_type::integer: integer = other.integer; break;
      case cli_type::unsigned_: unsigned_ = other.unsigned_; break;
      case cli_type::number: number = other.number; break;
      case cli_type::boolean: boolean = other.boolean; break;
      case cli_type::string: string_ = other.string_; break;
      case cli_type::array: array = other.array; break;
      case cli_type::object: object = other.object; break;
    }
  }
  cli_value(cli_value&& value) : cli_value() { _swap(value); }
  cli_value& operator=(cli_value other) {
    _swap(other);
    return *this;
  }

  cli_type get_type() const { return type; }
  void     set_type(cli_type type) {
    if (this->type == type) return;
    this->type = type;
  }

  void set_null() { _set(cli_type::none, integer, (int64_t)0); }
  void set_integer(int64_t value) { _set(cli_type::integer, integer, value); }
  void set_unsigned(uint64_t value) {
    _set(cli_type::unsigned_, unsigned_, value);
  }
  void set_number(double value) { _set(cli_type::number, number, value); }
  void set_boolean(bool value) { _set(cli_type::boolean, boolean, value); }
  void set_string(const char* value) {
    _set(cli_type::string, string_, string{value});
  }
  void set_string(const string& value) {
    _set(cli_type::string, string_, value);
  }
  void set_array() { _set(cli_type::array, array, {}); }
  void set_array(size_t size) { _set(cli_type::array, array, cli_array(size)); }
  void set_array(const cli_array& value) {
    _set(cli_type::array, array, value);
  }
  void set_object() { _set(cli_type::object, object, {}); }
  void set_object(const cli_object& value) {
    return _set(cli_type::object, object, value);
  }

  int64_t&    get_integer() { return _get(cli_type::integer, integer); }
  uint64_t&   get_unsigned() { return _get(cli_type::unsigned_, unsigned_); }
  double&     get_number() { return _get(cli_type::number, number); }
  bool&       get_boolean() { return _get(cli_type::boolean, boolean); }
  string&     get_string() { return _get(cli_type::string, string_); }
  cli_array&  get_array() { return _get(cli_type::array, array); }
  cli_object& get_object() { return _get(cli_type::object, object); }

  const int64_t& get_integer() const {
    return _get(cli_type::integer, integer);
  }
  const uint64_t& get_unsigned() const {
    return _get(cli_type::unsigned_, unsigned_);
  }
  const double& get_number() const { return _get(cli_type::number, number); }
  const bool&   get_boolean() const { return _get(cli_type::boolean, boolean); }
  const string& get_string() const { return _get(cli_type::string, string_); }
  const cli_array&  get_array() const { return _get(cli_type::array, array); }
  const cli_object& get_object() const {
    return _get(cli_type::object, object);
  }

  void set(std::nullptr_t) { set_null(); }
  void set(int32_t value) { set_integer(value); }
  void set(int64_t value) { set_integer(value); }
  void set(uint32_t value) { set_unsigned(value); }
  void set(uint64_t value) { set_unsigned(value); }
  void set(bool value) { set_boolean(value); }
  void set(const string& value) { set_string(value); }
  void set(const char* value) { set_string(value); }
  template <typename T, size_t N>
  void set(const std::array<T, N>& value) {
    set_array(value.size());
    for (auto idx = (size_t)0; idx < value.size(); idx++)
      get_array().at(idx).set(value.at(idx));
  }
  template <typename T>
  void set(const vector<T>& value) {
    set_array(value.size());
    for (auto idx = (size_t)0; idx < value.size(); idx++)
      get_array().at(idx).set(value.at(idx));
  }

 private:
  void _swap(cli_value& other) {
    std::swap(type, other.type);
    switch (type) {
      case cli_type::none: break;
      case cli_type::integer: std::swap(integer, other.integer); break;
      case cli_type::unsigned_: std::swap(unsigned_, other.unsigned_); break;
      case cli_type::number: std::swap(number, other.number); break;
      case cli_type::boolean: std::swap(boolean, other.boolean); break;
      case cli_type::string: std::swap(string_, other.string_); break;
      case cli_type::array: std::swap(array, other.array); break;
      case cli_type::object: std::swap(object, other.object); break;
    }
  }
  template <typename T>
  const T& _get(cli_type type, const T& value) const {
    if (this->type != type) throw std::invalid_argument{"bad json type"};
    return value;
  }
  template <typename T>
  T& _get(cli_type type, T& value) {
    if (this->type != type) throw std::invalid_argument{"bad json type"};
    return value;
  }
  template <typename T>
  void _set(cli_type type, T& value, const T& other) {
    if (this->type != type) set_type(type);
    value = other;
  }
};

// Command line schema.
struct cli_schema {
  cli_type           type            = cli_type::object;
  string             title           = "";
  string             description     = "";
  cli_value          default_        = {};
  cli_value          min             = {};
  cli_value          max             = {};
  vector<string>     enum_           = {};
  size_t             min_items       = 0;
  size_t             max_items       = std::numeric_limits<size_t>::max();
  vector<string>     required        = {};
  vector<string>     cli_positionals = {};
  string             cli_config      = "";
  vector<cli_schema> items           = {};
  ordered_map<string, cli_schema> properties = {};
};

// Command line setter.
using cli_setter = bool (*)(
    const cli_value&, void*, const vector<string>& choices);
// Command line variable.
struct cli_variable {
  void*                             value     = nullptr;
  cli_setter                        setter    = nullptr;
  vector<string>                    choices   = {};
  ordered_map<string, cli_variable> variables = {};
};
// Command line state.
struct cli_state {
  cli_value    defaults  = {};
  cli_schema   schema    = {};
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

#endif
