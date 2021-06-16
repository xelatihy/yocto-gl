//
// # Yocto/CLI: Utilities for writing command-line apps
//
// Yocto/CLI is a collection of utilities used in writing command-line
// applications, including parsing command line arguments, printing values,
// timers, progress bars, Json data, Json IO, file IO, path manipulation.
// Yocto/CLI is implemented in `yocto_cli.h` and `yocto_cli.cpp`, and
// depends on `json.hpp` for Json serialization.
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
// ORDERED MAP
// -----------------------------------------------------------------------------
namespace yocto {

// Simple ordered map
template <typename Key, typename Value>
struct ordered_map {
  size_t size() const { return _data.size(); }
  bool   empty() const { return _data.empty(); }

  bool contains(const Key& key) const { return (bool)_search(key); }

  Value& operator[](const Key& key) {
    if (auto it = _search(key); it) return it->second;
    return _data.emplace_back(key, Value{}).second;
  }
  const Value& operator[](const Key& key) const {
    if (auto it = _search(key); it) return it->second;
    throw std::out_of_range{"missing key for " + key};
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

  pair<Key, Value>* _search(const string& key) {
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

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON DATA
// -----------------------------------------------------------------------------
namespace yocto {

// Json forward declaration
struct json_value;

// Json error
struct json_error : std::logic_error {
  const json_value* _where = nullptr;
  json_error(const string& error, const json_value* where)
      : std::logic_error(error), _where(where) {}
};

// Json enum support
template <typename T>
struct json_enum;

// Json type
enum struct json_type {
  // clang-format off
  null, integer, uinteger, number, boolean, string, array, object
  // clang-format on
};

// Json typdefs
using json_array  = vector<json_value>;
using json_object = ordered_map<string, json_value>;

// Json value
struct json_value {
  // Json value
  json_value() : _type{json_type::null}, _integer{0} {}
  explicit json_value(json_type type) : _type{type} { _init(); }
  json_value(const json_value& other) : json_value() { _copy(other); }
  json_value(json_value&& value) : json_value() { _swap(value); }
  json_value& operator=(const json_value& other) {
    auto copy = json_value{other};
    return _swap(copy);
  }
  json_value& operator=(json_value&& other) { return _swap(other); }
  ~json_value() { _clear(); }

  // conversions, casts, assignments
  template <typename T>
  explicit json_value(const T& value) : json_value() {
    set(value);
  }
  template <typename T>
  explicit operator T() const {
    auto value = T{};
    get(value);
    return value;
  }
  template <typename T>
  json_value& operator=(const T& value) {
    set(value);
    return *this;
  }

  // type
  // clang-format off
  json_type type() const { return _type; }
  bool is_null() const { return _type == json_type::null; }
  bool is_integer() const { return _type == json_type::integer || _type == json_type::uinteger; }
  bool is_number() const { return _type == json_type::integer || _type == json_type::uinteger || _type == json_type::number; }
  bool is_boolean() const { return _type == json_type::boolean; }
  bool is_string() const { return _type == json_type::string; }
  bool is_array() const { return _type == json_type::array; }
  bool is_object() const { return _type == json_type::object; }
  // clang-format on

  // size
  bool   empty() const { return _empty(); }
  size_t size() const { return _size(); }

  // object creation
  static json_value array() { return json_value{json_array{}}; }
  static json_value array(size_t size) { return json_value{json_array(size)}; }
  static json_value object() { return json_value{json_object{}}; }

  // element access
  // clang-format off
  json_value&       operator[](size_t idx) { return _get_array().at(idx); }
  const json_value& operator[](size_t idx) const { return _get_array().at(idx); }
  json_value&       operator[](const string& key) { return _get_object()[key]; }
  const json_value& operator[](const string& key) const { return _get_object().at(key); }
  json_value&       at(size_t idx) { return _get_array().at(idx); }
  const json_value& at(size_t idx) const { return _get_array().at(idx); }
  json_value&       at(const string& key) { return _get_object().at(key); }
  const json_value& at(const string& key) const { return _get_object().at(key); }
  bool contains(const string& key) const { return _get_object().contains(key); }
  json_value& front() { return _get_array().front(); }
  const json_value& front() const { return _get_array().front(); }
  json_value& back() { return _get_array().back(); }
  const json_value& back() const { return _get_array().back(); }
  json_value& append() { return _get_array().emplace_back(); }
  template<typename T>
  void append(const T& value) { _get_array().push_back(json_value{value}); }
  json_value& emplace_back() { return _get_array().emplace_back(); }
  template<typename ... Args>
  json_value& emplace_back(Args&&... args) { return _get_array().emplace_back(std::forward<Args>(args)...); }
  void push_back(const json_value& item) { _get_array().push_back(item); }
  template<typename T>
  void push_back(const T& value) { _get_array().push_back(json_value{value}); }
  template<typename T>
  T value(const string& key, const T& default_) const {
    auto& object = _get_object();
    if (auto it = object.find(key); it != object.end()) return (T)(it->second);
    return default_;
  }
  // clang-format on

  // iteration
  // clang-format off
  json_value* begin() { return _get_array().data(); }
  const json_value* begin() const { return _get_array().data(); }
  json_value* end() { return _get_array().data() + _get_array().size(); }
  const json_value* end() const { return _get_array().data() + _get_array().size(); }
  json_object& items() { return _get_object(); }
  const json_object& items() const { return _get_object(); }
  // clang-format on

  // comparisons
  // clang-format off
  template<typename T>
  bool operator==(const T& value) const { return (T)(*this) == value; }
  bool operator==(const string& value) const { return _get_string() == value; }
  bool operator==(const char* value) const { return _get_string() == value; }
  template<typename T>
  bool operator!=(const T& value) const { return !(*this == value); }
  // clang-format on

  // setters
  void set(const json_value& other) {
    auto copy = json_value(other);
    _swap(copy);
  }
  void set(json_value&& other) { _swap(other); }
  void set(std::nullptr_t) { _set(json_type::null, _integer, 0); }
  void set(int32_t value) { _set(json_type::integer, _integer, value); }
  void set(int64_t value) { _set(json_type::integer, _integer, value); }
  void set(uint32_t value) { _set(json_type::uinteger, _uinteger, value); }
  void set(uint64_t value) { _set(json_type::uinteger, _uinteger, value); }
  void set(float value) { _set(json_type::number, _number, value); }
  void set(double value) { _set(json_type::number, _number, value); }
  void set(bool value) { _set(json_type::boolean, _boolean, value); }
  void set(const string& value) { _set(json_type::string, _string, value); }
  void set(const char* value) { _set(json_type::string, _string, value); }
  void set(const json_array& value) { _set(json_type::array, _array, value); }
  void set(const json_object& value) {
    _set(json_type::object, _object, value);
  }
  template <typename T, std::enable_if_t<std::is_enum_v<T>, bool> = true>
  void set(T value) {
    _set(json_type::string, _string, json_enum<T>::to_string(value));
  }
  template <typename T, size_t N>
  void set(const std::array<T, N>& value) {
    _set(json_type::array, _array, json_array{});
    _array->resize(value.size());
    for (auto idx = (size_t)0; idx < value.size(); idx++)
      (*_array)[idx].set(value.at(idx));
  }
  template <typename T>
  void set(const vector<T>& value) {
    _set(json_type::array, _array, json_array{});
    _array->resize(value.size());
    for (auto idx = (size_t)0; idx < value.size(); idx++)
      (*_array)[idx].set(value.at(idx));
  }
#ifdef __APPLE__
  void set(size_t value) {
    _set(json_type::uinteger, _uinteger, (uint64_t)value);
  }
#endif

  // generic get
  template <typename T>
  T get() const {
    auto value = T{};
    get(value);
    return value;
  }

  // getters
  void get(int32_t& value) const { value = _get_number<int32_t>(); }
  void get(int64_t& value) const { value = _get_number<int64_t>(); }
  void get(uint32_t& value) const { value = _get_number<uint32_t>(); }
  void get(uint64_t& value) const { value = _get_number<uint64_t>(); }
  void get(float& value) const { value = _get_number<float>(); }
  void get(double& value) const { value = _get_number<double>(); }
  void get(bool& value) const { value = _get_boolean(); }
  void get(string& value) const { value = _get_string(); }
  template <typename T, std::enable_if_t<std::is_enum_v<T>, bool> = true>
  void get(T& value) const {
    value = json_enum<T>::from_string(_get_string());
  }
  template <typename T, size_t N>
  void get(std::array<T, N>& value) const {
    auto& array = _get_array();
    if (array.size() != N)
      throw json_error{"array of fixed size expected", this};
    for (auto idx = (size_t)0; idx < value.size(); idx++)
      array[idx].get(value[idx]);
  }
  template <typename T>
  void get(vector<T>& value) const {
    auto& array = _get_array();
    value.resize(array.size());
    for (auto idx = (size_t)0; idx < value.size(); idx++)
      array[idx].get(value[idx]);
  }
#ifdef __APPLE__
  void get(size_t& value) const { value = _get_number<size_t>(); }
#endif

  // math types
  // clang-format off
  explicit json_value(const vec2i& value) : json_value((const std::array<int, 2>&)value) { }
  explicit json_value(const vec3i& value) : json_value((const std::array<int, 3>&)value) { }
  explicit json_value(const vec4i& value) : json_value((const std::array<int, 4>&)value) { }
  explicit json_value(const vec2f& value) : json_value((const std::array<float, 2>&)value) { }
  explicit json_value(const vec3f& value) : json_value((const std::array<float, 3>&)value) { }
  explicit json_value(const vec4f& value) : json_value((const std::array<float, 4>&)value) { }
  explicit json_value(const frame2f& value) : json_value((const std::array<float, 6>&)value) { }
  explicit json_value(const frame3f& value) : json_value((const std::array<float, 12>&)value) { }
  explicit json_value(const mat2f& value) : json_value((const std::array<float, 4>&)value) { }
  explicit json_value(const mat3f& value) : json_value((const std::array<float, 9>&)value) { }
  explicit json_value(const mat4f& value) : json_value((const std::array<float, 16>&)value) { }
  // clang-format on

  // math types
  // clang-format off
  void get(vec2i& value) const { get((std::array<int, 2>&)value); }
  void get(vec3i& value) const { get((std::array<int, 3>&)value); }
  void get(vec4i& value) const { get((std::array<int, 4>&)value); }
  void get(vec2f& value) const { get((std::array<float, 2>&)value); }
  void get(vec3f& value) const { get((std::array<float, 3>&)value); }
  void get(vec4f& value) const { get((std::array<float, 4>&)value); }
  void get(frame2f& value) const { get((std::array<float, 6>&)value); }
  void get(frame3f& value) const { get((std::array<float, 12>&)value); }
  void get(mat2f& value) const { get((std::array<float, 4>&)value); }
  void get(mat3f& value) const { get((std::array<float, 9>&)value); }
  void get(mat4f& value) const { get((std::array<float, 16>&)value); }
  // clang-format on

  // math types
  // clang-format off
  void set(const vec2i& value) { set((const std::array<int, 2>&)value); }
  void set(const vec3i& value) { set((const std::array<int, 3>&)value); }
  void set(const vec4i& value) { set((const std::array<int, 4>&)value); }
  void set(const vec2f& value) { set((const std::array<float, 2>&)value); }
  void set(const vec3f& value) { set((const std::array<float, 3>&)value); }
  void set(const vec4f& value) { set((const std::array<float, 4>&)value); }
  void set(const frame2f& value) { set((const std::array<float, 6>&)value); }
  void set(const frame3f& value) { set((const std::array<float, 12>&)value); }
  void set(const mat2f& value) { set((const std::array<float, 4>&)value); }
  void set(const mat3f& value) { set((const std::array<float, 9>&)value); }
  void set(const mat4f& value) { set((const std::array<float, 16>&)value); }
  // clang-format on

 private:
  json_type _type = json_type::null;
  union {
    int64_t      _integer = 0;
    uint64_t     _uinteger;
    double       _number;
    bool         _boolean;
    string*      _string;
    json_array*  _array;
    json_object* _object;
  };

  void _init() {
    switch (_type) {
      case json_type::string: _string = new string{}; break;
      case json_type::array: _array = new json_array{}; break;
      case json_type::object: _object = new json_object{}; break;
      default: break;
    }
  }

  void _copy(const json_value& other) {
    _type = other._type;
    switch (_type) {
      case json_type::null: break;
      case json_type::integer: _integer = other._integer; break;
      case json_type::uinteger: _uinteger = other._uinteger; break;
      case json_type::number: _number = other._number; break;
      case json_type::boolean: _boolean = other._boolean; break;
      case json_type::string: _string = new string{*other._string}; break;
      case json_type::array: _array = new json_array{*other._array}; break;
      case json_type::object: _object = new json_object{*other._object}; break;
    }
  }

  json_value& _swap(json_value& other) {
    std::swap(_type, other._type);
    std::swap(_integer, other._integer);
    return *this;
  }

  void _clear() {
    switch (_type) {
      case json_type::string: delete _string; break;
      case json_type::array: delete _array; break;
      case json_type::object: delete _object; break;
      default: break;
    }
    _type    = json_type::null;
    _integer = 0;
  }

  template <typename T>
  T _get_number() const {
    if (_type == json_type::integer) {
      return (T)_integer;
    } else if (_type == json_type::uinteger) {
      return (T)_uinteger;
    } else if (_type == json_type::number) {
      return (T)_number;
    } else {
      throw json_error{"number expected", this};
    }
  }

  bool _get_boolean() const {
    if (_type != json_type::boolean) throw json_error{"boolean expected", this};
    return _boolean;
  }
  const string& _get_string() const {
    if (_type != json_type::string) throw json_error{"string expected", this};
    return *_string;
  }
  string& _get_string() {
    if (_type != json_type::string) throw json_error{"string expected", this};
    return *_string;
  }
  const json_array& _get_array() const {
    if (_type != json_type::array) throw json_error{"array expected", this};
    return *_array;
  }
  json_array& _get_array() {
    if (_type != json_type::array) throw json_error{"array expected", this};
    return *_array;
  }
  const json_object& _get_object() const {
    if (_type != json_type::object) throw json_error{"object expected", this};
    return *_object;
  }
  json_object& _get_object() {
    if (_type != json_type::object) throw json_error{"object expected", this};
    return *_object;
  }

  template <typename T, typename V>
  void _set(json_type type, T& var, const V& value) {
    if (_type != type) _clear();
    _type = type;
    var   = (T)value;
  }
  template <typename T, typename V>
  void _set(json_type type, T*& var, const V& value) {
    if (_type != type) {
      _clear();
      _type = type;
      var   = new T(value);
    } else {
      *var = value;
    }
  }

  bool _empty() const {
    switch (_type) {
      case json_type::string: return _string->empty();
      case json_type::array: return _array->empty();
      case json_type::object: return _object->empty();
      default: throw json_error{"compound expected", this};
    }
  }
  size_t _size() const {
    switch (_type) {
      case json_type::string: return _string->size();
      case json_type::array: return _array->size();
      case json_type::object: return _object->size();
      default: throw json_error{"compound expected", this};
    }
  }
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON I/O
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save json
json_value load_json(const string& filename);
void       load_json(const string& filename, json_value& json);
void       save_json(const string& filename, json_value& json);

// Parse/dump json
json_value parse_json(const string& text);
void       parse_json(const string& text, json_value& json);
string     dump_json(const json_value& json);
void       dump_json(string& text, const json_value& json);

// Load json
bool load_json(const string& filename, json_value& json, string& error);
// Save json
bool save_json(const string& filename, json_value& json, string& error);

// Parse json
bool parse_json(const string& text, json_value& json, string& error);
// Dump json
bool dump_json(string& text, const json_value& json, string& error);

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

// Command line setter.
using cli_setter = bool (*)(
    const json_value&, void*, const vector<string>& choices);
// Command line variable.
struct cli_variable {
  void*                             value     = nullptr;
  cli_setter                        setter    = nullptr;
  vector<string>                    choices   = {};
  ordered_map<string, cli_variable> variables = {};
};
// Command line state.
struct cli_state {
  json_value   defaults  = {};
  json_value   schema    = {};
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
