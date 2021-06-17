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
// parse arguments, throw exceptions on error
void parse_cli(cli_state& cli, const vector<string>& args);
// parse arguments, checks for errors
bool parse_cli(cli_state& cli, const vector<string>& args, string& error);
// parse arguments, checks for errors, and exits on error or help
void parse_cli_and_handle_errors(cli_state& cli, int argc, const char** argv);
void parse_cli_and_handle_errors(cli_state& cli, const vector<string>& args);

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
template <typename K, typename T>
struct ordered_map {
 public:
  using container      = vector<pair<K, T>>;
  using const_iterator = typename container::const_iterator;
  using iterator       = typename container::iterator;
  using value_type     = typename container::value_type;

  ordered_map() {}
  ordered_map(const ordered_map& other) : _data{other._data} {}
  ordered_map(ordered_map&& other) : _data{std::move(other._data)} {};
  ordered_map& operator=(const ordered_map& other) {
    _data.operator=(other._data);
    return *this;
  }
  ordered_map& operator=(ordered_map&& other) {
    _data.operator=(std::move(other._data));
    return *this;
  }

  size_t size() const { return _data.size(); }
  bool   empty() const { return _data.empty(); }

  bool contains(const K& key) const { return _find(key) != _data.end(); }

  T& operator[](const K& key) {
    if (auto it = _find(key); it != _data.end()) return it->second;
    return _data.emplace_back(key, T{}).second;
  }
  const T& operator[](const K& key) const {
    if (auto it = _search(key); it != _data.end()) return it->second;
    throw std::out_of_range{"missing key for " + key};
  }
  T& at(const K& key) {
    if (auto it = _find(key); it != _data.end()) return it->second;
    throw std::out_of_range{"missing key for " + key};
  }
  const T& at(const K& key) const {
    if (auto it = _find(key); it != _data.end()) return it->second;
    throw std::out_of_range{"missing key for " + key};
  }

  iterator       find(const string& key) { return _find(key); }
  const_iterator find(const string& key) const { return _find(key); }

  iterator       begin() { return _data.begin(); }
  iterator       end() { return _data.end(); }
  const_iterator begin() const { return _data.begin(); }
  const_iterator end() const { return _data.end(); }

  void push_back(value_type&& item) { _data.push_back(std::move(item)); }
  void push_back(const value_type& item) { _data.push_back(item); }
  template <typename... Args>
  value_type& emplace_back(Args&&... args) {
    return _data.emplace_back(std::forward<Args>(args)...);
  }

 private:
  container _data;

  iterator _find(const string& key) {
    for (auto it = _data.begin(); it != _data.end(); ++it)
      if (key == it->first) return it;
    return _data.end();
  }
  const_iterator _find(const string& key) const {
    for (auto it = _data.begin(); it != _data.end(); ++it)
      if (key == it->first) return it;
    return _data.end();
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
// define  static const vector<string>& labels();
template <typename T>
struct json_enum_trait;

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
 public:
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
  explicit json_value(T&& value) : json_value() {
    set(std::forward<T>(value));
  }
  template <typename T>
  explicit operator T() const {
    auto value = T{};
    get(value);
    return value;
  }
  template <typename T>
  json_value& operator=(T&& value) {
    set(std::forward<T>(value));
    return *this;
  }

  // type
  json_type type() const { return _type; }
  bool      is_null() const { return _type == json_type::null; }
  bool      is_integer() const {
    return _type == json_type::integer || _type == json_type::uinteger;
  }
  bool is_number() const {
    return _type == json_type::integer || _type == json_type::uinteger ||
           _type == json_type::number;
  }
  bool is_boolean() const { return _type == json_type::boolean; }
  bool is_string() const { return _type == json_type::string; }
  bool is_array() const { return _type == json_type::array; }
  bool is_object() const { return _type == json_type::object; }

  // size
  bool   empty() const { return _empty(); }
  size_t size() const { return _size(); }

  // object creation
  static json_value array() { return json_value{json_array{}}; }
  static json_value array(size_t size) { return json_value{json_array(size)}; }
  static json_value object() { return json_value{json_object{}}; }

  // element access
  json_value&       operator[](size_t idx) { return _get_array().at(idx); }
  const json_value& operator[](size_t idx) const {
    return _get_array().at(idx);
  }
  json_value&       operator[](const string& key) { return _get_object()[key]; }
  const json_value& operator[](const string& key) const {
    return _get_object().at(key);
  }
  json_value&       at(size_t idx) { return _get_array().at(idx); }
  const json_value& at(size_t idx) const { return _get_array().at(idx); }
  json_value&       at(const string& key) { return _get_object().at(key); }
  const json_value& at(const string& key) const {
    return _get_object().at(key);
  }
  bool contains(const string& key) const { return _get_object().contains(key); }
  json_value&       front() { return _get_array().front(); }
  const json_value& front() const { return _get_array().front(); }
  json_value&       back() { return _get_array().back(); }
  const json_value& back() const { return _get_array().back(); }
  template <typename... Args>
  json_value& emplace_back(Args&&... args) {
    return _get_array().emplace_back(std::forward<Args>(args)...);
  }
  void push_back(json_value&& item) { _get_array().push_back(std::move(item)); }
  void push_back(const json_value& item) { _get_array().push_back(item); }
  template <typename T>
  void push_back(const T& value) {
    _get_array().push_back(json_value{value});
  }
  json_value& insert_back(string&& key) {
    return _get_object().emplace_back(std::move(key), json_value{}).second;
  }
  json_value& insert_back(const string& key) {
    return _get_object().emplace_back(key, json_value{}).second;
  }

  // get functions
  template <typename T>
  T get() const {
    auto value = T{};
    get(value);
    return value;
  }
  template <typename T>
  void try_get(const string& key, T& value) const {
    auto& object = _get_object();
    auto  it     = object.find(key);
    if (it != object.end()) it->second.get(value);
  }
  template <typename T>
  T get_at(const string& key) const {
    return _get_object().at(key).get<T>();
  }
  template <typename T>
  T get_or(const string& key, const T& default_) const {
    auto& object = _get_object();
    if (auto it = object.find(key); it != object.end())
      return it->second.get<T>();
    return default_;
  }
  template <typename T>
  T value(const string& key, const T& default_) const {
    auto& object = _get_object();
    if (auto it = object.find(key); it != object.end())
      return it->second.get<T>();
    return default_;
  }
  string value(const string& key, const char* default_) const {
    auto& object = _get_object();
    if (auto it = object.find(key); it != object.end())
      return it->second.get<string>();
    return default_;
  }

  // iteration
  json_value*       begin() { return _get_array().data(); }
  const json_value* begin() const { return _get_array().data(); }
  json_value*       end() { return _get_array().data() + _get_array().size(); }
  const json_value* end() const {
    return _get_array().data() + _get_array().size();
  }
  json_object&       items() { return _get_object(); }
  const json_object& items() const { return _get_object(); }

  // comparisons
  template <typename T>
  bool operator==(const T& value) const {
    return (T)(*this) == value;
  }
  template <typename T>
  bool operator!=(const T& value) const {
    return !(*this == value);
  }
  bool operator==(const string& value) const { return _get_string() == value; }
  bool operator==(const char* value) const { return _get_string() == value; }

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
  void set(json_array&& value) {
    _set(json_type::array, _array, std::move(value));
  }
  void set(const json_object& value) {
    _set(json_type::object, _object, value);
  }
  void set(json_object&& value) {
    _set(json_type::object, _object, std::move(value));
  }
  template <typename T, std::enable_if_t<std::is_enum_v<T>, bool> = true>
  void set(T value) {
    _set(json_type::string, _string,
        json_enum_trait<T>::labels().at((int)value));
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
    auto  values = get<string>();
    auto& labels = json_enum_trait<T>::labels();
    for (auto idx = 0; idx < (int)labels.size(); idx++)
      if (labels[idx] == values) {
        value = (T)idx;
        return;
      }
    throw json_error{"missing label", this};
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

  // math types
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
  template <typename T, typename V>
  void _set(json_type type, T*& var, V&& value) {
    if (_type != type) {
      _clear();
      _type = type;
      var   = new T(std::move(value));
    } else {
      *var = std::move(value);
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
void       save_json(const string& filename, const json_value& json);

// Parse/dump json
json_value parse_json(const string& text);
void       parse_json(const string& text, json_value& json);
string     format_json(const json_value& json);
void       format_json(string& text, const json_value& json);

// Load/save json
bool load_json(const string& filename, json_value& json, string& error);
bool save_json(const string& filename, const json_value& json, string& error);

// Parse/dump json
bool parse_json(const string& text, json_value& json, string& error);
bool format_json(string& text, const json_value& json, string& error);

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
using cli_setter = void (*)(
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
  json_value   value     = {};
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
  cli_error(const string& message) : std::runtime_error(message) {}
};
struct cli_help : std::runtime_error {
  cli_help(const string& message) : std::runtime_error(message) {}
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
