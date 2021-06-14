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

// Json error
struct json_error : std::logic_error {
  json_error(const string& error) : std::logic_error(error) {}
};

// Json type
enum struct json_type {
  // clang-format off
  null, integer, uinteger, number, boolean, string, array, object
  // clang-format on
};

// Json typdefs
struct json_value;
using json_array  = vector<json_value>;
using json_object = ordered_map<string, json_value>;

// Json value
struct json_value {
  // Json value
  json_value() : _type{json_type::null} {}
  json_value(json_type type) : _type{type} { _init(); }
  json_value(const json_value& other) { _copy(other); }
  json_value(json_value&& value) : json_value() { _swap(value); }
  json_value& operator=(json_value other) { return _swap(other); }
  ~json_value() { _clear(); }

  // conversions
  explicit json_value(std::nullptr_t) : _type{json_type::null}, _integer{0} {}
  explicit json_value(int32_t value)
      : _type{json_type::integer}, _integer{value} {}
  explicit json_value(int64_t value)
      : _type{json_type::integer}, _integer{value} {}
  explicit json_value(uint32_t value)
      : _type{json_type::uinteger}, _uinteger{value} {}
  explicit json_value(uint64_t value)
      : _type{json_type::uinteger}, _uinteger{value} {}
  explicit json_value(float value)  //
      : _type{json_type::number}, _number{value} {}
  explicit json_value(double value)
      : _type{json_type::number}, _number{value} {}
  explicit json_value(bool value)
      : _type{json_type::boolean}, _boolean{value} {}
  explicit json_value(const string& value)
      : _type{json_type::string}, _string{new string{value}} {}
  explicit json_value(const char* value)
      : _type{json_type::string}, _string{new string{value}} {}
  explicit json_value(const json_array& value)
      : _type{json_type::array}, _array{new json_array{value}} {}
  explicit json_value(const json_object& value)
      : _type{json_type::object}, _object{new json_object{value}} {}
  template <typename T, size_t N>
  explicit json_value(const array<T, N>& value)
      : _type{json_type::array}
      , _array{new json_array(value.data(), value.data() + value.size())} {}
  template <typename T>
  explicit json_value(const vector<T>& value)
      : _type{json_type::array}
      , _array{new json_array(value.data(), value.data() + value.size())} {}

  // casts
  explicit operator int32_t() const { return _get_number<int32_t>(); }
  explicit operator int64_t() const { return _get_number<int64_t>(); }
  explicit operator uint32_t() const { return _get_number<uint32_t>(); }
  explicit operator uint64_t() const { return _get_number<uint64_t>(); }
  explicit operator float() const { return _get_number<float>(); }
  explicit operator double() const { return _get_number<double>(); }
  explicit operator bool() const { return _get_boolean(); }
  explicit operator string() const { return _get_string(); }

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

  // set
  void set_null() { _set_type(json_type::null); }
  void set_integer(int64_t value) {
    _set_type(json_type::integer);
    _integer = value;
  }
  void set_uinteger(uint64_t value) {
    _set_type(json_type::uinteger);
    _uinteger = value;
  }
  void set_number(double value) {
    _set_type(json_type::number);
    _number = value;
  }
  void set_boolean(bool value) {
    _set_type(json_type::boolean);
    _boolean = value;
  }
  void set_string(const char* value) {
    _set_type(json_type::string);
    *_string = value;
  }
  void set_string(const string& value) {
    _set_type(json_type::string);
    *_string = value;
  }
  void set_array() {
    _set_type(json_type::array);
    *_array = {};
  }
  void set_array(size_t size) {
    _set_type(json_type::array);
    _array->resize(size);
  }
  void set_array(const json_array& value) {
    _set_type(json_type::array);
    *_array = value;
  }
  void set_object() {
    _set_type(json_type::object);
    *_object = {};
  }
  void set_object(const json_object& value) {
    _set_type(json_type::object);
    *_object = value;
  }

  int64_t& get_integer() {
    if (_type != json_type::integer) throw json_error{"integer expected"};
    return _integer;
  }
  uint64_t& get_uinteger() {
    if (_type != json_type::uinteger) throw json_error{"integer expected"};
    return _uinteger;
  }
  double& get_number() {
    if (_type != json_type::number) throw json_error{"number expected"};
    return _number;
  }
  bool& get_boolean() {
    if (_type != json_type::boolean) throw json_error{"boolean expected"};
    return _boolean;
  }
  string& get_string() {
    if (_type != json_type::string) throw json_error{"string expected"};
    return *_string;
  }
  json_array& get_array() {
    if (_type != json_type::array) throw json_error{"array expected"};
    return *_array;
  }
  json_object& get_object() {
    if (_type != json_type::object) throw json_error{"object expected"};
    return *_object;
  }

  const int64_t& get_integer() const {
    if (_type != json_type::integer) throw json_error{"integer expected"};
    return _integer;
  }
  const uint64_t& get_uinteger() const {
    if (_type != json_type::uinteger) throw json_error{"integer expected"};
    return _uinteger;
  }
  const double& get_number() const {
    if (_type != json_type::number) throw json_error{"number expected"};
    return _number;
  }
  const bool& get_boolean() const {
    if (_type != json_type::boolean) throw json_error{"boolean expected"};
    return _boolean;
  }
  const string& get_string() const {
    if (_type != json_type::string) throw json_error{"string expected"};
    return *_string;
  }
  const json_array& get_array() const {
    if (_type != json_type::array) throw json_error{"array expected"};
    return *_array;
  }
  const json_object& get_object() const {
    if (_type != json_type::object) throw json_error{"object expected"};
    return *_object;
  }

  // size
  bool empty() const {
    if (is_array()) return _array->empty();
    if (is_object()) return _object->empty();
    throw json_error{"array or object expected"};
  }
  size_t size() const {
    if (is_array()) return _array->size();
    if (is_object()) return _object->size();
    throw json_error{"array or object expected"};
  }

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
  json_value& insert(const string& key) { return _get_object()[key]; }
  json_value& append() { return _get_array().emplace_back(); }
  json_value& emplace_back() { return _get_array().emplace_back(); }
  void push_back(const json_value& item) { _get_array().push_back(item); }
  template<typename T>
  void push_back(const T& value) { _get_array().push_back(json_value{value}); }
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

  // setters
  // clang-format off
  void set(std::nullptr_t) { _set_value(json_type::null, _integer, 0); }
  void set(int32_t value) { _set_value(json_type::integer, _integer, value); }
  void set(int64_t value) { _set_value(json_type::integer, _integer, value); }
  void set(uint32_t value) { _set_value(json_type::uinteger, _uinteger, value); }
  void set(uint64_t value) { _set_value(json_type::uinteger, _uinteger, value); }
  void set(float value) { _set_value(json_type::number, _number, value); }
  void set(double value) { _set_value(json_type::number, _number, value); }
  void set(bool value) { _set_value(json_type::boolean, _boolean, value); }
  void set(const string& value) { _set_ptr(json_type::string, _string, value); }
  void set(const char* value) { _set_ptr(json_type::string, _string, value); }
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
  // clang-format on

  template <typename T>
  T get() const {
    auto value = T{};
    get(value);
    return value;
  }

  // conversions
  void get(int32_t& value) const { value = _get_number<int32_t>(); }
  void get(int64_t& value) const { value = _get_number<int64_t>(); }
  void get(uint32_t& value) const { value = _get_number<uint32_t>(); }
  void get(uint64_t& value) const { value = _get_number<uint64_t>(); }
  void get(float& value) const { value = _get_number<float>(); }
  void get(double& value) const { value = _get_number<double>(); }
  void get(bool& value) const { value = _get_boolean(); }
  void get(string& value) const { value = _get_string(); }
  template <typename T, size_t N>
  void get(std::array<T, N>& value) const {
    auto& array = _get_array();
    if (array.size() != N)
      throw json_error{"array of size " + std::to_string(N) + " expected"};
    for (auto idx = (size_t)0; idx < value.size(); idx++)
      array.at(idx).get(value.at(idx));
  }
  template <typename T>
  inline void get(vector<T>& value) const {
    auto& array = _get_array();
    value.resize(array.size());
    for (auto idx = (size_t)0; idx < value.size(); idx++)
      array.at(idx).get(value.at(idx));
  }

#ifdef __APPLE__
  void set(size_t value) { set((uint64_t)value); }
  void get(size_t& value) const { value = get<uint64_t>(); }
#endif

  template <typename T>
  json_value& insert(const string& key, const T& value) {
    auto& item = insert(key);
    item.set(value);
    return item;
  }
  template <typename T>
  json_value& append(const T& value) {
    auto& item = append();
    item.set(value);
    return item;
  }

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
  }

  void _set_type(json_type type) {
    if (_type == type) return;
    auto new_json = json_value{type};
    _swap(new_json);
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
      throw json_error{"number expected"};
    }
  }

  bool _get_boolean() const {
    if (_type != json_type::boolean) throw json_error{"boolean expected"};
    return _boolean;
  }
  const string& _get_string() const {
    if (_type != json_type::string) throw json_error{"string expected"};
    return *_string;
  }
  string& _get_string() {
    if (_type != json_type::string) throw json_error{"string expected"};
    return *_string;
  }
  const json_array& _get_array() const {
    if (_type == json_type::array) throw json_error{"array expected"};
    return *_array;
  }
  json_array& _get_array() {
    if (_type == json_type::array) throw json_error{"array expected"};
    return *_array;
  }
  const json_object& _get_object() const {
    if (_type != json_type::object) throw json_error{"object expected"};
    return *_object;
  }
  json_object& _get_object() {
    if (_type != json_type::object) throw json_error{"object expected"};
    return *_object;
  }

  template <typename T, typename V>
  void _set_value(json_type type, T& var, const V& value) {
    if (_type != type) _clear();
    _type = type;
    var   = (T)value;
  }
  template <typename T, typename V>
  void _set_ptr(json_type type, T* var, const V& value) {
    if (_type != type) {
      _clear();
      _type = type;
      var   = new T(value);
    } else {
      *var = value;
    }
  }
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON I/O
// -----------------------------------------------------------------------------
namespace yocto {

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
