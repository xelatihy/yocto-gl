//
// # Yocto/JSON: Utilities for manipulating JSON data
//
// Yocto/JSON is an implementation of a utilities for hanlding JSON data,
// including a Json variant data type, loading, saving and formatting JSON,
// and a parser for command line arguments to JSON.
// Compared to other libraries, it is more lightweight and provides a common
// implementation for the rest of Yocto/GL.
// Yocto/JSON is implemented in `yocto_json.h` and `yocto_json.cpp`, and
// depends on `json.hpp`.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
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

#ifndef _YOCTO_JSON_H_
#define _YOCTO_JSON_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <array>
#include <cstdio>
#include <functional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::function;
using std::pair;
using std::string;
using std::string_view;
using std::unordered_map;
using std::vector;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON DATA TYPE
// -----------------------------------------------------------------------------
namespace yocto {

// Json type
enum struct json_type {
  // clang-format off
  null, integer, unsigned_, real, boolean, string_, array, object, binary
  // clang-format on
};

// Json forward declarations
struct json_value;
using json_array      = vector<json_value>;
using json_object     = vector<pair<string, json_value>>;
using json_binary     = vector<uint8_t>;
using json_iterator   = json_value*;
using json_citerator  = const json_value*;
using json_oiterator  = pair<const string, json_value>*;
using json_ociterator = const pair<const string, json_value>*;

// Json type error
struct json_error : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// Json value
struct json_value {
  // constructors
  json_value();
  json_value(const json_value& other);
  json_value(json_value&& other);
  explicit json_value(std::nullptr_t);
  explicit json_value(int32_t);
  explicit json_value(int64_t);
  explicit json_value(uint32_t);
  explicit json_value(uint64_t);
  explicit json_value(float);
  explicit json_value(double);
  explicit json_value(bool);
  explicit json_value(const string&);
  explicit json_value(string_view);
  explicit json_value(const char* value);
  explicit json_value(const json_array&);
  explicit json_value(const json_object&);
  explicit json_value(const json_binary&);
  template <typename T>
  explicit json_value(const T& value);

  // assignments
  json_value& operator=(const json_value& other);
  json_value& operator=(json_value&& other);
  template <typename T>
  json_value& operator=(const T& value);

  // type
  json_type type() const;
  bool      is_null() const;
  bool      is_integer() const;
  bool      is_unsigned() const;
  bool      is_real() const;
  bool      is_integral() const;
  bool      is_number() const;
  bool      is_boolean() const;
  bool      is_string() const;
  bool      is_array() const;
  bool      is_object() const;
  bool      is_binary() const;

  // conversions (see get)
  explicit operator int32_t() const;
  explicit operator int64_t() const;
  explicit operator uint32_t() const;
  explicit operator uint64_t() const;
  explicit operator float() const;
  explicit operator double() const;
  explicit operator bool() const;
  explicit operator string() const;
  explicit operator string_view() const;
  template <typename T>
  explicit operator T() const;

  // get values via conversion
  template <typename T>
  T get() const;
  template <typename T>
  void get_to(T& value) const;

  // access references
  template <typename T>
  T& get_ref();
  template <typename T>
  const T& get_ref() const;

  // structure support
  bool   empty() const;
  size_t size() const;
  void   clear();
  void   resize(size_t size);
  void   reserve(size_t size);
  void   update(const json_value& other);

  // elemnt acceess
  bool              contains(const string& key) const;
  json_value&       operator[](size_t idx);
  json_value&       operator[](const string& key);
  json_value&       at(size_t idx);
  const json_value& at(size_t idx) const;
  json_value&       at(const string& key);
  const json_value& at(const string& key) const;
  json_value&       front();
  const json_value& front() const;
  json_value&       back();
  const json_value& back() const;
  void              push_back(const json_value& value);
  void              push_back(json_value&& value);
  template <typename T>
  void push_back(const T& value);
  template <typename... Args>
  json_value& emplace_back(Args&&... args);
  template <typename... Args>
  json_value& emplace(Args&&... args);

  // iteration
  json_iterator  begin();
  json_citerator begin() const;
  json_iterator  end();
  json_citerator end() const;
  struct json_object_it;
  struct json_object_cit;
  json_object_it  items();
  json_object_cit items() const;
  json_iterator   find(const string& key);
  json_citerator  find(const string& key) const;

  // get value at an object key
  template <typename T>
  T      value(const string& key, const T& default_) const;
  string value(const string& key, const char* default_) const;

  // array/object/binary creation
  static json_value array();
  static json_value object();
  static json_value binary();

  // swap
  void swap(json_value& other);

#ifdef __APPLE__
  explicit json_value(size_t);
  explicit operator size_t() const;
#endif

  // destructor
  ~json_value();

  json_type _type = json_type::null;
  union {
    int64_t      _integer;
    uint64_t     _unsigned;
    double       _real;
    bool         _boolean;
    string*      _string;
    json_array*  _array;
    json_object* _object;
    json_binary* _binary;
  };
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a json file
bool load_json(const string& filename, json_value& json, string& error);
bool save_json(const string& filename, const json_value& json, string& error);

// Formats/parse a Json to/from string
bool   parse_json(const string& text, json_value& json, string& error);
bool   format_json(string& text, const json_value& json, string& error);
string format_json(const json_value& json);

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

// Conversion shortcuts
template <typename T>
inline T from_json(const json_value& json);
template <typename T>
inline json_value to_json(const T& value);
template <typename T>
inline bool from_json(const json_value& json, T& value, json_error& error);
template <typename T>
inline bool to_json(const json_value& json, T& value, json_error& error);

// Conversion from json to values
inline void from_json(const json_value& json, int64_t& value);
inline void from_json(const json_value& json, int32_t& value);
inline void from_json(const json_value& json, uint64_t& value);
inline void from_json(const json_value& json, uint32_t& value);
inline void from_json(const json_value& json, double& value);
inline void from_json(const json_value& json, float& value);
inline void from_json(const json_value& json, bool& value);
inline void from_json(const json_value& json, string& value);
template <typename T>
inline void from_json(const json_value& json, vector<T>& value);
template <typename T, size_t N>
inline void from_json(const json_value& json, array<T, N>& value);

// Conversion to json from values
inline void to_json(json_value& json, int64_t value);
inline void to_json(json_value& json, int32_t value);
inline void to_json(json_value& json, uint64_t value);
inline void to_json(json_value& json, uint32_t value);
inline void to_json(json_value& json, double value);
inline void to_json(json_value& json, float value);
inline void to_json(json_value& json, bool value);
inline void to_json(json_value& json, const string& value);
template <typename T>
inline void to_json(json_value& json, const vector<T>& value);
template <typename T, size_t N>
inline void to_json(json_value& json, const array<T, N>& value);

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON SCHEMA
// -----------------------------------------------------------------------------
namespace yocto {

// Validate a value against a schema
bool validate_json(
    const json_value& value, const json_value& schema, string& error);
bool validate_json(const json_value& value, const json_value& schema,
    vector<string>& errors, size_t max_errors = 100);

}  // namespace yocto

// -----------------------------------------------------------------------------
// COMMAND LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// Parse the command line described by a JSON schema.
// Checks for errors, and exits on error or help.
struct json_value;
void parse_cli(
    json_value& value, const json_value& schema, int argc, const char** argv);
void parse_cli(
    json_value& value, const json_value& schema, const vector<string>& args);
// Parse the command line described by a schema.
bool parse_cli(json_value& value, const json_value& schema,
    const vector<string>& args, string& error, string& usage);
bool parse_cli(json_value& value, const json_value& schema, int argc,
    const char** argv, string& error, string& usage);

// Low-level parsing routine
bool parse_cli(json_value& value, const json_value& schema,
    const vector<string>& args, string& error, string& usage, string& command);

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON TREE DATA TYPE
// -----------------------------------------------------------------------------
namespace yocto {

// Declarations
struct json_view;
struct json_cview;

// Json tree
struct json_tree {
  union json_value {
    int64_t  _integer = 0;
    uint64_t _unsigned;
    double   _real;
    bool     _boolean;
    struct {
      uint32_t start;
      uint32_t length;
    } _string;
    struct {
      uint32_t length;
      uint32_t skip;
    } _array;
    struct {
      uint32_t length;
      uint32_t skip;
    } _object;
    struct {
      uint32_t start;
      uint32_t length;
    } _binary;
  };
  vector<json_type>  types            = {json_type::null};
  vector<json_value> values           = {json_value{}};
  vector<char>       strings          = {};
  vector<char>       keys             = {};
  vector<json_value> key_list         = {};
  vector<uint8_t>    binaries         = {};
  vector<uint32_t>   build_skip_stack = {};
  bool               valid            = true;
  string             error            = "";
};

// Load/save a json file
bool load_json(const string& filename, json_tree& json, string& error);
bool save_json(const string& filename, const json_tree& json, string& error);

// Get view from value
inline json_view  get_root(json_tree& json);
inline json_cview get_croot(json_tree& json);

// Error handling
inline void set_error(json_tree& json, string_view error);
inline void clear_error(json_tree& json);

// Json view
struct json_view {
  json_tree* root  = nullptr;
  uint32_t   index = 0;
  json_view(json_tree* root_) : root{root_}, index{(uint32_t)-1} {}
  json_view(json_tree* root_, uint32_t index_) : root{root_}, index{index_} {}
};
struct json_cview {
  json_tree* root  = nullptr;
  uint32_t   index = 0;
  json_cview(json_tree* root_) : root{root_}, index{(uint32_t)-1} {}
  json_cview(json_tree* root_, uint32_t index_) : root{root_}, index{index_} {}
  json_cview(json_view other) : root{other.root}, index{other.index} {}
};

// Error check
inline bool   is_valid(json_cview json);
inline bool   is_valid(json_view json);
inline string get_error(json_cview json);
inline string get_error(json_view json);
inline string compute_path(json_cview json);
inline bool   set_error(json_cview json, string_view error);

// Type
inline json_type get_type(json_cview json);
// Type
inline bool is_null(json_cview json);
inline bool is_integer(json_cview json);
inline bool is_unsigned(json_cview json);
inline bool is_real(json_cview json);
inline bool is_integral(json_cview json);
inline bool is_number(json_cview json);
inline bool is_boolean(json_cview json);
inline bool is_string(json_cview json);
inline bool is_array(json_cview json);
inline bool is_object(json_cview json);
inline bool is_binary(json_cview json);

// Initialization to basic types
inline bool set_null(json_view json);
inline bool set_integer(json_view json, int64_t value);
inline bool set_unsigned(json_view json, uint64_t value);
inline bool set_real(json_view json, double value);
inline bool set_boolean(json_view json, bool value);
inline bool set_string(json_view json, const string& value);
inline bool set_integral(json_view json, int64_t value);
inline bool set_integral(json_view json, int32_t value);
inline bool set_integral(json_view json, uint64_t value);
inline bool set_integral(json_view json, uint32_t value);
inline bool set_number(json_view json, double value);
inline bool set_number(json_view json, float value);

// Get basic values
inline bool get_integer(json_cview json, int64_t& value);
inline bool get_unsigned(json_cview json, uint64_t& value);
inline bool get_real(json_cview json, double& value);
inline bool get_boolean(json_cview json, bool& value);
inline bool get_string(json_cview json, string& value);
inline bool get_integral(json_cview json, int64_t& value);
inline bool get_integral(json_cview json, uint64_t& value);
inline bool get_number(json_cview json, double& value);
inline bool get_integral(json_cview json, int32_t& value);
inline bool get_integral(json_cview json, uint32_t& value);
inline bool get_number(json_cview json, float& value);

// Get basic values - ignore errors if present
inline int64_t  get_integer(json_cview json);
inline uint64_t get_unsigned(json_cview json);
inline double   get_real(json_cview json);
inline bool     get_boolean(json_cview json);
inline string   get_string(json_cview json);
inline int64_t  get_integral(json_cview json);
inline uint64_t get_uintegral(json_cview json);
inline double   get_number(json_cview json);

// Compound type
inline bool   is_empty(json_cview json);
inline size_t get_size(json_cview json);

// Array
inline bool       set_array(json_view json);
inline bool       set_array(json_view json, size_t size);
inline bool       array_size(json_cview json, size_t& size);
inline bool       has_element(json_view json, size_t idx);
inline bool       has_element(json_cview json, size_t idx);
inline json_view  get_element(json_view json, size_t idx);
inline json_cview get_element(json_cview json, size_t idx);
inline json_view  append_element(json_view json);
inline auto       iterate_array(json_view json);
inline auto       iterate_array(json_cview json);

// Object
inline bool       set_object(json_view json);
inline bool       object_size(json_cview json, size_t& size);
inline bool       has_element(json_view json, string_view key);
inline bool       has_element(json_cview json, string_view key);
inline json_view  get_element(json_view json, string_view key);
inline json_cview get_element(json_cview json, string_view key);
inline json_view  insert_element(json_view json, string_view key);
inline auto       iterate_object(json_view json);
inline auto       iterate_object(json_cview json);

// Binary
inline bool set_binary(json_view json, const json_binary& value);
inline bool get_binary(json_cview json, json_binary& value);

// Get the path of a json view
inline string compute_path(json_cview json);

// Conversion from json to values
template <typename T>
inline bool get_value(json_cview json, T& value);

// Conversion from json to values
inline bool get_value(json_cview json, int64_t& value);
inline bool get_value(json_cview json, int32_t& value);
inline bool get_value(json_cview json, uint64_t& value);
inline bool get_value(json_cview json, uint32_t& value);
inline bool get_value(json_cview json, double& value);
inline bool get_value(json_cview json, float& value);
inline bool get_value(json_cview json, bool& value);
inline bool get_value(json_cview json, string& value);
template <typename T>
inline bool get_value(json_cview json, vector<T>& value);
template <typename T, size_t N>
inline bool get_value(json_cview json, array<T, N>& value);

// Get value at a key or index
template <typename T>
inline bool get_value_at(json_cview json, string_view key, T& value);
template <typename T>
inline bool get_value_at(json_cview json, size_t idx, T& value);

// Get value at a key or nothing is key is not preesent
template <typename T>
inline bool get_value_if(json_cview json, string_view key, T& value);

// Conversion to json from values
template <typename T>
inline bool set_value(json_view json, const T& value);

// Conversion to json from values
inline bool set_value(json_view json, int64_t value);
inline bool set_value(json_view json, int32_t value);
inline bool set_value(json_view json, uint64_t value);
inline bool set_value(json_view json, uint32_t value);
inline bool set_value(json_view json, double value);
inline bool set_value(json_view json, float value);
inline bool set_value(json_view json, bool value);
inline bool set_value(json_view json, const string& value);
inline bool set_value(json_view json, const char* value);
template <typename T>
inline bool set_value(json_view json, const vector<T>& value);
template <typename T, size_t N>
inline bool set_value(json_view json, const array<T, N>& value);

// Helpers for user-defined types
inline bool check_array(json_cview json);
inline bool check_array(json_cview json, size_t size_);
inline bool check_object(json_cview json);

// Helpers for user-defined types
inline bool set_array(json_view json);
template <typename T>
inline bool set_value_at(json_view json, size_t idx, const T& value);
template <typename T>
inline bool      append_value(json_view json, const T& value);
inline json_view append_array(json_view json);
inline json_view append_object(json_view json);

// Helpers for user-defined types
inline bool set_object(json_view json);
template <typename T>
inline bool set_value_at(json_view json, string_view key, const T& value);
template <typename T>
inline bool insert_value(json_view json, string_view key, const T& value);
template <typename T>
inline bool insert_value_if(
    json_view json, string_view key, const T& value, const T& default_);
inline json_view insert_array(json_view json, string_view key);
inline json_view insert_object(json_view json, string_view key);

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// JSON DATA TYPE
// -----------------------------------------------------------------------------
namespace yocto {

// constructors
inline json_value::json_value() : _type{json_type::null}, _unsigned{0} {}
inline json_value::json_value(const json_value& other)
    : _type{json_type::null}, _unsigned{0} {
  switch (other._type) {
    case json_type::null:
      _type    = json_type::null;
      _integer = other._integer;
      break;
    case json_type::integer:
      _type    = json_type::integer;
      _integer = other._integer;
      break;
    case json_type::unsigned_:
      _type     = json_type::unsigned_;
      _unsigned = other._unsigned;
      break;
    case json_type::real:
      _type = json_type::real;
      _real = other._real;
      break;
    case json_type::boolean:
      _type    = json_type::boolean;
      _boolean = other._boolean;
      break;
    case json_type::string_:
      _type   = json_type::string_;
      _string = new string{*other._string};
      break;
    case json_type::array:
      _type  = json_type::array;
      _array = new json_array{*other._array};
      break;
    case json_type::object:
      _type   = json_type::object;
      _object = new json_object{*other._object};
      break;
    case json_type::binary:
      _type   = json_type::binary;
      _binary = new json_binary{*other._binary};
      break;
  }
}
inline json_value::json_value(json_value&& other)
    : _type{json_type::null}, _unsigned{0} {
  swap(other);
}
inline json_value::json_value(std::nullptr_t)
    : _type{json_type::null}, _unsigned{0} {}
inline json_value::json_value(int64_t value)
    : _type{json_type::integer}, _integer{value} {}
inline json_value::json_value(int32_t value)
    : _type{json_type::integer}, _integer{value} {}
inline json_value::json_value(uint64_t value)
    : _type{json_type::unsigned_}, _unsigned{value} {}
inline json_value::json_value(uint32_t value)
    : _type{json_type::unsigned_}, _unsigned{value} {}
inline json_value::json_value(double value)
    : _type{json_type::real}, _real{value} {}
inline json_value::json_value(float value)
    : _type{json_type::real}, _real{value} {}
inline json_value::json_value(bool value)
    : _type{json_type::boolean}, _boolean{value} {}
inline json_value::json_value(const string& value)
    : _type{json_type::string_}, _string{new string{value}} {}
inline json_value::json_value(string_view value)
    : _type{json_type::string_}, _string{new string{value}} {}
inline json_value::json_value(const char* value)
    : _type{json_type::string_}, _string{new string{value}} {}
inline json_value::json_value(const json_array& value)
    : _type{json_type::array}, _array{new json_array{value}} {}
inline json_value::json_value(const json_object& value)
    : _type{json_type::object}, _object{new json_object{value}} {}
inline json_value::json_value(const json_binary& value)
    : _type{json_type::binary}, _binary{new json_binary{value}} {}
template <typename T>
inline json_value::json_value(const T& value) {
  to_json(*this, value);
}
#ifdef __APPLE__
inline json_value::json_value(size_t value)
    : _type{json_type::unsigned_}, _unsigned{(uint64_t)value} {}
#endif

// assignments
inline json_value& json_value::operator=(const json_value& value) {
  auto json = json_value{value};
  this->swap(json);
  return *this;
}
inline json_value& json_value::operator=(json_value&& value) {
  this->swap(value);
  return *this;
}
template <typename T>
inline json_value& json_value::operator=(const T& value) {
  auto json = json_value{value};
  this->swap(json);
  return *this;
}

// type
inline json_type json_value::type() const { return _type; }
inline bool json_value::is_null() const { return _type == json_type::null; }
inline bool json_value::is_integer() const {
  return _type == json_type::integer;
}
inline bool json_value::is_unsigned() const {
  return _type == json_type::unsigned_;
}
inline bool json_value::is_real() const { return _type == json_type::real; }
inline bool json_value::is_integral() const {
  return _type == json_type::integer || _type == json_type::unsigned_;
}
inline bool json_value::is_number() const {
  return _type == json_type::real || _type == json_type::integer ||
         _type == json_type::unsigned_;
}
inline bool json_value::is_boolean() const {
  return _type == json_type::boolean;
}
inline bool json_value::is_string() const {
  return _type == json_type::string_;
}
inline bool json_value::is_array() const { return _type == json_type::array; }
inline bool json_value::is_object() const { return _type == json_type::object; }
inline bool json_value::is_binary() const { return _type == json_type::binary; }

// conversions
inline json_value::operator int64_t() const {
  if (_type != json_type::integer && _type != json_type::unsigned_)
    throw json_error{"integer expected"};
  return _type == json_type::integer ? (int64_t)_integer : (int64_t)_unsigned;
}
inline json_value::operator int32_t() const {
  if (_type != json_type::integer && _type != json_type::unsigned_)
    throw json_error{"integer expected"};
  return _type == json_type::integer ? (int32_t)_integer : (int32_t)_unsigned;
}
inline json_value::operator uint64_t() const {
  if (_type != json_type::integer && _type != json_type::unsigned_)
    throw json_error{"integer expected"};
  return _type == json_type::integer ? (uint64_t)_integer : (uint64_t)_unsigned;
}
inline json_value::operator uint32_t() const {
  if (_type != json_type::integer && _type != json_type::unsigned_)
    throw json_error{"integer expected"};
  return _type == json_type::integer ? (uint32_t)_integer : (uint32_t)_unsigned;
}
inline json_value::operator double() const {
  if (_type != json_type::real && _type != json_type::integer &&
      _type != json_type::unsigned_)
    throw json_error{"number expected"};
  return _type == json_type::real
             ? (double)_real
             : _type == json_type::integer ? (double)_integer
                                           : (double)_unsigned;
}
inline json_value::operator float() const {
  if (_type != json_type::real && _type != json_type::integer &&
      _type != json_type::unsigned_)
    throw json_error{"number expected"};
  return _type == json_type::real
             ? (float)_real
             : _type == json_type::integer ? (float)_integer : (float)_unsigned;
}
inline json_value::operator bool() const {
  if (_type != json_type::boolean) throw json_error{"boolean expected"};
  return _boolean;
}
inline json_value::operator string() const {
  if (_type != json_type::string_) throw json_error{"string expected"};
  return *_string;
}
inline json_value::operator string_view() const {
  if (_type != json_type::string_) throw json_error{"string expected"};
  return *_string;
}
template <typename T>
inline json_value::operator T() const {
  auto value = T{};
  from_json(*this, value);
  return value;
}
#ifdef __APPLE__
inline json_value::operator size_t() const {
  return (size_t) operator uint64_t();
}
#endif

// conversions
template <typename T>
inline T json_value::get() const {
  return operator T();
}
template <typename T>
inline void json_value::get_to(T& value) const {
  value = operator T();
}

// access
template <typename T>
inline T& json_value::get_ref() {
  static_assert(
      std::is_same_v<T, int64_t> || std::is_same_v<T, uint64_t> ||
          std::is_same_v<T, double> || std::is_same_v<T, bool> ||
          std::is_same_v<T, string> || std::is_same_v<T, json_array> ||
          std::is_same_v<T, json_object> || std::is_same_v<T, json_binary>,
      "type not in the json variant");
  if constexpr (std::is_same_v<T, int64_t>) {
    if (_type != json_type::integer) throw json_error{"integer expected"};
    return _integer;
  } else if constexpr (std::is_same_v<T, uint64_t>) {
    if (_type != json_type::unsigned_) throw json_error{"unsigned expected"};
    return _unsigned;
  } else if constexpr (std::is_same_v<T, double>) {
    if (_type != json_type::real) throw json_error{"real expected"};
    return _real;
  } else if constexpr (std::is_same_v<T, bool>) {
    if (_type != json_type::boolean) throw json_error{"boolean expected"};
    return _boolean;
  } else if constexpr (std::is_same_v<T, string>) {
    if (_type != json_type::string_) throw json_error{"string expected"};
    return *_string;
  } else if constexpr (std::is_same_v<T, json_array>) {
    if (_type != json_type::array) throw json_error{"array expected"};
    return *_array;
  } else if constexpr (std::is_same_v<T, json_object>) {
    if (_type != json_type::object) throw json_error{"object expected"};
    return *_object;
  } else if constexpr (std::is_same_v<T, json_binary>) {
    if (_type != json_type::binary) throw json_error{"binary expected"};
    return *_binary;
  } else {
    // will never get here
  }
}
// access
template <typename T>
inline const T& json_value::get_ref() const {
  return ((json_value*)this)->get_ref<T>();  // const cast
}

// structure support
inline bool json_value::empty() const {
  switch (_type) {
    case json_type::null: return true;
    case json_type::array: return _array->empty();
    case json_type::object: return _object->empty();
    default: return false;
  }
}
inline size_t json_value::size() const {
  switch (_type) {
    case json_type::null: return 0;
    case json_type::array: return _array->size();
    case json_type::object: return _object->size();
    default: return 1;
  }
}
inline void json_value::clear() {
  switch (_type) {
    case json_type::null: _unsigned = 0; break;
    case json_type::integer: _integer = 0; break;
    case json_type::unsigned_: _unsigned = 0; break;
    case json_type::real: _real = 0; break;
    case json_type::boolean: _boolean = false; break;
    case json_type::string_: _string->clear(); break;
    case json_type::array: _array->clear(); break;
    case json_type::object: _object->clear(); break;
    case json_type::binary: _binary->clear(); break;
    default: break;
  }
}
inline void json_value::resize(size_t size) {
  switch (_type) {
    case json_type::array: _array->resize(size); break;
    default: throw json_error{"array expected"};
  }
}
inline void json_value::reserve(size_t size) {
  switch (_type) {
    case json_type::array: _array->reserve(size); break;
    case json_type::object: _object->reserve(size); break;
    default: throw json_error{"structure expected"};
  }
}

// array support
inline json_value& json_value::operator[](size_t idx) {
  if (_type == json_type::null) *this = json_array{};
  if (_type != json_type::array) throw json_error{"array expected"};
  if (idx >= _array->size()) throw json_error{"index out of range"};
  return _array->operator[](idx);
}
inline json_value& json_value::operator[](const string& key) {
  if (_type == json_type::null) *this = json_object{};
  if (_type != json_type::object) throw json_error{"object expected"};
  if (auto elem = find(key); elem != nullptr) return *elem;
  return _object->emplace_back(key, json_value{}).second;
}
inline json_value& json_value::at(size_t idx) {
  if (_type != json_type::array) throw json_error{"array expected"};
  if (idx >= _array->size()) throw json_error{"index out of range"};
  return _array->operator[](idx);
}
inline const json_value& json_value::at(size_t idx) const {
  if (_type != json_type::array) throw json_error{"array expected"};
  if (idx >= _array->size()) throw json_error{"index out of range"};
  return _array->operator[](idx);
}
inline json_value& json_value::at(const string& key) {
  if (_type != json_type::object) throw json_error{"object expected"};
  if (auto elem = find(key); elem != nullptr) return *elem;
  throw json_error{"missing key " + key};
}
inline const json_value& json_value::at(const string& key) const {
  if (_type != json_type::object) throw json_error{"object expected"};
  if (auto elem = find(key); elem != nullptr) return *elem;
  throw json_error{"missing key " + key};
}
inline json_value& json_value::front() {
  if (_type != json_type::array) throw json_error{"array expected"};
  return _array->front();
}
inline const json_value& json_value::front() const {
  if (_type != json_type::array) throw json_error{"array expected"};
  return _array->front();
}
inline json_value& json_value::back() {
  if (_type != json_type::array) throw json_error{"array expected"};
  return _array->back();
}
inline const json_value& json_value::back() const {
  if (_type != json_type::array) throw json_error{"array expected"};
  return _array->back();
}
inline void json_value::push_back(const json_value& value) {
  if (_type == json_type::null) *this = json_array{};
  if (_type != json_type::array) throw json_error{"array expected"};
  _array->push_back(value);
}
inline void json_value::push_back(json_value&& value) {
  if (_type == json_type::null) *this = json_array{};
  if (_type != json_type::array) throw json_error{"array expected"};
  _array->push_back(std::move(value));
}
template <typename T>
inline void json_value::push_back(const T& value) {
  return push_back(json_value{value});
}
template <typename... Args>
inline json_value& json_value::emplace_back(Args&&... args) {
  if (_type == json_type::null) *this = json_array{};
  if (_type != json_type::array) throw json_error{"array expected"};
  return _array->emplace_back(std::forward<Args>(args)...);
}
template <typename... Args>
inline json_value& json_value::emplace(Args&&... args) {
  if (_type == json_type::null) *this = json_object{};
  if (_type != json_type::object) throw json_error{"object expected"};
  return _array->emplace_back(std::forward<Args>(args)...);
}
inline void json_value::update(const json_value& other) {
  if (_type == json_type::null) *this = json_object{};
  if (_type != json_type::object) throw json_error{"object expected"};
  if (other._type != json_type::object) throw json_error{"object expected"};
  for (auto& [key, value] : *other._object) this->operator[](key) = value;
}

// Iteration
inline json_value* json_value::begin() {
  if (_type != json_type::array) throw json_error{"array expected"};
  return _array->data();
}
inline const json_value* json_value::begin() const {
  if (_type != json_type::array) throw json_error{"array expected"};
  return _array->data();
}
inline json_value* json_value::end() {
  if (_type != json_type::array) throw json_error{"array expected"};
  return _array->data() + _array->size();
}
inline const json_value* json_value::end() const {
  if (_type != json_type::array) throw json_error{"array expected"};
  return _array->data() + _array->size();
}
struct json_value::json_object_it {
  pair<string, json_value>* _begin;
  pair<string, json_value>* _end;
  pair<string, json_value>* begin() { return _begin; }
  pair<string, json_value>* end() { return _end; }
};
struct json_value::json_object_cit {
  const pair<string, json_value>* _begin;
  const pair<string, json_value>* _end;
  const pair<string, json_value>* begin() { return _begin; }
  const pair<string, json_value>* end() { return _end; }
};
inline json_value::json_object_it json_value::items() {
  if (_type != json_type::object) throw json_error{"object expected"};
  return {_object->data(), _object->data() + _object->size()};
}
inline json_value::json_object_cit json_value::items() const {
  return {_object->data(), _object->data() + _object->size()};
}
inline json_value* json_value::find(const string& key) {
  if (_type != json_type::object) throw json_error{"object expected"};
  for (auto& [key_, value] : *_object) {
    if (key_ == key) return &value;
  }
  return nullptr;
}
inline const json_value* json_value::find(const string& key) const {
  if (_type != json_type::object) throw json_error{"object expected"};
  for (auto& [key_, value] : *_object) {
    if (key_ == key) return &value;
  }
  return nullptr;
}
inline bool json_value::contains(const string& key) const {
  return find(key) != nullptr;
}

// get value at an object key
template <typename T>
inline T json_value::value(const string& key, const T& default_) const {
  if (_type != json_type::object) throw json_error{"object expected"};
  auto element = find(key);
  return element ? element->get<T>() : default_;
}
inline string json_value::value(const string& key, const char* default_) const {
  return value<string>(key, default_);
}

// array/object/binary creation
inline json_value json_value::array() { return json_value{json_array{}}; }
inline json_value json_value::object() { return json_value{json_object{}}; }
inline json_value json_value::binary() { return json_value{json_binary{}}; }

// swap
inline void json_value::swap(json_value& other) {
  std::swap(_type, other._type);
  std::swap(_unsigned, other._unsigned);  // hask to swap bits
}

// destructor
inline json_value::~json_value() {
  switch (_type) {
    case json_type::string_: delete _string; break;
    case json_type::array: delete _array; break;
    case json_type::object: delete _object; break;
    case json_type::binary: delete _binary; break;
    default: break;
  }
  _type     = json_type::null;
  _unsigned = 0;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

// Conversion shortcuts
template <typename T>
inline T from_json(const json_value& json) {
  auto value = T{};
  from_json(json, value);
  return value;
}
template <typename T>
inline json_value to_json(const T& value) {
  auto json = json_value{};
  to_json(json, value);
  return json;
}
template <typename T>
inline bool from_json(const json_value& json, T& value, json_error& error) {
  try {
    from_json(json, value);
    return true;
  } catch (json_error& err) {
    error = err;
    return false;
  }
}
template <typename T>
inline bool to_json(const json_value& json, T& value, json_error& error) {
  try {
    to_json(json, value);
    return true;
  } catch (json_error& err) {
    error = err;
    return false;
  }
}

// Conversion from json to values
inline void from_json(const json_value& json, int64_t& value) {
  value = (int64_t)json;
}
inline void from_json(const json_value& json, int32_t& value) {
  value = (int32_t)json;
}
inline void from_json(const json_value& json, uint64_t& value) {
  value = (uint64_t)json;
}
inline void from_json(const json_value& json, uint32_t& value) {
  value = (uint32_t)json;
}
inline void from_json(const json_value& json, double& value) {
  value = (double)json;
}
inline void from_json(const json_value& json, float& value) {
  value = (float)json;
}
inline void from_json(const json_value& json, bool& value) {
  value = (bool)json;
}
inline void from_json(const json_value& json, string& value) {
  value = (string)json;
}
template <typename T>
inline void from_json(const json_value& json, vector<T>& value) {
  value.clear();
  value.reserve(json.size());
  for (auto& ejs : json) from_json(ejs, value.emplace_back());
}
template <typename T, size_t N>
inline void from_json(const json_value& json, array<T, N>& value) {
  if (json.size() != value.size()) throw json_error{"wrong array size"};
  for (auto idx = (size_t)0; idx < value.size(); idx++)
    from_json(json.at(idx), value.at(idx));
}

// Conversion to json from values
inline void to_json(json_value& json, int64_t value) { json = value; }
inline void to_json(json_value& json, int32_t value) { json = value; }
inline void to_json(json_value& json, uint64_t value) { json = value; }
inline void to_json(json_value& json, uint32_t value) { json = value; }
inline void to_json(json_value& json, double value) { json = value; }
inline void to_json(json_value& json, float value) { json = value; }
inline void to_json(json_value& json, bool value) { json = value; }
inline void to_json(json_value& json, const string& value) { json = value; }
template <typename T>
inline void to_json(json_value& json, const vector<T>& value) {
  json = json_array{};
  for (auto& v : value) to_json(json.emplace_back(), v);
}
template <typename T, size_t N>
inline void to_json(json_value& json, const array<T, N>& value) {
  json = json_array{};
  json.resize(value.size());
  for (auto idx = (size_t)0; idx < value.size(); idx++)
    to_json(json.at(idx), value.at(idx));
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON TREE DATA TYPE
// -----------------------------------------------------------------------------
namespace yocto {

// Get view from value
inline json_view  get_root(json_tree& json) { return {&json, 0}; }
inline json_cview get_croot(json_tree& json) { return {&json, 0}; }
inline bool       set_error(json_cview json, string_view error) {
  if (!is_valid(json)) return false;
  set_error(*json.root, string{error} + " at " + compute_path(json));
  return false;
}
inline json_view set_error_view(json_cview json, string_view error) {
  if (!is_valid(json)) return {json.root};
  set_error(*json.root, string{error} + " at " + compute_path(json));
  return {json.root};
}

// Error handling
inline void set_error(json_tree& json, string_view error) {
  if (!json.valid) return;
  json.valid = false;
  json.error = string{error};
}
inline void clear_error(json_tree& json) {
  json.valid = true;
  json.error = "";
}

// Helpers
inline json_type& _get_type(json_view json) {
  if (!is_valid(json)) throw std::invalid_argument{"bad json"};
  return json.root->types[json.index];
}
inline const json_type& _get_type(json_cview json) {
  if (!is_valid(json)) throw std::invalid_argument{"bad json"};
  return json.root->types[json.index];
}
inline json_tree::json_value& _get_value(json_view json) {
  if (!is_valid(json)) throw std::invalid_argument{"bad json"};
  return json.root->values[json.index];
}
inline const json_tree::json_value& _get_value(json_cview json) {
  if (!is_valid(json)) throw std::invalid_argument{"bad json"};
  return json.root->values[json.index];
}
inline uint32_t _get_capacity(uint32_t length) {
  if (length == 0) return 0;
  if (length <= 4) return 4;
  // TODO(fabio): faster pow2
  auto capacity = (uint32_t)4;
  while (capacity < length) capacity *= 2;
  return capacity;
}
inline string_view _get_key(json_cview json) {
  if (!is_valid(json)) throw std::invalid_argument{"bad tree"};
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::string_) throw std::invalid_argument{"bad key"};
  return {json.root->keys.data() + jsv._string.start, jsv._string.length};
}
inline void _find_path(json_view json, vector<json_view>& path);

// Error check
inline bool is_valid(json_cview json) {
  return json.root != nullptr && json.root->valid && json.index != (uint32_t)-1;
}
inline bool is_valid(json_view json) {
  return json.root != nullptr && json.root->valid && json.index != (uint32_t)-1;
}
inline string get_error(json_cview json) {
  if (json.root == nullptr) return "bad root";
  if (json.root->valid) return "";
  return json.root->error;
}
inline string get_error(json_view json) {
  if (json.root == nullptr) return "bad root";
  if (json.root->valid) return "";
  return json.root->error;
}

// Type
inline json_type get_type(json_cview json) {
  if (!is_valid(json)) return json_type::null;
  return _get_type(json);
}
inline bool is_null(json_cview json) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  return jst == json_type::null;
}
inline bool is_integer(json_cview json) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  return jst == json_type::integer;
}
inline bool is_unsigned(json_cview json) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  return jst == json_type::unsigned_;
}
inline bool is_real(json_cview json) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  return jst == json_type::real;
}
inline bool is_integral(json_cview json) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  return jst == json_type::integer || jst == json_type::unsigned_;
}
inline bool is_number(json_cview json) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  return jst == json_type::integer || jst == json_type::unsigned_ ||
         jst == json_type::real;
}
inline bool is_boolean(json_cview json) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  return jst == json_type::boolean;
}
inline bool is_string(json_cview json) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  return jst == json_type::string_;
}
inline bool is_array(json_cview json) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  return jst == json_type::array;
}
inline bool is_object(json_cview json) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  return jst == json_type::object;
}
inline bool is_binary(json_cview json) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  return jst == json_type::binary;
}

// Initialization to basic types
inline bool set_null(json_view json) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  jst       = json_type::null;
  return true;
}
inline bool set_integer(json_view json, int64_t value) {
  if (!is_valid(json)) return false;
  auto& jst    = _get_type(json);
  auto& jsv    = _get_value(json);
  jst          = json_type::integer;
  jsv._integer = value;
  return true;
}
inline bool set_unsigned(json_view json, uint64_t value) {
  if (!is_valid(json)) return false;
  auto& jst     = _get_type(json);
  auto& jsv     = _get_value(json);
  jst           = json_type::unsigned_;
  jsv._unsigned = value;
  return true;
}
inline bool set_real(json_view json, double value) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  jst       = json_type::real;
  jsv._real = value;
  return true;
}
inline bool set_boolean(json_view json, bool value) {
  if (!is_valid(json)) return false;
  auto& jst    = _get_type(json);
  auto& jsv    = _get_value(json);
  jst          = json_type::boolean;
  jsv._boolean = value;
  return true;
}
inline bool set_string(json_view json, const string& value) {
  if (!is_valid(json)) return false;
  auto& jst          = _get_type(json);
  auto& jsv          = _get_value(json);
  jst                = json_type::string_;
  jsv._string.start  = (uint32_t)json.root->strings.size();
  jsv._string.length = (uint32_t)value.size();
  json.root->strings.insert(
      json.root->strings.end(), value.begin(), value.end());
  json.root->strings.push_back(0);
  return true;
}
inline bool set_integral(json_view json, int64_t value) {
  return set_integer(json, value);
}
inline bool set_integral(json_view json, int32_t value) {
  return set_integer(json, value);
}
inline bool set_integral(json_view json, uint64_t value) {
  return set_unsigned(json, value);
}
inline bool set_integral(json_view json, uint32_t value) {
  return set_unsigned(json, value);
}
inline bool set_number(json_view json, double value) {
  return set_real(json, value);
}
inline bool set_number(json_view json, float value) {
  return set_real(json, value);
}

// Get basic values
inline bool get_integer(json_cview json, int64_t& value) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::integer) return set_error(json, "integer expected");
  value = jsv._integer;
  return true;
}
inline bool get_unsigned(json_cview json, uint64_t& value) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::unsigned_) return set_error(json, "unsigned expected");
  value = jsv._unsigned;
  return true;
}
inline bool get_real(json_cview json, double& value) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::real) return set_error(json, "real expected");
  value = jsv._real;
  return true;
}
inline bool get_boolean(json_cview json, bool& value) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::boolean) return set_error(json, "boolean expected");
  value = jsv._boolean;
  return true;
}
inline bool get_string(json_cview json, string& value) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::string_) return set_error(json, "string expected");
  value = string{json.root->strings.data() + jsv._string.start,
      json.root->strings.data() + jsv._string.start + jsv._string.length};
  return true;
}
inline bool get_integral(json_cview json, int64_t& value) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::integer && jst != json_type::unsigned_)
    return set_error(json, "integer expected");
  value = (jst == json_type::integer) ? (int64_t)jsv._integer
                                      : (int64_t)jsv._unsigned;
  return true;
}
inline bool get_integral(json_cview json, uint64_t& value) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::integer && jst != json_type::unsigned_)
    return set_error(json, "integer expected");
  value = (jst == json_type::integer) ? (uint64_t)jsv._integer
                                      : (uint64_t)jsv._unsigned;
  return true;
}
inline bool get_number(json_cview json, double& value) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::real && jst != json_type::integer &&
      jst != json_type::unsigned_)
    return set_error(json, "number expected");
  value = (jst == json_type::real)
              ? (double)jsv._real
              : (jst == json_type::integer) ? (double)jsv._integer
                                            : (double)jsv._unsigned;
  return true;
}
inline bool get_integral(json_cview json, int32_t& value) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::integer && jst != json_type::unsigned_)
    return set_error(json, "integer expected");
  value = (jst == json_type::integer) ? (int32_t)jsv._integer
                                      : (int32_t)jsv._unsigned;
  return true;
}
inline bool get_integral(json_cview json, uint32_t& value) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::integer && jst != json_type::unsigned_)
    return set_error(json, "integer expected");
  value = (jst == json_type::integer) ? (uint32_t)jsv._integer
                                      : (uint32_t)jsv._unsigned;
  return true;
}
inline bool get_number(json_cview json, float& value) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::real && jst != json_type::integer &&
      jst != json_type::unsigned_)
    return set_error(json, "number expected");
  value = (jst == json_type::real)
              ? (float)jsv._real
              : (jst == json_type::integer) ? (float)jsv._integer
                                            : (float)jsv._unsigned;
  return true;
}

// Get basic values
inline int64_t get_integer(json_cview json) {
  auto value = (int64_t)0;
  return get_integer(json, value) ? value : 0;
}
inline uint64_t get_unsigned(json_cview json) {
  auto value = (uint64_t)0;
  return get_unsigned(json, value) ? value : 0;
}
inline double get_real(json_cview json) {
  auto value = (double)0;
  return get_real(json, value) ? value : 0;
}
inline bool get_boolean(json_cview json) {
  auto value = false;
  return get_boolean(json, value) ? value : false;
}
inline string get_string(json_cview json) {
  auto value = string{};
  return get_string(json, value) ? value : string{};
}
inline int64_t get_integral(json_cview json) {
  auto value = (int64_t)0;
  return get_integral(json, value) ? value : 0;
}
inline double get_number(json_cview json) {
  auto value = (double)0;
  return get_number(json, value) ? value : 0;
}

// Compound type
inline bool   is_empty(json_cview json);
inline size_t get_size(json_cview json);

// Array iteeration
inline auto iterate_array(json_view json) {
  struct iterator {
    json_view json;
    bool      operator!=(const iterator& other) {
      return is_valid(json) && json.index != other.json.index;
    }
    iterator& operator++() {
      if (!is_valid(json)) return *this;
      auto& jst = _get_type(json);
      auto& jsv = _get_value(json);
      json.index += 1;
      if (jst == json_type::array) json.index += jsv._array.skip;
      if (jst == json_type::object) json.index += jsv._object.skip;
      return *this;
    }
    json_view operator*() const { return json; }
  };
  struct iterator_wrapper {
    json_view begin_;
    json_view end_;
    iterator  begin() { return {begin_}; }
    iterator  end() { return {end_}; }
  };
  if (!is_valid(json)) return iterator_wrapper{{json.root}, {json.root}};
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::array) {
    set_error_view(json, "array expected");
    return iterator_wrapper{{json.root}, {json.root}};
  }
  return iterator_wrapper{{json.root, json.index + 1},
      {json.root, json.index + 1 + jsv._array.skip}};
}
inline auto iterate_array(json_cview json) {
  struct iterator {
    json_cview json;
    bool       operator!=(const iterator& other) {
      return is_valid(json) && json.index != other.json.index;
    }
    iterator& operator++() {
      if (!is_valid(json)) return *this;
      auto& jst = _get_type(json);
      auto& jsv = _get_value(json);
      json.index += 1;
      if (jst == json_type::array) json.index += jsv._array.skip;
      if (jst == json_type::object) json.index += jsv._object.skip;
      return *this;
    }
    json_cview operator*() const { return json; }
  };
  struct iterator_wrapper {
    json_cview begin_;
    json_cview end_;
    iterator   begin() { return {begin_}; }
    iterator   end() { return {end_}; }
  };
  if (!is_valid(json)) return iterator_wrapper{{json.root}, {json.root}};
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::array) {
    set_error_view(json, "array expected");
    return iterator_wrapper{{json.root}, {json.root}};
  }
  return iterator_wrapper{
      {json.root, json.index + 1},
      {json.root, json.index + 1 + jsv._array.skip},
  };
}

// Array
inline bool set_array(json_view json) {
  if (!is_valid(json)) return false;
  if (json.index != json.root->values.size() - 1)
    throw std::out_of_range{"can only add at the end"};
  auto& jst  = _get_type(json);
  auto& jsv  = _get_value(json);
  jst        = json_type::array;
  jsv._array = {0, 0};
  return true;
}
inline bool array_size(json_cview json, size_t& size) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::array) return set_error(json, "array expected");
  size = (size_t)jsv._array.length;
  return true;
}
inline json_view get_element(json_view json, size_t idx) {
  if (!is_valid(json)) return {json.root};
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::array) return set_error_view(json, "array expected");
  if (idx >= jsv._array.length)
    return set_error_view(json, "index out of bounds");
  if (jsv._array.length == jsv._array.skip) {
    return {json.root, json.index + 1 + (uint32_t)idx};
  } else {
    auto count = 0;
    for (auto ejs : iterate_array(json)) {
      if (count++ == idx) return ejs;
    }
    return set_error_view(json, "index out of bounds");
  }
}
inline json_cview get_element(json_cview json, size_t idx) {
  if (!is_valid(json)) return {json.root};
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::array) return set_error_view(json, "array expected");
  if (idx >= jsv._array.length)
    return set_error_view(json, "index out of bounds");
  if (jsv._array.length == jsv._array.skip) {
    return {json.root, json.index + 1 + (uint32_t)idx};
  } else {
    auto count = 0;
    for (auto ejs : iterate_array(json)) {
      if (count++ == idx) return ejs;
    }
    return set_error_view(json, "index out of bounds");
  }
}
inline json_view append_element(json_view json) {
  if (!is_valid(json)) return {json.root};
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::array) return set_error_view(json, "array expected");
  if (json.index + 1 + jsv._array.skip != json.root->values.size())
    throw std::out_of_range{"can only add at the end"};
  jsv._array.length += 1;
  auto index = (uint32_t)json.root->values.size();
  json.root->types.emplace_back(json_type::null);
  json.root->values.emplace_back();
  auto stack = vector<json_view>{};
  _find_path(json, stack);
  for (auto jss : stack) {
    auto& jsst = _get_type(jss);
    auto& jssv = _get_value(jss);
    if (jsst == json_type::array) {
      jssv._array.skip += 1;
    } else if (jsst == json_type::object) {
      jssv._object.skip += 1;
    } else {
      throw std::runtime_error{"bad stack"};
    }
  }
  return {json.root, index};
}

// Object iteration
inline auto iterate_object(json_view json) {
  struct iterator {
    json_view json;
    bool      operator!=(const iterator& other) {
      return is_valid(json) && json.index != other.json.index;
    }
    iterator& operator++() {
      if (!is_valid(json)) return *this;
      auto& jst = _get_type(json_view{json.root, json.index + 1});
      auto& jsv = _get_value(json_view{json.root, json.index + 1});
      json.index += 2;
      if (jst == json_type::array) json.index += jsv._array.skip;
      if (jst == json_type::object) json.index += jsv._object.skip;
      return *this;
    }
    pair<string_view, json_view> operator*() const {
      return {_get_key(json), json_view{json.root, json.index + 1}};
    }
  };
  struct iterator_wrapper {
    json_view begin_;
    json_view end_;
    iterator  begin() { return {begin_}; }
    iterator  end() { return {end_}; }
  };
  if (!is_valid(json)) return iterator_wrapper{{json.root}, {json.root}};
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::object) {
    set_error_view(json, "object expected");
    return iterator_wrapper{{json.root}, {json.root}};
  }
  return iterator_wrapper{{json.root, json.index + 1},
      {json.root, json.index + 1 + jsv._object.skip}};
}
inline auto iterate_object(json_cview json) {
  struct iterator {
    json_cview json;
    bool       operator!=(const iterator& other) {
      return is_valid(json) && json.index != other.json.index;
    }
    iterator& operator++() {
      if (!is_valid(json)) return *this;
      auto& jst = _get_type(json_cview{json.root, json.index + 1});
      auto& jsv = _get_value(json_cview{json.root, json.index + 1});
      json.index += 2;
      if (jst == json_type::array) json.index += jsv._array.skip;
      if (jst == json_type::object) json.index += jsv._object.skip;
      return *this;
    }
    pair<string_view, json_cview> operator*() const {
      return {_get_key(json), json_cview{json.root, json.index + 1}};
    }
  };
  struct iterator_wrapper {
    json_cview begin_;
    json_cview end_;
    iterator   begin() { return {begin_}; }
    iterator   end() { return {end_}; }
  };
  if (!is_valid(json)) return iterator_wrapper{{json.root}, {json.root}};
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::object) {
    set_error_view(json, "object expected");
    return iterator_wrapper{{json.root}, {json.root}};
  }
  return iterator_wrapper{{json.root, json.index + 1},
      {json.root, json.index + 1 + jsv._object.skip}};
}

// Object
inline bool set_object(json_view json) {
  if (!is_valid(json)) return false;
  if (json.index != json.root->values.size() - 1)
    throw std::out_of_range{"can only add at the end"};
  auto& jst  = _get_type(json);
  auto& jsv  = _get_value(json);
  jst        = json_type::object;
  jsv._array = {0, 0};
  return true;
}
inline bool object_size(json_cview json, size_t& size) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::object) return set_error(json, "object expected");
  size = (size_t)jsv._object.length;
  return true;
}
inline json_view get_element(json_view json, string_view key) {
  if (!is_valid(json)) return {json.root};
  auto& jst = _get_type(json);
  if (jst != json_type::object) return set_error_view(json, "object expected");
  for (auto [okey, ejs] : iterate_object(json)) {
    if (okey == key) return ejs;
  }
  return set_error_view(json, "missing key " + string{key});
}
inline json_cview get_element(json_cview json, string_view key) {
  if (!is_valid(json)) return {json.root};
  auto& jst = _get_type(json);
  if (jst != json_type::object) return set_error_view(json, "object expected");
  for (auto [okey, value] : iterate_object(json)) {
    if (okey == key) return value;
  }
  return set_error_view(json, "missing key " + string{key});
}
inline bool has_element(json_view json, string_view key) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  if (jst != json_type::object) return set_error(json, "object expected");
  for (auto [okey, value] : iterate_object(json)) {
    if (okey == key) return true;
  }
  return false;
}
inline bool has_element(json_cview json, string_view key) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  if (jst != json_type::object) return set_error(json, "object expected");
  for (auto [okey, value] : iterate_object(json)) {
    if (okey == key) return true;
  }
  return false;
}
inline json_view insert_element(json_view json, string_view key) {
  if (!is_valid(json)) return {json.root};
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::object) return set_error_view(json, "object expected");
  if (json.index + 1 + jsv._object.skip != json.root->values.size())
    throw std::out_of_range{"can only add at the end"};
  for (auto [okey, ejs] : iterate_object(json)) {
    if (okey == key) return ejs;
  }
  jsv._object.length += 1;
  auto& jkt = json.root->types.emplace_back();
  auto& jkv = json.root->values.emplace_back();
  jkt       = json_type::string_;
  for (auto kv : json.root->key_list) {
    auto okey = string_view{
        json.root->keys.data() + kv._string.start, kv._string.length};
    if (okey == key) jkv = kv;
  }
  if (jkv._string.length == 0) {
    jkv._string.start  = (uint32_t)json.root->keys.size();
    jkv._string.length = (uint32_t)key.size();
    json.root->keys.insert(json.root->keys.end(), key.begin(), key.end());
    json.root->keys.push_back(0);
  }
  auto index = (uint32_t)json.root->values.size();
  json.root->types.emplace_back(json_type::null);
  json.root->values.emplace_back();
  auto stack = vector<json_view>{};
  _find_path(json, stack);
  for (auto jss : stack) {
    auto& jsst = _get_type(jss);
    auto& jssv = _get_value(jss);
    if (jsst == json_type::array) {
      jssv._array.skip += 2;
    } else if (jsst == json_type::object) {
      jssv._object.skip += 2;
    } else {
      throw std::runtime_error{"bad stack"};
    }
  }
  return {json.root, index};
}

// Binary
inline bool set_binary(json_view json, const json_binary& value) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  // TODO(fabio): implement reuse
  jst                = json_type::binary;
  jsv._binary.start  = (uint32_t)json.root->binaries.size();
  jsv._binary.length = (uint32_t)value.size();
  json.root->binaries.insert(
      json.root->binaries.end(), value.begin(), value.end());
  return true;
}
inline bool get_binary(json_cview json, json_binary& value) {
  if (!is_valid(json)) return false;
  auto& jst = _get_type(json);
  auto& jsv = _get_value(json);
  if (jst != json_type::binary) return set_error(json, "binary expected");
  value = json_binary{json.root->binaries.data() + jsv._binary.start,
      json.root->binaries.data() + jsv._binary.start + jsv._binary.length};
  return true;
}

// Get the path of a json view
inline bool _compute_path(json_cview json, json_cview jsv, string& path) {
  if (!is_valid(json) || !is_valid(jsv)) {
    return false;
  } else if (json.index == jsv.index) {
    path = "/";
    return true;
  } else if (is_array(json)) {
    auto idx = 0;
    for (auto ejs : iterate_array(json)) {
      if (!_compute_path(ejs, jsv, path)) continue;
      if (path.back() == '/') path.pop_back();
      path = "/"s + std::to_string(idx) + path;
      return true;
    }
    return false;
  } else if (is_object(json)) {
    for (auto [key, ejs] : iterate_object(json)) {
      if (!_compute_path(ejs, jsv, path)) continue;
      if (path.back() == '/') path.pop_back();
      path = "/" + string{key} + path;
      return true;
    }
    return false;
  } else {
    return false;
  }
}
inline string compute_path(json_cview json) {
  auto path = string{};
  if (_compute_path({json.root, 0}, json, path)) {
    return path;
  } else {
    return "";
  }
}

// Conversion from json to values
inline bool get_value(json_cview json, int64_t& value) {
  return get_integral(json, value);
}
inline bool get_value(json_cview json, int32_t& value) {
  return get_integral(json, value);
}
inline bool get_value(json_cview json, uint64_t& value) {
  return get_integral(json, value);
}
inline bool get_value(json_cview json, uint32_t& value) {
  return get_integral(json, value);
}
inline bool get_value(json_cview json, double& value) {
  return get_number(json, value);
}
inline bool get_value(json_cview json, float& value) {
  return get_number(json, value);
}
inline bool get_value(json_cview json, bool& value) {
  return get_boolean(json, value);
}
inline bool get_value(json_cview json, string& value) {
  return get_string(json, value);
}
template <typename T>
inline bool get_value(json_cview json, vector<T>& value) {
  if (!is_valid(json)) return false;
  if (!is_array(json)) return set_error(json, "array expected");
  value.clear();
  auto size = (size_t)0;
  array_size(json, size);
  value.reserve(size);
  for (auto ejs : iterate_array(json)) {
    if (!get_value(ejs, value.emplace_back())) return false;
  }
  return true;
}
template <typename T, size_t N>
inline bool get_value(json_cview json, array<T, N>& value) {
  if (!is_valid(json)) return false;
  if (!is_array(json)) return set_error(json, "array expected");
  auto size = (size_t)0;
  array_size(json, size);
  if (size != N) return set_error(json, "size mismatched");
  auto idx = 0;
  for (auto ejs : iterate_array(json)) {
    if (!get_value(ejs, value.at(idx++))) return false;
  }
  return true;
}

// Get value at a key or index
template <typename T>
inline bool get_value_at(json_cview json, string_view key, T& value) {
  if (!is_valid(json)) return false;
  auto element = get_element(json, key);
  return get_value(element, value);
}
template <typename T>
inline bool get_value_at(json_cview json, size_t idx, T& value) {
  if (!is_valid(json)) return false;
  auto element = get_element(json, idx);
  if (!is_valid(element)) return false;
  return get_value(element, value);
}

// Get value at a key or nothing is key is not preesent
template <typename T>
inline bool get_value_if(json_cview json, string_view key, T& value) {
  if (!is_valid(json)) return false;
  if (!has_element(json, key)) return true;
  auto element = get_element(json, key);
  if (!is_valid(element)) return false;
  return get_value(element, value);
}

// Conversion to json from values
inline bool set_value(json_view json, int64_t value) {
  return set_integral(json, value);
}
inline bool set_value(json_view json, int32_t value) {
  return set_integral(json, value);
}
inline bool set_value(json_view json, uint64_t value) {
  return set_integral(json, value);
}
inline bool set_value(json_view json, uint32_t value) {
  return set_integral(json, value);
}
inline bool set_value(json_view json, double value) {
  return set_number(json, value);
}
inline bool set_value(json_view json, float value) {
  return set_number(json, value);
}
inline bool set_value(json_view json, bool value) {
  return set_boolean(json, value);
}
inline bool set_value(json_view json, const string& value) {
  return set_string(json, value);
}
inline bool set_value(json_view json, const char* value) {
  return set_string(json, value);
}
template <typename T>
inline bool set_value(json_view json, const vector<T>& value) {
  if (!set_array(json)) return false;
  for (auto& v : value) {
    if (!set_value(append_element(json), v)) return false;
  }
  return true;
}
template <typename T, size_t N>
inline bool set_value(json_view json, const array<T, N>& value) {
  if (!set_array(json)) return false;
  for (auto& v : value) {
    if (!set_value(append_element(json), v)) return false;
  }
  return true;
}

// Helpers for user-defined types
inline bool check_array(json_cview json) {
  if (!is_valid(json)) return false;
  if (!is_array(json)) return set_error(json, "array expected");
  return true;
}
inline bool check_array(json_cview json, size_t size_) {
  if (!is_valid(json)) return false;
  if (!is_array(json)) return set_error(json, "array expected");
  auto size = (size_t)0;
  if (!array_size(json, size) || size != size_)
    return set_error(json, "mismatched size");
  return true;
}
inline bool check_object(json_cview json) {
  if (!is_valid(json)) return false;
  if (!is_object(json)) return set_error(json, "array expected");
  return true;
}

// Helpers for user-defined types
template <typename T>
inline bool set_value_at(json_view json, size_t idx, const T& value) {
  if (!is_valid(json)) return false;
  auto ejs = get_element(json, idx);
  if (!is_valid(ejs)) return false;
  return set_value(ejs, value);
}
template <typename T>
inline bool append_value(json_view json, const T& value) {
  if (!is_valid(json)) return false;
  auto ejs = append_element(json);
  if (!is_valid(ejs)) return false;
  return set_value(ejs, value);
}
inline json_view append_array(json_view json) {
  if (!is_valid(json)) return {json.root};
  auto ejs = append_element(json);
  if (!is_valid(ejs)) return {json.root};
  if (!set_array(ejs)) return {json.root};
  return ejs;
}
inline json_view append_object(json_view json) {
  if (!is_valid(json)) return {json.root};
  auto ejs = append_element(json);
  if (!is_valid(ejs)) return {json.root};
  if (!set_object(ejs)) return {json.root};
  return ejs;
}

// Helpers for user-defined types
template <typename T>
inline bool set_value_at(json_view json, string_view key, const T& value) {
  if (!is_valid(json)) return false;
  auto ejs = get_element(json, key);
  if (!is_valid(ejs)) return false;
  return set_value(ejs, value);
}
template <typename T>
inline bool insert_value(json_view json, string_view key, const T& value) {
  if (!is_valid(json)) return false;
  auto ejs = insert_element(json, key);
  if (!is_valid(ejs)) return false;
  return set_value(ejs, value);
}
template <typename T>
inline bool insert_value_if(
    json_view json, string_view key, const T& value, const T& default_) {
  if (!is_valid(json)) return false;
  if (value == default_) return true;
  auto ejs = insert_element(json, key);
  if (!is_valid(ejs)) return false;
  return set_value(ejs, value);
}
inline json_view insert_array(json_view json, string_view key) {
  if (!is_valid(json)) return {json.root};
  auto ejs = insert_element(json, key);
  if (!is_valid(ejs)) return {json.root};
  if (!set_array(ejs)) return {json.root};
  return ejs;
}
inline json_view insert_object(json_view json, string_view key) {
  if (!is_valid(json)) return {json.root};
  auto ejs = insert_element(json, key);
  if (!is_valid(ejs)) return {json.root};
  if (!set_object(ejs)) return {json.root};
  return ejs;
}

// Helpers that need to be declared here
inline bool _find_anchestors(
    json_view json, uint32_t index, vector<json_view>& path) {
  auto& jst = _get_type(json);
  if (jst == json_type::array) {
    if (json.index == index) {
      path.push_back(json);
      return true;
    }
    //    if (index <= json.index || index >= json.index + 1 + jsv._array.skip)
    //      return false;
    //    path.push_back(json.index);
    for (auto ejs : iterate_array(json)) {
      if (_find_anchestors(ejs, index, path)) {
        path.push_back(ejs);
        return true;
      }
    }
    return false;
  } else if (jst == json_type::object) {
    if (json.index == index) {
      path.push_back(json);
      return true;
    }
    //    if (index <= json.index || index >= json.index + 1 + jsv._object.skip)
    //      return false;
    //    path.push_back(json.index);
    for (auto [okey, ejs] : iterate_object(json)) {
      if (_find_anchestors(ejs, index, path)) {
        path.push_back(json);
        return true;
      }
    }
    return false;
  } else {
    return false;
  }
}
inline void _find_path(json_view json, vector<json_view>& path) {
  path.clear();
  _find_anchestors({json.root, 0}, json.index, path);
}

}  // namespace yocto

#endif
