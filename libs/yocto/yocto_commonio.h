//
// # Yocto/CommonIO: Utilities for writing command-line apps
//
// Yocto/CommonIO is a collection of utilities used in writing command-line
// applications, including parsing command line arguments, simple path
// manipulation, file lading and saving, and printing values, timers and
// progress bars.
// Yocto/CommonIO is implemented in `yocto_commonio.h` and `yocto_commonio.cpp`.
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

#ifndef _YOCTO_COMMONIO_H_
#define _YOCTO_COMMONIO_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <array>
#include <cstdio>
#include <functional>
#include <stdexcept>
#include <string>
#include <type_traits>
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
using std::vector;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// PRINT/FORMATTING UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Print a message to the console
void print_info(const string& msg);
// Prints a message to the console and exit with an error.
void print_fatal(const string& msg);

// Timer that prints as scope end. Create with `print_timed` and print with
// `print_elapsed`.
struct print_timer {
  int64_t start_time = -1;
  ~print_timer();  // print time if scope ends
};
// Print traces for timing and program debugging
print_timer print_timed(const string& msg);
int64_t     print_elapsed(print_timer& timer);

// Print progress
void print_progress(const string& message, int current, int total);

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
// parse arguments and checks for errors
bool parse_cli(cli_state& cli, int argc, const char** argv, string& error);
// gets usage message
string get_usage(const cli_state& cli);
// gets whether help was invoked
bool get_help(const cli_state& cli);

// Parses an optional or positional argument. Optional arguments' names start
// with "--" or "-", otherwise they are arguments. Supports strings, numbers,
// boolean flags and enums.
// Many names, separated by commas, can be used for each argument.
// Boolean flags are indicated with a pair of names "--name/--no-name", so that
// both options are explicitly specified.
template <typename T>
inline void add_option(cli_state& cli, const string& name, T& value,
    const string& usage, bool req = false);
// Parses an optional or positional argument where values can only be within a
// set of choices. Supports strings, integers and enums.
template <typename T>
inline void add_option(cli_state& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, bool req = false);
// Parse all arguments left on the command line. Can only be used as argument.
template <typename T>
inline void add_option(cli_state& cli, const string& name, vector<T>& value,
    const string& usage, bool req = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// PATH UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Utility to normalize a path
string normalize_path(const string& filename);

// Get directory name (not including '/').
string path_dirname(const string& filename);

// Get extension (including '.').
string path_extension(const string& filename);

// Get filename without directory.
string path_filename(const string& filename);

// Get filename without directory and extension.
string path_basename(const string& filename);

// Joins paths
string path_join(const string& patha, const string& pathb);
string path_join(const string& patha, const string& pathb, const string& pathc);

// Replaces extensions
string replace_extension(const string& filename, const string& ext);

// Check if a file can be opened for reading.
bool path_exists(const string& filename);

// Check if a file is a directory
bool path_isdir(const string& filename);

// Check if a file is a file
bool path_isfile(const string& filename);

// List the contents of a directory
vector<string> list_directory(const string& filename);

// Create a directory and all missing parent directories if needed
bool make_directory(const string& dirname, string& error);

// Get the current directory
string path_current();

}  // namespace yocto

// -----------------------------------------------------------------------------
// FILE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a text file
bool load_text(const string& filename, string& str, string& error);
bool save_text(const string& filename, const string& str, string& error);

// Using directive
using byte = unsigned char;

// Load/save a binary file
bool load_binary(const string& filename, vector<byte>& data, string& error);
bool save_binary(
    const string& filename, const vector<byte>& data, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON SUPPORT
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
using json_array  = vector<json_value>;
using json_object = vector<pair<string, json_value>>;
using json_binary = vector<uint8_t>;

// Json type error
struct json_type_error : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// Json value
struct json_value {
  // constructors
  json_value();
  json_value(const json_value& other);
  json_value(json_value&& other);
  explicit json_value(std::nullptr_t);
  explicit json_value(int64_t value);
  explicit json_value(int32_t value);
  explicit json_value(uint64_t value);
  explicit json_value(uint32_t value);
  explicit json_value(double value);
  explicit json_value(float value);
  explicit json_value(bool value);
  explicit json_value(const string& value);
  explicit json_value(const char* value);
  explicit json_value(const json_array& value);
  explicit json_value(const json_object& value);
  explicit json_value(const json_binary& value);

  // assignments
  json_value& operator=(const json_value& other);
  json_value& operator=(json_value&& other);
  json_value& operator=(std::nullptr_t);
  json_value& operator=(int64_t value);
  json_value& operator=(int32_t value);
  json_value& operator=(uint64_t value);
  json_value& operator=(uint32_t value);
  json_value& operator=(double value);
  json_value& operator=(float value);
  json_value& operator=(bool value);
  json_value& operator=(const string& value);
  json_value& operator=(const char* value);
  json_value& operator=(const json_array& value);
  json_value& operator=(const json_object& value);
  json_value& operator=(const json_binary& value);

  // type
  void      set_type(json_type type);
  json_type type() const;
  bool      is_null() const;
  bool      is_integer() const;
  bool      is_unsigned() const;
  bool      is_real() const;
  bool      is_bool() const;
  bool      is_string() const;
  bool      is_array() const;
  bool      is_object() const;
  bool      is_binary() const;

  // conversions
  explicit operator int64_t() const;
  explicit operator int32_t() const;
  explicit operator uint64_t() const;
  explicit operator uint32_t() const;
  explicit operator double() const;
  explicit operator float() const;
  explicit operator bool() const;
  explicit operator string() const;
  explicit operator json_array() const;
  explicit operator json_object() const;
  explicit operator json_binary() const;

// size_t fix
#ifdef __APPLE__
  explicit json_value(size_t value);
  explicit    operator size_t() const;
  json_value& operator=(size_t value);
#endif

  // access
  const int64_t&     get_integer() const;
  const uint64_t&    get_unsigned() const;
  const double&      get_real() const;
  const bool&        get_boolean() const;
  const string&      get_string() const;
  const json_array&  get_array() const;
  const json_object& get_object() const;
  const json_binary& get_binary() const;
  int64_t&           get_integer();
  uint64_t&          get_unsigned();
  double&            get_real();
  bool&              get_boolean();
  string&            get_string();
  json_array&        get_array();
  json_object&       get_object();
  json_binary&       get_binary();

  // structure support
  bool   empty() const;
  size_t size() const;
  void   resize(size_t size);
  void   reserve(size_t size);

  // array support
  static json_value array();
  json_value&       operator[](size_t idx);
  const json_value& operator[](size_t idx) const;
  json_value&       at(size_t idx);
  const json_value& at(size_t idx) const;
  json_value&       front();
  const json_value& front() const;
  json_value&       back();
  const json_value& back() const;
  void              push_back(const json_value& value);
  void              push_back(json_value&& value);
  template <typename... Args>
  json_value&       emplace_back(Args&&... args);
  json_value*       begin();
  const json_value* begin() const;
  json_value*       end();
  const json_value* end() const;

  // object support
  static json_value object();
  json_value&       operator[](const string& key);
  json_value&       at(const string& key);
  const json_value& at(const string& key) const;
  struct json_object_it;
  struct json_object_cit;
  json_object_it    items();
  json_object_cit   items() const;
  json_value*       find(const string& key);
  const json_value* find(const string& key) const;
  bool              contains(const string& key) const;

  // binary support
  static json_value binary();

  // swap
  void swap(json_value& other);

  // destructor
  ~json_value();

 private:
  json_type _type = json_type::null;
  union {
    int64_t      _integer;
    uint64_t     _unsigned;
    double       _real;
    bool         _boolean;
    string*      _string_;
    json_array*  _array;
    json_object* _object;
    json_binary* _binary;
  };
};

// Load/save a json file
bool load_json(const string& filename, json_value& js, string& error);
bool save_json(const string& filename, const json_value& js, string& error);

// Conversion shortcuts
template <typename T>
inline T& from_json(const json_value& js);
template <typename T>
inline json_value to_json(const T& value);

// Conversion from json to types
inline void from_json(const json_value& js, int64_t& value);
inline void from_json(const json_value& js, int32_t& value);
inline void from_json(const json_value& js, uint64_t& value);
inline void from_json(const json_value& js, uint32_t& value);
inline void from_json(const json_value& js, double& value);
inline void from_json(const json_value& js, float& value);
inline void from_json(const json_value& js, bool& value);
inline void from_json(const json_value& js, string& value);
template <typename T>
inline void from_json(const json_value& js, vector<T>& value);
template <typename T, size_t N>
inline void from_json(const json_value& js, array<T, N>& value);

// Conversion from types to json
inline void to_json(json_value& js, int64_t value);
inline void to_json(json_value& js, int32_t value);
inline void to_json(json_value& js, uint64_t value);
inline void to_json(json_value& js, uint32_t value);
inline void to_json(json_value& js, double value);
inline void to_json(json_value& js, float value);
inline void to_json(json_value& js, bool value);
inline void to_json(json_value& js, const string& value);
template <typename T>
inline void to_json(json_value& js, const vector<T>& value);
template <typename T, size_t N>
inline void to_json(json_value& js, const array<T, N>& value);

// Conversion from json to types
inline void from_json(const json_value& js, int64_t& value, string& error,
    const json_value& root);
inline void from_json(const json_value& js, int32_t& value, string& error,
    const json_value& root);
inline void from_json(const json_value& js, uint64_t& value, string& error,
    const json_value& root);
inline void from_json(const json_value& js, uint32_t& value, string& error,
    const json_value& root);
inline void from_json(const json_value& js, double& value, string& error);
inline void from_json(const json_value& js, float& value, string& error);
inline void from_json(const json_value& js, bool& value, string& error);
inline void from_json(const json_value& js, string& value, string& error);
template <typename T>
inline void from_json(const json_value& js, vector<T>& value, string& error,
    const json_value& root);
template <typename T, size_t N>
inline void from_json(const json_value& js, array<T, N>& value, string& error,
    const json_value& root);

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// FORMATTING
// -----------------------------------------------------------------------------
namespace yocto {

// This is a very crude replacement for `std::format()` that will be used when
// available on all platforms.
template <typename... Args>
inline string format(const string& fmt, Args&&... args);

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON SUPPORT
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
      _type    = json_type::string_;
      _string_ = new string{*other._string_};
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
    : _type{json_type::string_}, _string_{new string{value}} {}
inline json_value::json_value(const char* value)
    : _type{json_type::string_}, _string_{new string{value}} {}
inline json_value::json_value(const json_array& value)
    : _type{json_type::array}, _array{new json_array{value}} {}
inline json_value::json_value(const json_object& value)
    : _type{json_type::object}, _object{new json_object{value}} {}
inline json_value::json_value(const json_binary& value)
    : _type{json_type::binary}, _binary{new json_binary{value}} {}

// assignments
inline json_value& json_value::operator=(const json_value& value) {
  auto js = json_value{value};
  swap(js);
  return *this;
}
inline json_value& json_value::operator=(json_value&& value) {
  swap(value);
  return *this;
}
inline json_value& json_value::operator=(std::nullptr_t) {
  auto js = json_value{nullptr};
  swap(js);
  return *this;
}
inline json_value& json_value::operator=(int64_t value) {
  auto js = json_value{value};
  swap(js);
  return *this;
}
inline json_value& json_value::operator=(int32_t value) {
  auto js = json_value{value};
  swap(js);
  return *this;
}
inline json_value& json_value::operator=(uint64_t value) {
  auto js = json_value{value};
  swap(js);
  return *this;
}
inline json_value& json_value::operator=(uint32_t value) {
  auto js = json_value{value};
  swap(js);
  return *this;
}
inline json_value& json_value::operator=(double value) {
  auto js = json_value{value};
  swap(js);
  return *this;
}
inline json_value& json_value::operator=(float value) {
  auto js = json_value{value};
  swap(js);
  return *this;
}
inline json_value& json_value::operator=(bool value) {
  auto js = json_value{value};
  swap(js);
  return *this;
}
inline json_value& json_value::operator=(const string& value) {
  auto js = json_value{value};
  swap(js);
  return *this;
}
inline json_value& json_value::operator=(const char* value) {
  auto js = json_value{value};
  swap(js);
  return *this;
}
inline json_value& json_value::operator=(const json_array& value) {
  auto js = json_value{value};
  swap(js);
  return *this;
}
inline json_value& json_value::operator=(const json_object& value) {
  auto js = json_value{value};
  swap(js);
  return *this;
}
inline json_value& json_value::operator=(const json_binary& value) {
  auto js = json_value{value};
  swap(js);
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
inline bool json_value::is_bool() const { return _type == json_type::boolean; }
inline bool json_value::is_string() const {
  return _type == json_type::string_;
}
inline bool json_value::is_array() const { return _type == json_type::array; }
inline bool json_value::is_object() const { return _type == json_type::object; }
inline bool json_value::is_binary() const { return _type == json_type::binary; }

// conversions
inline json_value::operator int64_t() const {
  return is_unsigned() ? (int64_t)get_unsigned() : get_integer();
}
inline json_value::operator int32_t() const {
  return is_unsigned() ? (int32_t)get_unsigned() : (int32_t)get_integer();
}
inline json_value::operator uint64_t() const {
  return is_integer() ? (uint64_t)get_integer() : get_unsigned();
}
inline json_value::operator uint32_t() const {
  return is_integer() ? (uint32_t)get_integer() : (uint32_t)get_unsigned();
}
inline json_value::operator double() const {
  return is_integer() ? (double)get_integer()
                      : is_unsigned() ? (double)get_unsigned() : get_real();
}
inline json_value::operator float() const {
  return is_integer()
             ? (float)get_integer()
             : is_unsigned() ? (float)get_unsigned() : (float)get_real();
}
inline json_value::operator bool() const { return get_boolean(); }
inline json_value::operator string() const { return get_string(); }
inline json_value::operator json_array() const { return get_array(); }
inline json_value::operator json_object() const { return get_object(); }
inline json_value::operator json_binary() const { return get_binary(); }

// size_t fix
#ifdef __APPLE__
inline json_value::json_value(size_t value)
    : _type{json_type::unsigned_}, _unsigned{value} {}
inline json_value::operator size_t() const {
  return is_integer() ? (uint64_t)get_integer() : get_unsigned();
}
inline json_value& json_value::operator=(size_t value) {
  auto js = json_value{value};
  swap(js);
  return *this;
}
#endif

// access
inline const int64_t& json_value::get_integer() const {
  if (_type != json_type::integer) throw json_type_error{"integer expected"};
  return _integer;
}
inline const uint64_t& json_value::get_unsigned() const {
  if (_type != json_type::unsigned_) throw json_type_error{"unsigned expected"};
  return _unsigned;
}
inline const double& json_value::get_real() const {
  if (_type != json_type::real) throw json_type_error{"real expected"};
  return _real;
}
inline const bool& json_value::get_boolean() const {
  if (_type != json_type::boolean) throw json_type_error{"boolean expected"};
  return _boolean;
}
inline const string& json_value::get_string() const {
  if (_type != json_type::string_) throw json_type_error{"string expected"};
  return *_string_;
}
inline const json_array& json_value::get_array() const {
  if (_type != json_type::array) throw json_type_error{"array expected"};
  return *_array;
}
inline const json_object& json_value::get_object() const {
  if (_type != json_type::object) throw json_type_error{"object expected"};
  return *_object;
}
inline const json_binary& json_value::get_binary() const {
  if (_type != json_type::binary) throw json_type_error{"binary expected"};
  return *_binary;
}
inline int64_t& json_value::get_integer() {
  if (_type != json_type::integer) throw json_type_error{"integer expected"};
  return _integer;
}
inline uint64_t& json_value::get_unsigned() {
  if (_type != json_type::unsigned_) throw json_type_error{"unsigned expected"};
  return _unsigned;
}
inline double& json_value::get_real() {
  if (_type != json_type::real) throw json_type_error{"real expected"};
  return _real;
}
inline bool& json_value::get_boolean() {
  if (_type != json_type::boolean) throw json_type_error{"boolean expected"};
  return _boolean;
}
inline string& json_value::get_string() {
  if (_type != json_type::string_) throw json_type_error{"string expected"};
  return *_string_;
}
inline json_array& json_value::get_array() {
  if (_type != json_type::array) throw json_type_error{"array expected"};
  return *_array;
}
inline json_object& json_value::get_object() {
  if (_type != json_type::object) throw json_type_error{"object expected"};
  return *_object;
}
inline json_binary& json_value::get_binary() {
  if (_type != json_type::binary) throw json_type_error{"binary expected"};
  return *_binary;
}

// structure support
inline bool   json_value::empty() const { return size() == 0; }
inline size_t json_value::size() const {
  switch (_type) {
    case json_type::string_: return get_string().size();
    case json_type::array: return get_array().size();
    case json_type::object: return get_object().size();
    case json_type::binary: return get_binary().size();
    default: throw json_type_error{"bad json type"};
  }
}
inline void json_value::resize(size_t size) {
  switch (_type) {
    case json_type::string_: return get_string().resize(size);
    case json_type::array: return get_array().resize(size);
    case json_type::binary: return get_binary().resize(size, 0);
    default: throw json_type_error{"bad json type"};
  }
}
inline void json_value::reserve(size_t size) {
  switch (_type) {
    case json_type::string_: return get_string().reserve(size);
    case json_type::array: return get_array().reserve(size);
    case json_type::object: return get_object().reserve(size);
    case json_type::binary: return get_binary().reserve(size);
    default: throw json_type_error{"bad json type"};
  }
}

// array support
inline json_value  json_value::array() { return json_value{json_array{}}; }
inline json_value& json_value::operator[](size_t idx) {
  return get_array().at(idx);
}
inline const json_value& json_value::operator[](size_t idx) const {
  return get_array().at(idx);
}
inline json_value& json_value::at(size_t idx) { return get_array().at(idx); }
inline const json_value& json_value::at(size_t idx) const {
  return get_array().at(idx);
}
inline json_value&       json_value::front() { return get_array().front(); }
inline const json_value& json_value::front() const {
  return get_array().front();
}
inline json_value&       json_value::back() { return get_array().back(); }
inline const json_value& json_value::back() const { return get_array().back(); }
inline void              json_value::push_back(const json_value& value) {
  return get_array().push_back(value);
}
inline void json_value::push_back(json_value&& value) {
  return get_array().push_back(std::move(value));
}
template <typename... Args>
inline json_value& json_value::emplace_back(Args&&... args) {
  return get_array().emplace_back(std::forward(args)...);
}
inline json_value*       json_value::begin() { return get_array().data(); }
inline const json_value* json_value::begin() const {
  return get_array().data();
}
inline json_value* json_value::end() {
  return get_array().data() + get_array().size();
}
inline const json_value* json_value::end() const {
  return get_array().data() + get_array().size();
}

// object support
inline json_value  json_value::object() { return json_value{json_object{}}; }
inline json_value& json_value::operator[](const string& key) {
  if (auto ptr = find(key); ptr) {
    return *ptr;
  } else {
    return get_object().emplace_back(key, json_value{}).second;
  }
}
inline json_value& json_value::at(const string& key) {
  if (auto ptr = find(key); ptr) {
    return *ptr;
  } else {
    throw std::out_of_range{"missing key " + key};
  }
}
inline const json_value& json_value::at(const string& key) const {
  if (auto ptr = find(key); ptr) {
    return *ptr;
  } else {
    throw std::out_of_range{"missing key " + key};
  }
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
  return {get_object().data(), get_object().data() + get_object().size()};
}
inline json_value::json_object_cit json_value::items() const {
  return {get_object().data(), get_object().data() + get_object().size()};
}
inline json_value* json_value::find(const string& key) {
  for (auto& [key_, value] : json_value::get_object()) {
    if (key_ == key) return &value;
  }
  return nullptr;
}
inline const json_value* json_value::find(const string& key) const {
  for (auto& [key_, value] : json_value::get_object()) {
    if (key_ == key) return &value;
  }
  return nullptr;
}
inline bool json_value::contains(const string& key) const {
  return find(key) != nullptr;
}

// binary support
inline json_value json_value::binary() { return json_value{json_binary{}}; }

// swap
inline void json_value::swap(json_value& other) {
  std::swap(_type, other._type);
  std::swap(_unsigned, other._unsigned);  // hask to swap bits
}

// destructor
inline json_value::~json_value() {
  switch (_type) {
    case json_type::string_: delete _string_; break;
    case json_type::array: delete _array; break;
    case json_type::object: delete _object; break;
    case json_type::binary: delete _binary; break;
    default: break;
  }
  _type     = json_type::null;
  _unsigned = 0;
}

// set type
inline void json_value::set_type(json_type type) {
  switch (_type) {
    case json_type::string_: delete _string_; break;
    case json_type::array: delete _array; break;
    case json_type::object: delete _object; break;
    case json_type::binary: delete _binary; break;
    default: break;
  }
  _type = type;
  switch (type) {
    case json_type::null: _unsigned = 0; break;
    case json_type::integer: _integer = 0; break;
    case json_type::unsigned_: _unsigned = 0; break;
    case json_type::real: _real = 0; break;
    case json_type::boolean: _boolean = false; break;
    case json_type::string_: _string_ = new string{}; break;
    case json_type::array: _array = new json_array{}; break;
    case json_type::object: _object = new json_object{}; break;
    case json_type::binary: _binary = new json_binary{}; break;
  }
}

// Load/save a json file
bool load_json(const string& filename, json_value& js, string& error);
bool save_json(const string& filename, const json_value& js, string& error);

// Conversion shortcuts
template <typename T>
inline T& from_json(const json_value& js) {
  auto value = T{};
  from_json(js, value);
  return value;
}
template <typename T>
inline json_value to_json(const T& value) {
  auto js = json_value();
  to_json(js, value);
  return js;
}

// Conversion from json to types
inline void from_json(const json_value& js, int64_t& value) {
  value = (int64_t)js;
}
inline void from_json(const json_value& js, int32_t& value) {
  value = (int32_t)js;
}
inline void from_json(const json_value& js, uint64_t& value) {
  value = (uint64_t)js;
}
inline void from_json(const json_value& js, uint32_t& value) {
  value = (uint32_t)js;
}
inline void from_json(const json_value& js, double& value) {
  value = (double)js;
}
inline void from_json(const json_value& js, float& value) { value = (float)js; }
inline void from_json(const json_value& js, bool& value) { value = (bool)js; }
inline void from_json(const json_value& js, string& value) {
  value = (string)js;
}
template <typename T>
inline void from_json(const json_value& js, vector<T>& value) {
  value.clear();
  for (auto& ejs : js) from_json(ejs, value.emplace_back());
}
template <typename T, size_t N>
inline void from_json(const json_value& js, array<T, N>& value) {
  if (js.size() != N) throw std::out_of_range{"invalid size"};
  for (auto idx = 0; idx < N; idx++) from_json(js.at(idx), value.at(idx));
}

// Conversion from types to json
inline void to_json(json_value& js, int64_t value) { js = value; }
inline void to_json(json_value& js, int32_t value) { js = value; }
inline void to_json(json_value& js, uint64_t value) { js = value; }
inline void to_json(json_value& js, uint32_t value) { js = value; }
inline void to_json(json_value& js, double value) { js = value; }
inline void to_json(json_value& js, float value) { js = value; }
inline void to_json(json_value& js, bool value) { js = value; }
inline void to_json(json_value& js, const string& value) { js = value; }
template <typename T>
inline void to_json(json_value& js, const vector<T>& value) {
  js = json_array{};
  for (auto& v : value) to_json(js.emplace_back(), v);
}
template <typename T, size_t N>
inline void to_json(json_value& js, const array<T, N>& value) {
  js = json_array{};
  js.resize(N);
  for (auto idx = 0; idx < N; idx++) to_json(js.at(idx), value.at(idx));
}

// Type
inline json_type get_type(const json_value& js) { return js.type(); }
inline bool      set_type(json_value& js, json_type type) {
  js.set_type(type);
  return true;
}

// Access to internal representation
inline const int64_t* get_integer(const json_value& js) {
  return get_type(js) == json_type::integer ? &js.get_integer() : nullptr;
}
inline const uint64_t* get_unsigned(const json_value& js) {
  return get_type(js) == json_type::unsigned_ ? &js.get_unsigned() : nullptr;
}
inline const double* get_real(const json_value& js) {
  return get_type(js) == json_type::real ? &js.get_real() : nullptr;
}
inline const bool* get_boolean(const json_value& js) {
  return get_type(js) == json_type::boolean ? &js.get_boolean() : nullptr;
}
inline const string* get_string(const json_value& js) {
  return get_type(js) == json_type::string_ ? &js.get_string() : nullptr;
}
inline const json_array* get_array(const json_value& js) {
  return get_type(js) == json_type::array ? &js.get_array() : nullptr;
}
inline const json_object* get_object(const json_value& js) {
  return get_type(js) == json_type::object ? &js.get_object() : nullptr;
}
inline const json_binary* get_binary(const json_value& js) {
  return get_type(js) == json_type::binary ? &js.get_binary() : nullptr;
}

// Access to internal representation
inline int64_t* get_integer(json_value& js) {
  return get_type(js) == json_type::integer ? &js.get_integer() : nullptr;
}
inline uint64_t* get_unsigned(json_value& js) {
  return get_type(js) == json_type::unsigned_ ? &js.get_unsigned() : nullptr;
}
inline double* get_real(json_value& js) {
  return get_type(js) == json_type::real ? &js.get_real() : nullptr;
}
inline bool* get_boolean(json_value& js) {
  return get_type(js) == json_type::boolean ? &js.get_boolean() : nullptr;
}
inline string* get_string(json_value& js) {
  return get_type(js) == json_type::string_ ? &js.get_string() : nullptr;
}
inline json_array* get_array(json_value& js) {
  return get_type(js) == json_type::array ? &js.get_array() : nullptr;
}
inline json_object* get_object(json_value& js) {
  return get_type(js) == json_type::object ? &js.get_object() : nullptr;
}
inline json_binary* get_binary(json_value& js) {
  return get_type(js) == json_type::binary ? &js.get_binary() : nullptr;
}

// Initialization to basic types
inline bool set_integer(json_value& js, int64_t value) {
  if (!set_type(js, json_type::integer)) return false;
  *get_integer(js) = value;
  return true;
}
inline bool set_unsigned(json_value& js, uint64_t value) {
  if (!set_type(js, json_type::unsigned_)) return false;
  *get_unsigned(js) = value;
  return true;
}
inline bool set_real(json_value& js, double value) {
  if (!set_type(js, json_type::real)) return false;
  *get_real(js) = value;
  return true;
}
inline bool set_boolean(json_value& js, bool value) {
  if (!set_type(js, json_type::boolean)) return false;
  *get_boolean(js) = value;
  return true;
}
inline bool set_string(json_value& js, const string& value) {
  if (!set_type(js, json_type::string_)) return false;
  *get_string(js) = value;
  return true;
}
inline bool set_array(json_value& js, const json_array& value) {
  if (!set_type(js, json_type::array)) return false;
  *get_array(js) = value;
  return true;
}
inline bool set_object(json_value& js, const json_object& value) {
  if (!set_type(js, json_type::object)) return false;
  *get_object(js) = value;
  return true;
}
inline bool set_binary(json_value& js, const json_binary& value) {
  if (!set_type(js, json_type::binary)) return false;
  *get_binary(js) = value;
  return true;
}

// Compound type
inline bool empty(const json_value& js) {
  switch (get_type(js)) {
    case json_type::string_: return get_string(js)->empty();
    case json_type::array: return get_array(js)->empty();
    case json_type::object: return get_object(js)->empty();
    case json_type::binary: return get_binary(js)->empty();
    default: return true;
  }
}
inline size_t size(const json_value& js) {
  switch (get_type(js)) {
    case json_type::string_: return get_string(js)->size();
    case json_type::array: return get_array(js)->size();
    case json_type::object: return get_object(js)->size();
    case json_type::binary: return get_binary(js)->size();
    default: return 0;
  }
}

// Array
inline bool set_array(json_value& js) { return set_type(js, json_type::array); }
inline bool resize_array(json_value& js, size_t size) {
  if (auto array = get_array(js); array) {
    array->resize(size);
    return true;
  } else {
    return false;
  }
}
inline json_value* get_element(json_value& js, size_t idx) {
  if (auto array = get_array(js); array) {
    return idx < array->size() ? &array->at(idx) : nullptr;
  } else {
    return nullptr;
  }
}
inline const json_value* get_element(const json_value& js, size_t idx) {
  if (auto array = get_array(js); array) {
    return idx < array->size() ? &array->at(idx) : nullptr;
  } else {
    return nullptr;
  }
}
inline auto iterate_array(json_value& js) {
  struct iterator_wrapper {
    json_value* begin_;
    json_value* end_;
    json_value* begin() { return begin_; }
    json_value* end() { return end_; }
  };
  if (auto array = get_array(js); array) {
    return iterator_wrapper{array->data(), array->data() + array->size()};
  } else {
    return iterator_wrapper{nullptr, nullptr};
  }
}
inline auto iterate_array(const json_value& js) {
  struct iterator_wrapper {
    const json_value* begin_;
    const json_value* end_;
    const json_value* begin() { return begin_; }
    const json_value* end() { return end_; }
  };
  if (auto array = get_array(js); array) {
    return iterator_wrapper{array->data(), array->data() + array->size()};
  } else {
    return iterator_wrapper{nullptr, nullptr};
  }
}

// Conversion from json to values
inline bool get_value(const json_value& js, int64_t& value, string& error,
    const json_value& root) {
  if (auto integer = get_integer(js); integer) {
    value = (int64_t)*integer;
    return true;
  } else if (auto unsigned_ = get_unsigned(js); unsigned_) {
    value = (int64_t)*unsigned_;
    return true;
  } else {
    error = "integer expected";
    return false;
  }
}
inline bool get_value(const json_value& js, int32_t& value, string& error,
    const json_value& root) {
  auto value64 = (int64_t)0;
  if (!get_value(js, value, error)) return false;
  value = (int32_t)value64;
}
inline bool get_value(const json_value& js, uint64_t& value, string& error,
    const json_value& root) {
  if (auto integer = get_integer(js); integer) {
    value = (uint64_t)*integer;
    return true;
  } else if (auto unsigned_ = get_unsigned(js); unsigned_) {
    value = (uint64_t)*unsigned_;
    return true;
  } else {
    error = "integer expected";
    return false;
  }
}
inline bool get_value(const json_value& js, uint32_t& value, string& error,
    const json_value& root) {
  auto value64 = (uint64_t)0;
  if (!get_value(js, value, error)) return false;
  value = (uint32_t)value64;
}
inline bool get_value(const json_value& js, double& value, string& error,
    const json_value& root) {
  if (auto real = get_real(js); real) {
    value = (int64_t)*real;
    return true;
  } else if (auto integer = get_integer(js); integer) {
    value = (int64_t)*integer;
    return true;
  } else if (auto unsigned_ = get_unsigned(js); unsigned_) {
    value = (int64_t)*unsigned_;
    return true;
  } else {
    error = "number expected";
    return false;
  }
}
inline bool get_value(const json_value& js, float& value, string& error) {
  auto value64 = (double)0;
  if (!get_value(js, value, error)) return false;
  value = (float)value64;
}
inline bool get_value(const json_value& js, bool& value, string& error) {
  if (auto boolean = get_boolean(js); boolean) {
    value = *boolean;
    return true;
  } else {
    error = "boolean expected";
    return false;
  }
}
inline bool get_value(const json_value& js, string& value, string& error,
    const json_value& root) {
  if (auto string_ = get_string(js); string_) {
    value = *string_;
    return true;
  } else {
    error = "string expected";
    return false;
  }
}
template <typename T>
inline bool get_value(const json_value& js, vector<T>& value, string& error,
    const json_value& root) {
  if (get_type(js) == json_type::array) {
    value.clear();
    value.reserve(size(js));
    for (auto& ejs : iterate_array(js)) {
      if (!get_value(ejs, value.emplace_back())) return false;
    }
    return true;
  } else {
    error = "array expected";
    return false;
  }
}
template <typename T, size_t N>
inline bool get_value(const json_value& js, array<T, N>& value, string& error,
    const json_value& root) {
  if (get_type(js) == json_type::array && N == size(js)) {
    for (auto idx = (size_t)0; idx < size(js); idx++) {
      if (!get_value(*get_element(js, idx), value.at(idx))) return false;
    }
    return true;
  } else {
    error = get_type(js) == json_type::array ? "array expected"
                                             : "array size mismatched";
    return false;
  }
}

// Conversion to json from values
inline bool set_value(json_value& js, int64_t value, string& error) {
  if (!set_type(js, json_type::integer)) return false;
  *get_integer(js) = value;
  return true;
}
inline bool set_value(json_value& js, int32_t value, string& error) {
  return set_value(js, (int64_t)value, error);
}
inline bool set_value(json_value& js, uint64_t value, string& error) {
  if (!set_type(js, json_type::unsigned_)) return false;
  *get_unsigned(js) = value;
  return true;
}
inline bool set_value(json_value& js, uint32_t value, string& error) {
  return set_value(js, (uint64_t)value, error);
}
inline bool set_value(json_value& js, double value, string& error) {
  if (!set_type(js, json_type::real)) return false;
  *get_real(js) = value;
  return true;
}
inline bool set_value(json_value& js, float value, string& error) {
  return set_value(js, (double)value, error);
}
inline bool set_value(json_value& js, bool value, string& error) {
  if (!set_type(js, json_type::boolean)) return false;
  *get_boolean(js) = value;
  return true;
}
inline bool set_value(json_value& js, const string& value, string& error,
    const json_value& root) {
  if (!set_type(js, json_type::string_)) return false;
  return set_conversion_error("string expected", js, error);
  *get_string(js) = value;
  return true;
}
template <typename T>
inline bool set_value(json_value& js, const vector<T>& value, string& error,
    const json_value& root) {
  if (!set_array(js)) return false;
  if (!resize_array(js, value.size())) return false;
  auto idx = (size_t)0;
  for (auto& v : value) {
    if (!set_value(get_element(js, idx++), v)) return false;
  }
  return true;
}
template <typename T, size_t N>
inline bool set_value(json_value& js, const array<T, N>& value, string& error) {
  if (!set_array(js)) return false;
  if (!resize_array(js, value.size())) return false;
  auto idx = (size_t)0;
  for (auto& v : value) {
    if (!set_value(get_element(js, idx++), v, error)) return false;
  }
  return true;
}

// Helpers for user-defined types
inline bool check_array(const json_value& js, string& error) {
  if (auto array = get_array(js); array) {
    return true;
  } else {
    error = "array expected";
    return false;
  }
}
inline bool check_array(const json_value& js, size_t size_, string& error) {
  if (auto array = get_array(js); array) {
    if (size_ == array->size()) {
      return true;
    } else {
      error = "mismatchd array size";
      return false;
    }
  } else {
    error = "array expected";
    return false;
  }
}
inline bool check_object(const json_value& js, string& error) {
  if (auto array = get_array(js); array) {
    return true;
  } else {
    error = "object expected";
    return false;
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// FILE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Safe wrapper for FILE stream
struct file_stream {
  // file parameters
  string filename = "";
  FILE*  fs       = nullptr;
  bool   owned    = false;

  // move-only type
  file_stream(const file_stream&) = delete;
  file_stream& operator=(const file_stream&) = delete;
  ~file_stream();

  // operator bool to check for error
  explicit operator bool() const { return fs != nullptr; }
};

// Open a file
file_stream open_file(const string& filename, const string& mode);

// Close a file
void close_file(file_stream& fs);

// Read a line of text
bool read_line(file_stream& fs, char* buffer, size_t size);

// Read a line of text
template <size_t N>
inline bool read_line(file_stream& fs, array<char, N>& buffer) {
  return read_line(fs, buffer.data(), buffer.size());
}

// Write text to a file
bool write_text(file_stream& fs, const string& str);

// Read data from a file
bool read_data(file_stream& fs, void* buffer, size_t count);

// Write data from a file
bool write_data(file_stream& fs, const void* buffer, size_t count);

// Read data from a file
template <typename T>
inline bool read_value(file_stream& fs, T& buffer) {
  return read_data(fs, &buffer, sizeof(T));
}

// Write data from a file
template <typename T>
inline bool write_value(file_stream& fs, const T& buffer) {
  return write_data(fs, &buffer, sizeof(T));
}

// Read data from a file
template <typename T>
inline bool read_values(file_stream& fs, T* buffer, size_t count) {
  return read_data(fs, buffer, sizeof(T) * count);
}

// Write data from a file
template <typename T>
inline bool write_values(file_stream& fs, const T* buffer, size_t count) {
  return write_data(fs, buffer, sizeof(T) * count);
}

template <typename T>
inline T swap_endian(T value) {
  // https://stackoverflow.com/questions/105252/how-do-i-convert-between-big-endian-and-little-endian-values-in-c
  static_assert(sizeof(char) == 1, "sizeof(char) == 1");
  union {
    T             value;
    unsigned char bytes[sizeof(T)];
  } source, dest;
  source.value = value;
  for (auto k = (size_t)0; k < sizeof(T); k++)
    dest.bytes[k] = source.bytes[sizeof(T) - k - 1];
  return dest.value;
}

template <typename T>
inline bool read_value(file_stream& fs, T& value, bool big_endian) {
  if (!read_value(fs, value)) return false;
  if (big_endian) value = swap_endian(value);
  return true;
}

template <typename T>
inline bool write_value(file_stream& fs, const T& value_, bool big_endian) {
  auto value = big_endian ? swap_endian(value_) : value_;
  return write_value(fs, value);
}

// Opens a file with a utf8 file name
FILE* fopen_utf8(const char* filename, const char* mode);

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF FORMATTING
// -----------------------------------------------------------------------------
namespace yocto {

// Formats values to string
inline void format_value(string& str, const string& value) { str += value; }
inline void format_value(string& str, int8_t value) {
  str += std::to_string((int32_t)value);
}
inline void format_value(string& str, int16_t value) {
  str += std::to_string(value);
}
inline void format_value(string& str, int32_t value) {
  str += std::to_string(value);
}
inline void format_value(string& str, int64_t value) {
  str += std::to_string(value);
}
inline void format_value(string& str, uint8_t value) {
  str += std::to_string((uint32_t)value);
}
inline void format_value(string& str, uint16_t value) {
  str += std::to_string(value);
}
inline void format_value(string& str, uint32_t value) {
  str += std::to_string(value);
}
inline void format_value(string& str, uint64_t value) {
  str += std::to_string(value);
}
inline void format_value(string& str, float value) {
  auto buf = array<char, 256>{};
  snprintf(buf.data(), buf.size(), "%g", value);
  str += buf.data();
}
inline void format_value(string& str, double value) {
  auto buf = array<char, 256>{};
  snprintf(buf.data(), buf.size(), "%g", value);
  str += buf.data();
}

// Foramt to file
inline void format_values(string& str, const string& fmt) {
  auto pos = fmt.find("{}");
  if (pos != string::npos) throw std::invalid_argument("bad format string");
  str += fmt;
}
template <typename Arg, typename... Args>
inline void format_values(
    string& str, const string& fmt, const Arg& arg, const Args&... args) {
  auto pos = fmt.find("{}");
  if (pos == string::npos) throw std::invalid_argument("bad format string");
  str += fmt.substr(0, pos);
  format_value(str, arg);
  format_values(str, fmt.substr(pos + 2), args...);
}

template <typename... Args>
inline bool format_values(
    file_stream& fs, const string& fmt, const Args&... args) {
  auto str = ""s;
  format_values(str, fmt, args...);
  return write_text(fs, str);
}
template <typename T>
inline bool format_value(file_stream& fs, const T& value) {
  auto str = ""s;
  format_value(str, value);
  return write_text(fs, str);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF COMMAND-LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// Command line value type
enum struct cli_type { integer, uinteger, number, boolean, string };
// Command line value
struct cli_value {
  cli_type type     = cli_type::integer;
  int64_t  integer  = 0;
  uint64_t uinteger = 0;
  double   number   = 0;
  string   text     = "";
};
// Command line option. All data should be considered private.
struct cli_option {
  string                                          name    = "";
  cli_type                                        type    = cli_type::string;
  bool                                            req     = false;
  int                                             nargs   = 0;
  string                                          usage   = "";
  vector<cli_value>                               value   = {};
  vector<cli_value>                               def     = {};
  vector<string>                                  choices = {};
  bool                                            set     = false;
  function<void(const vector<cli_value>& values)> set_reference = {};
};
// Command line parser. All data should be considered private.
struct cli_state {
  string             name            = "";
  string             usage           = "";
  vector<cli_option> options         = {};
  string             usage_options   = "";
  string             usage_arguments = "";
  bool               help            = false;
};

template <typename T>
inline cli_type get_cli_type() {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  if constexpr (std::is_same_v<T, string>) {
    return cli_type::string;
  } else if constexpr (std::is_same_v<T, bool>) {
    return cli_type::boolean;
  } else if constexpr (std::is_enum_v<T>) {
    return cli_type::integer;
  } else if constexpr (std::is_integral_v<T> && !std::is_unsigned_v<T>) {
    return cli_type::integer;
  } else if constexpr (std::is_integral_v<T> && std::is_unsigned_v<T>) {
    return cli_type::uinteger;
  } else if constexpr (std::is_floating_point_v<T>) {
    return cli_type::number;
  } else {
    // probably should be an error
    return cli_type::string;
  }
}

template <typename T>
inline void set_value(cli_value& cvalue, const T& value) {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  cvalue.type = get_cli_type<T>();
  if constexpr (std::is_same_v<T, string>) {
    cvalue.text = value;
  } else if constexpr (std::is_same_v<T, bool>) {
    cvalue.integer = value ? 1 : 0;
  } else if constexpr (std::is_enum_v<T>) {
    cvalue.integer = (int64_t)value;
  } else if constexpr (std::is_integral_v<T> && !std::is_unsigned_v<T>) {
    cvalue.integer = value;
  } else if constexpr (std::is_integral_v<T> && std::is_unsigned_v<T>) {
    cvalue.uinteger = value;
  } else if constexpr (std::is_floating_point_v<T>) {
    cvalue.number = value;
  } else {
    // probably should be an error
    // pass
  }
}

template <typename T>
inline bool get_value(const cli_value& cvalue, T& value) {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  if constexpr (std::is_same_v<T, string>) {
    if (cvalue.type != cli_type::string) return false;
    value = cvalue.text;
    return true;
  } else if constexpr (std::is_same_v<T, bool>) {
    if (cvalue.type != cli_type::boolean) return false;
    value = cvalue.integer != 0;
    return true;
  } else if constexpr (std::is_enum_v<T>) {
    if (cvalue.type != cli_type::integer) return false;
    value = (T)cvalue.integer;
    return true;
  } else if constexpr (std::is_integral_v<T> && !std::is_unsigned_v<T>) {
    if (cvalue.type != cli_type::integer) return false;
    value = (T)cvalue.integer;
    return true;
  } else if constexpr (std::is_integral_v<T> && std::is_unsigned_v<T>) {
    if (cvalue.type != cli_type::uinteger) return false;
    value = (T)cvalue.uinteger;
    return true;
  } else if constexpr (std::is_floating_point_v<T>) {
    if (cvalue.type != cli_type::number) return false;
    value = (T)cvalue.number;
    return true;
  } else {
    return false;
  }
}

template <typename T>
inline void add_option(cli_state& cli, const string& name, T& value,
    const string& usage, bool req) {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  auto def = vector<cli_value>{};
  set_value(def.emplace_back(), value);
  auto& option         = cli.options.emplace_back();
  option.name          = name;
  option.type          = get_cli_type<T>();
  option.req           = req;
  option.nargs         = !std::is_same_v<T, bool> ? 1 : 0;
  option.usage         = usage;
  option.value         = def;
  option.def           = def;
  option.choices       = {};
  option.set_reference = [&value](const vector<cli_value>& cvalues) -> bool {
    if (cvalues.size() != 1) throw std::out_of_range{"invalid number of args"};
    return get_value(cvalues.front(), value);
  };
}

template <typename T>
inline void add_option(cli_state& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, bool req) {
  static_assert(
      std::is_same_v<T, string> || std::is_integral_v<T> || std::is_enum_v<T>,
      "unsupported type");
  auto def = vector<cli_value>{};
  set_value(def.emplace_back(), value);
  auto& option         = cli.options.emplace_back();
  option.name          = name;
  option.type          = get_cli_type<T>();
  option.req           = req;
  option.nargs         = 1;
  option.usage         = usage;
  option.value         = def;
  option.def           = def;
  option.choices       = choices;
  option.set_reference = [&value](const vector<cli_value>& cvalues) -> bool {
    if (cvalues.size() != 1) throw std::out_of_range{"invalid number of args"};
    return get_value(cvalues.front(), value);
  };
}

template <typename T>
inline void add_option(cli_state& cli, const string& name, vector<T>& values,
    const string& usage, bool req) {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  auto def = vector<cli_value>{};
  for (auto& value : values) set_value(def.emplace_back(), value);
  auto& option         = cli.options.emplace_back();
  option.name          = name;
  option.type          = get_cli_type<T>();
  option.req           = req;
  option.nargs         = -1;
  option.usage         = usage;
  option.value         = def;
  option.def           = def;
  option.choices       = {};
  option.set_reference = [&values](const vector<cli_value>& cvalues) -> bool {
    values.clear();
    for (auto& cvalue : cvalues) {
      if (!get_value(cvalue, values.emplace_back())) return false;
    }
    return true;
  };
}

}  // namespace yocto

#endif
