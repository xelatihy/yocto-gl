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

// Load/save a json file
bool load_json(const string& filename, json_value& js, string& error);
bool save_json(const string& filename, const json_value& js, string& error);

// Conversion shortcuts
template <typename T>
inline T from_json(const json_value& js);
template <typename T>
inline json_value to_json(const T& value);

// Conversion from json to values
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

// Conversion to json from values
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
bool load_json(const string& filename, json_tree& js, string& error);
bool save_json(const string& filename, const json_tree& js, string& error);

// Get view from value
inline json_view  get_root(json_tree& js);
inline json_cview get_croot(json_tree& js);

// Error handling
inline void set_error(json_tree& js, string_view error);
inline void clear_error(json_tree& js);

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
inline bool   is_valid(json_cview js);
inline bool   is_valid(json_view js);
inline string get_error(json_cview js);
inline string get_error(json_view js);
inline string compute_path(json_cview js);
inline bool   set_error(json_cview js, string_view error);

// Type
inline json_type get_type(json_cview js);
inline bool      set_type(json_view js, json_type type);
inline bool      is_null(json_cview js);
inline bool      is_integer(json_cview js);
inline bool      is_unsigned(json_cview js);
inline bool      is_real(json_cview js);
inline bool      is_integral(json_cview js);
inline bool      is_number(json_cview js);
inline bool      is_boolean(json_cview js);
inline bool      is_string(json_cview js);
inline bool      is_array(json_cview js);
inline bool      is_object(json_cview js);
inline bool      is_binary(json_cview js);

// Initialization to basic types
inline bool set_null(json_view js);
inline bool set_integer(json_view js, int64_t value);
inline bool set_unsigned(json_view js, uint64_t value);
inline bool set_real(json_view js, double value);
inline bool set_boolean(json_view js, bool value);
inline bool set_string(json_view js, const string& value);
inline bool set_integral(json_view js, int64_t value);
inline bool set_integral(json_view js, int32_t value);
inline bool set_integral(json_view js, uint64_t value);
inline bool set_integral(json_view js, uint32_t value);
inline bool set_number(json_view js, double value);
inline bool set_number(json_view js, float value);

// Get basic values
inline bool get_integer(json_cview js, int64_t& value);
inline bool get_unsigned(json_cview js, uint64_t& value);
inline bool get_real(json_cview js, double& value);
inline bool get_boolean(json_cview js, bool& value);
inline bool get_string(json_cview js, string& value);
inline bool get_integral(json_cview js, int64_t& value);
inline bool get_integral(json_cview js, uint64_t& value);
inline bool get_number(json_cview js, double& value);
inline bool get_integral(json_cview js, int32_t& value);
inline bool get_integral(json_cview js, uint32_t& value);
inline bool get_number(json_cview js, float& value);

// Get basic values - ignore errors if present
inline int64_t  get_integer(json_cview js);
inline uint64_t get_unsigned(json_cview js);
inline double   get_real(json_cview js);
inline bool     get_boolean(json_cview js);
inline string   get_string(json_cview js);
inline int64_t  get_integral(json_cview js);
inline uint64_t get_uintegral(json_cview js);
inline double   get_number(json_cview js);

// Compound type
inline bool   is_empty(json_cview js);
inline size_t get_size(json_cview js);

// Array
inline bool       set_array(json_view js);
inline bool       set_array(json_view js, size_t size);
inline bool       array_size(json_cview js, size_t& size);
inline bool       has_element(json_view js, size_t idx);
inline bool       has_element(json_cview js, size_t idx);
inline json_view  get_element(json_view js, size_t idx);
inline json_cview get_element(json_cview js, size_t idx);
inline json_view  append_element(json_view js);
inline auto       iterate_array(json_view js);
inline auto       iterate_array(json_cview js);

// Object
inline bool       set_object(json_view js);
inline bool       object_size(json_cview js, size_t& size);
inline bool       has_element(json_view js, string_view key);
inline bool       has_element(json_cview js, string_view key);
inline json_view  get_element(json_view js, string_view key);
inline json_cview get_element(json_cview js, string_view key);
inline json_view  insert_element(json_view js, string_view key);
inline auto       iterate_object(json_view js);
inline auto       iterate_object(json_cview js);

// Binary
inline bool set_binary(json_view js, const json_binary& value);
inline bool get_binary(json_cview js, json_binary& value);

// Get the path of a json view
inline string compute_path(json_cview js);

// Conversion from json to values
template <typename T>
inline bool get_value(json_cview js, T& value);

// Conversion from json to values
inline bool get_value(json_cview js, int64_t& value);
inline bool get_value(json_cview js, int32_t& value);
inline bool get_value(json_cview js, uint64_t& value);
inline bool get_value(json_cview js, uint32_t& value);
inline bool get_value(json_cview js, double& value);
inline bool get_value(json_cview js, float& value);
inline bool get_value(json_cview js, bool& value);
inline bool get_value(json_cview js, string& value);
template <typename T>
inline bool get_value(json_cview js, vector<T>& value);
template <typename T, size_t N>
inline bool get_value(json_cview js, array<T, N>& value);

// Get value at a key or index
template <typename T>
inline bool get_value_at(json_cview js, string_view key, T& value);
template <typename T>
inline bool get_value_at(json_cview js, size_t idx, T& value);

// Get value at a key or nothing is key is not preesent
template <typename T>
inline bool get_value_if(json_cview js, string_view key, T& value);

// Conversion to json from values
template <typename T>
inline bool set_value(json_view js, const T& value);

// Conversion to json from values
inline bool set_value(json_view js, int64_t value);
inline bool set_value(json_view js, int32_t value);
inline bool set_value(json_view js, uint64_t value);
inline bool set_value(json_view js, uint32_t value);
inline bool set_value(json_view js, double value);
inline bool set_value(json_view js, float value);
inline bool set_value(json_view js, bool value);
inline bool set_value(json_view js, const string& value);
inline bool set_value(json_view js, const char* value);
template <typename T>
inline bool set_value(json_view js, const vector<T>& value);
template <typename T, size_t N>
inline bool set_value(json_view js, const array<T, N>& value);

// Helpers for user-defined types
inline bool check_array(json_cview js);
inline bool check_array(json_cview js, size_t size_);
inline bool check_object(json_cview js);

// Helpers for user-defined types
inline bool set_array(json_view js);
template <typename T>
inline bool set_value_at(json_view js, size_t idx, const T& value);
template <typename T>
inline bool      append_value(json_view js, const T& value);
inline json_view append_array(json_view js);
inline json_view append_object(json_view js);

// Helpers for user-defined types
inline bool set_object(json_view js);
template <typename T>
inline bool set_value_at(json_view js, string_view key, const T& value);
template <typename T>
inline bool insert_value(json_view js, string_view key, const T& value);
template <typename T>
inline bool insert_value_if(
    json_view js, string_view key, const T& value, const T& default_);
inline json_view insert_array(json_view js, string_view key);
inline json_view insert_object(json_view js, string_view key);

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
inline json_value::json_value(const char* value)
    : _type{json_type::string_}, _string{new string{value}} {}
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
  if (_type != json_type::integer) throw json_error{"integer expected"};
  return _integer;
}
inline const uint64_t& json_value::get_unsigned() const {
  if (_type != json_type::unsigned_) throw json_error{"unsigned expected"};
  return _unsigned;
}
inline const double& json_value::get_real() const {
  if (_type != json_type::real) throw json_error{"real expected"};
  return _real;
}
inline const bool& json_value::get_boolean() const {
  if (_type != json_type::boolean) throw json_error{"boolean expected"};
  return _boolean;
}
inline const string& json_value::get_string() const {
  if (_type != json_type::string_) throw json_error{"string expected"};
  return *_string;
}
inline const json_array& json_value::get_array() const {
  if (_type != json_type::array) throw json_error{"array expected"};
  return *_array;
}
inline const json_object& json_value::get_object() const {
  if (_type != json_type::object) throw json_error{"object expected"};
  return *_object;
}
inline const json_binary& json_value::get_binary() const {
  if (_type != json_type::binary) throw json_error{"binary expected"};
  return *_binary;
}
inline int64_t& json_value::get_integer() {
  if (_type != json_type::integer) throw json_error{"integer expected"};
  return _integer;
}
inline uint64_t& json_value::get_unsigned() {
  if (_type != json_type::unsigned_) throw json_error{"unsigned expected"};
  return _unsigned;
}
inline double& json_value::get_real() {
  if (_type != json_type::real) throw json_error{"real expected"};
  return _real;
}
inline bool& json_value::get_boolean() {
  if (_type != json_type::boolean) throw json_error{"boolean expected"};
  return _boolean;
}
inline string& json_value::get_string() {
  if (_type != json_type::string_) throw json_error{"string expected"};
  return *_string;
}
inline json_array& json_value::get_array() {
  if (_type != json_type::array) throw json_error{"array expected"};
  return *_array;
}
inline json_object& json_value::get_object() {
  if (_type != json_type::object) throw json_error{"object expected"};
  return *_object;
}
inline json_binary& json_value::get_binary() {
  if (_type != json_type::binary) throw json_error{"binary expected"};
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
    default: throw json_error{"bad json type"};
  }
}
inline void json_value::resize(size_t size) {
  switch (_type) {
    case json_type::string_: return get_string().resize(size);
    case json_type::array: return get_array().resize(size);
    case json_type::binary: return get_binary().resize(size, 0);
    default: throw json_error{"bad json type"};
  }
}
inline void json_value::reserve(size_t size) {
  switch (_type) {
    case json_type::string_: return get_string().reserve(size);
    case json_type::array: return get_array().reserve(size);
    case json_type::object: return get_object().reserve(size);
    case json_type::binary: return get_binary().reserve(size);
    default: throw json_error{"bad json type"};
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
  for (auto& [key_, value] : get_object()) {
    if (key_ == key) return &value;
  }
  return nullptr;
}
inline const json_value* json_value::find(const string& key) const {
  for (auto& [key_, value] : get_object()) {
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
    case json_type::string_: delete _string; break;
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
    case json_type::string_: delete _string; break;
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
    case json_type::string_: _string = new string{}; break;
    case json_type::array: _array = new json_array{}; break;
    case json_type::object: _object = new json_object{}; break;
    case json_type::binary: _binary = new json_binary{}; break;
  }
}

// Conversion shortcuts
template <typename T>
inline T from_json(const json_value& js) {
  auto value = T{};
  from_json(js, value);
  return value;
}
template <typename T>
inline json_value to_json(const T& value) {
  auto js = json_value{};
  to_json(js, value);
  return js;
}

// Conversion from json to values
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
  value.reserve(js.size());
  for (auto& ejs : js) from_json(ejs, value.emplace_back());
}
template <typename T, size_t N>
inline void from_json(const json_value& js, array<T, N>& value) {
  if (js.size() != value.size()) throw json_error{"wrong array size"};
  for (auto idx = (size_t)0; idx < value.size(); idx++)
    from_json(js.at(idx), value.at(idx));
}

// Conversion to json from values
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
  js.resize(value.size());
  for (auto idx = (size_t)0; idx < value.size(); idx++)
    to_json(js.at(idx), value.at(idx));
}

// Get view from value
inline json_view  get_root(json_tree& js) { return {&js, 0}; }
inline json_cview get_croot(json_tree& js) { return {&js, 0}; }
inline bool       set_error(json_cview js, string_view error) {
  if (!is_valid(js)) return false;
  set_error(*js.root, string{error} + " at " + compute_path(js));
  return false;
}
inline json_view set_error_view(json_cview js, string_view error) {
  if (!is_valid(js)) return {js.root};
  set_error(*js.root, string{error} + " at " + compute_path(js));
  return {js.root};
}

// Error handling
inline void set_error(json_tree& js, string_view error) {
  if (!js.valid) return;
  js.valid = false;
  js.error = string{error};
}
inline void clear_error(json_tree& js) {
  js.valid = true;
  js.error = "";
}

// Helpers
inline json_type& _get_type(json_view js) {
  if (!is_valid(js)) throw std::invalid_argument{"bad json"};
  return js.root->types[js.index];
}
inline const json_type& _get_type(json_cview js) {
  if (!is_valid(js)) throw std::invalid_argument{"bad json"};
  return js.root->types[js.index];
}
inline json_tree::json_value& _get_value(json_view js) {
  if (!is_valid(js)) throw std::invalid_argument{"bad json"};
  return js.root->values[js.index];
}
inline const json_tree::json_value& _get_value(json_cview js) {
  if (!is_valid(js)) throw std::invalid_argument{"bad json"};
  return js.root->values[js.index];
}
inline uint32_t _get_capacity(uint32_t length) {
  if (length == 0) return 0;
  if (length <= 4) return 4;
  // TODO(fabio): faster pow2
  auto capacity = (uint32_t)4;
  while (capacity < length) capacity *= 2;
  return capacity;
}
inline string_view _get_key(json_cview js) {
  if (!is_valid(js)) throw std::invalid_argument{"bad tree"};
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::string_) throw std::invalid_argument{"bad key"};
  return {js.root->keys.data() + jsv._string.start, jsv._string.length};
}
inline void _find_path(json_view js, vector<json_view>& path);

// Error check
inline bool is_valid(json_cview js) {
  return js.root != nullptr && js.root->valid && js.index != (uint32_t)-1;
}
inline bool is_valid(json_view js) {
  return js.root != nullptr && js.root->valid && js.index != (uint32_t)-1;
}
inline string get_error(json_cview js) {
  if (js.root == nullptr) return "bad root";
  if (js.root->valid) return "";
  return js.root->error;
}
inline string get_error(json_view js) {
  if (js.root == nullptr) return "bad root";
  if (js.root->valid) return "";
  return js.root->error;
}

// Type
inline json_type get_type(json_cview js) {
  if (!is_valid(js)) return json_type::null;
  return _get_type(js);
}
inline bool is_null(json_cview js) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  return jst == json_type::null;
}
inline bool is_integer(json_cview js) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  return jst == json_type::integer;
}
inline bool is_unsigned(json_cview js) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  return jst == json_type::unsigned_;
}
inline bool is_real(json_cview js) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  return jst == json_type::real;
}
inline bool is_integral(json_cview js) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  return jst == json_type::integer || jst == json_type::unsigned_;
}
inline bool is_number(json_cview js) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  return jst == json_type::integer || jst == json_type::unsigned_ ||
         jst == json_type::real;
}
inline bool is_boolean(json_cview js) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  return jst == json_type::boolean;
}
inline bool is_string(json_cview js) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  return jst == json_type::string_;
}
inline bool is_array(json_cview js) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  return jst == json_type::array;
}
inline bool is_object(json_cview js) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  return jst == json_type::object;
}
inline bool is_binary(json_cview js) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  return jst == json_type::binary;
}

// Initialization to basic types
inline bool set_null(json_view js) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  jst       = json_type::null;
  return true;
}
inline bool set_integer(json_view js, int64_t value) {
  if (!is_valid(js)) return false;
  auto& jst    = _get_type(js);
  auto& jsv    = _get_value(js);
  jst          = json_type::integer;
  jsv._integer = value;
  return true;
}
inline bool set_unsigned(json_view js, uint64_t value) {
  if (!is_valid(js)) return false;
  auto& jst     = _get_type(js);
  auto& jsv     = _get_value(js);
  jst           = json_type::unsigned_;
  jsv._unsigned = value;
  return true;
}
inline bool set_real(json_view js, double value) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  jst       = json_type::real;
  jsv._real = value;
  return true;
}
inline bool set_boolean(json_view js, bool value) {
  if (!is_valid(js)) return false;
  auto& jst    = _get_type(js);
  auto& jsv    = _get_value(js);
  jst          = json_type::boolean;
  jsv._boolean = value;
  return true;
}
inline bool set_string(json_view js, const string& value) {
  if (!is_valid(js)) return false;
  auto& jst          = _get_type(js);
  auto& jsv          = _get_value(js);
  jst                = json_type::string_;
  jsv._string.start  = (uint32_t)js.root->strings.size();
  jsv._string.length = (uint32_t)value.size();
  js.root->strings.insert(js.root->strings.end(), value.begin(), value.end());
  js.root->strings.push_back(0);
  return true;
}
inline bool set_integral(json_view js, int64_t value) {
  return set_integer(js, value);
}
inline bool set_integral(json_view js, int32_t value) {
  return set_integer(js, value);
}
inline bool set_integral(json_view js, uint64_t value) {
  return set_unsigned(js, value);
}
inline bool set_integral(json_view js, uint32_t value) {
  return set_unsigned(js, value);
}
inline bool set_number(json_view js, double value) {
  return set_real(js, value);
}
inline bool set_number(json_view js, float value) {
  return set_real(js, value);
}

// Get basic values
inline bool get_integer(json_cview js, int64_t& value) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::integer) return set_error(js, "integer expected");
  value = jsv._integer;
  return true;
}
inline bool get_unsigned(json_cview js, uint64_t& value) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::unsigned_) return set_error(js, "unsigned expected");
  value = jsv._unsigned;
  return true;
}
inline bool get_real(json_cview js, double& value) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::real) return set_error(js, "real expected");
  value = jsv._real;
  return true;
}
inline bool get_boolean(json_cview js, bool& value) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::boolean) return set_error(js, "boolean expected");
  value = jsv._boolean;
  return true;
}
inline bool get_string(json_cview js, string& value) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::string_) return set_error(js, "string expected");
  value = string{js.root->strings.data() + jsv._string.start,
      js.root->strings.data() + jsv._string.start + jsv._string.length};
  return true;
}
inline bool get_integral(json_cview js, int64_t& value) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::integer && jst != json_type::unsigned_)
    return set_error(js, "integer expected");
  value = (jst == json_type::integer) ? (int64_t)jsv._integer
                                      : (int64_t)jsv._unsigned;
  return true;
}
inline bool get_integral(json_cview js, uint64_t& value) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::integer && jst != json_type::unsigned_)
    return set_error(js, "integer expected");
  value = (jst == json_type::integer) ? (uint64_t)jsv._integer
                                      : (uint64_t)jsv._unsigned;
  return true;
}
inline bool get_number(json_cview js, double& value) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::real && jst != json_type::integer &&
      jst != json_type::unsigned_)
    return set_error(js, "number expected");
  value = (jst == json_type::real)
              ? (double)jsv._real
              : (jst == json_type::integer) ? (double)jsv._integer
                                            : (double)jsv._unsigned;
  return true;
}
inline bool get_integral(json_cview js, int32_t& value) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::integer && jst != json_type::unsigned_)
    return set_error(js, "integer expected");
  value = (jst == json_type::integer) ? (int32_t)jsv._integer
                                      : (int32_t)jsv._unsigned;
  return true;
}
inline bool get_integral(json_cview js, uint32_t& value) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::integer && jst != json_type::unsigned_)
    return set_error(js, "integer expected");
  value = (jst == json_type::integer) ? (uint32_t)jsv._integer
                                      : (uint32_t)jsv._unsigned;
  return true;
}
inline bool get_number(json_cview js, float& value) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::real && jst != json_type::integer &&
      jst != json_type::unsigned_)
    return set_error(js, "number expected");
  value = (jst == json_type::real)
              ? (float)jsv._real
              : (jst == json_type::integer) ? (float)jsv._integer
                                            : (float)jsv._unsigned;
  return true;
}

// Get basic values
inline int64_t get_integer(json_cview js) {
  auto value = (int64_t)0;
  return get_integer(js, value) ? value : 0;
}
inline uint64_t get_unsigned(json_cview js) {
  auto value = (uint64_t)0;
  return get_unsigned(js, value) ? value : 0;
}
inline double get_real(json_cview js) {
  auto value = (double)0;
  return get_real(js, value) ? value : 0;
}
inline bool get_boolean(json_cview js) {
  auto value = false;
  return get_boolean(js, value) ? value : false;
}
inline string get_string(json_cview js) {
  auto value = string{};
  return get_string(js, value) ? value : string{};
}
inline int64_t get_integral(json_cview js) {
  auto value = (int64_t)0;
  return get_integral(js, value) ? value : 0;
}
inline double get_number(json_cview js) {
  auto value = (double)0;
  return get_number(js, value) ? value : 0;
}

// Compound type
inline bool   is_empty(json_cview js);
inline size_t get_size(json_cview js);

// Array iteeration
inline auto iterate_array(json_view js) {
  struct iterator {
    json_view js;
    bool      operator!=(const iterator& other) {
      return is_valid(js) && js.index != other.js.index;
    }
    iterator& operator++() {
      if (!is_valid(js)) return *this;
      auto& jst = _get_type(js);
      auto& jsv = _get_value(js);
      js.index += 1;
      if (jst == json_type::array) js.index += jsv._array.skip;
      if (jst == json_type::object) js.index += jsv._object.skip;
      return *this;
    }
    json_view operator*() const { return js; }
  };
  struct iterator_wrapper {
    json_view begin_;
    json_view end_;
    iterator  begin() { return {begin_}; }
    iterator  end() { return {end_}; }
  };
  if (!is_valid(js)) return iterator_wrapper{{js.root}, {js.root}};
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::array) {
    set_error_view(js, "array expected");
    return iterator_wrapper{{js.root}, {js.root}};
  }
  return iterator_wrapper{
      {js.root, js.index + 1}, {js.root, js.index + 1 + jsv._array.skip}};
}
inline auto iterate_array(json_cview js) {
  struct iterator {
    json_cview js;
    bool       operator!=(const iterator& other) {
      return is_valid(js) && js.index != other.js.index;
    }
    iterator& operator++() {
      if (!is_valid(js)) return *this;
      auto& jst = _get_type(js);
      auto& jsv = _get_value(js);
      js.index += 1;
      if (jst == json_type::array) js.index += jsv._array.skip;
      if (jst == json_type::object) js.index += jsv._object.skip;
      return *this;
    }
    json_cview operator*() const { return js; }
  };
  struct iterator_wrapper {
    json_cview begin_;
    json_cview end_;
    iterator   begin() { return {begin_}; }
    iterator   end() { return {end_}; }
  };
  if (!is_valid(js)) return iterator_wrapper{{js.root}, {js.root}};
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::array) {
    set_error_view(js, "array expected");
    return iterator_wrapper{{js.root}, {js.root}};
  }
  return iterator_wrapper{
      {js.root, js.index + 1},
      {js.root, js.index + 1 + jsv._array.skip},
  };
}

// Array
inline bool set_array(json_view js) {
  if (!is_valid(js)) return false;
  if (js.index != js.root->values.size() - 1)
    throw std::out_of_range{"can only add at the end"};
  auto& jst  = _get_type(js);
  auto& jsv  = _get_value(js);
  jst        = json_type::array;
  jsv._array = {0, 0};
  return true;
}
inline bool array_size(json_cview js, size_t& size) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::array) return set_error(js, "array expected");
  size = (size_t)jsv._array.length;
  return true;
}
inline json_view get_element(json_view js, size_t idx) {
  if (!is_valid(js)) return {js.root};
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::array) return set_error_view(js, "array expected");
  if (idx >= jsv._array.length)
    return set_error_view(js, "index out of bounds");
  if (jsv._array.length == jsv._array.skip) {
    return {js.root, js.index + 1 + (uint32_t)idx};
  } else {
    auto count = 0;
    for (auto ejs : iterate_array(js)) {
      if (count++ == idx) return ejs;
    }
    return set_error_view(js, "index out of bounds");
  }
}
inline json_cview get_element(json_cview js, size_t idx) {
  if (!is_valid(js)) return {js.root};
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::array) return set_error_view(js, "array expected");
  if (idx >= jsv._array.length)
    return set_error_view(js, "index out of bounds");
  if (jsv._array.length == jsv._array.skip) {
    return {js.root, js.index + 1 + (uint32_t)idx};
  } else {
    auto count = 0;
    for (auto ejs : iterate_array(js)) {
      if (count++ == idx) return ejs;
    }
    return set_error_view(js, "index out of bounds");
  }
}
inline json_view append_element(json_view js) {
  if (!is_valid(js)) return {js.root};
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::array) return set_error_view(js, "array expected");
  if (js.index + 1 + jsv._array.skip != js.root->values.size())
    throw std::out_of_range{"can only add at the end"};
  jsv._array.length += 1;
  auto index = (uint32_t)js.root->values.size();
  js.root->types.emplace_back(json_type::null);
  js.root->values.emplace_back();
  auto stack = vector<json_view>{};
  _find_path(js, stack);
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
  return {js.root, index};
}

// Object iteration
inline auto iterate_object(json_view js) {
  struct iterator {
    json_view js;
    bool      operator!=(const iterator& other) {
      return is_valid(js) && js.index != other.js.index;
    }
    iterator& operator++() {
      if (!is_valid(js)) return *this;
      auto& jst = _get_type(json_view{js.root, js.index + 1});
      auto& jsv = _get_value(json_view{js.root, js.index + 1});
      js.index += 2;
      if (jst == json_type::array) js.index += jsv._array.skip;
      if (jst == json_type::object) js.index += jsv._object.skip;
      return *this;
    }
    pair<string_view, json_view> operator*() const {
      return {_get_key(js), json_view{js.root, js.index + 1}};
    }
  };
  struct iterator_wrapper {
    json_view begin_;
    json_view end_;
    iterator  begin() { return {begin_}; }
    iterator  end() { return {end_}; }
  };
  if (!is_valid(js)) return iterator_wrapper{{js.root}, {js.root}};
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::object) {
    set_error_view(js, "object expected");
    return iterator_wrapper{{js.root}, {js.root}};
  }
  return iterator_wrapper{
      {js.root, js.index + 1}, {js.root, js.index + 1 + jsv._object.skip}};
}
inline auto iterate_object(json_cview js) {
  struct iterator {
    json_cview js;
    bool       operator!=(const iterator& other) {
      return is_valid(js) && js.index != other.js.index;
    }
    iterator& operator++() {
      if (!is_valid(js)) return *this;
      auto& jst = _get_type(json_cview{js.root, js.index + 1});
      auto& jsv = _get_value(json_cview{js.root, js.index + 1});
      js.index += 2;
      if (jst == json_type::array) js.index += jsv._array.skip;
      if (jst == json_type::object) js.index += jsv._object.skip;
      return *this;
    }
    pair<string_view, json_cview> operator*() const {
      return {_get_key(js), json_cview{js.root, js.index + 1}};
    }
  };
  struct iterator_wrapper {
    json_cview begin_;
    json_cview end_;
    iterator   begin() { return {begin_}; }
    iterator   end() { return {end_}; }
  };
  if (!is_valid(js)) return iterator_wrapper{{js.root}, {js.root}};
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::object) {
    set_error_view(js, "object expected");
    return iterator_wrapper{{js.root}, {js.root}};
  }
  return iterator_wrapper{
      {js.root, js.index + 1}, {js.root, js.index + 1 + jsv._object.skip}};
}

// Object
inline bool set_object(json_view js) {
  if (!is_valid(js)) return false;
  if (js.index != js.root->values.size() - 1)
    throw std::out_of_range{"can only add at the end"};
  auto& jst  = _get_type(js);
  auto& jsv  = _get_value(js);
  jst        = json_type::object;
  jsv._array = {0, 0};
  return true;
}
inline bool object_size(json_cview js, size_t& size) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::object) return set_error(js, "object expected");
  size = (size_t)jsv._object.length;
  return true;
}
inline json_view get_element(json_view js, string_view key) {
  if (!is_valid(js)) return {js.root};
  auto& jst = _get_type(js);
  if (jst != json_type::object) return set_error_view(js, "object expected");
  for (auto [okey, ejs] : iterate_object(js)) {
    if (okey == key) return ejs;
  }
  return set_error_view(js, "missing key " + string{key});
}
inline json_cview get_element(json_cview js, string_view key) {
  if (!is_valid(js)) return {js.root};
  auto& jst = _get_type(js);
  if (jst != json_type::object) return set_error_view(js, "object expected");
  for (auto [okey, value] : iterate_object(js)) {
    if (okey == key) return value;
  }
  return set_error_view(js, "missing key " + string{key});
}
inline bool has_element(json_view js, string_view key) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  if (jst != json_type::object) return set_error(js, "object expected");
  for (auto [okey, value] : iterate_object(js)) {
    if (okey == key) return true;
  }
  return false;
}
inline bool has_element(json_cview js, string_view key) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  if (jst != json_type::object) return set_error(js, "object expected");
  for (auto [okey, value] : iterate_object(js)) {
    if (okey == key) return true;
  }
  return false;
}
inline json_view insert_element(json_view js, string_view key) {
  if (!is_valid(js)) return {js.root};
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::object) return set_error_view(js, "object expected");
  if (js.index + 1 + jsv._object.skip != js.root->values.size())
    throw std::out_of_range{"can only add at the end"};
  for (auto [okey, ejs] : iterate_object(js)) {
    if (okey == key) return ejs;
  }
  jsv._object.length += 1;
  auto& jkt = js.root->types.emplace_back();
  auto& jkv = js.root->values.emplace_back();
  jkt       = json_type::string_;
  for (auto kv : js.root->key_list) {
    auto okey = string_view{
        js.root->keys.data() + kv._string.start, kv._string.length};
    if (okey == key) jkv = kv;
  }
  if (jkv._string.length == 0) {
    jkv._string.start  = (uint32_t)js.root->keys.size();
    jkv._string.length = (uint32_t)key.size();
    js.root->keys.insert(js.root->keys.end(), key.begin(), key.end());
    js.root->keys.push_back(0);
  }
  auto index = (uint32_t)js.root->values.size();
  js.root->types.emplace_back(json_type::null);
  js.root->values.emplace_back();
  auto stack = vector<json_view>{};
  _find_path(js, stack);
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
  return {js.root, index};
}

// Binary
inline bool set_binary(json_view js, const json_binary& value) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  // TODO(fabio): implement reuse
  jst                = json_type::binary;
  jsv._binary.start  = (uint32_t)js.root->binaries.size();
  jsv._binary.length = (uint32_t)value.size();
  js.root->binaries.insert(js.root->binaries.end(), value.begin(), value.end());
  return true;
}
inline bool get_binary(json_cview js, json_binary& value) {
  if (!is_valid(js)) return false;
  auto& jst = _get_type(js);
  auto& jsv = _get_value(js);
  if (jst != json_type::binary) return set_error(js, "binary expected");
  value = json_binary{js.root->binaries.data() + jsv._binary.start,
      js.root->binaries.data() + jsv._binary.start + jsv._binary.length};
  return true;
}

// Get the path of a json view
inline bool _compute_path(json_cview js, json_cview jsv, string& path) {
  if (!is_valid(js) || !is_valid(jsv)) {
    return false;
  } else if (js.index == jsv.index) {
    path = "/";
    return true;
  } else if (is_array(js)) {
    auto idx = 0;
    for (auto ejs : iterate_array(js)) {
      if (!_compute_path(ejs, jsv, path)) continue;
      if (path.back() == '/') path.pop_back();
      path = "/"s + std::to_string(idx) + path;
      return true;
    }
    return false;
  } else if (is_object(js)) {
    for (auto [key, ejs] : iterate_object(js)) {
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
inline string compute_path(json_cview js) {
  auto path = string{};
  if (_compute_path({js.root, 0}, js, path)) {
    return path;
  } else {
    return "";
  }
}

// Conversion from json to values
inline bool get_value(json_cview js, int64_t& value) {
  return get_integral(js, value);
}
inline bool get_value(json_cview js, int32_t& value) {
  return get_integral(js, value);
}
inline bool get_value(json_cview js, uint64_t& value) {
  return get_integral(js, value);
}
inline bool get_value(json_cview js, uint32_t& value) {
  return get_integral(js, value);
}
inline bool get_value(json_cview js, double& value) {
  return get_number(js, value);
}
inline bool get_value(json_cview js, float& value) {
  return get_number(js, value);
}
inline bool get_value(json_cview js, bool& value) {
  return get_boolean(js, value);
}
inline bool get_value(json_cview js, string& value) {
  return get_string(js, value);
}
template <typename T>
inline bool get_value(json_cview js, vector<T>& value) {
  if (!is_valid(js)) return false;
  if (!is_array(js)) return set_error(js, "array expected");
  value.clear();
  auto size = (size_t)0;
  array_size(js, size);
  value.reserve(size);
  for (auto ejs : iterate_array(js)) {
    if (!get_value(ejs, value.emplace_back())) return false;
  }
  return true;
}
template <typename T, size_t N>
inline bool get_value(json_cview js, array<T, N>& value) {
  if (!is_valid(js)) return false;
  if (!is_array(js)) return set_error(js, "array expected");
  auto size = (size_t)0;
  array_size(js, size);
  if (size != N) return set_error(js, "size mismatched");
  auto idx = 0;
  for (auto ejs : iterate_array(js)) {
    if (!get_value(ejs, value.at(idx++))) return false;
  }
  return true;
}

// Get value at a key or index
template <typename T>
inline bool get_value_at(json_cview js, string_view key, T& value) {
  if (!is_valid(js)) return false;
  auto element = get_element(js, key);
  return get_value(element, value);
}
template <typename T>
inline bool get_value_at(json_cview js, size_t idx, T& value) {
  if (!is_valid(js)) return false;
  auto element = get_element(js, idx);
  if (!is_valid(element)) return false;
  return get_value(element, value);
}

// Get value at a key or nothing is key is not preesent
template <typename T>
inline bool get_value_if(json_cview js, string_view key, T& value) {
  if (!is_valid(js)) return false;
  if (!has_element(js, key)) return true;
  auto element = get_element(js, key);
  if (!is_valid(element)) return false;
  return get_value(element, value);
}

// Conversion to json from values
inline bool set_value(json_view js, int64_t value) {
  return set_integral(js, value);
}
inline bool set_value(json_view js, int32_t value) {
  return set_integral(js, value);
}
inline bool set_value(json_view js, uint64_t value) {
  return set_integral(js, value);
}
inline bool set_value(json_view js, uint32_t value) {
  return set_integral(js, value);
}
inline bool set_value(json_view js, double value) {
  return set_number(js, value);
}
inline bool set_value(json_view js, float value) {
  return set_number(js, value);
}
inline bool set_value(json_view js, bool value) {
  return set_boolean(js, value);
}
inline bool set_value(json_view js, const string& value) {
  return set_string(js, value);
}
inline bool set_value(json_view js, const char* value) {
  return set_string(js, value);
}
template <typename T>
inline bool set_value(json_view js, const vector<T>& value) {
  if (!set_array(js)) return false;
  for (auto& v : value) {
    if (!set_value(append_element(js), v)) return false;
  }
  return true;
}
template <typename T, size_t N>
inline bool set_value(json_view js, const array<T, N>& value) {
  if (!set_array(js)) return false;
  for (auto& v : value) {
    if (!set_value(append_element(js), v)) return false;
  }
  return true;
}

// Helpers for user-defined types
inline bool check_array(json_cview js) {
  if (!is_valid(js)) return false;
  if (!is_array(js)) return set_error(js, "array expected");
  return true;
}
inline bool check_array(json_cview js, size_t size_) {
  if (!is_valid(js)) return false;
  if (!is_array(js)) return set_error(js, "array expected");
  auto size = (size_t)0;
  if (!array_size(js, size) || size != size_)
    return set_error(js, "mismatched size");
  return true;
}
inline bool check_object(json_cview js) {
  if (!is_valid(js)) return false;
  if (!is_object(js)) return set_error(js, "array expected");
  return true;
}

// Helpers for user-defined types
template <typename T>
inline bool set_value_at(json_view js, size_t idx, const T& value) {
  if (!is_valid(js)) return false;
  auto ejs = get_element(js, idx);
  if (!is_valid(ejs)) return false;
  return set_value(ejs, value);
}
template <typename T>
inline bool append_value(json_view js, const T& value) {
  if (!is_valid(js)) return false;
  auto ejs = append_element(js);
  if (!is_valid(ejs)) return false;
  return set_value(ejs, value);
}
inline json_view append_array(json_view js) {
  if (!is_valid(js)) return {js.root};
  auto ejs = append_element(js);
  if (!is_valid(ejs)) return {js.root};
  if (!set_array(ejs)) return {js.root};
  return ejs;
}
inline json_view append_object(json_view js) {
  if (!is_valid(js)) return {js.root};
  auto ejs = append_element(js);
  if (!is_valid(ejs)) return {js.root};
  if (!set_object(ejs)) return {js.root};
  return ejs;
}

// Helpers for user-defined types
template <typename T>
inline bool set_value_at(json_view js, string_view key, const T& value) {
  if (!is_valid(js)) return false;
  auto ejs = get_element(js, key);
  if (!is_valid(ejs)) return false;
  return set_value(ejs, value);
}
template <typename T>
inline bool insert_value(json_view js, string_view key, const T& value) {
  if (!is_valid(js)) return false;
  auto ejs = insert_element(js, key);
  if (!is_valid(ejs)) return false;
  return set_value(ejs, value);
}
template <typename T>
inline bool insert_value_if(
    json_view js, string_view key, const T& value, const T& default_) {
  if (!is_valid(js)) return false;
  if (value == default_) return true;
  auto ejs = insert_element(js, key);
  if (!is_valid(ejs)) return false;
  return set_value(ejs, value);
}
inline json_view insert_array(json_view js, string_view key) {
  if (!is_valid(js)) return {js.root};
  auto ejs = insert_element(js, key);
  if (!is_valid(ejs)) return {js.root};
  if (!set_array(ejs)) return {js.root};
  return ejs;
}
inline json_view insert_object(json_view js, string_view key) {
  if (!is_valid(js)) return {js.root};
  auto ejs = insert_element(js, key);
  if (!is_valid(ejs)) return {js.root};
  if (!set_object(ejs)) return {js.root};
  return ejs;
}

// Helpers that need to be declared here
inline bool _find_anchestors(
    json_view js, uint32_t index, vector<json_view>& path) {
  auto& jst = _get_type(js);
  if (jst == json_type::array) {
    if (js.index == index) {
      path.push_back(js);
      return true;
    }
    //    if (index <= js.index || index >= js.index + 1 + jsv._array.skip)
    //      return false;
    //    path.push_back(js.index);
    for (auto ejs : iterate_array(js)) {
      if (_find_anchestors(ejs, index, path)) {
        path.push_back(ejs);
        return true;
      }
    }
    return false;
  } else if (jst == json_type::object) {
    if (js.index == index) {
      path.push_back(js);
      return true;
    }
    //    if (index <= js.index || index >= js.index + 1 + jsv._object.skip)
    //      return false;
    //    path.push_back(js.index);
    for (auto [okey, ejs] : iterate_object(js)) {
      if (_find_anchestors(ejs, index, path)) {
        path.push_back(js);
        return true;
      }
    }
    return false;
  } else {
    return false;
  }
}
inline void _find_path(json_view js, vector<json_view>& path) {
  path.clear();
  _find_anchestors({js.root, 0}, js.index, path);
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
