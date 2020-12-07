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
  json_value() {}
  json_value(const json_value& other) { _copy(other); }
  json_value(json_value&& other) { _copy(std::move(other)); }
  explicit json_value(std::nullptr_t) : _type{json_type::null}, _unsigned{0} {}
  explicit json_value(int64_t value)
      : _type{json_type::integer}, _integer{value} {}
  explicit json_value(int32_t value)
      : _type{json_type::integer}, _integer{value} {}
  explicit json_value(uint64_t value)
      : _type{json_type::unsigned_}, _unsigned{value} {}
  explicit json_value(uint32_t value)
      : _type{json_type::unsigned_}, _unsigned{value} {}
  explicit json_value(double value) : _type{json_type::real}, _real{value} {}
  explicit json_value(float value) : _type{json_type::real}, _real{value} {}
  explicit json_value(bool value)
      : _type{json_type::boolean}, _boolean{value} {}
  explicit json_value(const string& value)
      : _type{json_type::string_}, _string_{new string{value}} {}
  explicit json_value(const char* value)
      : _type{json_type::string_}, _string_{new string{value}} {}
  explicit json_value(const json_array& value)
      : _type{json_type::array}, _array{new json_array{value}} {}
  explicit json_value(const json_object& value)
      : _type{json_type::object}, _object{new json_object{value}} {}
  explicit json_value(const json_binary& value)
      : _type{json_type::binary}, _binary{new json_binary{value}} {}

  // assignments
  json_value& operator=(const json_value& other) { return _copy(other); }
  json_value& operator=(json_value&& other) { return _copy(std::move(other)); }
  json_value& operator=(std::nullptr_t) { return _set(nullptr); }
  json_value& operator=(int64_t value) { return _set(value); }
  json_value& operator=(int32_t value) { return _set((int64_t)value); }
  json_value& operator=(uint64_t value) { return _set(value); }
  json_value& operator=(uint32_t value) { return _set((uint64_t)value); }
  json_value& operator=(double value) { return _set(value); }
  json_value& operator=(float value) { return _set((double)value); }
  json_value& operator=(bool value) { return _set(value); }
  json_value& operator=(const string& value) { return _set(value); }
  json_value& operator=(const char* value) { return _set(value); }
  json_value& operator=(const json_array& value) { return _set(value); }
  json_value& operator=(const json_object& value) { return _set(value); }
  json_value& operator=(const json_binary& value) { return _set(value); }

  // type
  json_type type() const { return _type; }
  bool      is_null() const { return _type == json_type::null; }
  bool      is_integer() const { return _type == json_type::integer; }
  bool      is_unsigned() const { return _type == json_type::unsigned_; }
  bool      is_real() const { return _type == json_type::real; }
  bool      is_bool() const { return _type == json_type::boolean; }
  bool      is_string() const { return _type == json_type::string_; }
  bool      is_array() const { return _type == json_type::array; }
  bool      is_object() const { return _type == json_type::object; }
  bool      is_binary() const { return _type == json_type::binary; }

  // conversions
  explicit operator int64_t() const {
    return is_unsigned() ? (int64_t)get_unsigned() : get_integer();
  }
  explicit operator int32_t() const {
    return is_unsigned() ? (int32_t)get_unsigned() : (int32_t)get_integer();
  }
  explicit operator uint64_t() const {
    return is_integer() ? (uint64_t)get_integer() : get_unsigned();
  }
  explicit operator uint32_t() const {
    return is_integer() ? (uint32_t)get_integer() : (uint32_t)get_unsigned();
  }
  explicit operator double() const {
    return is_integer() ? (double)get_integer()
                        : is_unsigned() ? (double)get_unsigned() : get_real();
  }
  explicit operator float() const {
    return is_integer()
               ? (float)get_integer()
               : is_unsigned() ? (float)get_unsigned() : (float)get_real();
  }
  explicit operator bool() const { return get_boolean(); }
  explicit operator string() const { return get_string(); }
  explicit operator json_array() const { return get_array(); }
  explicit operator json_object() const { return get_object(); }
  explicit operator json_binary() const { return get_binary(); }

  // access
  const int64_t& get_integer() const {
    _check_type(json_type::integer);
    return _integer;
  }
  const uint64_t& get_unsigned() const {
    _check_type(json_type::unsigned_);
    return _unsigned;
  }
  const double& get_real() const {
    _check_type(json_type::real);
    return _real;
  }
  const bool& get_boolean() const {
    _check_type(json_type::boolean);
    return _boolean;
  }
  const string& get_string() const {
    _check_type(json_type::string_);
    return *_string_;
  }
  const json_array& get_array() const {
    _check_type(json_type::array);
    return *_array;
  }
  const json_object& get_object() const {
    _check_type(json_type::object);
    return *_object;
  }
  const json_binary& get_binary() const {
    _check_type(json_type::binary);
    return *_binary;
  }
  int64_t& get_integer() {
    _check_type(json_type::integer);
    return _integer;
  }
  uint64_t& get_unsigned() {
    _check_type(json_type::unsigned_);
    return _unsigned;
  }
  double& get_real() {
    _check_type(json_type::real);
    return _real;
  }
  bool& get_boolean() {
    _check_type(json_type::boolean);
    return _boolean;
  }
  string& get_string() {
    _check_type(json_type::string_);
    return *_string_;
  }
  json_array& get_array() {
    _check_type(json_type::array);
    return *_array;
  }
  json_object& get_object() {
    _check_type(json_type::object);
    return *_object;
  }
  json_binary& get_binary() {
    _check_type(json_type::binary);
    return *_binary;
  }

  // structuree support
  bool   empty() const { return size() == 0; }
  size_t size() const { return _size(); }

  // array support
  json_value&       operator[](size_t idx) { return get_array().at(idx); }
  const json_value& operator[](size_t idx) const { return get_array().at(idx); }
  json_value&       at(size_t idx) { return get_array().at(idx); }
  const json_value& at(size_t idx) const { return get_array().at(idx); }
  void              push_back(const json_value& value) {
    return get_array().push_back(value);
  }
  void push_back(json_value&& value) {
    return get_array().push_back(std::move(value));
  }
  template <typename... Args>
  json_value& emplace_back(Args&&... args) {
    return get_array().emplace_back(std::forward(args)...);
  }
  json_value*       begin() { return get_array().data(); }
  const json_value* begin() const { return get_array().data(); }
  json_value*       end() { return get_array().data() + get_array().size(); }
  const json_value* end() const {
    return get_array().data() + get_array().size();
  }

  // object support
  json_value& operator[](const string& key) {
    if (auto ptr = find(key); ptr) {
      return *ptr;
    } else {
      return get_object().emplace_back(key, json_value{}).second;
    }
  }
  json_value& at(const string& key) {
    if (auto ptr = find(key); ptr) {
      return *ptr;
    } else {
      throw std::out_of_range{"missing key " + key};
    }
  }
  const json_value& at(const string& key) const {
    if (auto ptr = find(key); ptr) {
      return *ptr;
    } else {
      throw std::out_of_range{"missing key " + key};
    }
  }
  struct json_object_it {
    pair<string, json_value>* _begin;
    pair<string, json_value>* _end;
    pair<string, json_value>* begin() { return _begin; }
    pair<string, json_value>* end() { return _end; }
  };
  struct json_object_cit {
    const pair<string, json_value>* _begin;
    const pair<string, json_value>* _end;
    const pair<string, json_value>* begin() { return _begin; }
    const pair<string, json_value>* end() { return _end; }
  };
  json_object_it items() {
    return {get_object().data(), get_object().data() + get_object().size()};
  }
  json_object_cit items() const {
    return {get_object().data(), get_object().data() + get_object().size()};
  }
  json_value* find(const string& key) {
    for (auto& [key_, value] : get_object()) {
      if (key_ == key) return &value;
    }
    return nullptr;
  }
  const json_value* find(const string& key) const {
    for (auto& [key_, value] : get_object()) {
      if (key_ == key) return &value;
    }
    return nullptr;
  }

  // swap
  void swap(json_value& other) { _swap(other); }

  // destructor
  ~json_value() { _set_type(json_type::null); }

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

  // set type
  void _set_type(json_type type) {
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

  // destroy value
  void _destroy() {
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

  // swap operation
  json_value& _swap(json_value& other) {
    std::swap(_type, other._type);
    std::swap(_unsigned, other._unsigned);  // hask to swap bits
    return *this;
  }

  // assignment
  json_value& _copy(const json_value& other) {
    _set_type(other._type);
    switch (_type) {
      case json_type::null: _integer = other._integer; break;
      case json_type::integer: _integer = other._integer; break;
      case json_type::unsigned_: _unsigned = other._unsigned; break;
      case json_type::real: _real = other._real; break;
      case json_type::boolean: _boolean = other._boolean; break;
      case json_type::string_: _string_ = new string{*other._string_}; break;
      case json_type::array: _array = new json_array{*other._array}; break;
      case json_type::object: _object = new json_object{*other._object}; break;
      case json_type::binary: _binary = new json_binary{*other._binary}; break;
    }
    return *this;
  }
  json_value& _copy(json_value&& other) {
    _swap(other);
    return *this;
  }
  template <typename T>
  json_value& _set(const T& value) {
    auto js = json_value{value};
    _swap(js);
    return *this;
  }

  void _check_type(json_type type) const {
    if (_type != type) throw json_type_error{"bad json type"};
  }

  size_t _size() const {
    switch (_type) {
      case json_type::string_: return get_string().size();
      case json_type::array: return get_array().size();
      case json_type::object: return get_object().size();
      case json_type::binary: return get_binary().size();
      default: throw json_type_error{"bad json type"};
    }
  }
};

// Load/save a json file
bool load_json(const string& filename, json_value& js, string& error);
bool save_json(const string& filename, const json_value& js, string& error);

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
