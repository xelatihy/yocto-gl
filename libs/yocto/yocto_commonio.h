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
// Prints a message to the console and exit with an error. Returns error code.
int print_fatal(const string& msg);

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
struct cli_command;
cli_command make_cli(const string& cmd, const string& usage);
// parse arguments, checks for errors, and exits on error or help
void parse_cli(cli_command& cli, int argc, const char** argv);
// parse arguments and checks for errors
bool parse_cli(cli_command& cli, int argc, const char** argv, string& error);
// gets usage message
string get_usage(const cli_command& cli);
// gets whether help was invoked
bool get_help(const cli_command& cli);
// gets the set command
string get_command(const cli_command& cli);

// Add an optional argument. Supports strings, numbers, and boolean flags.
// Optional arguments will be parsed with name `--<name>` and `-<alt>`.
// Optional booleans will support both `--<name>` and `--no-<name>` to enabled
// and disable the flag.
template <typename T>
inline void add_optional(cli_command& cli, const string& name, T& value,
    const string& usage, const string& alt = "", bool req = false);
// Add a positional argument. Supports strings, numbers, and boolean flags.
template <typename T>
inline void add_positional(cli_command& cli, const string& name, T& value,
    const string& usage, bool req = true);
// Add an optional argument with values as labels. Supports integers, enums and
// strings.
template <typename T>
inline void add_optional(cli_command& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, const string& alt = "",
    bool req = false);
// Add a positional argument with values as labels. Supports string, integers
// and enums.
template <typename T>
inline void add_positional(cli_command& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, bool req = true);
// Add a positional argument that consumes all arguments left.
// Supports strings and enums.
template <typename T>
inline void add_positional(cli_command& cli, const string& name,
    vector<T>& value, const string& usage, bool req = true);
template <typename T>
inline void add_positional(cli_command& cli, const string& name,
    vector<T>& value, const string& usage, const vector<string>& choices,
    bool req = true);

// Add a subcommand
cli_command& add_command(
    cli_command& cli, const string& name, const string& usage);
void add_command_name(
    cli_command& cli, const string& name, string& value, const string& usage);

// Parses an optional or positional argument. Optional arguments' names start
// with "--" or "-", otherwise they are arguments. Supports strings, numbers,
// boolean flags and enums.
// Many names, separated by commas, can be used for each argument.
// Boolean flags are indicated with a pair of names "--name/--no-name", so that
// both options are explicitly specified.
template <typename T>
inline void add_option(cli_command& cli, const string& name, T& value,
    const string& usage, bool req = false);
// Parses an optional or positional argument where values can only be within a
// set of choices. Supports strings, integers and enums.
template <typename T>
inline void add_option(cli_command& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, bool req = false);
// Parse all arguments left on the command line. Can only be used as argument.
template <typename T>
inline void add_option(cli_command& cli, const string& name, vector<T>& value,
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
  int64_t  integer  = 0;
  uint64_t uinteger = 0;
  double   number   = 0;
  string   text     = "";
};
// Command line option. All data should be considered private.
struct cli_option {
  string                            name       = "";
  string                            alt        = "";
  bool                              positional = false;
  cli_type                          type       = cli_type::string;
  bool                              req        = false;
  int                               nargs      = 0;
  string                            usage      = "";
  vector<cli_value>                 value      = {};
  vector<cli_value>                 def        = {};
  vector<string>                    choices    = {};
  bool                              set        = false;
  function<bool(const cli_option&)> set_value  = {};
};
// Command line command. All data should be considered private.
struct cli_command {
  string                        name        = "";
  string                        usage       = "";
  vector<cli_command>           commands    = {};
  vector<cli_option>            options     = {};
  bool                          help        = false;
  string                        command     = "";
  function<void(const string&)> set_command = {};
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
inline vector<cli_value> make_cli_values(const vector<T>& values) {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  auto cvalues = vector<cli_value>{};
  for (auto value : values) {  // needs to copy for booleans
    auto& cvalue = cvalues.emplace_back();
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
  return cvalues;
}
template <typename T>
inline vector<cli_value> make_cli_values(const T& value) {
  return make_cli_values(vector<T>{value});
}

template <typename T>
inline bool get_value(const cli_option& option, T& value) {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  if (option.value.size() != 1) throw std::out_of_range{"bad option size"};
  auto& cvalue = option.value[0];
  if constexpr (std::is_same_v<T, string>) {
    if (option.type != cli_type::string) return false;
    value = cvalue.text;
  } else if constexpr (std::is_same_v<T, bool>) {
    if (option.type != cli_type::boolean) return false;
    value = cvalue.integer != 0;
  } else if constexpr (std::is_enum_v<T>) {
    if (option.type != cli_type::integer) return false;
    value = (T)cvalue.integer;
  } else if constexpr (std::is_integral_v<T> && !std::is_unsigned_v<T>) {
    if (option.type != cli_type::integer) return false;
    value = (T)cvalue.integer;
  } else if constexpr (std::is_integral_v<T> && std::is_unsigned_v<T>) {
    if (option.type != cli_type::uinteger) return false;
    value = (T)cvalue.uinteger;
  } else if constexpr (std::is_floating_point_v<T>) {
    if (option.type != cli_type::number) return false;
    value = (T)cvalue.number;
  } else {
    return false;
  }
  return true;
}

template <typename T>
inline bool get_value(const cli_option& option, vector<T>& values) {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  values.clear();
  for (auto& cvalue : option.value) {
    auto& value = values.emplace_back();
    if constexpr (std::is_same_v<T, string>) {
      if (option.type != cli_type::string) return false;
      value = cvalue.text;
    } else if constexpr (std::is_same_v<T, bool>) {
      if (option.type != cli_type::boolean) return false;
      value = cvalue.integer != 0;
    } else if constexpr (std::is_enum_v<T>) {
      if (option.type != cli_type::integer) return false;
      value = (T)cvalue.integer;
    } else if constexpr (std::is_integral_v<T> && !std::is_unsigned_v<T>) {
      if (option.type != cli_type::integer) return false;
      value = (T)cvalue.integer;
    } else if constexpr (std::is_integral_v<T> && std::is_unsigned_v<T>) {
      if (option.type != cli_type::uinteger) return false;
      value = (T)cvalue.uinteger;
    } else if constexpr (std::is_floating_point_v<T>) {
      if (option.type != cli_type::number) return false;
      value = (T)cvalue.number;
    } else {
      return false;
    }
  }
  return true;
}

template <typename T>
inline void add_optional(cli_command& cli, const string& name, T& value,
    const string& usage, const string& alt, bool req) {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  auto& option      = cli.options.emplace_back();
  option.name       = name;
  option.alt        = alt;
  option.positional = false;
  option.type       = get_cli_type<T>();
  option.req        = req;
  option.nargs      = !std::is_same_v<T, bool> ? 1 : 0;
  option.usage      = usage;
  option.value      = make_cli_values(value);
  option.def        = make_cli_values(value);
  option.choices    = {};
  option.set_value  = [&value](const cli_option& option) -> bool {
    return get_value(option, value);
  };
}

template <typename T>
inline void add_positional(cli_command& cli, const string& name, T& value,
    const string& usage, bool req) {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  auto& option      = cli.options.emplace_back();
  option.name       = name;
  option.alt        = "";
  option.positional = true;
  option.type       = get_cli_type<T>();
  option.req        = req;
  option.nargs      = 1;
  option.usage      = usage;
  option.value      = make_cli_values(value);
  option.def        = make_cli_values(value);
  option.choices    = {};
  option.set_value  = [&value](const cli_option& option) -> bool {
    return get_value(option, value);
  };
}

template <typename T>
inline void add_optional(cli_command& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, const string& alt,
    bool req) {
  static_assert(
      std::is_same_v<T, string> || std::is_integral_v<T> || std::is_enum_v<T>,
      "unsupported type");
  auto& option      = cli.options.emplace_back();
  option.name       = name;
  option.alt        = alt;
  option.positional = false;
  option.type       = get_cli_type<T>();
  option.req        = req;
  option.nargs      = 1;
  option.usage      = usage;
  option.value      = make_cli_values(value);
  option.def        = make_cli_values(value);
  option.choices    = choices;
  option.set_value  = [&value](const cli_option& option) -> bool {
    return get_value(option, value);
  };
}

template <typename T>
inline void add_positional(cli_command& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, bool req) {
  static_assert(
      std::is_same_v<T, string> || std::is_integral_v<T> || std::is_enum_v<T>,
      "unsupported type");
  auto& option      = cli.options.emplace_back();
  option.name       = name;
  option.alt        = "";
  option.positional = true;
  option.type       = get_cli_type<T>();
  option.req        = req;
  option.nargs      = 1;
  option.usage      = usage;
  option.value      = make_cli_values(value);
  option.def        = make_cli_values(value);
  option.choices    = choices;
  option.set_value  = [&value](const cli_option& option) -> bool {
    return get_value(option, value);
  };
}

template <typename T>
inline void add_positional(cli_command& cli, const string& name,
    vector<T>& value, const string& usage, bool req) {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  auto& option      = cli.options.emplace_back();
  option.name       = name;
  option.alt        = "";
  option.positional = true;
  option.type       = get_cli_type<T>();
  option.req        = req;
  option.nargs      = -1;
  option.usage      = usage;
  option.value      = make_cli_values(value);
  option.def        = make_cli_values(value);
  option.choices    = {};
  option.set_value  = [&value](const cli_option& option) -> bool {
    return get_value(option, value);
  };
}

template <typename T>
inline void add_positional(cli_command& cli, const string& name,
    vector<T>& value, const string& usage, const vector<string>& choices,
    bool req) {
  static_assert(
      std::is_same_v<T, string> || std::is_integral_v<T> || std::is_enum_v<T>,
      "unsupported type");
  auto& option      = cli.options.emplace_back();
  option.name       = name;
  option.alt        = "";
  option.positional = true;
  option.type       = get_cli_type<T>();
  option.req        = req;
  option.nargs      = -1;
  option.usage      = usage;
  option.value      = make_cli_values(value);
  option.def        = make_cli_values(value);
  option.choices    = choices;
  option.set_value  = [&value](const cli_option& option) -> bool {
    return get_value(option, value);
  };
}

inline vector<string> split_cli_names(const string& name_) {
  auto name  = name_;
  auto split = vector<string>{};
  if (name.empty()) throw std::invalid_argument("option name cannot be empty");
  if (name.find_first_of(" \t\r\n") != string::npos)
    throw std::invalid_argument("option name cannot contain whitespaces");
  while (name.find_first_of(",/") != string::npos) {
    auto pos = name.find_first_of(",/");
    if (pos > 0) split.push_back(name.substr(0, pos));
    name = name.substr(pos + 1);
  }
  if (!name.empty()) split.push_back(name);
  if (split.empty()) throw std::invalid_argument("option name cannot be empty");
  for (auto& name : split)
    if ((split[0][0] == '-') != (name[0] == '-'))
      throw std::invalid_argument("inconsistent option names for " + name);
  return split;
}

template <typename T>
inline void add_option(cli_command& cli, const string& name, T& value,
    const string& usage, bool req) {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  auto& option = cli.options.emplace_back();
  auto  names  = split_cli_names(name);
  if (names[0].find('-') == 0) {
    option.name = names[0].substr(1);
    if (option.name.find('-') == 0) option.name = option.name.substr(1);
    if (names.size() > 1 && names.at(1).find("--no-") != 0) {
      option.alt = names[1].substr(1);
      if (option.alt.find('-') == 0) option.alt = option.alt.substr(1);
    }
  } else {
    option.name       = name;
    option.alt        = "";
    option.positional = true;
  }
  option.type      = get_cli_type<T>();
  option.req       = req;
  option.nargs     = !std::is_same_v<T, bool> ? 1 : 0;
  option.usage     = usage;
  option.value     = make_cli_values(value);
  option.def       = make_cli_values(value);
  option.choices   = {};
  option.set_value = [&value](const cli_option& option) -> bool {
    return get_value(option, value);
  };
}

template <typename T>
inline void add_option(cli_command& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, bool req) {
  static_assert(
      std::is_same_v<T, string> || std::is_integral_v<T> || std::is_enum_v<T>,
      "unsupported type");
  auto& option = cli.options.emplace_back();
  auto  names  = split_cli_names(name);
  if (names[0].find('-') == 0) {
    option.name = names[0].substr(1);
    if (option.name.find('-') == 0) option.name = option.name.substr(1);
    if (names.size() > 1 && names.at(1).find("--no-") != 0) {
      option.alt = names[1].substr(1);
      if (option.alt.find('-') == 0) option.alt = option.alt.substr(1);
    }
  } else {
    option.name       = name;
    option.alt        = "";
    option.positional = true;
  }
  option.type      = get_cli_type<T>();
  option.req       = req;
  option.nargs     = 1;
  option.usage     = usage;
  option.value     = make_cli_values(value);
  option.def       = make_cli_values(value);
  option.choices   = choices;
  option.set_value = [&value](const cli_option& option) -> bool {
    return get_value(option, value);
  };
}

template <typename T>
inline void add_option(cli_command& cli, const string& name, vector<T>& value,
    const string& usage, bool req) {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  auto& option = cli.options.emplace_back();
  auto  names  = split_cli_names(name);
  if (names[0].find('-') == 0) {
    option.name = names[0].substr(1);
    if (option.name.find('-') == 0) option.name = option.name.substr(1);
    if (names.size() > 1 && names.at(1).find("--no-") != 0) {
      option.alt = names[1].substr(1);
      if (option.alt.find('-') == 0) option.alt = option.alt.substr(1);
    }
  } else {
    option.name       = name;
    option.alt        = "";
    option.positional = true;
  }
  option.type      = get_cli_type<T>();
  option.req       = req;
  option.nargs     = -1;
  option.usage     = usage;
  option.value     = make_cli_values(value);
  option.def       = make_cli_values(value);
  option.choices   = {};
  option.set_value = [&value](const cli_option& option) -> bool {
    return get_value(option, value);
  };
}

// Backward compatibility
using cli_state [[deprecated]] = cli_command;

}  // namespace yocto

#endif
