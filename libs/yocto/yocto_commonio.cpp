//
// Implementation for Yocto/CommonIO
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

#include "yocto_commonio.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <vector>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::unordered_set;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// PRINT/FORMATTING UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Print a message to the console
void print_info(const string& msg) { printf("%s\n", msg.c_str()); }
// Prints a messgae to the console and exit with an error.
int print_fatal(const string& msg) {
  printf("\n%s\n", msg.c_str());
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
print_timer print_timed(const string& msg) {
  printf("%s", msg.c_str());
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
// FILE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Cleanup
file_stream::~file_stream() {
  if (owned && fs) fclose(fs);
}

// Open a file
file_stream open_file(const string& filename, const string& mode) {
#ifdef _WIN32
  auto path8 = std::filesystem::u8path(filename);
  auto wmode = std::wstring(mode.begin(), mode.end());
  auto fs    = _wfopen(path8.c_str(), wmode.c_str());
#else
  auto fs = fopen(filename.c_str(), mode.c_str());
#endif
  return {filename, fs, true};
}

// Close a file
void close_file(file_stream& fs) {
  if (fs.owned && fs.fs) fclose(fs.fs);
  fs.filename = "";
  fs.fs       = nullptr;
  fs.owned    = false;
}

// Read a line of text
bool read_line(file_stream& fs, char* buffer, size_t size) {
  return fgets(buffer, (int)size, fs.fs);
}

// Write text to a file
bool write_text(file_stream& fs, const string& str) {
  return fprintf(fs.fs, "%s", str.c_str()) >= 0;
}

// Read data from a file
bool read_data(file_stream& fs, void* buffer, size_t count) {
  return fread(buffer, 1, count, fs.fs) == count;
}

// Write data from a file
bool write_data(file_stream& fs, const void* buffer, size_t count) {
  return fwrite(buffer, 1, count, fs.fs) == count;
}

// Opens a file with a utf8 file name
FILE* fopen_utf8(const char* filename, const char* mode) {
#ifdef _WIN32
  auto path8 = std::filesystem::u8path(filename);
  auto wmode = std::wstring(string{mode}.begin(), string{mode}.end());
  return _wfopen(path8.c_str(), wmode.c_str());
#else
  return fopen(filename, mode);
#endif
}

// Load a text file
bool load_text(const string& filename, string& str, string& error) {
  // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
  auto fs = open_file(filename, "rb");
  if (!fs) {
    error = filename + ": file not found";
    return false;
  }
  fseek(fs.fs, 0, SEEK_END);
  auto length = ftell(fs.fs);
  fseek(fs.fs, 0, SEEK_SET);
  str.resize(length);
  if (!read_values(fs, str.data(), length)) {
    error = filename + ": read error";
    return false;
  }
  return true;
}

// Save a text file
bool save_text(const string& filename, const string& str, string& error) {
  auto fs = open_file(filename, "wt");
  if (!fs) {
    error = filename + ": file not found";
    return false;
  }
  if (!write_text(fs, str)) {
    error = filename + ": write error";
    return false;
  }
  return true;
}

// Load a binary file
bool load_binary(const string& filename, vector<byte>& data, string& error) {
  // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
  auto fs = open_file(filename, "rb");
  if (!fs) {
    error = filename + ": file not found";
    return false;
  }
  fseek(fs.fs, 0, SEEK_END);
  auto length = ftell(fs.fs);
  fseek(fs.fs, 0, SEEK_SET);
  data.resize(length);
  if (!read_values(fs, data.data(), length)) {
    error = filename + ": read error";
    return false;
  }
  return true;
}

// Save a binary file
bool save_binary(
    const string& filename, const vector<byte>& data, string& error) {
  auto fs = open_file(filename, "wb");
  if (!fs) {
    error = filename + ": file not found";
    return false;
  }
  if (!write_values(fs, data.data(), data.size())) {
    error = filename + ": write error";
    return false;
  }
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PATH UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Make a path from a utf8 string
static std::filesystem::path make_path(const string& filename) {
  return std::filesystem::u8path(filename);
}

// Normalize path
string normalize_path(const string& filename) {
  return make_path(filename).generic_u8string();
}

// Get directory name (not including /)
string path_dirname(const string& filename) {
  return make_path(filename).parent_path().generic_u8string();
}

// Get extension (including .)
string path_extension(const string& filename) {
  return make_path(filename).extension().u8string();
}

// Get filename without directory.
string path_filename(const string& filename) {
  return make_path(filename).filename().u8string();
}

// Get filename without directory and extension.
string path_basename(const string& filename) {
  return make_path(filename).stem().u8string();
}

// Joins paths
string path_join(const string& patha, const string& pathb) {
  return (make_path(patha) / make_path(pathb)).generic_u8string();
}
string path_join(
    const string& patha, const string& pathb, const string& pathc) {
  return (make_path(patha) / make_path(pathb) / make_path(pathc))
      .generic_u8string();
}

// Replaces extensions
string replace_extension(const string& filename, const string& ext) {
  return make_path(filename).replace_extension(ext).u8string();
}

// Check if a file can be opened for reading.
bool path_exists(const string& filename) { return exists(make_path(filename)); }

// Check if a file is a directory
bool path_isdir(const string& filename) {
  return is_directory(make_path(filename));
}

// Check if a file is a file
bool path_isfile(const string& filename) {
  return is_regular_file(make_path(filename));
}

// List the contents of a directory
vector<string> list_directory(const string& filename) {
  auto entries = vector<string>{};
  for (auto entry : std::filesystem::directory_iterator(make_path(filename))) {
    entries.push_back(entry.path().generic_u8string());
  }
  return entries;
}

// Create a directory and all missing parent directories if needed
bool make_directory(const string& dirname, string& error) {
  if (path_exists(dirname)) return true;
  try {
    create_directories(make_path(dirname));
    return true;
  } catch (...) {
    error = dirname + ": cannot create directory";
    return false;
  }
}

// Get the current directory
string path_current() { return std::filesystem::current_path().u8string(); }

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF COMMAND-LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto {

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
inline void add_optional_impl(cli_command& cli, const string& name, T& value,
    const string& usage, const vector<T>& minmax, const vector<string>& choices,
    const string& alt, bool req) {
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
  option.nargs      = std::is_same_v<T, bool> ? 0 : 1;
  option.usage      = usage;
  option.minmax     = make_cli_values(minmax);
  option.choices    = choices;
  option.value      = make_cli_values(value);
  option.def        = make_cli_values(value);
  option.set_value  = [&value](const cli_option& option) -> bool {
    return get_value(option, value);
  };
}

template <typename T>
inline void add_positional_impl(cli_command& cli, const string& name, T& value,
    const string& usage, const vector<T>& minmax, const vector<string>& choices,
    bool req) {
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
  option.minmax     = make_cli_values(value);
  option.choices    = choices;
  option.value      = make_cli_values(value);
  option.def        = make_cli_values(value);
  option.set_value  = [&value](const cli_option& option) -> bool {
    return get_value(option, value);
  };
}

template <typename T>
inline void add_positionalv_impl(cli_command& cli, const string& name,
    vector<T>& value, const string& usage, const vector<T>& minmax,
    const vector<string>& choices, bool req) {
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
  option.minmax     = make_cli_values(value);
  option.choices    = choices;
  option.value      = make_cli_values(value);
  option.def        = make_cli_values(value);
  option.set_value  = [&value](const cli_option& option) -> bool {
    return get_value(option, value);
  };
}

// Add an optional argument. Supports strings, numbers, and boolean flags.
// Optional arguments will be parsed with name `--<name>` and `-<alt>`.
// Optional booleans will support both `--<name>` and `--no-<name>` to enabled
// and disable the flag.
void add_optional(cli_command& cli, const string& name, int& value,
    const string& usage, const vector<int>& minmax, const string& alt,
    bool req) {
  return add_optional_impl(
      cli, name, value, usage, {minmax[0], minmax[1]}, {}, alt, req);
}
void add_optional(cli_command& cli, const string& name, float& value,
    const string& usage, const vector<float>& minmax, const string& alt,
    bool req) {
  return add_optional_impl(cli, name, value, usage, minmax, {}, alt, req);
}
void add_optional(cli_command& cli, const string& name, bool& value,
    const string& usage, const vector<string>& choices, const string& alt,
    bool req) {
  return add_optional_impl(cli, name, value, usage, {}, choices, alt, req);
}
void add_optional(cli_command& cli, const string& name, string& value,
    const string& usage, const vector<string>& choices, const string& alt,
    bool req) {
  return add_optional_impl(cli, name, value, usage, {}, choices, alt, req);
}
void add_optional(cli_command& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, const string& alt,
    bool req) {
  return add_optional_impl(cli, name, value, usage, {}, choices, alt, req);
}
// Add a positional argument. Supports strings, numbers, and boolean flags.
void add_positional(cli_command& cli, const string& name, int& value,
    const string& usage, const vector<int>& minmax, bool req) {
  return add_positional_impl(cli, name, value, usage, minmax, {}, req);
}
void add_positional(cli_command& cli, const string& name, float& value,
    const string& usage, const vector<float>& minmax, bool req) {
  return add_positional_impl(cli, name, value, usage, minmax, {}, req);
}
void add_positional(cli_command& cli, const string& name, bool& value,
    const string& usage, const vector<string>& choices, bool req) {
  return add_positional_impl(cli, name, value, usage, {}, choices, req);
}
void add_positional(cli_command& cli, const string& name, string& value,
    const string& usage, const vector<string>& choices, bool req) {
  return add_positional_impl(cli, name, value, usage, {}, choices, req);
}
void add_positional(cli_command& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, bool req) {
  return add_positional_impl(cli, name, value, usage, {}, choices, req);
}
// Add a positional argument that consumes all arguments left.
// Supports strings and enums.
void add_positional(cli_command& cli, const string& name, vector<int>& value,
    const string& usage, const vector<int>& minmax, bool req) {
  return add_positionalv_impl(cli, name, value, usage, minmax, {}, req);
}
void add_positional(cli_command& cli, const string& name, vector<float>& value,
    const string& usage, const vector<float>& minmax, bool req) {
  return add_positionalv_impl(cli, name, value, usage, minmax, {}, req);
}
void add_positional(cli_command& cli, const string& name, vector<int>& value,
    const string& usage, const vector<string>& choices, bool req) {
  return add_positionalv_impl(cli, name, value, usage, {}, choices, req);
}
void add_positional(cli_command& cli, const string& name, vector<string>& value,
    const string& usage, const vector<string>& choices, bool req) {
  return add_positionalv_impl(cli, name, value, usage, {}, choices, req);
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
inline void add_option_impl(cli_command& cli, const string& name, T& value,
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
inline void add_option_impl(cli_command& cli, const string& name, T& value,
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
inline void add_option_impl(cli_command& cli, const string& name,
    vector<T>& value, const string& usage, bool req) {
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

// Parses an optional or positional argument. Optional arguments' names start
// with "--" or "-", otherwise they are arguments. Supports strings, numbers,
// boolean flags and enums.
// Many names, separated by commas, can be used for each argument.
// Boolean flags are indicated with a pair of names "--name/--no-name", so that
// both options are explicitly specified.
void add_option(cli_command& cli, const string& name, int& value,
    const string& usage, bool req) {
  return add_option_impl(cli, name, value, usage, req);
}
void add_option(cli_command& cli, const string& name, float& value,
    const string& usage, bool req) {
  return add_option_impl(cli, name, value, usage, req);
}
void add_option(cli_command& cli, const string& name, bool& value,
    const string& usage, bool req) {
  return add_option_impl(cli, name, value, usage, req);
}
void add_option(cli_command& cli, const string& name, string& value,
    const string& usage, bool req) {
  return add_option_impl(cli, name, value, usage, req);
}
// Parses an optional or positional argument where values can only be within a
// set of choices. Supports strings, integers and enums.
void add_option(cli_command& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, bool req) {
  return add_option_impl(cli, name, value, usage, choices, req);
}
void add_option(cli_command& cli, const string& name, string& value,
    const string& usage, const vector<string>& choices, bool req) {
  return add_option_impl(cli, name, value, usage, choices, req);
}
// Parse all arguments left on the command line. Can only be used as argument.
void add_option(cli_command& cli, const string& name, vector<string>& value,
    const string& usage, bool req) {
  return add_option_impl(cli, name, value, usage, req);
}

// initialize a command line parser
cli_command make_cli(const string& name, const string& usage) {
  auto cli  = cli_command{};
  cli.name  = name;
  cli.usage = usage;
  cli.commands.reserve(256);
  add_option(cli, "--help", cli.help, "Print usage.");
  return cli;
}

// add command
cli_command& add_command(
    cli_command& cli, const string& name, const string& usage) {
  for (auto& cmd : cli.commands) {
    if (cmd.name == name) {
      throw std::invalid_argument{"cannot add two commands with the same name"};
    }
  }
  auto& cmd = cli.commands.emplace_back();
  cmd.name  = name;
  cmd.usage = usage;
  cmd.commands.reserve(256);
  add_option(cmd, "--help", cmd.help, "Print usage.");
  return cmd;
}

void add_command_name(
    cli_command& cli, const string& name, string& value, const string& usage) {
  cli.set_command = [&value](const string& cvalue) { value = cvalue; };
}

static void validate_names(const cli_command& cmd) {
  // check for errors
  auto used = unordered_set<string>{};
  for (auto& option : cmd.options) {
    if (option.name.empty())
      throw std::invalid_argument("name cannot be empty");
    auto names = split_cli_names(option.name);
    if (names.empty()) throw std::invalid_argument("name cannot be empty");
    for (auto& name : names) {
      if (used.find(name) != used.end())
        throw std::invalid_argument("option name " + name + " already in use");
      used.insert(name);
      if ((name[0] == '-') != (option.name[0] == '-'))
        throw std::invalid_argument("inconsistent option type for " + name);
    }
  }
  for (auto& scmd : cmd.commands) validate_names(scmd);
}

bool get_help(const cli_command& cli) {
  if (cli.help) return true;
  for (auto& cmd : cli.commands) return get_help(cmd);
  return false;
}

static string get_usage(const cli_command& root, const cli_command& cli) {
  auto type_name = [](const cli_option& option) -> string {
    auto str = string{};
    str += "<";
    if (option.nargs < 0) str += "[";
    if (!option.choices.empty()) str += "string";
    switch (option.type) {
      case cli_type::integer: str += "integer"; break;
      case cli_type::uinteger: str += "uinteger"; break;
      case cli_type::number: str += "number"; break;
      case cli_type::string: str += "string"; break;
      case cli_type::boolean: str += "boolean"; break;
    }
    if (option.nargs < 0) str += "]";
    str += ">";
    return str;
  };
  auto def_string = [](const cli_option& option) -> string {
    if (option.req) return string{"[required]"};
    auto str = string{};
    str += "[";
    for (auto& value : option.def) {
      switch (option.type) {
        case cli_type::integer:
          str += option.choices.empty() ? std::to_string(value.integer)
                                        : option.choices[value.integer];
          break;
        case cli_type::uinteger:
          str += option.choices.empty() ? std::to_string(value.uinteger)
                                        : option.choices[value.uinteger];
          break;
        case cli_type::number: str += std::to_string(value.number); break;
        case cli_type::string: str += '\"' + value.text + '\"'; break;
        case cli_type::boolean: str += value.integer ? "true" : "false"; break;
      }
    }
    str += "]";
    return str;
  };

  if (!cli.command.empty()) {
    for (auto& subcommand : cli.commands)
      if (cli.command == subcommand.name) return get_usage(root, subcommand);
  }

  auto message      = string{};
  auto has_optional = false, has_positional = false, has_commands = false;
  auto usage_optional = string{}, usage_positional = string{},
       usage_command = string{};
  for (auto& option : cli.options) {
    auto line = string{};
    if (option.positional) {
      line += "  " + option.name;
    } else {
      line += "  --" + option.name;
      if (!option.alt.empty()) line += ", -" + option.alt;
    }
    if (option.nargs > 0) line += " " + type_name(option);
    while (line.size() < 32) line += " ";
    line += option.usage;
    line += " " + def_string(option) + "\n";
    if (!option.choices.empty()) {
      line += "    with choices: ";
      auto len = 16;
      for (auto& choice : option.choices) {
        if (len + choice.size() + 2 > 78) {
          line += "\n                  ";
          len = 16;
        }
        line += choice + ", ";
        len += choice.size() + 2;
      }
      line = line.substr(0, line.size() - 2);
      line += "\n";
    }
    if (option.positional) {
      has_positional = true;
      usage_positional += line;
    } else {
      has_optional = true;
      usage_optional += line;
    }
  }
  for (auto& scmd : cli.commands) {
    has_commands = true;
    auto line    = "  " + scmd.name;
    while (line.size() < 32) line += " ";
    line += scmd.usage + "\n";
    usage_command += line;
  }
  auto is_command = &cli != &root;
  message += "usage: " + root.name + (is_command ? " " + cli.name : "") +
             (has_commands ? " command" : "") +
             (has_optional ? " [options]" : "") +
             (has_positional ? " <arguments>" : "") + "\n";
  message += cli.usage + "\n\n";
  if (has_commands) {
    message += "commands:\n" + usage_command + "\n";
  }
  if (has_optional) {
    message += "options:\n" + usage_optional + "\n";
  }
  if (has_positional) {
    message += "arguments:\n" + usage_positional + "\n";
  }
  return message;
}

string get_usage(const cli_command& cli) { return get_usage(cli, cli); }

string get_command(const cli_command& cli) { return cli.command; }

static bool parse_value(
    cli_option& option, const vector<string>& args, size_t start) {
  option.value.resize(option.nargs > 0 ? option.nargs : (args.size() - start));
  for (auto idx = (size_t)0; idx < option.value.size(); idx++) {
    auto& value   = option.value.at(idx);
    auto& arg     = args.at(start + idx);
    auto& choices = option.choices;
    if (!choices.empty()) {
      if (std::find(choices.begin(), choices.end(), arg) == choices.end())
        return false;
    }
    switch (option.type) {
      case cli_type::string: {
        value.text = arg;
      } break;
      case cli_type::boolean: {
        if (arg == "true" || arg == "1") {
          value.integer = 1;
        } else if (arg == "false" || arg == "0") {
          value.integer = 0;
        } else {
          return false;
        }
      } break;
      case cli_type::integer: {
        if (choices.empty()) {
          auto end      = (char*)nullptr;
          value.integer = (int)strtol(arg.c_str(), &end, 10);
          if (end == nullptr) return false;
        } else {
          value.integer = (int64_t)(
              std::find(choices.begin(), choices.end(), arg) - choices.begin());
        }
      } break;
      case cli_type::uinteger: {
        if (choices.empty()) {
          auto end       = (char*)nullptr;
          value.uinteger = (int)strtoul(arg.c_str(), &end, 10);
          if (end == nullptr) return false;
        } else {
          value.uinteger = (uint64_t)(
              std::find(choices.begin(), choices.end(), arg) - choices.begin());
        }
      } break;
      case cli_type::number: {
        auto end     = (char*)nullptr;
        value.number = strtod(arg.c_str(), &end);
        if (end == nullptr) return false;
      } break;
    }
  }
  return true;
}

bool parse_cli(cli_command& cli, vector<string>& args, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  // current parsing state
  auto commands    = vector<cli_command*>{&cli};
  auto positionals = vector<int>{0};

  // parse arguments
  for (auto idx = (size_t)0; idx < args.size(); idx++) {
    auto& cmd = *commands.back();
    auto& arg = args[idx];
    if (arg == "--help") {
      cmd.help = true;
      break;
    }
    auto is_positional = args[idx].find('-') != 0;
    if (!cmd.commands.empty() && is_positional) {
      auto pos = std::find_if(cmd.commands.begin(), cmd.commands.end(),
          [&arg](auto& command) { return command.name == arg; });
      if (pos == cmd.commands.end()) return cli_error("unknown command " + arg);
      cmd.command = arg;
      commands.push_back(&(*pos));
      positionals.push_back(0);
      continue;
    } else if (is_positional) {
      auto pos   = cmd.options.end();
      auto count = 0;
      for (auto it = cmd.options.begin(); it != cmd.options.end(); ++it) {
        auto& option = *it;
        if (!option.positional) continue;
        if (count == positionals.back()) {
          pos = it;
          positionals.back()++;
          break;
        }
      }
      if (pos == cmd.options.end())
        return cli_error("too many positional arguments");
      auto& option = *pos;
      option.set   = true;
      if (option.nargs > 0) {
        if (idx + (size_t)option.nargs > args.size())
          return cli_error("missing value for " + option.name);
        if (!parse_value(option, args, idx))
          return cli_error("bad value for " + option.name);
        idx += option.nargs - 1;
      } else if (option.nargs < 0) {
        if (!parse_value(option, args, idx))
          return cli_error("bad value for " + option.name);
        idx += args.size();
      } else {
        throw std::invalid_argument{"unsupported number of arguments"};
      }
    } else {
      auto pos = std::find_if(
          cmd.options.begin(), cmd.options.end(), [&arg](auto& option) {
            if (option.positional) return false;
            return arg == "--" + option.name || arg == "-" + option.alt;
          });
      if (pos == cmd.options.end()) return cli_error("unknown option " + arg);
      auto& option = *pos;
      option.set   = true;
      if (option.nargs == 0) {
        if (option.type != cli_type::boolean)
          throw std::invalid_argument{"unsupported flag type"};
        option.value.resize(1);
        option.value[0].integer = 1;
      } else if (option.nargs > 0) {
        if (idx + (size_t)option.nargs >= args.size())
          return cli_error("missing value for " + option.name);
        if (!parse_value(option, args, idx + 1))
          return cli_error("bad value for " + option.name);
        idx += option.nargs;
      } else {
        throw std::invalid_argument{"unsupported number of arguments"};
      }
    }
  }

  // check for help
  for (auto command_ptr : commands) {
    auto& command = *command_ptr;
    if (command.help) return true;
  }

  // check for required, set defaults and set references
  for (auto command_ptr : commands) {
    auto& command = *command_ptr;
    if (!command.commands.empty() && command.command.empty())
      return cli_error("command not set for " + command.name);
    if (command.set_command) command.set_command(command.command);
    for (auto& option : command.options) {
      if (option.req && !option.set)
        return cli_error("missing value for " + option.name);
      if (!option.set) option.value = option.def;
      if (option.set_value) {
        if (!option.set_value(option)) {
          return cli_error("bad value for " + option.name);
        }
      }
    }
  }

  // done
  return true;
}

bool parse_cli(cli_command& cli, int argc, const char** argv, string& error) {
  // validate names
  validate_names(cli);
  // prepare args
  auto args = vector<string>{argv + 1, argv + argc};
  // parse
  return parse_cli(cli, args, error);
}

void parse_cli(cli_command& cli, int argc, const char** argv) {
  auto error = string{};
  if (!parse_cli(cli, argc, argv, error)) {
    print_info("error: " + error);
    print_info("");
    print_info(get_usage(cli));
    exit(1);
  } else if (get_help(cli)) {
    print_info(get_usage(cli));
    exit(0);
  }
}

}  // namespace yocto
