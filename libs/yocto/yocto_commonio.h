//
// # Yocto/CommonIO: Utilities for writing command-line apps
//
// Yocto/CommonIO is a collection of utilities used in writing command-line
// applications, including parsing command line arguments, simple path
// manipulation, file lading and saving, and printing values, timers and
// progress bars.

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
using std::function;
using std::string;
using std::unordered_set;
using std::vector;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// PRINT/FORMATTING UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Print a message to the console
inline void print_info(const string& msg);
// Prints a message to the console and exit with an error.
inline void print_fatal(const string& msg);

// Timer that prints as scope end. Create with `print_timed` and print with
// `print_elapsed`.
struct print_timer {
  int64_t start_time = -1;
  ~print_timer();  // print time if scope ends
};
// Print traces for timing and program debugging
inline print_timer print_timed(const string& msg);
inline void        print_elapsed(print_timer& timer);

// Print progress
inline void print_progress(const string& message, int current, int total);

// Format duration string from nanoseconds
inline string format_duration(int64_t duration);
// Format a large integer number in human readable form
inline string format_num(uint64_t num);

}  // namespace yocto

// -----------------------------------------------------------------------------
// COMMAND LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// Initialize a command line parser.
struct cli_state;
inline cli_state make_cli(const string& cmd, const string& usage);
// parse arguments, checks for errors, and exits on error or help
inline void parse_cli(cli_state& cli, int argc, const char** argv);
// parse arguments and checks for errors
inline bool parse_cli(
    cli_state& cli, int argc, const char** argv, string& error);
// gets usage message
inline string get_usage(const cli_state& cli);
// gets whether help was invoked
inline bool get_help(const cli_state& cli);

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
inline string normalize_path(const string& filename);

// Get directory name (not including '/').
inline string path_dirname(const string& filename);

// Get extension (including '.').
inline string path_extension(const string& filename);

// Get filename without directory.
inline string path_filename(const string& filename);

// Get filename without directory and extension.
inline string path_basename(const string& filename);

// Joins paths
inline string path_join(const string& patha, const string& pathb);
inline string path_join(
    const string& patha, const string& pathb, const string& pathc);

// Replaces extensions
inline string replace_extension(const string& filename, const string& ext);

// Check if a file can be opened for reading.
inline bool path_exists(const string& filename);

// Check if a file is a directory
inline bool path_isdir(const string& filename);

// Check if a file is a file
inline bool path_isfile(const string& filename);

// List the contents of a directory
inline vector<string> list_directory(const string& filename);

// Create a directory and all missing parent directories if needed
inline bool make_directory(const string& dirname, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// FILE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a text file
inline bool load_text(const string& filename, string& str, string& error);
inline bool save_text(const string& filename, const string& str, string& error);

// Using directive
using byte = unsigned char;

// Load/save a binary file
inline bool load_binary(
    const string& filename, vector<byte>& data, string& error);
inline bool save_binary(
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
// PRINT/FORMATTING UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Print a message to the console
inline void print_info(const string& msg) { printf("%s\n", msg.c_str()); }
// Prints a messgae to the console and exit with an error.
inline void print_fatal(const string& msg) {
  printf("%s\n", msg.c_str());
  exit(1);
}

// get time in nanoseconds - useful only to compute difference of times
inline int64_t get_time_() {
  return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}

// Format duration string from nanoseconds
inline string format_duration(int64_t duration) {
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
inline string format_num(uint64_t num) {
  auto rem = num % 1000;
  auto div = num / 1000;
  if (div > 0) return format_num(div) + "," + std::to_string(rem);
  return std::to_string(rem);
}

// Print traces for timing and program debugging
inline print_timer print_timed(const string& msg) {
  printf("%s", msg.c_str());
  fflush(stdout);
  // print_info(fmt + " [started]", args...);
  return print_timer{get_time_()};
}
inline void print_elapsed(print_timer& timer) {
  if (timer.start_time < 0) return;
  printf(" in %s\n", format_duration(get_time_() - timer.start_time).c_str());
  timer.start_time = -1;
}
inline print_timer::~print_timer() { print_elapsed(*this); }

// Print progress
inline void print_progress(const string& message, int current, int total) {
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
// PATH UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Make a path from a utf8 string
inline std::filesystem::path make_path(const string& filename) {
  return std::filesystem::u8path(filename);
}

// Normalize path
inline string normalize_path(const string& filename) {
  return make_path(filename).generic_u8string();
}

// Get directory name (not including /)
inline string path_dirname(const string& filename) {
  return make_path(filename).parent_path().generic_u8string();
}

// Get extension (including .)
inline string path_extension(const string& filename) {
  return make_path(filename).extension().u8string();
}

// Get filename without directory.
inline string path_filename(const string& filename) {
  return make_path(filename).filename().u8string();
}

// Get filename without directory and extension.
inline string path_basename(const string& filename) {
  return make_path(filename).stem().u8string();
}

// Joins paths
inline string path_join(const string& patha, const string& pathb) {
  return (make_path(patha) / make_path(pathb)).generic_u8string();
}
inline string path_join(
    const string& patha, const string& pathb, const string& pathc) {
  return (make_path(patha) / make_path(pathb) / make_path(pathc))
      .generic_u8string();
}

// Replaces extensions
inline string replace_extension(const string& filename, const string& ext) {
  return make_path(filename).replace_extension(ext).u8string();
}

// Check if a file can be opened for reading.
inline bool path_exists(const string& filename) {
  return exists(make_path(filename));
}

// Check if a file is a directory
inline bool path_isdir(const string& filename) {
  return is_directory(make_path(filename));
}

// Check if a file is a file
inline bool path_isfile(const string& filename) {
  return is_regular_file(make_path(filename));
}

// List the contents of a directory
inline vector<string> list_directory(const string& filename) {
  auto entries = vector<string>{};
  for (auto entry : std::filesystem::directory_iterator(make_path(filename))) {
    entries.push_back(entry.path().generic_u8string());
  }
  return entries;
}

// Create a directory and all missing parent directories if needed
inline bool make_directory(const string& dirname, string& error) {
  if (path_exists(dirname)) return true;
  try {
    create_directories(make_path(dirname));
    return true;
  } catch (...) {
    error = dirname + ": cannot create directory";
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

// Cleanup
inline file_stream::~file_stream() {
  if (owned && fs) fclose(fs);
}

// Open a file
inline file_stream open_file(const string& filename, const string& mode) {
#ifdef _Win32
  auto path8 = std::filesystem::u8path(filename);
  auto wmode = std::wstring(mode.begin(), mode.end());
  auto fs    = _wfopen(path.c_str(), wmode.c_str());
#else
  auto fs = fopen(filename.c_str(), mode.c_str());
#endif
  return {filename, fs, true};
}

// Close a file
inline void close_file(file_stream& fs) {
  if (fs.owned && fs.fs) fclose(fs.fs);
  fs.filename = "";
  fs.fs       = nullptr;
  fs.owned    = false;
}

// Read a line of text
inline bool read_line(file_stream& fs, char* buffer, size_t size) {
  return fgets(buffer, (int)size, fs.fs);
}

// Write text to a file
inline bool write_text(file_stream& fs, const string& str) {
  return fprintf(fs.fs, "%s", str.c_str()) >= 0;
}

// Read data from a file
template <typename T>
inline bool read_value(file_stream& fs, T& buffer) {
  return fread(&buffer, sizeof(T), 1, fs.fs) == 1;
}

// Write data from a file
template <typename T>
inline bool write_value(file_stream& fs, const T& buffer) {
  return fwrite(&buffer, sizeof(T), 1, fs.fs) == 1;
}

// Read data from a file
template <typename T>
inline bool read_values(file_stream& fs, T* buffer, size_t count) {
  return fread(buffer, sizeof(T), count, fs.fs) == count;
}

// Write data from a file
template <typename T>
inline bool write_values(file_stream& fs, const T* buffer, size_t count) {
  return fwrite(buffer, sizeof(T), count, fs.fs) == count;
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
  if (!write_value(fs, value)) return false;
  return true;
}

// Opens a file with a utf8 file name
inline FILE* fopen_utf8(const char* filename, const char* mode) {
#ifdef _Win32
  auto path8 = std::filesystem::u8path(filename);
  auto wmode = std::wstring(string{mode}.begin(), string{mode}.end());
  return _wfopen(path.c_str(), wmode.c_str());
#else
  return fopen(filename, mode);
#endif
}

// Load a text file
inline bool load_text(const string& filename, string& str, string& error) {
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
inline bool save_text(
    const string& filename, const string& str, string& error) {
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
inline bool load_binary(
    const string& filename, vector<byte>& data, string& error) {
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
inline bool save_binary(
    const string& filename, const vector<byte>& data, string& error) {
  auto fs = open_file(filename, "wb");
  if (!fs) {
    error = filename + ": file not found";
    return false;
  }
  if (!write_values(fs, data.data(), data.size()) != data.size()) {
    error = filename + ": write error";
    return false;
  }
  return true;
}

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
  char buf[256];
  snprintf(buf, sizeof(buf), "%g", value);
  str += buf;
}
inline void format_value(string& str, double value) {
  char buf[256];
  snprintf(buf, sizeof(buf), "%g", value);
  str += buf;
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
  if (!write_text(fs, str)) return false;
  return true;
}
template <typename T>
inline bool format_value(file_stream& fs, const T& value) {
  auto str = ""s;
  format_value(str, value);
  if (!write_text(fs, str)) return false;
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF COMMAND-LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// Command line parser data. All data should be considered private.
struct cli_option {
  string                                name          = "";
  bool                                  req           = false;
  int                                   nargs         = 0;
  function<bool(const vector<string>&)> parse_and_set = {};
};
struct cli_state {
  string             name            = "";
  string             usage           = "";
  vector<cli_option> options         = {};
  string             usage_options   = "";
  string             usage_arguments = "";
  bool               help            = false;
};

// initialize a command line parser
inline cli_state make_cli(const string& cmd, const string& usage) {
  auto cli  = cli_state{};
  cli.name  = cmd;
  cli.usage = usage;
  add_option(cli, "--help/--no-help", cli.help, "Print usage.");
  return cli;
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
struct cli_is_vector : std::false_type {};
template <typename T>
struct cli_is_vector<std::vector<T>> : std::true_type {};
template <typename T, typename A>
struct cli_is_vector<std::vector<T, A>> : std::true_type {};
template <class T>
inline constexpr bool cli_is_vector_v = cli_is_vector<T>::value;

template <typename T>
inline string cli_type_name() {
  if constexpr (std::is_same_v<T, string>) return "<string>";
  if constexpr (std::is_same_v<T, bool>) return "";
  if constexpr (std::is_integral_v<T>) return "<integer>";
  if constexpr (std::is_floating_point_v<T>) return "<number>";
  if constexpr (std::is_enum_v<T>) return "<enum>";
  if constexpr (cli_is_vector_v<T>)
    return "<[" + cli_type_name<typename T::value_type>() + "]>";
  return "<value>";
}

template <typename T>
inline string cli_to_string(const T& value, const vector<string>& choices) {
  if constexpr (std::is_same_v<T, string>) return value;
  if constexpr (std::is_same_v<T, bool>) return value ? "true" : "false";
  if constexpr (std::is_integral_v<T>)
    return choices.empty() ? std::to_string(value) : choices.at(value);
  if constexpr (std::is_floating_point_v<T>) return std::to_string(value);
  if constexpr (std::is_enum_v<T>)
    return choices.empty() ? std::to_string((int)value)
                           : choices.at((int)value);
  if constexpr (cli_is_vector_v<T>) {
    auto def = string{"["};
    for (auto i = 0; i < value.size(); i++)
      def += (i ? "," : "") + cli_to_string(value[i], choices);
    return def;
  }
  throw std::invalid_argument{"unsupported type"};
}

template <typename T>
inline int cli_nargs() {
  if constexpr (std::is_same_v<T, string>) return 1;
  if constexpr (std::is_same_v<T, bool>) return 0;
  if constexpr (std::is_integral_v<T>) return 1;
  if constexpr (std::is_floating_point_v<T>) return 1;
  if constexpr (std::is_enum_v<T>) return 1;
  if constexpr (cli_is_vector_v<T>) return -1;
  throw std::invalid_argument{"unsupported type"};
}

template <typename T>
inline bool parse_cli_value(
    const vector<string>& args, T& value, const vector<string>& choices) {
  if (!choices.empty()) {
    for (auto& arg : args) {
      if (std::find(choices.begin(), choices.end(), arg) == choices.end())
        return false;
    }
  }
  if constexpr (std::is_same_v<T, string>) {
    if (args.size() != 1) return false;
    value = args[0];
    return true;
  } else if constexpr (std::is_same_v<T, bool>) {
    if (args.size() != 1) return false;
    if (args[0] == "true" || args[0] == "1") {
      value = true;
      return true;
    } else if (args[0] == "false" || args[0] == "0") {
      value = false;
      return true;
    } else {
      return false;
    }
  } else if constexpr (std::is_integral_v<T>) {
    if (args.size() != 1) return false;
    if (choices.empty()) {
      auto end = (char*)nullptr;
      value    = (int)strtol(args[0].c_str(), &end, 10);
      return end != nullptr;
    } else {
      value = (T)(
          std::find(choices.begin(), choices.end(), args[0]) - choices.begin());
      return true;
    }
  } else if constexpr (std::is_same_v<T, float>) {
    if (args.size() != 1) return false;
    auto end = (char*)nullptr;
    value    = strtof(args[0].c_str(), &end);
    return end != nullptr;
  } else if constexpr (std::is_same_v<T, double>) {
    if (args.size() != 1) return false;
    auto end = (char*)nullptr;
    value    = strtod(args[0].c_str(), &end);
    return end != nullptr;
  } else if constexpr (std::is_enum_v<T>) {
    auto ivalue = 0;
    if (!parse_cli_value(args, ivalue, choices)) return false;
    value = (T)ivalue;
    return true;
  } else if constexpr (std::is_same_v<T, vector<string>>) {
    value = args;
    return true;
  } else if constexpr (cli_is_vector_v<T>) {
    value.clear();
    for (auto& arg : args)
      if (!parse_cli_value({arg}, value.emplace_back())) return false;
    return true;
  } else {
    throw std::invalid_argument{"unsupported type"};
  }
}

template <typename T>
inline void add_cli_option(cli_state& cli, const string& name, T& value,
    const string& usage, bool req, const vector<string>& choices) {
  // check for errors
  auto used = unordered_set<string>{};
  for (auto& option : cli.options) {
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

  // help message
  auto line = "  " + name + " " + cli_type_name<T>();
  while (line.size() < 32) line += " ";
  line += usage;
  line += !req ? " [" + cli_to_string(value, choices) + "]\n" : " [required]\n";
  if (!choices.empty()) {
    line += "    with choices: ";
    auto len = 16;
    for (auto& choice : choices) {
      if (len + choice.size() + 2 > 78) {
        line += "\n                 ";
        len = 16;
      }
      line += choice + ", ";
      len += choice.size() + 2;
    }
    line = line.substr(0, line.size() - 2);
    line += "\n";
  }
  if (name.find("-") == 0) {
    cli.usage_options += line;
  } else {
    cli.usage_arguments += line;
  }
  // add option
  cli.options.push_back({name, req, cli_nargs<T>(),
      [&value, choices](const vector<string>& args) -> bool {
        return parse_cli_value(args, value, choices);
      }});
}

template <typename T>
inline void add_option(cli_state& cli, const string& name, T& value,
    const string& usage, bool req) {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  return add_cli_option(cli, name, value, usage, req, {});
}
template <typename T>
inline void add_option(cli_state& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, bool req) {
  static_assert(
      std::is_same_v<T, string> || std::is_integral_v<T> || std::is_enum_v<T>,
      "unsupported type");
  return add_cli_option(cli, name, value, usage, req, choices);
}
template <typename T>
inline void add_option(cli_state& cli, const string& name, vector<T>& value,
    const string& usage, bool req) {
  static_assert(std::is_same_v<T, string> || std::is_same_v<T, bool> ||
                    std::is_integral_v<T> || std::is_floating_point_v<T> ||
                    std::is_enum_v<T>,
      "unsupported type");
  return add_cli_option(cli, name, value, usage, req, {});
}

inline bool get_help(const cli_state& cli) { return cli.help; }

inline string get_usage(const cli_state& cli) {
  auto message = string{};
  message +=
      "usage: " + cli.name + (cli.usage_options.empty() ? "" : " [options]") +
      (cli.usage_arguments.empty() ? "" : " <arguments>") + cli.usage + "\n\n";
  if (!cli.usage_options.empty())
    message += "options:\n" + cli.usage_options + "\n";
  if (!cli.usage_options.empty())
    message += "arguments:\n" + cli.usage_arguments + "\n";
  return message;
}

inline bool parse_cli(
    cli_state& cli, int argc, const char** argv, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  // prepare args
  auto args = vector<string>{argv + 1, argv + argc};
  // parse options
  for (auto& option : cli.options) {
    if (option.name[0] != '-') continue;
    auto set    = false;
    auto values = vector<string>{};
    for (auto& name : split_cli_names(option.name)) {
      if (std::find(args.begin(), args.end(), name) == args.end()) continue;
      auto pos = std::find(args.begin(), args.end(), name) - args.begin();
      args.erase(args.begin() + pos);
      if (option.nargs == 0) {
        values = {name.find("--no-") == string::npos ? "true" : "false"};
        set    = true;
      } else if (option.nargs > 0) {
        if (pos + option.nargs > args.size())
          return cli_error("missing value for " + name);
        values = {args.begin() + pos, args.begin() + pos + option.nargs};
        set    = true;
        args.erase(args.begin() + pos, args.begin() + pos + option.nargs);
      } else {
        throw std::invalid_argument{"unsupported number of arguments"};
      }
    }
    if (set) {
      if (!option.parse_and_set(values))
        return cli_error("bad value for " + option.name);
    } else {
      if (option.req) return cli_error("missing value for " + option.name);
    }
  }
  // check unknown options
  for (auto& arg : args) {
    if (arg.find("-") == 0) return cli_error("unknown option " + arg);
  }
  // parse positional
  for (auto& option : cli.options) {
    if (option.name[0] == '-') continue;
    auto set    = false;
    auto values = vector<string>{};
    if (args.empty()) {
      if (option.req) return cli_error("missing value for " + option.name);
    } else if (option.nargs < 0) {
      values = args;
      set    = true;
      args.clear();
    } else if (option.nargs > 0) {
      if (option.nargs > args.size())
        return cli_error("missing value for " + option.name);
      values = {args.begin(), args.begin() + option.nargs};
      args.erase(args.begin(), args.begin() + option.nargs);
      set = true;
    } else {
      throw std::invalid_argument{"unsupported number of arguments"};
    }
    if (set) {
      if (!option.parse_and_set(values))
        return cli_error("bad value for " + option.name);
    } else {
      if (option.req) return cli_error("missing value for " + option.name);
    }
  }
  // check remaining
  if (!args.empty()) return cli_error("mismatched value for " + args.front());
  // done
  return true;
}

inline void parse_cli(cli_state& cli, int argc, const char** argv) {
  auto error = string{};
  if (!parse_cli(cli, argc, argv, error)) {
    print_info("error: " + error);
    print_info("");
    print_info(get_usage(cli));
    exit(1);
  } else if (cli.help) {
    print_info(get_usage(cli));
    exit(0);
  }
}

}  // namespace yocto

#endif
