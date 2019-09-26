//
// # Yocto/Utils: Tiny collection of utilities to support Yocto/GL
//
//
// Yocto/Utils is a collection of utilities used in writing other Yocto/GL
// libraries and example applications. We support printing and parsing builting
// and Yocto/Math values, parsing command line arguments, simple path
// manipulation, file lading/saving and basic concurrency utilities.
//
//
// ## Printing values
//
// Use `print_info()` to print a message, `print_fatal()` to print and exit,
// and `print_timed()` to use a RIIA timer. Several overloads of `to_string()`
// are provided for both the basic types and Yocto/Math types.
//
//
// ## Python-like iterators and collection helpers
//
// This library includes a set of functions to help use C++ collections with
// more ease, inspired by Python. All functions and operators are defined in
// the yocto namespace so they will not affect the code outside. But within
// the Yocto/GL collection they are the best way to do this.
//
// 1. use `range()` to iterato over an integer sequence
// 2. use `enumerate()` to iteratare over a vector and number its elements
// 3. use opeartors + to either concatenate two vectors or a vector and an
//    element
// 4. use operators += to append an element or a vector to a given vector
//
//
// ## Command-Line Parsing
//
// We provide a simple, immediate-mode, command-line parser. The parser
// works in an immediate-mode manner since it reads each value as you call each
// function, rather than building a data structure and parsing offline. We
// support option and position arguments, automatic help generation, and
// error checking.
//
// 1. initialize the parser with `auto cli = make_cli(argc, argv, help)`
// 2. add options with `add_cli_option(cli, name, value, usage, req)`
//    - if name starts with '--' or '-' then it is an option
//    - otherwise it is a positional argument
//    - options and arguments may be intermixed
//    - the type of each option is determined by the passed reference `value`
//    - `req` indicates whether an option or argument is required or not
// 3. parse options with `parse_cli(cli, argc, argv)`
//    - if an error occurrs, the parser prints a usage message and returns false
//
//
// ## Path manipulation
//
// We define a few path manipulation utilities to split and join path
// components.
//
// 1. Get paths components with `get_dirname()`, `get_filename()` and
//   `get_extension()`
// 2. Replace the extension with `replace_path_extension()`
// 3. check if a file exists with `exists_file()`
//
//
// ## File IO
//
// 1. load and save text files with `load_text()` and `save_text()`
// 2. load and save binary files with `load_binary()` and `save_binary()`
//
//
// ## Concurrency utilities
//
// C++ has very basic supprt for concurrency and most of it is still platform
// dependent. We provide here very basic support for concurrency utlities
// built on top of C++ low-level threading and synchronization.
//
// 1. use `concurrent_queue()` for communicationing values between threads
// 2. use `parallel_for()` for basic parallel for loops
//
//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
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

#ifndef _YOCTO_UTILS_H_
#define _YOCTO_UTILS_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_math.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <deque>
#include <mutex>
#include <string>
#include <thread>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

using std::deque;
using std::thread;

}  // namespace yocto

// -----------------------------------------------------------------------------
// CONVERSION TO STRING
// -----------------------------------------------------------------------------
namespace yocto {

// Conversion to string for basic and Yocto/Math types
inline string to_string(int value) { return std::to_string(value); }
inline string to_string(float value) { return std::to_string(value); }
inline string to_string(double value) { return std::to_string(value); }
inline string to_string(const string& value) { return value; }
inline string to_string(const char* value) { return value; }
inline string to_string(const vec2f& value) {
  return to_string(value.x) + " " + to_string(value.y);
}
inline string to_string(const vec3f& value) {
  return to_string(value.x) + " " + to_string(value.y) + " " +
         to_string(value.z);
}
inline string to_string(const vec4f& value) {
  return to_string(value.x) + " " + to_string(value.y) + " " +
         to_string(value.z) + " " + to_string(value.w);
}
inline string to_string(const vec2i& value) {
  return to_string(value.x) + " " + to_string(value.y);
}
inline string to_string(const vec3i& value) {
  return to_string(value.x) + " " + to_string(value.y) + " " +
         to_string(value.z);
}
inline string to_string(const vec4i& value) {
  return to_string(value.x) + " " + to_string(value.y) + " " +
         to_string(value.z) + " " + to_string(value.w);
}
inline string to_string(const mat2f& value) {
  return to_string(value.x) + " " + to_string(value.y);
}
inline string to_string(const mat3f& value) {
  return to_string(value.x) + " " + to_string(value.y) + " " +
         to_string(value.z);
}
inline string to_string(const mat4f& value) {
  return to_string(value.x) + " " + to_string(value.y) + " " +
         to_string(value.z) + " " + to_string(value.w);
}
inline string to_string(const frame2f& value) {
  return to_string(value.x) + " " + to_string(value.y) + " " +
         to_string(value.o);
}
inline string to_string(const frame3f& value) {
  return to_string(value.x) + " " + to_string(value.y) + " " +
         to_string(value.z) + " " + to_string(value.o);
}
inline string to_string(const ray2f& value) {
  return to_string(value.o) + " " + to_string(value.d) + " " +
         to_string(value.tmin) + " " + to_string(value.tmax);
}
inline string to_string(const ray3f& value) {
  return to_string(value.o) + " " + to_string(value.d) + " " +
         to_string(value.tmin) + " " + to_string(value.tmax);
}
inline string to_string(const bbox2f& value) {
  return to_string(value.min) + " " + to_string(value.max);
}
inline string to_string(const bbox3f& value) {
  return to_string(value.min) + " " + to_string(value.max);
}

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
inline int64_t get_time() {
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
  sprintf(buffer, "%02d:%02d:%02d.%03d", hours, mins, secs, msecs);
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
inline auto print_timed(const string& msg) {
  struct scoped_timer {
    int64_t start_time = -1;
    ~scoped_timer() {
      printf(" in %s\n", format_duration(get_time() - start_time).c_str());
    }
  };
  printf("%s", msg.c_str());
  fflush(stdout);
  // print_info(fmt + " [started]", args...);
  return scoped_timer{get_time()};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PYTHON-LIKE ITERATORS
// -----------------------------------------------------------------------------
namespace yocto {

// Range object to support Python-like iteration. Use with `range()`.
struct range_helper {
  struct iterator {
    int       pos = 0;
    iterator& operator++() {
      pos++;
      return *this;
    }
    bool operator!=(const iterator& other) const { return pos != other.pos; }
    int  operator*() const { return pos; }
  };
  int      begin_ = 0, end_ = 0;
  iterator begin() const { return {begin_}; }
  iterator end() const { return {end_}; }
};

// Python `range()` equivalent. Construct an object to iterate over a sequence.
inline auto range(int min, int max) { return range_helper{min, max}; }
inline auto range(int max) { return range(0, max); }

// Enumerate object to support Python-like enumeration. Use with `enumerate()`.
template <typename T>
struct enumerate_helper {
  struct iterator {
    T*        data = nullptr;
    int       pos  = 0;
    iterator& operator++() {
      pos++;
      return *this;
    }
    bool operator!=(const iterator& other) const { return pos != other.pos; }
    pair<int&, T&> operator*() const { return {pos, *(data + pos)}; }
  };
  T*       data = nullptr;
  int      size = 0;
  iterator begin() const { return {data, 0}; }
  iterator end() const { return {data, size}; }
};

// Python `enumerate()` equivalent. Construct an object that iteraterates over a
// sequence of elements and numbers them.
template <typename T>
inline auto enumerate(const vector<T>& vals) {
  return enumerate_helper<const T>{vals.data(), vals.size()};
}
template <typename T>
inline auto enumerate(vector<T>& vals) {
  return enumerate_helper<T>{vals.data(), vals.size()};
}

// Vector append and concatenation
template <typename T>
inline vector<T>& operator+=(vector<T>& a, const vector<T>& b) {
  a.insert(a.end(), b.begin(), b.end());
  return a;
}
template <typename T>
inline vector<T>& operator+=(vector<T>& a, const T& b) {
  a.push_back(b);
  return a;
}
template <typename T>
inline vector<T> operator+(const vector<T>& a, const vector<T>& b) {
  auto c = a;
  return c += b;
}
template <typename T>
inline vector<T> operator+(const vector<T>& a, const T& b) {
  auto c = a;
  return c += b;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMMEDIATE-MODE COMMAND LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// Initialize a command line parser.
struct cli_state;
inline cli_state make_cli(const string& cmd, const string& usage);
// check if any error occurred and exit appropriately
inline bool parse_cli(cli_state& cli, int argc, const char** argv);

// Parse an int, float, string, and bool option or positional argument.
// Options's names starts with "--" or "-", otherwise they are arguments.
// The library support using many names for the same option/argument
// separate by commas. Boolean flags are indicated with a pair of names
// "--name/--no-name", so that we have both options available.
inline void add_cli_option(cli_state& cli, const string& name, string& value,
    const string& usage, bool req = false);
inline void add_cli_option(cli_state& cli, const string& name, int& value,
    const string& usage, bool req = false);
inline void add_cli_option(cli_state& cli, const string& name, float& value,
    const string& usage, bool req = false);
inline void add_cli_option(cli_state& cli, const string& name, bool& value,
    const string& usage, bool req = false);
// Parse an enum
inline void add_cli_option(cli_state& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, bool req = false);
// Parse all arguments left on the command line.
inline void add_cli_option(cli_state& cli, const string& name,
    vector<string>& value, const string& usage, bool req = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// PATH UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// These utilities are here only for backward compatibility. They should be
// considered deprecated.

// Utility to normalize a path
inline string normalize_path(const string& filename_) {
  auto filename = filename_;
  for (auto& c : filename)
    if (c == '\\') c = '/';
  if (filename.size() > 1 && filename[0] == '/' && filename[1] == '/') {
    throw std::invalid_argument("absolute paths are not supported");
    return filename_;
  }
  if (filename.size() > 3 && filename[1] == ':' && filename[2] == '/' &&
      filename[3] == '/') {
    throw std::invalid_argument("absolute paths are not supported");
    return filename_;
  }
  auto pos = (size_t)0;
  while ((pos = filename.find("//")) != filename.npos)
    filename = filename.substr(0, pos) + filename.substr(pos + 1);
  return filename;
}

// Get directory name (including '/').
inline string get_dirname(const string& filename_) {
  auto filename = normalize_path(filename_);
  auto pos      = filename.rfind('/');
  if (pos == string::npos) return "";
  return filename.substr(0, pos + 1);
}

// Get extension (not including '.').
inline string get_extension(const string& filename_) {
  auto filename = normalize_path(filename_);
  auto pos      = filename.rfind('.');
  if (pos == string::npos) return "";
  return filename.substr(pos);
}

// Get filename without directory.
inline string get_filename(const string& filename_) {
  auto filename = normalize_path(filename_);
  auto pos      = filename.rfind('/');
  if (pos == string::npos) return filename;
  return filename.substr(pos + 1);
}

// Get extension.
inline string get_noextension(const string& filename_) {
  auto filename = normalize_path(filename_);
  auto pos      = filename.rfind('.');
  if (pos == string::npos) return filename;
  return filename.substr(0, pos);
}

// Get filename without directory and extension.
inline string get_basename(const string& filename) {
  return get_noextension(get_filename(filename));
}

// Replaces extensions
inline string replace_extension(const string& filename, const string& ext) {
  return get_noextension(filename) + ext;
}

// Check if a file can be opened for reading.
inline bool exists_file(const string& filename) {
  auto fs = fopen(filename.c_str(), "r");
  if (fs) {
    fclose(fs);
    return true;
  } else {
    return false;
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// FILE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load a text file
inline void load_text(const string& filename, string& str) {
  // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
  auto fs = fopen(filename.c_str(), "rt");
  if (!fs) throw std::runtime_error("cannot open file " + filename);
  fseek(fs, 0, SEEK_END);
  auto length = ftell(fs);
  fseek(fs, 0, SEEK_SET);
  str.resize(length);
  if (fread(str.data(), 1, length, fs) != length) {
    fclose(fs);
    throw std::runtime_error("cannot read file " + filename);
  }
  fclose(fs);
}

// Save a text file
inline void save_text(const string& filename, const string& str) {
  auto fs = fopen(filename.c_str(), "wt");
  if (!fs) throw std::runtime_error("cannot open file " + filename);
  if (fprintf(fs, "%s", str.c_str()) < 0) {
    fclose(fs);
    throw std::runtime_error("cannot write file " + filename);
  }
  fclose(fs);
}

// Load a binary file
inline void load_binary(const string& filename, vector<byte>& data) {
  // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
  auto fs = fopen(filename.c_str(), "rb");
  if (!fs) throw std::runtime_error("cannot open file " + filename);
  fseek(fs, 0, SEEK_END);
  auto length = ftell(fs);
  fseek(fs, 0, SEEK_SET);
  data.resize(length);
  if (fread(data.data(), 1, length, fs) != length) {
    fclose(fs);
    throw std::runtime_error("cannot read file " + filename);
  }
  fclose(fs);
}

// Save a binary file
inline void save_binary(const string& filename, const vector<byte>& data) {
  auto fs = fopen(filename.c_str(), "wb");
  if (!fs) throw std::runtime_error("cannot open file " + filename);
  if (fwrite(data.data(), 1, data.size(), fs) != data.size()) {
    fclose(fs);
    throw std::runtime_error("cannot write file " + filename);
  }
  fclose(fs);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// CONCURRENCY UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// a simple concurrent queue that locks at every call
template <typename T>
struct concurrent_queue {
  concurrent_queue() {}
  concurrent_queue(const concurrent_queue& other) = delete;
  concurrent_queue& operator=(const concurrent_queue& other) = delete;

  bool empty() {
    std::lock_guard<std::mutex> lock(mutex);
    return queue.empty();
  }
  void clear() {
    std::lock_guard<std::mutex> lock(mutex);
    queue.clear();
  }
  void push(const T& value) {
    std::lock_guard<std::mutex> lock(mutex);
    queue.push_back(value);
  }
  bool try_pop(T& value) {
    std::lock_guard<std::mutex> lock(mutex);
    if (queue.empty()) return false;
    value = queue.front();
    queue.pop_front();
    return true;
  }

 private:
  std::mutex    mutex;
  std::deque<T> queue;
};

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename Func>
inline void parallel_for(int begin, int end, const Func& func,
    std::atomic<bool>* cancel = nullptr, bool serial = false) {
  if (serial) {
    for (auto idx = begin; idx < end; idx++) {
      if (cancel && *cancel) break;
      func(idx);
    }
  } else {
    auto             threads  = vector<thread>{};
    auto             nthreads = thread::hardware_concurrency();
    std::atomic<int> next_idx(begin);
    for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
      threads.emplace_back([&func, &next_idx, cancel, end]() {
        while (true) {
          if (cancel && *cancel) break;
          auto idx = next_idx.fetch_add(1);
          if (idx >= end) break;
          func(idx);
        }
      });
    }
    for (auto& t : threads) t.join();
  }
}

template <typename Func>
inline void parallel_for(int num, const Func& func,
    std::atomic<bool>* cancel = nullptr, bool serial = false) {
  parallel_for(0, num, func, cancel, serial);
}

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes a reference to a `T`.
template <typename T, typename Func>
inline void parallel_foreach(vector<T>& values, const Func& func,
    std::atomic<bool>* cancel = nullptr, bool serial = false) {
  parallel_for(
      0, (int)values.size(), [&func, &values](int idx) { func(values[idx]); },
      cancel, serial);
}
template <typename T, typename Func>
inline void parallel_foreach(const vector<T>& values, const Func& func,
    std::atomic<bool>* cancel = nullptr, bool serial = false) {
  parallel_for(
      0, (int)values.size(), [&func, &values](int idx) { func(values[idx]); },
      cancel, serial);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF COMMAND-LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// Command line parser data. All data should be considered private.
enum struct cli_type {
  // clang-format off
  string_, int_, float_, bool_, flag_, string_vector_, enum_
  // clang-format on
};
struct cmdline_option {
  string         name    = "";
  string         usage   = "";
  cli_type       type    = cli_type::string_;
  void*          value   = nullptr;
  bool           req     = false;
  bool           set     = false;
  vector<string> choices = {};
};
struct cli_state {
  string                 name            = "";
  string                 usage           = "";
  vector<cmdline_option> options         = {};
  string                 usage_options   = "";
  string                 usage_arguments = "";
};

// initialize a command line parser
inline cli_state make_cli(const string& cmd, const string& usage) {
  auto parser  = cli_state{};
  parser.name  = cmd;
  parser.usage = usage;
  return parser;
}

inline vector<string> split_cli_names(const string& name_) {
  auto name  = name_;
  auto split = vector<string>{};
  if (name.empty()) throw std::runtime_error("option name cannot be empty");
  if (name.find_first_of(" \t\r\n") != string::npos)
    throw std::runtime_error("option name cannot contain whitespaces");
  while (name.find_first_of(",/") != string::npos) {
    auto pos = name.find_first_of(",/");
    if (pos > 0) split.push_back(name.substr(0, pos));
    name = name.substr(pos + 1);
  }
  if (!name.empty()) split.push_back(name);
  if (split.empty()) throw std::runtime_error("option name cannot be empty");
  for (auto& name : split)
    if ((split[0][0] == '-') != (name[0] == '-'))
      throw std::runtime_error("inconsistent option names for " + name);
  return split;
}

inline void add_cli_option(cli_state& cli, const string& name, cli_type type,
    void* value, const string& usage, bool req, const vector<string>& choices) {
  static auto type_name = unordered_map<cli_type, string>{
      {cli_type::string_, "<string>"},
      {cli_type::int_, "<int>"},
      {cli_type::float_, "<float>"},
      {cli_type::bool_, "<true/false>"},
      {cli_type::flag_, ""},
      {cli_type::string_vector_, "<[string]>"},
      {cli_type::enum_, "<string>"},
  };
  // help message
  auto line = "  " + name + " " + type_name.at(type);
  while (line.size() < 32) line += " ";
  line += usage;
  if (!req) {
    line += " [";
    switch (type) {
      case cli_type::string_: line += *(string*)value; break;
      case cli_type::int_: line += std::to_string(*(int*)value); break;
      case cli_type::float_: line += std::to_string(*(float*)value); break;
      case cli_type::bool_: line += *(bool*)value ? "true" : "false"; break;
      case cli_type::flag_: line += *(bool*)value ? "true" : "false"; break;
      case cli_type::enum_: line += choices.at(*(int*)value); break;
      case cli_type::string_vector_: {
        for (auto i = 0; i < (*(vector<string>*)value).size(); i++) {
          if (i) line += ",";
          line += (*(vector<string>*)value)[i];
        }
      } break;
      default: throw std::runtime_error("unknown type");
    }
    line += "]";
  } else {
    line += " [required]";
  }
  line += "\n";
  if (name.find("-") == 0) {
    cli.usage_options += line;
  } else {
    cli.usage_arguments += line;
  }
  // add option
  cli.options.push_back(
      cmdline_option{name, usage, type, value, req, false, choices});
}

inline void add_cli_option(cli_state& cli, const string& name, string& value,
    const string& usage, bool req) {
  return add_cli_option(cli, name, cli_type::string_, &value, usage, req, {});
}
inline void add_cli_option(cli_state& cli, const string& name, int& value,
    const string& usage, bool req) {
  return add_cli_option(cli, name, cli_type::int_, &value, usage, req, {});
}
inline void add_cli_option(cli_state& cli, const string& name, float& value,
    const string& usage, bool req) {
  return add_cli_option(cli, name, cli_type::float_, &value, usage, req, {});
}
inline void add_cli_option(cli_state& cli, const string& name, bool& value,
    const string& usage, bool req) {
  return add_cli_option(cli, name, cli_type::flag_, &value, usage, req, {});
}
inline void add_cli_option(cli_state& cli, const string& name,
    vector<string>& value, const string& usage, bool req) {
  return add_cli_option(
      cli, name, cli_type::string_vector_, &value, usage, req, {});
}
inline void add_flag(cli_state& cli, const string& name, bool& value,
    const string& usage, bool req) {
  return add_cli_option(cli, name, cli_type::flag_, &value, usage, req, {});
}
inline void add_cli_option(cli_state& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, bool req) {
  return add_cli_option(
      cli, name, cli_type::enum_, &value, usage, req, choices);
}

inline bool print_cli_help(cli_state& cli, const string& error) {
  if (error != "") printf("error: %s\n\n", error.c_str());
  printf("usage: %s%s%s%s\n\n", cli.name.c_str(),
      cli.usage_options.empty() ? "" : " [options]",
      cli.usage_arguments.empty() ? "" : " <arguments>", cli.usage.c_str());
  if (!cli.usage_options.empty())
    printf("options:\n%s\n", cli.usage_options.c_str());
  if (!cli.usage_options.empty())
    printf("arguments:\n%s\n", cli.usage_arguments.c_str());
  return error.empty();
}

inline bool parse_cmdline_value(const string& str, int& value) {
  auto end = (char*)nullptr;
  value    = (int)strtol(str.c_str(), &end, 10);
  return end != nullptr;
}
inline bool parse_cmdline_value(const string& str, float& value) {
  auto end = (char*)nullptr;
  value    = strtof(str.c_str(), &end);
  return end != nullptr;
}
inline bool parse_cmdline_value(const string& str, bool& value) {
  if (str == "true" || str == "1") {
    value = true;
    return true;
  } else if (str == "false" || str == "0") {
    value = false;
    return true;
  } else {
    return false;
  }
}

inline bool parse_cli(cli_state& cli, int argc, const char** argv) {
  // check for errors
  auto used = unordered_map<string, int>{};
  for (auto& option : cli.options) {
    if (option.name.empty()) throw std::runtime_error("name cannot be empty");
    auto names = split_cli_names(option.name);
    if (names.empty()) throw std::runtime_error("name cannot be empty");
    for (auto& name : names) {
      if (used.find(name) != used.end())
        throw std::runtime_error("option name " + name + " already in use");
      used[name] = 1;
      if ((name[0] == '-') != (option.name[0] == '-'))
        throw std::runtime_error("incosisten option type for " + name);
    }
  }
  // prepare args
  auto args = vector<string>{argv + 1, argv + argc};
  // parse options
  for (auto& option : cli.options) {
    if (option.name[0] != '-') continue;
    for (auto& name : split_cli_names(option.name)) {
      auto pos = std::find(args.begin(), args.end(), name) - args.begin();
      if (pos >= args.size()) continue;
      if (option.type == cli_type::flag_) {
        *(bool*)option.value = name.find("--no-") == string::npos;
        option.set           = true;
        args.erase(args.begin() + pos);
      } else {
        if (pos + 1 >= args.size())
          return print_cli_help(cli, "missing value for " + name);
        auto value = args[pos + 1];
        args.erase(args.begin() + pos, args.begin() + pos + 2);
        if (option.type == cli_type::string_) {
          *(string*)option.value = value;
          option.set             = true;
        } else if (option.type == cli_type::int_) {
          if (!parse_cmdline_value(value, *(int*)option.value))
            return print_cli_help(cli, "incorrect value for " + name);
          option.set = true;
        } else if (option.type == cli_type::float_) {
          if (!parse_cmdline_value(value, *(float*)option.value))
            return print_cli_help(cli, "incorrect value for " + name);
          option.set = true;
        } else if (option.type == cli_type::bool_) {
          if (!parse_cmdline_value(value, *(bool*)option.value))
            return print_cli_help(cli, "incorrect value for " + name);
          option.set = true;
        } else if (option.type == cli_type::enum_) {
          auto pos = std::find(
              option.choices.begin(), option.choices.end(), value);
          if (pos == option.choices.end())
            return print_cli_help(cli, "incorrect value for " + name);
          else
            *(int*)option.value = (int)(pos - option.choices.begin());
          option.set = true;
        } else {
          throw std::runtime_error("unsupported type");
        }
      }
    }
    if (option.req && !option.set) {
      print_cli_help(cli, "missing value for " + option.name);
    }
  }
  // check unknown options
  for (auto& arg : args) {
    if (arg.find("-") == 0) return print_cli_help(cli, "unknown option " + arg);
  }
  // parse positional
  for (auto& option : cli.options) {
    if (option.name[0] == '-') continue;
    if (args.empty()) {
      if (option.req)
        return print_cli_help(cli, "missing value for " + option.name);
    } else if (option.type == cli_type::string_vector_) {
      *(vector<string>*)option.value = args;
      option.set                     = true;
      args.clear();
    } else {
      auto value = args.front();
      args.erase(args.begin());
      if (option.type == cli_type::string_) {
        *(string*)option.value = value;
        option.set             = true;
      } else if (option.type == cli_type::int_) {
        if (!parse_cmdline_value(value, *(int*)option.value))
          return print_cli_help(cli, "incorrect value for " + option.name);
        option.set = true;
      } else if (option.type == cli_type::float_) {
        if (!parse_cmdline_value(value, *(float*)option.value))
          return print_cli_help(cli, "incorrect value for " + option.name);
        option.set = true;
      } else if (option.type == cli_type::bool_) {
        if (!parse_cmdline_value(value, *(bool*)option.value))
          return print_cli_help(cli, "incorrect value for " + option.name);
        option.set = true;
      } else {
        throw std::runtime_error("unsupported type");
      }
    }
  }
  // check remaining
  if (!args.empty())
    return print_cli_help(cli, "mismatched value for " + args.front());
  return true;
}

}  // namespace yocto

#endif
