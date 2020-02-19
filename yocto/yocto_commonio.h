//
// # Yocto/CommonIO: Tiny collection of IO utilities
//
//
// Yocto/CommonIO is a collection of utilities used in writing Yocto/GL
// libraries and example applications. We support printing and parsing builtin
// types, parsing command line arguments, simple path manipulation, file
// lading/saving.
//
//
// ## Printing values
//
// Use `print_info()` to print a message, `print_fatal()` to print and exit.
// To time a block of code use `print_timed()` to use an RIIA timer or
// call `print_elapsed()` to print the elapsed time as needed.
// Several overloads of `to_string()` are provided for both the basic types
// and Yocto/Math types.
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
// 2. add options with `add_option(cli, name, value, usage, req)`
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
// 3. use `file` as a safe wrapper over C streams; use `open_file()`,
//  `close_file()`, `read_line()`, `read_value()`, `write_text()` and
//  `write_value()` to operate on the file.
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

#ifndef _YOCTO_COMMONIO_H_
#define _YOCTO_COMMONIO_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto::commonio {

using byte = unsigned char;
using std::string;
using std::unordered_map;
using std::vector;

}  // namespace yocto::commonio

// -----------------------------------------------------------------------------
// PRINT/FORMATTING UTILITIES
// -----------------------------------------------------------------------------
namespace yocto::commonio {

// Print a message to the console
inline void print_info(const string& msg);
// Prints a messgae to the console and exit with an error.
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

}  // namespace yocto::commonio

// -----------------------------------------------------------------------------
// COMMAND LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto::commonio {

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

// Parse an int, float, string, and bool option or positional argument.
// Options's names starts with "--" or "-", otherwise they are arguments.
// The library support using many names for the same option/argument
// separate by commas. Boolean flags are indicated with a pair of names
// "--name/--no-name", so that we have both options available.
inline void add_option(cli_state& cli, const string& name, string& value,
    const string& usage, bool req = false);
inline void add_option(cli_state& cli, const string& name, int& value,
    const string& usage, bool req = false);
inline void add_option(cli_state& cli, const string& name, float& value,
    const string& usage, bool req = false);
inline void add_option(cli_state& cli, const string& name, bool& value,
    const string& usage, bool req = false);
// Parse an enum
inline void add_option(cli_state& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, bool req = false);
template <typename T, typename = std::enable_if_t<std::is_enum_v<T>>>
inline void add_option(cli_state& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, bool req = false);
// Parse all arguments left on the command line.
inline void add_option(cli_state& cli, const string& name,
    vector<string>& value, const string& usage, bool req = false);

}  // namespace yocto::commonio

// -----------------------------------------------------------------------------
// PATH UTILITIES
// -----------------------------------------------------------------------------
namespace yocto::commonio {

// These utilities are here only for backward compatibility. They should be
// considered deprecated.

// Utility to normalize a path
inline string normalize_path(const string& filename);

// Get directory name (including '/').
inline string get_dirname(const string& filename);

// Get extension (not including '.').
inline string get_extension(const string& filename);

// Get filename without directory.
inline string get_filename(const string& filename);

// Get extension.
inline string get_noextension(const string& filename);

// Get filename without directory and extension.
inline string get_basename(const string& filename);

// Replaces extensions
inline string replace_extension(const string& filename, const string& ext);

// Check if a file can be opened for reading.
inline bool exists_file(const string& filename);
}  // namespace yocto::commonio

// -----------------------------------------------------------------------------
// FILE IO
// -----------------------------------------------------------------------------
namespace yocto::commonio {

// Load/save a text file
inline bool load_text(const string& filename, string& str, string& error);
inline bool save_text(const string& filename, const string& str, string& error);

// Load/save a binary file
inline bool load_binary(
    const string& filename, vector<byte>& data, string& error);
inline bool save_binary(
    const string& filename, const vector<byte>& data, string& error);

}  // namespace yocto::commonio

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
namespace yocto::commonio {

// This is a very crude replacement for `std::format()` that will be used when
// available on all platforms.
template <typename... Args>
inline string format(const string& fmt, Args&&... args);

}  // namespace yocto::commonio

// -----------------------------------------------------------------------------
// PRINT/FORMATTING UTILITIES
// -----------------------------------------------------------------------------
namespace yocto::commonio {

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
  using clock               = std::chrono::high_resolution_clock;
  static int64_t start_time = 0;
  if (current == 0) start_time = clock::now().time_since_epoch().count();
  auto elapsed = clock::now().time_since_epoch().count() - start_time;
  elapsed /= 1000000;  // millisecs
  auto mins  = pad(std::to_string(elapsed / 60000), 2);
  auto secs  = pad(std::to_string((elapsed % 60000) / 1000), 2);
  auto msecs = pad(std::to_string((elapsed % 60000) % 1000), 3);
  auto n     = (int)(30 * (float)current / (float)total);
  auto bar   = "[" + pade(string(n, '='), 30) + "]";
  auto line  = bar + " " + mins + ":" + secs + "." + msecs + " " +
              pade(message, 30);
  printf("\r%s\r", line.c_str());
  if (current == total) printf("\n");
  fflush(stdout);
}

}  // namespace yocto::commonio

// -----------------------------------------------------------------------------
// PATH UTILITIES
// -----------------------------------------------------------------------------
namespace yocto::commonio {

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

}  // namespace yocto::commonio

// -----------------------------------------------------------------------------
// FILE IO
// -----------------------------------------------------------------------------
namespace yocto::commonio {

// Load a text file
inline bool load_text(const string& filename, string& str, string& error) {
  // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
  auto fs = fopen(filename.c_str(), "rt");
  if (!fs) {
    error = filename + ": file not found";
    return false;
  }
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};
  fseek(fs, 0, SEEK_END);
  auto length = ftell(fs);
  fseek(fs, 0, SEEK_SET);
  str.resize(length);
  if (fread(str.data(), 1, length, fs) != length) {
    error = filename + ": read error";
    return false;
  }
  return true;
}

// Save a text file
inline bool save_text(
    const string& filename, const string& str, string& error) {
  auto fs = fopen(filename.c_str(), "wt");
  if (!fs) {
    error = filename + ": file not found";
    return false;
  }
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};
  if (fprintf(fs, "%s", str.c_str()) < 0) {
    error = filename + ": write error";
    return false;
  }
  fclose(fs);
  return true;
}

// Load a binary file
inline bool load_binary(
    const string& filename, vector<byte>& data, string& error) {
  // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
  auto fs = fopen(filename.c_str(), "rb");
  if (!fs) {
    error = filename + ": file not found";
    return false;
  }
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};
  fseek(fs, 0, SEEK_END);
  auto length = ftell(fs);
  fseek(fs, 0, SEEK_SET);
  data.resize(length);
  if (fread(data.data(), 1, length, fs) != length) {
    error = filename + ": read error";
    return false;
  }
  fclose(fs);
  return true;
}

// Save a binary file
inline bool save_binary(
    const string& filename, const vector<byte>& data, string& error) {
  auto fs = fopen(filename.c_str(), "wb");
  if (!fs) {
    error = filename + ": file not found";
    return false;
  }
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};
  if (fwrite(data.data(), 1, data.size(), fs) != data.size()) {
    error = filename + ": rewritead error";
    return false;
  }
  fclose(fs);
  return true;
}

}  // namespace yocto::commonio

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF COMMAND-LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto::commonio {

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
  bool                   help            = false;
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

inline void add_option(cli_state& cli, const string& name, cli_type type,
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

inline void add_option(cli_state& cli, const string& name, string& value,
    const string& usage, bool req) {
  return add_option(cli, name, cli_type::string_, &value, usage, req, {});
}
inline void add_option(cli_state& cli, const string& name, int& value,
    const string& usage, bool req) {
  return add_option(cli, name, cli_type::int_, &value, usage, req, {});
}
inline void add_option(cli_state& cli, const string& name, float& value,
    const string& usage, bool req) {
  return add_option(cli, name, cli_type::float_, &value, usage, req, {});
}
inline void add_option(cli_state& cli, const string& name, bool& value,
    const string& usage, bool req) {
  return add_option(cli, name, cli_type::flag_, &value, usage, req, {});
}
inline void add_option(cli_state& cli, const string& name,
    vector<string>& value, const string& usage, bool req) {
  return add_option(
      cli, name, cli_type::string_vector_, &value, usage, req, {});
}
inline void add_flag(cli_state& cli, const string& name, bool& value,
    const string& usage, bool req) {
  return add_option(cli, name, cli_type::flag_, &value, usage, req, {});
}
inline void add_option(cli_state& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, bool req) {
  return add_option(cli, name, cli_type::enum_, &value, usage, req, choices);
}
template <typename T, typename>
inline void add_option(cli_state& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, bool req) {
  return add_option(cli, name, (int&)value, usage, choices, req);
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

inline bool parse_cli_value(const string& str, int& value) {
  auto end = (char*)nullptr;
  value    = (int)strtol(str.c_str(), &end, 10);
  return end != nullptr;
}
inline bool parse_cli_value(const string& str, float& value) {
  auto end = (char*)nullptr;
  value    = strtof(str.c_str(), &end);
  return end != nullptr;
}
inline bool parse_cli_value(const string& str, bool& value) {
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

inline bool parse_cli(
    cli_state& cli, int argc, const char** argv, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

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
        throw std::runtime_error("inconsistent option type for " + name);
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
          return cli_error("missing value for " + name);
        auto value = args[pos + 1];
        args.erase(args.begin() + pos, args.begin() + pos + 2);
        if (option.type == cli_type::string_) {
          *(string*)option.value = value;
          option.set             = true;
        } else if (option.type == cli_type::int_) {
          if (!parse_cli_value(value, *(int*)option.value))
            return cli_error("incorrect value for " + name);
          option.set = true;
        } else if (option.type == cli_type::float_) {
          if (!parse_cli_value(value, *(float*)option.value))
            return cli_error("incorrect value for " + name);
          option.set = true;
        } else if (option.type == cli_type::bool_) {
          if (!parse_cli_value(value, *(bool*)option.value))
            return cli_error("incorrect value for " + name);
          option.set = true;
        } else if (option.type == cli_type::enum_) {
          auto pos = std::find(
              option.choices.begin(), option.choices.end(), value);
          if (pos == option.choices.end())
            return cli_error("incorrect value for " + name);
          else
            *(int*)option.value = (int)(pos - option.choices.begin());
          option.set = true;
        } else {
          throw std::runtime_error("unsupported type");
        }
      }
    }
    if (option.req && !option.set)
      return cli_error("missing value for " + option.name);
  }
  // check unknown options
  for (auto& arg : args) {
    if (arg.find("-") == 0) return cli_error("unknown option " + arg);
  }
  // parse positional
  for (auto& option : cli.options) {
    if (option.name[0] == '-') continue;
    if (args.empty()) {
      if (option.req) return cli_error("missing value for " + option.name);
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
        if (!parse_cli_value(value, *(int*)option.value))
          return cli_error("incorrect value for " + option.name);
        option.set = true;
      } else if (option.type == cli_type::float_) {
        if (!parse_cli_value(value, *(float*)option.value))
          return cli_error("incorrect value for " + option.name);
        option.set = true;
      } else if (option.type == cli_type::bool_) {
        if (!parse_cli_value(value, *(bool*)option.value))
          return cli_error("incorrect value for " + option.name);
        option.set = true;
      } else {
        throw std::runtime_error("unsupported type");
      }
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

}  // namespace yocto::commonio

#endif
