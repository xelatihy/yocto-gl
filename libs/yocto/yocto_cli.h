//
// # Yocto/CLI: Utilities for writing command-line apps
//
// Yocto/CLI is a collection of utilities used in writing command-line
// applications, including parsing command line arguments, printing values, and
// timers.
// Yocto/CLI is implemented in `yocto_cli.h`.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2022 Fabio Pellacini
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
#include <chrono>
#include <functional>
#include <sstream>
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
using std::pair;
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// COMMAND LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// init cli
struct cli_command;
inline cli_command make_cli(const string& name, const string& usage);

// add command
inline cli_command& add_command(
    cli_command& cli, const string& name, const string& usage);
// add command variable
template <typename T>
inline void add_command_var(cli_command& cli, T& value);
// add command by invoking add_options for the struct type passed
template <typename T>
inline cli_command& add_command(
    cli_command& cli, const string& name, T& value, const string& usage);

// add option
template <typename T>
inline void add_option(
    cli_command& cli, const string& name, T& value, const string& usage);
inline void add_option(
    cli_command& cli, const string& name, bool& value, const string& usage);
template <typename T>
inline void add_option(cli_command& cli, const string& name, vector<T>& value,
    const string& usage);
template <typename T, size_t N>
inline void add_option(cli_command& cli, const string& name, array<T, N>& value,
    const string& usage);

// add option with labels
template <typename T>
inline void add_option(cli_command& cli, const string& name, T& value,
    const string& usage, const vector<pair<T, string>>& labels);

// get usage
inline string get_usage(const cli_command& cli);

// parse cli
inline bool parse_cli(
    cli_command& cli, const vector<string>& args, string& error);
inline bool parse_cli(
    cli_command& cli, int argc, const char** argv, string& error);

// cli error
struct cli_error : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// parse cli, throws cli_error on error
inline void parse_cli(cli_command& cli, const vector<string>& args);
inline void parse_cli(cli_command& cli, int argc, const char** argv);

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

// Timer operations
inline void          start_timer(simple_timer& timer);
inline void          stop_timer(simple_timer& timer);
inline int64_t       elapsed_nanoseconds(const simple_timer& timer);
inline double        elapsed_seconds(const simple_timer& timer);
inline string        elapsed_formatted(const simple_timer& timer);
inline std::ostream& operator<<(
    std::ostream& stream, const simple_timer& timer);

}  // namespace yocto

// -----------------------------------------------------------------------------
// PRINTING VALUES
// -----------------------------------------------------------------------------
namespace yocto {

// Print a value to stdout. This is using stream for now.
// It will be moved to use std::format when it becomes available.
template <typename... Args>
inline void print(const string& format, const Args&... values);
template <typename... Args>
inline void println(const string& format, const Args&... values);

// Prints a message line.
template <typename... Args>
inline void print_info(const string& format, const Args&... values);
// Prints an error.
template <typename... Args>
inline void print_error(const string& format, const Args&... values);

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// ARGUMENT PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// Cli map. Cannot use map or unordered_map
template <typename Key, typename Value>
struct cli_map {
  // reserve to ensure pointer stability
  cli_map() { _data.reserve(256); }

  bool empty() const { return _data.empty(); }
  bool contains(const Key& key) const {
    for (auto& [key_, value] : _data)
      if (key_ == key) return true;
    return false;
  }

  Value& operator[](const Key& key) {
    for (auto& [key_, value] : _data)
      if (key_ == key) return value;
    return _data.emplace_back(key, Value{}).second;
  }
  const Value& operator[](const Key& key) const {
    for (auto& [key_, value] : _data)
      if (key_ == key) return value;
    throw std::out_of_range{"missing key"};
  }
  Value& at(const Key& key) {
    for (auto& [key_, value] : _data)
      if (key_ == key) return value;
    throw std::out_of_range{"missing key"};
  }
  const Value& at(const Key& key) const {
    for (auto& [key_, value] : _data)
      if (key_ == key) return value;
    throw std::out_of_range{"missing key"};
  }

  vector<pair<Key, Value>> _data = {};
};

// Cli data
using cli_setter = std::function<bool(const vector<string>&, string&)>;
struct cli_command {
  // options and commands
  cli_map<string, cli_setter>  options  = {};
  cli_map<string, cli_command> commands = {};
  // command
  string     command_sel = "";
  cli_setter command_var = {};
  // usage
  string usage_name     = "";
  string usage_descr    = "";
  string usage_options  = "";
  string usage_commands = "";
};

// Helpers
inline void _cli_check_option(const cli_command& cli, const string& name) {
  if (!cli.commands.empty())
    throw std::invalid_argument{"cannot add options and commands"};
  if (cli.options.contains(name))
    throw std::invalid_argument{"option already added " + name};
}
inline void _cli_check_command(const cli_command& cli, const string& name) {
  if (!cli.options.empty())
    throw std::invalid_argument{"cannot add options and commands"};
  if (cli.commands.contains(name))
    throw std::invalid_argument{"command already added " + name};
}

// parse helpers
inline bool _cli_parse_size(
    const vector<string>& args, size_t min, size_t max, string& error) {
  if (args.size() < min || args.size() > max) {
    error = "wrong number of arguments";
    return false;
  }
  return true;
}
template <typename T>
inline bool _cli_parse_value(const string& arg, T& value, string& error) {
  if constexpr (std::is_same_v<T, string>) {
    value = arg;
    return true;
  } else {
    auto stream = std::istringstream(arg);
    stream >> std::boolalpha >> value;
    if (!stream) {
      error = "parse error";
      return false;
    }
    return true;
  }
}
template <typename T>
inline bool _cli_parse_value(const string& arg, T& value, string& error,
    const vector<pair<T, string>>& labels) {
  for (auto& label : labels) {
    if (label.second == arg) {
      value = label.first;
      return true;
    }
  }
  error = "unknown value " + arg;
  return false;
}

// Usage helpers
inline string _cli_usage_command(const string& name, const string& descr) {
  auto usage = "  " + name;
  usage += string(std::max(22 - (int)usage.size(), 0), ' ');
  usage += descr + "\n";
  return usage;
}
inline string _cli_usage_option(const string& name, const string& var,
    const string& descr, const string& def, const vector<string>& labels = {}) {
  auto usage = "  --" + name + " " + var;
  usage += string(std::max(22 - (int)usage.size(), 0), ' ');
  usage += descr + "[" + def + "]\n";
  if (!labels.empty()) {
    usage += "    with labels: ";
    for (auto& label : labels) usage += label + ",";
    usage.back() = '\n';
  }
  return usage;
}
template <typename T>
inline string _cli_usage_option(
    const string& name, const T& value, const string& usage) {
  auto var = string{"var"};
  if constexpr (std::is_same_v<T, string>) var = "str";
  if constexpr (std::is_same_v<T, bool>) var = "";
  if constexpr (std::is_floating_point_v<T>) var = "num";
  if constexpr (std::is_integral_v<T> && !std::is_same_v<T, bool>) var = "int";
  auto stream = std::ostringstream();
  stream << std::boolalpha << value;
  auto def = stream.str();
  return _cli_usage_option(name, var, usage, def);
}
template <typename T>
inline string _cli_usage_option(
    const string& name, const vector<T>& value, const string& usage) {
  auto var = string{"var"};
  if constexpr (std::is_same_v<T, string>) var = "str";
  if constexpr (std::is_same_v<T, bool>) var = "bool";
  if constexpr (std::is_floating_point_v<T>) var = "num";
  if constexpr (std::is_integral_v<T> && !std::is_same_v<T, bool>) var = "int";
  var += "[]";
  auto stream = std::ostringstream();
  for (auto& item : value) {
    if (&item != &value.front()) stream << ",";
    stream << std::boolalpha << item;
  }
  auto def = stream.str();
  return _cli_usage_option(name, var, usage, def);
}
template <typename T, size_t N>
inline string _cli_usage_option(
    const string& name, const array<T, N>& value, const string& usage) {
  auto var = string{"var"};
  if constexpr (std::is_same_v<T, string>) var = "str";
  if constexpr (std::is_same_v<T, bool>) var = "bool";
  if constexpr (std::is_floating_point_v<T>) var = "num";
  if constexpr (std::is_integral_v<T> && !std::is_same_v<T, bool>) var = "int";
  var += "[" + std::to_string(N) + "]";
  auto stream = std::ostringstream();
  for (auto& item : value) {
    if (&item != &value.front()) stream << ",";
    stream << std::boolalpha << item;
  }
  auto def = stream.str();
  return _cli_usage_option(name, var, usage, def);
}
template <typename T>
inline string _cli_usage_option(const string& name, const T& value,
    const string& usage, const vector<pair<T, string>>& labels) {
  auto var     = string{"str"};
  auto def     = string{};
  auto labels_ = vector<string>{};
  for (auto& label : labels) {
    labels_.push_back(label.second);
    if (label.first == value) def = label.second;
  }
  if (def.empty()) throw std::invalid_argument{"undefined key"};
  return _cli_usage_option(name, var, usage, def, labels_);
}

// init cli
inline cli_command make_cli(const string& name, const string& usage) {
  auto cli        = cli_command{};
  cli.usage_name  = name;
  cli.usage_descr = usage;
  cli.command_var = [](const vector<string>&, string&) { return true; };
  return cli;
}

// add command
inline cli_command& add_command(
    cli_command& cli, const string& name, const string& usage) {
  _cli_check_command(cli, name);
  auto& cmd       = cli.commands[name];
  cmd.usage_name  = cli.usage_name + " " + name;
  cmd.usage_descr = usage;
  cli.usage_commands += _cli_usage_command(name, usage);
  cmd.command_var = [](const vector<string>&, string&) { return true; };
  return cmd;
}

// add command variable
template <typename T>
inline void add_command_var(cli_command& cli, T& value) {
  cli.command_var = [&value](const vector<string>& args, string& error) {
    if (!_cli_parse_size(args, 1, 1, error)) return false;
    if (!_cli_parse_value(args.front(), value, error)) return false;
    return true;
  };
}

// add command by invoking add_options for the struct type passed
template <typename T>
inline cli_command& add_command(
    cli_command& cli, const string& name, T& value, const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_options(cmd, value);
  return cmd;
}

// add option
template <typename T>
inline void add_option(
    cli_command& cli, const string& name, T& value, const string& usage) {
  _cli_check_option(cli, name);
  cli.usage_options += _cli_usage_option(name, value, usage);
  cli.options[name] = [&value](const vector<string>& args, string& error) {
    if (!_cli_parse_size(args, 1, 1, error)) return false;
    if (!_cli_parse_value(args.front(), value, error)) return false;
    return true;
  };
}

// add option
inline void add_option(
    cli_command& cli, const string& name, bool& value, const string& usage) {
  _cli_check_option(cli, name);
  cli.usage_options += _cli_usage_option(name, value, usage);
  cli.options[name] = [&value](const vector<string>& args, string& error) {
    if (!_cli_parse_size(args, 0, 1, error)) return false;
    if (!_cli_parse_value(args.empty() ? "true" : args.front(), value, error))
      return false;
    return true;
  };
}

// add option
template <typename T>
inline void add_option(cli_command& cli, const string& name, vector<T>& value,
    const string& usage) {
  _cli_check_option(cli, name);
  cli.usage_options += _cli_usage_option(name, value, usage);
  cli.options[name] = [&value](const vector<string>& args, string& error) {
    if (!_cli_parse_size(args, 1, 1024, error)) return false;
    value.resize(args.size());
    for (auto idx = 0; idx < (int)args.size(); idx++) {
      if (!_cli_parse_value(args[idx], value[idx], error)) return false;
    }
    return true;
  };
}

// add option
template <typename T, size_t N>
inline void add_option(cli_command& cli, const string& name, array<T, N>& value,
    const string& usage) {
  _cli_check_option(cli, name);
  cli.usage_options += _cli_usage_option(name, value, usage);
  cli.options[name] = [&value](const vector<string>& args, string& error) {
    if (!_cli_parse_size(args, N, N, error)) return false;
    for (auto idx = 0; idx < (int)args.size(); idx++) {
      if (!_cli_parse_value(args[idx], value[idx], error)) return false;
    }
    return true;
  };
}

// add option
template <typename T>
inline void add_option(cli_command& cli, const string& name, T& value,
    const string& usage, const vector<pair<T, string>>& labels) {
  _cli_check_option(cli, name);
  cli.usage_options += _cli_usage_option(name, value, usage, labels);
  cli.options[name] = [&value, labels](
                          const vector<string>& args, string& error) {
    if (!_cli_parse_size(args, 1, 1, error)) return false;
    if (!_cli_parse_value(args.front(), value, error, labels)) return false;
    return true;
  };
}

// parse cli
inline bool parse_cli(
    cli_command& cli, const vector<string>& args, size_t pos, string& error) {
  // parsee command or options
  if (!cli.commands.empty()) {
    // check command
    if (pos >= args.size()) {
      error = "missing command";
      return false;
    }
    // check help
    if (args[pos] == "--help") {
      error = "help invoked";
      return false;
    }
    // check option
    if (args[pos].find("--") == 0) {
      error = "missing command";
      return false;
    }
    // verify command
    if (!cli.commands.contains(args[pos])) {
      error = "unknown command " + args[pos];
      return false;
    }
    // get command
    auto name       = args[pos++];
    cli.command_sel = name;
    // set command
    if (!cli.command_var({name}, error)) {
      error += " for command " + name;
      return false;
    }
    // parse command recursively
    auto& command = cli.commands.at(name);
    return parse_cli(command, args, pos, error);
  } else {
    // check option
    if (pos < args.size() && args[pos].find("--") != 0) {
      error = "command should start with option";
      return false;
    }
    // parse options till done
    while (pos < args.size()) {
      // check for bugs
      if (args[pos].find("--") != 0) throw std::runtime_error{"parsing bug"};
      // check help
      if (args[pos] == "--help") {
        error = "help invoked";
        return false;
      }
      // verify command
      if (!cli.options.contains(args[pos].substr(2))) {
        error = "unknown option " + args[pos].substr(2);
        return false;
      }
      // get option
      auto  name   = args[pos++].substr(2);
      auto& option = cli.options.at(name);
      // get args
      auto values = vector<string>{};
      while (pos < args.size() && args[pos].find("--") != 0) {
        values.push_back(args[pos++]);
      }
      // set option
      if (!option(values, error)) {
        error += " for option " + name;
        return false;
      }
    }
    // done
    return true;
  }
}

inline string get_usage(const cli_command& cli) {
  if (!cli.command_sel.empty())
    return get_usage(cli.commands.at(cli.command_sel));
  auto usage = cli.usage_name + " [options]" +
               (cli.commands.empty() ? "" : " <command>") + "\n" +
               cli.usage_descr + "\n";
  if (!cli.usage_commands.empty())
    usage += "\ncommands:\n" + cli.usage_commands;
  if (!cli.usage_options.empty()) usage += "\noptions:\n" + cli.usage_options;
  usage += "\n";
  return usage;
}

// parse cli
inline bool parse_cli(
    cli_command& cli, const vector<string>& args, string& error) {
  return parse_cli(cli, args, 1, error);
}
inline bool parse_cli(
    cli_command& cli, int argc, const char** argv, string& error) {
  return parse_cli(cli, vector<string>{argv, argv + argc}, 1, error);
}

// parse cli
inline void parse_cli(cli_command& cli, const vector<string>& args) {
  auto error = string{};
  if (!parse_cli(cli, args, 1, error))
    throw cli_error{error + "\n\nusage: " + get_usage(cli)};
}
inline void parse_cli(cli_command& cli, int argc, const char** argv) {
  return parse_cli(cli, vector<string>{argv, argv + argc});
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIMPLE TIMER
// -----------------------------------------------------------------------------
namespace yocto {

// get time in nanoseconds - useful only to compute difference of times
inline int64_t _get_time() {
  return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}

// Format duration string from nanoseconds
inline string _format_duration(int64_t duration) {
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

// Simple timer
inline simple_timer::simple_timer() {
  start = _get_time();
  stop  = -1;
}

// Timer opreations
inline void start_timer(simple_timer& timer) {
  timer.start = _get_time();
  timer.stop  = -1;
}
inline void    stop_timer(simple_timer& timer) { timer.stop = _get_time(); }
inline int64_t elapsed_nanoseconds(const simple_timer& timer) {
  return (timer.stop < 0 ? _get_time() : timer.stop) - timer.start;
}
inline double elapsed_seconds(const simple_timer& timer) {
  return (double)elapsed_nanoseconds(timer) * 1e-9;
}
inline string elapsed_formatted(const simple_timer& timer) {
  return _format_duration(elapsed_nanoseconds(timer));
}
inline std::ostream& operator<<(
    std::ostream& stream, const simple_timer& timer) {
  return stream << elapsed_formatted(timer);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PRINTING VALUES
// -----------------------------------------------------------------------------
namespace yocto {

// Format to value
inline void format_to(std::stringstream& stream, const string& format) {
  auto pos = format.find("{}");
  if (pos != string::npos) throw std::invalid_argument("bad format string");
  stream << format;
}
template <typename Arg, typename... Args>
inline void format_to(std::stringstream& stream, const string& format,
    const Arg& arg, const Args&... args) {
  auto pos = format.find("{}");
  if (pos == string::npos) throw std::invalid_argument("bad format string");
  stream << format.substr(0, pos) << arg;
  format_to(stream, format.substr(pos + 2), args...);
}

// Print a value to stdout. This is using stream for now.
// It will be moved to use std::format when it becomes available.
template <typename... Args>
inline void print(const string& format, const Args&... values) {
  auto stream = std::stringstream{};
  format_to(stream, format, values...);
  printf("%s", stream.str().c_str());
}
template <typename... Args>
inline void println(const string& format, const Args&... values) {
  auto stream = std::stringstream{};
  format_to(stream, format, values...);
  printf("%s\n", stream.str().c_str());
}

// Prints a message line.
template <typename... Args>
inline void print_info(const string& format, const Args&... values) {
  auto stream = std::stringstream{};
  format_to(stream, format, values...);
  printf("%s\n", stream.str().c_str());
  fflush(stdout);
}
// Prints an error.
template <typename... Args>
inline void print_error(const string& format, const Args&... values) {
  auto stream = std::stringstream{};
  format_to(stream, format, values...);
  printf("error: %s\n", stream.str().c_str());
  fflush(stdout);
}

}  // namespace yocto

#endif
