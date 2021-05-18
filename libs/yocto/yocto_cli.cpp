//
// Implementation for Yocto/Bvh
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

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_cli.h"

#include <array>
#include <chrono>
#include <cstdio>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "ext/CLI11.hpp"
#include "ext/json.hpp"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::pair;
using std::string;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// PRINT/FORMATTING UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Print a message to the console
void print_info(const string& message) { printf("%s\n", message.c_str()); }
// Prints a message to the console and exit with an error.
int print_fatal(const string& message) {
  printf("\n%s\n", message.c_str());
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
print_timer print_timed(const string& message) {
  printf("%s", message.c_str());
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

static int    print_progress_current = 0;
static int    print_progress_total   = 0;
static string print_progress_message = "";
void          print_progress_begin(const string& message, int total) {
  print_progress_current = 0;
  print_progress_total   = total;
  print_progress_message = message;
  print_progress(
      print_progress_message, print_progress_current, print_progress_total);
}
void print_progress_end() {
  print_progress_current = print_progress_total;
  print_progress(
      print_progress_message, print_progress_current, print_progress_total);
}
void print_progress_next() {
  print_progress_current += 1;
  print_progress(
      print_progress_message, print_progress_current, print_progress_total);
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
// IMPLEMENTATION OF COMMAND-LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto {

template <typename T>
static vector<pair<string, T>> make_cli_map(const vector<string>& choices) {
  auto ret = vector<pair<string, T>>{};
  for (auto idx = 0; idx < (int)choices.size(); idx++) {
    ret.push_back({choices[idx], (T)idx});
  }
  return ret;
}
template <>
vector<pair<string, string>> make_cli_map<string>(
    const vector<string>& choices) {
  auto ret = vector<pair<string, string>>{};
  for (auto idx = 0; idx < (int)choices.size(); idx++) {
    ret.push_back({choices[idx], choices[idx]});
  }
  return ret;
}

template <typename T>
static void add_option_impl(const cli_command& cli, const string& name,
    T& value, const string& usage, const vector<T>& minmax,
    const vector<string>& choices, const string& alt, bool req) {
  auto cli11 = (CLI::App*)cli.state;
  if constexpr (std::is_same_v<T, bool>) {
    cli11->add_flag(
        "--" + name + (alt.empty() ? "" : (",-" + alt)), value, usage);
  } else {
    auto option = cli11->add_option(
        "--" + name + (alt.empty() ? "" : (",-" + alt)), value, usage);
    if (minmax.size() == 2) option->check(CLI::Bound(minmax[0], minmax[1]));
    if (!choices.empty())
      option->transform(CLI::CheckedTransformer(make_cli_map<T>(choices)));
  }
}

template <typename T, size_t N>
static void add_option_impl(const cli_command& cli, const string& name,
    array<T, N>& value, const string& usage, const vector<T>& minmax,
    const vector<string>& choices, const string& alt, bool req) {
  auto cli11  = (CLI::App*)cli.state;
  auto option = cli11->add_option(
      "--" + name + (alt.empty() ? "" : (",-" + alt)), value, usage);
  if constexpr (!std::is_same_v<T, bool>) {
    if (minmax.size() == 2) option->check(CLI::Bound(minmax[0], minmax[1]));
    if (!choices.empty())
      option->transform(CLI::CheckedTransformer(make_cli_map<T>(choices)));
  }
}
template <typename T>
static void add_argument_impl(const cli_command& cli, const string& name,
    T& value, const string& usage, const vector<T>& minmax,
    const vector<string>& choices, bool req) {
  auto cli11  = (CLI::App*)cli.state;
  auto option = cli11->add_option(name, value, usage);
  if constexpr (!std::is_same_v<T, bool>) {
    if (minmax.size() == 2) option->check(CLI::Bound(minmax[0], minmax[1]));
    if (!choices.empty())
      option->transform(CLI::CheckedTransformer(make_cli_map<T>(choices)));
  }
}

template <typename T, size_t N>
static void add_argument_impl(const cli_command& cli, const string& name,
    array<T, N>& value, const string& usage, const vector<T>& minmax,
    const vector<string>& choices, bool req) {
  auto cli11  = (CLI::App*)cli.state;
  auto option = cli11->add_option(name, value, usage);
  if constexpr (!std::is_same_v<T, bool>) {
    if (minmax.size() == 2) option->check(CLI::Bound(minmax[0], minmax[1]));
    if (!choices.empty())
      option->transform(CLI::CheckedTransformer(make_cli_map<T>(choices)));
  }
}

template <typename T>
static void add_argumentv_impl(const cli_command& cli, const string& name,
    vector<T>& value, const string& usage, const vector<T>& minmax,
    const vector<string>& choices, bool req) {
  auto cli11  = (CLI::App*)cli.state;
  auto option = cli11->add_option(name, value, usage);
  if constexpr (!std::is_same_v<T, bool>) {
    if (minmax.size() == 2) option->check(CLI::Bound(minmax[0], minmax[1]));
    if (!choices.empty())
      option->transform(CLI::CheckedTransformer(make_cli_map<T>(choices)));
  }
}

// Add an optional argument. Supports strings, numbers, and boolean flags.
// Optional arguments will be parsed with name `--<name>` and `-<alt>`.
// Optional booleans will support both `--<name>` and `--no-<name>` to enabled
// and disable the flag.
void add_option(const cli_command& cli, const string& name, int& value,
    const string& usage, const vector<int>& minmax, const string& alt,
    bool req) {
  return add_option_impl(cli, name, value, usage, minmax, {}, alt, req);
}
void add_option(const cli_command& cli, const string& name, float& value,
    const string& usage, const vector<float>& minmax, const string& alt,
    bool req) {
  return add_option_impl(cli, name, value, usage, minmax, {}, alt, req);
}
void add_option(const cli_command& cli, const string& name, bool& value,
    const string& usage, const vector<string>& choices, const string& alt,
    bool req) {
  return add_option_impl(cli, name, value, usage, {}, choices, alt, req);
}
void add_option(const cli_command& cli, const string& name, string& value,
    const string& usage, const vector<string>& choices, const string& alt,
    bool req) {
  return add_option_impl(cli, name, value, usage, {}, choices, alt, req);
}
void add_option(const cli_command& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, const string& alt,
    bool req) {
  return add_option_impl(cli, name, value, usage, {}, choices, alt, req);
}
// Add a positional argument. Supports strings, numbers, and boolean flags.
void add_argument(const cli_command& cli, const string& name, int& value,
    const string& usage, const vector<int>& minmax, bool req) {
  return add_argument_impl(cli, name, value, usage, minmax, {}, req);
}
void add_argument(const cli_command& cli, const string& name, float& value,
    const string& usage, const vector<float>& minmax, bool req) {
  return add_argument_impl(cli, name, value, usage, minmax, {}, req);
}
void add_argument(const cli_command& cli, const string& name, bool& value,
    const string& usage, const vector<string>& choices, bool req) {
  return add_argument_impl(cli, name, value, usage, {}, choices, req);
}
void add_argument(const cli_command& cli, const string& name, string& value,
    const string& usage, const vector<string>& choices, bool req) {
  return add_argument_impl(cli, name, value, usage, {}, choices, req);
}
void add_argument(const cli_command& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, bool req) {
  return add_argument_impl(cli, name, value, usage, {}, choices, req);
}
// Add a positional argument that consumes all arguments left.
// Supports strings and enums.
void add_argument(const cli_command& cli, const string& name,
    vector<int>& value, const string& usage, const vector<int>& minmax,
    bool req) {
  return add_argumentv_impl(cli, name, value, usage, minmax, {}, req);
}
void add_argument(const cli_command& cli, const string& name,
    vector<float>& value, const string& usage, const vector<float>& minmax,
    bool req) {
  return add_argumentv_impl(cli, name, value, usage, minmax, {}, req);
}
void add_argument(const cli_command& cli, const string& name,
    vector<int>& value, const string& usage, const vector<string>& choices,
    bool req) {
  return add_argumentv_impl(cli, name, value, usage, {}, choices, req);
}
void add_argument(const cli_command& cli, const string& name,
    vector<string>& value, const string& usage, const vector<string>& choices,
    bool req) {
  return add_argumentv_impl(cli, name, value, usage, {}, choices, req);
}

// Add an optional argument. Supports basic math types.
void add_option(const cli_command& cli, const string& name, vec2i& value,
    const string& usage, const vector<int>& minmax, const string& alt,
    bool req) {
  return add_option_impl(
      cli, name, (array<int, 2>&)value, usage, minmax, {}, alt, req);
}
void add_option(const cli_command& cli, const string& name, vec3i& value,
    const string& usage, const vector<int>& minmax, const string& alt,
    bool req) {
  return add_option_impl(
      cli, name, (array<int, 3>&)value, usage, minmax, {}, alt, req);
}
void add_option(const cli_command& cli, const string& name, vec4i& value,
    const string& usage, const vector<int>& minmax, const string& alt,
    bool req) {
  return add_option_impl(
      cli, name, (array<int, 4>&)value, usage, minmax, {}, alt, req);
}
void add_option(const cli_command& cli, const string& name, vec2f& value,
    const string& usage, const vector<float>& minmax, const string& alt,
    bool req) {
  return add_option_impl(
      cli, name, (array<float, 2>&)value, usage, minmax, {}, alt, req);
}
void add_option(const cli_command& cli, const string& name, vec3f& value,
    const string& usage, const vector<float>& minmax, const string& alt,
    bool req) {
  return add_option_impl(
      cli, name, (array<float, 3>&)value, usage, minmax, {}, alt, req);
}
void add_option(const cli_command& cli, const string& name, vec4f& value,
    const string& usage, const vector<float>& minmax, const string& alt,
    bool req) {
  return add_option_impl(
      cli, name, (array<float, 4>&)value, usage, minmax, {}, alt, req);
}

// initialize a command line parser
cli_state make_cli(const string& name, const string& usage) {
  auto cli  = cli_state{};
  cli.state = {
      new CLI::App(usage, name), [](void* state) { delete (CLI::App*)state; }};
  return cli;
}

// add command
cli_command add_command(
    const cli_command& cli, const string& name, const string& usage) {
  auto cli11 = (CLI::App*)cli.state;
  cli11->require_subcommand(1);
  return {cli11->add_subcommand(name, usage)};
}

void add_command_name(const cli_command& cli, const string& name, string& value,
    const string& usage) {
  auto cli11 = (CLI::App*)cli.state;
  cli11->final_callback([&value, cli11] {
    value = cli11->get_subcommands().front()->get_name();
  });
}

string get_command(const cli_state& cli) {
  auto cli11 = (CLI::App*)cli.state.get();
  return cli11->get_subcommands().front()->get_name();
}

void parse_cli(cli_state& cli, int argc, const char** argv) {
  auto cli11 = (CLI::App*)cli.state.get();
  try {
    cli11->parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    exit(cli11->exit(e));
  }
}

}  // namespace yocto
