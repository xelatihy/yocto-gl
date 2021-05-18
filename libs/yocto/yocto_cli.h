//
// # Yocto/CLI: Utilities for writing command-line apps
//
// Yocto/CLI is a collection of utilities used in writing command-line
// applications, including parsing command line arguments, printing values,
// timers and progress bars.
// Yocto/CLI is implemented in `yocto_cli.h`.
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

#ifndef _YOCTO_CLI_H_
#define _YOCTO_CLI_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <chrono>
#include <cstdio>
#include <functional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::function;
using std::string;
using std::vector;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// PRINT/LOG/FORMATTING UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Print a message to the console
void print_info(const string& message);
// Prints a message to the console and exit with an error. Returns error code.
int print_fatal(const string& message);

// Timer that prints as scope end. Create with `print_timed` and print with
// `print_elapsed`.
struct print_timer {
  int64_t start_time = -1;
  ~print_timer();  // print time if scope ends
};
// Print traces for timing and program debugging
print_timer print_timed(const string& message);
int64_t     print_elapsed(print_timer& timer);

// Print progress
void print_progress(const string& message, int current, int total);
void print_progress_begin(const string& message, int total = 1);
void print_progress_end();
void print_progress_next();

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
void add_option(cli_command& cli, const string& name, int& value,
    const string& usage, const vector<int>& minmax = {}, const string& alt = "",
    bool req = false);
void add_option(cli_command& cli, const string& name, float& value,
    const string& usage, const vector<float>& minmax = {},
    const string& alt = "", bool req = false);
void add_option(cli_command& cli, const string& name, bool& value,
    const string& usage, const vector<string>& choices = {},
    const string& alt = "", bool req = false);
void add_option(cli_command& cli, const string& name, string& value,
    const string& usage, const vector<string>& choices = {},
    const string& alt = "", bool req = false);
// Add a positional argument. Supports strings, numbers, and boolean flags.
void add_argument(cli_command& cli, const string& name, int& value,
    const string& usage, const vector<int>& minmax = {}, bool req = true);
void add_argument(cli_command& cli, const string& name, float& value,
    const string& usage, const vector<float>& minmax = {}, bool req = true);
void add_argument(cli_command& cli, const string& name, bool& value,
    const string& usage, const vector<string>& choices = {}, bool req = true);
void add_argument(cli_command& cli, const string& name, string& value,
    const string& usage, const vector<string>& choices = {}, bool req = true);
// Add an optional argument with values as labels. Supports integers, enums and
// strings.
void add_option(cli_command& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, const string& alt = "",
    bool req = false);
template <typename T, typename = std::enable_if_t<std::is_enum_v<T>>>
inline void add_option(cli_command& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, const string& alt = "",
    bool req = false);
// Add a positional argument with values as labels. Supports string, integers
// and enums.
void add_argument(cli_command& cli, const string& name, int& value,
    const string& usage, const vector<string>& choices, bool req = true);
template <typename T, typename = std::enable_if_t<std::is_enum_v<T>>>
inline void add_argument(cli_command& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, bool req = true);
// Add a positional argument that consumes all arguments left.
// Supports strings and enums.
void add_argument(cli_command& cli, const string& name, vector<int>& value,
    const string& usage, const vector<int>& minmax, bool req = true);
void add_argument(cli_command& cli, const string& name, vector<float>& value,
    const string& usage, const vector<float>& minmax, bool req = true);
void add_argument(cli_command& cli, const string& name, vector<int>& value,
    const string& usage, const vector<string>& choices = {}, bool req = true);
void add_argument(cli_command& cli, const string& name, vector<string>& value,
    const string& usage, const vector<string>& choices = {}, bool req = true);

// Add an optional argument. Supports basic math types.
void add_option(cli_command& cli, const string& name, vec2i& value,
    const string& usage, const vector<int>& minmax = {}, const string& alt = "",
    bool req = false);
void add_option(cli_command& cli, const string& name, vec3i& value,
    const string& usage, const vector<int>& minmax = {}, const string& alt = "",
    bool req = false);
void add_option(cli_command& cli, const string& name, vec4i& value,
    const string& usage, const vector<int>& minmax = {}, const string& alt = "",
    bool req = false);
void add_option(cli_command& cli, const string& name, vec2f& value,
    const string& usage, const vector<float>& minmax = {},
    const string& alt = "", bool req = false);
void add_option(cli_command& cli, const string& name, vec3f& value,
    const string& usage, const vector<float>& minmax = {},
    const string& alt = "", bool req = false);
void add_option(cli_command& cli, const string& name, vec4f& value,
    const string& usage, const vector<float>& minmax = {},
    const string& alt = "", bool req = false);

// Add a subcommand
cli_command& add_command(
    cli_command& cli, const string& name, const string& usage);
void add_command_name(
    cli_command& cli, const string& name, string& value, const string& usage);

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

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
  vector<cli_value>                 minmax     = {};
  vector<string>                    choices    = {};
  vector<cli_value>                 value      = {};
  vector<cli_value>                 def        = {};
  bool                              set        = false;
  function<bool(const cli_option&)> set_value  = {};
};
// Command line command. All data should be considered private.
struct cli_command {
  string                        name        = "";
  string                        usage       = "";
  vector<cli_command>           commands    = {};
  vector<cli_option>            options     = {};
  vector<cli_option>            arguments   = {};
  bool                          help        = false;
  string                        command     = "";
  function<void(const string&)> set_command = {};
};

template <typename T, typename>
inline void add_option(cli_command& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, const string& alt,
    bool req) {
  return add_option(
      cli, name, (std::underlying_type_t<T>&)value, usage, choices, alt, req);
}
template <typename T, typename>
inline void add_argument(cli_command& cli, const string& name, T& value,
    const string& usage, const vector<string>& choices, bool req) {
  return add_argument(
      cli, name, (std::underlying_type_t<T>&)value, usage, choices, req);
}

// Backward compatibility
using cli_state [[deprecated]] = cli_command;

}  // namespace yocto

#endif
