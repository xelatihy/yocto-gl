//
// Implementation for Yocto/CommonIO
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

#include "yocto_json.h"

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

// initialize a command line parser
cli_state make_cli(const string& name, const string& usage) {
  auto  cli             = cli_state{};
  auto& schema          = cli.schema;
  schema["title"]       = name;
  schema["description"] = usage;
  schema["type"]        = "object";
  return cli;
}

// add command
cli_command add_command(
    const cli_command& cli, const string& name, const string& usage) {
  auto& schema            = get_clischema(cli.cli.schema, cli.path);
  schema["cli_command"]   = "command";
  auto& property          = schema["properties"][name];
  property["title"]       = name;
  property["description"] = usage;
  property["type"]        = "object";
  property["properties"]  = json_object{};
  auto& setter            = get_clisetter(cli.cli.setter, cli.path);
  setter[name]            = cli_setter{};
  auto subcommand         = cli;
  subcommand.path.push_back(name);
  return subcommand;
}
cli_command add_command(
    cli_state& cli, const string& name, const string& usage) {
  return add_command({cli, {}}, name, usage);
}

static bool set_clivalues(
    const json_value& js, cli_setter& value, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  if (!value.is_object()) return cli_error("bad value");
  for (auto& [name, item] : value) {
    if (!js.contains(name)) continue;
    if (item.is_object()) {
      if (!set_clivalues(js.at(name), item, error)) return false;
    } else {
      if (!item.set(js.at(name))) return cli_error("bad value for " + name);
    }
  }
  return true;
}

static const char* cli_help_message = "Help invoked";

bool parse_cli(cli_state& cli, const vector<string>& args, string& error) {
  auto usage   = string{};
  auto command = string{};
  auto ok      = parse_cli(cli.value, cli.schema, args, error, usage, command);
  cli.usage    = usage;
  cli.command  = command;
  cli.help     = error == cli_help_message;
  if (!ok) return false;
  if (!set_clivalues(cli.value, cli.setter, error)) return false;
  return true;
}

void parse_cli(cli_state& cli, const vector<string>& args) {
  auto error = string{};
  if (!parse_cli(cli, args, error)) {
    if (error != cli_help_message) {
      print_info("error: " + error);
      print_info("");
    }
    print_info(get_usage(cli));
    exit(error != cli_help_message ? 1 : 0);
  }
}

bool parse_cli(cli_state& cli, int argc, const char** argv, string& error) {
  return parse_cli(cli, {argv, argv + argc}, error);
}

void parse_cli(cli_state& cli, int argc, const char** argv) {
  return parse_cli(cli, {argv, argv + argc});
}

}  // namespace yocto
