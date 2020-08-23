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

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF COMMAND-LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// initialize a command line parser
cli_state make_cli(const string& cmd, const string& usage) {
  auto cli  = cli_state{};
  cli.name  = cmd;
  cli.usage = usage;
  add_option(cli, "--help/--no-help", cli.help, "Print usage.");
  return cli;
}

static vector<string> split_cli_names(const string& name_) {
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

static void validate_names(const cli_state& cli) {
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
}

bool get_help(const cli_state& cli) { return cli.help; }

string get_usage(const cli_state& cli) {
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
      switch (value.type) {
        case cli_type::integer:
          str += option.choices.empty() ? std::to_string(value.number)
                                        : option.choices[value.integer];
          break;
        case cli_type::uinteger:
          str += option.choices.empty() ? std::to_string(value.number)
                                        : option.choices[value.uinteger];
          break;
        case cli_type::number: str += std::to_string(value.number); break;
        case cli_type::string: str += value.text; break;
        case cli_type::boolean: str += value.integer ? "true" : "false"; break;
      }
    }
    str += "]";
    return str;
  };
  auto message      = string{};
  auto has_optional = false, has_positional = false;
  auto usage_optional = string{}, usage_positional = string{};
  for (auto& option : cli.options) {
    auto line = "  " + option.name + " " + type_name(option);
    while (line.size() < 32) line += " ";
    line += option.usage;
    line += " " + def_string(option) + "\n";
    if (!option.choices.empty()) {
      line += "    with choices: ";
      auto len = 16;
      for (auto& choice : option.choices) {
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
    if (option.name.find("-") == 0) {
      has_optional = true;
      usage_optional += line;
    } else {
      has_positional = true;
      usage_positional += line;
    }
  }
  message += "usage: " + cli.name + (has_optional ? "" : " [options]") +
             (has_positional ? "" : " <arguments>") + cli.usage + "\n\n";
  if (has_optional) {
    message += "options:\n" + usage_optional + "\n";
  }
  if (has_positional) {
    message += "arguments:\n" + usage_positional + "\n";
  }
  return message;
}

static bool parse_value(
    cli_value& value, const string& arg, const vector<string>& choices) {
  if (!choices.empty()) {
    if (std::find(choices.begin(), choices.end(), arg) == choices.end())
      return false;
  }
  switch (value.type) {
    case cli_type::string: {
      value.text = arg;
      return true;
    } break;
    case cli_type::boolean: {
      if (arg == "true" || arg == "1") {
        value.integer = 1;
        return true;
      } else if (arg == "false" || arg == "0") {
        value.integer = 0;
        return true;
      } else {
        return false;
      }
    } break;
    case cli_type::integer: {
      if (choices.empty()) {
        auto end      = (char*)nullptr;
        value.integer = (int)strtol(arg.c_str(), &end, 10);
        return end != nullptr;
      } else {
        value.integer = (int64_t)(
            std::find(choices.begin(), choices.end(), arg) - choices.begin());
        return true;
      }
    } break;
    case cli_type::uinteger: {
      if (choices.empty()) {
        auto end       = (char*)nullptr;
        value.uinteger = (int)strtoul(arg.c_str(), &end, 10);
        return end != nullptr;
      } else {
        value.uinteger = (uint64_t)(
            std::find(choices.begin(), choices.end(), arg) - choices.begin());
        return true;
      }
    } break;
    case cli_type::number: {
      auto end     = (char*)nullptr;
      value.number = strtod(arg.c_str(), &end);
      return end != nullptr;
      return true;
    } break;
  }
  return false;
}

bool parse_cli(cli_state& cli, int argc, const char** argv, string& error) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  // validate names
  validate_names(cli);

  // prepare args
  auto args = vector<string>{argv + 1, argv + argc};
  // parse options
  for (auto& option : cli.options) {
    if (option.name[0] != '-') continue;
    option.value = option.def;
    option.set   = false;
    auto values  = vector<string>{};
    for (auto& name : split_cli_names(option.name)) {
      if (std::find(args.begin(), args.end(), name) == args.end()) continue;
      auto pos = std::find(args.begin(), args.end(), name) - args.begin();
      args.erase(args.begin() + pos);
      if (option.nargs == 0) {
        values     = {name.find("--no-") == string::npos ? "true" : "false"};
        option.set = true;
      } else if (option.nargs > 0) {
        if (pos + option.nargs > args.size())
          return cli_error("missing value for " + name);
        values     = {args.begin() + pos, args.begin() + pos + option.nargs};
        option.set = true;
        args.erase(args.begin() + pos, args.begin() + pos + option.nargs);
      } else {
        throw std::invalid_argument{"unsupported number of arguments"};
      }
    }
    if (option.set) {
      option.value.clear();
      for (auto& value : values) {
        option.value.emplace_back();
        option.value.back().type = option.type;
        if (!parse_value(option.value.back(), value, option.choices))
          return cli_error("bad value for " + option.name);
      }
      option.set_reference(option.value);
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
    option.value = option.def;
    option.set   = false;
    auto values  = vector<string>{};
    if (args.empty()) {
      if (option.req) return cli_error("missing value for " + option.name);
    } else if (option.nargs < 0) {
      values     = args;
      option.set = true;
      args.clear();
    } else if (option.nargs > 0) {
      if (option.nargs > args.size())
        return cli_error("missing value for " + option.name);
      values = {args.begin(), args.begin() + option.nargs};
      args.erase(args.begin(), args.begin() + option.nargs);
      option.set = true;
    } else {
      throw std::invalid_argument{"unsupported number of arguments"};
    }
    if (option.set) {
      option.value.clear();
      for (auto& value : values) {
        option.value.emplace_back();
        option.value.back().type = option.type;
        if (!parse_value(option.value.back(), value, option.choices))
          return cli_error("bad value for " + option.name);
      }
      option.set_reference(option.value);
    } else {
      if (option.req) return cli_error("missing value for " + option.name);
    }
  }
  // check remaining
  if (!args.empty()) return cli_error("mismatched value for " + args.front());
  // done
  return true;
}

void parse_cli(cli_state& cli, int argc, const char** argv) {
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
