//
// # Yocto/SceneIO: Scene serialization
//
// Yocto/SceneIO supports loading and saving scenes from Ply, Obj, Pbrt, glTF
// and a custom Json format.
// Yocto/SceneIO is implemented in `yocto_sceneio.h` and `yocto_sceneio.cpp`
// and depends on `stb_image.h`, `stb_image_write.h`, `tinyexr.h`.
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

#ifndef _YOCTO_SCENEIO_H_
#define _YOCTO_SCENEIO_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <functional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "yocto_scene.h"

#define JSON_USE_IMPLICIT_CONVERSIONS 0
#include "ext/json.hpp"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IO ERROR
// -----------------------------------------------------------------------------
namespace yocto {

// Result object modeled on std::expected
struct io_error : std::runtime_error {
  using std::runtime_error::runtime_error;
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Check if an image is HDR or LDR based on filename.
bool is_hdr_filename(const string& filename);
bool is_ldr_filename(const string& filename);

// Loads/saves a 4 channels float/byte image in linear/srgb color space.
bool load_image(const string& filename, image_data& img, string& error);
bool save_image(const string& filename, const image_data& img, string& error);

// Loads/saves a 4 channels float/byte image in linear/srgb color space.
image_data load_image(const string& filename);
void       load_image(const string& filename, image_data& image);
void       save_image(const string& filename, const image_data& image);

// Make presets. Supported mostly in IO.
image_data make_image_preset(const string& type);

// Make presets. Supported mostly in IO.
bool make_image_preset(
    const string& filename, image_data& image, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// TEXTURE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a texture in the supported formats.
bool load_texture(const string& filename, texture_data& texture, string& error);
bool save_texture(
    const string& filename, const texture_data& texture, string& error);

// Load/save a texture in the supported formats.
texture_data load_texture(const string& filename);
void         load_texture(const string& filename, texture_data& texture);
void         save_texture(const string& filename, const texture_data& texture);

// Make presets. Supported mostly in IO.
texture_data make_texture_preset(const string& type);

// Make presets. Supported mostly in IO.
bool make_texture_preset(
    const string& filname, texture_data& texture, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a shape
bool load_shape(const string& filename, shape_data& shape, string& error,
    bool flip_texcoords = true);
bool save_shape(const string& filename, const shape_data& shape, string& error,
    bool flip_texcoords = true, bool ascii = false);

// Load/save a shape
shape_data load_shape(const string& filename, bool flip_texcoords = true);
void       load_shape(
          const string& filename, shape_data& shape, bool flip_texcoords = true);
void save_shape(const string& filename, const shape_data& shape,
    bool flip_texcoords = true, bool ascii = false);

// Load/save a subdiv
bool load_fvshape(const string& filename, fvshape_data& shape, string& error,
    bool flip_texcoords = true);
bool save_fvshape(const string& filename, const fvshape_data& shape,
    string& error, bool flip_texcoords = true, bool ascii = false);

// Load/save a subdiv
fvshape_data load_fvshape(const string& filename, bool flip_texcoords = true);
void         load_fvshape(
            const string& filename, fvshape_data& shape, bool flip_texcoords = true);
void save_fvshape(const string& filename, const fvshape_data& shape,
    bool flip_texcoords = true, bool ascii = false);

// Make presets. Supported mostly in IO.
shape_data   make_shape_preset(const string& type);
fvshape_data make_fvshape_preset(const string& type);

// Make presets. Supported mostly in IO.
bool make_shape_preset(const string& filname, shape_data& shape, string& error);
bool make_fvshape_preset(
    const string& filname, fvshape_data& shape, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SUBDIV IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a subdiv in the supported formats.
bool load_subdiv(const string& filename, subdiv_data& subdiv, string& error);
bool save_subdiv(
    const string& filename, const subdiv_data& subdiv, string& error);

// Load/save a subdiv in the supported formats.
subdiv_data load_subdiv(const string& filename);
void        load_subdiv(const string& filename, subdiv_data& subdiv);
void        save_subdiv(const string& filename, const subdiv_data& subdiv);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a scene in the supported formats.
bool load_scene(const string& filename, scene_data& scene, string& error,
    bool noparallel = false);
bool save_scene(const string& filename, const scene_data& scene, string& error,
    bool noparallel = false);

// Make missing scene directories
bool make_scene_directories(
    const string& filename, const scene_data& scene, string& error);

// Add environment
bool add_environment(scene_data& scene, const string& filename, string& error);

// Load/save a scene in the supported formats.
scene_data load_scene(const string& filename, bool noparallel = false);
void       load_scene(
          const string& filename, scene_data& scene, bool noparallel = false);
void save_scene(
    const string& filename, const scene_data& scene, bool noparallel = false);

// Add environment
void add_environment(scene_data& scene, const string& filename);

// Make missing scene directories
void make_scene_directories(const string& filename, const scene_data& scene);

// Scene presets used for testing.
scene_data make_scene_preset(const string& type);

// Scene presets used for testing.
bool make_scene_preset(
    const string& filename, scene_data& scene, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// FILE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Using directive
using byte = unsigned char;

// Load/save a text file
bool load_text(const string& filename, string& str, string& error);
bool save_text(const string& filename, const string& str, string& error);

// Load/save a binary file
bool load_binary(const string& filename, vector<byte>& data, string& error);
bool save_binary(
    const string& filename, const vector<byte>& data, string& error);

// Load/save a text file
string load_text(const string& filename);
void   load_text(const string& filename, string& str);
void   save_text(const string& filename, const string& str);

// Load/save a binary file
vector<byte> load_binary(const string& filename);
void         load_binary(const string& filename, vector<byte>& data);
void         save_binary(const string& filename, const vector<byte>& data);

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON SUPPORT
// -----------------------------------------------------------------------------
namespace yocto {

// Json values
using json_value = nlohmann::ordered_json;

// Load/save a json file
bool load_json(const string& filename, json_value& json, string& error);
bool save_json(const string& filename, const json_value& json, string& error);

// Load/save a json file
json_value load_json(const string& filename);
void       load_json(const string& filename, json_value& json);
void       save_json(const string& filename, const json_value& json);

// Json conversions
inline void to_json(json_value& json, const vec2f& value);
inline void to_json(json_value& json, const vec3f& value);
inline void to_json(json_value& json, const vec4f& value);
inline void to_json(json_value& json, const frame2f& value);
inline void to_json(json_value& json, const frame3f& value);
inline void to_json(json_value& json, const mat2f& value);
inline void to_json(json_value& json, const mat3f& value);
inline void to_json(json_value& json, const mat4f& value);
inline void from_json(const json_value& json, vec2f& value);
inline void from_json(const json_value& json, vec3f& value);
inline void from_json(const json_value& json, vec4f& value);
inline void from_json(const json_value& json, frame2f& value);
inline void from_json(const json_value& json, frame3f& value);
inline void from_json(const json_value& json, mat2f& value);
inline void from_json(const json_value& json, mat3f& value);
inline void from_json(const json_value& json, mat4f& value);

}  // namespace yocto

// -----------------------------------------------------------------------------
// ARGUMENT PARSING
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

// add option
struct cli_option;
template <typename T>
inline cli_option& add_option(
    cli_command& cli, const string& name, T& value, const string& usage);
inline cli_option& add_option(
    cli_command& cli, const string& name, bool& value, const string& usage);
template <typename T>
inline cli_option& add_option(cli_command& cli, const string& name,
    vector<T>& value, const string& usage);
template <typename T, size_t N>
inline cli_option& add_option(cli_command& cli, const string& name,
    array<T, N>& value, const string& usage);

// add option with labels
template <typename T>
inline cli_option& add_option(cli_command& cli, const string& name, T& value,
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
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// JSON MANIPULATION
// -----------------------------------------------------------------------------
namespace yocto {

// Json conversions
inline void to_json(json_value& json, const vec2f& value) {
  nlohmann::to_json(json, (const array<float, 2>&)value);
}
inline void to_json(json_value& json, const vec3f& value) {
  nlohmann::to_json(json, (const array<float, 3>&)value);
}
inline void to_json(json_value& json, const vec4f& value) {
  nlohmann::to_json(json, (const array<float, 4>&)value);
}
inline void to_json(json_value& json, const frame2f& value) {
  nlohmann::to_json(json, (const array<float, 6>&)value);
}
inline void to_json(json_value& json, const frame3f& value) {
  nlohmann::to_json(json, (const array<float, 12>&)value);
}
inline void to_json(json_value& json, const mat2f& value) {
  nlohmann::to_json(json, (const array<float, 4>&)value);
}
inline void to_json(json_value& json, const mat3f& value) {
  nlohmann::to_json(json, (const array<float, 9>&)value);
}
inline void to_json(json_value& json, const mat4f& value) {
  nlohmann::to_json(json, (const array<float, 16>&)value);
}
inline void from_json(const json_value& json, vec2f& value) {
  nlohmann::from_json(json, (array<float, 2>&)value);
}
inline void from_json(const json_value& json, vec3f& value) {
  nlohmann::from_json(json, (array<float, 3>&)value);
}
inline void from_json(const json_value& json, vec4f& value) {
  nlohmann::from_json(json, (array<float, 4>&)value);
}
inline void from_json(const json_value& json, frame2f& value) {
  nlohmann::from_json(json, (array<float, 6>&)value);
}
inline void from_json(const json_value& json, frame3f& value) {
  nlohmann::from_json(json, (array<float, 12>&)value);
}
inline void from_json(const json_value& json, mat2f& value) {
  nlohmann::from_json(json, (array<float, 4>&)value);
}
inline void from_json(const json_value& json, mat3f& value) {
  nlohmann::from_json(json, (array<float, 9>&)value);
}
inline void from_json(const json_value& json, mat4f& value) {
  nlohmann::from_json(json, (array<float, 16>&)value);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// ARGUMENT PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// Cli data
using cli_setter = std::function<bool(const vector<string>&, string&)>;
struct cli_option {
  // name
  string name = "";
  // setter
  cli_setter setter = {};
};
struct cli_command {
  // name
  string name = "";
  // options and commands
  unordered_map<string, cli_option>  options  = {};
  unordered_map<string, cli_command> commands = {};
  // command
  string     command = "";
  cli_setter setter  = {};
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
  if (cli.options.find(name) != cli.options.end())
    throw std::invalid_argument{"option already added " + name};
}
inline void _cli_check_command(const cli_command& cli, const string& name) {
  if (!cli.options.empty())
    throw std::invalid_argument{"cannot add options and commands"};
  if (cli.commands.find(name) != cli.commands.end())
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
  cli.name        = name;
  cli.usage_name  = name;
  cli.usage_descr = usage;
  cli.options.reserve(256);
  cli.commands.reserve(256);
  cli.setter = [](const vector<string>&, string&) { return true; };
  return cli;
}

// add command
inline cli_command& add_command(
    cli_command& cli, const string& name, const string& usage) {
  _cli_check_command(cli, name);
  auto& cmd       = cli.commands[name];
  cmd.name        = name;
  cmd.usage_name  = cli.usage_name + " " + name;
  cmd.usage_descr = usage;
  cli.usage_commands += _cli_usage_command(name, usage);
  cmd.options.reserve(256);
  cmd.commands.reserve(256);
  cmd.setter = [](const vector<string>&, string&) { return true; };
  return cmd;
}

// add command variable
template <typename T>
inline void add_command_var(cli_command& cli, T& value) {
  cli.setter = [&value](const vector<string>& args, string& error) {
    if (!_cli_parse_size(args, 1, 1, error)) return false;
    if (!_cli_parse_value(args.front(), value, error)) return false;
    return true;
  };
}

// add option
template <typename T>
inline cli_option& add_option(
    cli_command& cli, const string& name, T& value, const string& usage) {
  _cli_check_option(cli, name);
  auto& opt = cli.options[name];
  opt.name  = name;
  cli.usage_options += _cli_usage_option(name, value, usage);
  opt.setter = [&value](const vector<string>& args, string& error) {
    if (!_cli_parse_size(args, 1, 1, error)) return false;
    if (!_cli_parse_value(args.front(), value, error)) return false;
    return true;
  };
  return opt;
}

// add option
inline cli_option& add_option(
    cli_command& cli, const string& name, bool& value, const string& usage) {
  _cli_check_option(cli, name);
  auto& opt = cli.options[name];
  opt.name  = name;
  cli.usage_options += _cli_usage_option(name, value, usage);
  opt.setter = [&value](const vector<string>& args, string& error) {
    if (!_cli_parse_size(args, 0, 1, error)) return false;
    if (!_cli_parse_value(args.empty() ? "true" : args.front(), value, error))
      return false;
    return true;
  };
  return opt;
}

// add option
template <typename T>
inline cli_option& add_option(cli_command& cli, const string& name,
    vector<T>& value, const string& usage) {
  _cli_check_option(cli, name);
  auto& opt = cli.options[name];
  opt.name  = name;
  cli.usage_options += _cli_usage_option(name, value, usage);
  opt.setter = [&value](const vector<string>& args, string& error) {
    if (!_cli_parse_size(args, 1, 1024, error)) return false;
    value.resize(args.size());
    for (auto idx = (size_t)0; idx < args.size(); idx++) {
      if (!_cli_parse_value(args[idx], value[idx], error)) return false;
    }
    return true;
  };
  return opt;
}

// add option
template <typename T, size_t N>
inline cli_option& add_option(cli_command& cli, const string& name,
    array<T, N>& value, const string& usage) {
  _cli_check_option(cli, name);
  auto& opt = cli.options[name];
  opt.name  = name;
  cli.usage_options += _cli_usage_option(name, value, usage);
  opt.setter = [&value](const vector<string>& args, string& error) {
    if (!_cli_parse_size(args, N, N, error)) return false;
    for (auto idx = (size_t)0; idx < args.size(); idx++) {
      if (!_cli_parse_value(args[idx], value[idx], error)) return false;
    }
    return true;
  };
  return opt;
}

// add option
template <typename T>
inline cli_option& add_option(cli_command& cli, const string& name, T& value,
    const string& usage, const vector<pair<T, string>>& labels) {
  _cli_check_option(cli, name);
  auto& opt = cli.options[name];
  opt.name  = name;
  cli.usage_options += _cli_usage_option(name, value, usage, labels);
  opt.setter = [&value, labels](const vector<string>& args, string& error) {
    if (!_cli_parse_size(args, 1, 1, error)) return false;
    if (!_cli_parse_value(args.front(), value, error, labels)) return false;
    return true;
  };
  return opt;
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
    if (cli.commands.find(args[pos]) == cli.commands.end()) {
      error = "unknown command " + args[pos];
      return false;
    }
    // get command
    auto name   = args[pos++];
    cli.command = name;
    // set command
    if (!cli.setter({name}, error)) {
      error += " for command " + name;
      return false;
    }
    // parse command recursively
    auto& command = cli.commands.at(name);
    return parse_cli(command, args, pos, error);
  } else {
    // check option
    if (pos < args.size() - 1 && args[pos].find("--") != 0) {
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
      if (cli.options.find(args[pos].substr(2)) == cli.options.end()) {
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
      if (!option.setter(values, error)) {
        error += " for option " + name;
        return false;
      }
    }
    // done
    return true;
  }
}

inline string get_usage(const cli_command& cli) {
  if (!cli.command.empty()) return get_usage(cli.commands.at(cli.command));
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

#endif
