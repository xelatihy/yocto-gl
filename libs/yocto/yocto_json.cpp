//
// Implementation for Yocto/JSON
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

#include "yocto_json.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <string>

#include "ext/json.hpp"
#include "yocto_commonio.h"

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF JSON IO
// -----------------------------------------------------------------------------
namespace yocto {

#define YOCTO_JSON_SAX 1

using njson = nlohmann::ordered_json;

// load/save json
bool load_json(const string& filename, njson& json, string& error) {
  // error helpers
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error in json";
    return false;
  };
  auto text = ""s;
  if (!load_text(filename, text, error)) return false;
  try {
    json = njson::parse(text);
    return true;
  } catch (std::exception&) {
    return parse_error();
  }
}

bool save_json(const string& filename, const njson& json, string& error) {
  return save_text(filename, json.dump(2), error);
}

// convert json
void to_json(json_value& json, const njson& njs) {
  switch (njs.type()) {
    case njson::value_t::null: json = json_value{}; break;
    case njson::value_t::number_integer: json = (int64_t)njs; break;
    case njson::value_t::number_unsigned: json = (uint64_t)njs; break;
    case njson::value_t::number_float: json = (double)njs; break;
    case njson::value_t::boolean: json = (bool)njs; break;
    case njson::value_t::string: json = (string)njs; break;
    case njson::value_t::array:
      json = json_array();
      for (auto& ejs : njs) to_json(json.emplace_back(), ejs);
      break;
    case njson::value_t::object:
      json = json_object();
      for (auto& [key, ejs] : njs.items()) to_json(json[key], ejs);
      break;
    case njson::value_t::binary:
      json                        = json_binary();
      json.get_ref<json_binary>() = njs.get_binary();
      break;
    case njson::value_t::discarded: json = json_value{}; break;
  }
}

// convert json
void from_json(const json_value& json, njson& njs) {
  switch (json.type()) {
    case json_type::null: njs = {}; break;
    case json_type::ninteger: njs = json.get_ref<int64_t>(); break;
    case json_type::nunsigned: njs = json.get_ref<uint64_t>(); break;
    case json_type::nfloat: njs = json.get_ref<double>(); break;
    case json_type::boolean: njs = json.get_ref<bool>(); break;
    case json_type::string: njs = json.get_ref<string>(); break;
    case json_type::array:
      njs = njson::array();
      for (auto& ejs : json) from_json(ejs, njs.emplace_back());
      break;
    case json_type::object:
      njs = njson::object();
      for (auto& [key, ejs] : json.items()) from_json(ejs, njs[key]);
      break;
    case json_type::binary:
      njs              = njson::binary({});
      njs.get_binary() = json.get_ref<json_binary>();
      break;
  }
}

#if YOCTO_JSON_SAX == 1

// load json
bool load_json(const string& filename, json_value& json, string& error) {
  // error helpers
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error in json";
    return false;
  };

  // sax handler
  struct sax_handler {
    // stack
    yocto::json_value*              root  = nullptr;
    std::vector<yocto::json_value*> stack = {};
    std::string                     current_key;
    explicit sax_handler(yocto::json_value* root_) {
      *root_ = yocto::json_value{};
      root   = root_;
      stack.push_back(root);
    }

    // get current value
    yocto::json_value& next_value() {
      if (stack.size() == 1) return *root;
      if (stack.back()->is_array()) return (*stack.back()).emplace_back();
      if (stack.back()->is_object()) return (*stack.back())[current_key];
      throw yocto::json_error{"bad json type"};
    }

    // values
    bool null() {
      next_value() = yocto::json_value{};
      return true;
    }
    bool boolean(bool value) {
      next_value() = value;
      return true;
    }
    bool number_integer(int64_t value) {
      next_value() = value;
      return true;
    }
    bool number_unsigned(uint64_t value) {
      next_value() = value;
      return true;
    }
    bool number_float(double value, const std::string&) {
      next_value() = value;
      return true;
    }
    bool string(std::string& value) {
      next_value() = value;
      return true;
    }
    bool binary(std::vector<uint8_t>& value) {
      next_value() = value;
      return true;
    }

    // objects
    bool start_object(size_t elements) {
      next_value() = yocto::json_object{};
      stack.push_back(&next_value());
      return true;
    }
    bool end_object() {
      stack.pop_back();
      return true;
    }
    bool key(std::string& value) {
      current_key = value;
      return true;
    }

    // arrays
    bool start_array(size_t elements) {
      next_value() = yocto::json_array{};
      stack.push_back(&next_value());
      return true;
    }
    bool end_array() {
      stack.pop_back();
      return true;
    }

    bool parse_error(size_t position, const std::string& last_token,
        const nlohmann::detail::exception&) {
      return false;
    }
  };

  // set up parsing
  json         = json_value{};
  auto handler = sax_handler{&json};

  // load text
  auto text = ""s;
  if (!load_text(filename, text, error)) return false;

  // parse json
  if (!njson::sax_parse(text, &handler)) return parse_error();
  return true;
}

#else

// load json
bool load_json(const string& filename, json_value& json, string& error) {
  // parse json
  auto njs = njson{};
  if (!load_json(filename, njs, error)) return false;

  // convert
  to_json(json, njs);
  return true;
}

#endif

// save json
bool save_json(const string& filename, const json_value& json, string& error) {
  // convert
  auto njs = njson{};
  from_json(json, njs);

  // save
  return save_json(filename, njs, error);
}

// Formats a Json to string
bool format_json(string& text, const json_value& json, string& error) {
  // convert
  auto njs = njson{};
  from_json(json, njs);

  // save
  text = njs.dump(2);
  return true;
}
string format_json(const json_value& json) {
  auto text  = string{};
  auto error = string{};
  if (!format_json(text, json, error)) return "";
  return text;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF JSON VALIDATION
// -----------------------------------------------------------------------------
namespace yocto {

// Validate a value against a schema
static bool validate_json(const json_value& value, const string& path,
    const json_value& schema, vector<string>& errors, size_t max_error) {
  // error handling
  auto emit_error = [&errors, max_error, &path](const string& message) {
    errors.push_back(message + (path.empty() ? ""s : ("at " + path)));
    return errors.size() >= max_error;
  };

  // early exit
  if (schema.is_boolean() && schema.get<bool>()) return true;
  if (schema.is_object() && schema.empty()) return true;

  // validate type
  if (schema.contains("type") && schema.at("type").is_string()) {
    auto& type    = schema.at("type").get_ref<string>();
    auto  type_ok = (type == "null" && value.is_null()) ||
                   (type == "integer" && value.is_integer()) ||
                   (type == "number" && value.is_number()) ||
                   (type == "boolean" && value.is_boolean()) ||
                   (type == "string" && value.is_string()) ||
                   (type == "array" && value.is_array()) ||
                   (type == "object" && value.is_object());
    if (!type_ok) {
      if (!emit_error(type + " expected")) return false;
    }
  }
  if (schema.contains("type") && schema.at("type").is_array()) {
    auto type_ok = false;
    for (auto& tschema : schema) {
      if (type_ok) break;
      auto& type = tschema.get_ref<string>();
      type_ok    = (type == "null" && value.is_null()) ||
                (type == "integer" && value.is_integer()) ||
                (type == "number" && value.is_number()) ||
                (type == "boolean" && value.is_boolean()) ||
                (type == "string" && value.is_string()) ||
                (type == "array" && value.is_array()) ||
                (type == "object" && value.is_object());
    }
    if (!type_ok) {
      auto types = ""s;
      for (auto& tschema : schema)
        types += (types.empty() ? "" : " or") + tschema.get_ref<string>();
      if (!emit_error(types + " expected")) return false;
    }
  }

  // check range
  // TODO(fabio): fix number precision
  if (schema.contains("minimum") && value.is_number()) {
    if (schema.at("minimum").get<double>() > value.get<double>()) {
      if (!emit_error("value out of range")) return false;
    }
  }
  if (schema.contains("maximum") && value.is_number()) {
    if (schema.at("maximum").get<double>() > value.get<double>()) {
      if (!emit_error("value out of range")) return false;
    }
  }
  if (schema.contains("exclusiveMinimum") && value.is_number()) {
    if (schema.at("exclusiveMinimum").get<double>() >= value.get<double>()) {
      if (!emit_error("value out of range")) return false;
    }
  }
  if (schema.contains("exclusiveMaximum") && value.is_number()) {
    if (schema.at("exclusiveMaximum").get<double>() <= value.get<double>()) {
      if (!emit_error("value out of range")) return false;
    }
  }

  // enum checks
  if (schema.contains("enum") && schema.at("enum").is_array()) {
    auto found = false;
    for (auto& item : schema.at("enum")) {
      if (found) break;
      if (item.is_string() && value.is_string() &&
          item.get_ref<string>() == value.get_ref<string>())
        found = true;
      if (item.is_integer() && value.is_integer() &&
          item.get<int64_t>() == value.get<int64_t>())
        found = true;
      if (item.is_number() && value.is_number() &&
          item.get<double>() == value.get<double>())
        found = true;
    }
    if (!found) {
      if (!emit_error("invalid enum")) return false;
    }
  }

  // size checks
  if (schema.contains("minLength") && value.is_string()) {
    if (schema.at("minLength").get<size_t>() > value.get_ref<string>().size()) {
      if (!emit_error("size out of range")) return false;
    }
  }
  if (schema.contains("maxLength") && value.is_string()) {
    if (schema.at("maxLength").get<size_t>() < value.get_ref<string>().size()) {
      if (!emit_error("size out of range")) return false;
    }
  }
  if (schema.contains("minItems") && value.is_array()) {
    if (schema.at("minItems").get<size_t>() >
        value.get_ref<json_array>().size()) {
      if (!emit_error("size out of range")) return false;
    }
  }
  if (schema.contains("maxItems") && value.is_array()) {
    if (schema.at("maxItems").get<size_t>() <
        value.get_ref<json_array>().size()) {
      if (!emit_error("size out of range")) return false;
    }
  }
  if (schema.contains("minProperties") && value.is_object()) {
    if (schema.at("minProperties").get<size_t>() >
        value.get_ref<json_object>().size()) {
      if (!emit_error("size out of range")) return false;
    }
  }
  if (schema.contains("maxProperties") && value.is_object()) {
    if (schema.at("maxProperties").get<size_t>() <
        value.get_ref<json_object>().size()) {
      if (!emit_error("size out of range")) return false;
    }
  }

  // check array items
  if (schema.contains("items") && value.is_object() &&
      schema.at("items").is_object()) {
    auto& items = schema.at("items");
    for (auto idx = (size_t)0; idx < value.size(); idx++) {
      if (!validate_json(value.at(idx), path + "/" + std::to_string(idx), items,
              errors, max_error)) {
        if (errors.size() > max_error) break;
      }
    }
  }
  if (schema.contains("items") && value.is_array() &&
      schema.at("items").is_array()) {
    auto& items = schema.at("items").get_ref<json_array>();
    for (auto idx = (size_t)0; idx < std::min(items.size(), value.size());
         idx++) {
      if (!validate_json(value.at(idx), path + "/" + std::to_string(idx),
              items.at(idx), errors, max_error)) {
        if (errors.size() > max_error) break;
      }
    }
  }

  // check object properties
  if (schema.contains("properties") && value.is_object() &&
      schema.at("properties").is_object()) {
    auto& properties = schema.at("properties").get_ref<json_object>();
    for (auto& [name, property] : properties) {
      if (!value.contains(name)) continue;
      if (!validate_json(
              value.at(name), path + "/" + name, property, errors, max_error)) {
        if (errors.size() > max_error) break;
      }
    }
  }
  if (schema.contains("additionalProperties") && value.is_object() &&
      schema.contains("properties") &&
      schema.at("additionalProperties").is_boolean() &&
      schema.at("additionalProperties").get<bool>() == false) {
    auto& properties = schema.at("properties");
    for (auto& [name, item] : value.get_ref<json_object>()) {
      if (properties.contains(name)) {
        if (!emit_error("unknown property " + name)) return false;
      }
    }
  }
  if (schema.contains("additionalProperties") && value.is_object() &&
      schema.contains("properties") &&
      schema.at("additionalProperties").is_object()) {
    auto& properties = schema.at("properties");
    for (auto& [name, item] : value.get_ref<json_object>()) {
      if (properties.contains(name)) continue;
      if (!validate_json(
              item, path + "/" + name, properties, errors, max_error)) {
        if (errors.size() > max_error) break;
      }
    }
  }
  if (schema.contains("required") && value.is_object() &&
      schema.at("required").is_array()) {
    auto& required = schema.at("required").get_ref<json_array>();
    for (auto& name_ : required) {
      auto& name = name_.get_ref<string>();
      if (!value.contains(name)) {
        if (emit_error("missing value for " + name)) return false;
      }
    }
  }

  // done
  return false;
}
bool validate_json(
    const json_value& value, const json_value& schema, string& error) {
  auto errors = vector<string>{};
  if (validate_json(value, "", schema, errors, 1)) return true;
  error = errors.at(0);
  return false;
}
bool validate_json(const json_value& value, const json_value& schema,
    vector<string>& errors, size_t max_errors) {
  return validate_json(value, "", schema, errors, max_errors);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF COMMAND-LINE PARSING
// -----------------------------------------------------------------------------
namespace yocto {

static json_value fix_cli_schema(const json_value& schema) { return schema; }

static string get_cliusage(
    const json_value& schema, const string& app_name, const string& command) {
  // helper
  auto is_positional = [](const json_value& schema,
                           const string&    name) -> bool {
    if (!schema.contains("cli_positional")) return false;
    if (!schema.at("cli_positional").is_array()) return false;
    for (auto& pname : schema.at("cli_positional")) {
      if (pname.is_string() && pname.get_ref<string>() == name) return true;
    }
    return false;
  };
  auto is_required = [](const json_value& schema, const string& name) -> bool {
    if (!schema.contains("required")) return false;
    if (!schema.at("required").is_array()) return false;
    for (auto& pname : schema.at("required")) {
      if (pname.is_string() && pname.get_ref<string>() == name) return true;
    }
    return false;
  };
  auto has_commands = [](const json_value& schema) -> bool {
    for (auto& [name, property] : schema.at("properties").items()) {
      if (property.value("type", "") == "object") return true;
    }
    return false;
  };
  auto get_alternate = [](const json_value& schema,
                           const string&    alt) -> string {
    if (!schema.contains("cli_alternate")) return "";
    if (!schema.at("cli_alternate").is_object()) return "";
    if (!schema.at("cli_alternate").contains(alt)) return "";
    if (!schema.at("cli_alternate").at(alt).is_string()) return "";
    return schema.at("cli_alternate").at(alt).get<string>();
  };

  auto message        = string{};
  auto usage_optional = string{}, usage_positional = string{},
       usage_command = string{};
  for (auto& [name, property] : schema.at("properties").items()) {
    if (property.value("type", "") == "object") continue;
    auto decorated_name = name;
    auto positional     = is_positional(schema, name);
    if (!positional) {
      decorated_name = "--" + name;
      if (auto alt = get_alternate(schema, name); !alt.empty())
        decorated_name += ", -" + alt;
    }
    auto line = "  " + decorated_name;
    if (property.value("type", "") != "boolean") {
      line += " <" + property.value("type", "") + ">";
    }
    while (line.size() < 32) line += " ";
    line += property.value("description", "");
    if (is_required(schema, name)) {
      line += " [req]\n";
    } else if (property.contains("default")) {
      auto& default_ = property.at("default");
      if (default_.is_array()) {
        line += " [";
        for (auto& item : default_) {
          if (&item != &default_.at(0)) line += ", ";
          line += format_json(item);
        }
        line += "]\n";
      } else if (default_.is_object()) {
        line += " [object]\n";
      } else {
        line += " [" + format_json(property.at("default")) + "]\n";
      }
    } else {
      line += "\n";
    }
    if (property.contains("enum")) {
      line += "    with choices: ";
      auto len = 16;
      for (auto& choice_ : property.at("enum")) {
        auto choice = format_json(choice_);
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
    if (positional) {
      usage_positional += line;
    } else {
      usage_optional += line;
    }
  }
  if (has_commands(schema)) {
    for (auto& [name, property] : schema.at("properties").items()) {
      if (property.value("type", "") != "object") continue;
      auto line = "  " + name;
      while (line.size() < 32) line += " ";
      line += property.value("description", "") + "\n";
      usage_command += line;
    }
  }

  {
    auto line = string{};
    line += "  --help";
    while (line.size() < 32) line += " ";
    line += "Prints an help message\n";
    usage_optional += line;
  }

  message += "usage: " + path_basename(app_name);
  if (!command.empty()) message += " " + command;
  if (!usage_command.empty()) message += " command";
  if (!usage_optional.empty()) message += " [options]";
  if (!usage_positional.empty()) message += " <arguments>";
  message += "\n";
  message += schema.value("description", "") + "\n\n";
  if (!usage_command.empty()) {
    message += "commands:\n" + usage_command + "\n";
  }
  if (!usage_optional.empty()) {
    message += "options:\n" + usage_optional + "\n";
  }
  if (!usage_positional.empty()) {
    message += "arguments:\n" + usage_positional + "\n";
  }
  return message;
}

string get_command(const cli_state& cli) { return cli.command; }
bool   get_help(const cli_state& cli) { return cli.help; }
string get_usage(const cli_state& cli) { return cli.usage; }

static bool parse_clivalue(
    json_value& value, const string& arg, const json_value& schema) {
  // if (!choices.empty()) {
  //   if (std::find(choices.begin(), choices.end(), arg) == choices.end())
  //     return false;
  // }
  auto type = schema.value("type", "string");
  if (type == "string") {
    value = arg;
    return true;
  } else if (type == "integer") {
    auto end = (char*)nullptr;
    if (arg.find('-') == 0) {
      value = (int64_t)strtol(arg.c_str(), &end, 10);
    } else {
      value = (uint64_t)strtoul(arg.c_str(), &end, 10);
    }
    return end != nullptr;
  } else if (type == "number") {
    auto end = (char*)nullptr;
    value    = strtod(arg.c_str(), &end);
    return end != nullptr;
    return true;
  } else if (type == "boolean") {
    if (arg == "true" || arg == "1") {
      value = true;
      return true;
    } else if (arg == "false" || arg == "0") {
      value = false;
      return true;
    } else {
      return false;
    }
  }
  return false;
}

static bool parse_clivalue(
    json_value& value, const vector<string>& args, const json_value& schema) {
  auto type = schema.value("type", "string");
  if (type == "array") {
    value = json_array{};
    for (auto& arg : args) {
      if (!parse_clivalue(value.emplace_back(), arg, schema.at("items")))
        return false;
    }
    return true;
  }
  return false;
}

static const char* cli_help_message = "Help invoked";

bool parse_cli(json_value& value, const json_value& schema_,
    const vector<string>& args, string& error, string& usage, string& command) {
  auto cli_error = [&error](const string& message) {
    error = message;
    return false;
  };

  // helpers
  auto advance_positional = [](const json_value& schema,
                                size_t&          last_positional) -> string {
    if (!schema.contains("cli_positional")) return "";
    auto& positionals = schema.at("cli_positional");
    if (!positionals.is_array()) return "";
    if (positionals.size() == last_positional) return "";
    if (!positionals.at(last_positional).is_string()) return "";
    return positionals.at(last_positional++).get<string>();
  };
  auto is_positional = [](const json_value& schema,
                           const string&    name) -> bool {
    if (!schema.contains("cli_positional")) return false;
    if (!schema.at("cli_positional").is_array()) return false;
    for (auto& pname : schema.at("cli_positional")) {
      if (pname.is_string() && pname.get_ref<string>() == name) return true;
    }
    return false;
  };
  auto is_required = [](const json_value& schema, const string& name) -> bool {
    if (!schema.contains("required")) return false;
    if (!schema.at("required").is_array()) return false;
    for (auto& pname : schema.at("required")) {
      if (pname.is_string() && pname.get_ref<string>() == name) return true;
    }
    return false;
  };
  auto get_alternate = [](const json_value& schema,
                           const string&    alt) -> string {
    if (!schema.contains("cli_alternate")) return "";
    if (!schema.at("cli_alternate").is_object()) return "";
    if (!schema.at("cli_alternate").contains(alt)) return "";
    if (!schema.at("cli_alternate").at(alt).is_string()) return "";
    return schema.at("cli_alternate").at(alt).get<string>();
  };
  auto get_command = [](const json_value& schema) -> string {
    if (!schema.contains("cli_command")) return "$command";
    if (!schema.at("cli_command").is_string()) return "$command";
    return schema.at("cli_command").get<string>();
  };
  auto has_commands = [](const json_value& schema) -> bool {
    for (auto& [name, property] : schema.at("properties").items()) {
      if (property.value("type", "") == "object") return true;
    }
    return false;
  };

  // parsing stack
  struct stack_elem {
    string      name = "";
    json_value& schema;
    json_value& value;
    size_t      positional = 0;
  };

  // initialize parsing
  auto schema = fix_cli_schema(schema_);
  value       = json_object{};
  auto stack  = vector<stack_elem>{{"", schema, value, 0}};
  command     = "";
  usage       = get_cliusage(schema, args[0], command);

  // parse the command line
  for (auto idx = (size_t)1; idx < args.size(); idx++) {
    auto& [_, schema, value, cpositional] = stack.back();
    auto arg                              = args.at(idx);
    auto positional                       = arg.find('-') != 0;
    if (positional && has_commands(schema)) {
      auto name = string{};
      for (auto& [pname, property] : schema.at("properties").items()) {
        if (property.value("type", "string") != "object") continue;
        if (pname != arg) continue;
        name = arg;
      }
      if (name.empty()) return cli_error("missing value for command");
      value[get_command(schema)] = name;
      value[name]                = json_object{};
      stack.push_back(
          {name, schema.at("properties").at(name), value.at(name), 0});
      command += (command.empty() ? "" : " ") + name;
      usage = get_cliusage(stack.back().schema, args[0], command);
      continue;
    } else if (positional) {
      auto name            = string{};
      auto next_positional = advance_positional(schema, cpositional);
      for (auto& [pname, property] : schema.at("properties").items()) {
        if (property.value("type", "string") == "object") continue;
        if (pname != next_positional) continue;
        name = pname;
      }
      if (name.empty()) return cli_error("too many positional arguments");
      auto& property = schema.at("properties").at(name);
      if (property.value("type", "string") == "array") {
        auto array_args = vector<string>(args.begin() + idx, args.end());
        if (!parse_clivalue(value[name], array_args, property))
          return cli_error("bad value for " + name);
        idx = args.size();
      } else if (property.value("type", "string") != "object") {
        if (!parse_clivalue(value[name], args[idx], property))
          return cli_error("bad value for " + name);
      }
    } else {
      if (arg == "--help" || arg == "-?") {
        return cli_error(cli_help_message);
      }
      arg = arg.substr(1);
      if (arg.find('-') == 0) arg = arg.substr(1);
      auto name = string{};
      for (auto& [pname, property] : schema.at("properties").items()) {
        if (property.value("type", "string") == "object") continue;
        if (is_positional(schema, pname)) continue;
        if (pname != arg && get_alternate(schema, pname) != arg) continue;
        name = pname;
        break;
      }
      if (name.empty()) return cli_error("unknown option " + args[idx]);
      if (value.contains(name)) return cli_error("option already set " + name);
      auto& property = schema.at("properties").at(name);
      if (property.value("type", "string") == "boolean") {
        if (!parse_clivalue(value[name], "true", property))
          return cli_error("bad value for " + name);
      } else if (property.value("type", "string") == "array") {
        auto min_items = (int)property.value("minItems", 0);
        auto max_items = (int)property.value("maxItems", 0);
        if (min_items != max_items || max_items == 0)
          return cli_error("bad array for " + name);
        if (idx + max_items >= args.size())
          return cli_error("missing value for " + name);
        auto array_args = vector<string>(
            args.begin() + idx + 1, args.begin() + idx + 1 + max_items);
        if (!parse_clivalue(value[name], array_args, property))
          return cli_error("bad value for " + name);
        idx += 1 + max_items;
      } else {
        if (idx + 1 >= args.size())
          return cli_error("missing value for " + name);
        if (!parse_clivalue(value[name], args[idx + 1], property))
          return cli_error("bad value for " + name);
        idx += 1;
      }
    }
  }

  // check for required and apply defaults
  for (auto& [_, schema, value, __] : stack) {
    if (has_commands(schema) && !value.contains(get_command(schema)))
      return cli_error("missing value for " + get_command(schema));
    for (auto& [name, property] : schema.at("properties").items()) {
      if (property.value("type", "string") == "object") continue;
      if (is_required(schema, name) && !value.contains(name))
        return cli_error("missing value for " + name);
      if (property.contains("default") && !value.contains(name))
        value[name] = property.at("default");
    }
  }

  // done
  return true;
}

bool parse_cli(json_value& value, const json_value& schema,
    const vector<string>& args, string& error, string& usage) {
  auto command = string{};
  return parse_cli(value, schema, args, error, usage, command);
}

bool parse_cli(json_value& value, const json_value& schema, int argc,
    const char** argv, string& error, string& usage) {
  return parse_cli(value, schema, {argv, argv + argc}, error, usage);
}

void parse_cli(
    json_value& value, const json_value& schema, const vector<string>& args) {
  auto error = string{};
  auto usage = string{};
  if (!parse_cli(value, schema, args, error, usage)) {
    if (error != cli_help_message) {
      print_info("error: " + error);
      print_info("");
    }
    print_info(usage);
    exit(error != cli_help_message ? 1 : 0);
  }
}

void parse_cli(
    json_value& value, const json_value& schema, int argc, const char** argv) {
  return parse_cli(value, schema, {argv, argv + argc});
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF JSON TREE DATA STRUCTURE
// -----------------------------------------------------------------------------
namespace yocto {

// convert json
void to_json(njson& njs, json_ctview json) {
  switch (get_type(json)) {
    case json_type::null: njs = nullptr; break;
    case json_type::ninteger: njs = get_number_integer(json); break;
    case json_type::nunsigned: njs = get_number_unsigned(json); break;
    case json_type::nfloat: njs = get_number_real(json); break;
    case json_type::boolean: njs = get_boolean(json); break;
    case json_type::string: njs = get_string(json); break;
    case json_type::array:
      njs = njson::array();
      for (auto ejs : iterate_array(json)) to_json(njs.emplace_back(), ejs);
      break;
    case json_type::object:
      njs = njson::object();
      for (auto [key, ejs] : iterate_object(json))
        to_json(njs[string{key}], ejs);
      break;
    case json_type::binary:
      njs = njson::binary({});
      get_binary(json, njs.get_binary());
      break;
  }
}

// convert json
void from_json(const njson& njs, json_tview json) {
  switch (njs.type()) {
    case njson::value_t::null: set_null(json); break;
    case njson::value_t::number_integer: set_integer(json, (int64_t)njs); break;
    case njson::value_t::number_unsigned: set_unsigned(json, njs); break;
    case njson::value_t::number_float: set_real(json, njs); break;
    case njson::value_t::boolean: set_boolean(json, (bool)njs); break;
    case njson::value_t::string: set_string(json, (string)njs); break;
    case njson::value_t::array:
      set_array(json);
      for (auto& ejs : njs) from_json(ejs, append_element(json));
      break;
    case njson::value_t::object:
      set_object(json);
      for (auto& [key, ejs] : njs.items())
        from_json(ejs, insert_element(json, key));
      break;
    case njson::value_t::binary: set_binary(json, njs.get_binary()); break;
    case njson::value_t::discarded: set_null(json); break;
  }
}

#if YOCTO_JSON_SAX == 1

// load json
bool load_json(const string& filename, json_tree& json, string& error) {
  // error helpers
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error in json";
    return false;
  };

  // sax handler
  struct sax_handler {
    // stack
    json_tview              root;
    std::vector<json_tview> stack = {};
    std::string             current_key;
    explicit sax_handler(json_tview root_) : root{root_}, stack{root_} {}

    // get current value
    json_tview next_value() {
      if (stack.size() == 1) return root;
      auto& jst = _get_type(stack.back());
      if (jst == json_type::array) return append_element(stack.back());
      if (jst == json_type::object)
        return insert_element(stack.back(), current_key);
      throw yocto::json_error{"bad json type"};
    }

    // values
    bool null() {
      set_null(next_value());
      return true;
    }
    bool boolean(bool value) {
      set_boolean(next_value(), value);
      return true;
    }
    bool number_integer(int64_t value) {
      set_integer(next_value(), value);
      return true;
    }
    bool number_unsigned(uint64_t value) {
      set_unsigned(next_value(), value);
      return true;
    }
    bool number_float(double value, const std::string&) {
      set_real(next_value(), value);
      return true;
    }
    bool string(std::string& value) {
      set_string(next_value(), value);
      return true;
    }
    bool binary(std::vector<uint8_t>& value) {
      set_binary(next_value(), value);
      return true;
    }

    // objects
    bool start_object(size_t elements) {
      set_object(next_value());
      stack.push_back(next_value());
      return true;
    }
    bool end_object() {
      stack.pop_back();
      return true;
    }
    bool key(std::string& value) {
      current_key = value;
      return true;
    }

    // arrays
    bool start_array(size_t elements) {
      set_array(next_value());
      stack.push_back(next_value());
      return true;
    }
    bool end_array() {
      stack.pop_back();
      return true;
    }

    bool parse_error(size_t position, const std::string& last_token,
        const nlohmann::detail::exception&) {
      return false;
    }
  };

  // set up parsing
  json         = json_tree{};
  auto handler = sax_handler{get_root(json)};

  // load text
  auto text = ""s;
  if (!load_text(filename, text, error)) return false;

  // parse json
  if (!njson::sax_parse(text, &handler)) return parse_error();
  return true;
}

#else

// load json
bool load_json(const string& filename, json_tree& json, string& error) {
  // parse json
  auto njs = njson{};
  if (!load_json(filename, njs, error)) return false;

  // convert
  from_json(njs, get_root(json));
  return true;
}

#endif

// save json
bool save_json(const string& filename, const json_tree& json, string& error) {
  // convert
  auto njs = njson{};
  to_json(njs, get_root((json_tree&)json));

  // save
  return save_json(filename, njs, error);
}

}  // namespace yocto
