//
// Implementation for Yocto/Json
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

#include "yocto_json.h"

#include <fstream>

#include "ext/json.hpp"

// -----------------------------------------------------------------------------
// FILE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Opens a file with a utf8 file name
static FILE* fopen_utf8(const char* filename, const char* mode) {
#ifdef _WIN32
  auto path8    = std::filesystem::u8path(filename);
  auto str_mode = string{mode};
  auto wmode    = std::wstring(str_mode.begin(), str_mode.end());
  return _wfopen(path8.c_str(), wmode.c_str());
#else
  return fopen(filename, mode);
#endif
}

// Load a text file
static bool load_text(const string& filename, string& str, string& error) {
  // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
  auto fs = fopen_utf8(filename.c_str(), "rb");
  if (!fs) {
    error = filename + ": file not found";
    return false;
  }
  fseek(fs, 0, SEEK_END);
  auto length = ftell(fs);
  fseek(fs, 0, SEEK_SET);
  str.resize(length);
  if (fread(str.data(), 1, length, fs) != length) {
    fclose(fs);
    error = filename + ": read error";
    return false;
  }
  fclose(fs);
  return true;
}

// Save a text file
static bool save_text(
    const string& filename, const string& str, string& error) {
  auto fs = fopen_utf8(filename.c_str(), "wt");
  if (!fs) {
    error = filename + ": file not found";
    return false;
  }
  if (fprintf(fs, "%s", str.c_str()) < 0) {
    fclose(fs);
    error = filename + ": write error";
    return false;
  }
  fclose(fs);
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON I/O
// -----------------------------------------------------------------------------
namespace yocto {

// Load json
bool load_json(const string& filename, json_value& json, string& error) {
  // load file
  auto text = string{};
  if (!load_text(filename, text, error)) return false;
  // parse
  return parse_json(text, json, error);
}

static void from_json(const nlohmann::json& js, json_value& value) {
  switch (js.type()) {
    case nlohmann::json::value_t::null: {
      value.set_null();
    } break;
    case nlohmann::json::value_t::number_integer: {
      value.set_integer((int64_t)js);
    } break;
    case nlohmann::json::value_t::number_unsigned: {
      value.set_uinteger((uint64_t)js);
    } break;
    case nlohmann::json::value_t::number_float: {
      value.set_number((double)js);
    } break;
    case nlohmann::json::value_t::boolean: {
      value.set_boolean((bool)js);
    } break;
    case nlohmann::json::value_t::string: {
      value.set_string((string)js);
    } break;
    case nlohmann::json::value_t::array: {
      value.set_array();
      for (auto& jitem : js) from_json(jitem, value.get_array().emplace_back());
    } break;
    case nlohmann::json::value_t::object: {
      value.set_object();
      for (auto& [key, jitem] : js.items())
        from_json(jitem, value.get_object()[key]);
    } break;
    case nlohmann::json::value_t::binary: {
      value.set_null();
    } break;
    case nlohmann::json::value_t::discarded: {
      value.set_null();
    } break;
  }
}

// Pars json
bool parse_json(const string& text, json_value& json, string& error) {
  try {
    auto njson = nlohmann::json::parse(text);
    from_json(njson, json);
    return true;
  } catch (...) {
    error = "error parsing json";
    return false;
  }
}

// Save json
bool save_json(const string& filename, json_value& json, string& error) {
  // convert to string
  auto text = string{};
  if (!dump_json(text, json, error)) return false;
  // save file
  return save_text(filename, text, error);
}

static void to_json(nlohmann::json& js, const json_value& value) {
  switch (value.get_type()) {
    case json_type::null: {
      js = {};
    } break;
    case json_type::integer: {
      js = value.get_integer();
    } break;
    case json_type::uinteger: {
      js = value.get_uinteger();
    } break;
    case json_type::number: {
      js = value.get_number();
    } break;
    case json_type::boolean: {
      js = value.get_boolean();
    } break;
    case json_type::string: {
      js = value.get_string();
    } break;
    case json_type::array: {
      js = nlohmann::json::array();
      for (auto& item : value.get_array()) to_json(js.emplace_back(), item);
    } break;
    case json_type::object: {
      js = nlohmann::json::object();
      for (auto& [key, item] : value.get_object()) to_json(js[key], item);
    } break;
  }
}

// Dump json
bool dump_json(string& text, const json_value& json, string& error) {
  try {
    auto njson = nlohmann::json{};
    to_json(njson, json);
    text = njson.dump(2);
    return true;
  } catch (...) {
    error = "error dumping json";
    return false;
  }
}

}  // namespace yocto