//
// Implementation for Yocto/CommonIO
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

#include "yocto_commonio.h"

#include <charconv>
#include <cstdio>
#include <filesystem>
#include <limits>

#include "ext/fast_float/fast_float.h"
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
vector<string> list_directory(const string& dirname) {
  try {
    auto entries = vector<string>{};
    for (auto entry : std::filesystem::directory_iterator(make_path(dirname))) {
      entries.push_back(entry.path().generic_u8string());
    }
    return entries;
  } catch (...) {
    throw io_error{dirname, "cannot list directory"};
  }
}

// Create a directory and all missing parent directories if needed
void make_directory(const string& dirname) {
  if (path_exists(dirname)) return;
  try {
    create_directories(make_path(dirname));
  } catch (...) {
    throw io_error{dirname, "cannot create directory"};
  }
}

// Get the current directory
string path_current() { return std::filesystem::current_path().u8string(); }

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

// Opens a file with utf8 filename
FILE* fopen_utf8(const string& filename, const string& mode) {
#ifdef _WIN32
  auto path8 = std::filesystem::u8path(filename);
  auto wmode = std::wstring(mode.begin(), mode.end());
  return _wfopen(path8.c_str(), wmode.c_str());
#else
  return fopen(filename.c_str(), mode.c_str());
#endif
}

// Load a text file
string load_text(const string& filename) {
  auto str = string{};
  load_text(filename, str);
  return str;
}
void load_text(const string& filename, string& text) {
  // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
  auto fs = fopen_utf8(filename.c_str(), "rb");
  if (!fs) throw io_error::open_error(filename);
  fseek(fs, 0, SEEK_END);
  auto length = ftell(fs);
  fseek(fs, 0, SEEK_SET);
  text.resize(length);
  if (fread(text.data(), 1, length, fs) != length) {
    fclose(fs);
    throw io_error::read_error(filename);
  }
  fclose(fs);
}

// Save a text file
void save_text(const string& filename, const string& text) {
  auto fs = fopen_utf8(filename.c_str(), "wt");
  if (!fs) throw io_error::open_error(filename);
  if (fprintf(fs, "%s", text.c_str()) < 0) {
    fclose(fs);
    throw io_error::write_error(filename);
  }
  fclose(fs);
}

// Load a binary file
vector<byte> load_binary(const string& filename) {
  auto data = vector<byte>{};
  load_binary(filename, data);
  return data;
}
void load_binary(const string& filename, vector<byte>& data) {
  // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
  auto fs = fopen_utf8(filename.c_str(), "rb");
  if (!fs) throw io_error::open_error(filename);
  fseek(fs, 0, SEEK_END);
  auto length = ftell(fs);
  fseek(fs, 0, SEEK_SET);
  data.resize(length);
  if (fread(data.data(), 1, length, fs) != length) {
    fclose(fs);
    throw io_error::read_error(filename);
  }
  fclose(fs);
}

// Save a binary file
void save_binary(const string& filename, const vector<byte>& data) {
  auto fs = fopen_utf8(filename.c_str(), "wb");
  if (!fs) throw io_error::open_error(filename);
  if (fwrite(data.data(), 1, data.size(), fs) != data.size()) {
    fclose(fs);
    throw io_error::write_error(filename);
  }
  fclose(fs);
}

// Load a text file
bool load_text(const string& filename, string& str, string& error) {
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
bool save_text(const string& filename, const string& str, string& error) {
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

// Load a binary file
bool load_binary(const string& filename, vector<byte>& data, string& error) {
  // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
  auto fs = fopen_utf8(filename.c_str(), "rb");
  if (!fs) {
    error = filename + ": file not found";
    return false;
  }
  fseek(fs, 0, SEEK_END);
  auto length = ftell(fs);
  fseek(fs, 0, SEEK_SET);
  data.resize(length);
  if (fread(data.data(), 1, length, fs) != length) {
    fclose(fs);
    error = filename + ": read error";
    return false;
  }
  fclose(fs);
  return true;
}

// Save a binary file
bool save_binary(
    const string& filename, const vector<byte>& data, string& error) {
  auto fs = fopen_utf8(filename.c_str(), "wb");
  if (!fs) {
    error = filename + ": file not found";
    return false;
  }
  if (fwrite(data.data(), 1, data.size(), fs) != data.size()) {
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
json_value load_json(const string& filename) {
  auto json = json_value{};
  load_json(filename, json);
  return json;
}
void load_json(const string& filename, json_value& json) {
  auto text = load_text(filename);
  parse_json(text, json);
}
void save_json(const string& filename, const json_value& json) {
  auto text = format_json(json);
  save_text(filename, text);
}

struct json_sax {
  explicit json_sax(json_value* root, bool allow_exceptions = true)
      : _root(root), _allow_exceptions(allow_exceptions) {}

  bool null() {
    handle_value(nullptr);
    return true;
  }

  bool boolean(bool val) {
    handle_value(val);
    return true;
  }

  bool number_integer(int64_t val) {
    handle_value(val);
    return true;
  }

  bool number_unsigned(uint64_t val) {
    handle_value(val);
    return true;
  }

  bool number_float(double val, const std::string& /*unused*/) {
    handle_value(val);
    return true;
  }

  bool string(std::string& val) {
    handle_value(val);
    return true;
  }

  bool binary(vector<byte>& val) {
    handle_value(std::move(val));
    return true;
  }

  bool start_object(std::size_t len) {
    _stack.push_back(handle_value(json_object{}));
    return true;
  }

  bool key(std::string& val) {
    // add null at given key and store the reference for later
    // _object_item = &(_stack.back()->operator[](val));
    _object_item = &(_stack.back()->insert_back(std::move(val)));
    return true;
  }

  bool end_object() {
    _stack.pop_back();
    return true;
  }

  bool start_array(std::size_t len) {
    _stack.push_back(handle_value(json_array{}));
    return true;
  }

  bool end_array() {
    _stack.pop_back();
    return true;
  }

  template <class Exception>
  bool parse_error(std::size_t /*unused*/, const std::string& /*unused*/,
      const Exception& ex) {
    _errored = true;
    static_cast<void>(ex);
    if (_allow_exceptions) {
      if (_stack.empty()) throw json_error{"parse error", _root};
      if (_stack.back()->is_array())
        throw json_error{"parse error", _stack.back()};
      if (_stack.back()->is_object())
        throw json_error{"parse error", _stack.back()};
    }
    return false;
  }

  bool is_errored() const { return _errored; }

 private:
  template <typename Value>
  json_value* handle_value(Value&& v) {
    if (_stack.empty()) {
      *_root = json_value(std::forward<Value>(v));
      return _root;
    }

    assert(_stack.back()->is_array() || _stack.back()->is_object());

    if (_stack.back()->is_array()) {
      _stack.back()->emplace_back(std::forward<Value>(v));
      return &(_stack.back()->back());
    }

    assert(_stack.back()->is_object());
    assert(_object_item);
    *_object_item = json_value(std::forward<Value>(v));
    return _object_item;
  }

  json_value*         _root;
  vector<json_value*> _stack            = {};
  json_value*         _object_item      = nullptr;
  bool                _errored          = false;
  bool                _allow_exceptions = true;
};

// Parse json
json_value parse_json(const string& text) {
  auto json = json_value{};
  parse_json(text, json);
  return json;
}
void parse_json(const string& text, json_value& json) {
  json     = json_value{};
  auto sax = json_sax{&json};
  nlohmann::json::sax_parse(text, &sax);
}

static void _json_write_char(string& text, char value) {
  text.push_back(value);
}
static void _json_write_chars(string& text, const char* value, int size) {
  text.append(value, size);
}
static void _json_write_indent(string& text, int size) {
  text.append(size, ' ');
}
static void _json_format_value(string& text, int64_t value) {
  char buffer[64];
  auto result = std::to_chars(buffer, buffer + 64, value);
  text.append(buffer, result.ptr - buffer);
}
static void _json_format_value(string& text, uint64_t value) {
  char buffer[64];
  auto result = std::to_chars(buffer, buffer + 64, value);
  text.append(buffer, result.ptr - buffer);
}
static void _json_format_value(string& text, double value) {
  char  buffer[128];
  char* end = ::nlohmann::detail::to_chars(buffer, buffer + 128, value);
  text.append(buffer, end - buffer);
}
static void _json_format_value(string& text, bool value) {
  if (value) {
    text.append("true", 4);
  } else {
    text.append("false", 5);
  }
}
static void _json_format_value(string& text, const string& value) {
  text.push_back('\"');
  for (auto c : value) {
    switch (c) {
      case '\"': text.append("\\\"", 2); break;
      case '\t': text.append("\\t", 2); break;
      case '\n': text.append("\\n", 2); break;
      case '\r': text.append("\\r", 2); break;
      case '\b': text.append("\\b", 2); break;
      case '\f': text.append("\\f", 2); break;
      case '\\': text.append("\\\\", 2); break;
      default: text.push_back(c); break;
    }
  }
  text.push_back('\"');
}

static void _json_format_element(
    string& text, const json_value& json, int indent) {
  switch (json.type()) {
    case json_type::null: {
      _json_write_chars(text, "null", 4);
    } break;
    case json_type::integer: {
      _json_format_value(text, json.get<int64_t>());
    } break;
    case json_type::uinteger: {
      _json_format_value(text, json.get<uint64_t>());
    } break;
    case json_type::number: {
      _json_format_value(text, json.get<double>());
    } break;
    case json_type::boolean: {
      _json_format_value(text, json.get<bool>());
    } break;
    case json_type::string: {
      _json_format_value(text, json.get<string>());
    } break;
    case json_type::array: {
      if (json.empty()) {
        _json_write_chars(text, "[]", 2);
      } else if (!json.at(0).is_array() && !json.at(0).is_object()) {
        _json_write_chars(text, "[ ", 2);
        auto count = (size_t)0, size = json.size();
        for (auto& item : json) {
          _json_format_element(text, item, indent + 2);
          count += 1;
          if (count < size) _json_write_chars(text, ", ", 2);
        }
        _json_write_chars(text, " ]", 2);
      } else {
        _json_write_chars(text, "[\n", 2);
        auto count = (size_t)0, size = json.size();
        for (auto& item : json) {
          _json_write_indent(text, indent + 2);
          _json_format_element(text, item, indent + 2);
          count += 1;
          if (count < size) _json_write_chars(text, ",\n", 2);
        }
        _json_write_char(text, '\n');
        _json_write_indent(text, indent);
        _json_write_char(text, ']');
      }
    } break;
    case json_type::object: {
      if (json.empty()) {
        _json_write_chars(text, "{}", 2);
      } else {
        _json_write_chars(text, "{\n", 2);
        auto count = (size_t)0, size = json.size();
        for (auto& [key, item] : json.items()) {
          _json_write_indent(text, indent + 2);
          _json_format_value(text, key);
          _json_write_chars(text, ": ", 2);
          _json_format_element(text, item, indent + 2);
          count += 1;
          if (count < size) _json_write_chars(text, ",\n", 2);
        }
        _json_write_char(text, '\n');
        _json_write_indent(text, indent);
        _json_write_char(text, '}');
      }
    }
  }
}

// Dump json
string format_json(const json_value& json) {
  auto text = string{};
  format_json(text, json);
  return text;
}
void format_json(string& text, const json_value& json) {
  text.clear();
  _json_format_element(text, json, 0);
}

// Parse json
bool load_json(const string& filename, json_value& json, string& error) {
  try {
    load_json(filename, json);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// Save json
bool save_json(const string& filename, json_value& json, string& error) {
  try {
    save_json(filename, json);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// Parse json
bool parse_json(const string& text, json_value& json, string& error) {
  try {
    parse_json(text, json);
    return true;
  } catch (const json_error& exception) {
    error = exception.what();
    return false;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// Dump json
bool format_json(string& text, const json_value& json, string& error) {
  try {
    format_json(text, json);
    return true;
  } catch (const json_error& exception) {
    error = exception.what();
    return false;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// FAST CONVERSIONS FROM/TO CHARS
// -----------------------------------------------------------------------------
namespace yocto {

from_chars_result from_chars(
    const char* first, const char* last, float& value, chars_format fmt_) {
  auto fmt    = fmt_ == std::chars_format::general
                    ? fast_float::chars_format::general
                : fmt_ == std::chars_format::fixed
                    ? fast_float::chars_format::fixed
                    : fast_float::chars_format::scientific;
  auto result = fast_float::from_chars(first, last, value, fmt);
  return {result.ptr, result.ec};
}
from_chars_result from_chars(
    const char* first, const char* last, double& value, chars_format fmt_) {
  auto fmt    = fmt_ == std::chars_format::general
                    ? fast_float::chars_format::general
                : fmt_ == std::chars_format::fixed
                    ? fast_float::chars_format::fixed
                    : fast_float::chars_format::scientific;
  auto result = fast_float::from_chars(first, last, value, fmt);
  return {result.ptr, result.ec};
}

to_chars_result to_chars(
    char* first, char* last, float value, chars_format fmt) {
  if (last - first >= std::numeric_limits<float>::max_digits10)
    return {last, std::errc::value_too_large};
  auto ptr = ::nlohmann::detail::to_chars(first, last, value);
  return {ptr, std::errc()};
}
to_chars_result to_chars(
    char* first, char* last, double value, chars_format fmt) {
  if (last - first >= std::numeric_limits<double>::max_digits10)
    return {last, std::errc::value_too_large};
  auto ptr = ::nlohmann::detail::to_chars(first, last, value);
  return {ptr, std::errc()};
}

}  // namespace yocto
