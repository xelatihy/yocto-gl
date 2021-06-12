//
// # Yocto/Json: Utilities for handling Json
//
// Yocto/Json is a collection of utilities used in hadling Json.
// Yocto/Json is implemented in `yocto_json.h` and `yocto_json.cpp`, and
// depends on `json.hpp` for Json serialization.
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

#ifndef _YOCTO_JSON_H_
#define _YOCTO_JSON_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <array>
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
using std::array;
using std::string;
using std::vector;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON DATA
// -----------------------------------------------------------------------------
namespace yocto {

// Simple ordered map
template <typename Key, typename Value>
struct ordered_map {
  size_t size() const { return _data.size(); }
  bool   empty() const { return _data.empty(); }

  bool contains(const Key& key) const { return (bool)_search(key); }

  Value& operator[](const Key& key) {
    if (auto it = _search(key); it) return it->second;
    _data.emplace_back(key, Value{});
    return _data.back().second;
  }
  Value& at(const Key& key) {
    if (auto it = _search(key); it) return it->second;
    throw std::out_of_range{"missing key for " + key};
  }
  const Value& at(const Key& key) const {
    if (auto it = _search(key); it) return it->second;
    throw std::out_of_range{"missing key for " + key};
  }

  pair<Key, Value>* find(const string& key) {
    if (auto it = _search(key); it) return it;
    return end();
  }
  const pair<Key, Value>* find(const string& key) const {
    if (auto it = _search(key); it) return it;
    return end();
  }

  pair<Key, Value>*       begin() { return _data.data(); }
  pair<Key, Value>*       end() { return _data.data() + _data.size(); }
  const pair<Key, Value>* begin() const { return _data.data(); }
  const pair<Key, Value>* end() const { return _data.data() + _data.size(); }

 private:
  vector<pair<Key, Value>> _data;
  pair<Key, Value>*        _search(const string& key) {
    for (auto& item : _data)
      if (key == item.first) return &item;
    return nullptr;
  }
  const pair<Key, Value>* _search(const string& key) const {
    for (auto& item : _data)
      if (key == item.first) return &item;
    return nullptr;
  }
};

struct json_error : std::logic_error {
  json_error(const string& error) : std::logic_error(error) {}
};

// Json type
enum struct json_type {
  // clang-format off
  null, integer, uinteger, number, boolean, string, array, object
  // clang-format on
};

// Json typdefs
struct json_value;
using json_array  = vector<json_value>;
using json_object = ordered_map<string, json_value>;

// Json value
struct json_value {
  json_value() : _type{json_type::null} {}
  json_value(json_type type) : _type{type} {
    switch (_type) {
      case json_type::string: _string = new string{}; break;
      case json_type::array: _array = new json_array{}; break;
      case json_type::object: _object = new json_object{}; break;
      default: break;
    }
  }
  json_value(const json_value& other) {
    _type = other._type;
    switch (_type) {
      case json_type::null: break;
      case json_type::integer: _integer = other._integer; break;
      case json_type::uinteger: _uinteger = other._uinteger; break;
      case json_type::number: _number = other._number; break;
      case json_type::boolean: _boolean = other._boolean; break;
      case json_type::string: _string = new string{*other._string}; break;
      case json_type::array: _array = new json_array{*other._array}; break;
      case json_type::object: _object = new json_object{*other._object}; break;
    }
  }
  json_value(json_value&& value) : json_value() { _swap(value); }
  json_value& operator=(json_value other) {
    _swap(other);
    return *this;
  }
  ~json_value() {
    switch (_type) {
      case json_type::string: delete _string; break;
      case json_type::array: delete _array; break;
      case json_type::object: delete _object; break;
      default: break;
    }
  }

  json_type get_type() const { return _type; }
  void      set_type(json_type type) {
    if (_type == type) return;
    auto new_json = json_value{type};
    _swap(new_json);
  }

  bool is_null() const { return _type == json_type::null; }
  bool is_integer() const {
    return _type == json_type::integer || _type == json_type::uinteger;
  }
  bool is_number() const {
    return _type == json_type::integer || _type == json_type::uinteger ||
           _type == json_type::number;
  }
  bool is_boolean() const { return _type == json_type::boolean; }
  bool is_string() const { return _type == json_type::string; }
  bool is_array() const { return _type == json_type::array; }
  bool is_object() const { return _type == json_type::object; }

  void set_null() { set_type(json_type::null); }
  void set_integer(int64_t value) {
    set_type(json_type::integer);
    _integer = value;
  }
  void set_uinteger(uint64_t value) {
    set_type(json_type::uinteger);
    _uinteger = value;
  }
  void set_number(double value) {
    set_type(json_type::number);
    _number = value;
  }
  void set_boolean(bool value) {
    set_type(json_type::boolean);
    _boolean = value;
  }
  void set_string(const char* value) {
    set_type(json_type::string);
    *_string = value;
  }
  void set_string(const string& value) {
    set_type(json_type::string);
    *_string = value;
  }
  void set_array() {
    set_type(json_type::array);
    *_array = {};
  }
  void set_array(size_t size) {
    set_type(json_type::array);
    _array->resize(size);
  }
  void set_array(const json_array& value) {
    set_type(json_type::array);
    *_array = value;
  }
  void set_object() {
    set_type(json_type::object);
    *_object = {};
  }
  void set_object(const json_object& value) {
    set_type(json_type::object);
    *_object = value;
  }

  int64_t& get_integer() {
    if (_type != json_type::integer) throw json_error{"integer expected"};
    return _integer;
  }
  uint64_t& get_uinteger() {
    if (_type != json_type::uinteger) throw json_error{"integer expected"};
    return _uinteger;
  }
  double& get_number() {
    if (_type != json_type::number) throw json_error{"number expected"};
    return _number;
  }
  bool& get_boolean() {
    if (_type != json_type::boolean) throw json_error{"boolean expected"};
    return _boolean;
  }
  string& get_string() {
    if (_type != json_type::string) throw json_error{"string expected"};
    return *_string;
  }
  json_array& get_array() {
    if (_type != json_type::array) throw json_error{"array expected"};
    return *_array;
  }
  json_object& get_object() {
    if (_type != json_type::object) throw json_error{"object expected"};
    return *_object;
  }

  const int64_t& get_integer() const {
    if (_type != json_type::integer) throw json_error{"integer expected"};
    return _integer;
  }
  const uint64_t& get_uinteger() const {
    if (_type != json_type::uinteger) throw json_error{"integer expected"};
    return _uinteger;
  }
  const double& get_number() const {
    if (_type != json_type::number) throw json_error{"number expected"};
    return _number;
  }
  const bool& get_boolean() const {
    if (_type != json_type::boolean) throw json_error{"boolean expected"};
    return _boolean;
  }
  const string& get_string() const {
    if (_type != json_type::string) throw json_error{"string expected"};
    return *_string;
  }
  const json_array& get_array() const {
    if (_type != json_type::array) throw json_error{"array expected"};
    return *_array;
  }
  const json_object& get_object() const {
    if (_type != json_type::object) throw json_error{"object expected"};
    return *_object;
  }

  bool empty() const {
    if (is_array()) return _array->empty();
    if (is_object()) return _object->empty();
    throw json_error{"array or object expected"};
  }
  size_t size() const {
    if (is_array()) return _array->size();
    if (is_object()) return _object->size();
    throw json_error{"array or object expected"};
  }

  json_value& operator[](size_t idx) {
    if (is_array()) return _array->at(idx);
    throw json_error{"array expected"};
  }
  const json_value& operator[](size_t idx) const {
    if (is_array()) return _array->at(idx);
    throw json_error{"array expected"};
  }
  json_value& operator[](const string& key) {
    if (is_object()) return _object->operator[](key);
    throw json_error{"object expected"};
  }
  json_value& at(size_t idx) {
    if (is_array()) return _array->at(idx);
    throw json_error{"array expected"};
  }
  const json_value& at(size_t idx) const {
    if (is_array()) return _array->at(idx);
    throw json_error{"array expected"};
  }
  json_value& at(const string& key) {
    if (is_object()) return _object->at(key);
    throw json_error{"object expected"};
  }
  const json_value& at(const string& key) const {
    if (is_object()) return _object->at(key);
    throw json_error{"object expected"};
  }

  template <typename... Args>
  json_value& emplace_back(Args&&... args) {
    if (is_array()) return _array->emplace_back(std::forward(args)...);
    throw json_error{"array expected"};
  }

  void set(std::nullptr_t) { set_null(); }
  void set(int32_t value) { set_integer(value); }
  void set(int64_t value) { set_integer(value); }
  void set(uint32_t value) { set_uinteger(value); }
  void set(uint64_t value) { set_uinteger(value); }
  void set(float value) { set_number(value); }
  void set(double value) { set_number(value); }
  void set(bool value) { set_boolean(value); }
  void set(const string& value) { set_string(value); }
  void set(const char* value) { set_string(value); }
  template <typename T, size_t N>
  void set(const std::array<T, N>& value) {
    set_array(value.size());
    for (auto idx = (size_t)0; idx < value.size(); idx++)
      get_array().at(idx).set(value.at(idx));
  }
  template <typename T>
  void set(const vector<T>& value) {
    set_array(value.size());
    for (auto idx = (size_t)0; idx < value.size(); idx++)
      get_array().at(idx).set(value.at(idx));
  }

  template <typename T>
  T get() const {
    auto value = T{};
    get(value);
    return value;
  }

  void get(int32_t& value) const { value = (int32_t)get<int64_t>(); }
  void get(int64_t& value) const {
    if (get_type() == json_type::integer) {
      value = (int64_t)get_integer();
    } else if (get_type() == json_type::uinteger) {
      value = (int64_t)get_uinteger();
    } else {
      throw json_error{"integer expected"};
    }
  }
  void get(uint32_t& value) const { value = (uint32_t)get<uint64_t>(); }
  void get(uint64_t& value) const {
    if (get_type() == json_type::integer) {
      value = (uint64_t)get_integer();
    } else if (get_type() == json_type::uinteger) {
      value = (uint64_t)get_uinteger();
    } else {
      throw json_error{"integer expected"};
    }
  }
  void get(float& value) const { value = (double)get<double>(); }
  void get(double& value) const {
    if (get_type() == json_type::integer) {
      value = (double)get_integer();
    } else if (get_type() == json_type::uinteger) {
      value = (double)get_uinteger();
    } else if (get_type() == json_type::number) {
      value = (double)get_number();
    } else {
      throw json_error{"number expected"};
    }
  }
  void get(bool& value) const {
    if (get_type() == json_type::boolean) {
      value = get_boolean();
    } else {
      throw json_error{"boolean expected"};
    }
  }
  void get(string& value) const {
    if (get_type() == json_type::string) {
      value = get_string();
    } else {
      throw json_error{"string expected"};
    }
  }
  template <typename T, size_t N>
  void get(std::array<T, N>& value) const {
    if (get_type() == json_type::array) {
      if (get_array().size() != N)
        throw json_error{"array of size " + std::to_string(N) + " expected"};
      for (auto idx = (size_t)0; idx < value.size(); idx++)
        get_array().at(idx).get(value.at(idx));
    } else {
      throw json_error{"array expected"};
    }
  }
  template <typename T>
  void get(vector<T>& value) const {
    if (get_type() == json_type::array) {
      value.resize(get_array().size());
      for (auto idx = (size_t)0; idx < value.size(); idx++)
        get_array().at(idx).get(value.at(idx));
    } else {
      throw json_error{"array expected"};
    }
  }

 private:
  json_type _type = json_type::null;
  union {
    int64_t      _integer = 0;
    uint64_t     _uinteger;
    double       _number;
    bool         _boolean;
    string*      _string;
    json_array*  _array;
    json_object* _object;
  };

  void _swap(json_value& other) {
    std::swap(_type, other._type);
    std::swap(_integer, other._integer);
  }
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON I/O
// -----------------------------------------------------------------------------
namespace yocto {

// Load json
bool load_json(const string& filename, json_value& json, string& error);
// Save json
bool save_json(const string& filename, json_value& json, string& error);

// Parse json
bool parse_json(const string& text, json_value& json, string& error);
// Dump json
bool dump_json(string& text, const json_value& json, string& error);

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON SCHEMA
// -----------------------------------------------------------------------------
namespace yocto {

// Json schema
struct json_schema {
  void      set_type(json_type type) { _type = type; }
  json_type get_type() const { return _type; }

  bool is_null() const { return _type == json_type::null; }
  bool is_integer() const {
    return _type == json_type::integer || _type == json_type::uinteger;
  }
  bool is_number() const {
    return _type == json_type::integer || _type == json_type::uinteger ||
           _type == json_type::number;
  }
  bool is_boolean() const { return _type == json_type::boolean; }
  bool is_string() const { return _type == json_type::string; }
  bool is_array() const { return _type == json_type::array; }
  bool is_object() const { return _type == json_type::object; }

  void set_integer() { set_type(json_type::integer); }
  void set_uinteger() { set_type(json_type::uinteger); }
  void set_number() { set_type(json_type::number); }
  void set_boolean() { set_type(json_type::boolean); }
  void set_string() { set_type(json_type::string); }
  void set_array() { set_type(json_type::array); }
  void set_object() { set_type(json_type::object); }

  void          set_title(const string& title) { _title = title; }
  const string& get_title() const { return _title; }
  void          set_description(const string& description) {
    _description = description;
  }
  const string& get_description() const { return _description; }

  template <typename T>
  void set_default(const T& default_) {
    _default.set(default_);
  }
  const json_value& get_default() const { return _default; }
  json_value&       get_default() { return _default; }

  template <typename T>
  void set_min(const T& min) {
    _min.set(min);
  }
  const json_value& get_min() const { return _min; }
  json_value&       get_min() { return _min; }
  template <typename T>
  void set_max(const T& max) {
    _max.set(max);
  }
  const json_value& get_max() const { return _max; }
  json_value&       get_max() { return _max; }

  void                  set_enum(const vector<string>& enum_) { _enum = enum_; }
  const vector<string>& get_enum() const { return _enum; }

  void   set_minitems(size_t minitems) { _minitems = minitems; }
  size_t get_minitems() const { return _minitems; }
  void   set_maxitems(size_t maxitems) { _maxitems = maxitems; }
  size_t get_maxitems() const { return _maxitems; }

  void set_required(const vector<string>& required) { _required = required; }
  const vector<string>& get_required() const { return _required; }
  void add_required(const string& required) { _required.push_back(required); }

  void set_clipositional(const vector<string>& clipositional) {
    _clipositional = clipositional;
  }
  const vector<string>& get_clipositional() const { return _clipositional; }
  void                  add_clipositional(const string& clipositional) {
    _clipositional.push_back(clipositional);
  }

  void set_cliconfig(const string& cliconfig) { _cliconfig = cliconfig; }
  const string& get_cliconfig() const { return _cliconfig; }

  void set_items(const vector<json_schema>& items) { _items = items; }
  const vector<json_schema>& get_items() const { return _items; }

  json_schema&       get_item(size_t idx) { return _items.at(idx); }
  const json_schema& get_item(size_t idx) const { return _items.at(idx); }
  json_schema&       add_item() { return _items.emplace_back(); }
  void add_item(const json_schema& item) { _items.push_back(item); }

  void set_properties(const ordered_map<string, json_schema>& properties) {
    _properties = properties;
  }
  const ordered_map<string, json_schema>& get_properties() const {
    return _properties;
  }

  bool has_property(const string& name) const {
    return _properties.contains(name);
  }
  void set_property(const string& name, const json_schema& item) {
    _properties[name] = item;
  }
  json_schema& get_property(const string& name) { return _properties.at(name); }
  const json_schema& get_property(const string& name) const {
    return _properties.at(name);
  }
  json_schema& add_property(const string& name) { return _properties[name]; }
  void         add_property(const string& name, const json_schema& item) {
    _properties[name] = item;
  }

 private:
  json_type           _type          = json_type::object;
  string              _title         = "";
  string              _description   = "";
  json_value          _default       = {};
  json_value          _min           = {};
  json_value          _max           = {};
  vector<string>      _enum          = {};
  size_t              _minitems      = 0;
  size_t              _maxitems      = std::numeric_limits<size_t>::max();
  vector<string>      _required      = {};
  vector<string>      _clipositional = {};
  string              _cliconfig     = "";
  vector<json_schema> _items         = {};
  ordered_map<string, json_schema> _properties = {};
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

inline size_t json_empty(const json_value& json) {
  if (json.get_type() == json_type::array) {
    return json.get_array().empty();
  } else if (json.get_type() == json_type::object) {
    return json.get_object().empty();
  } else {
    throw json_error{"array or object expected"};
  }
}
inline size_t json_size(const json_value& json) {
  if (json.get_type() == json_type::array) {
    return json.get_array().size();
  } else if (json.get_type() == json_type::object) {
    return json.get_object().size();
  } else {
    throw json_error{"array or object expected"};
  }
}
inline const json_value& json_at(const json_value& json, size_t idx) {
  if (json.get_type() == json_type::array) {
    if (idx >= json.get_array().size()) throw json_error{"out of bounds"};
    return json.get_array().at(idx);
  } else {
    throw json_error{"array expected"};
  }
}
inline json_value& json_at(json_value& json, size_t idx) {
  if (json.get_type() == json_type::array) {
    if (idx >= json.get_array().size()) throw json_error{"out of bounds"};
    return json.get_array().at(idx);
  } else {
    throw json_error{"array expected"};
  }
}
inline const json_value& json_at(const json_value& json, const string& key) {
  if (json.get_type() == json_type::array) {
    auto it = json.get_object().find(key);
    if (it == json.get_object().end()) throw json_error{"missing key " + key};
    return it->second;
  } else {
    throw json_error{"object expected"};
  }
}
inline json_value& json_at(json_value& json, const string& key) {
  if (json.get_type() == json_type::array) {
    auto it = json.get_object().find(key);
    if (it == json.get_object().end()) throw json_error{"missing key " + key};
    return it->second;
  } else {
    throw json_error{"object expected"};
  }
}
inline json_value& json_append(json_value& json) {
  if (json.get_type() == json_type::array) {
    return json.get_array().emplace_back();
  } else {
    throw json_error{"array expected"};
  }
}
inline json_value& json_insert(json_value& json, const string& key) {
  if (json.get_type() == json_type::object) {
    return json.get_object()[key];
  } else {
    throw json_error{"object expected"};
  }
}

template <typename T>
inline json_value to_json(const T& value) {
  auto json = json_value{};
  to_json(json, value);
  return json;
}

inline void to_json(json_value& json, int64_t value) {
  json.set_integer(value);
}
inline void to_json(json_value& json, int32_t value) {
  json.set_integer(value);
}
inline void to_json(json_value& json, uint64_t value) {
  json.set_uinteger(value);
}
inline void to_json(json_value& json, uint32_t value) {
  json.set_uinteger(value);
}
inline void to_json(json_value& json, double value) { json.set_number(value); }
inline void to_json(json_value& json, float value) { json.set_number(value); }
inline void to_json(json_value& json, bool value) { json.set_boolean(value); }
inline void to_json(json_value& json, const char* value) {
  json.set_string(value);
}
inline void to_json(json_value& json, const string& value) {
  json.set_string(value);
}
inline void to_json(json_value& json, const json_array& value) {
  json.set_array(value);
}
inline void to_json(json_value& json, const json_object& value) {
  json.set_object(value);
}
template <typename T, size_t N>
inline void to_json(json_value& json, const array<T, N>& value) {
  json.set_array(value.size());
  for (auto idx = 0; idx < value.size(); idx++)
    to_json(json.get_array().at(idx), value.at(idx));
}
template <typename T>
inline void to_json(json_value& json, const vector<T>& value) {
  json.set_array(value.size());
  for (auto idx = 0; idx < value.size(); idx++)
    to_json(json.get_array().at(idx), value.at(idx));
}

inline void to_json_array(json_value& json) { json.set_array(); }
inline void to_json_array(json_value& json, size_t size) {
  json.set_array(size);
}
inline json_value& to_json_append(json_value& json) {
  if (json.get_type() != json_type::array) json.set_array();
  return json.get_array().emplace_back();
}
inline void        to_json_object(json_value& json) { json.set_object(); }
inline json_value& to_json_insert(json_value& json, const string& key) {
  if (json.get_type() != json_type::object) json.set_object();
  return json.get_object()[key];
}

template <typename T>
inline T from_json(const json_value& json) {
  auto value = T{};
  from_json(json, value);
  return value;
}

inline void from_json(const json_value& json, int64_t& value) {
  if (json.get_type() == json_type::integer) {
    value = (int64_t)json.get_integer();
  } else if (json.get_type() == json_type::uinteger) {
    value = (int64_t)json.get_uinteger();
  } else {
    throw json_error{"integer expected"};
  }
}
inline void from_json(const json_value& json, int32_t& value) {
  value = (int32_t)from_json<int64_t>(json);
}
inline void from_json(const json_value& json, uint64_t& value) {
  if (json.get_type() == json_type::integer) {
    value = (uint64_t)json.get_integer();
  } else if (json.get_type() == json_type::uinteger) {
    value = (uint64_t)json.get_uinteger();
  } else {
    throw json_error{"integer expected"};
  }
}
inline void from_json(const json_value& json, uint32_t& value) {
  value = (uint32_t)from_json<uint64_t>(json);
}
inline void from_json(const json_value& json, double& value) {
  if (json.get_type() == json_type::integer) {
    value = (double)json.get_integer();
  } else if (json.get_type() == json_type::uinteger) {
    value = (double)json.get_uinteger();
  } else if (json.get_type() == json_type::number) {
    value = (double)json.get_number();
  } else {
    throw json_error{"number expected"};
  }
}
inline void from_json(const json_value& json, float& value) {
  value = (double)from_json<double>(json);
}
inline void from_json(const json_value& json, bool& value) {
  if (json.get_type() == json_type::boolean) {
    value = json.get_boolean();
  } else {
    throw json_error{"boolean expected"};
  }
}
inline void from_json(const json_value& json, string& value) {
  if (json.get_type() == json_type::string) {
    value = json.get_string();
  } else {
    throw json_error{"string expected"};
  }
}
template <typename T, size_t N>
inline void from_json(const json_value& json, array<T, N>& value) {
  if (json.get_type() == json_type::array) {
    if (json.get_array().size() != N)
      throw json_error{"array of size " + std::to_string(N) + " expected"};
    for (auto idx = (size_t)0; idx < value.size(); idx++)
      from_json(json.get_array().at(idx), value.at(idx));
  } else {
    throw json_error{"array expected"};
  }
}
template <typename T>
inline void from_json(const json_value& json, vector<T>& value) {
  if (json.get_type() == json_type::array) {
    value.resize(json.get_array().size());
    for (auto idx = (size_t)0; idx < value.size(); idx++)
      from_json(json.get_array().at(idx), value.at(idx));
  } else {
    throw json_error{"array expected"};
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON SCHEMA
// -----------------------------------------------------------------------------
namespace yocto {

inline void to_schema(json_schema& schema, int64_t value, const string& title,
    const string& description) {
  schema.set_integer();
  schema.set_title(title);
  schema.set_description(description);
  schema.set_default(value);
}
inline void to_schema(json_schema& schema, int64_t value, const string& title,
    const string& description, int64_t min, int64_t max) {
  schema.set_integer();
  schema.set_title(title);
  schema.set_description(description);
  schema.set_default(value);
  schema.set_min(min);
  schema.set_max(max);
}
inline void to_schema(json_schema& schema, int32_t value, const string& title,
    const string& description) {
  to_schema(schema, (int64_t)value, title, description);
}
inline void to_schema(json_schema& schema, int32_t value, const string& title,
    const string& description, int32_t min, int32_t max) {
  to_schema(
      schema, (int64_t)value, title, description, (int64_t)min, (int64_t)max);
}
inline void to_schema(json_schema& schema, uint64_t value, const string& title,
    const string& description) {
  schema.set_uinteger();
  schema.set_title(title);
  schema.set_description(description);
  schema.set_default(value);
}
inline void to_schema(json_schema& schema, uint64_t value, const string& title,
    const string& description, uint64_t min, uint64_t max) {
  schema.set_uinteger();
  schema.set_title(title);
  schema.set_description(description);
  schema.set_default(value);
  schema.set_min(min);
  schema.set_max(max);
}
inline void to_schema(json_schema& schema, uint32_t value, const string& title,
    const string& description) {
  to_schema(schema, (uint64_t)value, title, description);
}
inline void to_schema(json_schema& schema, uint32_t value, const string& title,
    const string& description, uint32_t min, uint32_t max) {
  to_schema(schema, (uint64_t)value, title, description, (uint64_t)min,
      (uint64_t)max);
}
inline void to_schema(json_schema& schema, double value, const string& title,
    const string& description) {
  schema.set_number();
  schema.set_title(title);
  schema.set_description(description);
  schema.set_default(value);
}
inline void to_schema(json_schema& schema, double value, const string& title,
    const string& description, double min, double max) {
  schema.set_number();
  schema.set_title(title);
  schema.set_description(description);
  schema.set_default(value);
  schema.set_min(min);
  schema.set_max(max);
}
inline void to_schema(json_schema& schema, float value, const string& title,
    const string& description) {
  to_schema(schema, (double)value, title, description);
}
inline void to_schema(json_schema& schema, float value, const string& title,
    const string& description, float min, float max) {
  to_schema(
      schema, (double)value, title, description, (double)min, (double)max);
}
inline void to_schema(json_schema& schema, bool value, const string& title,
    const string& description) {
  schema.set_boolean();
  schema.set_title(title);
  schema.set_description(description);
  schema.set_default(value);
}
inline void to_schema(json_schema& schema, const string& value,
    const string& title, const string& description) {
  schema.set_string();
  schema.set_title(title);
  schema.set_description(description);
  schema.set_default(value);
}
inline void to_schema(json_schema& schema, const string& value,
    const string& title, const string& description,
    const vector<string>& enum_) {
  schema.set_string();
  schema.set_title(title);
  schema.set_description(description);
  schema.set_default(value);
  schema.set_enum(enum_);
}
template <typename T, size_t N>
inline void to_schema(json_schema& schema, const array<T, N>& value,
    const string& title, const string& description) {
  schema.set_array();
  schema.set_title(title);
  schema.set_description(description);
  schema.set_default(value);
  schema.set_minitems(N);
  schema.set_maxitems(N);
  to_schema(schema.add_item(), T{}, "item", "");
}
template <typename T>
inline void to_schema(json_schema& schema, const vector<T>& value,
    const string& title, const string& description) {
  schema.set_array();
  schema.set_title(title);
  schema.set_description(description);
  schema.set_default(value);
  schema.set_minitems(0);
  schema.set_maxitems(std::numeric_limits<size_t>::max());
  to_schema(schema.add_item(), T{}, "item", "");
}

}  // namespace yocto

#endif
