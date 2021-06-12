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
  none,
  integer,
  unsigned_,
  number,
  boolean,
  string,
  array,
  object
};

// Json typdefs
struct json_value;
using json_array  = vector<json_value>;
using json_object = ordered_map<string, json_value>;

// Json value
struct json_value {
  json_value() : type{json_type::none} {}
  json_value(json_type type) : type{type} {}
  json_value(const json_value& other) {
    type = other.type;
    switch (type) {
      case json_type::none: break;
      case json_type::integer: integer = other.integer; break;
      case json_type::unsigned_: unsigned_ = other.unsigned_; break;
      case json_type::number: number = other.number; break;
      case json_type::boolean: boolean = other.boolean; break;
      case json_type::string: string_ = other.string_; break;
      case json_type::array: array = other.array; break;
      case json_type::object: object = other.object; break;
    }
  }
  json_value(json_value&& value) : json_value() { _swap(value); }
  json_value& operator=(json_value other) {
    _swap(other);
    return *this;
  }

  json_type get_type() const { return type; }
  void      set_type(json_type type) {
    if (this->type == type) return;
    this->type = type;
  }

  void set_null() { _set(json_type::none, integer, (int64_t)0); }
  void set_integer(int64_t value) { _set(json_type::integer, integer, value); }
  void set_unsigned(uint64_t value) {
    _set(json_type::unsigned_, unsigned_, value);
  }
  void set_number(double value) { _set(json_type::number, number, value); }
  void set_boolean(bool value) { _set(json_type::boolean, boolean, value); }
  void set_string(const char* value) {
    _set(json_type::string, string_, string{value});
  }
  void set_string(const string& value) {
    _set(json_type::string, string_, value);
  }
  void set_array() { _set(json_type::array, array, {}); }
  void set_array(size_t size) {
    _set(json_type::array, array, json_array(size));
  }
  void set_array(const json_array& value) {
    _set(json_type::array, array, value);
  }
  void set_object() { _set(json_type::object, object, {}); }
  void set_object(const json_object& value) {
    return _set(json_type::object, object, value);
  }

  int64_t&     get_integer() { return _get(json_type::integer, integer); }
  uint64_t&    get_unsigned() { return _get(json_type::unsigned_, unsigned_); }
  double&      get_number() { return _get(json_type::number, number); }
  bool&        get_boolean() { return _get(json_type::boolean, boolean); }
  string&      get_string() { return _get(json_type::string, string_); }
  json_array&  get_array() { return _get(json_type::array, array); }
  json_object& get_object() { return _get(json_type::object, object); }

  const int64_t& get_integer() const {
    return _get(json_type::integer, integer);
  }
  const uint64_t& get_unsigned() const {
    return _get(json_type::unsigned_, unsigned_);
  }
  const double& get_number() const { return _get(json_type::number, number); }
  const bool& get_boolean() const { return _get(json_type::boolean, boolean); }
  const string& get_string() const { return _get(json_type::string, string_); }
  const json_array&  get_array() const { return _get(json_type::array, array); }
  const json_object& get_object() const {
    return _get(json_type::object, object);
  }

  void set(std::nullptr_t) { set_null(); }
  void set(int32_t value) { set_integer(value); }
  void set(int64_t value) { set_integer(value); }
  void set(uint32_t value) { set_unsigned(value); }
  void set(uint64_t value) { set_unsigned(value); }
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
    } else if (get_type() == json_type::unsigned_) {
      value = (int64_t)get_unsigned();
    } else {
      throw json_error{"integer expected"};
    }
  }
  void get(uint32_t& value) const { value = (uint32_t)get<uint64_t>(); }
  void get(uint64_t& value) const {
    if (get_type() == json_type::integer) {
      value = (uint64_t)get_integer();
    } else if (get_type() == json_type::unsigned_) {
      value = (uint64_t)get_unsigned();
    } else {
      throw json_error{"integer expected"};
    }
  }
  void get(float& value) const { value = (double)get<double>(); }
  void get(double& value) const {
    if (get_type() == json_type::integer) {
      value = (double)get_integer();
    } else if (get_type() == json_type::unsigned_) {
      value = (double)get_unsigned();
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
  json_type   type      = json_type::none;
  int64_t     integer   = 0;
  uint64_t    unsigned_ = 0;
  double      number    = 0;
  bool        boolean   = false;
  string      string_   = {};
  json_array  array     = {};
  json_object object    = {};

  void _swap(json_value& other) {
    std::swap(type, other.type);
    switch (type) {
      case json_type::none: break;
      case json_type::integer: std::swap(integer, other.integer); break;
      case json_type::unsigned_: std::swap(unsigned_, other.unsigned_); break;
      case json_type::number: std::swap(number, other.number); break;
      case json_type::boolean: std::swap(boolean, other.boolean); break;
      case json_type::string: std::swap(string_, other.string_); break;
      case json_type::array: std::swap(array, other.array); break;
      case json_type::object: std::swap(object, other.object); break;
    }
  }
  template <typename T>
  const T& _get(json_type type, const T& value) const {
    if (this->type != type) throw std::invalid_argument{"bad json type"};
    return value;
  }
  template <typename T>
  T& _get(json_type type, T& value) {
    if (this->type != type) throw std::invalid_argument{"bad json type"};
    return value;
  }
  template <typename T>
  void _set(json_type type, T& value, const T& other) {
    if (this->type != type) set_type(type);
    value = other;
  }
};

// Json schema
struct json_schema {
  json_type           type            = json_type::object;
  string              title           = "";
  string              description     = "";
  json_value          default_        = {};
  json_value          min             = {};
  json_value          max             = {};
  vector<string>      enum_           = {};
  size_t              min_items       = 0;
  size_t              max_items       = std::numeric_limits<size_t>::max();
  vector<string>      required        = {};
  vector<string>      cli_positionals = {};
  string              cli_config      = "";
  vector<json_schema> items           = {};
  ordered_map<string, json_schema> properties = {};
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
  json.set_unsigned(value);
}
inline void to_json(json_value& json, uint32_t value) {
  json.set_unsigned(value);
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
  } else if (json.get_type() == json_type::unsigned_) {
    value = (int64_t)json.get_unsigned();
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
  } else if (json.get_type() == json_type::unsigned_) {
    value = (uint64_t)json.get_unsigned();
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
  } else if (json.get_type() == json_type::unsigned_) {
    value = (double)json.get_unsigned();
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
  schema.type        = to_json(value).get_type();
  schema.title       = title;
  schema.description = description;
  schema.default_    = to_json(value);
}
inline void to_schema(json_schema& schema, int64_t value, const string& title,
    const string& description, int64_t min, int64_t max) {
  schema.type        = to_json(value).get_type();
  schema.title       = title;
  schema.description = description;
  schema.default_    = to_json(value);
  schema.min         = to_json(min);
  schema.max         = to_json(max);
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
  schema.type        = to_json(value).get_type();
  schema.title       = title;
  schema.description = description;
  schema.default_    = to_json(value);
}
inline void to_schema(json_schema& schema, uint64_t value, const string& title,
    const string& description, uint64_t min, uint64_t max) {
  schema.type        = to_json(value).get_type();
  schema.title       = title;
  schema.description = description;
  schema.default_    = to_json(value);
  schema.min         = to_json(min);
  schema.max         = to_json(max);
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
  schema.type        = to_json(value).get_type();
  schema.title       = title;
  schema.description = description;
  schema.default_    = to_json(value);
}
inline void to_schema(json_schema& schema, double value, const string& title,
    const string& description, double min, double max) {
  schema.type        = to_json(value).get_type();
  schema.title       = title;
  schema.description = description;
  schema.default_    = to_json(value);
  schema.min         = to_json(min);
  schema.max         = to_json(max);
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
  schema.type        = to_json(value).get_type();
  schema.title       = title;
  schema.description = description;
  schema.default_    = to_json(value);
}
inline void to_schema(json_schema& schema, const string& value,
    const string& title, const string& description) {
  schema.type        = to_json(value).get_type();
  schema.title       = title;
  schema.description = description;
  schema.default_    = to_json(value);
}
inline void to_schema(json_schema& schema, const string& value,
    const string& title, const string& description,
    const vector<string>& enum_) {
  schema.type        = to_json(value).get_type();
  schema.title       = title;
  schema.description = description;
  schema.default_    = to_json(value);
  schema.enum_       = enum_;
}
template <typename T, size_t N>
inline void to_schema(json_schema& schema, const array<T, N>& value,
    const string& title, const string& description) {
  schema.type        = to_json(value).get_type();
  schema.title       = title;
  schema.description = description;
  schema.default_    = to_json(value);
  schema.min_items   = N;
  schema.max_items   = N;
  to_schema(schema.items.emplace_back(), T{}, "item", "");
}
template <typename T>
inline void to_schema(json_schema& schema, const vector<T>& value,
    const string& title, const string& description) {
  schema.type        = to_json(value).type;
  schema.title       = title;
  schema.description = description;
  schema.default_    = to_json(value);
  schema.min_items   = 0;
  schema.max_items   = std::numeric_limits<size_t>::max();
  to_schema(schema.items.emplace_back(), T{}, "item", "");
}

}  // namespace yocto

#endif
