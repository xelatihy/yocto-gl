//
// # Yocto/Json: Tiny collection of utilities to support JSON in Yocto/GL
//
//
// Yocto/Json is a collection of utilities used in writing other Yocto/GL
// libraries that use the JSON file format. We support converting basic Yocto/GL
// types to JSON and load/saving JSON files. Internally we use nlohmann::json as
// the Json implementation. This might change in the future. You should consider
// this library as just a helper for the rest of Yocto/GL.
//
//
// ## JSON conversion
//
// We adopt nlohmann::json (https://github.com/nlohmann/json) conventions to
// convert to and from JSON. Please see the documentation of the library. That
// implementation makes heavy use of exceptions. We prefer to use exception-free
// code and provide the functions `serialize_json_value()` for it.
//
//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#include "yocto_math.h"
#include "yocto_utils.h"

#include "ext/json.hpp"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

using nlohmann::json;

}  // namespace yocto

// -----------------------------------------------------------------------------
// FILE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a text file
inline bool load_json(const string& filename, json& js);
inline bool save_json(const string& filename, const json& js);

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON SERIALIZATION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Serialize/deserialize basic types
inline bool serialize_json_value(json& js, int& value, bool save);
inline bool serialize_json_value(json& js, bool& value, bool save);
inline bool serialize_json_value(json& js, unsigned char& value, bool save);
inline bool serialize_json_value(json& js, float& value, bool save);
inline bool serialize_json_value(json& js, double& value, bool save);
inline bool serialize_json_value(json& js, string& value, bool save);
inline bool serialize_json_value(json& js, vec2f& value, bool save);
inline bool serialize_json_value(json& js, vec3f& value, bool save);
inline bool serialize_json_value(json& js, vec4f& value, bool save);
inline bool serialize_json_value(json& js, vec2i& value, bool save);
inline bool serialize_json_value(json& js, vec3i& value, bool save);
inline bool serialize_json_value(json& js, vec4i& value, bool save);
inline bool serialize_json_value(json& js, vec4b& value, bool save);
inline bool serialize_json_value(json& js, mat2f& value, bool save);
inline bool serialize_json_value(json& js, mat3f& value, bool save);
inline bool serialize_json_value(json& js, mat4f& value, bool save);
inline bool serialize_json_value(json& js, frame2f& value, bool save);
inline bool serialize_json_value(json& js, frame3f& value, bool save);
inline bool serialize_json_value(json& js, bbox3f& value, bool save);

// Serialize/deserialize compound types
template <typename T>
inline bool serialize_json_value(json& js, vector<T>& value, bool save);

// Serialize/deserialize values as keys in a JSON object
template <typename T>
inline bool serialize_json_value(
    json& js, T& value, const char* name, const T& def, bool save);

// Check if a JSON value has a key
inline bool has_json_key(const json& js, const char* key);

// Get a value from a JSON key or a default value if any error occurs
template <typename T>
inline T get_json_value(
    const json& js, const char* key, const T& default_value);

}  // namespace yocto

// ---------------------------------------------------------------------------//
//                                                                            //
//                             IMPLEMENTATION                                 //
//                                                                            //
// ---------------------------------------------------------------------------//

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR FILE READING
// -----------------------------------------------------------------------------
namespace yocto {

// Load a JSON object
inline bool load_json(const string& filename, json& js) {
    auto text = ""s;
    if (!load_text(filename, text)) return false;
    try {
        js = json::parse(text.begin(), text.end());
    } catch (...) {
        log_io_error("could not parse json {}", filename);
        return false;
    }
    return true;
}

// Save a JSON object
inline bool save_json(const string& filename, const json& js) {
    auto str = ""s;
    try {
        str = js.dump(4);
    } catch (...) {
        log_io_error("could not dump json {}", filename);
        return false;
    }
    return save_text(filename, str);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF CONVERSION TO/FROM JSON
// -----------------------------------------------------------------------------
namespace yocto {

inline void to_json_array(json& js, const float* vals, int N) {
    js = json::array();
    for (auto i = 0; i < N; i++) js.push_back(vals[i]);
}
inline void to_json_array(json& js, const int* vals, int N) {
    js = json::array();
    for (auto i = 0; i < N; i++) js.push_back(vals[i]);
}
inline void to_json_array(json& js, const byte* vals, int N) {
    js = json::array();
    for (auto i = 0; i < N; i++) js.push_back((int)vals[i]);
}
inline void from_json_array(const json& js, float* vals, int N) {
    if (js.is_array()) throw std::runtime_error("canot convert json");
    if (js.size() != N) throw std::runtime_error("canot convert json");
    for (auto i = 0; i < N; i++) vals[i] = js.at(i).get<float>();
}
inline void from_json_array(const json& js, int* vals, int N) {
    if (js.is_array()) throw std::runtime_error("canot convert json");
    if (js.size() != N) throw std::runtime_error("canot convert json");
    for (auto i = 0; i < N; i++) vals[i] = js.at(i).get<int>();
}
inline void from_json_array(const json& js, byte* vals, int N) {
    if (js.is_array()) throw std::runtime_error("canot convert json");
    if (js.size() != N) throw std::runtime_error("canot convert json");
    for (auto i = 0; i < N; i++) vals[i] = js.at(i).get<int>();
}

inline void to_json(json& js, const vec2f& val) {
    to_json_array(js, &val.x, 2);
}
inline void from_json(const json& js, vec2f& val) {
    from_json_array(js, &val.x, 2);
}
inline void to_json(json& js, const vec3f& val) {
    to_json_array(js, &val.x, 3);
}
inline void from_json(const json& js, vec3f& val) {
    from_json_array(js, &val.x, 3);
}
inline void to_json(json& js, const vec4f& val) {
    to_json_array(js, &val.x, 4);
}
inline void from_json(const json& js, vec4f& val) {
    from_json_array(js, &val.x, 4);
}

inline void to_json(json& js, const vec2i& val) {
    to_json_array(js, &val.x, 2);
}
inline void from_json(const json& js, vec2i& val) {
    from_json_array(js, &val.x, 2);
}
inline void to_json(json& js, const vec3i& val) {
    to_json_array(js, &val.x, 3);
}
inline void from_json(const json& js, vec3i& val) {
    from_json_array(js, &val.x, 3);
}
inline void to_json(json& js, const vec4i& val) {
    to_json_array(js, &val.x, 4);
}
inline void from_json(const json& js, vec4i& val) {
    from_json_array(js, &val.x, 4);
}
inline void to_json(json& js, const vec4b& val) {
    to_json_array(js, &val.x, 4);
}
inline void from_json(const json& js, vec4b& val) {
    from_json_array(js, &val.x, 4);
}

inline void to_json(json& js, const frame3f& val) {
    to_json_array(js, &val.x.x, 12);
}
inline void from_json(const json& js, frame3f& val) {
    from_json_array(js, &val.x.x, 12);
}

inline void to_json(json& js, const mat4f& val) {
    to_json_array(js, &val.x.x, 16);
}
inline void from_json(const json& js, mat4f& val) {
    from_json_array(js, &val.x.x, 16);
}

inline void to_json(json& js, const bbox3f& val) {
    to_json_array(js, &val.min.x, 6);
}
inline void from_json(const json& js, bbox3f& val) {
    from_json_array(js, &val.min.x, 6);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR JSON SERIALIZATION
// -----------------------------------------------------------------------------
namespace yocto {

// Dumps a json value
inline bool serialize_json_value(json& js, int& value, bool save) {
    if (save) {
        js = value;
        return true;
    } else {
        if (!js.is_number_integer()) return false;
        value = js.get<int>();
        return true;
    }
}
inline bool serialize_json_value(json& js, bool& value, bool save) {
    if (save) {
        js = value;
        return true;
    } else {
        if (!js.is_boolean()) return false;
        value = js.get<bool>();
        return true;
    }
}
inline bool serialize_json_value(json& js, unsigned char& value, bool save) {
    if (save) {
        js = (int)value;
        return true;
    } else {
        if (!js.is_number_integer()) return false;
        value = (unsigned char)js.get<int>();
        return true;
    }
}
inline bool serialize_json_value(json& js, float& value, bool save) {
    if (save) {
        js = value;
        return true;
    } else {
        if (!js.is_number()) return false;
        value = js.get<float>();
        return true;
    }
}
inline bool serialize_json_value(json& js, double& value, bool save) {
    if (save) {
        js = value;
        return true;
    } else {
        if (!js.is_number()) return false;
        value = js.get<float>();
        return true;
    }
}
inline bool serialize_json_value(json& js, string& value, bool save) {
    if (save) {
        js = value;
        return true;
    } else {
        if (!js.is_string()) return false;
        value = js.get<string>();
        return true;
    }
}

template <typename T>
inline bool serialize_json_values(json& js, T* values, int num, bool save) {
    if (save) {
        js = json::array();
        for (auto i = 0; i < num; i++) {
            js.push_back({});
            if (!serialize_json_value(js.back(), values[i], save)) return false;
        }
        return true;
    } else {
        if (!js.is_array()) return false;
        if (js.size() != num) return false;
        for (auto i = 0; i < num; i++)
            if (!serialize_json_value(js.at(i), values[i], save)) return false;
        return true;
    }
}

template <typename T>
inline bool serialize_json_value(json& js, vector<T>& value, bool save) {
    if (save) {
        js = json::array();
        for (auto i = 0; i < value.size(); i++) {
            js.push_back({});
            if (!serialize_json_value(js.back(), value[i], save)) return false;
        }
        return true;
    } else {
        if (!js.is_array()) return false;
        value.resize(js.size());
        for (auto i = 0; i < value.size(); i++)
            if (!serialize_json_value(js.at(i), value[i], save)) return false;
        return true;
    }
}

inline bool serialize_json_value(json& js, vec1f& value, bool save) {
    return serialize_json_values(js, &value.x, 1, save);
}
inline bool serialize_json_value(json& js, vec2f& value, bool save) {
    return serialize_json_values(js, &value.x, 2, save);
}
inline bool serialize_json_value(json& js, vec3f& value, bool save) {
    return serialize_json_values(js, &value.x, 3, save);
}
inline bool serialize_json_value(json& js, vec4f& value, bool save) {
    return serialize_json_values(js, &value.x, 4, save);
}
inline bool serialize_json_value(json& js, vec2i& value, bool save) {
    return serialize_json_values(js, &value.x, 2, save);
}
inline bool serialize_json_value(json& js, vec3i& value, bool save) {
    return serialize_json_values(js, &value.x, 3, save);
}
inline bool serialize_json_value(json& js, vec4i& value, bool save) {
    return serialize_json_values(js, &value.x, 4, save);
}
inline bool serialize_json_value(json& js, vec4b& value, bool save) {
    return serialize_json_values(js, &value.x, 4, save);
}
inline bool serialize_json_value(json& js, mat2f& value, bool save) {
    return serialize_json_values(js, &value.x.x, 4, save);
}
inline bool serialize_json_value(json& js, mat3f& value, bool save) {
    return serialize_json_values(js, &value.x.x, 9, save);
}
inline bool serialize_json_value(json& js, mat4f& value, bool save) {
    return serialize_json_values(js, &value.x.x, 16, save);
}
inline bool serialize_json_value(json& js, frame2f& value, bool save) {
    return serialize_json_values(js, &value.x.x, 6, save);
}
inline bool serialize_json_value(json& js, frame3f& value, bool save) {
    return serialize_json_values(js, &value.x.x, 12, save);
}
inline bool serialize_json_value(json& js, bbox3f& value, bool save) {
    return serialize_json_values(js, &value.min.x, 6, save);
}

// Dumps a json value
template <typename T>
inline bool serialize_json_value(
    json& js, T& value, const char* name, const T& def, bool save) {
    if (save) {
        if (value == def) return true;
        return serialize_json_value(js[name], value, save);
    } else {
        if (!js.count(name)) return true;
        value = def;
        return serialize_json_value(js.at(name), value, save);
    }
}

// Check if a JSON value has a key
inline bool has_json_key(const json& js, const char* key) {
    return js.is_object() && js.count(key) > 0;
}

// Get a value from a JSON key or a default value if any error occurs
template <typename T>
inline T get_json_value(
    const json& js, const char* key, const T& default_value) {
    if (!js.is_object()) return default_value;
    if (js.count(key) <= 0) return default_value;
    auto value = default_value;
    if (!serialize_json_value((json&)js.at(key), value, false))
        return default_value;
    return value;
}

}  // namespace yocto

#endif
