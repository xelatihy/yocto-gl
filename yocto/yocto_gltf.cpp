//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
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
// LICENSE OF INCLUDED CODE
//
//
// base64.cpp and base64.h
//
// Copyright (C) 2004-2008 René Nyffenegger
//
// This source code is provided 'as-is', without any express or implied
// warranty. In no event will the author be held liable for any damages
// arising from the use of this software.
//
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
//
// 1. The origin of this source code must not be misrepresented; you must not
// claim that you wrote the original source code. If you use this source code
// in a product, an acknowledgment in the product documentation would be
// appreciated but is not required.
//
// 2. Altered source versions must be plainly marked as such, and must not be
// misrepresented as being the original source code.
//
// 3. This notice may not be removed or altered from any source distribution.
//
// René Nyffenegger rene.nyffenegger@adp-gmbh.ch
//

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR YOCTO_GLTF
// -----------------------------------------------------------------------------

//
// TODO: avoid streams
// TODO: unflatten: error for pos.size()
//

#include "yocto_gltf.h"

#include <cstdio>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <vector>

#ifndef YGLTF_NO_IMAGE
#include "yocto_img.h"
#endif

namespace ygltf {

// #codegen begin func ---------------------------------------------------------

// Parse error
struct _parse_stack {
    std::vector<std::string> path;
    std::string pathname() {
        auto p = std::string();
        for (auto n : path) p += '/' + n;
        return p;
    }
};

// Parse support function.
static inline bool _parse(std::string& val, const json& js, _parse_stack& err) {
    if (!js.is_string()) return false;
    val = js;
    return true;
}

// Parse support function.
static inline bool _parse(int& val, const json& js, _parse_stack& err) {
    if (!js.is_number_integer() && !js.is_number_unsigned()) return false;
    val = js;
    return true;
}

// Parse support function.
static inline bool _parse(bool& val, const json& js, _parse_stack& err) {
    if (!js.is_boolean()) return false;
    val = js;
    return true;
}

// Parse support function.
static inline bool _parse(float& val, const json& js, _parse_stack& err) {
    if (!js.is_number()) return false;
    val = js;
    return true;
}

// Parse support function.
static inline bool _parse(double& val, const json& js, _parse_stack& err) {
    if (!js.is_number()) return false;
    val = js;
    return true;
}

// Parse support function.
static inline bool _parse(json& val, const json& js, _parse_stack& err) {
    val = js;
    return true;
}

// Parse support function.
template <typename T>
static inline bool _parse(
    std::vector<T>& vals, const json& js, _parse_stack& err) {
    if (!js.is_array()) return false;
    vals.resize(js.size());
    for (auto i = 0; i < js.size(); i++) {
        // this is contrived to support for vector<bool>
        auto v = T();
        _parse(v, js[i], err);
        vals[i] = v;
    }
    return true;
}

// Parse support function.
template <typename T, std::size_t N>
static inline bool _parse(
    std::array<T, N>& vals, const json& js, _parse_stack& err) {
    if (!js.is_array()) return false;
    if (N != js.size()) return false;
    for (auto i = 0; i < N; i++) { _parse(vals[i], js[i], err); }
    return true;
}

// Parse support function.
template <typename T>
static inline bool _parse(
    std::map<std::string, T>& vals, const json& js, _parse_stack& err) {
    if (!js.is_object()) return false;
    for (auto kv = js.begin(); kv != js.end(); ++kv) {
        _parse(vals[kv.key()], kv.value(), err);
    }
    return true;
}

// Parse support function.
template <typename T>
static inline bool _parse(optional<T>& val, const json& js, _parse_stack& err) {
    val.valid = true;
    return _parse(val.value, js, err);
}

// Parse support function.
template <typename T>
static inline bool _parse_attr(T& val, const char* name, bool required,
    const json& js, _parse_stack& err) {
    auto exists = js.find(name) != js.end();
    err.path.push_back(name);
    if (required && !exists) return false;
    if (exists) _parse(val, js[name], err);
    err.path.pop_back();
    return true;
}

// Parse support function.
static inline bool _parse_begin_obj(const json& js, _parse_stack& err) {
    if (!js.is_object()) return false;
    return true;
}

// Parse support function.
static inline bool _parse_end_obj(const json& js, _parse_stack& err) {
    return true;
}

// Dump support function.
static inline bool _dump(const std::string& val, json& js, _parse_stack& err) {
    js = val;
    return true;
}

// Parse support function.
static inline bool _dump(const int& val, json& js, _parse_stack& err) {
    js = val;
    return true;
}

// Parse support function.
static inline bool _dump(const bool& val, json& js, _parse_stack& err) {
    js = val;
    return true;
}

// Parse support function.
static inline bool _dump(const float& val, json& js, _parse_stack& err) {
    js = val;
    return true;
}

// Parse support function.
static inline bool _dump(const double& val, json& js, _parse_stack& err) {
    js = val;
    return true;
}

// Parse support function.
static inline bool _dump(const json& val, json& js, _parse_stack& err) {
    js = val;
    return true;
}

// Parse support function.
template <typename T>
static inline bool _dump(
    const std::vector<T>& vals, json& js, _parse_stack& err) {
    js = json::array();
    for (auto i = 0; i < vals.size(); i++) { _dump(vals[i], js[i], err); }
    return true;
}

// Parse support function.
template <typename T, std::size_t N>
static inline bool _dump(
    const std::array<T, N>& vals, json& js, _parse_stack& err) {
    js = json::array();
    for (auto i = 0; i < N; i++) { _dump(vals[i], js[i], err); }
    return true;
}

// Parse support function.
template <typename T>
static inline bool _dump(
    const std::map<std::string, T>& vals, json& js, _parse_stack& err) {
    js = json::object();
    for (auto&& kv : vals) { _dump(kv.second, js[kv.first], err); }
    return true;
}

// Parse support function.
template <typename T>
static inline bool _dump(const optional<T>& val, json& js, _parse_stack& err) {
    if (!val) return true;
    _dump(val.value, js, err);
    return true;
}

// Parse support function.
template <typename T>
static inline bool _dump_attr(const T& val, const char* name, const T& defval,
    bool required, json& js, _parse_stack& err) {
    err.path.push_back(name);
    if (required || !(val == defval)) { _dump(val, js[name], err); }
    err.path.pop_back();
    return true;
}

// Parse support function.
static inline bool _dump_begin_obj(json& js, _parse_stack& err) {
    if (!js.is_object()) js = json::object();
    return true;
}

// Parse support function.
static inline bool _dump_end_obj(json& js, _parse_stack& err) { return true; }

// Equality check for defaults (might go away in the future)
static inline bool operator==(
    const glTFProperty_t& a, const glTFProperty_t& b) {
    if (!(a.extensions == b.extensions)) return false;
    if (!(a.extras == b.extras)) return false;
    return true;
}

// Parses a glTFProperty_t object
static inline bool _parse(
    glTFProperty_t& val, const json& js, _parse_stack& err) {
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.extensions, "extensions", false, js, err))
        return false;
    if (!_parse_attr(val.extras, "extras", false, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a glTFProperty_t object to JSON
static inline bool _dump(
    const glTFProperty_t& val, json& js, _parse_stack& err) {
    static const auto defval = glTFProperty_t();
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(
            val.extensions, "extensions", defval.extensions, false, js, err))
        return false;
    if (!_dump_attr(val.extras, "extras", defval.extras, false, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(
    const glTFChildOfRootProperty_t& a, const glTFChildOfRootProperty_t& b) {
    if (!((glTFProperty_t&)a == (glTFProperty_t&)b)) return false;
    if (!(a.name == b.name)) return false;
    return true;
}

// Parses a glTFChildOfRootProperty_t object
static inline bool _parse(
    glTFChildOfRootProperty_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.name, "name", false, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a glTFChildOfRootProperty_t object to JSON
static inline bool _dump(
    const glTFChildOfRootProperty_t& val, json& js, _parse_stack& err) {
    static const auto defval = glTFChildOfRootProperty_t();
    if (!_dump((glTFProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(val.name, "name", defval.name, false, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(
    const accessor_sparse_indices_t& a, const accessor_sparse_indices_t& b) {
    if (!((glTFProperty_t&)a == (glTFProperty_t&)b)) return false;
    if (!(a.bufferView == b.bufferView)) return false;
    if (!(a.byteOffset == b.byteOffset)) return false;
    if (!(a.componentType == b.componentType)) return false;
    return true;
}

// Parse a componentType_t enum
static inline bool _parse(accessor_sparse_indices_t::componentType_t& val,
    const json& js, _parse_stack& err) {
    static std::map<int, accessor_sparse_indices_t::componentType_t> table = {
        {5121, accessor_sparse_indices_t::componentType_t::unsigned_byte_t},
        {5123, accessor_sparse_indices_t::componentType_t::unsigned_short_t},
        {5125, accessor_sparse_indices_t::componentType_t::unsigned_int_t},
    };
    auto v = 0;
    if (!_parse(v, js, err)) return false;
    if (table.find(v) == table.end()) return false;
    val = table[v];
    return true;
}

// Parses a accessor_sparse_indices_t object
static inline bool _parse(
    accessor_sparse_indices_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.bufferView, "bufferView", true, js, err)) return false;
    if (!_parse_attr(val.byteOffset, "byteOffset", false, js, err))
        return false;
    if (!_parse_attr(val.componentType, "componentType", true, js, err))
        return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a componentType_t enum to JSON
static inline bool _dump(const accessor_sparse_indices_t::componentType_t& val,
    json& js, _parse_stack& err) {
    static std::map<accessor_sparse_indices_t::componentType_t, int> table = {
        {accessor_sparse_indices_t::componentType_t::unsigned_byte_t, 5121},
        {accessor_sparse_indices_t::componentType_t::unsigned_short_t, 5123},
        {accessor_sparse_indices_t::componentType_t::unsigned_int_t, 5125},
    };
    auto v = table[val];
    if (!_dump(v, js, err)) return false;
    return true;
}

// Converts a accessor_sparse_indices_t object to JSON
static inline bool _dump(
    const accessor_sparse_indices_t& val, json& js, _parse_stack& err) {
    static const auto defval = accessor_sparse_indices_t();
    if (!_dump((glTFProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(
            val.bufferView, "bufferView", defval.bufferView, true, js, err))
        return false;
    if (!_dump_attr(
            val.byteOffset, "byteOffset", defval.byteOffset, false, js, err))
        return false;
    if (!_dump_attr(val.componentType, "componentType", defval.componentType,
            true, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(
    const accessor_sparse_values_t& a, const accessor_sparse_values_t& b) {
    if (!((glTFProperty_t&)a == (glTFProperty_t&)b)) return false;
    if (!(a.bufferView == b.bufferView)) return false;
    if (!(a.byteOffset == b.byteOffset)) return false;
    return true;
}

// Parses a accessor_sparse_values_t object
static inline bool _parse(
    accessor_sparse_values_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.bufferView, "bufferView", true, js, err)) return false;
    if (!_parse_attr(val.byteOffset, "byteOffset", false, js, err))
        return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a accessor_sparse_values_t object to JSON
static inline bool _dump(
    const accessor_sparse_values_t& val, json& js, _parse_stack& err) {
    static const auto defval = accessor_sparse_values_t();
    if (!_dump((glTFProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(
            val.bufferView, "bufferView", defval.bufferView, true, js, err))
        return false;
    if (!_dump_attr(
            val.byteOffset, "byteOffset", defval.byteOffset, false, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(
    const accessor_sparse_t& a, const accessor_sparse_t& b) {
    if (!((glTFProperty_t&)a == (glTFProperty_t&)b)) return false;
    if (!(a.count == b.count)) return false;
    if (!(a.indices == b.indices)) return false;
    if (!(a.values == b.values)) return false;
    return true;
}

// Parses a accessor_sparse_t object
static inline bool _parse(
    accessor_sparse_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.count, "count", true, js, err)) return false;
    if (!_parse_attr(val.indices, "indices", true, js, err)) return false;
    if (!_parse_attr(val.values, "values", true, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a accessor_sparse_t object to JSON
static inline bool _dump(
    const accessor_sparse_t& val, json& js, _parse_stack& err) {
    static const auto defval = accessor_sparse_t();
    if (!_dump((glTFProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(val.count, "count", defval.count, true, js, err))
        return false;
    if (!_dump_attr(val.indices, "indices", defval.indices, true, js, err))
        return false;
    if (!_dump_attr(val.values, "values", defval.values, true, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(const accessor_t& a, const accessor_t& b) {
    if (!((glTFChildOfRootProperty_t&)a == (glTFChildOfRootProperty_t&)b))
        return false;
    if (!(a.bufferView == b.bufferView)) return false;
    if (!(a.byteOffset == b.byteOffset)) return false;
    if (!(a.componentType == b.componentType)) return false;
    if (!(a.count == b.count)) return false;
    if (!(a.max == b.max)) return false;
    if (!(a.min == b.min)) return false;
    if (!(a.normalized == b.normalized)) return false;
    if (!(a.sparse == b.sparse)) return false;
    if (!(a.type == b.type)) return false;
    return true;
}

// Parse a componentType_t enum
static inline bool _parse(
    accessor_t::componentType_t& val, const json& js, _parse_stack& err) {
    static std::map<int, accessor_t::componentType_t> table = {
        {5120, accessor_t::componentType_t::byte_t},
        {5121, accessor_t::componentType_t::unsigned_byte_t},
        {5122, accessor_t::componentType_t::short_t},
        {5123, accessor_t::componentType_t::unsigned_short_t},
        {5125, accessor_t::componentType_t::unsigned_int_t},
        {5126, accessor_t::componentType_t::float_t},
    };
    auto v = 0;
    if (!_parse(v, js, err)) return false;
    if (table.find(v) == table.end()) return false;
    val = table[v];
    return true;
}

// Parse a type_t enum
static inline bool _parse(
    accessor_t::type_t& val, const json& js, _parse_stack& err) {
    static std::map<std::string, accessor_t::type_t> table = {
        {"SCALAR", accessor_t::type_t::scalar_t},
        {"VEC2", accessor_t::type_t::vec2_t},
        {"VEC3", accessor_t::type_t::vec3_t},
        {"VEC4", accessor_t::type_t::vec4_t},
        {"MAT2", accessor_t::type_t::mat2_t},
        {"MAT3", accessor_t::type_t::mat3_t},
        {"MAT4", accessor_t::type_t::mat4_t},
    };
    auto v = std::string();
    if (!_parse(v, js, err)) return false;
    if (table.find(v) == table.end()) return false;
    val = table[v];
    return true;
}

// Parses a accessor_t object
static inline bool _parse(accessor_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.bufferView, "bufferView", false, js, err))
        return false;
    if (!_parse_attr(val.byteOffset, "byteOffset", false, js, err))
        return false;
    if (!_parse_attr(val.componentType, "componentType", true, js, err))
        return false;
    if (!_parse_attr(val.count, "count", true, js, err)) return false;
    if (!_parse_attr(val.max, "max", true, js, err)) return false;
    if (!_parse_attr(val.min, "min", true, js, err)) return false;
    if (!_parse_attr(val.normalized, "normalized", false, js, err))
        return false;
    if (!_parse_attr(val.sparse, "sparse", false, js, err)) return false;
    if (!_parse_attr(val.type, "type", true, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a componentType_t enum to JSON
static inline bool _dump(
    const accessor_t::componentType_t& val, json& js, _parse_stack& err) {
    static std::map<accessor_t::componentType_t, int> table = {
        {accessor_t::componentType_t::byte_t, 5120},
        {accessor_t::componentType_t::unsigned_byte_t, 5121},
        {accessor_t::componentType_t::short_t, 5122},
        {accessor_t::componentType_t::unsigned_short_t, 5123},
        {accessor_t::componentType_t::unsigned_int_t, 5125},
        {accessor_t::componentType_t::float_t, 5126},
    };
    auto v = table[val];
    if (!_dump(v, js, err)) return false;
    return true;
}

// Converts a type_t enum to JSON
static inline bool _dump(
    const accessor_t::type_t& val, json& js, _parse_stack& err) {
    static std::map<accessor_t::type_t, std::string> table = {
        {accessor_t::type_t::scalar_t, "SCALAR"},
        {accessor_t::type_t::vec2_t, "VEC2"},
        {accessor_t::type_t::vec3_t, "VEC3"},
        {accessor_t::type_t::vec4_t, "VEC4"},
        {accessor_t::type_t::mat2_t, "MAT2"},
        {accessor_t::type_t::mat3_t, "MAT3"},
        {accessor_t::type_t::mat4_t, "MAT4"},
    };
    auto v = table[val];
    if (!_dump(v, js, err)) return false;
    return true;
}

// Converts a accessor_t object to JSON
static inline bool _dump(const accessor_t& val, json& js, _parse_stack& err) {
    static const auto defval = accessor_t();
    if (!_dump((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(
            val.bufferView, "bufferView", defval.bufferView, false, js, err))
        return false;
    if (!_dump_attr(
            val.byteOffset, "byteOffset", defval.byteOffset, false, js, err))
        return false;
    if (!_dump_attr(val.componentType, "componentType", defval.componentType,
            true, js, err))
        return false;
    if (!_dump_attr(val.count, "count", defval.count, true, js, err))
        return false;
    if (!_dump_attr(val.max, "max", defval.max, true, js, err)) return false;
    if (!_dump_attr(val.min, "min", defval.min, true, js, err)) return false;
    if (!_dump_attr(
            val.normalized, "normalized", defval.normalized, false, js, err))
        return false;
    if (!_dump_attr(val.sparse, "sparse", defval.sparse, false, js, err))
        return false;
    if (!_dump_attr(val.type, "type", defval.type, true, js, err)) return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(
    const animation_channel_target_t& a, const animation_channel_target_t& b) {
    if (!((glTFProperty_t&)a == (glTFProperty_t&)b)) return false;
    if (!(a.node == b.node)) return false;
    if (!(a.path == b.path)) return false;
    return true;
}

// Parse a path_t enum
static inline bool _parse(animation_channel_target_t::path_t& val,
    const json& js, _parse_stack& err) {
    static std::map<std::string, animation_channel_target_t::path_t> table = {
        {"translation", animation_channel_target_t::path_t::translation_t},
        {"rotation", animation_channel_target_t::path_t::rotation_t},
        {"scale", animation_channel_target_t::path_t::scale_t},
    };
    auto v = std::string();
    if (!_parse(v, js, err)) return false;
    if (table.find(v) == table.end()) return false;
    val = table[v];
    return true;
}

// Parses a animation_channel_target_t object
static inline bool _parse(
    animation_channel_target_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.node, "node", true, js, err)) return false;
    if (!_parse_attr(val.path, "path", true, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a path_t enum to JSON
static inline bool _dump(const animation_channel_target_t::path_t& val,
    json& js, _parse_stack& err) {
    static std::map<animation_channel_target_t::path_t, std::string> table = {
        {animation_channel_target_t::path_t::translation_t, "translation"},
        {animation_channel_target_t::path_t::rotation_t, "rotation"},
        {animation_channel_target_t::path_t::scale_t, "scale"},
    };
    auto v = table[val];
    if (!_dump(v, js, err)) return false;
    return true;
}

// Converts a animation_channel_target_t object to JSON
static inline bool _dump(
    const animation_channel_target_t& val, json& js, _parse_stack& err) {
    static const auto defval = animation_channel_target_t();
    if (!_dump((glTFProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(val.node, "node", defval.node, true, js, err)) return false;
    if (!_dump_attr(val.path, "path", defval.path, true, js, err)) return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(
    const animation_channel_t& a, const animation_channel_t& b) {
    if (!((glTFProperty_t&)a == (glTFProperty_t&)b)) return false;
    if (!(a.sampler == b.sampler)) return false;
    if (!(a.target == b.target)) return false;
    return true;
}

// Parses a animation_channel_t object
static inline bool _parse(
    animation_channel_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.sampler, "sampler", true, js, err)) return false;
    if (!_parse_attr(val.target, "target", true, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a animation_channel_t object to JSON
static inline bool _dump(
    const animation_channel_t& val, json& js, _parse_stack& err) {
    static const auto defval = animation_channel_t();
    if (!_dump((glTFProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(val.sampler, "sampler", defval.sampler, true, js, err))
        return false;
    if (!_dump_attr(val.target, "target", defval.target, true, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(
    const animation_sampler_t& a, const animation_sampler_t& b) {
    if (!((glTFProperty_t&)a == (glTFProperty_t&)b)) return false;
    if (!(a.input == b.input)) return false;
    if (!(a.interpolation == b.interpolation)) return false;
    if (!(a.output == b.output)) return false;
    return true;
}

// Parse a interpolation_t enum
static inline bool _parse(animation_sampler_t::interpolation_t& val,
    const json& js, _parse_stack& err) {
    static std::map<std::string, animation_sampler_t::interpolation_t> table = {
        {"LINEAR", animation_sampler_t::interpolation_t::linear_t},
        {"STEP", animation_sampler_t::interpolation_t::step_t},
    };
    auto v = std::string();
    if (!_parse(v, js, err)) return false;
    if (table.find(v) == table.end()) return false;
    val = table[v];
    return true;
}

// Parses a animation_sampler_t object
static inline bool _parse(
    animation_sampler_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.input, "input", true, js, err)) return false;
    if (!_parse_attr(val.interpolation, "interpolation", false, js, err))
        return false;
    if (!_parse_attr(val.output, "output", true, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a interpolation_t enum to JSON
static inline bool _dump(const animation_sampler_t::interpolation_t& val,
    json& js, _parse_stack& err) {
    static std::map<animation_sampler_t::interpolation_t, std::string> table = {
        {animation_sampler_t::interpolation_t::linear_t, "LINEAR"},
        {animation_sampler_t::interpolation_t::step_t, "STEP"},
    };
    auto v = table[val];
    if (!_dump(v, js, err)) return false;
    return true;
}

// Converts a animation_sampler_t object to JSON
static inline bool _dump(
    const animation_sampler_t& val, json& js, _parse_stack& err) {
    static const auto defval = animation_sampler_t();
    if (!_dump((glTFProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(val.input, "input", defval.input, true, js, err))
        return false;
    if (!_dump_attr(val.interpolation, "interpolation", defval.interpolation,
            false, js, err))
        return false;
    if (!_dump_attr(val.output, "output", defval.output, true, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(const animation_t& a, const animation_t& b) {
    if (!((glTFChildOfRootProperty_t&)a == (glTFChildOfRootProperty_t&)b))
        return false;
    if (!(a.channels == b.channels)) return false;
    if (!(a.samplers == b.samplers)) return false;
    return true;
}

// Parses a animation_t object
static inline bool _parse(animation_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.channels, "channels", true, js, err)) return false;
    if (!_parse_attr(val.samplers, "samplers", true, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a animation_t object to JSON
static inline bool _dump(const animation_t& val, json& js, _parse_stack& err) {
    static const auto defval = animation_t();
    if (!_dump((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(val.channels, "channels", defval.channels, true, js, err))
        return false;
    if (!_dump_attr(val.samplers, "samplers", defval.samplers, true, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(const asset_t& a, const asset_t& b) {
    if (!((glTFProperty_t&)a == (glTFProperty_t&)b)) return false;
    if (!(a.copyright == b.copyright)) return false;
    if (!(a.generator == b.generator)) return false;
    if (!(a.version == b.version)) return false;
    return true;
}

// Parse a version_t enum
static inline bool _parse(
    asset_t::version_t& val, const json& js, _parse_stack& err) {
    static std::map<std::string, asset_t::version_t> table = {
        {"2.0", asset_t::version_t::_2_0_t},
    };
    auto v = std::string();
    if (!_parse(v, js, err)) return false;
    if (table.find(v) == table.end()) return false;
    val = table[v];
    return true;
}

// Parses a asset_t object
static inline bool _parse(asset_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.copyright, "copyright", false, js, err)) return false;
    if (!_parse_attr(val.generator, "generator", false, js, err)) return false;
    if (!_parse_attr(val.version, "version", true, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a version_t enum to JSON
static inline bool _dump(
    const asset_t::version_t& val, json& js, _parse_stack& err) {
    static std::map<asset_t::version_t, std::string> table = {
        {asset_t::version_t::_2_0_t, "2.0"},
    };
    auto v = table[val];
    if (!_dump(v, js, err)) return false;
    return true;
}

// Converts a asset_t object to JSON
static inline bool _dump(const asset_t& val, json& js, _parse_stack& err) {
    static const auto defval = asset_t();
    if (!_dump((glTFProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(
            val.copyright, "copyright", defval.copyright, false, js, err))
        return false;
    if (!_dump_attr(
            val.generator, "generator", defval.generator, false, js, err))
        return false;
    if (!_dump_attr(val.version, "version", defval.version, true, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(const buffer_t& a, const buffer_t& b) {
    if (!((glTFChildOfRootProperty_t&)a == (glTFChildOfRootProperty_t&)b))
        return false;
    if (!(a.byteLength == b.byteLength)) return false;
    if (!(a.uri == b.uri)) return false;
    return true;
}

// Parses a buffer_t object
static inline bool _parse(buffer_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.byteLength, "byteLength", true, js, err)) return false;
    if (!_parse_attr(val.uri, "uri", false, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a buffer_t object to JSON
static inline bool _dump(const buffer_t& val, json& js, _parse_stack& err) {
    static const auto defval = buffer_t();
    if (!_dump((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(
            val.byteLength, "byteLength", defval.byteLength, true, js, err))
        return false;
    if (!_dump_attr(val.uri, "uri", defval.uri, false, js, err)) return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(const bufferView_t& a, const bufferView_t& b) {
    if (!((glTFChildOfRootProperty_t&)a == (glTFChildOfRootProperty_t&)b))
        return false;
    if (!(a.buffer == b.buffer)) return false;
    if (!(a.byteLength == b.byteLength)) return false;
    if (!(a.byteOffset == b.byteOffset)) return false;
    if (!(a.byteStride == b.byteStride)) return false;
    if (!(a.target == b.target)) return false;
    return true;
}

// Parse a target_t enum
static inline bool _parse(
    bufferView_t::target_t& val, const json& js, _parse_stack& err) {
    static std::map<int, bufferView_t::target_t> table = {
        {34962, bufferView_t::target_t::array_buffer_t},
        {34963, bufferView_t::target_t::element_array_buffer_t},
    };
    auto v = 0;
    if (!_parse(v, js, err)) return false;
    if (table.find(v) == table.end()) return false;
    val = table[v];
    return true;
}

// Parses a bufferView_t object
static inline bool _parse(
    bufferView_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.buffer, "buffer", true, js, err)) return false;
    if (!_parse_attr(val.byteLength, "byteLength", true, js, err)) return false;
    if (!_parse_attr(val.byteOffset, "byteOffset", true, js, err)) return false;
    if (!_parse_attr(val.byteStride, "byteStride", false, js, err))
        return false;
    if (!_parse_attr(val.target, "target", false, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a target_t enum to JSON
static inline bool _dump(
    const bufferView_t::target_t& val, json& js, _parse_stack& err) {
    static std::map<bufferView_t::target_t, int> table = {
        {bufferView_t::target_t::array_buffer_t, 34962},
        {bufferView_t::target_t::element_array_buffer_t, 34963},
    };
    auto v = table[val];
    if (!_dump(v, js, err)) return false;
    return true;
}

// Converts a bufferView_t object to JSON
static inline bool _dump(const bufferView_t& val, json& js, _parse_stack& err) {
    static const auto defval = bufferView_t();
    if (!_dump((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(val.buffer, "buffer", defval.buffer, true, js, err))
        return false;
    if (!_dump_attr(
            val.byteLength, "byteLength", defval.byteLength, true, js, err))
        return false;
    if (!_dump_attr(
            val.byteOffset, "byteOffset", defval.byteOffset, true, js, err))
        return false;
    if (!_dump_attr(
            val.byteStride, "byteStride", defval.byteStride, false, js, err))
        return false;
    if (!_dump_attr(val.target, "target", defval.target, false, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(
    const camera_orthographic_t& a, const camera_orthographic_t& b) {
    if (!((glTFProperty_t&)a == (glTFProperty_t&)b)) return false;
    if (!(a.xmag == b.xmag)) return false;
    if (!(a.ymag == b.ymag)) return false;
    if (!(a.zfar == b.zfar)) return false;
    if (!(a.znear == b.znear)) return false;
    return true;
}

// Parses a camera_orthographic_t object
static inline bool _parse(
    camera_orthographic_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.xmag, "xmag", true, js, err)) return false;
    if (!_parse_attr(val.ymag, "ymag", true, js, err)) return false;
    if (!_parse_attr(val.zfar, "zfar", true, js, err)) return false;
    if (!_parse_attr(val.znear, "znear", true, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a camera_orthographic_t object to JSON
static inline bool _dump(
    const camera_orthographic_t& val, json& js, _parse_stack& err) {
    static const auto defval = camera_orthographic_t();
    if (!_dump((glTFProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(val.xmag, "xmag", defval.xmag, true, js, err)) return false;
    if (!_dump_attr(val.ymag, "ymag", defval.ymag, true, js, err)) return false;
    if (!_dump_attr(val.zfar, "zfar", defval.zfar, true, js, err)) return false;
    if (!_dump_attr(val.znear, "znear", defval.znear, true, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(
    const camera_perspective_t& a, const camera_perspective_t& b) {
    if (!((glTFProperty_t&)a == (glTFProperty_t&)b)) return false;
    if (!(a.aspectRatio == b.aspectRatio)) return false;
    if (!(a.yfov == b.yfov)) return false;
    if (!(a.zfar == b.zfar)) return false;
    if (!(a.znear == b.znear)) return false;
    return true;
}

// Parses a camera_perspective_t object
static inline bool _parse(
    camera_perspective_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.aspectRatio, "aspectRatio", false, js, err))
        return false;
    if (!_parse_attr(val.yfov, "yfov", true, js, err)) return false;
    if (!_parse_attr(val.zfar, "zfar", false, js, err)) return false;
    if (!_parse_attr(val.znear, "znear", true, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a camera_perspective_t object to JSON
static inline bool _dump(
    const camera_perspective_t& val, json& js, _parse_stack& err) {
    static const auto defval = camera_perspective_t();
    if (!_dump((glTFProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(
            val.aspectRatio, "aspectRatio", defval.aspectRatio, false, js, err))
        return false;
    if (!_dump_attr(val.yfov, "yfov", defval.yfov, true, js, err)) return false;
    if (!_dump_attr(val.zfar, "zfar", defval.zfar, false, js, err))
        return false;
    if (!_dump_attr(val.znear, "znear", defval.znear, true, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(const camera_t& a, const camera_t& b) {
    if (!((glTFChildOfRootProperty_t&)a == (glTFChildOfRootProperty_t&)b))
        return false;
    if (!(a.orthographic == b.orthographic)) return false;
    if (!(a.perspective == b.perspective)) return false;
    if (!(a.type == b.type)) return false;
    return true;
}

// Parse a type_t enum
static inline bool _parse(
    camera_t::type_t& val, const json& js, _parse_stack& err) {
    static std::map<std::string, camera_t::type_t> table = {
        {"perspective", camera_t::type_t::perspective_t},
        {"orthographic", camera_t::type_t::orthographic_t},
    };
    auto v = std::string();
    if (!_parse(v, js, err)) return false;
    if (table.find(v) == table.end()) return false;
    val = table[v];
    return true;
}

// Parses a camera_t object
static inline bool _parse(camera_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.orthographic, "orthographic", false, js, err))
        return false;
    if (!_parse_attr(val.perspective, "perspective", false, js, err))
        return false;
    if (!_parse_attr(val.type, "type", true, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a type_t enum to JSON
static inline bool _dump(
    const camera_t::type_t& val, json& js, _parse_stack& err) {
    static std::map<camera_t::type_t, std::string> table = {
        {camera_t::type_t::perspective_t, "perspective"},
        {camera_t::type_t::orthographic_t, "orthographic"},
    };
    auto v = table[val];
    if (!_dump(v, js, err)) return false;
    return true;
}

// Converts a camera_t object to JSON
static inline bool _dump(const camera_t& val, json& js, _parse_stack& err) {
    static const auto defval = camera_t();
    if (!_dump((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(val.orthographic, "orthographic", defval.orthographic,
            false, js, err))
        return false;
    if (!_dump_attr(
            val.perspective, "perspective", defval.perspective, false, js, err))
        return false;
    if (!_dump_attr(val.type, "type", defval.type, true, js, err)) return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(const image_t& a, const image_t& b) {
    if (!((glTFChildOfRootProperty_t&)a == (glTFChildOfRootProperty_t&)b))
        return false;
    if (!(a.bufferView == b.bufferView)) return false;
    if (!(a.mimeType == b.mimeType)) return false;
    if (!(a.uri == b.uri)) return false;
    return true;
}

// Parses a image_t object
static inline bool _parse(image_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.bufferView, "bufferView", false, js, err))
        return false;
    if (!_parse_attr(val.mimeType, "mimeType", false, js, err)) return false;
    if (!_parse_attr(val.uri, "uri", false, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a image_t object to JSON
static inline bool _dump(const image_t& val, json& js, _parse_stack& err) {
    static const auto defval = image_t();
    if (!_dump((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(
            val.bufferView, "bufferView", defval.bufferView, false, js, err))
        return false;
    if (!_dump_attr(val.mimeType, "mimeType", defval.mimeType, false, js, err))
        return false;
    if (!_dump_attr(val.uri, "uri", defval.uri, false, js, err)) return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(
    const textureInfoBase_t& a, const textureInfoBase_t& b) {
    if (!(a.index == b.index)) return false;
    if (!(a.texCoord == b.texCoord)) return false;
    return true;
}

// Parses a textureInfoBase_t object
static inline bool _parse(
    textureInfoBase_t& val, const json& js, _parse_stack& err) {
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.index, "index", true, js, err)) return false;
    if (!_parse_attr(val.texCoord, "texCoord", false, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a textureInfoBase_t object to JSON
static inline bool _dump(
    const textureInfoBase_t& val, json& js, _parse_stack& err) {
    static const auto defval = textureInfoBase_t();
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(val.index, "index", defval.index, true, js, err))
        return false;
    if (!_dump_attr(val.texCoord, "texCoord", defval.texCoord, false, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(const textureInfo_t& a, const textureInfo_t& b) {
    if (!((textureInfoBase_t&)a == (textureInfoBase_t&)b)) return false;
    return true;
}

// Parses a textureInfo_t object
static inline bool _parse(
    textureInfo_t& val, const json& js, _parse_stack& err) {
    if (!_parse((textureInfoBase_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a textureInfo_t object to JSON
static inline bool _dump(
    const textureInfo_t& val, json& js, _parse_stack& err) {
    static const auto defval = textureInfo_t();
    if (!_dump((textureInfoBase_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(const texture_t& a, const texture_t& b) {
    if (!((glTFChildOfRootProperty_t&)a == (glTFChildOfRootProperty_t&)b))
        return false;
    if (!(a.format == b.format)) return false;
    if (!(a.internalFormat == b.internalFormat)) return false;
    if (!(a.sampler == b.sampler)) return false;
    if (!(a.source == b.source)) return false;
    if (!(a.target == b.target)) return false;
    if (!(a.type == b.type)) return false;
    return true;
}

// Parse a format_t enum
static inline bool _parse(
    texture_t::format_t& val, const json& js, _parse_stack& err) {
    static std::map<int, texture_t::format_t> table = {
        {6406, texture_t::format_t::alpha_t},
        {6407, texture_t::format_t::rgb_t}, {6408, texture_t::format_t::rgba_t},
        {6409, texture_t::format_t::luminance_t},
        {6410, texture_t::format_t::luminance_alpha_t},
    };
    auto v = 0;
    if (!_parse(v, js, err)) return false;
    if (table.find(v) == table.end()) return false;
    val = table[v];
    return true;
}

// Parse a internalFormat_t enum
static inline bool _parse(
    texture_t::internalFormat_t& val, const json& js, _parse_stack& err) {
    static std::map<int, texture_t::internalFormat_t> table = {
        {6406, texture_t::internalFormat_t::alpha_t},
        {6407, texture_t::internalFormat_t::rgb_t},
        {6408, texture_t::internalFormat_t::rgba_t},
        {6409, texture_t::internalFormat_t::luminance_t},
        {6410, texture_t::internalFormat_t::luminance_alpha_t},
    };
    auto v = 0;
    if (!_parse(v, js, err)) return false;
    if (table.find(v) == table.end()) return false;
    val = table[v];
    return true;
}

// Parse a target_t enum
static inline bool _parse(
    texture_t::target_t& val, const json& js, _parse_stack& err) {
    static std::map<int, texture_t::target_t> table = {
        {3553, texture_t::target_t::texture_2d_t},
    };
    auto v = 0;
    if (!_parse(v, js, err)) return false;
    if (table.find(v) == table.end()) return false;
    val = table[v];
    return true;
}

// Parse a type_t enum
static inline bool _parse(
    texture_t::type_t& val, const json& js, _parse_stack& err) {
    static std::map<int, texture_t::type_t> table = {
        {5121, texture_t::type_t::unsigned_byte_t},
        {33635, texture_t::type_t::unsigned_short_5_6_5_t},
        {32819, texture_t::type_t::unsigned_short_4_4_4_4_t},
        {32820, texture_t::type_t::unsigned_short_5_5_5_1_t},
    };
    auto v = 0;
    if (!_parse(v, js, err)) return false;
    if (table.find(v) == table.end()) return false;
    val = table[v];
    return true;
}

// Parses a texture_t object
static inline bool _parse(texture_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.format, "format", false, js, err)) return false;
    if (!_parse_attr(val.internalFormat, "internalFormat", false, js, err))
        return false;
    if (!_parse_attr(val.sampler, "sampler", true, js, err)) return false;
    if (!_parse_attr(val.source, "source", true, js, err)) return false;
    if (!_parse_attr(val.target, "target", false, js, err)) return false;
    if (!_parse_attr(val.type, "type", false, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a format_t enum to JSON
static inline bool _dump(
    const texture_t::format_t& val, json& js, _parse_stack& err) {
    static std::map<texture_t::format_t, int> table = {
        {texture_t::format_t::alpha_t, 6406},
        {texture_t::format_t::rgb_t, 6407}, {texture_t::format_t::rgba_t, 6408},
        {texture_t::format_t::luminance_t, 6409},
        {texture_t::format_t::luminance_alpha_t, 6410},
    };
    auto v = table[val];
    if (!_dump(v, js, err)) return false;
    return true;
}

// Converts a internalFormat_t enum to JSON
static inline bool _dump(
    const texture_t::internalFormat_t& val, json& js, _parse_stack& err) {
    static std::map<texture_t::internalFormat_t, int> table = {
        {texture_t::internalFormat_t::alpha_t, 6406},
        {texture_t::internalFormat_t::rgb_t, 6407},
        {texture_t::internalFormat_t::rgba_t, 6408},
        {texture_t::internalFormat_t::luminance_t, 6409},
        {texture_t::internalFormat_t::luminance_alpha_t, 6410},
    };
    auto v = table[val];
    if (!_dump(v, js, err)) return false;
    return true;
}

// Converts a target_t enum to JSON
static inline bool _dump(
    const texture_t::target_t& val, json& js, _parse_stack& err) {
    static std::map<texture_t::target_t, int> table = {
        {texture_t::target_t::texture_2d_t, 3553},
    };
    auto v = table[val];
    if (!_dump(v, js, err)) return false;
    return true;
}

// Converts a type_t enum to JSON
static inline bool _dump(
    const texture_t::type_t& val, json& js, _parse_stack& err) {
    static std::map<texture_t::type_t, int> table = {
        {texture_t::type_t::unsigned_byte_t, 5121},
        {texture_t::type_t::unsigned_short_5_6_5_t, 33635},
        {texture_t::type_t::unsigned_short_4_4_4_4_t, 32819},
        {texture_t::type_t::unsigned_short_5_5_5_1_t, 32820},
    };
    auto v = table[val];
    if (!_dump(v, js, err)) return false;
    return true;
}

// Converts a texture_t object to JSON
static inline bool _dump(const texture_t& val, json& js, _parse_stack& err) {
    static const auto defval = texture_t();
    if (!_dump((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(val.format, "format", defval.format, false, js, err))
        return false;
    if (!_dump_attr(val.internalFormat, "internalFormat", defval.internalFormat,
            false, js, err))
        return false;
    if (!_dump_attr(val.sampler, "sampler", defval.sampler, true, js, err))
        return false;
    if (!_dump_attr(val.source, "source", defval.source, true, js, err))
        return false;
    if (!_dump_attr(val.target, "target", defval.target, false, js, err))
        return false;
    if (!_dump_attr(val.type, "type", defval.type, false, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(const material_normalTextureInfo_t& a,
    const material_normalTextureInfo_t& b) {
    if (!((textureInfoBase_t&)a == (textureInfoBase_t&)b)) return false;
    if (!(a.scale == b.scale)) return false;
    return true;
}

// Parses a material_normalTextureInfo_t object
static inline bool _parse(
    material_normalTextureInfo_t& val, const json& js, _parse_stack& err) {
    if (!_parse((textureInfoBase_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.scale, "scale", false, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a material_normalTextureInfo_t object to JSON
static inline bool _dump(
    const material_normalTextureInfo_t& val, json& js, _parse_stack& err) {
    static const auto defval = material_normalTextureInfo_t();
    if (!_dump((textureInfoBase_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(val.scale, "scale", defval.scale, false, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(const material_occlusionTextureInfo_t& a,
    const material_occlusionTextureInfo_t& b) {
    if (!((textureInfoBase_t&)a == (textureInfoBase_t&)b)) return false;
    if (!(a.strength == b.strength)) return false;
    return true;
}

// Parses a material_occlusionTextureInfo_t object
static inline bool _parse(
    material_occlusionTextureInfo_t& val, const json& js, _parse_stack& err) {
    if (!_parse((textureInfoBase_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.strength, "strength", false, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a material_occlusionTextureInfo_t object to JSON
static inline bool _dump(
    const material_occlusionTextureInfo_t& val, json& js, _parse_stack& err) {
    static const auto defval = material_occlusionTextureInfo_t();
    if (!_dump((textureInfoBase_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(val.strength, "strength", defval.strength, false, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(const material_pbrMetallicRoughness_t& a,
    const material_pbrMetallicRoughness_t& b) {
    if (!(a.baseColorFactor == b.baseColorFactor)) return false;
    if (!(a.baseColorTexture == b.baseColorTexture)) return false;
    if (!(a.metallicFactor == b.metallicFactor)) return false;
    if (!(a.metallicRoughnessTexture == b.metallicRoughnessTexture))
        return false;
    if (!(a.roughnessFactor == b.roughnessFactor)) return false;
    return true;
}

// Parses a material_pbrMetallicRoughness_t object
static inline bool _parse(
    material_pbrMetallicRoughness_t& val, const json& js, _parse_stack& err) {
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.baseColorFactor, "baseColorFactor", false, js, err))
        return false;
    if (!_parse_attr(val.baseColorTexture, "baseColorTexture", false, js, err))
        return false;
    if (!_parse_attr(val.metallicFactor, "metallicFactor", false, js, err))
        return false;
    if (!_parse_attr(val.metallicRoughnessTexture, "metallicRoughnessTexture",
            false, js, err))
        return false;
    if (!_parse_attr(val.roughnessFactor, "roughnessFactor", false, js, err))
        return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a material_pbrMetallicRoughness_t object to JSON
static inline bool _dump(
    const material_pbrMetallicRoughness_t& val, json& js, _parse_stack& err) {
    static const auto defval = material_pbrMetallicRoughness_t();
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(val.baseColorFactor, "baseColorFactor",
            defval.baseColorFactor, false, js, err))
        return false;
    if (!_dump_attr(val.baseColorTexture, "baseColorTexture",
            defval.baseColorTexture, false, js, err))
        return false;
    if (!_dump_attr(val.metallicFactor, "metallicFactor", defval.metallicFactor,
            false, js, err))
        return false;
    if (!_dump_attr(val.metallicRoughnessTexture, "metallicRoughnessTexture",
            defval.metallicRoughnessTexture, false, js, err))
        return false;
    if (!_dump_attr(val.roughnessFactor, "roughnessFactor",
            defval.roughnessFactor, false, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Parses a arrayValues object

static inline bool _parse(
    arrayValues_t& val, const json& js, _parse_stack& err) {
    // alsp support constants
    if (_parse(val.items_number, js, err)) return true;
    if (_parse(val.items_string, js, err)) return true;
    if (_parse(val.items_boolean, js, err)) return true;
    {
        float v;
        if (_parse(v, js, err)) {
            val.items_number.push_back(v);
            return true;
        }
    }
    {
        std::string v;
        if (_parse(v, js, err)) {
            val.items_string.push_back(v);
            return true;
        }
    }
    {
        bool v;
        if (_parse(v, js, err)) {
            val.items_boolean.push_back(v);
            return true;
        }
    }
    return false;
}

// Equality check

static inline bool operator==(const arrayValues_t& a, const arrayValues_t& b) {
    return a.items_number == b.items_number &&
           a.items_string == b.items_string &&
           a.items_boolean == b.items_boolean;
}

// Converts a arrayValues object to JSON

static inline bool _dump(
    const arrayValues_t& val, json& js, _parse_stack& err) {
    if (!val.items_number.empty()) {
        if (_dump(val.items_number, js, err)) return false;
    } else if (!val.items_string.empty()) {
        if (_dump(val.items_string, js, err)) return false;
    } else if (!val.items_boolean.empty()) {
        if (_dump(val.items_boolean, js, err)) return false;
    } else {
        if (_dump(std::vector<int>(), js, err)) return false;
    }
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(const material_t& a, const material_t& b) {
    if (!((glTFChildOfRootProperty_t&)a == (glTFChildOfRootProperty_t&)b))
        return false;
    if (!(a.emissiveFactor == b.emissiveFactor)) return false;
    if (!(a.emissiveTexture == b.emissiveTexture)) return false;
    if (!(a.normalTexture == b.normalTexture)) return false;
    if (!(a.occlusionTexture == b.occlusionTexture)) return false;
    if (!(a.pbrMetallicRoughness == b.pbrMetallicRoughness)) return false;
    return true;
}

// Parses a material_t object
static inline bool _parse(material_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.emissiveFactor, "emissiveFactor", false, js, err))
        return false;
    if (!_parse_attr(val.emissiveTexture, "emissiveTexture", false, js, err))
        return false;
    if (!_parse_attr(val.normalTexture, "normalTexture", false, js, err))
        return false;
    if (!_parse_attr(val.occlusionTexture, "occlusionTexture", false, js, err))
        return false;
    if (!_parse_attr(
            val.pbrMetallicRoughness, "pbrMetallicRoughness", false, js, err))
        return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a material_t object to JSON
static inline bool _dump(const material_t& val, json& js, _parse_stack& err) {
    static const auto defval = material_t();
    if (!_dump((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(val.emissiveFactor, "emissiveFactor", defval.emissiveFactor,
            false, js, err))
        return false;
    if (!_dump_attr(val.emissiveTexture, "emissiveTexture",
            defval.emissiveTexture, false, js, err))
        return false;
    if (!_dump_attr(val.normalTexture, "normalTexture", defval.normalTexture,
            false, js, err))
        return false;
    if (!_dump_attr(val.occlusionTexture, "occlusionTexture",
            defval.occlusionTexture, false, js, err))
        return false;
    if (!_dump_attr(val.pbrMetallicRoughness, "pbrMetallicRoughness",
            defval.pbrMetallicRoughness, false, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(
    const mesh_primitive_t& a, const mesh_primitive_t& b) {
    if (!((glTFProperty_t&)a == (glTFProperty_t&)b)) return false;
    if (!(a.attributes == b.attributes)) return false;
    if (!(a.indices == b.indices)) return false;
    if (!(a.material == b.material)) return false;
    if (!(a.mode == b.mode)) return false;
    if (!(a.targets == b.targets)) return false;
    return true;
}

// Parse a mode_t enum
static inline bool _parse(
    mesh_primitive_t::mode_t& val, const json& js, _parse_stack& err) {
    static std::map<int, mesh_primitive_t::mode_t> table = {
        {0, mesh_primitive_t::mode_t::points_t},
        {1, mesh_primitive_t::mode_t::lines_t},
        {2, mesh_primitive_t::mode_t::line_loop_t},
        {3, mesh_primitive_t::mode_t::line_strip_t},
        {4, mesh_primitive_t::mode_t::triangles_t},
        {5, mesh_primitive_t::mode_t::triangle_strip_t},
        {6, mesh_primitive_t::mode_t::triangle_fan_t},
    };
    auto v = 0;
    if (!_parse(v, js, err)) return false;
    if (table.find(v) == table.end()) return false;
    val = table[v];
    return true;
}

// Parses a mesh_primitive_t object
static inline bool _parse(
    mesh_primitive_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.attributes, "attributes", true, js, err)) return false;
    if (!_parse_attr(val.indices, "indices", false, js, err)) return false;
    if (!_parse_attr(val.material, "material", false, js, err)) return false;
    if (!_parse_attr(val.mode, "mode", false, js, err)) return false;
    if (!_parse_attr(val.targets, "targets", false, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a mode_t enum to JSON
static inline bool _dump(
    const mesh_primitive_t::mode_t& val, json& js, _parse_stack& err) {
    static std::map<mesh_primitive_t::mode_t, int> table = {
        {mesh_primitive_t::mode_t::points_t, 0},
        {mesh_primitive_t::mode_t::lines_t, 1},
        {mesh_primitive_t::mode_t::line_loop_t, 2},
        {mesh_primitive_t::mode_t::line_strip_t, 3},
        {mesh_primitive_t::mode_t::triangles_t, 4},
        {mesh_primitive_t::mode_t::triangle_strip_t, 5},
        {mesh_primitive_t::mode_t::triangle_fan_t, 6},
    };
    auto v = table[val];
    if (!_dump(v, js, err)) return false;
    return true;
}

// Converts a mesh_primitive_t object to JSON
static inline bool _dump(
    const mesh_primitive_t& val, json& js, _parse_stack& err) {
    static const auto defval = mesh_primitive_t();
    if (!_dump((glTFProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(
            val.attributes, "attributes", defval.attributes, true, js, err))
        return false;
    if (!_dump_attr(val.indices, "indices", defval.indices, false, js, err))
        return false;
    if (!_dump_attr(val.material, "material", defval.material, false, js, err))
        return false;
    if (!_dump_attr(val.mode, "mode", defval.mode, false, js, err))
        return false;
    if (!_dump_attr(val.targets, "targets", defval.targets, false, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(const mesh_t& a, const mesh_t& b) {
    if (!((glTFChildOfRootProperty_t&)a == (glTFChildOfRootProperty_t&)b))
        return false;
    if (!(a.primitives == b.primitives)) return false;
    if (!(a.weights == b.weights)) return false;
    return true;
}

// Parses a mesh_t object
static inline bool _parse(mesh_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.primitives, "primitives", true, js, err)) return false;
    if (!_parse_attr(val.weights, "weights", false, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a mesh_t object to JSON
static inline bool _dump(const mesh_t& val, json& js, _parse_stack& err) {
    static const auto defval = mesh_t();
    if (!_dump((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(
            val.primitives, "primitives", defval.primitives, true, js, err))
        return false;
    if (!_dump_attr(val.weights, "weights", defval.weights, false, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(const node_t& a, const node_t& b) {
    if (!((glTFChildOfRootProperty_t&)a == (glTFChildOfRootProperty_t&)b))
        return false;
    if (!(a.camera == b.camera)) return false;
    if (!(a.children == b.children)) return false;
    if (!(a.matrix == b.matrix)) return false;
    if (!(a.mesh == b.mesh)) return false;
    if (!(a.rotation == b.rotation)) return false;
    if (!(a.scale == b.scale)) return false;
    if (!(a.skin == b.skin)) return false;
    if (!(a.translation == b.translation)) return false;
    if (!(a.weights == b.weights)) return false;
    return true;
}

// Parses a node_t object
static inline bool _parse(node_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.camera, "camera", false, js, err)) return false;
    if (!_parse_attr(val.children, "children", false, js, err)) return false;
    if (!_parse_attr(val.matrix, "matrix", false, js, err)) return false;
    if (!_parse_attr(val.mesh, "mesh", false, js, err)) return false;
    if (!_parse_attr(val.rotation, "rotation", false, js, err)) return false;
    if (!_parse_attr(val.scale, "scale", false, js, err)) return false;
    if (!_parse_attr(val.skin, "skin", false, js, err)) return false;
    if (!_parse_attr(val.translation, "translation", false, js, err))
        return false;
    if (!_parse_attr(val.weights, "weights", false, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a node_t object to JSON
static inline bool _dump(const node_t& val, json& js, _parse_stack& err) {
    static const auto defval = node_t();
    if (!_dump((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(val.camera, "camera", defval.camera, false, js, err))
        return false;
    if (!_dump_attr(val.children, "children", defval.children, false, js, err))
        return false;
    if (!_dump_attr(val.matrix, "matrix", defval.matrix, false, js, err))
        return false;
    if (!_dump_attr(val.mesh, "mesh", defval.mesh, false, js, err))
        return false;
    if (!_dump_attr(val.rotation, "rotation", defval.rotation, false, js, err))
        return false;
    if (!_dump_attr(val.scale, "scale", defval.scale, false, js, err))
        return false;
    if (!_dump_attr(val.skin, "skin", defval.skin, false, js, err))
        return false;
    if (!_dump_attr(
            val.translation, "translation", defval.translation, false, js, err))
        return false;
    if (!_dump_attr(val.weights, "weights", defval.weights, false, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(const sampler_t& a, const sampler_t& b) {
    if (!((glTFChildOfRootProperty_t&)a == (glTFChildOfRootProperty_t&)b))
        return false;
    if (!(a.magFilter == b.magFilter)) return false;
    if (!(a.minFilter == b.minFilter)) return false;
    if (!(a.wrapS == b.wrapS)) return false;
    if (!(a.wrapT == b.wrapT)) return false;
    return true;
}

// Parse a magFilter_t enum
static inline bool _parse(
    sampler_t::magFilter_t& val, const json& js, _parse_stack& err) {
    static std::map<int, sampler_t::magFilter_t> table = {
        {9728, sampler_t::magFilter_t::nearest_t},
        {9729, sampler_t::magFilter_t::linear_t},
    };
    auto v = 0;
    if (!_parse(v, js, err)) return false;
    if (table.find(v) == table.end()) return false;
    val = table[v];
    return true;
}

// Parse a minFilter_t enum
static inline bool _parse(
    sampler_t::minFilter_t& val, const json& js, _parse_stack& err) {
    static std::map<int, sampler_t::minFilter_t> table = {
        {9728, sampler_t::minFilter_t::nearest_t},
        {9729, sampler_t::minFilter_t::linear_t},
        {9984, sampler_t::minFilter_t::nearest_mipmap_nearest_t},
        {9985, sampler_t::minFilter_t::linear_mipmap_nearest_t},
        {9986, sampler_t::minFilter_t::nearest_mipmap_linear_t},
        {9987, sampler_t::minFilter_t::linear_mipmap_linear_t},
    };
    auto v = 0;
    if (!_parse(v, js, err)) return false;
    if (table.find(v) == table.end()) return false;
    val = table[v];
    return true;
}

// Parse a wrapS_t enum
static inline bool _parse(
    sampler_t::wrapS_t& val, const json& js, _parse_stack& err) {
    static std::map<int, sampler_t::wrapS_t> table = {
        {33071, sampler_t::wrapS_t::clamp_to_edge_t},
        {33648, sampler_t::wrapS_t::mirrored_repeat_t},
        {10497, sampler_t::wrapS_t::repeat_t},
    };
    auto v = 0;
    if (!_parse(v, js, err)) return false;
    if (table.find(v) == table.end()) return false;
    val = table[v];
    return true;
}

// Parse a wrapT_t enum
static inline bool _parse(
    sampler_t::wrapT_t& val, const json& js, _parse_stack& err) {
    static std::map<int, sampler_t::wrapT_t> table = {
        {33071, sampler_t::wrapT_t::clamp_to_edge_t},
        {33648, sampler_t::wrapT_t::mirrored_repeat_t},
        {10497, sampler_t::wrapT_t::repeat_t},
    };
    auto v = 0;
    if (!_parse(v, js, err)) return false;
    if (table.find(v) == table.end()) return false;
    val = table[v];
    return true;
}

// Parses a sampler_t object
static inline bool _parse(sampler_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.magFilter, "magFilter", false, js, err)) return false;
    if (!_parse_attr(val.minFilter, "minFilter", false, js, err)) return false;
    if (!_parse_attr(val.wrapS, "wrapS", false, js, err)) return false;
    if (!_parse_attr(val.wrapT, "wrapT", false, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a magFilter_t enum to JSON
static inline bool _dump(
    const sampler_t::magFilter_t& val, json& js, _parse_stack& err) {
    static std::map<sampler_t::magFilter_t, int> table = {
        {sampler_t::magFilter_t::nearest_t, 9728},
        {sampler_t::magFilter_t::linear_t, 9729},
    };
    auto v = table[val];
    if (!_dump(v, js, err)) return false;
    return true;
}

// Converts a minFilter_t enum to JSON
static inline bool _dump(
    const sampler_t::minFilter_t& val, json& js, _parse_stack& err) {
    static std::map<sampler_t::minFilter_t, int> table = {
        {sampler_t::minFilter_t::nearest_t, 9728},
        {sampler_t::minFilter_t::linear_t, 9729},
        {sampler_t::minFilter_t::nearest_mipmap_nearest_t, 9984},
        {sampler_t::minFilter_t::linear_mipmap_nearest_t, 9985},
        {sampler_t::minFilter_t::nearest_mipmap_linear_t, 9986},
        {sampler_t::minFilter_t::linear_mipmap_linear_t, 9987},
    };
    auto v = table[val];
    if (!_dump(v, js, err)) return false;
    return true;
}

// Converts a wrapS_t enum to JSON
static inline bool _dump(
    const sampler_t::wrapS_t& val, json& js, _parse_stack& err) {
    static std::map<sampler_t::wrapS_t, int> table = {
        {sampler_t::wrapS_t::clamp_to_edge_t, 33071},
        {sampler_t::wrapS_t::mirrored_repeat_t, 33648},
        {sampler_t::wrapS_t::repeat_t, 10497},
    };
    auto v = table[val];
    if (!_dump(v, js, err)) return false;
    return true;
}

// Converts a wrapT_t enum to JSON
static inline bool _dump(
    const sampler_t::wrapT_t& val, json& js, _parse_stack& err) {
    static std::map<sampler_t::wrapT_t, int> table = {
        {sampler_t::wrapT_t::clamp_to_edge_t, 33071},
        {sampler_t::wrapT_t::mirrored_repeat_t, 33648},
        {sampler_t::wrapT_t::repeat_t, 10497},
    };
    auto v = table[val];
    if (!_dump(v, js, err)) return false;
    return true;
}

// Converts a sampler_t object to JSON
static inline bool _dump(const sampler_t& val, json& js, _parse_stack& err) {
    static const auto defval = sampler_t();
    if (!_dump((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(
            val.magFilter, "magFilter", defval.magFilter, false, js, err))
        return false;
    if (!_dump_attr(
            val.minFilter, "minFilter", defval.minFilter, false, js, err))
        return false;
    if (!_dump_attr(val.wrapS, "wrapS", defval.wrapS, false, js, err))
        return false;
    if (!_dump_attr(val.wrapT, "wrapT", defval.wrapT, false, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(const scene_t& a, const scene_t& b) {
    if (!((glTFChildOfRootProperty_t&)a == (glTFChildOfRootProperty_t&)b))
        return false;
    if (!(a.nodes == b.nodes)) return false;
    return true;
}

// Parses a scene_t object
static inline bool _parse(scene_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.nodes, "nodes", false, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a scene_t object to JSON
static inline bool _dump(const scene_t& val, json& js, _parse_stack& err) {
    static const auto defval = scene_t();
    if (!_dump((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(val.nodes, "nodes", defval.nodes, false, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(const skin_t& a, const skin_t& b) {
    if (!((glTFChildOfRootProperty_t&)a == (glTFChildOfRootProperty_t&)b))
        return false;
    if (!(a.inverseBindMatrices == b.inverseBindMatrices)) return false;
    if (!(a.joints == b.joints)) return false;
    if (!(a.skeleton == b.skeleton)) return false;
    return true;
}

// Parses a skin_t object
static inline bool _parse(skin_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(
            val.inverseBindMatrices, "inverseBindMatrices", false, js, err))
        return false;
    if (!_parse_attr(val.joints, "joints", true, js, err)) return false;
    if (!_parse_attr(val.skeleton, "skeleton", false, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a skin_t object to JSON
static inline bool _dump(const skin_t& val, json& js, _parse_stack& err) {
    static const auto defval = skin_t();
    if (!_dump((glTFChildOfRootProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(val.inverseBindMatrices, "inverseBindMatrices",
            defval.inverseBindMatrices, false, js, err))
        return false;
    if (!_dump_attr(val.joints, "joints", defval.joints, true, js, err))
        return false;
    if (!_dump_attr(val.skeleton, "skeleton", defval.skeleton, false, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// Equality check for defaults (might go away in the future)
static inline bool operator==(const glTF_t& a, const glTF_t& b) {
    if (!((glTFProperty_t&)a == (glTFProperty_t&)b)) return false;
    if (!(a.accessors == b.accessors)) return false;
    if (!(a.animations == b.animations)) return false;
    if (!(a.asset == b.asset)) return false;
    if (!(a.bufferViews == b.bufferViews)) return false;
    if (!(a.buffers == b.buffers)) return false;
    if (!(a.cameras == b.cameras)) return false;
    if (!(a.extensionsRequired == b.extensionsRequired)) return false;
    if (!(a.extensionsUsed == b.extensionsUsed)) return false;
    if (!(a.images == b.images)) return false;
    if (!(a.materials == b.materials)) return false;
    if (!(a.meshes == b.meshes)) return false;
    if (!(a.nodes == b.nodes)) return false;
    if (!(a.samplers == b.samplers)) return false;
    if (!(a.scene == b.scene)) return false;
    if (!(a.scenes == b.scenes)) return false;
    if (!(a.skins == b.skins)) return false;
    if (!(a.textures == b.textures)) return false;
    return true;
}

// Parses a glTF_t object
static inline bool _parse(glTF_t& val, const json& js, _parse_stack& err) {
    if (!_parse((glTFProperty_t&)val, js, err)) return false;
    if (!_parse_begin_obj(js, err)) return false;
    if (!_parse_attr(val.accessors, "accessors", false, js, err)) return false;
    if (!_parse_attr(val.animations, "animations", false, js, err))
        return false;
    if (!_parse_attr(val.asset, "asset", true, js, err)) return false;
    if (!_parse_attr(val.bufferViews, "bufferViews", false, js, err))
        return false;
    if (!_parse_attr(val.buffers, "buffers", false, js, err)) return false;
    if (!_parse_attr(val.cameras, "cameras", false, js, err)) return false;
    if (!_parse_attr(
            val.extensionsRequired, "extensionsRequired", false, js, err))
        return false;
    if (!_parse_attr(val.extensionsUsed, "extensionsUsed", false, js, err))
        return false;
    if (!_parse_attr(val.images, "images", false, js, err)) return false;
    if (!_parse_attr(val.materials, "materials", false, js, err)) return false;
    if (!_parse_attr(val.meshes, "meshes", false, js, err)) return false;
    if (!_parse_attr(val.nodes, "nodes", false, js, err)) return false;
    if (!_parse_attr(val.samplers, "samplers", false, js, err)) return false;
    if (!_parse_attr(val.scene, "scene", false, js, err)) return false;
    if (!_parse_attr(val.scenes, "scenes", false, js, err)) return false;
    if (!_parse_attr(val.skins, "skins", false, js, err)) return false;
    if (!_parse_attr(val.textures, "textures", false, js, err)) return false;
    if (!_parse_end_obj(js, err)) return false;
    return true;
}

// Converts a glTF_t object to JSON
static inline bool _dump(const glTF_t& val, json& js, _parse_stack& err) {
    static const auto defval = glTF_t();
    if (!_dump((glTFProperty_t&)val, js, err)) return false;
    if (!_dump_begin_obj(js, err)) return false;
    if (!_dump_attr(
            val.accessors, "accessors", defval.accessors, false, js, err))
        return false;
    if (!_dump_attr(
            val.animations, "animations", defval.animations, false, js, err))
        return false;
    if (!_dump_attr(val.asset, "asset", defval.asset, true, js, err))
        return false;
    if (!_dump_attr(
            val.bufferViews, "bufferViews", defval.bufferViews, false, js, err))
        return false;
    if (!_dump_attr(val.buffers, "buffers", defval.buffers, false, js, err))
        return false;
    if (!_dump_attr(val.cameras, "cameras", defval.cameras, false, js, err))
        return false;
    if (!_dump_attr(val.extensionsRequired, "extensionsRequired",
            defval.extensionsRequired, false, js, err))
        return false;
    if (!_dump_attr(val.extensionsUsed, "extensionsUsed", defval.extensionsUsed,
            false, js, err))
        return false;
    if (!_dump_attr(val.images, "images", defval.images, false, js, err))
        return false;
    if (!_dump_attr(
            val.materials, "materials", defval.materials, false, js, err))
        return false;
    if (!_dump_attr(val.meshes, "meshes", defval.meshes, false, js, err))
        return false;
    if (!_dump_attr(val.nodes, "nodes", defval.nodes, false, js, err))
        return false;
    if (!_dump_attr(val.samplers, "samplers", defval.samplers, false, js, err))
        return false;
    if (!_dump_attr(val.scene, "scene", defval.scene, false, js, err))
        return false;
    if (!_dump_attr(val.scenes, "scenes", defval.scenes, false, js, err))
        return false;
    if (!_dump_attr(val.skins, "skins", defval.skins, false, js, err))
        return false;
    if (!_dump_attr(val.textures, "textures", defval.textures, false, js, err))
        return false;
    if (!_dump_end_obj(js, err)) return false;
    return true;
}

// #codegen end func

//
// Get directory name (including '/').
//
inline std::string _get_dirname(const std::string& filename) {
    auto pos = filename.rfind('/');
    if (pos == std::string::npos) pos = filename.rfind('\\');
    if (pos == std::string::npos) return "";
    return filename.substr(0, pos + 1);
}

//
// Fix path
//
inline std::string _fix_path(const std::string& path_) {
    auto path = path_;
    for (auto& c : path)
        if (c == '\\') c = '/';
    return path;
}

//
// Load a binary file in memory
// http://stackoverflow.com/questions/116038/what-is-the-best-way-to-read-an-entire-file-into-a-stdstring-in-c
//
static inline std::vector<unsigned char> _load_binfile(
    const std::string& filename, bool skip_missing) {
    std::ifstream ifs(
        filename.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
    if (!ifs) {
        if (skip_missing) return {};
        throw gltf_exception("could not open file " + filename);
    }
    std::ifstream::pos_type fileSize = ifs.tellg();
    ifs.seekg(0, std::ios::beg);
    std::vector<unsigned char> bytes(fileSize);
    ifs.read((char*)&bytes[0], fileSize);
    return bytes;
}

//
// Load a text file in memory
// http://stackoverflow.com/questions/116038/what-is-the-best-way-to-read-an-entire-file-into-a-stdstring-in-c
//
static inline std::string _load_textfile(
    const std::string& filename, bool skip_missing) {
    std::ifstream stream(filename);
    if (!stream) {
        if (skip_missing) return "";
        throw gltf_exception("could not open file " + filename);
    }
    std::stringstream sstr;
    sstr << stream.rdbuf();
    return sstr.str();
}

//
// Saves text.
//
static inline bool _save_textfile(
    const std::string& filename, const std::string& txt, std::string& errmsg) {
    auto f = fopen(filename.c_str(), "wt");
    if (!f) {
        errmsg = "cannot write file " + filename;
        return false;
    }
    fwrite(txt.c_str(), 1, (int)txt.size(), f);
    fclose(f);
    return true;
}

//
// Saves text.
//
static inline void _save_textfile(
    const std::string& filename, const std::string& txt) {
    std::string errmsg;
    auto ok = _save_textfile(filename, txt, errmsg);
    if (!ok) throw gltf_exception(errmsg);
}

//
// Saves binary.
//
static inline bool _save_binfile(const std::string& filename,
    const std::vector<unsigned char>& bin, std::string& errmsg) {
    auto f = fopen(filename.c_str(), "wb");
    if (!f) {
        errmsg = "cannot write file " + filename;
        return false;
    }
    fwrite(bin.data(), 1, (int)bin.size(), f);
    fclose(f);
    return true;
}

//
// Saves binary.
//
static inline void _save_binfile(
    const std::string& filename, const std::vector<unsigned char>& bin) {
    std::string errmsg;
    auto ok = _save_binfile(filename, bin, errmsg);
    if (!ok) throw gltf_exception(errmsg);
}
//
// Get extension name
//
static inline std::string _get_extname(const std::string& filename) {
    auto pos = filename.rfind(".");
    if (pos == std::string::npos) return "";
    return filename.substr(pos);
}

//
// Loads a gltf.
//
YGLTF_API glTF_t* load_gltf(const std::string& filename, bool load_bin,
    bool load_shader, bool load_image, bool skip_missing) {
    // clear data
    auto gltf = std::unique_ptr<glTF_t>(new glTF_t());

    // load json
    auto js = json();
    try {
        std::ifstream stream(filename.c_str());
        if (!stream) throw gltf_exception("could not load json");
        stream >> js;
    } catch (const std::exception&) {
        throw gltf_exception("could not load json");
    }

    // parse json
    auto err = _parse_stack();
    if (!_parse(*gltf.get(), js, err))
        throw(gltf_exception("error parsing gltf at " + err.pathname()));

    // load external resources
    auto dirname = _get_dirname(filename);
    if (load_bin) load_buffers(gltf.get(), dirname, skip_missing);
    if (load_shader) load_shaders(gltf.get(), dirname, skip_missing);
    if (load_image) load_images(gltf.get(), dirname, skip_missing);
    // done
    return gltf.release();
}

//
// Saves a gltf.
//
YGLTF_API void save_gltf(const std::string& filename, const glTF_t* gltf,
    bool save_bin, bool save_shader, bool save_image) {
    // dumps json
    auto js = json();
    auto err = _parse_stack();
    if (!_dump(*gltf, js, err))
        throw(gltf_exception("error dumping gltf at " + err.pathname()));

    // save json
    _save_textfile(filename, js.dump(2));

    // save external resources
    auto dirname = _get_dirname(filename);
    if (save_bin) save_buffers(gltf, dirname);
    if (save_shader) save_shaders(gltf, dirname);
    if (save_image) save_images(gltf, dirname);
}

//
// Binary gltf header
//
struct _gltf_binary_header {
    // the ASCII string `'glTF'`
    char magic[4];
    // version of the Binary glTF container format
    uint32_t version;
    // total length of the Binary glTF
    uint32_t length;
    // length, in bytes, of the glTF content
    uint32_t contentLength;
    // format of the glTF content
    uint32_t contentFormat;
};

//
// reading shortcut
//
template <typename T>
static inline void _fread(FILE* f, T* v, int count) {
    if (fread(v, sizeof(T), count, f) != count)
        throw gltf_exception("could not read binary file");
}

//
// writing shortcut
//
template <typename T>
static inline void _fwrite(FILE* f, const T* v, int count) {
    if (fwrite(v, sizeof(T), count, f) != count)
        throw gltf_exception("could not write binary file");
}

//
// Loads a binary gltf.
//
YGLTF_API glTF_t* load_binary_gltf(const std::string& filename, bool load_bin,
    bool load_shader, bool load_image, bool skip_missing) {
    // clear data
    auto gltf = std::unique_ptr<glTF_t>(new glTF_t());

    // opens binary file
    auto f = std::fopen(filename.c_str(), "rb");
    if (!f) throw gltf_exception("could not load binary file");

    // read header
    auto header = _gltf_binary_header();
    _fread(f, &header, 1);

    // validate header
    if (header.magic[0] != 'g' || header.magic[1] != 'l' ||
        header.magic[2] != 'T' || header.magic[3] != 'F')
        throw gltf_exception("invalid binary gltf format");

    // read json string
    auto json_vec = std::vector<char>(header.contentLength + 1);
    json_vec.back() = 0;
    _fread(f, &json_vec[0], header.contentLength);

    // load json
    auto js = json();
    try {
        js = json::parse(json_vec.data());
    } catch (const std::exception&) {
        throw gltf_exception("could not load json");
    }

    // parse json
    auto err = _parse_stack();
    if (!_parse(*gltf.get(), js, err))
        throw gltf_exception("error parsing gltf at " + err.pathname());

    // load external resources
    auto dirname = _get_dirname(filename);
    if (load_bin) load_buffers(gltf.get(), dirname, skip_missing);
    if (load_shader) load_shaders(gltf.get(), dirname, skip_missing);
    if (load_image) load_images(gltf.get(), dirname, skip_missing);

// load internal buffer
#if 0
    if (load_bin) {
        auto buffer = &gltf->buffers["binary_glTF"];
        buffer->data.resize(header.length - header.contentLength -
                            sizeof(header));
        if (buffer->data.size() != buffer->byteLength)
            throw gltf_exception("corrupt binary gltf");
        _fread(f, buffer->data.data(), (int)buffer->data.size());
    }
#endif

    // close
    fclose(f);

    // done
    return gltf.release();
}

//
// Saves a binary gltf.
//
YGLTF_API void save_binary_gltf(const std::string& filename, const glTF_t* gltf,
    bool save_bin, bool save_shader, bool save_image) {
    // opens binary file
    auto f = std::fopen(filename.c_str(), "wb");
    if (!f) throw gltf_exception("could not write binary file");

    // dumps json
    auto js = json();
    auto err = _parse_stack();
    if (!_dump(*gltf, js, err))
        throw gltf_exception("error parsing gltf at " + err.pathname());

    // fix string
    auto js_str = js.dump(2);
    if (js_str.length() % 4) js_str += 4 - js_str.length() % 4;

    // internal buffer
    auto buffer = &gltf->buffers.at(0);

    // prepare header
    auto header = _gltf_binary_header();
    header.magic[0] = 'g';
    header.magic[1] = 'l';
    header.magic[2] = 'T';
    header.magic[3] = 'F';
    header.version = 1;
    header.contentFormat = 0;
    header.contentLength = (int)js_str.length();
    header.length =
        sizeof(header) + (int)js_str.length() + (int)buffer->data.size();

    // write header
    _fwrite(f, &header, 1);

    // write internal buffer
    _fwrite(f, buffer->data.data(), (int)buffer->data.size());

    // close
    fclose(f);

    // save external resources
    auto dirname = _get_dirname(filename);
    if (save_bin) save_buffers(gltf, dirname);
    if (save_shader) save_shaders(gltf, dirname);
    if (save_image) save_images(gltf, dirname);
}

//
// Base 64 support
//
namespace _base64 {

static const std::string base64_chars =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz"
    "0123456789+/";

static inline bool is_base64(unsigned char c) {
    return (isalnum(c) || (c == '+') || (c == '/'));
}

static inline std::string base64_encode(
    unsigned char const* bytes_to_encode, unsigned int in_len) {
    std::string ret;
    int i = 0;
    int j = 0;
    unsigned char char_array_3[3];
    unsigned char char_array_4[4];

    while (in_len--) {
        char_array_3[i++] = *(bytes_to_encode++);
        if (i == 3) {
            char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
            char_array_4[1] = ((char_array_3[0] & 0x03) << 4) +
                              ((char_array_3[1] & 0xf0) >> 4);
            char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) +
                              ((char_array_3[2] & 0xc0) >> 6);
            char_array_4[3] = char_array_3[2] & 0x3f;

            for (i = 0; (i < 4); i++) ret += base64_chars[char_array_4[i]];
            i = 0;
        }
    }

    if (i) {
        for (j = i; j < 3; j++) char_array_3[j] = '\0';

        char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
        char_array_4[1] =
            ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
        char_array_4[2] =
            ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
        char_array_4[3] = char_array_3[2] & 0x3f;

        for (j = 0; (j < i + 1); j++) ret += base64_chars[char_array_4[j]];

        while ((i++ < 3)) ret += '=';
    }

    return ret;
}

static inline std::string base64_decode(std::string const& encoded_string) {
    int in_len = (int)encoded_string.size();
    int i = 0;
    int j = 0;
    int in_ = 0;
    unsigned char char_array_4[4], char_array_3[3];
    std::string ret;

    while (in_len-- && (encoded_string[in_] != '=') &&
           is_base64(encoded_string[in_])) {
        char_array_4[i++] = encoded_string[in_];
        in_++;
        if (i == 4) {
            for (i = 0; i < 4; i++)
                char_array_4[i] = base64_chars.find(char_array_4[i]);

            char_array_3[0] =
                (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
            char_array_3[1] = ((char_array_4[1] & 0xf) << 4) +
                              ((char_array_4[2] & 0x3c) >> 2);
            char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

            for (i = 0; (i < 3); i++) ret += char_array_3[i];
            i = 0;
        }
    }

    if (i) {
        for (j = i; j < 4; j++) char_array_4[j] = 0;

        for (j = 0; j < 4; j++)
            char_array_4[j] = base64_chars.find(char_array_4[j]);

        char_array_3[0] =
            (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
        char_array_3[1] =
            ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
        char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

        for (j = 0; (j < i - 1); j++) ret += char_array_3[j];
    }

    return ret;
}

}  // namespace _base64

//
// Check if a string starts with a prefix
//
static inline bool _startsiwith(
    const std::string& str, const std::string& prefix) {
    if (str.length() < prefix.length()) return false;
    return str.substr(0, prefix.length()) == prefix;
}

//
// Load buffer data.
//
YGLTF_API void load_buffers(
    glTF_t* gltf, const std::string& dirname, bool skip_missing) {
    for (auto& buffer_ : gltf->buffers) {
        auto buffer = &buffer_;
        if (_startsiwith(buffer->uri, "data:")) {
            // assume it is base64 and find ','
            auto pos = buffer->uri.find(',');
            if (pos == buffer->uri.npos)
                throw gltf_exception("could not decode base64 data");
            // decode
            auto data = _base64::base64_decode(buffer->uri.substr(pos + 1));
            buffer->data =
                std::vector<unsigned char>((unsigned char*)data.c_str(),
                    (unsigned char*)data.c_str() + data.length());
        } else {
            buffer->data =
                _load_binfile(_fix_path(dirname + buffer->uri), skip_missing);
        }
    }
}

//
// Load shaders data.
//
YGLTF_API void load_shaders(
    glTF_t* gltf, const std::string& dirname, bool skip_missing) {
#if 0
    for (auto&& kv : gltf->shaders) {
        auto shader = &kv.second;
        if (_startsiwith(shader->uri, "data:")) {
            // assume it is base64 and find ','
            auto pos = shader->uri.find(',');
            if (pos == shader->uri.npos)
                throw gltf_exception("could not decode base64 data");
            // decode
            shader->data = _base64::base64_decode(shader->uri.substr(pos + 1));
        } else {
            shader->data =
                _load_textfile(_fix_path(dirname + shader->uri), skip_missing);
        }
    }
#endif
}

//
// Loads images.
//
YGLTF_API void load_images(
    glTF_t* gltf, const std::string& dirname, bool skip_missing) {
#ifndef YGL_NO_STBIMAGE

    for (auto& image_ : gltf->images) {
        auto image = &image_;
        image->data = image_data_t();
        auto img = (yimg::simage*)nullptr;
        if (_startsiwith(image->uri, "data:")) {
            // assume it is base64 and find ','
            auto pos = image->uri.find(',');
            auto ext = std::string();
            if (pos == image->uri.npos)
                throw gltf_exception("could not decode base64 data");
            auto header = image->uri.substr(0, pos);
            for (auto format : {"png", "jpg", "jpeg", "tga", "ppm", "hdr"})
                if (header.find(format) != header.npos) ext = format;
            if (ext.empty())
                throw gltf_exception("unsupported embedded image format " +
                                     header.substr(0, pos));
            // decode
            auto data = _base64::base64_decode(image->uri.substr(pos + 1));
            img = yimg::load_image_from_memory(
                ext, (unsigned char*)data.c_str(), (int)data.length());
        } else {
            img = yimg::load_image(_fix_path(dirname + image->uri));
        }
        if (!img) throw gltf_exception("problem loading image");
        image->data.width = img->width;
        image->data.height = img->height;
        image->data.ncomp = img->ncomp;
        if (img->hdr)
            image->data.dataf = std::vector<float>(
                img->hdr, img->hdr + img->width * img->height * img->ncomp);
        if (img->ldr)
            image->data.datab = std::vector<unsigned char>(
                img->ldr, img->ldr + img->width * img->height * img->ncomp);
        delete img;
    }

#endif
}

//
// Save buffer data.
//
YGLTF_API void save_buffers(const glTF_t* gltf, const std::string& dirname) {
    for (auto& buffer_ : gltf->buffers) {
        auto buffer = &buffer_;
        if (_startsiwith(buffer->uri, "data:"))
            throw gltf_exception("saving of embedded data not supported");
        _save_binfile(dirname + buffer->uri, buffer->data);
    }
}

//
// Save shaders data.
//
YGLTF_API void save_shaders(const glTF_t* gltf, const std::string& dirname) {
#if 0
    for (auto&& kv : gltf->shaders) {
        auto shader = &kv.second;
        if (_startsiwith(shader->uri, "data:"))
            throw gltf_exception("saving of embedded data not supported");
        _save_textfile(dirname + shader->uri, shader->data);
    }
#endif
}

//
// Save images.
//
YGLTF_API void save_images(const glTF_t* gltf, const std::string& dirname) {
#ifndef YGL_NO_STBIMAGE

    for (auto& image_ : gltf->images) {
        auto image = &image_;
        if (_startsiwith(image->uri, "data:"))
            throw gltf_exception("saving of embedded data not supported");
        yimg::save_image(dirname + image->uri, image->data.width,
            image->data.height, image->data.ncomp, image->data.dataf.data(),
            image->data.datab.data());
    }

#endif
}

//
// Math support
//
static inline std::array<float, 16> _float4x4_mul(
    const std::array<float, 16>& a, const std::array<float, 16>& b) {
    auto c = std::array<float, 16>();
    for (auto i = 0; i < 4; i++) {
        for (auto j = 0; j < 4; j++) {
            c[j * 4 + i] = 0;
            for (auto k = 0; k < 4; k++)
                c[j * 4 + i] += a[k * 4 + i] * b[j * 4 + k];
        }
    }
    return c;
}

inline vec_array_view::vec_array_view(
    const glTF_t* gltf, const accessor_t& accessor) {
    _size = accessor.count;
    _ncomp = _num_components(accessor.type);
    _ctype = accessor.componentType;
    _normalize = accessor.normalized;
    auto buffer_view = &gltf->bufferViews.at(accessor.bufferView);
    _stride = (buffer_view->byteStride) ? buffer_view->byteStride :
                                          (_ctype_size(_ctype) * _ncomp);
    auto buffer = &gltf->buffers.at(buffer_view->buffer);
    _data = buffer->data.data() + accessor.byteOffset + buffer_view->byteOffset;
    assert(buffer_view->target == bufferView_t::target_t::array_buffer_t);
}

inline std::array<float, 4> vec_array_view::operator[](int idx) const {
    auto v = std::array<float, 4>{0, 0, 0, 1};
    for (auto i = 0; i < _ncomp; i++) {
        auto valb = _data + _stride * idx + i * _ctype_size(_ctype);
        // use double for integer conversion to attempt to maintain precision
        switch (_ctype) {
            case accessor_t::componentType_t::float_t:
                v[i] = (float)(*(float*)valb);
                break;
            case accessor_t::componentType_t::byte_t:
                if (_normalize)
                    v[i] = (float)((double)(*(char*)valb) /
                                   (double)std::numeric_limits<char>::max());
                else
                    v[i] = (float)(*(char*)valb);
                break;
            case accessor_t::componentType_t::unsigned_byte_t:
                if (_normalize)
                    v[i] =
                        (float)((double)(*(unsigned char*)valb) /
                                (double)
                                    std::numeric_limits<unsigned char>::max());
                else
                    v[i] = (float)(*(unsigned char*)valb);
                break;
            case accessor_t::componentType_t::short_t:
                if (_normalize)
                    v[i] = (float)((double)(*(short*)valb) /
                                   (double)std::numeric_limits<short>::max());
                else
                    v[i] = (float)(*(short*)valb);
                break;
            case accessor_t::componentType_t::unsigned_short_t:
                if (_normalize)
                    v[i] =
                        (float)((double)(*(unsigned short*)valb) /
                                (double)
                                    std::numeric_limits<unsigned short>::max());
                else
                    v[i] = (float)(*(unsigned short*)valb);
                break;
            case accessor_t::componentType_t::unsigned_int_t:
                if (_normalize)
                    v[i] =
                        (float)((double)(*(unsigned int*)valb) /
                                (double)
                                    std::numeric_limits<unsigned int>::max());
                else
                    v[i] = (float)(*(unsigned int*)valb);
                break;
        }
    }
    return v;
}

inline int vec_array_view::_num_components(accessor_t::type_t type) {
    switch (type) {
        case accessor_t::type_t::scalar_t: return 1;
        case accessor_t::type_t::vec2_t: return 2;
        case accessor_t::type_t::vec3_t: return 3;
        case accessor_t::type_t::vec4_t: return 4;
        default: assert(false); return 0;
    }
}

inline int vec_array_view::_ctype_size(
    accessor_t::componentType_t componentType) {
    switch (componentType) {
        case accessor_t::componentType_t::byte_t: return 1;
        case accessor_t::componentType_t::unsigned_byte_t: return 1;
        case accessor_t::componentType_t::short_t: return 2;
        case accessor_t::componentType_t::unsigned_short_t: return 2;
        case accessor_t::componentType_t::unsigned_int_t: return 4;
        case accessor_t::componentType_t::float_t: return 4;
        default: assert(false); return 0;
    }
}

//
// element_attay_view implementation
//
inline element_array_view::element_array_view(
    const glTF_t* gltf, const accessor_t& accessor) {
    _size = accessor.count;
    _ctype = accessor.componentType;
    auto buffer_view = &gltf->bufferViews.at(accessor.bufferView);
    _stride = (buffer_view->byteStride) ? buffer_view->byteStride :
                                          _ctype_size(_ctype);
    auto buffer = &gltf->buffers.at(buffer_view->buffer);
    _data = buffer->data.data() + accessor.byteOffset + buffer_view->byteOffset;
    assert(
        buffer_view->target == bufferView_t::target_t::element_array_buffer_t);
    assert(accessor.type == accessor_t::type_t::scalar_t);
}

//
// element_attay_view implementation
//
inline int element_array_view::operator[](int idx) const {
    auto valb = _data + _stride * idx;
    switch (_ctype) {
        case accessor_t::componentType_t::byte_t: return int(*(char*)valb);
        case accessor_t::componentType_t::unsigned_byte_t:
            return int(*(unsigned char*)valb);
        case accessor_t::componentType_t::short_t: return int(*(short*)valb);
        case accessor_t::componentType_t::unsigned_short_t:
            return int(*(unsigned short*)valb);
        case accessor_t::componentType_t::unsigned_int_t:
            return int(*(unsigned int*)valb);
        default: assert(false); return 0;
    }
}

//
// element_attay_view implementation
//
inline int element_array_view::_ctype_size(
    accessor_t::componentType_t componentType) {
    switch (componentType) {
        case accessor_t::componentType_t::byte_t: return 1;
        case accessor_t::componentType_t::unsigned_byte_t: return 1;
        case accessor_t::componentType_t::short_t: return 2;
        case accessor_t::componentType_t::unsigned_short_t: return 2;
        case accessor_t::componentType_t::unsigned_int_t: return 4;
        case accessor_t::componentType_t::float_t: assert(false); return 0;
        default: assert(false); return 0;
    }
}

//
// Identity matrix
//
const std::array<float, 16> _identity_float4x4 = {
    1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};

//
// Math support
//
YGLTF_API std::array<float, 16> node_transform(const node_t* node) {
    auto xf = _identity_float4x4;

    // matrix
    if (node->matrix != _identity_float4x4) {
        xf = _float4x4_mul(node->matrix, xf);
    }

    // scale
    if (node->scale != std::array<float, 3>{1, 1, 1}) {
        xf = _float4x4_mul(
            std::array<float, 16>{node->scale[0], 0, 0, 0, 0, node->scale[1], 0,
                0, 0, 0, node->scale[2], 0, 0, 0, 0, 1},
            xf);
    }

    // rotation
    if (node->rotation != std::array<float, 4>{0, 0, 0, 1}) {
        // from https://github.com/sgorsten/linalg/blob/master/linalg.h
        // public domain code
        auto q = node->rotation;
        std::array<float, 3> qx = {
            q[3] * q[3] + q[0] * q[0] - q[1] * q[1] - q[2] * q[2],
            (q[0] * q[1] + q[2] * q[3]) * 2, (q[2] * q[0] - q[1] * q[3]) * 2};
        std::array<float, 3> qy = {(q[0] * q[1] - q[2] * q[3]) * 2,
            q[3] * q[3] - q[0] * q[0] + q[1] * q[1] - q[2] * q[2],
            (q[1] * q[2] + q[0] * q[3]) * 2};
        std::array<float, 3> qz = {(q[2] * q[0] + q[1] * q[3]) * 2,
            (q[1] * q[2] - q[0] * q[3]) * 2,
            q[3] * q[3] - q[0] * q[0] - q[1] * q[1] + q[2] * q[2]};
        std::array<float, 16> r = {qx[0], qx[1], qx[2], 0, qy[0], qy[1], qy[2],
            0, qz[0], qz[1], qz[2], 0, 0, 0, 0, 1};
        xf = _float4x4_mul(r, xf);
    }

    // translation
    if (node->translation != std::array<float, 3>{0, 0, 0}) {
        xf = _float4x4_mul(std::array<float, 16>{1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                               1, 0, node->translation[0], node->translation[1],
                               node->translation[2], 1},
            xf);
    }

    return xf;
}

//
// Flattens a gltf file into a flattened asset.
//
YGLTF_API fl_gltf* flatten_gltf(const glTF_t* gltf, int scene_idx) {
    // clear asset
    auto fl_gltf = std::unique_ptr<ygltf::fl_gltf>(new ygltf::fl_gltf());

    // get scene names
    auto scenes = std::vector<int>();
    if (scene_idx < 0) {
        for (auto i = 0; i < gltf->scenes.size(); i++) scenes.push_back(i);
    } else {
        scenes.push_back(scene_idx);
    }

    // convert texture
    for (auto& tf_txt_ : gltf->textures) {
        auto tf_txt = &tf_txt_;
        auto txt = new fl_texture();
        txt->name = tf_txt->name;
        auto tf_img = &gltf->images.at(tf_txt->source);
        txt->path = (_startsiwith(tf_img->uri, "data:")) ?
                        std::string("inlines") :
                        tf_img->uri;
        txt->width = tf_img->data.width;
        txt->height = tf_img->data.height;
        txt->ncomp = tf_img->data.ncomp;
        txt->datab = tf_img->data.datab;
        txt->dataf = tf_img->data.dataf;
        fl_gltf->textures.push_back(txt);
    }

    // convert materials
    for (auto& tf_mat_ : gltf->materials) {
        auto tf_mat = &tf_mat_;
        auto mat = new fl_material();
        mat->name = tf_mat->name;
        mat->ke = tf_mat->emissiveFactor;
        mat->ke_txt = tf_mat->emissiveTexture.index;
        auto color = tf_mat->pbrMetallicRoughness.baseColorFactor;
        auto metallic = tf_mat->pbrMetallicRoughness.metallicFactor;
        mat->kd = {color[0] * (1 - metallic), color[1] * (1 - metallic),
            color[2] * (1 - metallic)};
        mat->ks = {
            color[0] * metallic, color[1] * metallic, color[2] * metallic};
        mat->rs = tf_mat->pbrMetallicRoughness.roughnessFactor;
        if (metallic < 0.5f) {
            mat->kd_txt = tf_mat->pbrMetallicRoughness.baseColorTexture.index;
            mat->ks_txt = -1;
        } else {
            mat->kd_txt = -1;
            mat->ks_txt = tf_mat->pbrMetallicRoughness.baseColorTexture.index;
        }
        fl_gltf->materials.push_back(mat);
    }

    // convert meshes
    auto meshes = std::vector<std::vector<int>>();
    for (auto tf_mesh_ : gltf->meshes) {
        auto tf_mesh = &tf_mesh_;
        meshes.push_back({});
        auto mesh = &meshes.back();
        // primitives
        for (auto& tf_primitives : tf_mesh->primitives) {
            auto prim = new fl_primitives();
            mesh->push_back((int)fl_gltf->primitives.size());
            prim->material = tf_primitives.material;
            // vertex data
            for (auto&& tf_attribute : tf_primitives.attributes) {
                auto&& semantic = tf_attribute.first;
                auto&& vals = vec_array_view(
                    gltf, gltf->accessors.at(tf_attribute.second));
                if (semantic == "POSITION") {
                    prim->pos.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++) {
                        auto v = vals[i];
                        prim->pos.push_back({v[0], v[1], v[2]});
                    }
                } else if (semantic == "NORMAL") {
                    prim->norm.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++) {
                        auto v = vals[i];
                        prim->norm.push_back({v[0], v[1], v[2]});
                    }
                } else if (semantic == "TEXCOORD" || semantic == "TEXCOORD_0") {
                    prim->texcoord.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++) {
                        auto v = vals[i];
                        prim->texcoord.push_back({v[0], v[1]});
                    }
                } else if (semantic == "COLOR" || semantic == "COLOR_0") {
                    prim->color.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++) {
                        auto v = vals[i];
                        prim->color.push_back({v[0], v[1], v[2]});
                    }
                } else {
                    // ignore
                }
            }
            // indices
            if (tf_primitives.indices < 0) {
                switch (tf_primitives.mode) {
                    case mesh_primitive_t::mode_t::triangles_t: {
                        prim->triangles.reserve(prim->pos.size() / 3);
                        for (auto i = 0; i < prim->pos.size() / 3; i++) {
                            prim->triangles.push_back(
                                {i * 3 + 0, i * 3 + 1, i * 3 + 2});
                        }
                    } break;
                    case mesh_primitive_t::mode_t::triangle_fan_t: {
                        prim->triangles.reserve(prim->pos.size() - 2);
                        for (auto i = 2; i < prim->pos.size(); i++) {
                            prim->triangles.push_back({0, i - 1, i});
                        }
                    } break;
                    case mesh_primitive_t::mode_t::triangle_strip_t: {
                        prim->triangles.reserve(prim->pos.size() - 2);
                        for (auto i = 2; i < prim->pos.size(); i++) {
                            prim->triangles.push_back({i - 2, i - 1, i});
                        }
                    } break;
                    case mesh_primitive_t::mode_t::lines_t: {
                        prim->lines.reserve(prim->pos.size() / 2);
                        for (auto i = 0; i < prim->pos.size() / 2; i++) {
                            prim->lines.push_back({i * 2 + 0, i * 2 + 1});
                        }
                    } break;
                    case mesh_primitive_t::mode_t::line_loop_t: {
                        prim->lines.reserve(prim->pos.size());
                        for (auto i = 1; i < prim->pos.size(); i++) {
                            prim->lines.push_back({i - 1, i});
                        }
                        prim->lines.back() = {(int)prim->pos.size() - 1, 0};
                    } break;
                    case mesh_primitive_t::mode_t::line_strip_t: {
                        prim->lines.reserve(prim->pos.size() - 1);
                        for (auto i = 1; i < prim->pos.size(); i++) {
                            prim->lines.push_back({i - 1, i});
                        }
                    } break;
                    case mesh_primitive_t::mode_t::points_t: {
                        prim->points.reserve(prim->pos.size());
                        for (auto i = 0; i < prim->pos.size(); i++) {
                            prim->points.push_back(i);
                        }
                    } break;
                }
            } else {
                auto indices = element_array_view(
                    gltf, gltf->accessors.at(tf_primitives.indices));
                switch (tf_primitives.mode) {
                    case mesh_primitive_t::mode_t::triangles_t: {
                        prim->triangles.reserve(indices.size());
                        for (auto i = 0; i < indices.size() / 3; i++) {
                            prim->triangles.push_back({indices[i * 3 + 0],
                                indices[i * 3 + 1], indices[i * 3 + 2]});
                        }
                    } break;
                    case mesh_primitive_t::mode_t::triangle_fan_t: {
                        prim->triangles.reserve(indices.size() - 2);
                        for (auto i = 2; i < indices.size(); i++) {
                            prim->triangles.push_back(
                                {indices[0], indices[i - 1], indices[i]});
                        }
                    } break;
                    case mesh_primitive_t::mode_t::triangle_strip_t: {
                        prim->triangles.reserve(indices.size() - 2);
                        for (auto i = 2; i < indices.size(); i++) {
                            prim->triangles.push_back(
                                {indices[i - 2], indices[i - 1], indices[i]});
                        }
                    } break;
                    case mesh_primitive_t::mode_t::lines_t: {
                        prim->lines.reserve(indices.size() / 2);
                        for (auto i = 0; i < indices.size() / 2; i++) {
                            prim->lines.push_back(
                                {indices[i * 2 + 0], indices[i * 2 + 1]});
                        }
                    } break;
                    case mesh_primitive_t::mode_t::line_loop_t: {
                        prim->lines.reserve(indices.size());
                        for (auto i = 1; i < indices.size(); i++) {
                            prim->lines.push_back({indices[i - 1], indices[i]});
                        }
                        prim->lines.back() = {
                            indices[indices.size() - 1], indices[0]};
                    } break;
                    case mesh_primitive_t::mode_t::line_strip_t: {
                        prim->lines.reserve(indices.size() - 1);
                        for (auto i = 1; i < indices.size(); i++) {
                            prim->lines.push_back({indices[i - 1], indices[i]});
                        }
                    } break;
                    case mesh_primitive_t::mode_t::points_t: {
                        prim->points.reserve(indices.size());
                        for (auto i = 0; i < indices.size(); i++) {
                            prim->points.push_back(indices[i]);
                        }
                    } break;
                }
            }
            fl_gltf->primitives.push_back(prim);
        }
    }

    // convert cameras
    auto cameras = std::vector<std::unique_ptr<fl_camera>>();
    for (auto& tf_cam_ : gltf->cameras) {
        auto tf_cam = &tf_cam_;
        auto cam = new fl_camera();
        cam->name = tf_cam->name;
        cam->ortho = tf_cam->type == camera_t::type_t::orthographic_t;
        if (cam->ortho) {
            cam->yfov = tf_cam->orthographic.ymag;
            cam->aspect = tf_cam->orthographic.xmag / tf_cam->orthographic.ymag;
        } else {
            cam->yfov = tf_cam->perspective.yfov;
            cam->aspect = tf_cam->perspective.aspectRatio;
        }
        cameras.push_back(std::unique_ptr<fl_camera>(cam));
    }

    // walk the scenes and add objects
    for (auto scn_id : scenes) {
        auto scn = &gltf->scenes.at(scn_id);
        auto fl_scn = new fl_scene();
        auto stack = std::vector<std::tuple<int, std::array<float, 16>>>();
        for (auto node_id : scn->nodes) {
            stack.push_back(std::make_tuple(node_id, _identity_float4x4));
        }
        while (!stack.empty()) {
            int node_id;
            std::array<float, 16> xf;
            std::tie(node_id, xf) = stack.back();
            stack.pop_back();
            auto node = &gltf->nodes.at(node_id);
            xf = _float4x4_mul(xf, node_transform(node));
            if (node->camera >= 0) {
                // BUG: initialization
                fl_gltf->cameras.push_back(
                    new fl_camera(*cameras.at(node->camera)));
                fl_gltf->cameras.back()->xform = xf;
                fl_scn->cameras.push_back((int)fl_gltf->cameras.size() - 1);
            }
            if (node->mesh >= 0) {
// BUG: initialization
#ifdef _WIN32
                auto fm = new fl_mesh();
                fm->name = gltf->meshes.at(node->mesh).name;
                fm->xform = xf;
                fm->primitives = meshes.at(node->mesh);
                fl_gltf->meshes.push_back(fm);
#else
                fl_gltf->meshes.push_back(
                    new fl_mesh{gltf->meshes.at(node->mesh).name, xf,
                        meshes.at(node->mesh)});
#endif
                fl_scn->meshes.push_back((int)fl_gltf->meshes.size() - 1);
            }
            for (auto child : node->children) { stack.push_back({child, xf}); }
            fl_gltf->scenes.push_back(fl_scn);
        }
    }

    return fl_gltf.release();
}

//
// Unflattnes gltf
//
YGLTF_API glTF_t* unflatten_gltf(
    const fl_gltf* fl_gltf, const std::string& buffer_uri) {
    auto gltf = std::unique_ptr<glTF_t>(new glTF_t());

    // convert cameras
    for (auto cid = 0; cid < fl_gltf->cameras.size(); cid++) {
        auto fl_cam = fl_gltf->cameras[cid];
        gltf->cameras.push_back({});
        auto cam = &gltf->cameras.back();
        cam->name = fl_cam->name;
        cam->type = (fl_cam->ortho) ? camera_t::type_t::orthographic_t :
                                      camera_t::type_t::perspective_t;
        if (fl_cam->ortho) {
            cam->orthographic.ymag = fl_cam->yfov;
            cam->orthographic.xmag = fl_cam->aspect * fl_cam->yfov;
            cam->orthographic.znear = 0.001;
            cam->orthographic.znear = 100000;
        } else {
            cam->perspective.yfov = fl_cam->yfov;
            cam->perspective.aspectRatio = fl_cam->aspect;
            cam->perspective.znear = 0.001;
            cam->perspective.znear = 100000;
        }
        gltf->nodes.push_back({});
        auto node = &gltf->nodes.back();
        node->camera = cid;
        node->matrix = fl_cam->xform;
    }

    // init buffers
    auto bid = 0;
    gltf->buffers.push_back({});
    auto buffer = &gltf->buffers.back();
    buffer->uri = buffer_uri;

    // attribute handling
    auto add_array_accessor = [&](const std::string& name, int count, int nc,
                                  const float* df) {
        gltf->bufferViews.push_back({});
        auto bufferView = &gltf->bufferViews.back();
        bufferView->buffer = bid;
        bufferView->byteOffset = (int)buffer->data.size();
        bufferView->byteStride = 0;
        bufferView->byteLength = count * nc * sizeof(float);
        buffer->data.resize(buffer->data.size() + bufferView->byteLength);
        buffer->byteLength += bufferView->byteLength;
        auto ptr =
            buffer->data.data() + buffer->data.size() - bufferView->byteLength;
        bufferView->target = bufferView_t::target_t::array_buffer_t;
        memcpy(ptr, df, bufferView->byteLength);
        gltf->accessors.push_back({});
        auto accessor = &gltf->accessors.back();
        accessor->bufferView = (int)gltf->bufferViews.size() - 1;
        accessor->byteOffset = 0;
        accessor->componentType = accessor_t::componentType_t::float_t;
        accessor->count = count;
        switch (nc) {
            case 1: accessor->type = accessor_t::type_t::scalar_t; break;
            case 2: accessor->type = accessor_t::type_t::vec2_t; break;
            case 3: accessor->type = accessor_t::type_t::vec3_t; break;
            default: assert(false);
        }
        return (int)gltf->accessors.size() - 1;
    };

    // attribute handling
    auto add_element_accessor = [&](const std::string& name, int count,
                                    const int* di, bool use_uint) {
        gltf->bufferViews.push_back({});
        auto bufferView = &gltf->bufferViews.back();
        bufferView->buffer = bid;
        bufferView->byteOffset = (int)buffer->data.size();
        bufferView->byteStride = 0;
        bufferView->byteLength = (use_uint) ? count * sizeof(unsigned int) :
                                              count * sizeof(unsigned short);
        buffer->data.resize(buffer->data.size() + bufferView->byteLength);
        buffer->byteLength += bufferView->byteLength;
        auto ptr =
            buffer->data.data() + buffer->data.size() - bufferView->byteLength;
        bufferView->target = bufferView_t::target_t::element_array_buffer_t;
        if (use_uint) {
            memcpy(ptr, di, bufferView->byteLength);
        } else {
            for (auto i = 0; i < count; i++) {
                assert(di[i] < std::numeric_limits<unsigned short>::max());
                auto s = (unsigned short)di[i];
                memcpy(ptr + i * sizeof(unsigned short), &s,
                    sizeof(unsigned short));
            }
        }
        gltf->accessors.push_back({});
        auto accessor = &gltf->accessors.back();
        accessor->bufferView = (int)gltf->bufferViews.size() - 1;
        accessor->byteOffset = 0;
        accessor->componentType =
            (use_uint) ? accessor_t::componentType_t::unsigned_int_t :
                         accessor_t::componentType_t::unsigned_short_t;
        accessor->count = count;
        accessor->type = accessor_t::type_t::scalar_t;
        return (int)gltf->accessors.size() - 1;
    };

    // convert shapes
    for (auto i = 0; i < fl_gltf->meshes.size(); i++) {
        auto fl_mesh = fl_gltf->meshes[i];
        auto gid = "mesh" + std::to_string(i);
        gltf->meshes.push_back({});
        auto mesh = &gltf->meshes.back();
        mesh->name = fl_mesh->name;
        for (auto j = 0; j < fl_mesh->primitives.size(); j++) {
            auto fl_prim = fl_gltf->primitives[fl_mesh->primitives[j]];
            auto pid = std::to_string(j);
            mesh->primitives.emplace_back();
            auto prim = &mesh->primitives.back();
            prim->material = fl_prim->material;
            if (!fl_prim->pos.empty())
                prim->attributes["POSITION"] = add_array_accessor(
                    gid + pid + "_pos", (int)fl_prim->pos.size(), 3,
                    (float*)fl_prim->pos.data());
            if (!fl_prim->norm.empty())
                prim->attributes["NORMAL"] = add_array_accessor(
                    gid + pid + "_norm", (int)fl_prim->norm.size(), 3,
                    (float*)fl_prim->norm.data());
            if (!fl_prim->texcoord.empty())
                prim->attributes["TEXCOORD_0"] = add_array_accessor(
                    gid + pid + "_texcoord", (int)fl_prim->texcoord.size(), 2,
                    (float*)fl_prim->texcoord.data());
            if (!fl_prim->color.empty())
                prim->attributes["COLOR"] = add_array_accessor(
                    gid + pid + "_color", (int)fl_prim->color.size(), 3,
                    (float*)fl_prim->color.data());
            auto elem_as_uint = fl_prim->pos.size() >
                                std::numeric_limits<unsigned short>::max();
            if (!fl_prim->points.empty()) {
                prim->indices = add_element_accessor(gid + pid + "_points",
                    (int)fl_prim->points.size(), (int*)fl_prim->points.data(),
                    elem_as_uint);
                prim->mode = mesh_primitive_t::mode_t::points_t;
            } else if (!fl_prim->lines.empty()) {
                prim->indices = add_element_accessor(gid + pid + "_lines",
                    (int)fl_prim->lines.size() * 2, (int*)fl_prim->lines.data(),
                    elem_as_uint);
                prim->mode = mesh_primitive_t::mode_t::lines_t;
            } else if (!fl_prim->triangles.empty()) {
                prim->indices = add_element_accessor(gid + pid + "_triangles",
                    (int)fl_prim->triangles.size() * 3,
                    (int*)fl_prim->triangles.data(), elem_as_uint);
                prim->mode = mesh_primitive_t::mode_t::triangles_t;
            } else
                assert(false);
        }
        gltf->nodes.push_back({});
        auto node = &gltf->nodes.back();
        node->mesh = (int)gltf->meshes.size() - 1;
        node->matrix = fl_mesh->xform;
    }

    // convert materials
    for (auto i = 0; i < fl_gltf->materials.size(); i++) {
        auto fl_mat = fl_gltf->materials[i];
        gltf->materials.push_back({});
        auto mat = &gltf->materials.back();
        mat->name = fl_mat->name;
        mat->emissiveFactor = fl_mat->ke;
        mat->emissiveTexture.index = fl_mat->ke_txt;
        mat->pbrMetallicRoughness.baseColorFactor = {
            fl_mat->kd[0] + fl_mat->ks[0], fl_mat->kd[1] + fl_mat->ks[1],
            fl_mat->kd[2] + fl_mat->ks[2], 1.0f};
        auto wall = (fl_mat->kd[0] + fl_mat->kd[1] + fl_mat->kd[2]) / 3 +
                    (fl_mat->ks[0] + fl_mat->ks[1] + fl_mat->ks[2]) / 3;
        if (wall) {
            mat->pbrMetallicRoughness.metallicFactor =
                (fl_mat->ks[0] + fl_mat->ks[1] + fl_mat->ks[2]) / wall;
        }
        if (mat->pbrMetallicRoughness.metallicFactor < 0.5f) {
            mat->pbrMetallicRoughness.baseColorTexture.index = fl_mat->kd_txt;
        } else {
            mat->pbrMetallicRoughness.baseColorTexture.index = fl_mat->ks_txt;
        }
    }

    // convert textures
    for (auto tid = 0; tid < fl_gltf->textures.size(); tid++) {
        auto&& fl_txt = fl_gltf->textures[tid];
        gltf->textures.push_back({});
        auto txt = &gltf->textures.back();
        txt->name = fl_txt->name;
        txt->type = texture_t::type_t::unsigned_byte_t;
        switch (fl_txt->ncomp) {
            case 1:
                txt->format = texture_t::format_t::luminance_t;
                txt->internalFormat = texture_t::internalFormat_t::luminance_t;
                break;
            case 2:
                txt->format = texture_t::format_t::luminance_alpha_t;
                txt->internalFormat =
                    texture_t::internalFormat_t::luminance_alpha_t;
                break;
            case 3:
                txt->format = texture_t::format_t::rgb_t;
                txt->internalFormat = texture_t::internalFormat_t::rgb_t;
                break;
            case 0:
            case 4:
                txt->format = texture_t::format_t::rgba_t;
                txt->internalFormat = texture_t::internalFormat_t::rgba_t;
                break;
            default: assert(false);
        }
        txt->source = gltf->images.size();
        txt->sampler = gltf->samplers.size();
        gltf->samplers.push_back({});
        auto sampler = &gltf->samplers.back();
        gltf->images.push_back({});
        auto image = &gltf->images.back();
        sampler->minFilter = sampler_t::minFilter_t::linear_mipmap_linear_t;
        sampler->magFilter = sampler_t::magFilter_t::linear_t;
        sampler->wrapS = sampler_t::wrapS_t::repeat_t;
        sampler->wrapT = sampler_t::wrapT_t::repeat_t;
        image->uri = fl_txt->path;
        image->data.width = fl_txt->width;
        image->data.height = fl_txt->height;
        image->data.ncomp = fl_txt->ncomp;
        image->data.datab = fl_txt->datab;
        image->data.dataf = fl_txt->dataf;
    }

    // convert scenes
    gltf->scene = 0;
    gltf->scenes.push_back({});
    auto scene = &gltf->scenes.back();
    gltf->nodes.push_back({});
    auto node = &gltf->nodes.back();
    scene->nodes.push_back((int)gltf->nodes.size() - 1);
    for (auto nid = 0; nid < gltf->nodes.size(); nid++) {
        if (nid != scene->nodes[0]) node->children.push_back(nid);
    }

    // done
    return gltf.release();
}

}  // namespace ygltf
