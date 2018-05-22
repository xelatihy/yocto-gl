//
// Implementation for Yocto/glTF.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
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
//
// LICENSE OF INCLUDED CODE FOR BASE64 (base64.h, base64.cpp)
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
//

#include "yocto_gltf.h"

#include "yocto_image.h"

#include <fstream>
#include "ext/json.hpp"

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR KHRONOS GLTF
// -----------------------------------------------------------------------------
namespace ygl {

// Json alias
using json = nlohmann::json;

// Parse int function.
void serialize(int& val, json& js, bool reading) {
    if (reading) {
        if (!js.is_number_integer())
            throw std::runtime_error("integer expected");
        val = js;
    } else {
        js = val;
    }
}

// Parse float function.
void serialize(float& val, json& js, bool reading) {
    if (reading) {
        if (!js.is_number()) throw std::runtime_error("number expected");
        val = js;
    } else {
        js = val;
    }
}

// Parse bool function.
void serialize(bool& val, json& js, bool reading) {
    if (reading) {
        if (!js.is_boolean()) throw std::runtime_error("bool expected");
        val = js;
    } else {
        js = val;
    }
}

// Parse string function.
void serialize(std::string& val, json& js, bool reading) {
    if (reading) {
        if (!js.is_string()) throw std::runtime_error("string expected");
        val = js;
    } else {
        js = val;
    }
}

// Parse json function.
void serialize(json& val, json& js, bool reading) {
    if (reading) {
        val = js;
    } else {
        js = val;
    }
}

// Parse support function.
template <typename T>
void serialize(T*& val, json& js, bool reading) {
    if (reading) {
        if (js.is_null()) {
            if (val) delete val;
            val = nullptr;
            return;
        }
        if (!js.is_object()) throw std::runtime_error("object expected");
        if (!val) val = new T();
        serialize(*val, js, reading);
    } else {
        if (!val) {
            js = nullptr;
            return;
        }
        if (!js.is_object()) js = json::object();
        serialize(*val, js, reading);
    }
}

// Parse support function.
template <typename T>
void serialize(std::vector<T>& vals, json& js, bool reading) {
    if (reading) {
        if (!js.is_array()) throw std::runtime_error("array expected");
        vals.resize(js.size());
        for (auto i = 0; i < js.size(); i++) {
            // this is contrived to support for std::vector<bool>
            auto v = T();
            serialize(v, js[i], reading);
            vals[i] = v;
        }
    } else {
        js = json::array();
        for (auto i = 0; i < vals.size(); i++)
            serialize(vals[i], js[i], reading);
    }
}

// Parse support function.
template <typename T, size_t N>
void serialize(std::array<T, N>& vals, json& js, bool reading) {
    if (reading) {
        if (!js.is_array()) throw std::runtime_error("array expected");
        if (N != js.size()) throw std::runtime_error("wrong array size");
        for (auto i = 0; i < N; i++) serialize(vals[i], js.at(i), reading);
    } else {
        js = json::array();
        for (auto i = 0; i < N; i++) serialize(vals[i], js[i], reading);
    }
}

// Parse support function.
template <typename T>
void serialize(std::map<std::string, T>& vals, json& js, bool reading) {
    if (reading) {
        if (!js.is_object()) throw std::runtime_error("object expected");
        for (auto kv = js.begin(); kv != js.end(); ++kv) {
            serialize(vals[kv.key()], kv.value(), reading);
        }
    } else {
        js = json::object();
        for (auto& kv : vals) serialize(kv.second, js[kv.first], reading);
    }
}

// Parse support function.
template <typename T, typename T1>
void serialize(T& val, json& js, bool reading,
    const std::vector<std::pair<T1, T>>& table) {
    if (reading) {
        auto v = T1();
        serialize(v, js, reading);
        auto found = false;
        for (auto& kv : table) {
            if (kv.first == v) {
                val = kv.second;
                found = true;
                break;
            }
        }
        if (!found) throw std::runtime_error("bad enum value");
    } else {
        auto found = false;
        auto v = T1();
        for (auto& kv : table) {
            if (kv.second == val) {
                v = kv.first;
                found = true;
                break;
            }
        }
        if (!found) throw std::runtime_error("invalid value");
        serialize(v, js, reading);
    }
}

// Parse support function.
template <typename T, typename T1>
void serialize(T& val, json& js, bool reading, const std::map<T, T1>& table) {
    if (reading) {
        auto v = T1();
        serialize(v, js, reading);
        auto found = false;
        for (auto& kv : table) {
            if (kv.second == v) {
                val = kv.first;
                found = true;
                break;
            }
        }
        if (!found) throw std::runtime_error("bad enum value");
    } else {
        auto found = false;
        auto v = T1();
        for (auto& kv : table) {
            if (kv.first == val) {
                v = kv.second;
                found = true;
                break;
            }
        }
        if (!found) throw std::runtime_error("invalid value");
        serialize(v, js, reading);
    }
}

// Parse support function.
void serialize(vec2f& vals, json& js, bool reading) {
    serialize((std::array<float, 2>&)vals, js, reading);
}
void serialize(vec3f& vals, json& js, bool reading) {
    serialize((std::array<float, 3>&)vals, js, reading);
}
void serialize(vec4f& vals, json& js, bool reading) {
    serialize((std::array<float, 4>&)vals, js, reading);
}
void serialize(vec2i& vals, json& js, bool reading) {
    serialize((std::array<int, 2>&)vals, js, reading);
}
void serialize(vec3i& vals, json& js, bool reading) {
    serialize((std::array<int, 3>&)vals, js, reading);
}
void serialize(vec4i& vals, json& js, bool reading) {
    serialize((std::array<int, 4>&)vals, js, reading);
}
void serialize(mat3f& vals, json& js, bool reading) {
    serialize((std::array<float, 9>&)vals, js, reading);
}
void serialize(mat4f& vals, json& js, bool reading) {
    serialize((std::array<float, 16>&)vals, js, reading);
}
void serialize(frame3f& vals, json& js, bool reading) {
    serialize((std::array<float, 12>&)vals, js, reading);
}

// Parse support function.
void serialize_obj(json& js, bool reading) {
    if (reading) {
        if (!js.is_object()) throw std::runtime_error("object expected");
    } else {
        if (!js.is_object()) js = json::object();
    }
}

// Parse support function.
template <typename T>
void serialize_attr(T& val, json& js, const char* name, bool reading,
    bool required = true, const T& def = {}) {
    if (reading) {
        if (required) {
            if (!js.count(name)) throw std::runtime_error("missing value");
            serialize(val, js.at(name), reading);
        } else {
            if (!js.count(name))
                val = def;
            else
                serialize(val, js.at(name), reading);
        }
    } else {
        if (required || val != def) serialize(val, js[name], reading);
    }
}

// Dump support function.
template <typename T>
void serialize_attr(std::vector<T>& val, json& js, const char* name,
    bool reading, bool required = true, const std::vector<T>& def = {}) {
    if (reading) {
        if (required) {
            if (!js.count(name)) throw std::runtime_error("missing value");
            serialize(val, js.at(name), reading);
        } else {
            if (!js.count(name))
                val = def;
            else
                serialize(val, js.at(name), reading);
        }
    } else {
        if (required || !val.empty()) serialize(val, js[name], reading);
    }
}

// Dump support function.
template <typename T>
void serialize_attr(std::map<std::string, T>& val, json& js, const char* name,
    bool reading, bool required = true,
    const std::map<std::string, T>& def = {}) {
    if (reading) {
        if (required) {
            if (!js.count(name)) throw std::runtime_error("missing value");
            serialize(val, js.at(name), reading);
        } else {
            if (!js.count(name))
                val = def;
            else
                serialize(val, js.at(name), reading);
        }
    } else {
        if (required || !val.empty()) serialize(val, js[name], reading);
    }
}

// #codegen begin gltf-func

// Check for default value
template <typename T>
bool operator==(const glTFid<T>& a, const glTFid<T>& b) {
    return (int)a == (int)b;
}

// Check for default value
template <typename T>
bool operator!=(const glTFid<T>& a, const glTFid<T>& b) {
    return (int)a != (int)b;
}

// Parse id function.
template <typename T>
void serialize(glTFid<T>& val, json& js, bool reading) {
    if (reading) {
        if (!js.is_number_integer()) throw std::runtime_error("int expected");
        val = glTFid<T>((int)js);
    } else {
        js = (int)val;
    }
}

// Parses a glTFProperty object
void serialize(glTFProperty& val, json& js, bool reading) {
    if (reading) {
        if (!js.is_object()) throw std::runtime_error("object expected");
#if YGL_GLTFJSON
        if (js.count("extensions"))
            serialize(val.extensions, js.at("extensions"), reading);
        if (js.count("extras")) serialize(val.extras, js.at("extras"), reading);
#endif
    } else {
        if (!js.is_object()) js = json::object();
#if YGL_GLTFJSON
        if (!val.extensions.empty())
            serialize(val.extensions, js["extensions"], reading);
        if (!val.extras.is_null()) serialize(val.extras, js["extras"], reading);
#endif
    }
}

// Parses a glTFChildOfRootProperty object
void serialize(glTFChildOfRootProperty& val, json& js, bool reading) {
    static auto def = glTFChildOfRootProperty();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(val.name, js, "name", reading, false, def.name);
}
// Parse a glTFAccessorSparseIndicesComponentType enum
void serialize(
    glTFAccessorSparseIndicesComponentType& val, json& js, bool reading) {
    static std::vector<std::pair<int, glTFAccessorSparseIndicesComponentType>>
        table = {
            {5121, glTFAccessorSparseIndicesComponentType::UnsignedByte},
            {5123, glTFAccessorSparseIndicesComponentType::UnsignedShort},
            {5125, glTFAccessorSparseIndicesComponentType::UnsignedInt},
        };
    serialize(val, js, reading, table);
}

// Parses a glTFAccessorSparseIndices object
void serialize(glTFAccessorSparseIndices& val, json& js, bool reading) {
    static auto def = glTFAccessorSparseIndices();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(
        val.bufferView, js, "bufferView", reading, true, def.bufferView);
    serialize_attr(
        val.byteOffset, js, "byteOffset", reading, false, def.byteOffset);
    serialize_attr(val.componentType, js, "componentType", reading, true,
        def.componentType);
}

// Parses a glTFAccessorSparseValues object
void serialize(glTFAccessorSparseValues& val, json& js, bool reading) {
    static auto def = glTFAccessorSparseValues();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(
        val.bufferView, js, "bufferView", reading, true, def.bufferView);
    serialize_attr(
        val.byteOffset, js, "byteOffset", reading, false, def.byteOffset);
}

// Parses a glTFAccessorSparse object
void serialize(glTFAccessorSparse& val, json& js, bool reading) {
    static auto def = glTFAccessorSparse();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(val.count, js, "count", reading, true, def.count);
    serialize_attr(val.indices, js, "indices", reading, true, def.indices);
    serialize_attr(val.values, js, "values", reading, true, def.values);
}
// Parse a glTFAccessorComponentType enum
void serialize(glTFAccessorComponentType& val, json& js, bool reading) {
    static std::vector<std::pair<int, glTFAccessorComponentType>> table = {
        {5120, glTFAccessorComponentType::Byte},
        {5121, glTFAccessorComponentType::UnsignedByte},
        {5122, glTFAccessorComponentType::Short},
        {5123, glTFAccessorComponentType::UnsignedShort},
        {5125, glTFAccessorComponentType::UnsignedInt},
        {5126, glTFAccessorComponentType::Float},
    };
    serialize(val, js, reading, table);
}

// Parse a glTFAccessorType enum
void serialize(glTFAccessorType& val, json& js, bool reading) {
    static std::vector<std::pair<std::string, glTFAccessorType>> table = {
        {"SCALAR", glTFAccessorType::Scalar},
        {"VEC2", glTFAccessorType::Vec2},
        {"VEC3", glTFAccessorType::Vec3},
        {"VEC4", glTFAccessorType::Vec4},
        {"MAT2", glTFAccessorType::Mat2},
        {"MAT3", glTFAccessorType::Mat3},
        {"MAT4", glTFAccessorType::Mat4},
    };
    serialize(val, js, reading, table);
}

// Parses a glTFAccessor object
void serialize(glTFAccessor& val, json& js, bool reading) {
    static auto def = glTFAccessor();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(
        val.bufferView, js, "bufferView", reading, false, def.bufferView);
    serialize_attr(
        val.byteOffset, js, "byteOffset", reading, false, def.byteOffset);
    serialize_attr(val.componentType, js, "componentType", reading, true,
        def.componentType);
    serialize_attr(
        val.normalized, js, "normalized", reading, false, def.normalized);
    serialize_attr(val.count, js, "count", reading, true, def.count);
    serialize_attr(val.type, js, "type", reading, true, def.type);
    serialize_attr(val.max, js, "max", reading, false, def.max);
    serialize_attr(val.min, js, "min", reading, false, def.min);
    serialize_attr(val.sparse, js, "sparse", reading, false, def.sparse);
}
// Parse a glTFAnimationChannelTargetPath enum
void serialize(glTFAnimationChannelTargetPath& val, json& js, bool reading) {
    static std::vector<std::pair<std::string, glTFAnimationChannelTargetPath>>
        table = {
            {"translation", glTFAnimationChannelTargetPath::Translation},
            {"rotation", glTFAnimationChannelTargetPath::Rotation},
            {"scale", glTFAnimationChannelTargetPath::Scale},
            {"weights", glTFAnimationChannelTargetPath::Weights},
        };
    serialize(val, js, reading, table);
}

// Parses a glTFAnimationChannelTarget object
void serialize(glTFAnimationChannelTarget& val, json& js, bool reading) {
    static auto def = glTFAnimationChannelTarget();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(val.node, js, "node", reading, true, def.node);
    serialize_attr(val.path, js, "path", reading, true, def.path);
}

// Parses a glTFAnimationChannel object
void serialize(glTFAnimationChannel& val, json& js, bool reading) {
    static auto def = glTFAnimationChannel();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(val.sampler, js, "sampler", reading, true, def.sampler);
    serialize_attr(val.target, js, "target", reading, true, def.target);
}
// Parse a glTFAnimationSamplerInterpolation enum
void serialize(glTFAnimationSamplerInterpolation& val, json& js, bool reading) {
    static std::vector<
        std::pair<std::string, glTFAnimationSamplerInterpolation>>
        table = {
            {"LINEAR", glTFAnimationSamplerInterpolation::Linear},
            {"STEP", glTFAnimationSamplerInterpolation::Step},
            {"CUBICSPLINE", glTFAnimationSamplerInterpolation::CubicSpline},
        };
    serialize(val, js, reading, table);
}

// Parses a glTFAnimationSampler object
void serialize(glTFAnimationSampler& val, json& js, bool reading) {
    static auto def = glTFAnimationSampler();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(val.input, js, "input", reading, true, def.input);
    serialize_attr(val.interpolation, js, "interpolation", reading, false,
        def.interpolation);
    serialize_attr(val.output, js, "output", reading, true, def.output);
}

// Parses a glTFAnimation object
void serialize(glTFAnimation& val, json& js, bool reading) {
    static auto def = glTFAnimation();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(val.channels, js, "channels", reading, true, def.channels);
    serialize_attr(val.samplers, js, "samplers", reading, true, def.samplers);
}

// Parses a glTFAsset object
void serialize(glTFAsset& val, json& js, bool reading) {
    static auto def = glTFAsset();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(
        val.copyright, js, "copyright", reading, false, def.copyright);
    serialize_attr(
        val.generator, js, "generator", reading, false, def.generator);
    serialize_attr(val.version, js, "version", reading, true, def.version);
    serialize_attr(
        val.minVersion, js, "minVersion", reading, false, def.minVersion);
}

// Parses a glTFBuffer object
void serialize(glTFBuffer& val, json& js, bool reading) {
    static auto def = glTFBuffer();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(val.uri, js, "uri", reading, false, def.uri);
    serialize_attr(
        val.byteLength, js, "byteLength", reading, true, def.byteLength);
}
// Parse a glTFBufferViewTarget enum
void serialize(glTFBufferViewTarget& val, json& js, bool reading) {
    static std::vector<std::pair<int, glTFBufferViewTarget>> table = {
        {34962, glTFBufferViewTarget::ArrayBuffer},
        {34963, glTFBufferViewTarget::ElementArrayBuffer},
    };
    serialize(val, js, reading, table);
}

// Parses a glTFBufferView object
void serialize(glTFBufferView& val, json& js, bool reading) {
    static auto def = glTFBufferView();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(val.buffer, js, "buffer", reading, true, def.buffer);
    serialize_attr(
        val.byteOffset, js, "byteOffset", reading, false, def.byteOffset);
    serialize_attr(
        val.byteLength, js, "byteLength", reading, true, def.byteLength);
    serialize_attr(
        val.byteStride, js, "byteStride", reading, false, def.byteStride);
    serialize_attr(val.target, js, "target", reading, false, def.target);
}

// Parses a glTFCameraOrthographic object
void serialize(glTFCameraOrthographic& val, json& js, bool reading) {
    static auto def = glTFCameraOrthographic();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(val.xmag, js, "xmag", reading, true, def.xmag);
    serialize_attr(val.ymag, js, "ymag", reading, true, def.ymag);
    serialize_attr(val.zfar, js, "zfar", reading, true, def.zfar);
    serialize_attr(val.znear, js, "znear", reading, true, def.znear);
}

// Parses a glTFCameraPerspective object
void serialize(glTFCameraPerspective& val, json& js, bool reading) {
    static auto def = glTFCameraPerspective();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(
        val.aspectRatio, js, "aspectRatio", reading, false, def.aspectRatio);
    serialize_attr(val.yfov, js, "yfov", reading, true, def.yfov);
    serialize_attr(val.zfar, js, "zfar", reading, false, def.zfar);
    serialize_attr(val.znear, js, "znear", reading, true, def.znear);
}
// Parse a glTFCameraType enum
void serialize(glTFCameraType& val, json& js, bool reading) {
    static std::vector<std::pair<std::string, glTFCameraType>> table = {
        {"perspective", glTFCameraType::Perspective},
        {"orthographic", glTFCameraType::Orthographic},
    };
    serialize(val, js, reading, table);
}

// Parses a glTFCamera object
void serialize(glTFCamera& val, json& js, bool reading) {
    static auto def = glTFCamera();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(
        val.orthographic, js, "orthographic", reading, false, def.orthographic);
    serialize_attr(
        val.perspective, js, "perspective", reading, false, def.perspective);
    serialize_attr(val.type, js, "type", reading, true, def.type);
}
// Parse a glTFImageMimeType enum
void serialize(glTFImageMimeType& val, json& js, bool reading) {
    static std::vector<std::pair<std::string, glTFImageMimeType>> table = {
        {"image/jpeg", glTFImageMimeType::ImageJpeg},
        {"image/png", glTFImageMimeType::ImagePng},
    };
    serialize(val, js, reading, table);
}

// Parses a glTFImage object
void serialize(glTFImage& val, json& js, bool reading) {
    static auto def = glTFImage();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(val.uri, js, "uri", reading, false, def.uri);
    serialize_attr(val.mimeType, js, "mimeType", reading, false, def.mimeType);
    serialize_attr(
        val.bufferView, js, "bufferView", reading, false, def.bufferView);
}

// Parses a glTFTextureInfo object
void serialize(glTFTextureInfo& val, json& js, bool reading) {
    static auto def = glTFTextureInfo();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(val.index, js, "index", reading, true, def.index);
    serialize_attr(val.texCoord, js, "texCoord", reading, false, def.texCoord);
}

// Parses a glTFTexture object
void serialize(glTFTexture& val, json& js, bool reading) {
    static auto def = glTFTexture();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(val.sampler, js, "sampler", reading, false, def.sampler);
    serialize_attr(val.source, js, "source", reading, false, def.source);
}

// Parses a glTFMaterialNormalTextureInfo object
void serialize(glTFMaterialNormalTextureInfo& val, json& js, bool reading) {
    static auto def = glTFMaterialNormalTextureInfo();
    serialize_obj(js, reading);
    serialize((glTFTextureInfo&)val, js, reading);
    serialize_attr(val.scale, js, "scale", reading, false, def.scale);
}

// Parses a glTFMaterialOcclusionTextureInfo object
void serialize(glTFMaterialOcclusionTextureInfo& val, json& js, bool reading) {
    static auto def = glTFMaterialOcclusionTextureInfo();
    serialize_obj(js, reading);
    serialize((glTFTextureInfo&)val, js, reading);
    serialize_attr(val.strength, js, "strength", reading, false, def.strength);
}

// Parses a glTFMaterialPbrMetallicRoughness object
void serialize(glTFMaterialPbrMetallicRoughness& val, json& js, bool reading) {
    static auto def = glTFMaterialPbrMetallicRoughness();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(val.baseColorFactor, js, "baseColorFactor", reading, false,
        def.baseColorFactor);
    serialize_attr(val.baseColorTexture, js, "baseColorTexture", reading, false,
        def.baseColorTexture);
    serialize_attr(val.metallicFactor, js, "metallicFactor", reading, false,
        def.metallicFactor);
    serialize_attr(val.roughnessFactor, js, "roughnessFactor", reading, false,
        def.roughnessFactor);
    serialize_attr(val.metallicRoughnessTexture, js, "metallicRoughnessTexture",
        reading, false, def.metallicRoughnessTexture);
}

// Parses a glTFMaterialPbrSpecularGlossiness object
void serialize(glTFMaterialPbrSpecularGlossiness& val, json& js, bool reading) {
    static auto def = glTFMaterialPbrSpecularGlossiness();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(val.diffuseFactor, js, "diffuseFactor", reading, false,
        def.diffuseFactor);
    serialize_attr(val.diffuseTexture, js, "diffuseTexture", reading, false,
        def.diffuseTexture);
    serialize_attr(val.specularFactor, js, "specularFactor", reading, false,
        def.specularFactor);
    serialize_attr(val.glossinessFactor, js, "glossinessFactor", reading, false,
        def.glossinessFactor);
    serialize_attr(val.specularGlossinessTexture, js,
        "specularGlossinessTexture", reading, false,
        def.specularGlossinessTexture);
}
// Parse a glTFMaterialAlphaMode enum
void serialize(glTFMaterialAlphaMode& val, json& js, bool reading) {
    static std::vector<std::pair<std::string, glTFMaterialAlphaMode>> table = {
        {"OPAQUE", glTFMaterialAlphaMode::Opaque},
        {"MASK", glTFMaterialAlphaMode::Mask},
        {"BLEND", glTFMaterialAlphaMode::Blend},
    };
    serialize(val, js, reading, table);
}

// Parses a glTFMaterial object
void serialize(glTFMaterial& val, json& js, bool reading) {
    static auto def = glTFMaterial();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(val.pbrMetallicRoughness, js, "pbrMetallicRoughness",
        reading, false, def.pbrMetallicRoughness);
    serialize_attr(val.normalTexture, js, "normalTexture", reading, false,
        def.normalTexture);
    serialize_attr(val.occlusionTexture, js, "occlusionTexture", reading, false,
        def.occlusionTexture);
    serialize_attr(val.emissiveTexture, js, "emissiveTexture", reading, false,
        def.emissiveTexture);
    serialize_attr(val.emissiveFactor, js, "emissiveFactor", reading, false,
        def.emissiveFactor);
    serialize_attr(
        val.alphaMode, js, "alphaMode", reading, false, def.alphaMode);
    serialize_attr(
        val.alphaCutoff, js, "alphaCutoff", reading, false, def.alphaCutoff);
    serialize_attr(
        val.doubleSided, js, "doubleSided", reading, false, def.doubleSided);
    if (reading) {
        if (js.count("extensions")) {
            auto& js_ext = js["extensions"];
            serialize_attr(val.pbrSpecularGlossiness, js_ext,
                "KHR_materials_pbrSpecularGlossiness", reading, false,
                def.pbrSpecularGlossiness);
        }
    } else {
        if (val.pbrSpecularGlossiness != nullptr) {
            auto& js_ext = js["extensions"];
            serialize_attr(val.pbrSpecularGlossiness, js_ext,
                "KHR_materials_pbrSpecularGlossiness", reading, false,
                def.pbrSpecularGlossiness);
        }
    }
}
// Parse a glTFMeshPrimitiveMode enum
void serialize(glTFMeshPrimitiveMode& val, json& js, bool reading) {
    static std::vector<std::pair<int, glTFMeshPrimitiveMode>> table = {
        {0, glTFMeshPrimitiveMode::Points},
        {1, glTFMeshPrimitiveMode::Lines},
        {2, glTFMeshPrimitiveMode::LineLoop},
        {3, glTFMeshPrimitiveMode::LineStrip},
        {4, glTFMeshPrimitiveMode::Triangles},
        {5, glTFMeshPrimitiveMode::TriangleStrip},
        {6, glTFMeshPrimitiveMode::TriangleFan},
    };
    serialize(val, js, reading, table);
}

// Parses a glTFMeshPrimitive object
void serialize(glTFMeshPrimitive& val, json& js, bool reading) {
    static auto def = glTFMeshPrimitive();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(
        val.attributes, js, "attributes", reading, true, def.attributes);
    serialize_attr(val.indices, js, "indices", reading, false, def.indices);
    serialize_attr(val.material, js, "material", reading, false, def.material);
    serialize_attr(val.mode, js, "mode", reading, false, def.mode);
    serialize_attr(val.targets, js, "targets", reading, false, def.targets);
}

// Parses a glTFMesh object
void serialize(glTFMesh& val, json& js, bool reading) {
    static auto def = glTFMesh();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(
        val.primitives, js, "primitives", reading, true, def.primitives);
    serialize_attr(val.weights, js, "weights", reading, false, def.weights);
}

// Parses a glTFNode object
void serialize(glTFNode& val, json& js, bool reading) {
    static auto def = glTFNode();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(val.camera, js, "camera", reading, false, def.camera);
    serialize_attr(val.children, js, "children", reading, false, def.children);
    serialize_attr(val.skin, js, "skin", reading, false, def.skin);
    serialize_attr(val.matrix, js, "matrix", reading, false, def.matrix);
    serialize_attr(val.mesh, js, "mesh", reading, false, def.mesh);
    serialize_attr(val.rotation, js, "rotation", reading, false, def.rotation);
    serialize_attr(val.scale, js, "scale", reading, false, def.scale);
    serialize_attr(
        val.translation, js, "translation", reading, false, def.translation);
    serialize_attr(val.weights, js, "weights", reading, false, def.weights);
}
// Parse a glTFSamplerMagFilter enum
void serialize(glTFSamplerMagFilter& val, json& js, bool reading) {
    static std::vector<std::pair<int, glTFSamplerMagFilter>> table = {
        {9728, glTFSamplerMagFilter::Nearest},
        {9729, glTFSamplerMagFilter::Linear},
    };
    serialize(val, js, reading, table);
}

// Parse a glTFSamplerMinFilter enum
void serialize(glTFSamplerMinFilter& val, json& js, bool reading) {
    static std::vector<std::pair<int, glTFSamplerMinFilter>> table = {
        {9728, glTFSamplerMinFilter::Nearest},
        {9729, glTFSamplerMinFilter::Linear},
        {9984, glTFSamplerMinFilter::NearestMipmapNearest},
        {9985, glTFSamplerMinFilter::LinearMipmapNearest},
        {9986, glTFSamplerMinFilter::NearestMipmapLinear},
        {9987, glTFSamplerMinFilter::LinearMipmapLinear},
    };
    serialize(val, js, reading, table);
}

// Parse a glTFSamplerWrapS enum
void serialize(glTFSamplerWrapS& val, json& js, bool reading) {
    static std::vector<std::pair<int, glTFSamplerWrapS>> table = {
        {33071, glTFSamplerWrapS::ClampToEdge},
        {33648, glTFSamplerWrapS::MirroredRepeat},
        {10497, glTFSamplerWrapS::Repeat},
    };
    serialize(val, js, reading, table);
}

// Parse a glTFSamplerWrapT enum
void serialize(glTFSamplerWrapT& val, json& js, bool reading) {
    static std::vector<std::pair<int, glTFSamplerWrapT>> table = {
        {33071, glTFSamplerWrapT::ClampToEdge},
        {33648, glTFSamplerWrapT::MirroredRepeat},
        {10497, glTFSamplerWrapT::Repeat},
    };
    serialize(val, js, reading, table);
}

// Parses a glTFSampler object
void serialize(glTFSampler& val, json& js, bool reading) {
    static auto def = glTFSampler();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(
        val.magFilter, js, "magFilter", reading, false, def.magFilter);
    serialize_attr(
        val.minFilter, js, "minFilter", reading, false, def.minFilter);
    serialize_attr(val.wrapS, js, "wrapS", reading, false, def.wrapS);
    serialize_attr(val.wrapT, js, "wrapT", reading, false, def.wrapT);
}

// Parses a glTFScene object
void serialize(glTFScene& val, json& js, bool reading) {
    static auto def = glTFScene();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(val.nodes, js, "nodes", reading, false, def.nodes);
}

// Parses a glTFSkin object
void serialize(glTFSkin& val, json& js, bool reading) {
    static auto def = glTFSkin();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(val.inverseBindMatrices, js, "inverseBindMatrices", reading,
        false, def.inverseBindMatrices);
    serialize_attr(val.skeleton, js, "skeleton", reading, false, def.skeleton);
    serialize_attr(val.joints, js, "joints", reading, true, def.joints);
}

// Parses a glTF object
void serialize(glTF& val, json& js, bool reading) {
    static auto def = glTF();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(val.extensionsUsed, js, "extensionsUsed", reading, false,
        def.extensionsUsed);
    serialize_attr(val.extensionsRequired, js, "extensionsRequired", reading,
        false, def.extensionsRequired);
    serialize_attr(
        val.accessors, js, "accessors", reading, false, def.accessors);
    serialize_attr(
        val.animations, js, "animations", reading, false, def.animations);
    serialize_attr(val.asset, js, "asset", reading, true, def.asset);
    serialize_attr(val.buffers, js, "buffers", reading, false, def.buffers);
    serialize_attr(
        val.bufferViews, js, "bufferViews", reading, false, def.bufferViews);
    serialize_attr(val.cameras, js, "cameras", reading, false, def.cameras);
    serialize_attr(val.images, js, "images", reading, false, def.images);
    serialize_attr(
        val.materials, js, "materials", reading, false, def.materials);
    serialize_attr(val.meshes, js, "meshes", reading, false, def.meshes);
    serialize_attr(val.nodes, js, "nodes", reading, false, def.nodes);
    serialize_attr(val.samplers, js, "samplers", reading, false, def.samplers);
    serialize_attr(val.scene, js, "scene", reading, false, def.scene);
    serialize_attr(val.scenes, js, "scenes", reading, false, def.scenes);
    serialize_attr(val.skins, js, "skins", reading, false, def.skins);
    serialize_attr(val.textures, js, "textures", reading, false, def.textures);
}

// #codegen end gltf-func

// Encode in base64
std::string base64_encode(
    unsigned char const* bytes_to_encode, unsigned int in_len) {
    static const std::string base64_chars =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz"
        "0123456789+/";

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

// Decode from base64
std::string base64_decode(std::string const& encoded_string) {
    static const std::string base64_chars =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz"
        "0123456789+/";

    auto is_base64 = [](unsigned char c) -> bool {
        return (isalnum(c) || (c == '+') || (c == '/'));
    };

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

static bool startswith(const std::string& str, const std::string& substr) {
    if (str.length() < substr.length()) return false;
    for (auto i = 0; i < substr.length(); i++)
        if (str[i] != substr[i]) return false;
    return true;
}

// Load buffer data.
void load_buffers(const glTF* gltf, const std::string& dirname) {
    auto fix_path = [](const std::string& path_) {
        auto path = path_;
        for (auto& c : path)
            if (c == '\\') c = '/';
        return path;
    };
    auto load_binary = [](const std::string& filename) {
        // https://stackoverflow.com/questions/174531/easiest-way-to-get-files-contents-in-c
        auto f = fopen(filename.c_str(), "rb");
        if (!f) throw std::runtime_error("cannot read file " + filename);
        fseek(f, 0, SEEK_END);
        auto len = ftell(f);
        fseek(f, 0, SEEK_SET);
        auto buf = std::vector<unsigned char>(len);
        if (fread(buf.data(), 1, len, f) != len)
            throw std::runtime_error("cannot read file " + filename);
        fclose(f);
        return buf;
    };

    for (auto buffer : gltf->buffers) {
        if (buffer->uri == "") continue;
        if (startswith(buffer->uri, "data:")) {
            // assume it is base64 and find ','
            auto pos = buffer->uri.find(',');
            if (pos == buffer->uri.npos) {
                throw std::runtime_error("could not decode base64 data");
            }
            // decode
            auto data = base64_decode(buffer->uri.substr(pos + 1));
            buffer->data =
                std::vector<unsigned char>((unsigned char*)data.c_str(),
                    (unsigned char*)data.c_str() + data.length());
        } else {
            buffer->data = load_binary(fix_path(dirname + buffer->uri));
            if (buffer->data.empty()) {
                throw std::runtime_error("could not load binary file " +
                                         fix_path(dirname + buffer->uri));
            }
        }
        if (buffer->byteLength != buffer->data.size()) {
            throw std::runtime_error("mismatched buffer size");
        }
    }
}

// Path dirname
static inline std::string path_dirname(const std::string& filename) {
    auto pos = filename.rfind('/');
    if (pos == std::string::npos) pos = filename.rfind('\\');
    if (pos == std::string::npos) return "";
    return filename.substr(0, pos + 1);
}

// Loads a gltf.
glTF* load_gltf(const std::string& filename, bool load_bin) {
    // clear data
    auto gltf = new glTF();

    // load json
    std::ifstream stream(filename.c_str());
    if (!stream) throw std::runtime_error("could not load json " + filename);
    auto js = json();
    try {
        stream >> js;
    } catch (const std::exception& e) {
        throw std::runtime_error(
            std::string("could not load json with error ") + e.what());
    }

    // parse json
    try {
        serialize(gltf, js, true);
    } catch (const std::exception& e) {
        throw std::runtime_error("error parsing gltf " + std::string(e.what()));
    }

    // load external resources
    auto dirname = path_dirname(filename);
    if (load_bin) load_buffers(gltf, dirname);

    // done
    return gltf;
}

// Save buffer data.
void save_buffers(const glTF* gltf, const std::string& dirname) {
    auto save_binary = [](const std::string& filename,
                           const std::vector<unsigned char>& data) {
        auto f = fopen(filename.c_str(), "wb");
        if (!f) throw std::runtime_error("cannot write file " + filename);
        auto num = fwrite(data.data(), 1, data.size(), f);
        if (num != data.size())
            throw std::runtime_error("cannot write file " + filename);
        fclose(f);
    };

    for (auto buffer : gltf->buffers) {
        if (startswith(buffer->uri, "data:")) {
            throw std::runtime_error("saving of embedded data not supported");
        }
        save_binary(dirname + buffer->uri, buffer->data);
    }
}

// Saves a gltf.
void save_gltf(const std::string& filename, const glTF* gltf, bool save_bin) {
    auto save_text = [](const std::string& filename, const std::string& str) {
        auto f = fopen(filename.c_str(), "wb");
        if (!f) throw std::runtime_error("cannot write file " + filename);
        auto num = fwrite(str.c_str(), 1, str.size(), f);
        if (num != str.size())
            throw std::runtime_error("cannot write file " + filename);
        fclose(f);
    };

    // dumps json
    auto js = json();
    serialize((glTF*&)gltf, js, false);

    // save json
    save_text(filename, js.dump(2));

    // save external resources
    auto dirname = path_dirname(filename);
    if (save_bin) save_buffers(gltf, dirname);
}

// reading shortcut
template <typename T>
void gltf_fread(FILE* f, T* v, int count) {
    if (fread(v, sizeof(T), count, f) != count)
        throw std::runtime_error("could not read binary file");
}

// writing shortcut
template <typename T>
void gltf_fwrite(FILE* f, const T* v, int count) {
    if (fwrite(v, sizeof(T), count, f) != count)
        std::runtime_error("could not write binary file");
}

// Loads a binary gltf.
glTF* load_binary_gltf(const std::string& filename, bool load_bin) {
    // clear data
    auto gltf = new glTF();

    // opens binary file
    auto f = fopen(filename.c_str(), "rb");
    if (!f) throw std::runtime_error("could not load binary file " + filename);

    // read magic
    uint32_t magic;
    gltf_fread(f, &magic, 1);
    if (magic != 0x46546C67) throw std::runtime_error("corrupted glb format");

    // read version
    uint32_t version;
    gltf_fread(f, &version, 1);
    if (version != 1 && version != 2)
        throw std::runtime_error("unsupported glb version");

    // read length
    uint32_t length;
    gltf_fread(f, &length, 1);

    // data
    auto json_bytes = std::vector<char>();
    auto buffer_bytes = std::vector<unsigned char>();
    uint32_t buffer_length = 0;

    if (version == 1) {
        // read content length and format
        uint32_t json_length, json_format;
        gltf_fread(f, &json_length, 1);
        gltf_fread(f, &json_format, 1);

        // read json bytes
        json_bytes.resize(json_length);
        gltf_fread(f, json_bytes.data(), json_length);

        // read buffer bytes
        if (load_bin) {
            buffer_bytes.resize(length - json_length - 20);
            gltf_fread(f, buffer_bytes.data(), (int)buffer_bytes.size());
            buffer_length = (int)buffer_bytes.size();
        }
    }

    if (version == 2) {
        // read content length and format
        uint32_t json_length, json_format;
        gltf_fread(f, &json_length, 1);
        gltf_fread(f, &json_format, 1);
        if (json_format != 0x4E4F534A) {
            throw std::runtime_error("corrupt binary format");
            return nullptr;
        }

        // read json bytes
        json_bytes.resize(json_length);
        gltf_fread(f, json_bytes.data(), (int)json_bytes.size());

        // read content length and format
        uint32_t buffer_format;
        gltf_fread(f, &buffer_length, 1);
        gltf_fread(f, &buffer_format, 1);
        if (buffer_format != 0x004E4942)
            throw std::runtime_error("corrupt binary format");

        // read buffer bytes
        if (load_bin) {
            buffer_bytes.resize(buffer_length);
            gltf_fread(f, buffer_bytes.data(), (int)buffer_bytes.size());
        }
    }

    // load json
    auto js = json();
    try {
        json_bytes.push_back(0);
        js = json::parse(json_bytes.data());
    } catch (const std::exception& e) {
        throw std::runtime_error(
            std::string("could not load json with error ") + e.what());
    }

    // parse json
    try {
        serialize(gltf, js, true);
    } catch (const std::exception& e) {
        throw std::runtime_error(
            "cannot parse gltf json " + std::string(e.what()));
        return nullptr;
    }

    // fix internal buffer
    auto buffer = gltf->buffers.at(0);
    buffer->byteLength = buffer_length;
    if (version == 2) buffer->uri = "";
    if (load_bin) { buffer->data = buffer_bytes; }

    // load external resources
    auto dirname = path_dirname(filename);
    if (load_bin) load_buffers(gltf, dirname);

    // close
    fclose(f);

    // done
    return gltf;
}

// Saves a binary gltf.
void save_binary_gltf(
    const std::string& filename, const glTF* gltf, bool save_bin) {
    // opens binary file
    auto f = fopen(filename.c_str(), "wb");
    if (!f) throw std::runtime_error("could not write binary file");

    // dumps json
    auto js = json();
    serialize((glTF*&)gltf, js, false);

    // fix string
    auto js_str = js.dump(2);
    while (js_str.length() % 4) js_str += " ";
    uint32_t json_length = (uint32_t)js_str.size();

    // internal buffer
    auto buffer = gltf->buffers.at(0);
    uint32_t buffer_length = buffer->byteLength;
    if (buffer_length % 4) buffer_length += 4 - buffer_length % 4;

    // write header
    uint32_t magic = 0x46546C67;
    gltf_fwrite(f, &magic, 1);
    uint32_t version = 2;
    gltf_fwrite(f, &version, 1);
    uint32_t length = 12 + 8 + json_length + 8 + buffer_length;
    gltf_fwrite(f, &length, 1);

    // write json
    uint32_t json_type = 0x4E4F534A;
    gltf_fwrite(f, &json_length, 1);
    gltf_fwrite(f, &json_type, 1);
    gltf_fwrite(f, js_str.data(), (int)json_length);

    if (save_bin) {
        uint32_t buffer_type = 0x004E4942;
        gltf_fwrite(f, &buffer_length, 1);
        gltf_fwrite(f, &buffer_type, 1);
        gltf_fwrite(f, buffer->data.data(), (int)buffer->data.size());
        char pad = 0;
        for (auto i = 0; i < buffer_length - buffer->data.size(); i++)
            gltf_fwrite(f, &pad, 1);
    }

    // close
    fclose(f);

    // save external resources
    auto dirname = path_dirname(filename);
    if (save_bin) save_buffers(gltf, dirname);
}

// Load glTF texture images.
void load_gltf_textures(
    glTF* gltf, const std::string& dirname, bool skip_missing) {
#if YGL_IMAGEIO
    for (auto image : gltf->images) {
        auto filename = std::string();
        image->data = gltf_image_data();
        if (image->bufferView || startswith(image->uri, "data:")) {
            auto buffer = std::string();
            auto data = (unsigned char*)nullptr;
            auto data_size = 0;
            if (image->bufferView) {
                auto view = gltf->get(image->bufferView);
                auto buffer = gltf->get(view->buffer);
                if (!view || !buffer || view->byteStride) {
                    if (skip_missing) continue;
                    throw std::runtime_error("invalid image buffer view");
                }
                if (image->mimeType == glTFImageMimeType::ImagePng)
                    filename = "internal_data.png";
                else if (image->mimeType == glTFImageMimeType::ImageJpeg)
                    filename = "internal_data.jpg";
                else {
                    if (skip_missing) continue;
                    throw std::runtime_error("unsupported image format");
                }
                data = buffer->data.data() + view->byteOffset;
                data_size = view->byteLength;
            } else {
                // assume it is base64 and find ','
                auto pos = image->uri.find(',');
                if (pos == image->uri.npos) {
                    if (skip_missing) continue;
                    throw std::runtime_error("could not decode base64 data");
                }
                auto header = image->uri.substr(0, pos);
                for (auto format : {"png", "jpg", "jpeg", "tga", "ppm", "hdr"})
                    if (header.find(format) != header.npos)
                        filename = std::string("fake.") + format;
                if (is_hdr_filename(filename)) {
                    if (skip_missing) continue;
                    throw std::runtime_error(
                        "unsupported embedded image format " +
                        header.substr(0, pos));
                }
                // decode
                buffer = base64_decode(image->uri.substr(pos + 1));
                data_size = (int)buffer.size();
                data = (unsigned char*)buffer.data();
            }
            if (is_hdr_filename(filename)) {
                image->data.hdr = load_image4f_from_memory(
                    data, data_size, image->data.width, image->data.height);
            } else {
                image->data.ldr = load_image4b_from_memory(
                    data, data_size, image->data.width, image->data.height);
            }
        } else {
            filename = dirname + image->uri;
            for (auto& c : filename)
                if (c == '\\') c = '/';
            if (is_hdr_filename(filename)) {
                image->data.hdr = load_image4f(
                    filename, image->data.width, image->data.height);
            } else {
                image->data.ldr = load_image4b(
                    filename, image->data.width, image->data.height);
            }
        }
        if (image->data.hdr.empty() && image->data.ldr.empty()) {
            if (skip_missing) continue;
            throw std::runtime_error("cannot load image " + filename);
        }
    }
#else
    throw std::runtime_error("cannot load images");
#endif
}

// Save glTF texture images.
void save_gltf_textures(
    const glTF* gltf, const std::string& dirname, bool skip_missing) {
#if YGL_IMAGEIO
    for (auto image : gltf->images) {
        if (image->data.ldr.empty() && image->data.hdr.empty()) continue;
        if (startswith(image->uri, "data:")) {
            if (skip_missing) continue;
            throw std::runtime_error("saving of embedded data not supported");
        }
        auto filename = dirname + image->uri;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        auto ok = false;
        if (!image->data.ldr.empty()) {
            ok = save_image4b(filename, image->data.width, image->data.height,
                image->data.ldr);
        }
        if (!image->data.hdr.empty()) {
            ok = save_image4f(filename, image->data.width, image->data.height,
                image->data.hdr);
        }
        if (!ok) {
            if (skip_missing) continue;
            throw std::runtime_error("cannot save image " + filename);
        }
    }
#else
    throw std::runtime_error("cannot save images");
#endif
}

accessor_view::accessor_view(const glTF* gltf, const glTFAccessor* accessor) {
    _size = accessor->count;
    _ncomp = _num_components(accessor->type);
    _ctype = accessor->componentType;
    _normalize = accessor->normalized;
    auto buffer_view = gltf->get(accessor->bufferView);
    _stride = (buffer_view->byteStride) ? buffer_view->byteStride :
                                          (_ctype_size(_ctype) * _ncomp);
    auto buffer = gltf->get(buffer_view->buffer);
    _data =
        buffer->data.data() + accessor->byteOffset + buffer_view->byteOffset;
    auto remaining_buffer_bytes =
        buffer->data.size() - (_data - buffer->data.data());
    auto view_bytes = _size * _stride;
    _valid = remaining_buffer_bytes >= view_bytes;
    if (!_valid) throw std::runtime_error("corrupted glTF accessor view");
}

float accessor_view::get(int idx, int c) const {
    auto i = min(max(c, 0), ncomp() - 1);
    auto valb = _data + _stride * idx + i * _ctype_size(_ctype);
    // use double for integer conversion to attempt to maintain
    // precision
    if (!_normalize) {
        switch (_ctype) {
            case glTFAccessorComponentType::Float:
                return (float)(*(float*)valb);
            case glTFAccessorComponentType::Byte: return (float)(*(char*)valb);
            case glTFAccessorComponentType::UnsignedByte:
                return (float)(*(unsigned char*)valb);
            case glTFAccessorComponentType::Short:
                return (float)(*(short*)valb);
            case glTFAccessorComponentType::UnsignedShort:
                return (float)(*(unsigned short*)valb);
            case glTFAccessorComponentType::UnsignedInt:
                return (float)(*(unsigned int*)valb);
            case glTFAccessorComponentType::NotSet:
                throw std::runtime_error("bad enum value");
                break;
        }

    } else {
        switch (_ctype) {
            case glTFAccessorComponentType::Float:
                return (float)(*(float*)valb);
            case glTFAccessorComponentType::Byte:
                return (float)max((float)(c / 127.0), -1.0f);
            case glTFAccessorComponentType::UnsignedByte:
                return (float)(c / 255.0);
            case glTFAccessorComponentType::Short:
                return (float)(max((float)(c / 32767.0), -1.0f));
            case glTFAccessorComponentType::UnsignedShort:
                return (float)(c / 65535.0);
            case glTFAccessorComponentType::UnsignedInt:
                return (float)(max((float)(c / 2147483647.0), -1.0f));
            case glTFAccessorComponentType::NotSet:
                throw std::runtime_error("bad enum value");
                break;
        }
    }
    return 0;
}

int accessor_view::geti(int idx, int c) const {
    auto i = min(max(c, 0), ncomp() - 1);
    auto valb = _data + _stride * idx + i * _ctype_size(_ctype);
    // use double for integer conversion to attempt to maintain
    // precision
    switch (_ctype) {
        case glTFAccessorComponentType::Float: return (int)(*(float*)valb);
        case glTFAccessorComponentType::Byte: return (int)(*(char*)valb);
        case glTFAccessorComponentType::UnsignedByte:
            return (int)(*(unsigned char*)valb);
        case glTFAccessorComponentType::Short: return (int)(*(short*)valb);
        case glTFAccessorComponentType::UnsignedShort:
            return (int)(*(unsigned short*)valb);
        case glTFAccessorComponentType::UnsignedInt:
            return (int)(*(unsigned int*)valb);
        case glTFAccessorComponentType::NotSet:
            throw std::runtime_error("bad enum value");
            break;
    }
    return 0;
}

int accessor_view::_num_components(glTFAccessorType type) {
    switch (type) {
        case glTFAccessorType::NotSet: return 1;
        case glTFAccessorType::Scalar: return 1;
        case glTFAccessorType::Vec2: return 2;
        case glTFAccessorType::Vec3: return 3;
        case glTFAccessorType::Vec4: return 4;
        case glTFAccessorType::Mat2: return 4;
        case glTFAccessorType::Mat3: return 9;
        case glTFAccessorType::Mat4: return 16;
    }
}

int accessor_view::_ctype_size(glTFAccessorComponentType componentType) {
    switch (componentType) {
        case glTFAccessorComponentType::NotSet: return 1;
        case glTFAccessorComponentType::Byte: return 1;
        case glTFAccessorComponentType::UnsignedByte: return 1;
        case glTFAccessorComponentType::Short: return 2;
        case glTFAccessorComponentType::UnsignedShort: return 2;
        case glTFAccessorComponentType::UnsignedInt: return 4;
        case glTFAccessorComponentType::Float: return 4;
    }
}

}  // namespace ygl
