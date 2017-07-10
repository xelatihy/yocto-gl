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
// BUG: bounding boxes are incorrect
//

//
// TODO: spline animation
//

#include "yocto_gltf.h"

#include <cfloat>
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
struct parse_stack {
    std::vector<std::string> path = {"glTF"};
    std::string pathname() {
        auto p = std::string();
        for (auto n : path) p += '/' + n;
        return p;
    }
};

// Parse support function.
template <typename T>
static void parse(std::vector<T>& vals, const json& js, parse_stack& err) {
    if (!js.is_array())
        throw gltf_exception("json array expected at " + err.pathname());
    vals.resize(js.size());
    for (auto i = 0; i < js.size(); i++) {
        // this is contrived to support for vector<bool>
        auto v = T();
        parse(v, js[i], err);
        vals[i] = v;
    }
}

// Parse support function.
template <typename T, int N>
static void parse(ym::vec<T, N>& vals, const json& js, parse_stack& err) {
    if (!js.is_array())
        throw gltf_exception("json array expected at " + err.pathname());
    if (N != js.size())
        throw gltf_exception("json array expected at " + err.pathname());
    for (auto i = 0; i < N; i++) { parse(vals[i], js[i], err); }
}

// Parse support function.
template <typename T, int N>
static void parse(ym::quat<T, N>& vals, const json& js, parse_stack& err) {
    if (!js.is_array())
        throw gltf_exception("json array expected at " + err.pathname());
    if (N != js.size())
        throw gltf_exception("json array expected at " + err.pathname());
    for (auto i = 0; i < N; i++) { parse(vals[i], js[i], err); }
}

// Parse support function.
template <typename T, int N, int M>
static void parse(ym::mat<T, N, M>& vals, const json& js, parse_stack& err) {
    if (!js.is_array())
        throw gltf_exception("json array expected at " + err.pathname());
    if (N * M != js.size())
        throw gltf_exception("json array expected at " + err.pathname());
    for (auto j = 0; j < M; j++) {
        for (auto i = 0; i < N; i++) { parse(vals[j][i], js[j * N + i], err); }
    }
}

// Parse support function.
template <typename T>
static void parse(
    std::map<std::string, T>& vals, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected at " + err.pathname());
    for (auto kv = js.begin(); kv != js.end(); ++kv) {
        parse(vals[kv.key()], kv.value(), err);
    }
}

// Parse support function.
template <typename T>
static void parse_attr(
    T& val, const char* name, const json& js, parse_stack& err) {
    auto iter = js.find(name);
    if (iter == js.end()) return;
    err.path.push_back(name);
    parse(val, *iter, err);
    err.path.pop_back();
}

// Parse int function.
static void parse(int& val, const json& js, parse_stack& err) {
    if (!js.is_number_integer()) throw gltf_exception(err.pathname());
    val = js;
}

// Parse float function.
static void parse(float& val, const json& js, parse_stack& err) {
    if (!js.is_number())
        throw gltf_exception("json number expected at " + err.pathname());
    val = js;
}

// Parse bool function.
static void parse(bool& val, const json& js, parse_stack& err) {
    if (!js.is_boolean())
        throw gltf_exception("json bool expected at " + err.pathname());
    val = js;
}

// Parse std::string function.
static void parse(std::string& val, const json& js, parse_stack& err) {
    if (!js.is_string())
        throw gltf_exception("json string expected at " + err.pathname());
    val = js;
}

// Parse json function.
static void parse(json& val, const json& js, parse_stack& err) { val = js; }

// Parse id function.
template <typename T>
static void parse(glTFid<T>& val, const json& js, parse_stack& err) {
    if (!js.is_number_integer())
        throw gltf_exception(
            "json unsigned integer expected at " + err.pathname());
    val = glTFid<T>((int)js);
}
// Parses a glTFProperty object
static void parse(glTFProperty*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFProperty;
    parse_attr(val->extensions, "extensions", js, err);
    parse_attr(val->extras, "extras", js, err);
}

// Parses a glTFChildOfRootProperty object
static void parse(
    glTFChildOfRootProperty*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFChildOfRootProperty;
    parse((glTFProperty*&)val, js, err);
    parse_attr(val->name, "name", js, err);
}

// Parse a glTFAccessorSparseIndicesComponentType enum
static void parse(glTFAccessorSparseIndicesComponentType& val, const json& js,
    parse_stack& err) {
    static std::map<int, glTFAccessorSparseIndicesComponentType> table = {
        {5121, glTFAccessorSparseIndicesComponentType::UnsignedByte},
        {5123, glTFAccessorSparseIndicesComponentType::UnsignedShort},
        {5125, glTFAccessorSparseIndicesComponentType::UnsignedInt},
    };
    auto v = int();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw gltf_exception(err.pathname());
    val = table[v];
}

// Parses a glTFAccessorSparseIndices object
static void parse(
    glTFAccessorSparseIndices*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFAccessorSparseIndices;
    parse((glTFProperty*&)val, js, err);
    if (!js.count("bufferView"))
        throw gltf_exception("missing required json value");
    parse_attr(val->bufferView, "bufferView", js, err);
    parse_attr(val->byteOffset, "byteOffset", js, err);
    if (!js.count("componentType"))
        throw gltf_exception("missing required json value");
    parse_attr(val->componentType, "componentType", js, err);
}

// Parses a glTFAccessorSparseValues object
static void parse(
    glTFAccessorSparseValues*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFAccessorSparseValues;
    parse((glTFProperty*&)val, js, err);
    if (!js.count("bufferView"))
        throw gltf_exception("missing required json value");
    parse_attr(val->bufferView, "bufferView", js, err);
    parse_attr(val->byteOffset, "byteOffset", js, err);
}

// Parses a glTFAccessorSparse object
static void parse(glTFAccessorSparse*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFAccessorSparse;
    parse((glTFProperty*&)val, js, err);
    if (!js.count("count")) throw gltf_exception("missing required json value");
    parse_attr(val->count, "count", js, err);
    if (!js.count("indices"))
        throw gltf_exception("missing required json value");
    parse_attr(val->indices, "indices", js, err);
    if (!js.count("values"))
        throw gltf_exception("missing required json value");
    parse_attr(val->values, "values", js, err);
}

// Parse a glTFAccessorComponentType enum
static void parse(
    glTFAccessorComponentType& val, const json& js, parse_stack& err) {
    static std::map<int, glTFAccessorComponentType> table = {
        {5120, glTFAccessorComponentType::Byte},
        {5121, glTFAccessorComponentType::UnsignedByte},
        {5122, glTFAccessorComponentType::Short},
        {5123, glTFAccessorComponentType::UnsignedShort},
        {5125, glTFAccessorComponentType::UnsignedInt},
        {5126, glTFAccessorComponentType::Float},
    };
    auto v = int();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw gltf_exception(err.pathname());
    val = table[v];
}

// Parse a glTFAccessorType enum
static void parse(glTFAccessorType& val, const json& js, parse_stack& err) {
    static std::map<std::string, glTFAccessorType> table = {
        {"SCALAR", glTFAccessorType::Scalar}, {"VEC2", glTFAccessorType::Vec2},
        {"VEC3", glTFAccessorType::Vec3}, {"VEC4", glTFAccessorType::Vec4},
        {"MAT2", glTFAccessorType::Mat2}, {"MAT3", glTFAccessorType::Mat3},
        {"MAT4", glTFAccessorType::Mat4},
    };
    auto v = std::string();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw gltf_exception(err.pathname());
    val = table[v];
}

// Parses a glTFAccessor object
static void parse(glTFAccessor*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFAccessor;
    parse((glTFChildOfRootProperty*&)val, js, err);
    parse_attr(val->bufferView, "bufferView", js, err);
    parse_attr(val->byteOffset, "byteOffset", js, err);
    if (!js.count("componentType"))
        throw gltf_exception("missing required json value");
    parse_attr(val->componentType, "componentType", js, err);
    if (!js.count("count")) throw gltf_exception("missing required json value");
    parse_attr(val->count, "count", js, err);
    parse_attr(val->max, "max", js, err);
    parse_attr(val->min, "min", js, err);
    parse_attr(val->normalized, "normalized", js, err);
    parse_attr(val->sparse, "sparse", js, err);
    if (!js.count("type")) throw gltf_exception("missing required json value");
    parse_attr(val->type, "type", js, err);
}

// Parse a glTFAnimationChannelTargetPath enum
static void parse(
    glTFAnimationChannelTargetPath& val, const json& js, parse_stack& err) {
    static std::map<std::string, glTFAnimationChannelTargetPath> table = {
        {"translation", glTFAnimationChannelTargetPath::Translation},
        {"rotation", glTFAnimationChannelTargetPath::Rotation},
        {"scale", glTFAnimationChannelTargetPath::Scale},
        {"weights", glTFAnimationChannelTargetPath::Weights},
    };
    auto v = std::string();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw gltf_exception(err.pathname());
    val = table[v];
}

// Parses a glTFAnimationChannelTarget object
static void parse(
    glTFAnimationChannelTarget*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFAnimationChannelTarget;
    parse((glTFProperty*&)val, js, err);
    if (!js.count("node")) throw gltf_exception("missing required json value");
    parse_attr(val->node, "node", js, err);
    if (!js.count("path")) throw gltf_exception("missing required json value");
    parse_attr(val->path, "path", js, err);
}

// Parses a glTFAnimationChannel object
static void parse(
    glTFAnimationChannel*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFAnimationChannel;
    parse((glTFProperty*&)val, js, err);
    if (!js.count("sampler"))
        throw gltf_exception("missing required json value");
    parse_attr(val->sampler, "sampler", js, err);
    if (!js.count("target"))
        throw gltf_exception("missing required json value");
    parse_attr(val->target, "target", js, err);
}

// Parse a glTFAnimationSamplerInterpolation enum
static void parse(
    glTFAnimationSamplerInterpolation& val, const json& js, parse_stack& err) {
    static std::map<std::string, glTFAnimationSamplerInterpolation> table = {
        {"LINEAR", glTFAnimationSamplerInterpolation::Linear},
        {"STEP", glTFAnimationSamplerInterpolation::Step},
        {"CATMULLROMSPLINE",
            glTFAnimationSamplerInterpolation::Catmullromspline},
        {"CUBICSPLINE", glTFAnimationSamplerInterpolation::Cubicspline},
    };
    auto v = std::string();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw gltf_exception(err.pathname());
    val = table[v];
}

// Parses a glTFAnimationSampler object
static void parse(
    glTFAnimationSampler*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFAnimationSampler;
    parse((glTFProperty*&)val, js, err);
    if (!js.count("input")) throw gltf_exception("missing required json value");
    parse_attr(val->input, "input", js, err);
    parse_attr(val->interpolation, "interpolation", js, err);
    if (!js.count("output"))
        throw gltf_exception("missing required json value");
    parse_attr(val->output, "output", js, err);
}

// Parses a glTFAnimation object
static void parse(glTFAnimation*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFAnimation;
    parse((glTFChildOfRootProperty*&)val, js, err);
    if (!js.count("channels"))
        throw gltf_exception("missing required json value");
    parse_attr(val->channels, "channels", js, err);
    if (!js.count("samplers"))
        throw gltf_exception("missing required json value");
    parse_attr(val->samplers, "samplers", js, err);
}

// Parses a glTFAsset object
static void parse(glTFAsset*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFAsset;
    parse((glTFProperty*&)val, js, err);
    parse_attr(val->copyright, "copyright", js, err);
    parse_attr(val->generator, "generator", js, err);
    parse_attr(val->minVersion, "minVersion", js, err);
    if (!js.count("version"))
        throw gltf_exception("missing required json value");
    parse_attr(val->version, "version", js, err);
}

// Parses a glTFBuffer object
static void parse(glTFBuffer*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFBuffer;
    parse((glTFChildOfRootProperty*&)val, js, err);
    if (!js.count("byteLength"))
        throw gltf_exception("missing required json value");
    parse_attr(val->byteLength, "byteLength", js, err);
    parse_attr(val->uri, "uri", js, err);
}

// Parse a glTFBufferViewTarget enum
static void parse(glTFBufferViewTarget& val, const json& js, parse_stack& err) {
    static std::map<int, glTFBufferViewTarget> table = {
        {34962, glTFBufferViewTarget::ArrayBuffer},
        {34963, glTFBufferViewTarget::ElementArrayBuffer},
    };
    auto v = int();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw gltf_exception(err.pathname());
    val = table[v];
}

// Parses a glTFBufferView object
static void parse(glTFBufferView*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFBufferView;
    parse((glTFChildOfRootProperty*&)val, js, err);
    if (!js.count("buffer"))
        throw gltf_exception("missing required json value");
    parse_attr(val->buffer, "buffer", js, err);
    if (!js.count("byteLength"))
        throw gltf_exception("missing required json value");
    parse_attr(val->byteLength, "byteLength", js, err);
    parse_attr(val->byteOffset, "byteOffset", js, err);
    parse_attr(val->byteStride, "byteStride", js, err);
    parse_attr(val->target, "target", js, err);
}

// Parses a glTFCameraOrthographic object
static void parse(
    glTFCameraOrthographic*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFCameraOrthographic;
    parse((glTFProperty*&)val, js, err);
    if (!js.count("xmag")) throw gltf_exception("missing required json value");
    parse_attr(val->xmag, "xmag", js, err);
    if (!js.count("ymag")) throw gltf_exception("missing required json value");
    parse_attr(val->ymag, "ymag", js, err);
    if (!js.count("zfar")) throw gltf_exception("missing required json value");
    parse_attr(val->zfar, "zfar", js, err);
    if (!js.count("znear")) throw gltf_exception("missing required json value");
    parse_attr(val->znear, "znear", js, err);
}

// Parses a glTFCameraPerspective object
static void parse(
    glTFCameraPerspective*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFCameraPerspective;
    parse((glTFProperty*&)val, js, err);
    parse_attr(val->aspectRatio, "aspectRatio", js, err);
    if (!js.count("yfov")) throw gltf_exception("missing required json value");
    parse_attr(val->yfov, "yfov", js, err);
    parse_attr(val->zfar, "zfar", js, err);
    if (!js.count("znear")) throw gltf_exception("missing required json value");
    parse_attr(val->znear, "znear", js, err);
}

// Parse a glTFCameraType enum
static void parse(glTFCameraType& val, const json& js, parse_stack& err) {
    static std::map<std::string, glTFCameraType> table = {
        {"perspective", glTFCameraType::Perspective},
        {"orthographic", glTFCameraType::Orthographic},
    };
    auto v = std::string();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw gltf_exception(err.pathname());
    val = table[v];
}

// Parses a glTFCamera object
static void parse(glTFCamera*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFCamera;
    parse((glTFChildOfRootProperty*&)val, js, err);
    parse_attr(val->orthographic, "orthographic", js, err);
    parse_attr(val->perspective, "perspective", js, err);
    if (!js.count("type")) throw gltf_exception("missing required json value");
    parse_attr(val->type, "type", js, err);
}

// Parse a glTFImageMimeType enum
static void parse(glTFImageMimeType& val, const json& js, parse_stack& err) {
    static std::map<std::string, glTFImageMimeType> table = {
        {"image/jpeg", glTFImageMimeType::ImageJpeg},
        {"image/png", glTFImageMimeType::ImagePng},
    };
    auto v = std::string();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw gltf_exception(err.pathname());
    val = table[v];
}

// Parses a glTFImage object
static void parse(glTFImage*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFImage;
    parse((glTFChildOfRootProperty*&)val, js, err);
    parse_attr(val->bufferView, "bufferView", js, err);
    parse_attr(val->mimeType, "mimeType", js, err);
    parse_attr(val->uri, "uri", js, err);
}

// Parses a glTFTextureInfo object
static void parse(glTFTextureInfo*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFTextureInfo;
    parse((glTFProperty*&)val, js, err);
    if (!js.count("index")) throw gltf_exception("missing required json value");
    parse_attr(val->index, "index", js, err);
    parse_attr(val->texCoord, "texCoord", js, err);
}

// Parses a glTFTexture object
static void parse(glTFTexture*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFTexture;
    parse((glTFChildOfRootProperty*&)val, js, err);
    parse_attr(val->sampler, "sampler", js, err);
    parse_attr(val->source, "source", js, err);
}

// Parses a glTFMaterialNormalTextureInfo object
static void parse(
    glTFMaterialNormalTextureInfo*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFMaterialNormalTextureInfo;
    parse((glTFTextureInfo*&)val, js, err);
    parse_attr(val->scale, "scale", js, err);
}

// Parses a glTFMaterialOcclusionTextureInfo object
static void parse(
    glTFMaterialOcclusionTextureInfo*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFMaterialOcclusionTextureInfo;
    parse((glTFTextureInfo*&)val, js, err);
    parse_attr(val->strength, "strength", js, err);
}

// Parses a glTFMaterialPbrMetallicRoughness object
static void parse(
    glTFMaterialPbrMetallicRoughness*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFMaterialPbrMetallicRoughness;
    parse((glTFProperty*&)val, js, err);
    parse_attr(val->baseColorFactor, "baseColorFactor", js, err);
    parse_attr(val->baseColorTexture, "baseColorTexture", js, err);
    parse_attr(val->metallicFactor, "metallicFactor", js, err);
    parse_attr(
        val->metallicRoughnessTexture, "metallicRoughnessTexture", js, err);
    parse_attr(val->roughnessFactor, "roughnessFactor", js, err);
}

// Parses a glTFMaterialPbrSpecularGlossiness object
static void parse(
    glTFMaterialPbrSpecularGlossiness*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFMaterialPbrSpecularGlossiness;
    parse((glTFProperty*&)val, js, err);
    parse_attr(val->diffuseFactor, "diffuseFactor", js, err);
    parse_attr(val->diffuseTexture, "diffuseTexture", js, err);
    parse_attr(val->glossinessFactor, "glossinessFactor", js, err);
    parse_attr(val->specularFactor, "specularFactor", js, err);
    parse_attr(
        val->specularGlossinessTexture, "specularGlossinessTexture", js, err);
}

// Parse a glTFMaterialAlphaMode enum
static void parse(
    glTFMaterialAlphaMode& val, const json& js, parse_stack& err) {
    static std::map<std::string, glTFMaterialAlphaMode> table = {
        {"OPAQUE", glTFMaterialAlphaMode::Opaque},
        {"MASK", glTFMaterialAlphaMode::Mask},
        {"BLEND", glTFMaterialAlphaMode::Blend},
    };
    auto v = std::string();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw gltf_exception(err.pathname());
    val = table[v];
}

// Parses a glTFMaterial object
static void parse(glTFMaterial*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFMaterial;
    parse((glTFChildOfRootProperty*&)val, js, err);
    parse_attr(val->alphaCutoff, "alphaCutoff", js, err);
    parse_attr(val->alphaMode, "alphaMode", js, err);
    parse_attr(val->doubleSided, "doubleSided", js, err);
    parse_attr(val->emissiveFactor, "emissiveFactor", js, err);
    parse_attr(val->emissiveTexture, "emissiveTexture", js, err);
    parse_attr(val->normalTexture, "normalTexture", js, err);
    parse_attr(val->occlusionTexture, "occlusionTexture", js, err);
    parse_attr(val->pbrMetallicRoughness, "pbrMetallicRoughness", js, err);
    if (js.count("extensions")) {
        auto& js_ext = js["extensions"];
        parse_attr(val->pbrSpecularGlossiness,
            "KHR_materials_pbrSpecularGlossiness", js_ext, err);
    }
}

// Parse a glTFMeshPrimitiveMode enum
static void parse(
    glTFMeshPrimitiveMode& val, const json& js, parse_stack& err) {
    static std::map<int, glTFMeshPrimitiveMode> table = {
        {0, glTFMeshPrimitiveMode::Points}, {1, glTFMeshPrimitiveMode::Lines},
        {2, glTFMeshPrimitiveMode::LineLoop},
        {3, glTFMeshPrimitiveMode::LineStrip},
        {4, glTFMeshPrimitiveMode::Triangles},
        {5, glTFMeshPrimitiveMode::TriangleStrip},
        {6, glTFMeshPrimitiveMode::TriangleFan},
    };
    auto v = int();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw gltf_exception(err.pathname());
    val = table[v];
}

// Parses a glTFMeshPrimitive object
static void parse(glTFMeshPrimitive*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFMeshPrimitive;
    parse((glTFProperty*&)val, js, err);
    if (!js.count("attributes"))
        throw gltf_exception("missing required json value");
    parse_attr(val->attributes, "attributes", js, err);
    parse_attr(val->indices, "indices", js, err);
    parse_attr(val->material, "material", js, err);
    parse_attr(val->mode, "mode", js, err);
    parse_attr(val->targets, "targets", js, err);
}

// Parses a glTFMesh object
static void parse(glTFMesh*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFMesh;
    parse((glTFChildOfRootProperty*&)val, js, err);
    if (!js.count("primitives"))
        throw gltf_exception("missing required json value");
    parse_attr(val->primitives, "primitives", js, err);
    parse_attr(val->weights, "weights", js, err);
}

// Parses a glTFNode object
static void parse(glTFNode*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFNode;
    parse((glTFChildOfRootProperty*&)val, js, err);
    parse_attr(val->camera, "camera", js, err);
    parse_attr(val->children, "children", js, err);
    parse_attr(val->matrix, "matrix", js, err);
    parse_attr(val->mesh, "mesh", js, err);
    parse_attr(val->rotation, "rotation", js, err);
    parse_attr(val->scale, "scale", js, err);
    parse_attr(val->skin, "skin", js, err);
    parse_attr(val->translation, "translation", js, err);
    parse_attr(val->weights, "weights", js, err);
}

// Parse a glTFSamplerMagFilter enum
static void parse(glTFSamplerMagFilter& val, const json& js, parse_stack& err) {
    static std::map<int, glTFSamplerMagFilter> table = {
        {9728, glTFSamplerMagFilter::Nearest},
        {9729, glTFSamplerMagFilter::Linear},
    };
    auto v = int();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw gltf_exception(err.pathname());
    val = table[v];
}

// Parse a glTFSamplerMinFilter enum
static void parse(glTFSamplerMinFilter& val, const json& js, parse_stack& err) {
    static std::map<int, glTFSamplerMinFilter> table = {
        {9728, glTFSamplerMinFilter::Nearest},
        {9729, glTFSamplerMinFilter::Linear},
        {9984, glTFSamplerMinFilter::NearestMipmapNearest},
        {9985, glTFSamplerMinFilter::LinearMipmapNearest},
        {9986, glTFSamplerMinFilter::NearestMipmapLinear},
        {9987, glTFSamplerMinFilter::LinearMipmapLinear},
    };
    auto v = int();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw gltf_exception(err.pathname());
    val = table[v];
}

// Parse a glTFSamplerWrapS enum
static void parse(glTFSamplerWrapS& val, const json& js, parse_stack& err) {
    static std::map<int, glTFSamplerWrapS> table = {
        {33071, glTFSamplerWrapS::ClampToEdge},
        {33648, glTFSamplerWrapS::MirroredRepeat},
        {10497, glTFSamplerWrapS::Repeat},
    };
    auto v = int();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw gltf_exception(err.pathname());
    val = table[v];
}

// Parse a glTFSamplerWrapT enum
static void parse(glTFSamplerWrapT& val, const json& js, parse_stack& err) {
    static std::map<int, glTFSamplerWrapT> table = {
        {33071, glTFSamplerWrapT::ClampToEdge},
        {33648, glTFSamplerWrapT::MirroredRepeat},
        {10497, glTFSamplerWrapT::Repeat},
    };
    auto v = int();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw gltf_exception(err.pathname());
    val = table[v];
}

// Parses a glTFSampler object
static void parse(glTFSampler*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFSampler;
    parse((glTFChildOfRootProperty*&)val, js, err);
    parse_attr(val->magFilter, "magFilter", js, err);
    parse_attr(val->minFilter, "minFilter", js, err);
    parse_attr(val->wrapS, "wrapS", js, err);
    parse_attr(val->wrapT, "wrapT", js, err);
}

// Parses a glTFScene object
static void parse(glTFScene*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFScene;
    parse((glTFChildOfRootProperty*&)val, js, err);
    parse_attr(val->nodes, "nodes", js, err);
}

// Parses a glTFSkin object
static void parse(glTFSkin*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTFSkin;
    parse((glTFChildOfRootProperty*&)val, js, err);
    parse_attr(val->inverseBindMatrices, "inverseBindMatrices", js, err);
    if (!js.count("joints"))
        throw gltf_exception("missing required json value");
    parse_attr(val->joints, "joints", js, err);
    parse_attr(val->skeleton, "skeleton", js, err);
}

// Parses a glTF object
static void parse(glTF*& val, const json& js, parse_stack& err) {
    if (!js.is_object())
        throw gltf_exception("json object expected" + err.pathname());
    if (!val) val = new glTF;
    parse((glTFProperty*&)val, js, err);
    parse_attr(val->accessors, "accessors", js, err);
    parse_attr(val->animations, "animations", js, err);
    if (!js.count("asset")) throw gltf_exception("missing required json value");
    parse_attr(val->asset, "asset", js, err);
    parse_attr(val->bufferViews, "bufferViews", js, err);
    parse_attr(val->buffers, "buffers", js, err);
    parse_attr(val->cameras, "cameras", js, err);
    parse_attr(val->extensionsRequired, "extensionsRequired", js, err);
    parse_attr(val->extensionsUsed, "extensionsUsed", js, err);
    parse_attr(val->images, "images", js, err);
    parse_attr(val->materials, "materials", js, err);
    parse_attr(val->meshes, "meshes", js, err);
    parse_attr(val->nodes, "nodes", js, err);
    parse_attr(val->samplers, "samplers", js, err);
    parse_attr(val->scene, "scene", js, err);
    parse_attr(val->scenes, "scenes", js, err);
    parse_attr(val->skins, "skins", js, err);
    parse_attr(val->textures, "textures", js, err);
}

// Dump support function.
template <typename T>
static void dump(const std::vector<T>& vals, json& js, parse_stack& err) {
    js = json::array();
    for (auto i = 0; i < vals.size(); i++) { dump(vals[i], js[i], err); }
}

// Dump support function.
template <typename T, int N>
static void dump(const ym::vec<T, N>& vals, json& js, parse_stack& err) {
    js = json::array();
    for (auto i = 0; i < N; i++) { dump(vals[i], js[i], err); }
}

// Dump support function.
template <typename T, int N>
static void dump(const ym::quat<T, N>& vals, json& js, parse_stack& err) {
    js = json::array();
    for (auto i = 0; i < N; i++) { dump(vals[i], js[i], err); }
}

// Dump support function.
template <typename T, int N, int M>
static void dump(const ym::mat<T, N, M>& vals, json& js, parse_stack& err) {
    js = json::array();
    for (auto j = 0; j < M; j++) {
        for (auto i = 0; i < N; i++) { dump(vals[j][i], js[j * N + i], err); }
    }
}

// Dump support function.
template <typename T>
static void dump(
    const std::map<std::string, T>& vals, json& js, parse_stack& err) {
    js = json::object();
    for (auto&& kv : vals) { dump(kv.second, js[kv.first], err); }
}

// Dump support function.
template <typename T>
static void dump_attr(
    const T& val, const char* name, json& js, parse_stack& err) {
    err.path.push_back(name);
    dump(val, js[name], err);
    err.path.pop_back();
}

// Converts int to json.
static void dump(const int& val, json& js, parse_stack& err) { js = val; }

// Converts float to json.
static void dump(const float& val, json& js, parse_stack& err) { js = val; }

// Converts bool to json.
static void dump(const bool& val, json& js, parse_stack& err) { js = val; }

// Converts std::string to json.
static void dump(const std::string& val, json& js, parse_stack& err) {
    js = val;
}

// Converts json to json.
static void dump(const json& val, json& js, parse_stack& err) { js = val; }

// Converts __TYPE__ to json.
template <typename T>
static void dump(const glTFid<T>& val, json& js, parse_stack& err) {
    js = (int)val;
}
// Converts a glTFProperty object to JSON
static void dump(const glTFProperty* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    if (!val->extensions.empty())
        dump_attr(val->extensions, "extensions", js, err);
    if (!val->extras.empty()) dump_attr(val->extras, "extras", js, err);
}

// Converts a glTFChildOfRootProperty object to JSON
static void dump(
    const glTFChildOfRootProperty* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    if (!val->name.empty()) dump_attr(val->name, "name", js, err);
}

// Converts a glTFAccessorSparseIndicesComponentType enum to JSON
static void dump(const glTFAccessorSparseIndicesComponentType& val, json& js,
    parse_stack& err) {
    static std::map<glTFAccessorSparseIndicesComponentType, int> table = {
        {glTFAccessorSparseIndicesComponentType::UnsignedByte, 5121},
        {glTFAccessorSparseIndicesComponentType::UnsignedShort, 5123},
        {glTFAccessorSparseIndicesComponentType::UnsignedInt, 5125},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFAccessorSparseIndices object to JSON
static void dump(
    const glTFAccessorSparseIndices* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    dump_attr(val->bufferView, "bufferView", js, err);
    if (val->byteOffset != 0) dump_attr(val->byteOffset, "byteOffset", js, err);
    dump_attr(val->componentType, "componentType", js, err);
}

// Converts a glTFAccessorSparseValues object to JSON
static void dump(
    const glTFAccessorSparseValues* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    dump_attr(val->bufferView, "bufferView", js, err);
    if (val->byteOffset != 0) dump_attr(val->byteOffset, "byteOffset", js, err);
}

// Converts a glTFAccessorSparse object to JSON
static void dump(const glTFAccessorSparse* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    dump_attr(val->count, "count", js, err);
    dump_attr(val->indices, "indices", js, err);
    dump_attr(val->values, "values", js, err);
}

// Converts a glTFAccessorComponentType enum to JSON
static void dump(
    const glTFAccessorComponentType& val, json& js, parse_stack& err) {
    static std::map<glTFAccessorComponentType, int> table = {
        {glTFAccessorComponentType::Byte, 5120},
        {glTFAccessorComponentType::UnsignedByte, 5121},
        {glTFAccessorComponentType::Short, 5122},
        {glTFAccessorComponentType::UnsignedShort, 5123},
        {glTFAccessorComponentType::UnsignedInt, 5125},
        {glTFAccessorComponentType::Float, 5126},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFAccessorType enum to JSON
static void dump(const glTFAccessorType& val, json& js, parse_stack& err) {
    static std::map<glTFAccessorType, std::string> table = {
        {glTFAccessorType::Scalar, "SCALAR"}, {glTFAccessorType::Vec2, "VEC2"},
        {glTFAccessorType::Vec3, "VEC3"}, {glTFAccessorType::Vec4, "VEC4"},
        {glTFAccessorType::Mat2, "MAT2"}, {glTFAccessorType::Mat3, "MAT3"},
        {glTFAccessorType::Mat4, "MAT4"},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFAccessor object to JSON
static void dump(const glTFAccessor* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    if ((bool)val->bufferView)
        dump_attr(val->bufferView, "bufferView", js, err);
    if (val->byteOffset != 0) dump_attr(val->byteOffset, "byteOffset", js, err);
    dump_attr(val->componentType, "componentType", js, err);
    dump_attr(val->count, "count", js, err);
    if (!val->max.empty()) dump_attr(val->max, "max", js, err);
    if (!val->min.empty()) dump_attr(val->min, "min", js, err);
    if (val->normalized) dump_attr(val->normalized, "normalized", js, err);
    if (val->sparse) dump_attr(val->sparse, "sparse", js, err);
    dump_attr(val->type, "type", js, err);
}

// Converts a glTFAnimationChannelTargetPath enum to JSON
static void dump(
    const glTFAnimationChannelTargetPath& val, json& js, parse_stack& err) {
    static std::map<glTFAnimationChannelTargetPath, std::string> table = {
        {glTFAnimationChannelTargetPath::Translation, "translation"},
        {glTFAnimationChannelTargetPath::Rotation, "rotation"},
        {glTFAnimationChannelTargetPath::Scale, "scale"},
        {glTFAnimationChannelTargetPath::Weights, "weights"},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFAnimationChannelTarget object to JSON
static void dump(
    const glTFAnimationChannelTarget* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    dump_attr(val->node, "node", js, err);
    dump_attr(val->path, "path", js, err);
}

// Converts a glTFAnimationChannel object to JSON
static void dump(const glTFAnimationChannel* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    dump_attr(val->sampler, "sampler", js, err);
    dump_attr(val->target, "target", js, err);
}

// Converts a glTFAnimationSamplerInterpolation enum to JSON
static void dump(
    const glTFAnimationSamplerInterpolation& val, json& js, parse_stack& err) {
    static std::map<glTFAnimationSamplerInterpolation, std::string> table = {
        {glTFAnimationSamplerInterpolation::Linear, "LINEAR"},
        {glTFAnimationSamplerInterpolation::Step, "STEP"},
        {glTFAnimationSamplerInterpolation::Catmullromspline,
            "CATMULLROMSPLINE"},
        {glTFAnimationSamplerInterpolation::Cubicspline, "CUBICSPLINE"},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFAnimationSampler object to JSON
static void dump(const glTFAnimationSampler* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    dump_attr(val->input, "input", js, err);
    if (val->interpolation != glTFAnimationSamplerInterpolation::NotSet)
        dump_attr(val->interpolation, "interpolation", js, err);
    dump_attr(val->output, "output", js, err);
}

// Converts a glTFAnimation object to JSON
static void dump(const glTFAnimation* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    dump_attr(val->channels, "channels", js, err);
    dump_attr(val->samplers, "samplers", js, err);
}

// Converts a glTFAsset object to JSON
static void dump(const glTFAsset* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    if (!val->copyright.empty())
        dump_attr(val->copyright, "copyright", js, err);
    if (!val->generator.empty())
        dump_attr(val->generator, "generator", js, err);
    if (!val->minVersion.empty())
        dump_attr(val->minVersion, "minVersion", js, err);
    dump_attr(val->version, "version", js, err);
}

// Converts a glTFBuffer object to JSON
static void dump(const glTFBuffer* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    dump_attr(val->byteLength, "byteLength", js, err);
    if (!val->uri.empty()) dump_attr(val->uri, "uri", js, err);
}

// Converts a glTFBufferViewTarget enum to JSON
static void dump(const glTFBufferViewTarget& val, json& js, parse_stack& err) {
    static std::map<glTFBufferViewTarget, int> table = {
        {glTFBufferViewTarget::ArrayBuffer, 34962},
        {glTFBufferViewTarget::ElementArrayBuffer, 34963},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFBufferView object to JSON
static void dump(const glTFBufferView* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    dump_attr(val->buffer, "buffer", js, err);
    dump_attr(val->byteLength, "byteLength", js, err);
    if (val->byteOffset != 0) dump_attr(val->byteOffset, "byteOffset", js, err);
    if (val->byteStride != 0) dump_attr(val->byteStride, "byteStride", js, err);
    if (val->target != glTFBufferViewTarget::NotSet)
        dump_attr(val->target, "target", js, err);
}

// Converts a glTFCameraOrthographic object to JSON
static void dump(
    const glTFCameraOrthographic* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    dump_attr(val->xmag, "xmag", js, err);
    dump_attr(val->ymag, "ymag", js, err);
    dump_attr(val->zfar, "zfar", js, err);
    dump_attr(val->znear, "znear", js, err);
}

// Converts a glTFCameraPerspective object to JSON
static void dump(const glTFCameraPerspective* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    if (val->aspectRatio != 0)
        dump_attr(val->aspectRatio, "aspectRatio", js, err);
    dump_attr(val->yfov, "yfov", js, err);
    if (val->zfar != 0) dump_attr(val->zfar, "zfar", js, err);
    dump_attr(val->znear, "znear", js, err);
}

// Converts a glTFCameraType enum to JSON
static void dump(const glTFCameraType& val, json& js, parse_stack& err) {
    static std::map<glTFCameraType, std::string> table = {
        {glTFCameraType::Perspective, "perspective"},
        {glTFCameraType::Orthographic, "orthographic"},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFCamera object to JSON
static void dump(const glTFCamera* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    if (val->orthographic)
        dump_attr(val->orthographic, "orthographic", js, err);
    if (val->perspective) dump_attr(val->perspective, "perspective", js, err);
    dump_attr(val->type, "type", js, err);
}

// Converts a glTFImageMimeType enum to JSON
static void dump(const glTFImageMimeType& val, json& js, parse_stack& err) {
    static std::map<glTFImageMimeType, std::string> table = {
        {glTFImageMimeType::ImageJpeg, "image/jpeg"},
        {glTFImageMimeType::ImagePng, "image/png"},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFImage object to JSON
static void dump(const glTFImage* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    if ((bool)val->bufferView)
        dump_attr(val->bufferView, "bufferView", js, err);
    if (val->mimeType != glTFImageMimeType::NotSet)
        dump_attr(val->mimeType, "mimeType", js, err);
    if (!val->uri.empty()) dump_attr(val->uri, "uri", js, err);
}

// Converts a glTFTextureInfo object to JSON
static void dump(const glTFTextureInfo* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    dump_attr(val->index, "index", js, err);
    if (val->texCoord != 0) dump_attr(val->texCoord, "texCoord", js, err);
}

// Converts a glTFTexture object to JSON
static void dump(const glTFTexture* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    if ((bool)val->sampler) dump_attr(val->sampler, "sampler", js, err);
    if ((bool)val->source) dump_attr(val->source, "source", js, err);
}

// Converts a glTFMaterialNormalTextureInfo object to JSON
static void dump(
    const glTFMaterialNormalTextureInfo* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFTextureInfo*)val, js, err);
    if (val->scale != 1) dump_attr(val->scale, "scale", js, err);
}

// Converts a glTFMaterialOcclusionTextureInfo object to JSON
static void dump(
    const glTFMaterialOcclusionTextureInfo* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFTextureInfo*)val, js, err);
    if (val->strength != 1) dump_attr(val->strength, "strength", js, err);
}

// Converts a glTFMaterialPbrMetallicRoughness object to JSON
static void dump(
    const glTFMaterialPbrMetallicRoughness* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    if (val->baseColorFactor != ym::vec4f{1, 1, 1, 1})
        dump_attr(val->baseColorFactor, "baseColorFactor", js, err);
    if (val->baseColorTexture)
        dump_attr(val->baseColorTexture, "baseColorTexture", js, err);
    if (val->metallicFactor != 1)
        dump_attr(val->metallicFactor, "metallicFactor", js, err);
    if (val->metallicRoughnessTexture)
        dump_attr(
            val->metallicRoughnessTexture, "metallicRoughnessTexture", js, err);
    if (val->roughnessFactor != 1)
        dump_attr(val->roughnessFactor, "roughnessFactor", js, err);
}

// Converts a glTFMaterialPbrSpecularGlossiness object to JSON
static void dump(
    const glTFMaterialPbrSpecularGlossiness* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    if (val->diffuseFactor != ym::vec4f{1, 1, 1, 1})
        dump_attr(val->diffuseFactor, "diffuseFactor", js, err);
    if (val->diffuseTexture)
        dump_attr(val->diffuseTexture, "diffuseTexture", js, err);
    if (val->glossinessFactor != 1)
        dump_attr(val->glossinessFactor, "glossinessFactor", js, err);
    if (val->specularFactor != ym::vec3f{1, 1, 1})
        dump_attr(val->specularFactor, "specularFactor", js, err);
    if (val->specularGlossinessTexture)
        dump_attr(val->specularGlossinessTexture, "specularGlossinessTexture",
            js, err);
}

// Converts a glTFMaterialAlphaMode enum to JSON
static void dump(const glTFMaterialAlphaMode& val, json& js, parse_stack& err) {
    static std::map<glTFMaterialAlphaMode, std::string> table = {
        {glTFMaterialAlphaMode::Opaque, "OPAQUE"},
        {glTFMaterialAlphaMode::Mask, "MASK"},
        {glTFMaterialAlphaMode::Blend, "BLEND"},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFMaterial object to JSON
static void dump(const glTFMaterial* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    if (val->alphaCutoff != 0.5)
        dump_attr(val->alphaCutoff, "alphaCutoff", js, err);
    if (val->alphaMode != glTFMaterialAlphaMode::NotSet)
        dump_attr(val->alphaMode, "alphaMode", js, err);
    if (val->doubleSided) dump_attr(val->doubleSided, "doubleSided", js, err);
    if (val->emissiveFactor != ym::vec3f{0, 0, 0})
        dump_attr(val->emissiveFactor, "emissiveFactor", js, err);
    if (val->emissiveTexture)
        dump_attr(val->emissiveTexture, "emissiveTexture", js, err);
    if (val->normalTexture)
        dump_attr(val->normalTexture, "normalTexture", js, err);
    if (val->occlusionTexture)
        dump_attr(val->occlusionTexture, "occlusionTexture", js, err);
    if (val->pbrMetallicRoughness)
        dump_attr(val->pbrMetallicRoughness, "pbrMetallicRoughness", js, err);
    if (val->pbrSpecularGlossiness) {
        auto& js_ext = js["extensions"];
        if (val->pbrSpecularGlossiness)
            dump_attr(val->pbrSpecularGlossiness,
                "KHR_materials_pbrSpecularGlossiness", js_ext, err);
    }
}

// Converts a glTFMeshPrimitiveMode enum to JSON
static void dump(const glTFMeshPrimitiveMode& val, json& js, parse_stack& err) {
    static std::map<glTFMeshPrimitiveMode, int> table = {
        {glTFMeshPrimitiveMode::Points, 0}, {glTFMeshPrimitiveMode::Lines, 1},
        {glTFMeshPrimitiveMode::LineLoop, 2},
        {glTFMeshPrimitiveMode::LineStrip, 3},
        {glTFMeshPrimitiveMode::Triangles, 4},
        {glTFMeshPrimitiveMode::TriangleStrip, 5},
        {glTFMeshPrimitiveMode::TriangleFan, 6},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFMeshPrimitive object to JSON
static void dump(const glTFMeshPrimitive* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    dump_attr(val->attributes, "attributes", js, err);
    if ((bool)val->indices) dump_attr(val->indices, "indices", js, err);
    if ((bool)val->material) dump_attr(val->material, "material", js, err);
    if (val->mode != glTFMeshPrimitiveMode::NotSet)
        dump_attr(val->mode, "mode", js, err);
    if (!val->targets.empty()) dump_attr(val->targets, "targets", js, err);
}

// Converts a glTFMesh object to JSON
static void dump(const glTFMesh* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    dump_attr(val->primitives, "primitives", js, err);
    if (!val->weights.empty()) dump_attr(val->weights, "weights", js, err);
}

// Converts a glTFNode object to JSON
static void dump(const glTFNode* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    if ((bool)val->camera) dump_attr(val->camera, "camera", js, err);
    if (!val->children.empty()) dump_attr(val->children, "children", js, err);
    if (val->matrix !=
        ym::mat4f{{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}})
        dump_attr(val->matrix, "matrix", js, err);
    if ((bool)val->mesh) dump_attr(val->mesh, "mesh", js, err);
    if (val->rotation != ym::quat4f{0, 0, 0, 1})
        dump_attr(val->rotation, "rotation", js, err);
    if (val->scale != ym::vec3f{1, 1, 1})
        dump_attr(val->scale, "scale", js, err);
    if ((bool)val->skin) dump_attr(val->skin, "skin", js, err);
    if (val->translation != ym::vec3f{0, 0, 0})
        dump_attr(val->translation, "translation", js, err);
    if (!val->weights.empty()) dump_attr(val->weights, "weights", js, err);
}

// Converts a glTFSamplerMagFilter enum to JSON
static void dump(const glTFSamplerMagFilter& val, json& js, parse_stack& err) {
    static std::map<glTFSamplerMagFilter, int> table = {
        {glTFSamplerMagFilter::Nearest, 9728},
        {glTFSamplerMagFilter::Linear, 9729},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFSamplerMinFilter enum to JSON
static void dump(const glTFSamplerMinFilter& val, json& js, parse_stack& err) {
    static std::map<glTFSamplerMinFilter, int> table = {
        {glTFSamplerMinFilter::Nearest, 9728},
        {glTFSamplerMinFilter::Linear, 9729},
        {glTFSamplerMinFilter::NearestMipmapNearest, 9984},
        {glTFSamplerMinFilter::LinearMipmapNearest, 9985},
        {glTFSamplerMinFilter::NearestMipmapLinear, 9986},
        {glTFSamplerMinFilter::LinearMipmapLinear, 9987},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFSamplerWrapS enum to JSON
static void dump(const glTFSamplerWrapS& val, json& js, parse_stack& err) {
    static std::map<glTFSamplerWrapS, int> table = {
        {glTFSamplerWrapS::ClampToEdge, 33071},
        {glTFSamplerWrapS::MirroredRepeat, 33648},
        {glTFSamplerWrapS::Repeat, 10497},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFSamplerWrapT enum to JSON
static void dump(const glTFSamplerWrapT& val, json& js, parse_stack& err) {
    static std::map<glTFSamplerWrapT, int> table = {
        {glTFSamplerWrapT::ClampToEdge, 33071},
        {glTFSamplerWrapT::MirroredRepeat, 33648},
        {glTFSamplerWrapT::Repeat, 10497},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFSampler object to JSON
static void dump(const glTFSampler* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    if (val->magFilter != glTFSamplerMagFilter::NotSet)
        dump_attr(val->magFilter, "magFilter", js, err);
    if (val->minFilter != glTFSamplerMinFilter::NotSet)
        dump_attr(val->minFilter, "minFilter", js, err);
    if (val->wrapS != glTFSamplerWrapS::NotSet)
        dump_attr(val->wrapS, "wrapS", js, err);
    if (val->wrapT != glTFSamplerWrapT::NotSet)
        dump_attr(val->wrapT, "wrapT", js, err);
}

// Converts a glTFScene object to JSON
static void dump(const glTFScene* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    if (!val->nodes.empty()) dump_attr(val->nodes, "nodes", js, err);
}

// Converts a glTFSkin object to JSON
static void dump(const glTFSkin* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    if ((bool)val->inverseBindMatrices)
        dump_attr(val->inverseBindMatrices, "inverseBindMatrices", js, err);
    dump_attr(val->joints, "joints", js, err);
    if ((bool)val->skeleton) dump_attr(val->skeleton, "skeleton", js, err);
}

// Converts a glTF object to JSON
static void dump(const glTF* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    if (!val->accessors.empty())
        dump_attr(val->accessors, "accessors", js, err);
    if (!val->animations.empty())
        dump_attr(val->animations, "animations", js, err);
    dump_attr(val->asset, "asset", js, err);
    if (!val->bufferViews.empty())
        dump_attr(val->bufferViews, "bufferViews", js, err);
    if (!val->buffers.empty()) dump_attr(val->buffers, "buffers", js, err);
    if (!val->cameras.empty()) dump_attr(val->cameras, "cameras", js, err);
    if (!val->extensionsRequired.empty())
        dump_attr(val->extensionsRequired, "extensionsRequired", js, err);
    if (!val->extensionsUsed.empty())
        dump_attr(val->extensionsUsed, "extensionsUsed", js, err);
    if (!val->images.empty()) dump_attr(val->images, "images", js, err);
    if (!val->materials.empty())
        dump_attr(val->materials, "materials", js, err);
    if (!val->meshes.empty()) dump_attr(val->meshes, "meshes", js, err);
    if (!val->nodes.empty()) dump_attr(val->nodes, "nodes", js, err);
    if (!val->samplers.empty()) dump_attr(val->samplers, "samplers", js, err);
    if ((bool)val->scene) dump_attr(val->scene, "scene", js, err);
    if (!val->scenes.empty()) dump_attr(val->scenes, "scenes", js, err);
    if (!val->skins.empty()) dump_attr(val->skins, "skins", js, err);
    if (!val->textures.empty()) dump_attr(val->textures, "textures", js, err);
}

// Validate support function.
template <typename T>
static void validate(const std::vector<T>& vals, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    for (auto i = 0; i < vals.size(); i++) { validate(vals[i], err, errs); }
}

// Validate support function.
template <typename T, int N>
static void validate(const ym::vec<T, N>& vals, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {}

// Validate support function.
template <typename T, int N>
static void validate(const ym::quat<T, N>& vals, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {}

// Validate support function.
template <typename T, int N, int M>
static void validate(const ym::mat<T, N, M>& vals, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {}

// Validate support function.
template <typename T>
static void validate(const std::map<std::string, T>& vals, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    for (auto&& kv : vals) {
        if (kv.first == "")
            errs.push_back({"missing dictionary key", err.pathname()});
        validate(kv.second, err, errs);
    }
}

// Validate support function.
template <typename T>
static void validate_attr(const T& val, const char* name, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    err.path.push_back(name);
    validate(val, err, errs);
    err.path.pop_back();
}

// Validates int (placeholder).
static void validate(const int& val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {}

// Validates float (placeholder).
static void validate(const float& val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {}

// Validates bool (placeholder).
static void validate(const bool& val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {}

// Validates std::string (placeholder).
static void validate(const std::string& val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {}

// Validates json (placeholder).
static void validate(const json& val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {}

// Converts __TYPE__ to json.
template <typename T>
static void validate(const glTFid<T>& val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {}
// Validates a glTFProperty object
static void validate(const glTFProperty* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate_attr(val->extensions, "extensions", err, errs);
    validate_attr(val->extras, "extras", err, errs);
}

// Validates a glTFChildOfRootProperty object
static void validate(const glTFChildOfRootProperty* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFProperty*)val, err, errs);
    validate_attr(val->name, "name", err, errs);
}

// Validates a glTFAccessorSparseIndices object
static void validate(const glTFAccessorSparseIndices* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFProperty*)val, err, errs);
    if (!(bool)val->bufferView)
        errs.push_back({"missing requried value bufferView", err.pathname()});
    validate_attr(val->bufferView, "bufferView", err, errs);
    validate_attr(val->byteOffset, "byteOffset", err, errs);
    validate_attr(val->componentType, "componentType", err, errs);
}

// Validates a glTFAccessorSparseValues object
static void validate(const glTFAccessorSparseValues* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFProperty*)val, err, errs);
    if (!(bool)val->bufferView)
        errs.push_back({"missing requried value bufferView", err.pathname()});
    validate_attr(val->bufferView, "bufferView", err, errs);
    validate_attr(val->byteOffset, "byteOffset", err, errs);
}

// Validates a glTFAccessorSparse object
static void validate(const glTFAccessorSparse* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFProperty*)val, err, errs);
    validate_attr(val->count, "count", err, errs);
    if (!val->indices)
        errs.push_back({"missing requried value indices", err.pathname()});
    validate_attr(val->indices, "indices", err, errs);
    if (!val->values)
        errs.push_back({"missing requried value values", err.pathname()});
    validate_attr(val->values, "values", err, errs);
}

// Validates a glTFAccessor object
static void validate(const glTFAccessor* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFChildOfRootProperty*)val, err, errs);
    validate_attr(val->bufferView, "bufferView", err, errs);
    validate_attr(val->byteOffset, "byteOffset", err, errs);
    validate_attr(val->componentType, "componentType", err, errs);
    validate_attr(val->count, "count", err, errs);
    validate_attr(val->max, "max", err, errs);
    validate_attr(val->min, "min", err, errs);
    validate_attr(val->normalized, "normalized", err, errs);
    validate_attr(val->sparse, "sparse", err, errs);
    validate_attr(val->type, "type", err, errs);
}

// Validates a glTFAnimationChannelTarget object
static void validate(const glTFAnimationChannelTarget* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFProperty*)val, err, errs);
    if (!(bool)val->node)
        errs.push_back({"missing requried value node", err.pathname()});
    validate_attr(val->node, "node", err, errs);
    validate_attr(val->path, "path", err, errs);
}

// Validates a glTFAnimationChannel object
static void validate(const glTFAnimationChannel* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFProperty*)val, err, errs);
    if (!(bool)val->sampler)
        errs.push_back({"missing requried value sampler", err.pathname()});
    validate_attr(val->sampler, "sampler", err, errs);
    if (!val->target)
        errs.push_back({"missing requried value target", err.pathname()});
    validate_attr(val->target, "target", err, errs);
}

// Validates a glTFAnimationSampler object
static void validate(const glTFAnimationSampler* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFProperty*)val, err, errs);
    if (!(bool)val->input)
        errs.push_back({"missing requried value input", err.pathname()});
    validate_attr(val->input, "input", err, errs);
    validate_attr(val->interpolation, "interpolation", err, errs);
    if (!(bool)val->output)
        errs.push_back({"missing requried value output", err.pathname()});
    validate_attr(val->output, "output", err, errs);
}

// Validates a glTFAnimation object
static void validate(const glTFAnimation* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFChildOfRootProperty*)val, err, errs);
    if (val->channels.empty())
        errs.push_back({"missing requried value channels", err.pathname()});
    validate_attr(val->channels, "channels", err, errs);
    if (val->samplers.empty())
        errs.push_back({"missing requried value samplers", err.pathname()});
    validate_attr(val->samplers, "samplers", err, errs);
}

// Validates a glTFAsset object
static void validate(const glTFAsset* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFProperty*)val, err, errs);
    validate_attr(val->copyright, "copyright", err, errs);
    validate_attr(val->generator, "generator", err, errs);
    validate_attr(val->minVersion, "minVersion", err, errs);
    if (val->version.empty())
        errs.push_back({"missing requried value version", err.pathname()});
    validate_attr(val->version, "version", err, errs);
}

// Validates a glTFBuffer object
static void validate(const glTFBuffer* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFChildOfRootProperty*)val, err, errs);
    validate_attr(val->byteLength, "byteLength", err, errs);
    validate_attr(val->uri, "uri", err, errs);
}

// Validates a glTFBufferView object
static void validate(const glTFBufferView* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFChildOfRootProperty*)val, err, errs);
    if (!(bool)val->buffer)
        errs.push_back({"missing requried value buffer", err.pathname()});
    validate_attr(val->buffer, "buffer", err, errs);
    validate_attr(val->byteLength, "byteLength", err, errs);
    validate_attr(val->byteOffset, "byteOffset", err, errs);
    validate_attr(val->byteStride, "byteStride", err, errs);
    validate_attr(val->target, "target", err, errs);
}

// Validates a glTFCameraOrthographic object
static void validate(const glTFCameraOrthographic* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFProperty*)val, err, errs);
    validate_attr(val->xmag, "xmag", err, errs);
    validate_attr(val->ymag, "ymag", err, errs);
    validate_attr(val->zfar, "zfar", err, errs);
    validate_attr(val->znear, "znear", err, errs);
}

// Validates a glTFCameraPerspective object
static void validate(const glTFCameraPerspective* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFProperty*)val, err, errs);
    validate_attr(val->aspectRatio, "aspectRatio", err, errs);
    validate_attr(val->yfov, "yfov", err, errs);
    validate_attr(val->zfar, "zfar", err, errs);
    validate_attr(val->znear, "znear", err, errs);
}

// Validates a glTFCamera object
static void validate(const glTFCamera* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFChildOfRootProperty*)val, err, errs);
    validate_attr(val->orthographic, "orthographic", err, errs);
    validate_attr(val->perspective, "perspective", err, errs);
    validate_attr(val->type, "type", err, errs);
}

// Validates a glTFImage object
static void validate(const glTFImage* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFChildOfRootProperty*)val, err, errs);
    validate_attr(val->bufferView, "bufferView", err, errs);
    validate_attr(val->mimeType, "mimeType", err, errs);
    validate_attr(val->uri, "uri", err, errs);
}

// Validates a glTFTextureInfo object
static void validate(const glTFTextureInfo* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFProperty*)val, err, errs);
    if (!(bool)val->index)
        errs.push_back({"missing requried value index", err.pathname()});
    validate_attr(val->index, "index", err, errs);
    validate_attr(val->texCoord, "texCoord", err, errs);
}

// Validates a glTFTexture object
static void validate(const glTFTexture* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFChildOfRootProperty*)val, err, errs);
    validate_attr(val->sampler, "sampler", err, errs);
    validate_attr(val->source, "source", err, errs);
}

// Validates a glTFMaterialNormalTextureInfo object
static void validate(const glTFMaterialNormalTextureInfo* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFTextureInfo*)val, err, errs);
    validate_attr(val->scale, "scale", err, errs);
}

// Validates a glTFMaterialOcclusionTextureInfo object
static void validate(const glTFMaterialOcclusionTextureInfo* val,
    parse_stack& err, std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFTextureInfo*)val, err, errs);
    validate_attr(val->strength, "strength", err, errs);
}

// Validates a glTFMaterialPbrMetallicRoughness object
static void validate(const glTFMaterialPbrMetallicRoughness* val,
    parse_stack& err, std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFProperty*)val, err, errs);
    validate_attr(val->baseColorFactor, "baseColorFactor", err, errs);
    validate_attr(val->baseColorTexture, "baseColorTexture", err, errs);
    validate_attr(val->metallicFactor, "metallicFactor", err, errs);
    validate_attr(
        val->metallicRoughnessTexture, "metallicRoughnessTexture", err, errs);
    validate_attr(val->roughnessFactor, "roughnessFactor", err, errs);
}

// Validates a glTFMaterialPbrSpecularGlossiness object
static void validate(const glTFMaterialPbrSpecularGlossiness* val,
    parse_stack& err, std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFProperty*)val, err, errs);
    validate_attr(val->diffuseFactor, "diffuseFactor", err, errs);
    validate_attr(val->diffuseTexture, "diffuseTexture", err, errs);
    validate_attr(val->glossinessFactor, "glossinessFactor", err, errs);
    validate_attr(val->specularFactor, "specularFactor", err, errs);
    validate_attr(
        val->specularGlossinessTexture, "specularGlossinessTexture", err, errs);
}

// Validates a glTFMaterial object
static void validate(const glTFMaterial* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFChildOfRootProperty*)val, err, errs);
    validate_attr(val->alphaCutoff, "alphaCutoff", err, errs);
    validate_attr(val->alphaMode, "alphaMode", err, errs);
    validate_attr(val->doubleSided, "doubleSided", err, errs);
    validate_attr(val->emissiveFactor, "emissiveFactor", err, errs);
    validate_attr(val->emissiveTexture, "emissiveTexture", err, errs);
    validate_attr(val->normalTexture, "normalTexture", err, errs);
    validate_attr(val->occlusionTexture, "occlusionTexture", err, errs);
    validate_attr(val->pbrMetallicRoughness, "pbrMetallicRoughness", err, errs);
    validate_attr(
        val->pbrSpecularGlossiness, "pbrSpecularGlossiness", err, errs);
}

// Validates a glTFMeshPrimitive object
static void validate(const glTFMeshPrimitive* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFProperty*)val, err, errs);
    if (val->attributes.empty())
        errs.push_back({"missing requried value attributes", err.pathname()});
    validate_attr(val->attributes, "attributes", err, errs);
    validate_attr(val->indices, "indices", err, errs);
    validate_attr(val->material, "material", err, errs);
    validate_attr(val->mode, "mode", err, errs);
    validate_attr(val->targets, "targets", err, errs);
}

// Validates a glTFMesh object
static void validate(const glTFMesh* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFChildOfRootProperty*)val, err, errs);
    if (val->primitives.empty())
        errs.push_back({"missing requried value primitives", err.pathname()});
    validate_attr(val->primitives, "primitives", err, errs);
    validate_attr(val->weights, "weights", err, errs);
}

// Validates a glTFNode object
static void validate(const glTFNode* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFChildOfRootProperty*)val, err, errs);
    validate_attr(val->camera, "camera", err, errs);
    validate_attr(val->children, "children", err, errs);
    validate_attr(val->matrix, "matrix", err, errs);
    validate_attr(val->mesh, "mesh", err, errs);
    validate_attr(val->rotation, "rotation", err, errs);
    validate_attr(val->scale, "scale", err, errs);
    validate_attr(val->skin, "skin", err, errs);
    validate_attr(val->translation, "translation", err, errs);
    validate_attr(val->weights, "weights", err, errs);
}

// Validates a glTFSampler object
static void validate(const glTFSampler* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFChildOfRootProperty*)val, err, errs);
    validate_attr(val->magFilter, "magFilter", err, errs);
    validate_attr(val->minFilter, "minFilter", err, errs);
    validate_attr(val->wrapS, "wrapS", err, errs);
    validate_attr(val->wrapT, "wrapT", err, errs);
}

// Validates a glTFScene object
static void validate(const glTFScene* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFChildOfRootProperty*)val, err, errs);
    validate_attr(val->nodes, "nodes", err, errs);
}

// Validates a glTFSkin object
static void validate(const glTFSkin* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFChildOfRootProperty*)val, err, errs);
    validate_attr(val->inverseBindMatrices, "inverseBindMatrices", err, errs);
    if (val->joints.empty())
        errs.push_back({"missing requried value joints", err.pathname()});
    validate_attr(val->joints, "joints", err, errs);
    validate_attr(val->skeleton, "skeleton", err, errs);
}

// Validates a glTF object
static void validate(const glTF* val, parse_stack& err,
    std::vector<std::pair<std::string, std::string>>& errs) {
    if (!val) return;
    validate((const glTFProperty*)val, err, errs);
    validate_attr(val->accessors, "accessors", err, errs);
    validate_attr(val->animations, "animations", err, errs);
    if (!val->asset)
        errs.push_back({"missing requried value asset", err.pathname()});
    validate_attr(val->asset, "asset", err, errs);
    validate_attr(val->bufferViews, "bufferViews", err, errs);
    validate_attr(val->buffers, "buffers", err, errs);
    validate_attr(val->cameras, "cameras", err, errs);
    validate_attr(val->extensionsRequired, "extensionsRequired", err, errs);
    validate_attr(val->extensionsUsed, "extensionsUsed", err, errs);
    validate_attr(val->images, "images", err, errs);
    validate_attr(val->materials, "materials", err, errs);
    validate_attr(val->meshes, "meshes", err, errs);
    validate_attr(val->nodes, "nodes", err, errs);
    validate_attr(val->samplers, "samplers", err, errs);
    validate_attr(val->scene, "scene", err, errs);
    validate_attr(val->scenes, "scenes", err, errs);
    validate_attr(val->skins, "skins", err, errs);
    validate_attr(val->textures, "textures", err, errs);
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
// Get extension name
//
static inline std::string _get_extension(const std::string& filename) {
    auto pos = filename.rfind(".");
    if (pos == std::string::npos) return "";
    return filename.substr(pos);
}

//
// Get base name.
//
inline std::string _get_basename(const std::string& filename) {
    auto dirname = _get_dirname(filename);
    auto extension = _get_extension(filename);
    return filename.substr(
        dirname.size(), filename.size() - dirname.size() - extension.size());
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
// Loads a gltf.
//
glTF* load_gltf(const std::string& filename, bool load_bin, bool load_image,
    bool skip_missing) {
    // clear data
    auto gltf = std::unique_ptr<glTF>(new glTF());

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
    auto err = parse_stack();
    auto gltf_ = gltf.get();
    parse(gltf_, js, err);

    // load external resources
    auto dirname = _get_dirname(filename);
    if (load_bin) load_buffers(gltf.get(), dirname, skip_missing);
    if (load_image) load_images(gltf.get(), dirname, skip_missing);
    // done
    return gltf.release();
}

//
// Saves a gltf.
//
void save_gltf(const std::string& filename, const glTF* gltf, bool save_bin,
    bool save_image) {
    // dumps json
    auto js = json();
    auto err = parse_stack();
    dump(gltf, js, err);

    // save json
    _save_textfile(filename, js.dump(2));

    // save external resources
    auto dirname = _get_dirname(filename);
    if (save_bin) save_buffers(gltf, dirname);
    if (save_image) save_images(gltf, dirname);
}

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
glTF* load_binary_gltf(const std::string& filename, bool load_bin,
    bool load_image, bool skip_missing) {
    // clear data
    auto gltf = std::unique_ptr<glTF>(new glTF());

    // opens binary file
    auto f = std::fopen(filename.c_str(), "rb");
    if (!f) throw gltf_exception("could not load binary file");

    // read magic
    uint32_t magic;
    _fread(f, &magic, 1);
    if (magic != 0x46546C67) throw gltf_exception("corrupted glb format");

    // read version
    uint32_t version;
    _fread(f, &version, 1);
    if (version != 1 && version != 2)
        throw gltf_exception("unsupported glb version");

    // read length
    uint32_t length;
    _fread(f, &length, 1);

    // data
    auto json_bytes = std::vector<char>();
    auto buffer_bytes = std::vector<unsigned char>();
    uint32_t buffer_length = 0;

    if (version == 1) {
        // read content length and format
        uint32_t json_length, json_format;
        _fread(f, &json_length, 1);
        _fread(f, &json_format, 1);

        // read json bytes
        json_bytes.resize(json_length);
        _fread(f, json_bytes.data(), json_length);

        // read buffer bytes
        if (load_bin) {
            buffer_bytes.resize(length - json_length - 20);
            _fread(f, buffer_bytes.data(), (int)buffer_bytes.size());
            buffer_length = buffer_bytes.size();
        }
    }

    if (version == 2) {
        // read content length and format
        uint32_t json_length, json_format;
        _fread(f, &json_length, 1);
        _fread(f, &json_format, 1);
        if (json_format != 0x4E4F534A)
            throw gltf_exception("corrupt binary format");

        // read json bytes
        json_bytes.resize(json_length);
        _fread(f, json_bytes.data(), (int)json_bytes.size());

        // read content length and format
        uint32_t buffer_format;
        _fread(f, &buffer_length, 1);
        _fread(f, &buffer_format, 1);
        if (buffer_format != 0x004E4942)
            throw gltf_exception("corrupt binary format");

        // read buffer bytes
        if (load_bin) {
            buffer_bytes.resize(buffer_length);
            _fread(f, buffer_bytes.data(), (int)buffer_bytes.size());
        }
    }

    // load json
    auto js = json();
    try {
        json_bytes.push_back(0);
        js = json::parse(json_bytes.data());
    } catch (const std::exception&) {
        throw gltf_exception("could not load json");
    }

    // parse json
    auto err = parse_stack();
    auto gltf_ = gltf.get();
    parse(gltf_, js, err);

    // fix internal buffer
    auto buffer = gltf->buffers.at(0);
    buffer->byteLength = buffer_length;
    if (version == 2) buffer->uri = "";
    if (load_bin) { buffer->data = buffer_bytes; }

    // load external resources
    auto dirname = _get_dirname(filename);
    if (load_bin) load_buffers(gltf.get(), dirname, skip_missing);
    if (load_image) load_images(gltf.get(), dirname, skip_missing);

    // close
    fclose(f);

    // done
    return gltf.release();
}

//
// Saves a binary gltf.
//
void save_binary_gltf(const std::string& filename, const glTF* gltf,
    bool save_bin, bool save_image) {
    // opens binary file
    auto f = std::fopen(filename.c_str(), "wb");
    if (!f) throw gltf_exception("could not write binary file");

    // dumps json
    auto js = json();
    auto err = parse_stack();
    dump(gltf, js, err);

    // fix string
    auto js_str = js.dump(2);
    if (js_str.length() % 4) {
        auto count = js_str.length() % 4;
        for (auto c = 0; c < count; c++) js_str += " ";
    }
    uint32_t json_length = js_str.size();

    // internal buffer
    auto buffer = gltf->buffers.at(0);
    uint32_t buffer_length = buffer->byteLength;
    if (buffer_length % 4) buffer_length += 4 - buffer_length % 4;

    // write header
    uint32_t magic = 0x46546C67;
    _fwrite(f, &magic, 1);
    uint32_t version = 2;
    _fwrite(f, &version, 1);
    uint32_t length = 12 + 8 + json_length + 8 + buffer_length;
    _fread(f, &length, 1);

    // write json
    uint32_t json_type = 0x4E4F534A;
    _fwrite(f, &json_length, 1);
    _fwrite(f, &json_type, 1);
    _fwrite(f, js_str.data(), (int)json_length);

    if (save_bin) {
        uint32_t buffer_type = 0x004E4942;
        _fwrite(f, &buffer_length, 1);
        _fwrite(f, &buffer_type, 1);
        _fwrite(f, buffer->data.data(), (int)buffer->data.size());
        char pad = 0;
        for (auto i = 0; i < buffer_length - buffer->data.size(); i++) {
            _fwrite(f, &pad, 1);
        }
    }

    // close
    fclose(f);

    // save external resources
    auto dirname = _get_dirname(filename);
    if (save_bin) save_buffers(gltf, dirname);
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

#if 0
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
#endif

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
void load_buffers(glTF* gltf, const std::string& dirname, bool skip_missing) {
    for (auto buffer : gltf->buffers) {
        if (buffer->uri == "") continue;
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
void load_shaders(glTF* gltf, const std::string& dirname, bool skip_missing) {
#if 0
    for (auto kv : gltf->shaders) {
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
void load_images(glTF* gltf, const std::string& dirname, bool skip_missing) {
#ifndef YGL_NO_STBIMAGE

    for (auto image : gltf->images) {
        image->data = image_data();
        if (image->bufferView) {
            auto view = gltf->get(image->bufferView);
            auto buffer = gltf->get(view->buffer);
            if (!view || !buffer || view->byteStride) {
                if (skip_missing) continue;
                throw gltf_exception("invalid image buffer view");
            }
            auto ext = std::string();
            if (image->mimeType == glTFImageMimeType::ImagePng)
                ext = "png";
            else if (image->mimeType == glTFImageMimeType::ImageJpeg)
                ext = "jpg";
            else {
                if (skip_missing) continue;
                throw gltf_exception("unsupported image format");
            }
            yimg::load_image_from_memory(ext,
                buffer->data.data() + view->byteOffset, view->byteLength,
                image->data.width, image->data.height, image->data.ncomp,
                image->data.dataf, image->data.datab);
        } else if (_startsiwith(image->uri, "data:")) {
            // assume it is base64 and find ','
            auto pos = image->uri.find(',');
            auto ext = std::string();
            if (pos == image->uri.npos)
                throw gltf_exception("could not decode base64 data");
            auto header = image->uri.substr(0, pos);
            for (auto format : {"png", "jpg", "jpeg", "tga", "ppm", "hdr"})
                if (header.find(format) != header.npos) ext = format;
            if (ext.empty()) {
                if (skip_missing) continue;
                throw gltf_exception("unsupported embedded image format " +
                                     header.substr(0, pos));
            }
            // decode
            auto data = _base64::base64_decode(image->uri.substr(pos + 1));
            yimg::load_image_from_memory(ext, (unsigned char*)data.c_str(),
                (int)data.length(), image->data.width, image->data.height,
                image->data.ncomp, image->data.dataf, image->data.datab);
        } else {
            try {
                yimg::load_image(_fix_path(dirname + image->uri),
                    image->data.width, image->data.height, image->data.ncomp,
                    image->data.dataf, image->data.datab);
            } catch (...) {
                if (!skip_missing) throw;
            }
        }
    }

#endif
}

//
// Save buffer data.
//
void save_buffers(const glTF* gltf, const std::string& dirname) {
    for (auto buffer : gltf->buffers) {
        if (_startsiwith(buffer->uri, "data:"))
            throw gltf_exception("saving of embedded data not supported");
        _save_binfile(dirname + buffer->uri, buffer->data);
    }
}

//
// Save images.
//
void save_images(const glTF* gltf, const std::string& dirname) {
#ifndef YGL_NO_STBIMAGE

    for (auto image : gltf->images) {
        if (_startsiwith(image->uri, "data:"))
            throw gltf_exception("saving of embedded data not supported");
        if (!image->data.dataf.empty()) {
            yimg::save_image(dirname + image->uri, image->data.width,
                image->data.height, image->data.ncomp,
                image->data.dataf.data());
        }
        if (!image->data.datab.empty()) {
            yimg::save_image(dirname + image->uri, image->data.width,
                image->data.height, image->data.ncomp,
                image->data.datab.data());
        }
    }

#endif
}

inline vec_array_view::vec_array_view(
    const glTF* gltf, const glTFAccessor* accessor) {
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
}

inline float vec_array_view::get(int idx, int c) const {
    auto i = std::min(std::max(c, 0), ncomp() - 1);
    auto valb = _data + _stride * idx + i * _ctype_size(_ctype);
    // use double for integer conversion to attempt to maintain precision
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
                return (float)std::max(c / 127.0, -1.0);
            case glTFAccessorComponentType::UnsignedByte:
                return (float)(c / 255.0);
            case glTFAccessorComponentType::Short:
                return (float)(std::max(c / 32767.0, -1.0));
            case glTFAccessorComponentType::UnsignedShort:
                return (float)(c / 65535.0);
            case glTFAccessorComponentType::UnsignedInt:
                return (float)(std::max(c / 2147483647.0, -1.0));
            case glTFAccessorComponentType::NotSet:
                throw std::runtime_error("bad enum value");
                break;
        }
    }
    return 0;
}

template <int N>
inline ym::vec<float, N> vec_array_view::get(int idx) const {
    auto def = ym::vec<float, N>();
    for (auto i = 0; i < N; i++) def[i] = 0;
    return get<N>(idx, def);
}

template <int N>
inline ym::vec<float, N> vec_array_view::get(
    int idx, const ym::vec<float, N>& def) const {
    auto v = def;
    for (auto i = 0; i < std::min(_ncomp, N); i++) v[i] = get(idx, i);
    return v;
}

template <int N, int M>
inline ym::mat<float, N, M> vec_array_view::get(int idx) const {
    auto v = ym::mat<float, N, M>();
    if (_ncomp != N * M) throw gltf_exception("bad array view access");
    for (auto j = 0; j < M; j++)
        for (auto i = 0; i < N; i++) v[j][i] = get(idx, j * N + i);
    return v;
}

inline int vec_array_view::geti(int idx, int c) const {
    auto i = std::min(std::max(c, 0), ncomp() - 1);
    auto valb = _data + _stride * idx + i * _ctype_size(_ctype);
    // use double for integer conversion to attempt to maintain precision
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

template <int N>
inline ym::vec<int, N> vec_array_view::geti(int idx) const {
    auto v = ym::vec<int, N>();
    for (auto i = 0; i < std::min(_ncomp, N); i++) { v[i] = geti(idx, i); }
    for (auto i = std::min(_ncomp, N); i < N; i++) v[i] = 0;
    return v;
}

inline int vec_array_view::_num_components(glTFAccessorType type) {
    switch (type) {
        case glTFAccessorType::Scalar: return 1;
        case glTFAccessorType::Vec2: return 2;
        case glTFAccessorType::Vec3: return 3;
        case glTFAccessorType::Vec4: return 4;
        case glTFAccessorType::Mat2: return 4;
        case glTFAccessorType::Mat3: return 9;
        case glTFAccessorType::Mat4: return 16;
        default: assert(false); return 0;
    }
}

inline int vec_array_view::_ctype_size(
    glTFAccessorComponentType componentType) {
    switch (componentType) {
        case glTFAccessorComponentType::Byte: return 1;
        case glTFAccessorComponentType::UnsignedByte: return 1;
        case glTFAccessorComponentType::Short: return 2;
        case glTFAccessorComponentType::UnsignedShort: return 2;
        case glTFAccessorComponentType::UnsignedInt: return 4;
        case glTFAccessorComponentType::Float: return 4;
        default: assert(false); return 0;
    }
}

//
// element_attay_view implementation
//
inline element_array_view::element_array_view(
    const glTF* gltf, const glTFAccessor* accessor) {
    _size = accessor->count;
    _ctype = accessor->componentType;
    auto buffer_view = gltf->get(accessor->bufferView);
    _stride = (buffer_view->byteStride) ? buffer_view->byteStride :
                                          _ctype_size(_ctype);
    auto buffer = gltf->get(buffer_view->buffer);
    _data =
        buffer->data.data() + accessor->byteOffset + buffer_view->byteOffset;
    assert(accessor->type == glTFAccessorType::Scalar);
}

//
// element_attay_view implementation
//
inline int element_array_view::operator[](int idx) const {
    auto valb = _data + _stride * idx;
    switch (_ctype) {
        case glTFAccessorComponentType::Byte: return int(*(char*)valb);
        case glTFAccessorComponentType::UnsignedByte:
            return int(*(unsigned char*)valb);
        case glTFAccessorComponentType::Short: return int(*(short*)valb);
        case glTFAccessorComponentType::UnsignedShort:
            return int(*(unsigned short*)valb);
        case glTFAccessorComponentType::UnsignedInt:
            return int(*(unsigned int*)valb);
        default: assert(false); return 0;
    }
}

//
// element_attay_view implementation
//
inline int element_array_view::_ctype_size(
    glTFAccessorComponentType componentType) {
    switch (componentType) {
        case glTFAccessorComponentType::Byte: return 1;
        case glTFAccessorComponentType::UnsignedByte: return 1;
        case glTFAccessorComponentType::Short: return 2;
        case glTFAccessorComponentType::UnsignedShort: return 2;
        case glTFAccessorComponentType::UnsignedInt: return 4;
        case glTFAccessorComponentType::Float: assert(false); return 0;
        default: assert(false); return 0;
    }
}

//
// Math support
//
ym::mat4f node_transform(const glTFNode* node) {
    return ym::translation_mat4(node->translation) *
           ym::rotation_mat4(node->rotation) * ym::scaling_mat4(node->scale) *
           node->matrix;
}

//
// Math support
//
ym::mat4f node_transform(const node* node) {
    return ym::translation_mat4(node->translation) *
           ym::rotation_mat4(node->rotation) * ym::scaling_mat4(node->scale) *
           node->matrix;
}

//
// cleanup
//
mesh::~mesh() {
    for (auto e : shapes)
        if (e) delete e;
}

//
// cleanup
//
animation_group::~animation_group() {
    for (auto e : animations)
        if (e) delete e;
}

//
// cleanup
//
scene_group::~scene_group() {
    for (auto e : cameras)
        if (e) delete e;
    for (auto e : materials)
        if (e) delete e;
    for (auto e : meshes)
        if (e) delete e;
    for (auto e : textures)
        if (e) delete e;
    for (auto e : nodes)
        if (e) delete e;
    for (auto e : scenes)
        if (e) delete e;
    for (auto e : animations)
        if (e) delete e;
    for (auto e : skins)
        if (e) delete e;
}

//
// Flattens a gltf file into a flattened asset.
//
scene_group* gltf_to_scenes(const glTF* gltf, int scene_idx) {
    // clear asset
    auto scns = new ygltf::scene_group();

    // convert images
    for (auto gtxt : gltf->images) {
        auto txt = new texture();
        txt->name = gtxt->name;
        txt->path = (_startsiwith(gtxt->uri, "data:")) ?
                        std::string("inlines") :
                        gtxt->uri;
        txt->width = gtxt->data.width;
        txt->height = gtxt->data.height;
        txt->ncomp = gtxt->data.ncomp;
        txt->datab = gtxt->data.datab;
        txt->dataf = gtxt->data.dataf;
        scns->textures.push_back(txt);
    }

    // maps for translation
    static const auto filter_min_map =
        std::map<glTFSamplerMinFilter, texture_filter>{
            {glTFSamplerMinFilter::NotSet,
                texture_filter::linear_mipmap_linear},
            {glTFSamplerMinFilter::Linear, texture_filter::linear},
            {glTFSamplerMinFilter::Nearest, texture_filter::nearest},
            {glTFSamplerMinFilter::LinearMipmapLinear,
                texture_filter::linear_mipmap_linear},
            {glTFSamplerMinFilter::LinearMipmapNearest,
                texture_filter::linear_mipmap_nearest},
            {glTFSamplerMinFilter::NearestMipmapLinear,
                texture_filter::nearest_mipmap_linear},
            {glTFSamplerMinFilter::NearestMipmapNearest,
                texture_filter::nearest_mipmap_nearest},
        };
    static const auto filter_mag_map =
        std::map<glTFSamplerMagFilter, texture_filter>{
            {glTFSamplerMagFilter::NotSet, texture_filter::linear},
            {glTFSamplerMagFilter::Linear, texture_filter::linear},
            {glTFSamplerMagFilter::Nearest, texture_filter::nearest},
        };
    static const auto wrap_s_map = std::map<glTFSamplerWrapS, texture_wrap>{
        {glTFSamplerWrapS::NotSet, texture_wrap::repeat},
        {glTFSamplerWrapS::Repeat, texture_wrap::repeat},
        {glTFSamplerWrapS::ClampToEdge, texture_wrap::clamp},
        {glTFSamplerWrapS::MirroredRepeat, texture_wrap::mirror},
    };
    static const auto wrap_t_map = std::map<glTFSamplerWrapT, texture_wrap>{
        {glTFSamplerWrapT::NotSet, texture_wrap::repeat},
        {glTFSamplerWrapT::Repeat, texture_wrap::repeat},
        {glTFSamplerWrapT::ClampToEdge, texture_wrap::clamp},
        {glTFSamplerWrapT::MirroredRepeat, texture_wrap::mirror},
    };

    // add a texture
    auto add_texture = [gltf, scns](glTFTextureInfo* ginfo, texture*& txt,
                           texture_info*& info, bool normal = false,
                           bool occlusion = false) {
        auto gtxt = gltf->get(ginfo->index);
        if (!gtxt) return;
        txt = (!gtxt->source) ? nullptr : scns->textures[(int)gtxt->source];
        if (!txt) return;
        info = new texture_info();
        auto gsmp = gltf->get(gtxt->sampler);
        if (gsmp) {
            info->filter_mag = filter_mag_map.at(gsmp->magFilter);
            info->filter_min = filter_min_map.at(gsmp->minFilter);
            info->wrap_s = wrap_s_map.at(gsmp->wrapS);
            info->wrap_t = wrap_t_map.at(gsmp->wrapT);
        }
        if (normal) {
            auto ninfo = (glTFMaterialNormalTextureInfo*)ginfo;
            info->scale = ninfo->scale;
        }
        if (occlusion) {
            auto ninfo = (glTFMaterialOcclusionTextureInfo*)ginfo;
            info->scale = ninfo->strength;
        }
    };

    // convert materials
    for (auto gmat : gltf->materials) {
        auto mat = new material();
        mat->name = gmat->name;
        mat->emission = gmat->emissiveFactor;
        if (gmat->emissiveTexture) {
            add_texture(gmat->emissiveTexture, mat->emission_txt,
                mat->emission_txt_info);
        }
        if (gmat->pbrMetallicRoughness) {
            auto gmr = gmat->pbrMetallicRoughness;
            mat->metallic_roughness = new material_metallic_rooughness();
            auto mr = mat->metallic_roughness;
            mr->base = {gmr->baseColorFactor[0], gmr->baseColorFactor[1],
                gmr->baseColorFactor[2]};
            mr->opacity = gmr->baseColorFactor[3];
            mr->metallic = gmr->metallicFactor;
            mr->roughness = gmr->roughnessFactor;
            if (gmr->baseColorTexture) {
                add_texture(
                    gmr->baseColorTexture, mr->base_txt, mr->base_txt_info);
            }
            if (gmr->metallicRoughnessTexture) {
                add_texture(gmr->metallicRoughnessTexture, mr->metallic_txt,
                    mr->metallic_txt_info);
            }
        }
        if (gmat->pbrSpecularGlossiness) {
            auto gsg = gmat->pbrSpecularGlossiness;
            mat->specular_glossiness = new material_specular_glossiness();
            auto sg = mat->specular_glossiness;
            sg->diffuse = {gsg->diffuseFactor[0], gsg->diffuseFactor[1],
                gsg->diffuseFactor[2]};
            sg->opacity = gsg->diffuseFactor[3];
            sg->specular = gsg->specularFactor;
            sg->glossiness = gsg->glossinessFactor;
            if (gsg->diffuseTexture) {
                add_texture(
                    gsg->diffuseTexture, sg->diffuse_txt, sg->diffuse_txt_info);
            }
            if (gsg->specularGlossinessTexture) {
                add_texture(gsg->specularGlossinessTexture, sg->specular_txt,
                    sg->specular_txt_info);
            }
        }
        if (gmat->normalTexture) {
            add_texture(gmat->normalTexture, mat->normal_txt,
                mat->normal_txt_info, true, false);
        }
        if (gmat->occlusionTexture) {
            add_texture(gmat->occlusionTexture, mat->occlusion_txt,
                mat->occlusion_txt_info, false, true);
        }
        mat->double_sided = gmat->doubleSided;
        scns->materials.push_back(mat);
    }

    // convert meshes
    auto meshes = std::vector<std::vector<shape*>>();
    for (auto gmesh : gltf->meshes) {
        auto msh = new mesh();
        // primitives
        for (auto gprim : gmesh->primitives) {
            auto prim = new shape();
            if (gprim->material) {
                prim->mat = scns->materials[(int)gprim->material];
            }
            // vertex data
            for (auto gattr : gprim->attributes) {
                auto semantic = gattr.first;
                auto vals = vec_array_view(gltf, gltf->get(gattr.second));
                if (semantic == "POSITION") {
                    prim->pos.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        prim->pos.push_back(vals.get<3>(i));
                } else if (semantic == "NORMAL") {
                    prim->norm.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        prim->norm.push_back(vals.get<3>(i));
                } else if (semantic == "TEXCOORD" || semantic == "TEXCOORD_0") {
                    prim->texcoord.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        prim->texcoord.push_back(vals.get<2>(i));
                } else if (semantic == "TEXCOORD_1") {
                    prim->texcoord1.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        prim->texcoord1.push_back(vals.get<2>(i));
                } else if (semantic == "COLOR" || semantic == "COLOR_0") {
                    prim->color.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        prim->color.push_back(vals.get<4>(i, {0, 0, 0, 1}));
                } else if (semantic == "TANGENT") {
                    prim->tangsp.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        prim->tangsp.push_back(vals.get<4>(i));
                } else if (semantic == "WEIGHTS_0") {
                    prim->skin_weights.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        prim->skin_weights.push_back(vals.get<4>(i));
                } else if (semantic == "JOINTS_0") {
                    prim->skin_joints.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        prim->skin_joints.push_back(vals.geti<4>(i));
                } else {
                    // ignore
                }
            }
            // indices
            if (!gprim->indices) {
                switch (gprim->mode) {
                    case glTFMeshPrimitiveMode::Triangles: {
                        prim->triangles.reserve(prim->pos.size() / 3);
                        for (auto i = 0; i < prim->pos.size() / 3; i++) {
                            prim->triangles.push_back(
                                {i * 3 + 0, i * 3 + 1, i * 3 + 2});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::TriangleFan: {
                        prim->triangles.reserve(prim->pos.size() - 2);
                        for (auto i = 2; i < prim->pos.size(); i++) {
                            prim->triangles.push_back({0, i - 1, i});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::TriangleStrip: {
                        prim->triangles.reserve(prim->pos.size() - 2);
                        for (auto i = 2; i < prim->pos.size(); i++) {
                            prim->triangles.push_back({i - 2, i - 1, i});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::Lines: {
                        prim->lines.reserve(prim->pos.size() / 2);
                        for (auto i = 0; i < prim->pos.size() / 2; i++) {
                            prim->lines.push_back({i * 2 + 0, i * 2 + 1});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::LineLoop: {
                        prim->lines.reserve(prim->pos.size());
                        for (auto i = 1; i < prim->pos.size(); i++) {
                            prim->lines.push_back({i - 1, i});
                        }
                        prim->lines.back() = {(int)prim->pos.size() - 1, 0};
                    } break;
                    case glTFMeshPrimitiveMode::LineStrip: {
                        prim->lines.reserve(prim->pos.size() - 1);
                        for (auto i = 1; i < prim->pos.size(); i++) {
                            prim->lines.push_back({i - 1, i});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::NotSet:
                    case glTFMeshPrimitiveMode::Points: {
                        prim->points.reserve(prim->pos.size());
                        for (auto i = 0; i < prim->pos.size(); i++) {
                            prim->points.push_back(i);
                        }
                    } break;
                }
            } else {
                auto indices =
                    element_array_view(gltf, gltf->get(gprim->indices));
                switch (gprim->mode) {
                    case glTFMeshPrimitiveMode::Triangles: {
                        prim->triangles.reserve(indices.size());
                        for (auto i = 0; i < indices.size() / 3; i++) {
                            prim->triangles.push_back({indices[i * 3 + 0],
                                indices[i * 3 + 1], indices[i * 3 + 2]});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::TriangleFan: {
                        prim->triangles.reserve(indices.size() - 2);
                        for (auto i = 2; i < indices.size(); i++) {
                            prim->triangles.push_back(
                                {indices[0], indices[i - 1], indices[i]});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::TriangleStrip: {
                        prim->triangles.reserve(indices.size() - 2);
                        for (auto i = 2; i < indices.size(); i++) {
                            prim->triangles.push_back(
                                {indices[i - 2], indices[i - 1], indices[i]});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::Lines: {
                        prim->lines.reserve(indices.size() / 2);
                        for (auto i = 0; i < indices.size() / 2; i++) {
                            prim->lines.push_back(
                                {indices[i * 2 + 0], indices[i * 2 + 1]});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::LineLoop: {
                        prim->lines.reserve(indices.size());
                        for (auto i = 1; i < indices.size(); i++) {
                            prim->lines.push_back({indices[i - 1], indices[i]});
                        }
                        prim->lines.back() = {
                            indices[indices.size() - 1], indices[0]};
                    } break;
                    case glTFMeshPrimitiveMode::LineStrip: {
                        prim->lines.reserve(indices.size() - 1);
                        for (auto i = 1; i < indices.size(); i++) {
                            prim->lines.push_back({indices[i - 1], indices[i]});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::NotSet:
                    case glTFMeshPrimitiveMode::Points: {
                        prim->points.reserve(indices.size());
                        for (auto i = 0; i < indices.size(); i++) {
                            prim->points.push_back(indices[i]);
                        }
                    } break;
                }
            }

            // morph targets
            int target_index = 0;
            for (auto& gtarget : gprim->targets) {
                auto target = new shape_morph();
                for (auto gattr : gtarget) {
                    auto semantic = gattr.first;
                    auto vals = vec_array_view(gltf, gltf->get(gattr.second));
                    if (semantic == "POSITION") {
                        target->pos.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            target->pos.push_back(vals.get<3>(i));
                    } else if (semantic == "NORMAL") {
                        target->norm.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            target->norm.push_back(vals.get<3>(i));
                    } else if (semantic == "TANGENT") {
                        target->tangsp.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            target->tangsp.push_back(vals.get<3>(i));
                    } else {
                        // ignore
                    }
                }
                if (target_index < (int)gmesh->weights.size() - 1)
                    target->weight = gmesh->weights[target_index];
                target_index++;
                prim->morph_targets.push_back(target);
            }
            msh->shapes.push_back(prim);
        }
        scns->meshes.push_back(msh);
    }

    // convert cameras
    for (auto gcam : gltf->cameras) {
        auto cam = new camera();
        cam->name = gcam->name;
        cam->ortho = gcam->type == glTFCameraType::Orthographic;
        if (cam->ortho) {
            auto ortho = gcam->orthographic;
            cam->yfov = ortho->ymag;
            cam->aspect = ortho->xmag / ortho->ymag;
            cam->near = ortho->znear;
            cam->far = ortho->zfar;
        } else {
            auto persp = gcam->perspective;
            cam->yfov = persp->yfov;
            cam->aspect = persp->aspectRatio;
            if (!cam->aspect) cam->aspect = 16.0f / 9.0f;
            cam->near = persp->znear;
            cam->far = persp->zfar;
        }
        scns->cameras.push_back(cam);
    }

    // convert nodes
    for (auto gnode : gltf->nodes) {
        auto node = new ygltf::node();
        node->name = gnode->name;
        node->local_xform = node_transform(gnode);
        node->cam =
            (!gnode->camera) ? nullptr : scns->cameras[(int)gnode->camera];
        node->msh = (!gnode->mesh) ? nullptr : scns->meshes[(int)gnode->mesh];
        node->translation = gnode->translation;
        node->rotation = gnode->rotation;
        node->scale = gnode->scale;
        node->matrix = gnode->matrix;
        node->morph_weights = gnode->weights;
        scns->nodes.push_back(node);
    }

    // set up children pointers
    auto is_root = std::vector<bool>(gltf->nodes.size(), true);
    for (auto nid = 0; nid < gltf->nodes.size(); nid++) {
        auto gnode = gltf->nodes[nid];
        auto node = scns->nodes[nid];
        for (auto n : gnode->children)
            node->children.push_back(scns->nodes[(int)n]);
        for (auto n : gnode->children) is_root[(int)n] = false;
    }

    // fix node morph weights
    for (auto node : scns->nodes) {
        if (!node->msh) continue;
        for (auto shp : node->msh->shapes) {
            if (node->morph_weights.size() < shp->morph_targets.size()) {
                node->morph_weights.resize(shp->morph_targets.size());
            }
        }
    }

    // convert animations
    for (auto ganim : gltf->animations) {
        auto anim_group = new animation_group();
        anim_group->name = ganim->name;
        std::map<ym::vec2i, animation*> sampler_map;
        for (auto gchannel : ganim->channels) {
            if (sampler_map.find({(int)gchannel->sampler,
                    (int)gchannel->target->path}) == sampler_map.end()) {
                auto gsampler = ganim->get(gchannel->sampler);
                auto keyframes = new animation();
                auto input_view =
                    vec_array_view(gltf, gltf->get(gsampler->input));
                keyframes->time.resize(input_view.size());
                for (auto i = 0; i < input_view.size(); i++)
                    keyframes->time[i] = input_view.get<1>(i)[0];
                keyframes->interp =
                    (animation_interpolation)gsampler->interpolation;
                auto output_view =
                    vec_array_view(gltf, gltf->get(gsampler->output));
                switch (gchannel->target->path) {
                    case glTFAnimationChannelTargetPath::Translation: {
                        keyframes->translation.reserve(output_view.size());
                        for (auto i = 0; i < output_view.size(); i++)
                            keyframes->translation.push_back(
                                output_view.get<3>(i));
                    } break;
                    case glTFAnimationChannelTargetPath::Rotation: {
                        keyframes->rotation.reserve(output_view.size());
                        for (auto i = 0; i < output_view.size(); i++)
                            keyframes->rotation.push_back(
                                (ym::quat4f)output_view.get<4>(i));
                    } break;
                    case glTFAnimationChannelTargetPath::Scale: {
                        keyframes->scale.reserve(output_view.size());
                        for (auto i = 0; i < output_view.size(); i++)
                            keyframes->scale.push_back(output_view.get<3>(i));
                    } break;
                    case glTFAnimationChannelTargetPath::Weights: {
                        // get a node that it refers to
                        auto ncomp = 0;
                        auto gnode = gltf->get(gchannel->target->node);
                        auto gmesh = gltf->get(gnode->mesh);
                        if (gmesh) {
                            for (auto gshp : gmesh->primitives) {
                                ncomp =
                                    ym::max((int)gshp->targets.size(), ncomp);
                            }
                        }
                        if (ncomp) {
                            auto values = std::vector<float>();
                            values.reserve(output_view.size());
                            for (auto i = 0; i < output_view.size(); i++)
                                values.push_back(output_view.get(i, 0));
                            keyframes->morph_weights.resize(
                                values.size() / ncomp);
                            for (auto i = 0;
                                 i < keyframes->morph_weights.size(); i++) {
                                keyframes->morph_weights[i].resize(ncomp);
                                for (auto j = 0; j < ncomp; j++)
                                    keyframes->morph_weights[i][j] =
                                        values[i * ncomp + j];
                            }
                        }
                    } break;
                    default: {
                        // skip
                    }
                }
                sampler_map[{(int)gchannel->sampler,
                    (int)gchannel->target->path}] = keyframes;
                anim_group->animations.push_back(keyframes);
            }
            sampler_map
                .at({(int)gchannel->sampler, (int)gchannel->target->path})
                ->nodes.push_back(scns->nodes[(int)gchannel->target->node]);
        }
        scns->animations.push_back(anim_group);
    }

    // convert skins
    for (auto gskin : gltf->skins) {
        auto skin = new ygltf::skin();
        skin->name = gskin->name;
        for (auto gnode : gskin->joints)
            skin->joints.push_back(scns->nodes[(int)gnode]);
        skin->root = scns->nodes[(int)gskin->skeleton];
        if (!gskin->inverseBindMatrices) {
            skin->pose_matrices.assign(skin->joints.size(), ym::identity_mat4f);
        } else {
            auto pose_matrix_view =
                vec_array_view(gltf, gltf->get(gskin->inverseBindMatrices));
            skin->pose_matrices.resize(skin->joints.size());
            assert(pose_matrix_view.size() == skin->joints.size());
            assert(pose_matrix_view.ncomp() == 16);
            for (auto i = 0; i < pose_matrix_view.size(); i++) {
                skin->pose_matrices[i] = pose_matrix_view.get<4, 4>(i);
            }
        }
        scns->skins.push_back(skin);
    }

    // set skin pointers
    for (auto nid = 0; nid < gltf->nodes.size(); nid++) {
        if (!gltf->nodes[nid]->skin) continue;
        scns->nodes[nid]->skn = scns->skins[(int)gltf->nodes[nid]->skin];
    }

    // convert scenes
    for (auto gscn : gltf->scenes) {
        auto scn = new scene();
        scn->name = gscn->name;
        for (auto n : gscn->nodes) scn->nodes.push_back(scns->nodes[(int)n]);
        scns->scenes.push_back(scn);
    }
    if (gltf->scene) { scns->default_scene = scns->scenes[(int)gltf->scene]; }

    // update transforms
    update_transforms(scns);

    return scns;
}

//
// helper
//
template <typename T>
static inline int index(const std::vector<T*>& vec, T* val) {
    auto pos = std::find(vec.begin(), vec.end(), val);
    if (pos == vec.end()) return -1;
    return (int)(pos - vec.begin());
}

//
// Unflattnes gltf
//
glTF* scenes_to_gltf(const scene_group* scns, const std::string& buffer_uri) {
    auto gltf = std::unique_ptr<glTF>(new glTF());

    // add asset info
    gltf->asset = new glTFAsset();
    gltf->asset->generator = "Yocto/gltf";
    gltf->asset->version = "2.0";

    // convert cameras
    for (auto cam : scns->cameras) {
        auto gcam = new glTFCamera();
        gcam->name = cam->name;
        gcam->type = (cam->ortho) ? glTFCameraType::Orthographic :
                                    glTFCameraType::Perspective;
        if (cam->ortho) {
            auto ortho = new glTFCameraOrthographic();
            ortho->ymag = cam->yfov;
            ortho->xmag = cam->aspect * cam->yfov;
            ortho->znear = cam->near;
            ortho->znear = cam->far;
            gcam->orthographic = ortho;
        } else {
            auto persp = new glTFCameraPerspective();
            persp->yfov = cam->yfov;
            persp->aspectRatio = cam->aspect;
            persp->znear = cam->near;
            persp->zfar = cam->far;
            gcam->perspective = persp;
        }
        gltf->cameras.push_back(gcam);
    }

    // convert images
    for (auto txt : scns->textures) {
        auto gimg = new glTFImage();
        gimg->uri = txt->path;
        gimg->data.width = txt->width;
        gimg->data.height = txt->height;
        gimg->data.ncomp = txt->ncomp;
        gimg->data.datab = txt->datab;
        gimg->data.dataf = txt->dataf;
        gltf->images.push_back(gimg);
    }

    // conversion maps
    static const auto wrap_s_map = std::map<texture_wrap, glTFSamplerWrapS>{
        {texture_wrap::repeat, glTFSamplerWrapS::Repeat},
        {texture_wrap::clamp, glTFSamplerWrapS::ClampToEdge},
        {texture_wrap::mirror, glTFSamplerWrapS::MirroredRepeat},
    };
    static const auto wrap_t_map = std::map<texture_wrap, glTFSamplerWrapT>{
        {texture_wrap::repeat, glTFSamplerWrapT::Repeat},
        {texture_wrap::clamp, glTFSamplerWrapT::ClampToEdge},
        {texture_wrap::mirror, glTFSamplerWrapT::MirroredRepeat},
    };
    static const auto texture_min_map =
        std::map<texture_filter, glTFSamplerMinFilter>{
            {texture_filter::linear, glTFSamplerMinFilter::Linear},
            {texture_filter::nearest, glTFSamplerMinFilter::Nearest},
            {texture_filter::linear_mipmap_linear,
                glTFSamplerMinFilter::LinearMipmapLinear},
            {texture_filter::linear_mipmap_nearest,
                glTFSamplerMinFilter::LinearMipmapNearest},
            {texture_filter::nearest_mipmap_linear,
                glTFSamplerMinFilter::NearestMipmapNearest},
            {texture_filter::nearest_mipmap_nearest,
                glTFSamplerMinFilter::NearestMipmapNearest},
        };
    static const auto texture_mag_map =
        std::map<texture_filter, glTFSamplerMagFilter>{
            {texture_filter::linear, glTFSamplerMagFilter::Linear},
            {texture_filter::nearest, glTFSamplerMagFilter::Nearest},
        };

    // add a texture and sampler
    auto add_texture = [&gltf, scns](texture* txt, texture_info* info) {
        if (!txt) return glTFid<glTFTexture>{-1};
        auto gtxt = new glTFTexture();
        gtxt->name = txt->name;
        gtxt->source = glTFid<glTFImage>(index(scns->textures, txt));
        if (info && !info->is_default()) {
            auto gsmp = new glTFSampler();
            gsmp->wrapS = wrap_s_map.at(info->wrap_s);
            gsmp->wrapT = wrap_t_map.at(info->wrap_t);
            gsmp->minFilter = texture_min_map.at(info->filter_min);
            gsmp->magFilter = texture_mag_map.at(info->filter_mag);
            gtxt->sampler = glTFid<glTFSampler>(gltf->samplers.size());
            gltf->samplers.push_back(gsmp);
        }
        gltf->textures.push_back(gtxt);
        return glTFid<glTFTexture>{(int)gltf->textures.size() - 1};
    };

    // convert materials
    for (auto mat : scns->materials) {
        auto gmat = new glTFMaterial();
        gmat->name = mat->name;
        gmat->emissiveFactor = mat->emission;
        if (mat->emission_txt) {
            gmat->emissiveTexture = new glTFTextureInfo();
            gmat->emissiveTexture->index =
                add_texture(mat->emission_txt, mat->emission_txt_info);
        }
        if (mat->metallic_roughness) {
            gmat->pbrMetallicRoughness = new glTFMaterialPbrMetallicRoughness();
            auto gmr = gmat->pbrMetallicRoughness;
            auto mr = mat->metallic_roughness;
            gmr->baseColorFactor = {
                mr->base[0], mr->base[1], mr->base[2], mr->opacity};
            gmr->metallicFactor = mr->metallic;
            gmr->roughnessFactor = mr->roughness;
            if (mr->base_txt) {
                gmr->baseColorTexture = new glTFTextureInfo();
                gmr->baseColorTexture->index =
                    add_texture(mr->base_txt, mr->base_txt_info);
            }
            if (mr->metallic_txt) {
                gmr->metallicRoughnessTexture = new glTFTextureInfo();
                gmr->metallicRoughnessTexture->index =
                    add_texture(mr->metallic_txt, mr->metallic_txt_info);
            }
        }
        if (mat->specular_glossiness) {
            gmat->pbrSpecularGlossiness =
                new glTFMaterialPbrSpecularGlossiness();
            auto gsg = gmat->pbrSpecularGlossiness;
            auto sg = mat->specular_glossiness;
            gsg->diffuseFactor = {
                sg->diffuse[0], sg->diffuse[1], sg->diffuse[2], sg->opacity};
            gsg->specularFactor = sg->specular;
            gsg->glossinessFactor = sg->glossiness;
            if (sg->diffuse_txt) {
                gsg->diffuseTexture = new glTFTextureInfo();
                gsg->diffuseTexture->index =
                    add_texture(sg->diffuse_txt, sg->diffuse_txt_info);
            }
            if (sg->specular_txt) {
                gsg->specularGlossinessTexture = new glTFTextureInfo();
                gsg->specularGlossinessTexture->index =
                    add_texture(sg->specular_txt, sg->specular_txt_info);
            }
        }
        if (mat->normal_txt) {
            gmat->normalTexture = new glTFMaterialNormalTextureInfo();
            gmat->normalTexture->index =
                add_texture(mat->normal_txt, mat->normal_txt_info);
            if (gmat->normalTexture && mat->normal_txt_info)
                gmat->normalTexture->scale = mat->normal_txt_info->scale;
        }
        if (mat->occlusion_txt) {
            gmat->occlusionTexture = new glTFMaterialOcclusionTextureInfo();
            gmat->occlusionTexture->index =
                add_texture(mat->occlusion_txt, mat->occlusion_txt_info);
            if (gmat->occlusionTexture && mat->occlusion_txt_info)
                gmat->occlusionTexture->strength =
                    mat->occlusion_txt_info->scale;
        }
        gmat->doubleSided = mat->double_sided;
        gltf->materials.push_back(gmat);
    }

    // init buffers
    auto gbuffer = new glTFBuffer();
    gltf->buffers.push_back(gbuffer);
    gbuffer->uri = buffer_uri;

    // attribute handling
    auto add_accessor = [&gltf, gbuffer](const std::string& name,
                            glTFAccessorType type,
                            glTFAccessorComponentType ctype, int count,
                            int csize, const void* data, bool save_min_max) {
        gltf->bufferViews.push_back(new glTFBufferView());
        auto bufferView = gltf->bufferViews.back();
        bufferView->buffer = glTFid<glTFBuffer>(0);
        bufferView->byteOffset = (int)gbuffer->data.size();
        bufferView->byteStride = 0;
        bufferView->byteLength = count * csize;
        gbuffer->data.resize(gbuffer->data.size() + bufferView->byteLength);
        gbuffer->byteLength += bufferView->byteLength;
        auto ptr = gbuffer->data.data() + gbuffer->data.size() -
                   bufferView->byteLength;
        bufferView->target = glTFBufferViewTarget::ArrayBuffer;
        memcpy(ptr, data, bufferView->byteLength);
        gltf->accessors.push_back(new glTFAccessor());
        auto accessor = gltf->accessors.back();
        accessor->bufferView =
            glTFid<glTFBufferView>((int)gltf->bufferViews.size() - 1);
        accessor->byteOffset = 0;
        accessor->componentType = ctype;
        accessor->count = count;
        accessor->type = type;
        if (save_min_max && count &&
            ctype == glTFAccessorComponentType::Float) {
            switch (type) {
                case glTFAccessorType::Scalar: {
                    auto bbox = ym::make_bbox(count, (ym::vec1f*)data);
                    accessor->min = {bbox.min.x};
                    accessor->max = {bbox.max.x};
                } break;
                case glTFAccessorType::Vec2: {
                    auto bbox = ym::make_bbox(count, (ym::vec2f*)data);
                    accessor->min = {bbox.min.x, bbox.min.y};
                    accessor->max = {bbox.max.x, bbox.max.y};
                } break;
                case glTFAccessorType::Vec3: {
                    auto bbox = ym::make_bbox(count, (ym::vec3f*)data);
                    accessor->min = {bbox.min.x, bbox.min.y, bbox.min.z};
                    accessor->max = {bbox.max.x, bbox.max.y, bbox.max.z};
                } break;
                case glTFAccessorType::Vec4: {
                    auto bbox = ym::make_bbox(count, (ym::vec4f*)data);
                    accessor->min = {
                        bbox.min.x, bbox.min.y, bbox.min.z, bbox.min.w};
                    accessor->max = {
                        bbox.max.x, bbox.max.y, bbox.max.z, bbox.max.w};
                } break;
                default: break;
            }
        }
        return glTFid<glTFAccessor>((int)gltf->accessors.size() - 1);
    };

    // convert meshes
    for (auto i = 0; i < scns->meshes.size(); i++) {
        auto msh = scns->meshes[i];
        auto gid = "mesh" + std::to_string(i);
        auto gmesh = new glTFMesh();
        gmesh->name = msh->name;
        for (auto j = 0; j < msh->shapes.size(); j++) {
            auto gprim = msh->shapes[j];
            auto pid = std::to_string(j);
            auto prim = new glTFMeshPrimitive();
            prim->material =
                glTFid<glTFMaterial>(index(scns->materials, gprim->mat));
            if (!gprim->pos.empty())
                prim->attributes["POSITION"] = add_accessor(gid + pid + "_pos",
                    glTFAccessorType::Vec3, glTFAccessorComponentType::Float,
                    (int)gprim->pos.size(), sizeof(ym::vec3f),
                    gprim->pos.data(), true);
            if (!gprim->norm.empty())
                prim->attributes["NORMAL"] = add_accessor(gid + pid + "_norm",
                    glTFAccessorType::Vec3, glTFAccessorComponentType::Float,
                    (int)gprim->norm.size(), sizeof(ym::vec3f),
                    gprim->norm.data(), false);
            if (!gprim->texcoord.empty())
                prim->attributes["TEXCOORD_0"] = add_accessor(
                    gid + pid + "_texcoord", glTFAccessorType::Vec2,
                    glTFAccessorComponentType::Float,
                    (int)gprim->texcoord.size(), sizeof(ym::vec2f),
                    gprim->texcoord.data(), false);
            if (!gprim->texcoord1.empty())
                prim->attributes["TEXCOORD_1"] = add_accessor(
                    gid + pid + "_texcoord1", glTFAccessorType::Vec2,
                    glTFAccessorComponentType::Float,
                    (int)gprim->texcoord1.size(), sizeof(ym::vec2f),
                    gprim->texcoord1.data(), false);
            if (!gprim->color.empty())
                prim->attributes["COLOR_0"] = add_accessor(gid + pid + "_color",
                    glTFAccessorType::Vec4, glTFAccessorComponentType::Float,
                    (int)gprim->color.size(), sizeof(ym::vec4f),
                    gprim->color.data(), false);
            if (!gprim->skin_weights.empty())
                prim->attributes["WEIGHTS_0"] = add_accessor(
                    gid + pid + "_skin_weights", glTFAccessorType::Vec4,
                    glTFAccessorComponentType::Float,
                    (int)gprim->skin_weights.size(), sizeof(ym::vec4f),
                    gprim->skin_weights.data(), false);
            if (!gprim->skin_joints.empty()) {
                using ushort = unsigned short;
                auto joints_short = std::vector<ym::vec<ushort, 4>>();
                joints_short.reserve(gprim->skin_joints.size());
                for (auto&& j : gprim->skin_joints)
                    joints_short.push_back(
                        {(ushort)j.x, (ushort)j.y, (ushort)j.z, (ushort)j.w});
                prim->attributes["JOINTS_0"] = add_accessor(
                    gid + pid + "_skin_joints", glTFAccessorType::Vec4,
                    glTFAccessorComponentType::UnsignedShort,
                    (int)joints_short.size(), sizeof(ushort) * 4,
                    joints_short.data(), false);
            }
            // auto elem_as_uint = gprim->pos.size() >
            // std::numeric_limits<unsigned short>::max();
            if (!gprim->points.empty()) {
                prim->indices = add_accessor(gid + pid + "_points",
                    glTFAccessorType::Scalar,
                    glTFAccessorComponentType::UnsignedInt,
                    (int)gprim->points.size(), sizeof(int),
                    (int*)gprim->points.data(), false);
                prim->mode = glTFMeshPrimitiveMode::Points;
            } else if (!gprim->lines.empty()) {
                prim->indices =
                    add_accessor(gid + pid + "_lines", glTFAccessorType::Scalar,
                        glTFAccessorComponentType::UnsignedInt,
                        (int)gprim->lines.size() * 2, sizeof(int),
                        (int*)gprim->lines.data(), false);
                prim->mode = glTFMeshPrimitiveMode::Lines;
            } else if (!gprim->triangles.empty()) {
                prim->indices = add_accessor(gid + pid + "_triangles",
                    glTFAccessorType::Scalar,
                    glTFAccessorComponentType::UnsignedInt,
                    (int)gprim->triangles.size() * 3, sizeof(int),
                    (int*)gprim->triangles.data(), false);
                prim->mode = glTFMeshPrimitiveMode::Triangles;
            } else {
                assert(false);
            }
            auto target_index = 0;
            for (auto target : gprim->morph_targets) {
                auto mid = std::to_string(target_index++);
                prim->targets.push_back({});
                if (!target->pos.empty()) {
                    prim->targets.back()["POSITION"] = add_accessor(
                        gid + pid + "_" + mid + "_pos", glTFAccessorType::Vec3,
                        glTFAccessorComponentType::Float,
                        (int)target->pos.size(), sizeof(ym::vec3f),
                        target->pos.data(), true);
                } else if (!target->norm.empty()) {
                    prim->targets.back()["NORMAL"] = add_accessor(
                        gid + pid + "_" + mid + "_norm", glTFAccessorType::Vec3,
                        glTFAccessorComponentType::Float,
                        (int)target->norm.size(), sizeof(ym::vec3f),
                        target->norm.data(), true);
                } else if (!target->tangsp.empty()) {
                    prim->targets.back()["TANGENT"] = add_accessor(
                        gid + pid + "_" + mid + "_tang", glTFAccessorType::Vec3,
                        glTFAccessorComponentType::Float,
                        (int)target->tangsp.size(), sizeof(ym::vec3f),
                        target->tangsp.data(), true);
                } else {
                    assert(false);
                }
            }
            gmesh->primitives.push_back(prim);
        }
        gltf->meshes.push_back(gmesh);
    }

    // nodes
    for (auto node : scns->nodes) {
        auto gnode = new glTFNode();
        gnode->name = node->name;
        gnode->camera = glTFid<glTFCamera>(index(scns->cameras, node->cam));
        gnode->mesh = glTFid<glTFMesh>(index(scns->meshes, node->msh));
        gnode->matrix = node->matrix;
        gnode->translation = node->translation;
        gnode->rotation = node->rotation;
        gnode->scale = node->scale;
        gnode->weights = node->morph_weights;
        gltf->nodes.push_back(gnode);
    }

    // fix children
    auto nid = 0;
    for (auto node : scns->nodes) {
        auto gnode = gltf->nodes[nid++];
        for (auto child : node->children)
            gnode->children.push_back(
                glTFid<glTFNode>(index(scns->nodes, child)));
    }

    // skins
    for (auto sk : scns->skins) {
        auto gsk = new glTFSkin();
        gsk->name = sk->name;
        gsk->skeleton = glTFid<glTFNode>(index(scns->nodes, sk->root));
        for (auto joint : sk->joints) {
            gsk->joints.push_back(glTFid<glTFNode>(index(scns->nodes, joint)));
        }
        if (!sk->pose_matrices.empty()) {
            gsk->inverseBindMatrices = add_accessor(sk->name + "_pose",
                glTFAccessorType::Mat4, glTFAccessorComponentType::Float,
                (int)sk->pose_matrices.size(), sizeof(ym::mat4f),
                sk->pose_matrices.data(), false);
        }
        gltf->skins.push_back(gsk);
    }

    // fix skin references
    nid = 0;
    for (auto node : scns->nodes) {
        auto gnode = gltf->nodes[nid++];
        if (node->skn) {
            gnode->skin = glTFid<glTFSkin>{index(scns->skins, node->skn)};
        }
    }

    // interpolation map
    static const auto interpolation_map =
        std::map<animation_interpolation, glTFAnimationSamplerInterpolation>{
            {animation_interpolation::step,
                glTFAnimationSamplerInterpolation::Step},
            {animation_interpolation::linear,
                glTFAnimationSamplerInterpolation::Linear},
            {animation_interpolation::cubic,
                glTFAnimationSamplerInterpolation::Cubicspline},
            {animation_interpolation::catmull_rom,
                glTFAnimationSamplerInterpolation::Catmullromspline},
        };

    // animation
    for (auto anims : scns->animations) {
        auto ganim = new glTFAnimation();
        ganim->name = anims->name;
        auto count = 0;
        for (auto anim : anims->animations) {
            auto aid = ganim->name + std::to_string(count++);
            auto gsmp = new glTFAnimationSampler();
            gsmp->input = add_accessor(aid + "_time", glTFAccessorType::Scalar,
                glTFAccessorComponentType::Float, (int)anim->time.size(),
                sizeof(float), anim->time.data(), false);
            auto path = glTFAnimationChannelTargetPath::NotSet;
            if (!anim->translation.empty()) {
                gsmp->output = add_accessor(aid + "_translation",
                    glTFAccessorType::Vec3, glTFAccessorComponentType::Float,
                    (int)anim->translation.size(), sizeof(ym::vec3f),
                    anim->translation.data(), false);
                path = glTFAnimationChannelTargetPath::Translation;
            } else if (!anim->rotation.empty()) {
                gsmp->output = add_accessor(aid + "_rotation",
                    glTFAccessorType::Vec4, glTFAccessorComponentType::Float,
                    (int)anim->rotation.size(), sizeof(ym::vec4f),
                    anim->rotation.data(), false);
                path = glTFAnimationChannelTargetPath::Rotation;
            } else if (!anim->scale.empty()) {
                gsmp->output = add_accessor(aid + "_scale",
                    glTFAccessorType::Vec3, glTFAccessorComponentType::Float,
                    (int)anim->scale.size(), sizeof(ym::vec3f),
                    anim->scale.data(), false);
                path = glTFAnimationChannelTargetPath::Scale;
            } else if (!anim->morph_weights.empty()) {
                auto values = std::vector<float>();
                values.reserve(
                    anim->morph_weights.size() * anim->morph_weights[0].size());
                for (auto i = 0; i < anim->morph_weights.size(); i++) {
                    values.insert(values.end(), anim->morph_weights[i].begin(),
                        anim->morph_weights[i].end());
                }
                gsmp->output = add_accessor(aid + "_weights",
                    glTFAccessorType::Scalar, glTFAccessorComponentType::Float,
                    (int)values.size(), sizeof(float), values.data(), false);
                path = glTFAnimationChannelTargetPath::Weights;
            } else {
            }
            gsmp->interpolation = interpolation_map.at(anim->interp);
            for (auto node : anim->nodes) {
                auto gchan = new glTFAnimationChannel();
                gchan->sampler =
                    glTFid<glTFAnimationSampler>{(int)ganim->samplers.size()};
                gchan->target = new glTFAnimationChannelTarget();
                gchan->target->node =
                    glTFid<glTFNode>{index(scns->nodes, node)};
                gchan->target->path = path;
                ganim->channels.push_back(gchan);
            }
            ganim->samplers.push_back(gsmp);
        }
        gltf->animations.push_back(ganim);
    }

    // scenes
    for (auto scn : scns->scenes) {
        auto gscn = new glTFScene();
        gscn->name = scn->name;
        for (auto child : scn->nodes)
            gscn->nodes.push_back(glTFid<glTFNode>(index(scns->nodes, child)));
        gltf->scenes.push_back(gscn);
    }
    gltf->scene = glTFid<glTFScene>(index(scns->scenes, scns->default_scene));

    // done
    return gltf.release();
}

//
// Validate a gltf. Missing many validation as of this version.
//
std::vector<std::pair<std::string, std::string>> validate_gltf(
    const glTF* gltf) {
    auto errs = std::vector<std::pair<std::string, std::string>>();
    parse_stack err;
    validate(gltf, err, errs);
    return errs;
}

//
// Load scene
//
scene_group* load_scenes(
    const std::string& filename, bool load_textures, bool skip_missing) {
    auto ext = _get_extension(filename);
    auto gltf = std::unique_ptr<glTF>();
    if (ext != ".glb") {
        gltf = std::unique_ptr<glTF>(
            load_gltf(filename, true, load_textures, skip_missing));
    } else {
        gltf = std::unique_ptr<glTF>(
            load_binary_gltf(filename, true, load_textures, skip_missing));
    }
    return gltf_to_scenes(gltf.get());
}

///
/// Save scene
///
void save_scenes(
    const std::string& filename, const scene_group* scn, bool save_textures) {
    auto buffer_uri = _get_basename(filename) + ".bin";
    auto gltf = std::unique_ptr<glTF>(scenes_to_gltf(scn, buffer_uri));
    save_gltf(filename, gltf.get(), true, save_textures);
}

#ifndef YOBJ_NO_SHAPE

//
// Computes a scene bounding box
//
ym::bbox3f compute_scene_bounds(const scene_group* scn) {
    auto bbox_meshes =
        std::map<mesh*, ym::bbox3f>{{nullptr, ym::invalid_bbox3f}};
    for (auto mesh : scn->meshes) {
        bbox_meshes[mesh] = ym::invalid_bbox3f;
        auto& bbox = bbox_meshes[mesh];
        for (auto shp : mesh->shapes)
            for (auto& p : shp->pos) bbox += p;
    }
    auto bbox = ym::invalid_bbox3f;
    if (!scn->nodes.empty()) {
        for (auto ist : scn->nodes) {
            if (ist->msh)
                bbox += ym::transform_bbox(
                    ym::mat4f(ist->xform), bbox_meshes.at(ist->msh));
        }
    } else {
        for (auto mesh : scn->meshes) bbox += bbox_meshes[mesh];
    }
    return bbox;
}

//
// Add missing data to the scene.
//
void add_normals(scene_group* scn) {
    for (auto msh : scn->meshes) {
        for (auto shp : msh->shapes) {
            if (!shp->norm.empty()) continue;
            shp->norm.resize(shp->pos.size());
            if (!shp->points.empty()) {
                shp->norm.assign(shp->pos.size(), {0, 0, 1});
            } else if (!shp->lines.empty()) {
                ym::compute_tangents(shp->lines, shp->pos, shp->norm);
            } else if (!shp->triangles.empty()) {
                ym::compute_normals(shp->triangles, shp->pos, shp->norm);
            }
        }
    }
}

//
// Add missing data to the scene.
//
void add_tangent_space(scene_group* scn) {
    for (auto msh : scn->meshes) {
        for (auto shp : msh->shapes) {
            if (shp->triangles.empty()) continue;
            if (!shp->tangsp.empty() || shp->texcoord.empty() ||
                !shp->mat->normal_txt)
                continue;
            shp->tangsp.resize(shp->pos.size());
            ym::compute_tangent_frame(shp->triangles, shp->pos, shp->norm,
                shp->texcoord, shp->tangsp);
        }
    }
}

//
// Add missing data to the scene.
//
void add_radius(scene_group* scn, float radius) {
    for (auto msh : scn->meshes) {
        for (auto shp : msh->shapes) {
            if (shp->points.empty() && shp->lines.empty()) continue;
            if (!shp->radius.empty()) continue;
            shp->radius.resize(shp->pos.size(), radius);
        }
    }
}

//
// Add missing data to the scene.
//
void add_texture_data(scene_group* scn) {
    for (auto txt : scn->textures) {
        if (txt->dataf.empty() && txt->datab.empty()) {
            printf("unable to load texture %s\n", txt->path.c_str());
            txt->width = 1;
            txt->height = 1;
            txt->ncomp = 4;
            txt->datab = {{255, 255, 255, 255}};
        }
    }
}

//
// Add missing data to the scene.
//
void add_nodes(scene_group* scn) {
    if (!scn->nodes.empty()) return;
    for (auto mesh : scn->meshes) {
        auto ist = new node();
        ist->name = mesh->name;
        ist->msh = mesh;
        scn->nodes.push_back(ist);
    }
}

//
// Add missing data to the scene.
//
void add_scene(scene_group* scn) {
    if (!scn->scenes.empty()) return;
    auto s = new scene();
    s->name = "scene";
    update_transforms(scn);
    for (auto n : scn->nodes) {
        if (!n->parent) s->nodes.push_back(n);
    }
}

//
// Add missing data to the scene.
//
void add_names(scene_group* scn) {
    auto cid = 0;
    for (auto cam : scn->cameras) {
        if (cam->name.empty())
            cam->name = "<camera " + std::to_string(cid) + ">";
        cid++;
    }

    auto mid = 0;
    for (auto mat : scn->materials) {
        if (mat->name.empty())
            mat->name = "<material " + std::to_string(mid) + ">";
        mid++;
    }

    auto mmid = 0;
    for (auto mesh : scn->meshes) {
        if (mesh->name.empty())
            mesh->name = "<mesh " + std::to_string(mmid) + ">";
        mmid++;
        auto sid = 0;
        for (auto shp : mesh->shapes) {
            if (shp->name.empty())
                shp->name = "<shape " + std::to_string(sid) + ">";
            sid++;
        }
    }

    auto nid = 0;
    for (auto node : scn->nodes) {
        if (node->name.empty())
            node->name = "<node " + std::to_string(nid) + ">";
        nid++;
    }

    auto sid = 0;
    for (auto scene : scn->scenes) {
        if (scene->name.empty())
            scene->name = "<scene " + std::to_string(sid) + ">";
        sid++;
    }
}

//
// Add a default camera that views the entire scene.
//
void add_default_cameras(scene_group* scns) {
    for (auto scn : scns->scenes) {
        auto cams = get_camera_nodes(scn);
        if (cams.empty()) {
            // TODO: scene bounds
            auto bbox = ym::bbox3f{compute_scene_bounds(scns)};
            auto center = ym::center(bbox);
            auto bbox_size = ym::diagonal(bbox);
            auto bbox_msize =
                ym::max(bbox_size[0], ym::max(bbox_size[1], bbox_size[2]));
            // set up camera
            auto cam = new camera();
            cam->name = scn->name + " camera";
            auto camera_dir = ym::vec3f{1, 0.4f, 1};
            auto from = camera_dir * bbox_msize + center;
            auto to = center;
            auto up = ym::vec3f{0, 1, 0};
            cam->ortho = false;
            cam->aspect = 16.0f / 9.0f;
            cam->yfov = 2 * atanf(0.5f);
            cam->aperture = 0;
            cam->focus = ym::length(to - from);
            auto node = new ygltf::node();
            node->matrix = ym::to_mat(ym::lookat_frame3(from, to, up));
            node->cam = cam;
            node->name = cam->name;
            scns->cameras.push_back(cam);
            scns->nodes.push_back(node);
            scn->nodes.push_back(node);
        }
    }
}

#endif

//
// Update node trasforms
//
void update_transforms(node* ist) {
    ist->local_xform = node_transform(ist);
    if (ist->parent) {
        ist->xform = ist->parent->xform * ist->local_xform;
    } else {
        ist->xform = ist->local_xform;
    }
    for (auto child : ist->children) update_transforms(child);
}

//
// Update node trasforms
//
void update_transforms(scene_group* scns) {
    for (auto node : scns->nodes) node->parent = nullptr;
    for (auto node : scns->nodes) {
        for (auto child : node->children) child->parent = node;
    }
    for (auto node : scns->nodes) { update_transforms(node); }
}

//
// Evalute interpolated values
//
void update_animated_node_transforms(const animation* anim, float time) {
    time = ym::clamp(time, anim->time.front(), anim->time.back() - 0.001f);
    // get time slice
    auto i1 = 0, i2 = 0;
    auto t = 0.0f;
    auto interp = anim->interp;
    if (time <= anim->time.front()) {
        interp = animation_interpolation::step;
    } else if (time >= anim->time.back()) {
        i1 = anim->time.size() - 1;
        i2 = anim->time.size() - 2;
        interp = animation_interpolation::step;
    } else {
        for (i2 = 0; i2 < anim->time.size() && anim->time[i2] < time; i2++)
            ;
        i1 = i2 - 1;
        t = (time - anim->time[i1]) / (anim->time[i2] - anim->time[i1]);
    }

    // apply transforms
    if (!anim->translation.empty()) {
        auto trans = ym::vec3f{0, 0, 0};
        switch (interp) {
            case animation_interpolation::step: {
                trans = anim->translation[i1];
            } break;
            case animation_interpolation::linear: {
                trans =
                    anim->translation[i1] * (1 - t) + anim->translation[i2] * t;
            } break;
            case animation_interpolation::catmull_rom: {
            } break;
            case animation_interpolation::cubic: {
            } break;
        }
        for (auto node : anim->nodes) node->translation = trans;
    } else if (!anim->rotation.empty()) {
        auto rot = ym::quat4f{0, 0, 0, 1};
        switch (interp) {
            case animation_interpolation::step: {
                rot = anim->rotation[i1];
            } break;
            case animation_interpolation::linear: {
                // rot = ym::slerp(anim->rotation[i1], anim->rotation[i2], t);
                auto rot_ = (ym::quat4f)normalize(
                    ((ym::vec4f)anim->rotation[i1]) * (1 - t) +
                    ((ym::vec4f)anim->rotation[i2]) * t);
                rot = ym::quat4f{rot_.x, rot_.y, rot_.z, rot_.w};

            } break;
            case animation_interpolation::catmull_rom: {
            } break;
            case animation_interpolation::cubic: {
            } break;
        }
        for (auto node : anim->nodes) node->rotation = rot;
    } else if (!anim->scale.empty()) {
        auto scale = ym::vec3f{1, 1, 1};
        switch (interp) {
            case animation_interpolation::step: {
                scale = anim->scale[i1];
            } break;
            case animation_interpolation::linear: {
                scale = anim->scale[i1] * (1 - t) + anim->scale[i2] * t;
            } break;
            case animation_interpolation::catmull_rom: {
            } break;
            case animation_interpolation::cubic: {
            } break;
        }
        for (auto node : anim->nodes) node->scale = scale;
    } else if (!anim->morph_weights.empty()) {
        auto weights = std::vector<float>(anim->morph_weights[0].size());
        switch (interp) {
            case animation_interpolation::step: {
                weights = anim->morph_weights[i1];
            } break;
            case animation_interpolation::linear: {
                for (auto i = 0; i < weights.size(); i++) {
                    weights[i] = anim->morph_weights[i1][i] * (1 - t) +
                                 anim->morph_weights[i2][i] * t;
                }
            } break;
            case animation_interpolation::catmull_rom: {
            } break;
            case animation_interpolation::cubic: {
            } break;
        }
        for (auto node : anim->nodes) node->morph_weights = weights;
    } else {
    }
}

//
// Compute animation
//
void update_animated_transforms(scene_group* scns, float time) {
    for (auto anim_group : scns->animations) {
        for (auto anim : anim_group->animations) {
            update_animated_node_transforms(anim, time);
        }
    }
}

//
// Update skin trasforms
//
void update_skin_transforms(node* ist, node* parent) {
    ist->skin_xform = node_transform(ist);
    if (parent) ist->skin_xform = parent->skin_xform * ist->skin_xform;
    for (auto child : ist->children) update_skin_transforms(child, ist);
}

//
// Skin transforms (local-to-object)
//
std::vector<ym::mat4f> get_skin_transforms(
    const skin* sk, const ym::mat4f& xform) {
    auto ret = std::vector<ym::mat4f>(sk->joints.size());
    update_skin_transforms(sk->root, nullptr);
    auto inv_root = ym::inverse(xform);
    for (auto i = 0; i < sk->joints.size(); i++) {
        if (!sk->pose_matrices.empty()) {
            ret[i] = inv_root * sk->joints[i]->xform * sk->pose_matrices[i];
        } else {
            ret[i] = inv_root * sk->joints[i]->xform;
        }
    }
    return ret;
}

//
// Compute shape morphing
//
void compute_morphing_deformation(const shape* shp,
    const std::vector<float>& weights, std::vector<ym::vec3f>& pos,
    std::vector<ym::vec3f>& norm, std::vector<ym::vec4f>& tangsp) {
    pos = shp->pos;
    norm = shp->norm;
    tangsp = shp->tangsp;
    for (auto idx = 0; idx < shp->morph_targets.size(); idx++) {
        auto morph = shp->morph_targets[idx];
        auto weight = (idx < weights.size()) ? weights[idx] : morph->weight;
        if (weight == 0) continue;
        if (!morph->pos.empty()) {
            for (auto i = 0; i < pos.size(); i++) {
                pos[i] += weight * morph->pos[i];
            }
        }
        if (!morph->norm.empty()) {
            for (auto i = 0; i < pos.size(); i++) {
                norm[i] += weight * morph->norm[i];
            }
        }
        if (!morph->tangsp.empty()) {
            for (auto i = 0; i < tangsp.size(); i++) {
                *(ym::vec3f*)(&tangsp[i]) += weight * morph->tangsp[i];
            }
        }
    }
}

//
// Animation times
//
ym::vec2f get_animation_bounds(const scene_group* scns) {
    auto range = ym::vec2f{0, 0};
    for (auto anim_group : scns->animations) {
        for (auto anim : anim_group->animations) {
            range[0] = std::min(anim->time.front(), range[0]);
            range[1] = std::max(anim->time.back(), range[1]);
        }
    }
    return range;
}

//
// Helpder to filter nodes
//
void get_nodes(const std::vector<node*> nodes, bool has_camera, bool has_mesh,
    std::vector<node*>* filtered) {
    for (auto node : nodes) {
        if (has_camera && node->cam) filtered->push_back(node);
        if (has_mesh && node->msh) filtered->push_back(node);
        get_nodes(node->children, has_camera, has_mesh, filtered);
    }
}

//
// Helpder to filter nodes
//
std::vector<node*> get_nodes(
    const std::vector<node*> nodes, bool has_camera, bool has_mesh) {
    std::vector<node*> filtered;
    get_nodes(nodes, has_camera, has_mesh, &filtered);
    return filtered;
}

//
// Get a list of nodes with meshes
//
std::vector<node*> get_mesh_nodes(const scene* scn) {
    if (!scn) return {};
    return get_nodes(scn->nodes, false, true);
}

//
// Get a list of nodes with cameras
//
std::vector<node*> get_camera_nodes(const scene* scn) {
    if (!scn) return {};
    return get_nodes(scn->nodes, true, false);
}

//
// Support for buffer sections. Skip refcount on purpose.
//
static bool operator==(const buffer_section& a, const buffer_section& b) {
    return a.start == b.start && a.size == b.size && a.stride == b.stride &&
           a.count == b.count && a.type == b.type && a.ctype == b.ctype &&
           a.ncomp == b.ncomp && a.csize == b.csize;
}

//
// Generate buffer descriptions.
//
std::vector<buffer_descr*> gen_buffer_descriptors(const glTF* gltf) {
    std::vector<buffer_descr*> descriptors;

    auto bid = 0;
    for (auto buffer : gltf->buffers) {
        descriptors.push_back(new buffer_descr());
        auto descr = descriptors.back();
        descr->buffer = bid++;
        if (buffer->uri.substr(0, 5) == "data:") {
            descr->uri = "data:";
        } else {
            descr->uri = buffer->uri;
        }
        descr->size = buffer->byteLength;
        descr->name = buffer->name;
    }

    static const auto csize_map = std::map<glTFAccessorComponentType, int>{
        {glTFAccessorComponentType::NotSet, 0},
        {glTFAccessorComponentType::Byte, 1},
        {glTFAccessorComponentType::UnsignedByte, 1},
        {glTFAccessorComponentType::Short, 2},
        {glTFAccessorComponentType::UnsignedShort, 2},
        {glTFAccessorComponentType::UnsignedInt, 4},
        {glTFAccessorComponentType::Float, 4}};

    static const auto ncomp_map = std::map<glTFAccessorType, int>{
        {glTFAccessorType::NotSet, 0}, {glTFAccessorType::Scalar, 1},
        {glTFAccessorType::Vec2, 2}, {glTFAccessorType::Vec3, 3},
        {glTFAccessorType::Vec4, 4}, {glTFAccessorType::Mat2, 4},
        {glTFAccessorType::Mat3, 9}, {glTFAccessorType::Mat4, 16},
    };

    auto add_section = [gltf, &descriptors](glTFAccessor* accessor) {
        if (!accessor) return;
        auto bufferView = gltf->get(accessor->bufferView);
        if (!bufferView) return;
        auto buffer = gltf->get(bufferView->buffer);
        if (!buffer) return;
        auto descr = descriptors.at((int)bufferView->buffer);
        auto sect = buffer_section();
        sect.refcount = 1;
        sect.start = accessor->byteOffset + bufferView->byteOffset;
        sect.size = bufferView->byteLength - accessor->byteOffset;
        sect.stride = bufferView->byteStride;
        sect.count = accessor->count;
        sect.type = accessor->type;
        sect.ctype = accessor->componentType;
        sect.ncomp = ncomp_map.at(accessor->type);
        sect.csize = csize_map.at(accessor->componentType);
        auto found = false;
        for (auto dsect : descr->sections) {
            if (*dsect == sect) {
                dsect->refcount++;
                found = true;
                break;
            }
        }
        if (!found) descr->sections.push_back(new buffer_section(sect));
    };

    for (auto mesh : gltf->meshes) {
        for (auto prim : mesh->primitives) {
            for (auto attr : prim->attributes) {
                add_section(gltf->get(attr.second));
            }
            add_section(gltf->get(prim->indices));
        }
    }

    for (auto skin : gltf->skins) {
        add_section(gltf->get(skin->inverseBindMatrices));
    }

    for (auto anim : gltf->animations) {
        for (auto smp : anim->samplers) {
            add_section(gltf->get(smp->input));
            add_section(gltf->get(smp->output));
        }
    }

    for (auto descr : descriptors) {
        std::sort(descr->sections.data(),
            descr->sections.data() + descr->sections.size(),
            [](ygltf::buffer_section* a, ygltf::buffer_section* b) {
                return a->start < b->start;
            });
    }

    return descriptors;
}

}  // namespace ygltf
