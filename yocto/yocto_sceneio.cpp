//
// Implementation for Yocto/GL Input and Output functions.
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

#include "yocto_sceneio.h"
#include "yocto_imageio.h"
#include "yocto_random.h"
#include "yocto_shape.h"
#include "yocto_utils.h"

#include <cstdlib>
#include <regex>

#include <array>
#include <climits>

#include "ext/json.hpp"

#if YGL_HAPPLY
#include "ext/happly.h"
#endif

// -----------------------------------------------------------------------------
// JSON UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Json alias
using json = nlohmann::json;

// Helper for printing JSON values
bool print_value(string& str, const json& js) {
    str += js.dump();
    return true;
}

// Load a JSON object
bool load_json(const string& filename, json& js) {
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
bool save_json(const string& filename, const json& js) {
    auto str = ""s;
    try {
        str = js.dump(4);
    } catch (...) {
        log_io_error("could not dump json {}", filename);
        return false;
    }
    return save_text(filename, str);
}

inline void to_json(json& js, const vec2f& val) {
    js = std::array<float, 2>{{val.x, val.y}};
}
inline void from_json(const json& js, vec2f& val) {
    auto vala = js.get<std::array<float, 2>>();
    val       = {vala[0], vala[1]};
}
inline void to_json(json& js, const vec3f& val) {
    js = std::array<float, 3>{{val.x, val.y, val.z}};
}
inline void from_json(const json& js, vec3f& val) {
    auto vala = js.get<std::array<float, 3>>();
    val       = {vala[0], vala[1], vala[2]};
}
inline void to_json(json& js, const vec4f& val) {
    js = std::array<float, 4>{{val.x, val.y, val.z, val.w}};
}
inline void from_json(const json& js, vec4f& val) {
    auto vala = js.get<std::array<float, 4>>();
    val       = {vala[0], vala[1], vala[2], vala[3]};
}

inline void to_json(json& js, const vec2i& val) {
    js = std::array<int, 2>{{val.x, val.y}};
}
inline void from_json(const json& js, vec2i& val) {
    auto vala = js.get<std::array<int, 2>>();
    val       = {vala[0], vala[1]};
}
inline void to_json(json& js, const vec3i& val) {
    js = std::array<int, 3>{{val.x, val.y, val.z}};
}
inline void from_json(const json& js, vec3i& val) {
    auto vala = js.get<std::array<int, 3>>();
    val       = {vala[0], vala[1], vala[2]};
}
inline void to_json(json& js, const vec4i& val) {
    js = std::array<int, 4>{{val.x, val.y, val.z, val.w}};
}
inline void from_json(const json& js, vec4b& val) {
    auto vala = js.get<std::array<byte, 4>>();
    val       = {vala[0], vala[1], vala[2], vala[3]};
}
inline void to_json(json& js, const vec4b& val) {
    js = std::array<byte, 4>{{val.x, val.y, val.z, val.w}};
}
inline void from_json(const json& js, vec4i& val) {
    auto vala = js.get<std::array<int, 4>>();
    val       = {vala[0], vala[1], vala[2], vala[3]};
}

inline void to_json(json& js, const frame3f& val) {
    js = std::array<vec3f, 4>{{val.x, val.y, val.z, val.o}};
}
inline void from_json(const json& js, frame3f& val) {
    auto vala = js.get<std::array<vec3f, 4>>();
    val       = {vala[0], vala[1], vala[2], vala[3]};
}

inline void to_json(json& js, const mat4f& val) {
    js = std::array<vec4f, 4>{{val.x, val.y, val.z, val.w}};
}
inline void from_json(const json& js, mat4f& val) {
    auto vala = js.get<std::array<vec4f, 4>>();
    val       = {vala[0], vala[1], vala[2], vala[3]};
}

inline void to_json(json& js, const bbox3f& val) {
    js = std::array<vec3f, 2>{{val.min, val.max}};
}
inline void from_json(const json& js, bbox3f& val) {
    auto vala = js.get<std::array<vec3f, 2>>();
    val       = {vala[0], vala[1]};
}

inline void to_json(json& js, const image4f& value) {
    js           = json::object();
    js["width"]  = value.width;
    js["height"] = value.height;
    js["pixels"] = vector<vec4f>{
        data(value), data(value) + value.width * value.height};
}
inline void to_json(json& js, const image4b& value) {
    js           = json::object();
    js["width"]  = value.width;
    js["height"] = value.height;
    js["pixels"] = vector<vec4b>{
        data(value), data(value) + value.width * value.height};
}
inline void from_json(const json& js, image4f& value) {
    auto width  = js.at("width").get<int>();
    auto height = js.at("height").get<int>();
    auto pixels = js.at("pixels").get<vector<vec4f>>();
    value       = make_image(width, height, data(pixels));
}
inline void from_json(const json& js, image4b& value) {
    auto width  = js.at("width").get<int>();
    auto height = js.at("height").get<int>();
    auto pixels = js.at("pixels").get<vector<vec4b>>();
    value       = make_image(width, height, data(pixels));
}
inline void to_json(json& js, const volume1f& value) {
    js           = json::object();
    js["width"]  = value.width;
    js["height"] = value.height;
    js["depth"]  = value.width;
    js["voxels"] = vector<float>{
        data(value), data(value) + value.width * value.height * value.depth};
}
inline void from_json(const json& js, volume1f& value) {
    auto width  = js.at("width").get<int>();
    auto height = js.at("height").get<int>();
    auto depth  = js.at("depth").get<int>();
    auto voxels = js.at("voxels").get<vector<float>>();
    value       = make_volume(width, height, depth, data(voxels));
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GENERIC SCENE LOADING
// -----------------------------------------------------------------------------
namespace yocto {

// Load a scene
bool load_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options) {
    auto ext = get_extension(filename);
    if (ext == "json" || ext == "JSON") {
        return load_json_scene(filename, scene, options);
    } else if (ext == "obj" || ext == "OBJ") {
        return load_obj_scene(filename, scene, options);
    } else if (ext == "gltf" || ext == "GLTF") {
        return load_gltf_scene(filename, scene, options);
    } else if (ext == "pbrt" || ext == "PBRT") {
        return load_pbrt_scene(filename, scene, options);
    } else if (ext == "ybin" || ext == "YBIN") {
        return load_ybin_scene(filename, scene, options);
    } else {
        scene = {};
        log_io_error("unsupported scene format {}", ext);
        return false;
    }
}

// Save a scene
bool save_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options) {
    auto ext = get_extension(filename);
    if (ext == "json" || ext == "JSON") {
        return save_json_scene(filename, scene, options);
    } else if (ext == "obj" || ext == "OBJ") {
        return save_obj_scene(filename, scene, options);
    } else if (ext == "gltf" || ext == "GLTF") {
        return save_gltf_scene(filename, scene, options);
    } else if (ext == "pbrt" || ext == "PBRT") {
        return save_pbrt_scene(filename, scene, options);
    } else if (ext == "ybin" || ext == "YBIN") {
        return save_ybin_scene(filename, scene, options);
    } else {
        log_io_error("unsupported scene format {}", ext);
        return false;
    }
}

bool load_image_nolog(const string& filename, image4f& img);
bool load_image_nolog(const string& filename, image4b& img);
bool load_volume_nolog(const string& filename, volume1f& vol);

bool load_scene_textures(yocto_scene& scene, const string& dirname,
    const load_scene_options& options) {
    if (options.skip_textures) return true;

    // load images
    atomic<bool> exit_error(false);
    parallel_foreach(scene.textures,
        [&exit_error, &options, &dirname](yocto_texture& texture) {
            if (exit_error) return;
            if (texture.filename == "" || !empty(texture.hdr_image) ||
                !empty(texture.ldr_image))
                return;
            auto filename = normalize_path(dirname + texture.filename);
            if (is_hdr_filename(filename)) {
                if (!load_image_nolog(filename, texture.hdr_image)) {
                    if (options.exit_on_error) {
                        exit_error = true;
                        return;
                    }
                }
            } else {
                if (!load_image_nolog(filename, texture.ldr_image)) {
                    if (options.exit_on_error) {
                        exit_error = true;
                        return;
                    }
                }
            }
        },
        options.cancel_flag, options.run_serially);
    if (exit_error) return false;

    // load volumes
    parallel_foreach(scene.voltextures,
        [&exit_error, &options, &dirname](yocto_voltexture& texture) {
            if (exit_error) return;
            if (texture.filename == "" || !empty(texture.volume_data)) return;
            auto filename = normalize_path(dirname + texture.filename);
            if (!load_volume_nolog(filename, texture.volume_data)) {
                if (options.exit_on_error) {
                    exit_error = true;
                    return;
                }
            }
        },
        options.cancel_flag, options.run_serially);
    if (exit_error) return false;

    // assign opacity texture if needed
    if (options.assign_texture_opacity) {
        auto has_opacity = vector<int>(scene.textures.size());
        parallel_for((int)scene.textures.size(),
            [&scene, &has_opacity](int texture_id) {
                auto& texture           = scene.textures[texture_id];
                has_opacity[texture_id] = 0;
                for (auto& p : texture.hdr_image) {
                    if (p.w < 0.999f) {
                        has_opacity[texture_id] = 1;
                        break;
                    }
                }
                for (auto& p : texture.ldr_image) {
                    if (p.w < 255) {
                        has_opacity[texture_id] = 1;
                        break;
                    }
                }
            },
            options.cancel_flag, options.run_serially);
        for (auto& material : scene.materials) {
            if (material.diffuse_texture >= 0 && material.opacity_texture < 0 &&
                has_opacity[material.diffuse_texture])
                material.opacity_texture = material.diffuse_texture;
        }
    }

    // done
    return true;
}

bool save_image_nolog(const string& filename, const image4f& img);
bool save_image_nolog(const string& filename, const image4b& img);
bool save_volume_nolog(const string& filename, const volume1f& vol);

// helper to save textures
bool save_scene_textures(const yocto_scene& scene, const string& dirname,
    const save_scene_options& options) {
    if (options.skip_textures) return true;

    // save images
    atomic<bool> exit_error(false);
    parallel_foreach(scene.textures,
        [&exit_error, &options, &dirname](const yocto_texture& texture) {
            if (exit_error) return;
            if (empty(texture.hdr_image) && empty(texture.ldr_image)) return;
            auto filename = normalize_path(dirname + texture.filename);
            if (is_hdr_filename(filename)) {
                if (!save_image_nolog(filename, texture.hdr_image)) {
                    if (options.exit_on_error) {
                        exit_error = true;
                        return;
                    }
                }
            } else {
                if (!save_image_nolog(filename, texture.ldr_image)) {
                    if (options.exit_on_error) {
                        exit_error = true;
                        return;
                    }
                }
            }
        },
        options.cancel_flag, options.run_serially);
    if (exit_error) return false;

    // save volumes
    parallel_foreach(scene.voltextures,
        [&exit_error, &options, &dirname](const yocto_voltexture& texture) {
            if (exit_error) return;
            if (empty(texture.volume_data)) return;
            auto filename = normalize_path(dirname + texture.filename);
            if (!save_volume_nolog(filename, texture.volume_data)) {
                if (options.exit_on_error) {
                    exit_error = true;
                    return;
                }
            }
        },
        options.cancel_flag, options.run_serially);
    if (exit_error) return false;

    // done
    return true;
}

// merge quads and triangles
void merge_triangles_and_quads(
    vector<vec3i>& triangles, vector<vec4i>& quads, bool force_triangles) {
    if (empty(quads)) return;
    if (force_triangles) {
        auto qtriangles = convert_quads_to_triangles(quads);
        triangles.insert(triangles.end(), qtriangles.begin(), qtriangles.end());
        quads = {};
    } else {
        auto tquads = convert_triangles_to_quads(triangles);
        quads.insert(quads.end(), tquads.begin(), tquads.end());
        triangles = {};
    }
}

// check if it is really face varying
bool is_face_varying(const vector<vec4i>& quads_positions,
    const vector<vec4i>& quads_normals, const vector<vec4i>& quads_texcoords) {
    if (empty(quads_positions)) return false;
    if (!empty(quads_normals)) {
        for (auto i = 0; i < quads_positions.size(); i++)
            if (quads_positions[i] != quads_normals[i]) return true;
    }
    if (!empty(quads_texcoords)) {
        for (auto i = 0; i < quads_positions.size(); i++)
            if (quads_positions[i] != quads_texcoords[i]) return true;
    }
    return false;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IO UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Encode in base64
string base64_encode(unsigned char const* bytes_to_encode, unsigned int in_len) {
    static const string base64_chars =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz"
        "0123456789+/";

    string        ret;
    int           i = 0;
    int           j = 0;
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
        char_array_4[1] = ((char_array_3[0] & 0x03) << 4) +
                          ((char_array_3[1] & 0xf0) >> 4);
        char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) +
                          ((char_array_3[2] & 0xc0) >> 6);
        char_array_4[3] = char_array_3[2] & 0x3f;

        for (j = 0; (j < i + 1); j++) ret += base64_chars[char_array_4[j]];

        while ((i++ < 3)) ret += '=';
    }

    return ret;
}

// Decode from base64
string base64_decode(string const& encoded_string) {
    static const string base64_chars =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz"
        "0123456789+/";

    auto is_base64 = [](unsigned char c) -> bool {
        return (isalnum(c) || (c == '+') || (c == '/'));
    };

    int           in_len = (int)encoded_string.size();
    int           i      = 0;
    int           j      = 0;
    int           in_    = 0;
    unsigned char char_array_4[4], char_array_3[3];
    string        ret;

    while (in_len-- && (encoded_string[in_] != '=') &&
           is_base64(encoded_string[in_])) {
        char_array_4[i++] = encoded_string[in_];
        in_++;
        if (i == 4) {
            for (i = 0; i < 4; i++)
                char_array_4[i] = base64_chars.find(char_array_4[i]);

            char_array_3[0] = (char_array_4[0] << 2) +
                              ((char_array_4[1] & 0x30) >> 4);
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

        char_array_3[0] = (char_array_4[0] << 2) +
                          ((char_array_4[1] & 0x30) >> 4);
        char_array_3[1] = ((char_array_4[1] & 0xf) << 4) +
                          ((char_array_4[2] & 0x3c) >> 2);
        char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

        for (j = 0; (j < i - 1); j++) ret += char_array_3[j];
    }

    return ret;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// BUILTIN JSON FORMAT
// -----------------------------------------------------------------------------
namespace yocto {

bool operator==(const image4f& a, const image4f& b) {
    return a.width == b.width && a.height == b.height && a.pixels == b.pixels;
}
bool operator==(const image4b& a, const image4b& b) {
    return a.width == b.width && a.height == b.height && a.pixels == b.pixels;
}
bool operator==(const volume1f& a, const volume1f& b) {
    return a.width == b.width && a.height == b.height && a.depth == b.depth &&
           a.voxels == b.voxels;
}

// Dumps a json value
bool serialize_json_value(json& js, int& value, bool save) {
    if (save) {
        js = value;
        return true;
    } else {
        if (!js.is_number_integer()) return false;
        value = js.get<int>();
        return true;
    }
}
bool serialize_json_value(json& js, bool& value, bool save) {
    if (save) {
        js = value;
        return true;
    } else {
        if (!js.is_boolean()) return false;
        value = js.get<bool>();
        return true;
    }
}
bool serialize_json_value(json& js, unsigned char& value, bool save) {
    if (save) {
        js = (int)value;
        return true;
    } else {
        if (!js.is_number_integer()) return false;
        value = (unsigned char)js.get<int>();
        return true;
    }
}
bool serialize_json_value(json& js, float& value, bool save) {
    if (save) {
        js = value;
        return true;
    } else {
        if (!js.is_number()) return false;
        value = js.get<float>();
        return true;
    }
}
bool serialize_json_value(json& js, double& value, bool save) {
    if (save) {
        js = value;
        return true;
    } else {
        if (!js.is_number()) return false;
        value = js.get<float>();
        return true;
    }
}
bool serialize_json_value(json& js, string& value, bool save) {
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
bool serialize_json_values(json& js, T* values, int num, bool save) {
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
bool serialize_json_value(json& js, vector<T>& value, bool save) {
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

bool serialize_json_value(json& js, vec2f& value, bool save) {
    return serialize_json_values(js, &value.x, 2, save);
}
bool serialize_json_value(json& js, vec3f& value, bool save) {
    return serialize_json_values(js, &value.x, 3, save);
}
bool serialize_json_value(json& js, vec4f& value, bool save) {
    return serialize_json_values(js, &value.x, 4, save);
}
bool serialize_json_value(json& js, vec2i& value, bool save) {
    return serialize_json_values(js, &value.x, 2, save);
}
bool serialize_json_value(json& js, vec3i& value, bool save) {
    return serialize_json_values(js, &value.x, 3, save);
}
bool serialize_json_value(json& js, vec4i& value, bool save) {
    return serialize_json_values(js, &value.x, 4, save);
}
bool serialize_json_value(json& js, vec4b& value, bool save) {
    return serialize_json_values(js, &value.x, 4, save);
}
bool serialize_json_value(json& js, mat2f& value, bool save) {
    return serialize_json_values(js, &value.x.x, 4, save);
}
bool serialize_json_value(json& js, mat3f& value, bool save) {
    return serialize_json_values(js, &value.x.x, 9, save);
}
bool serialize_json_value(json& js, mat4f& value, bool save) {
    return serialize_json_values(js, &value.x.x, 16, save);
}
bool serialize_json_value(json& js, frame2f& value, bool save) {
    return serialize_json_values(js, &value.x.x, 6, save);
}
bool serialize_json_value(json& js, frame3f& value, bool save) {
    return serialize_json_values(js, &value.x.x, 12, save);
}
bool serialize_json_value(json& js, bbox3f& value, bool save) {
    return serialize_json_values(js, &value.min.x, 6, save);
}

// Dumps a json value
template <typename T>
bool serialize_json_value(
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

// Dumps a json value
template <typename T>
bool serialize_json_values(
    json& js, T* values, int num, const char* name, bool save) {
    if (save) {
        if (!values || num == 0) return true;
        return serialize_json_values(js[name], values, num, save);
    } else {
        if (!js.count(name)) return true;
        return serialize_json_values(js.at(name), values, num, save);
    }
}

bool serialize_json_value(json& js, image4f& value, bool save) {
    auto width = 0, height = 0;
    if (!serialize_json_value(js, width, "width", -1, save)) return false;
    if (!serialize_json_value(js, height, "height", -1, save)) return false;
    if (!save) value = make_image(width, height, zero4f);
    if (!serialize_json_values(js, data(value), width * height, "pixels", save))
        return false;
    return true;
}
bool serialize_json_value(json& js, image4b& value, bool save) {
    auto width = 0, height = 0;
    if (!serialize_json_value(js, width, "width", -1, save)) return false;
    if (!serialize_json_value(js, height, "height", -1, save)) return false;
    if (!save) value = make_image(width, height, zero4b);
    if (!serialize_json_values(js, data(value), width * height, "pixels", save))
        return false;
    return true;
}
template <typename T>
bool serialize_json_value(json& js, volume1f& value, bool save) {
    auto width = 0, height = 0, depth = 0;
    if (!serialize_json_value(js, width, "width", -1, save)) return false;
    if (!serialize_json_value(js, height, "height", -1, save)) return false;
    if (!serialize_json_value(js, depth, "heidepthght", -1, save)) return false;
    if (!save) value = make_volume(width, height, depth, 0.0f);
    if (!serialize_json_values(
            js, data(value), width * height * depth, "voxels", save))
        return false;
    return true;
}

// Dumps a json value
template <typename T>
bool serialize_json_objref(
    json& js, int& value, const vector<T>& refs, bool save) {
    if (save) {
        auto vals = value >= 0 ? refs[value].name : ""s;
        return serialize_json_value(js, vals, save);
    } else {
        if (!js.is_string()) return false;
        auto name = ""s;
        if (!serialize_json_value(js, name, save)) return false;
        value = -1;
        for (auto index = 0; index < refs.size(); index++) {
            if (refs[index].name == name) {
                value = index;
                break;
            }
        }
        if (value < 0) return false;
        return true;
    }
}

// Dumps a json value
template <typename T>
bool serialize_json_objref(
    json& js, vector<int>& value, const vector<T>& refs, bool save) {
    if (save) {
        js = json::array();
        for (auto v : value) {
            js.push_back({});
            if (!serialize_json_objref(js.back(), v, refs, save)) return false;
        }
        return true;
    } else {
        if (!js.is_array()) return false;
        for (auto& j : js) {
            value.push_back(-1);
            if (!serialize_json_objref(j, value.back(), refs, save))
                return false;
        }
        return true;
    }
}

// Dumps a json value
template <typename T>
bool serialize_json_objref(
    json& js, int& value, const char* name, const vector<T>& refs, bool save) {
    if (save) {
        if (value < 0) return true;
        return serialize_json_objref(js[name], value, refs, save);
    } else {
        if (!js.count(name)) return true;
        value = -1;
        return serialize_json_objref(js.at(name), value, refs, save);
    }
}

// Dumps a json value
template <typename T>
bool serialize_json_objref(json& js, vector<int>& value, const char* name,
    const vector<T>& refs, bool save) {
    if (save) {
        if (empty(value)) return true;
        return serialize_json_objref(js[name], value, refs, save);
    } else {
        if (!js.count(name)) return true;
        value = {};
        return serialize_json_objref(js.at(name), value, refs, save);
    }
}

// Starts a json object
bool serialize_json_objbegin(json& js, bool save) {
    if (save) {
        js = json::object();
        return true;
    } else {
        return js.is_object();
    }
}

// Dumps a json value
template <typename T>
bool serialize_json_objarray(
    json& js, vector<T>& value, const yocto_scene& scene, bool save) {
    if (save) {
        js = json::array();
        for (auto& v : value) {
            js.push_back({});
            if (!serialize_json_object(js.back(), v, scene, save)) return false;
        }
        return true;
    } else {
        if (!js.is_array()) return false;
        for (auto& j : js) {
            value.push_back(T{});
            if (!serialize_json_object(j, value.back(), scene, save))
                return false;
        }
        return true;
    }
}

// Dumps a json value
template <typename T>
bool serialize_json_objarray(json& js, vector<T>& value, const char* name,
    const yocto_scene& scene, bool save) {
    if (save) {
        if (empty(value)) return true;
        return serialize_json_objarray(js[name], value, scene, save);
    } else {
        if (!js.count(name)) return true;
        value = {};
        return serialize_json_objarray(js.at(name), value, scene, save);
    }
}

// Parses and applied a JSON procedural
template <typename T>
bool serialize_json_procedural(const json& js, T& value, const char* name,
    const yocto_scene& scene, bool save) {
    if (save) {
        return true;
    } else {
        if (!js.count(name)) return true;
        return apply_json_procedural(js.at(name), value, scene);
    }
}

// Procedural commands for cameras
bool apply_json_procedural(
    const json& js, yocto_camera& value, const yocto_scene& scene) {
    if (!serialize_json_objbegin((json&)js, false)) return false;
    if (js.count("from") || js.count("to")) {
        auto from            = js.value("from", zero3f);
        auto to              = js.value("to", zero3f);
        auto up              = js.value("up", vec3f{0, 1, 0});
        value.frame          = make_lookat_frame(from, to, up);
        value.focus_distance = length(from - to);
    }
    return true;
}

// Serialize struct
bool serialize_json_object(
    json& js, yocto_camera& value, const yocto_scene& scene, bool save) {
    static const auto def = yocto_camera();
    if (!serialize_json_objbegin(js, save)) return false;
    if (!serialize_json_value(js, value.name, "name", def.name, save))
        return false;
    if (!serialize_json_value(js, value.frame, "frame", def.frame, save))
        return false;
    if (!serialize_json_value(
            js, value.orthographic, "orthographic", def.orthographic, save))
        return false;
    if (!serialize_json_value(
            js, value.film_width, "film_width", def.film_width, save))
        return false;
    if (!serialize_json_value(
            js, value.film_height, "film_height", def.film_height, save))
        return false;
    if (!serialize_json_value(
            js, value.focal_length, "focal_length", def.focal_length, save))
        return false;
    if (!serialize_json_value(js, value.focus_distance, "focus_distance",
            def.focus_distance, save))
        return false;
    if (!serialize_json_value(
            js, value.lens_aperture, "lens_aperture", def.lens_aperture, save))
        return false;
    if (!serialize_json_procedural(js, value, "!!proc", scene, save))
        return false;
    return true;
}

// Procedural commands for textures
bool apply_json_procedural(
    const json& js, yocto_texture& value, const yocto_scene& scene) {
    if (!serialize_json_objbegin((json&)js, false)) return false;
    auto type = js.value("type", ""s);
    if (type == "") return true;
    auto is_hdr = false;
    auto width  = js.value("width", 512);
    auto height = js.value("height", 512);
    if (type == "grid") {
        value.hdr_image = make_grid_image(width, height, js.value("tile", 8),
            js.value("c0", vec4f{0.2f, 0.2f, 0.2f, 1}),
            js.value("c1", vec4f{0.8f, 0.8f, 0.8f, 1}));
    } else if (type == "checker") {
        value.hdr_image = make_checker_image(width, height, js.value("tile", 8),
            js.value("c0", vec4f{0.2f, 0.2f, 0.2f, 1}),
            js.value("c1", vec4f{0.8f, 0.8f, 0.8f, 1}));
    } else if (type == "bump") {
        value.hdr_image = make_bumpdimple_image(
            width, height, js.value("tile", 8));
    } else if (type == "uvramp") {
        value.hdr_image = make_uvramp_image(width, height);
    } else if (type == "gammaramp") {
        value.hdr_image = make_gammaramp_image(width, height);
    } else if (type == "blackbodyramp") {
        value.hdr_image = make_blackbodyramp_image(width, height);
    } else if (type == "uvgrid") {
        value.hdr_image = make_uvgrid_image(width, height);
    } else if (type == "sky") {
        if (width < height * 2) width = height * 2;
        value.hdr_image = make_sunsky_image(width, height,
            js.value("sun_angle", pif / 4), js.value("turbidity", 3.0f),
            js.value("has_sun", false), js.value("sun_angle_scale", 1.0f),
            js.value("sun_emission_scale", 1.0f),
            js.value("ground_albedo", vec3f{0.7f, 0.7f, 0.7f}));
        is_hdr          = true;
    } else if (type == "noise") {
        value.hdr_image = make_noise_image(
            width, height, js.value("scale", 1.0f), js.value("wrap", true));
    } else if (type == "fbm") {
        value.hdr_image = make_fbm_image(width, height, js.value("scale", 1.0f),
            js.value("lacunarity", 2.0f), js.value("gain", 0.5f),
            js.value("octaves", 6), js.value("wrap", true));
    } else if (type == "ridge") {
        value.hdr_image = make_ridge_image(width, height,
            js.value("scale", 1.0f), js.value("lacunarity", 2.0f),
            js.value("gain", 0.5f), js.value("offset", 1.0f),
            js.value("octaves", 6), js.value("wrap", true));
    } else if (type == "turbulence") {
        value.hdr_image = make_turbulence_image(width, height,
            js.value("scale", 1.0f), js.value("lacunarity", 2.0f),
            js.value("gain", 0.5f), js.value("octaves", 6),
            js.value("wrap", true));
    } else {
        log_error("unknown texture type {}", type);
        return false;
    }
    if (js.value("bump_to_normal", false)) {
        value.hdr_image = bump_to_normal_map(
            value.hdr_image, js.value("bump_scale", 1.0f));
        value.ldr_as_linear = true;
    }
    if (!is_hdr) {
        if (!value.ldr_as_linear) {
            value.ldr_image = float_to_byte(linear_to_srgb(value.hdr_image));
        } else {
            value.ldr_image = float_to_byte(value.hdr_image);
        }
        value.hdr_image = {};
    }
    if (value.filename == "") {
        auto ext       = (is_hdr) ? string("hdr") : string("png");
        value.filename = "textures/" + value.name + "." + ext;
    }
    return true;
}

// Serialize struct
bool serialize_json_object(
    json& js, yocto_texture& value, const yocto_scene& scene, bool save) {
    static const auto def = yocto_texture();
    if (!serialize_json_objbegin(js, save)) return false;
    if (!serialize_json_value(js, value.name, "name", def.name, save))
        return false;
    if (!serialize_json_value(js, value.filename, "filename", def.filename, save))
        return false;
    if (!serialize_json_value(
            js, value.clamp_to_edge, "clamp_to_edge", def.clamp_to_edge, save))
        return false;
    if (!serialize_json_value(
            js, value.height_scale, "height_scale", def.height_scale, save))
        return false;
    if (!serialize_json_value(js, value.no_interpolation, "no_interpolation",
            def.no_interpolation, save))
        return false;
    if (!serialize_json_value(
            js, value.ldr_as_linear, "ldr_as_linear", def.ldr_as_linear, save))
        return false;
    if (value.filename == "" || !save) {
        if (!serialize_json_value(
                js, value.hdr_image, "hdr_image", def.hdr_image, save))
            return false;
        if (!serialize_json_value(
                js, value.ldr_image, "ldr_image", def.ldr_image, save))
            return false;
    }
    if (!serialize_json_procedural(js, value, "!!proc", scene, save))
        return false;
    return true;
}

// Procedural commands for textures
bool apply_json_procedural(
    const json& js, yocto_voltexture& value, const yocto_scene& scene) {
    if (!serialize_json_objbegin((json&)js, false)) return false;
    auto type = js.value("type", ""s);
    if (type == "") return true;
    auto width  = js.value("width", 512);
    auto height = js.value("height", 512);
    auto depth  = js.value("depth", 512);
    if (type == "test_volume") {
        value.volume_data = make_test_volume(width, height, depth,
            js.value("scale", 10.0f), js.value("exponent", 6.0f));
    } else {
        log_error("unknown texture type {}", type);
        return false;
    }
    if (value.filename == "") {
        auto ext       = string("vol");
        value.filename = "textures/" + value.name + "." + ext;
    }
    return true;
}

// Serialize struct
bool serialize_json_object(
    json& js, yocto_voltexture& value, const yocto_scene& scene, bool save) {
    static const auto def = yocto_voltexture();
    if (!serialize_json_objbegin(js, save)) return false;
    if (!serialize_json_value(js, value.name, "name", def.name, save))
        return false;
    if (!serialize_json_value(js, value.filename, "filename", def.filename, save))
        return false;
    if (!serialize_json_value(
            js, value.clamp_to_edge, "clamp_to_edge", def.clamp_to_edge, save))
        return false;
    if (!serialize_json_value(js, value.no_interpolation, "no_interpolation",
            def.no_interpolation, save))
        return false;
    if (value.filename == "") {
        if (!empty(value.volume_data)) js["vol"] = value.volume_data;
    }
    if (!serialize_json_procedural(js, value, "!!proc", scene, save))
        return false;
    return true;
}

// Procedural commands for materials
bool apply_json_procedural(
    const json& js, yocto_material& value, const yocto_scene& scene) {
    if (!serialize_json_objbegin((json&)js, false)) return false;
    return true;
}

// Serialize struct
bool serialize_json_object(
    json& js, yocto_material& value, const yocto_scene& scene, bool save) {
    static const auto def = yocto_material();
    if (!serialize_json_objbegin(js, save)) return false;
    if (!serialize_json_value(js, value.name, "name", def.name, save))
        return false;
    if (value.base_metallic != def.base_metallic)
        js["base_metallic"] = value.base_metallic;
    if (value.gltf_textures != def.gltf_textures)
        js["gltf_textures"] = value.gltf_textures;
    if (value.double_sided != def.double_sided)
        js["double_sided"] = value.double_sided;
    if (!serialize_json_value(js, value.emission, "emission", def.emission, save))
        return false;
    if (!serialize_json_value(js, value.diffuse, "diffuse", def.diffuse, save))
        return false;
    if (!serialize_json_value(js, value.specular, "specular", def.specular, save))
        return false;
    if (!serialize_json_value(
            js, value.transmission, "transmission", def.transmission, save))
        return false;
    if (!serialize_json_value(
            js, value.roughness, "roughness", def.roughness, save))
        return false;
    if (!serialize_json_value(js, value.opacity, "opacity", def.opacity, save))
        return false;
    if (!serialize_json_value(js, value.fresnel, "fresnel", def.fresnel, save))
        return false;
    if (!serialize_json_value(js, value.refract, "refract", def.refract, save))
        return false;
    if (!serialize_json_objref(js, value.emission_texture, "emission_texture",
            scene.textures, save))
        return false;
    if (!serialize_json_objref(
            js, value.diffuse_texture, "diffuse_texture", scene.textures, save))
        return false;
    if (!serialize_json_objref(js, value.specular_texture, "specular_texture",
            scene.textures, save))
        return false;
    if (!serialize_json_objref(js, value.transmission_texture,
            "transmission_texture", scene.textures, save))
        return false;
    if (!serialize_json_objref(js, value.roughness_texture, "roughness_texture",
            scene.textures, save))
        return false;
    if (!serialize_json_objref(
            js, value.opacity_texture, "opacity_texture", scene.textures, save))
        return false;
    if (!serialize_json_objref(js, value.occlusion_texture, "occlusion_texture",
            scene.textures, save))
        return false;
    if (!serialize_json_objref(
            js, value.bump_texture, "bump_texture", scene.textures, save))
        return false;
    if (!serialize_json_objref(js, value.displacement_texture,
            "displacement_texture", scene.textures, save))
        return false;
    if (!serialize_json_objref(
            js, value.normal_texture, "normal_texture", scene.textures, save))
        return false;
    if (!serialize_json_objref(js, value.volume_density_texture,
            "volume_density_texture", scene.voltextures, save))
        return false;
    if (!serialize_json_procedural(js, value, "!!proc", scene, save))
        return false;
    return true;
}

// Procedural commands for materials
bool apply_json_procedural(
    const json& js, yocto_shape& value, const yocto_scene& scene) {
    if (!serialize_json_objbegin((json&)js, false)) return false;
    auto type = js.value("type", ""s);
    if (type == "") return true;
    value.points        = {};
    value.lines         = {};
    value.triangles     = {};
    value.quads         = {};
    value.positions     = {};
    value.normals       = {};
    value.texturecoords = {};
    value.radius        = {};
    if (type == "quad") {
        tie(value.quads, value.positions, value.normals,
            value.texturecoords) = make_quad_shape(js.value("steps", vec2i{1, 1}),
            js.value("size", vec2f{2, 2}), js.value("uvsize", vec2f{1, 1}));
    } else if (type == "quady") {
        tie(value.quads, value.positions, value.normals,
            value.texturecoords) = make_quad_shape(js.value("steps", vec2i{1, 1}),
            js.value("size", vec2f{2, 2}), js.value("uvsize", vec2f{1, 1}));
    } else if (type == "quad_stack") {
        tie(value.quads, value.positions, value.normals,
            value.texturecoords) = make_quad_stack_shape(js.value("steps",
                                                             vec3i{1, 1, 1}),
            js.value("size", vec3f{2, 2, 2}), js.value("uvsize", vec2f{1, 1}));
    } else if (type == "cube") {
        tie(value.quads, value.positions, value.normals,
            value.texturecoords) = make_cube_shape(js.value("steps",
                                                       vec3i{1, 1, 1}),
            js.value("size", vec3f{2, 2, 2}), js.value("uvsize", vec3f{1, 1, 1}));
    } else if (type == "cube_rounded") {
        tie(value.quads, value.positions, value.normals,
            value.texturecoords) = make_cube_rounded_shape(js.value("steps",
                                                               vec3i{32, 32, 32}),
            js.value("size", vec3f{2, 2, 2}),
            js.value("uvsize", vec3f{1, 1, 1}), js.value("radius", 0.3f));
    } else if (type == "sphere") {
        tie(value.quads, value.positions, value.normals,
            value.texturecoords) = make_sphere_shape(js.value("steps",
                                                         vec2i{64, 32}),
            js.value("size", 2.0f), js.value("uvsize", vec2f{1, 1}));
    } else if (type == "sphere_cube") {
        tie(value.quads, value.positions, value.normals,
            value.texturecoords) = make_sphere_cube_shape(js.value("steps", 32),
            js.value("size", 2.0f), js.value("uvsize", 1.0f));
    } else if (type == "sphere_flipcap") {
        tie(value.quads, value.positions, value.normals,
            value.texturecoords) = make_sphere_flipcap_shape(js.value("steps",
                                                                 vec2i{64, 32}),
            js.value("size", 2.0f), js.value("uvsize", vec2f{1, 1}),
            js.value("zflip", vec2f{-0.75f, +0.75f}));
    } else if (type == "disk") {
        tie(value.quads, value.positions, value.normals,
            value.texturecoords) = make_disk_shape(js.value("steps",
                                                       vec2i{32, 16}),
            js.value("size", 2.0f), js.value("uvsize", vec2f{1, 1}));
    } else if (type == "disk_quad") {
        tie(value.quads, value.positions, value.normals,
            value.texturecoords) = make_disk_quad_shape(js.value("steps", 32),
            js.value("size", 2.0f), js.value("uvsize", 1.0f));
    } else if (type == "disk_bulged") {
        tie(value.quads, value.positions, value.normals,
            value.texturecoords) = make_disk_bulged_shape(js.value("steps", 32),
            js.value("size", 2.0f), js.value("uvsize", 1.0f),
            js.value("height", 0.25f));
    } else if (type == "cylinder_side") {
        tie(value.quads, value.positions, value.normals,
            value.texturecoords) = make_cylinder_side_shape(js.value("steps",
                                                                vec2i{64, 32}),
            js.value("size", vec2f{2.0f, 2.0f}), js.value("uvsize", vec2f{1, 1}));
    } else if (type == "cylinder") {
        tie(value.quads, value.positions, value.normals,
            value.texturecoords) = make_cylinder_shape(js.value("steps",
                                                           vec3i{64, 32, 16}),
            js.value("size", vec2f{2.0f, 2.0f}),
            js.value("uvsize", vec3f{1, 1, 1}));
    } else if (type == "cylinder_rounded") {
        tie(value.quads, value.positions, value.normals,
            value.texturecoords) = make_cylinder_rounded_shape(js.value("steps",
                                                                   vec3i{64, 32,
                                                                       16}),
            js.value("size", vec2f{2.0f, 2.0f}),
            js.value("uvsize", vec3f{1, 1, 1}), js.value("radius", 0.15f));
    } else if (type == "sphere_geodesic") {
        tie(value.triangles, value.positions,
            value.normals) = make_geodesic_sphere_shape(js.value("tesselation", 4),
            js.value("size", 2.0f));
    } else if (type == "floor") {
        tie(value.quads, value.positions, value.normals,
            value.texturecoords) = make_floor_shape(js.value("steps", vec2i{1, 1}),
            js.value("size", vec2f{40, 40}), js.value("uvsize", vec2f{20, 20}));
    } else if (type == "matball") {
        tie(value.quads, value.positions, value.normals,
            value.texturecoords) = make_sphere_shape(js.value("steps",
                                                         vec2i{64, 32}),
            js.value("size", 2.0f), js.value("uvsize", vec2f{1, 1}));
    } else if (type == "hairball") {
        auto [base_quads, base_positions, base_normals,
            base_texturecoords] = make_sphere_cube_shape(32,
            js.value("size", 2.0f) * 0.8f, 1);
        tie(value.lines, value.positions, value.normals, value.texturecoords,
            value.radius) = make_hair_shape(js.value("steps", vec2i{4, 65536}),
            {}, base_quads, base_positions, base_normals, base_texturecoords,
            js.value("length", vec2f{0.2f, 0.2f}),
            js.value("radius", vec2f{0.001f, 0.001f}),
            js.value("noise", vec2f{0, 0}), js.value("clump", vec2f{0, 0}));
    } else if (type == "hairball_interior") {
        tie(value.quads, value.positions, value.normals,
            value.texturecoords) = make_sphere_cube_shape(32,
            js.value("size", 2.0f) * 0.8f, 1);
    } else if (type == "suzanne") {
        tie(value.quads, value.positions) = make_suzanne_shape(
            js.value("size", 2.0f));
    } else if (type == "cube_posonly") {
        tie(value.quads, value.positions) = make_cube_posonly_shape(
            js.value("steps", vec3i{1, 1, 1}), js.value("size", vec3f{2, 2, 2}),
            js.value("uvsize", vec3f{1, 1, 1}));
    } else {
        log_error("unknown shape type {}", type);
        return false;
    }
    if (!value.quads.empty() && js.value("shell_thickness", 0.0f) > 0) {
        tie(value.quads, value.positions, value.normals,
            value.texturecoords) = make_shell_shape(value.quads,
            value.positions, value.normals, value.texturecoords,
            js.value("shell_thickness", 0.0f));
    }
    if (!value.quads.empty() && js.value("as_triangles", false)) {
        value.triangles = convert_quads_to_triangles(value.quads);
        value.quads     = {};
    }
    if (js.value("flipyz", false)) {
        for (auto& p : value.positions) p = {p.x, p.z, p.y};
        for (auto& n : value.normals) n = {n.x, n.z, n.y};
    }
    return true;
}

// Serialize struct
bool serialize_json_object(
    json& js, yocto_shape& value, const yocto_scene& scene, bool save) {
    static const auto def = yocto_shape();
    if (!serialize_json_objbegin(js, save)) return false;
    if (!serialize_json_value(js, value.name, "name", def.name, save))
        return false;
    if (!serialize_json_value(js, value.filename, "filename", def.filename, save))
        return false;
    if (!serialize_json_objref(
            js, value.material, "material", scene.materials, save))
        return false;
    if (!serialize_json_value(js, value.subdivision_level, "subdivision_level",
            def.subdivision_level, save))
        return false;
    if (!serialize_json_value(
            js, value.catmull_clark, "catmull_clark", def.catmull_clark, save))
        return false;
    if (!serialize_json_value(js, value.compute_vertex_normals,
            "compute_vertex_normals", def.compute_vertex_normals, save))
        return false;
    if (value.filename == "" || !save) {
        if (!serialize_json_value(js, value.points, "points", def.points, save))
            return false;
        if (!serialize_json_value(js, value.lines, "lines", def.lines, save))
            return false;
        if (!serialize_json_value(
                js, value.triangles, "triangles", def.triangles, save))
            return false;
        if (!serialize_json_value(js, value.quads, "quads", def.quads, save))
            return false;
        if (!serialize_json_value(
                js, value.positions, "positions", def.positions, save))
            return false;
        if (!serialize_json_value(js, value.normals, "normals", def.normals, save))
            return false;
        if (!serialize_json_value(js, value.texturecoords, "texturecoords",
                def.texturecoords, save))
            return false;
        if (!serialize_json_value(js, value.colors, "colors", def.colors, save))
            return false;
        if (!serialize_json_value(js, value.radius, "radius", def.radius, save))
            return false;
        if (!serialize_json_value(js, value.tangentspaces, "tangent_spaces",
                def.tangentspaces, save))
            return false;
    }
    if (!serialize_json_procedural(js, value, "!!proc", scene, save))
        return false;
    return true;
}

// Procedural commands for materials
bool apply_json_procedural(
    const json& js, yocto_surface& value, const yocto_scene& scene) {
    if (!serialize_json_objbegin((json&)js, false)) return false;
    auto type = js.value("type", ""s);
    if (type == "") return true;
    value.quads_positions     = {};
    value.quads_normals       = {};
    value.quads_texturecoords = {};
    value.quads_materials     = {};
    value.positions           = {};
    value.normals             = {};
    value.texturecoords       = {};
    if (type == "quad") {
        tie(value.quads_positions, value.positions, value.normals,
            value.texturecoords) = make_quad_shape(js.value("steps", vec2i{1, 1}),
            js.value("size", vec2f{2, 2}), js.value("uvsize", vec2f{1, 1}));
    } else if (type == "quady") {
        tie(value.quads_positions, value.positions, value.normals,
            value.texturecoords) = make_quad_shape(js.value("steps", vec2i{1, 1}),
            js.value("size", vec2f{2, 2}), js.value("uvsize", vec2f{1, 1}));
    } else if (type == "quad_stack") {
        tie(value.quads_positions, value.positions, value.normals,
            value.texturecoords) = make_quad_stack_shape(js.value("steps",
                                                             vec3i{1, 1, 1}),
            js.value("size", vec3f{2, 2, 2}), js.value("uvsize", vec2f{1, 1}));
    } else if (type == "cube") {
        tie(value.quads_positions, value.positions, value.normals,
            value.texturecoords) = make_cube_shape(js.value("steps",
                                                       vec3i{1, 1, 1}),
            js.value("size", vec3f{2, 2, 2}), js.value("uvsize", vec3f{1, 1, 1}));
    } else if (type == "cube_rounded") {
        tie(value.quads_positions, value.positions, value.normals,
            value.texturecoords) = make_cube_rounded_shape(js.value("steps",
                                                               vec3i{32, 32, 32}),
            js.value("size", vec3f{2, 2, 2}),
            js.value("uvsize", vec3f{1, 1, 1}), js.value("radius", 0.3f));
    } else if (type == "sphere") {
        tie(value.quads_positions, value.positions, value.normals,
            value.texturecoords) = make_sphere_shape(js.value("steps",
                                                         vec2i{64, 32}),
            js.value("size", 2.0f), js.value("uvsize", vec2f{1, 1}));
    } else if (type == "sphere_cube") {
        tie(value.quads_positions, value.positions, value.normals,
            value.texturecoords) = make_sphere_cube_shape(js.value("steps", 32),
            js.value("size", 2.0f), js.value("uvsize", 1.0f));
    } else if (type == "sphere_flipcap") {
        tie(value.quads_positions, value.positions, value.normals,
            value.texturecoords) = make_sphere_flipcap_shape(js.value("steps",
                                                                 vec2i{64, 32}),
            js.value("size", 2.0f), js.value("uvsize", vec2f{1, 1}),
            js.value("zflip", vec2f{-0.75f, +0.75f}));
    } else if (type == "disk") {
        tie(value.quads_positions, value.positions, value.normals,
            value.texturecoords) = make_disk_shape(js.value("steps",
                                                       vec2i{32, 16}),
            js.value("size", 2.0f), js.value("uvsize", vec2f{1, 1}));
    } else if (type == "disk_quad") {
        tie(value.quads_positions, value.positions, value.normals,
            value.texturecoords) = make_disk_quad_shape(js.value("steps", 32),
            js.value("size", 2.0f), js.value("uvsize", 1.0f));
    } else if (type == "disk_bulged") {
        tie(value.quads_positions, value.positions, value.normals,
            value.texturecoords) = make_disk_bulged_shape(js.value("steps", 32),
            js.value("size", 2.0f), js.value("uvsize", 1.0f),
            js.value("height", 0.25f));
    } else if (type == "cylinder_side") {
        tie(value.quads_positions, value.positions, value.normals,
            value.texturecoords) = make_cylinder_side_shape(js.value("steps",
                                                                vec2i{64, 32}),
            js.value("size", vec2f{2.0f, 2.0f}), js.value("uvsize", vec2f{1, 1}));
    } else if (type == "cylinder") {
        tie(value.quads_positions, value.positions, value.normals,
            value.texturecoords) = make_cylinder_shape(js.value("steps",
                                                           vec3i{64, 32, 16}),
            js.value("size", vec2f{2.0f, 2.0f}),
            js.value("uvsize", vec3f{1, 1, 1}));
    } else if (type == "cylinder_rounded") {
        tie(value.quads_positions, value.positions, value.normals,
            value.texturecoords) = make_cylinder_rounded_shape(js.value("steps",
                                                                   vec3i{64, 32,
                                                                       16}),
            js.value("size", vec2f{2.0f, 2.0f}),
            js.value("uvsize", vec3f{1, 1, 1}), js.value("radius", 0.15f));
    } else if (type == "floor") {
        tie(value.quads_positions, value.positions, value.normals,
            value.texturecoords) = make_floor_shape(js.value("steps", vec2i{1, 1}),
            js.value("size", vec2f{40, 40}), js.value("uvsize", vec2f{20, 20}));
    } else if (type == "matball") {
        tie(value.quads_positions, value.positions, value.normals,
            value.texturecoords) = make_sphere_shape(js.value("steps",
                                                         vec2i{64, 32}),
            js.value("size", 2.0f), js.value("uvsize", vec2f{1, 1}));
    } else if (type == "hairball_interior") {
        tie(value.quads_positions, value.positions, value.normals,
            value.texturecoords) = make_sphere_cube_shape(32,
            js.value("size", 2.0f) * 0.8f, 1);
    } else if (type == "suzanne") {
        tie(value.quads_positions, value.positions) = make_suzanne_shape(
            js.value("size", 2.0f));
    } else if (type == "cube_facevarying") {
        tie(value.quads_positions, value.quads_normals,
            value.quads_texturecoords, value.positions, value.normals,
            value.texturecoords) = make_cube_facevarying_shape(js.value("steps",
                                                                   vec3i{1, 1, 1}),
            js.value("size", vec3f{2, 2, 2}), js.value("uvsize", vec3f{1, 1, 1}));
    } else if (type == "cube_multiplematerials") {
        tie(value.quads_positions, value.quads_normals, value.quads_texturecoords,
            value.quads_materials, value.positions, value.normals,
            value.texturecoords) = make_cube_multiplematerials_shape(js.value("steps",
                                                                         vec3i{1, 1,
                                                                             1}),
            js.value("size", vec3f{2, 2, 2}), js.value("uvsize", vec3f{1, 1, 1}));
    } else if (type == "cube_posonly") {
        tie(value.quads_positions, value.positions) = make_cube_posonly_shape(
            js.value("steps", vec3i{1, 1, 1}), js.value("size", vec3f{2, 2, 2}),
            js.value("uvsize", vec3f{1, 1, 1}));
    } else {
        log_error("unknown shape type {}", type);
        return false;
    }
    if (empty(value.quads_normals) && !empty(value.normals))
        value.quads_normals = value.quads_positions;
    if (empty(value.quads_texturecoords) && !empty(value.texturecoords))
        value.quads_texturecoords = value.quads_positions;
    if (js.value("flipyz", false)) {
        for (auto& p : value.positions) p = {p.x, p.z, p.y};
        for (auto& n : value.normals) n = {n.x, n.z, n.y};
    }
    return true;
}

// Serialize struct
bool serialize_json_object(
    json& js, yocto_surface& value, const yocto_scene& scene, bool save) {
    static const auto def = yocto_surface();
    if (!serialize_json_objbegin(js, save)) return false;
    if (!serialize_json_value(js, value.name, "name", def.name, save))
        return false;
    if (!serialize_json_value(js, value.filename, "filename", def.filename, save))
        return false;
    if (!serialize_json_objref(
            js, value.materials, "materials", scene.materials, save))
        return false;
    if (!serialize_json_value(js, value.subdivision_level, "subdivision_level",
            def.subdivision_level, save))
        return false;
    if (!serialize_json_value(
            js, value.catmull_clark, "catmull_clark", def.catmull_clark, save))
        return false;
    if (!serialize_json_value(js, value.compute_vertex_normals,
            "compute_vertex_normals", def.compute_vertex_normals, save))
        return false;
    if (value.filename == "" || !save) {
        if (!serialize_json_value(js, value.quads_positions, "quads_positions",
                def.quads_positions, save))
            return false;
        if (!serialize_json_value(js, value.quads_normals, "quads_normals",
                def.quads_normals, save))
            return false;
        if (!serialize_json_value(js, value.quads_texturecoords,
                "quads_texturecoords", def.quads_texturecoords, save))
            return false;
        if (!serialize_json_value(js, value.quads_materials, "quads_materials",
                def.quads_materials, save))
            return false;
        if (!serialize_json_value(
                js, value.positions, "positions", def.positions, save))
            return false;
        if (!serialize_json_value(js, value.normals, "normals", def.normals, save))
            return false;
        if (!serialize_json_value(js, value.texturecoords, "texturecoords",
                def.texturecoords, save))
            return false;
    }
    if (!serialize_json_procedural(js, value, "!!proc", scene, save))
        return false;
    return true;
}

// Procedural commands for instances
bool apply_json_procedural(
    const json& js, yocto_instance& value, const yocto_scene& scene) {
    if (!serialize_json_objbegin((json&)js, false)) return false;
    if (js.count("from")) {
        auto from   = js.value("from", zero3f);
        auto to     = js.value("to", zero3f);
        auto up     = js.value("up", vec3f{0, 1, 0});
        value.frame = make_lookat_frame(from, to, up, true);
    }
    if (js.count("translation") || js.count("rotation") || js.count("scale")) {
        auto translation = js.value("translation", zero3f);
        auto rotation    = js.value("rotation", zero4f);
        auto scaling     = js.value("scale", vec3f{1, 1, 1});
        value.frame      = make_translation_frame(translation) *
                      make_scaling_frame(scaling) *
                      make_rotation_frame(xyz(rotation), rotation.w);
    }
    return true;
}

// Serialize struct
bool serialize_json_object(
    json& js, yocto_instance& value, const yocto_scene& scene, bool save) {
    static const auto def = yocto_instance();
    if (!serialize_json_objbegin(js, save)) return false;
    if (!serialize_json_value(js, value.name, "name", def.name, save))
        return false;
    if (!serialize_json_value(js, value.frame, "frame", def.frame, save))
        return false;
    if (!serialize_json_objref(js, value.shape, "shape", scene.shapes, save))
        return false;
    if (!serialize_json_objref(js, value.surface, "surface", scene.surfaces, save))
        return false;
    if (!serialize_json_procedural(js, value, "!!proc", scene, save))
        return false;
    return true;
}

// Procedural commands for materials
bool apply_json_procedural(
    const json& js, yocto_environment& value, const yocto_scene& scene) {
    if (!serialize_json_objbegin((json&)js, false)) return false;
    if (js.count("rotation")) {
        auto rotation = js.value("rotation", zero4f);
        value.frame   = make_rotation_frame(xyz(rotation), rotation.w);
    }
    return true;
}

// Serialize struct
bool serialize_json_object(
    json& js, yocto_environment& value, const yocto_scene& scene, bool save) {
    static const auto def = yocto_environment();
    if (!serialize_json_objbegin(js, save)) return false;
    if (!serialize_json_value(js, value.name, "name", def.name, save))
        return false;
    if (!serialize_json_value(js, value.frame, "frame", def.frame, save))
        return false;
    if (!serialize_json_value(js, value.emission, "emission", def.emission, save))
        return false;
    if (!serialize_json_objref(js, value.emission_texture, "emission_texture",
            scene.textures, save))
        return false;
    if (!serialize_json_procedural(js, value, "!!proc", scene, save))
        return false;
    return true;
}

// Procedural commands for nodes
bool apply_json_procedural(
    const json& js, yocto_scene_node& value, const yocto_scene& scene) {
    if (!serialize_json_objbegin((json&)js, false)) return false;
    if (js.count("from")) {
        auto from   = js.value("from", zero3f);
        auto to     = js.value("to", zero3f);
        auto up     = js.value("up", vec3f{0, 1, 0});
        value.local = make_lookat_frame(from, to, up, true);
    }
    return true;
}

// Serialize struct
bool serialize_json_object(
    json& js, yocto_scene_node& value, const yocto_scene& scene, bool save) {
    static const auto def = yocto_scene_node();
    if (!serialize_json_objbegin(js, save)) return false;
    if (!serialize_json_value(js, value.name, "name", def.name, save))
        return false;
    if (!serialize_json_value(js, value.local, "local", def.local, save))
        return false;
    if (!serialize_json_value(
            js, value.translation, "translation", def.translation, save))
        return false;
    if (!serialize_json_value(js, value.rotation, "rotation", def.rotation, save))
        return false;
    if (!serialize_json_value(js, value.scale, "scale", def.scale, save))
        return false;
    if (!serialize_json_value(js, value.weights, "weights", def.weights, save))
        return false;
    if (!serialize_json_objref(js, value.parent, "parent", scene.nodes, save))
        return false;
    if (!serialize_json_objref(js, value.camera, "camera", scene.cameras, save))
        return false;
    if (!serialize_json_objref(
            js, value.instance, "instance", scene.instances, save))
        return false;
    if (!serialize_json_objref(
            js, value.environment, "environment", scene.environments, save))
        return false;
    if (!serialize_json_procedural(js, value, "!!proc", scene, save))
        return false;
    return true;
}

// Serialize enum
bool serialize_json_value(json& js, yocto_interpolation_type& value, bool save) {
    if (save) {
        static auto names = unordered_map<int, string>{
            {(int)yocto_interpolation_type::linear, "linear"},
            {(int)yocto_interpolation_type::step, "step"},
            {(int)yocto_interpolation_type::bezier, "bezier"},
        };
        auto vals = names.at((int)value);
        return serialize_json_value(js, vals, save);
    } else {
        static auto names = unordered_map<string, int>{
            {"linear", (int)yocto_interpolation_type::linear},
            {"step", (int)yocto_interpolation_type::step},
            {"bezier", (int)yocto_interpolation_type::bezier},
        };
        auto vals = ""s;
        if (!serialize_json_value(js, vals, save)) return false;
        try {
            value = (yocto_interpolation_type)names.at(vals);
        } catch (...) {
            return false;
        }
        return true;
    }
}

// Serialize struct
bool serialize_json_object(
    json& js, yocto_animation& value, const yocto_scene& scene, bool save) {
    static const auto def = yocto_animation();
    if (!serialize_json_objbegin(js, save)) return false;
    if (!serialize_json_value(js, value.name, "name", def.name, save))
        return false;
    if (!serialize_json_value(js, value.filename, "filename", def.filename, save))
        return false;
    if (!serialize_json_value(js, value.animation_group, "animation_group",
            def.animation_group, save))
        return false;
    if (!serialize_json_value(
            js, value.interpolation_type, "type", def.interpolation_type, save))
        return false;
    if (value.filename == "" || !save) {
        if (!serialize_json_value(js, value.keyframes_times, "keyframes_times",
                def.keyframes_times, save))
            return false;
        if (!serialize_json_value(js, value.translation_keyframes,
                "translation_keyframes", def.translation_keyframes, save))
            return false;
        if (!serialize_json_value(js, value.rotation_keyframes,
                "rotation_keyframes", def.rotation_keyframes, save))
            return false;
        if (!serialize_json_value(js, value.scale_keyframes, "scale_keyframes",
                def.scale_keyframes, save))
            return false;
    }
    if (!serialize_json_objref(
            js, value.node_targets, "node_targets", scene.nodes, save))
        return false;
    if (!serialize_json_procedural(js, value, "!!proc", scene, save))
        return false;
    return true;
}

// Procedural commands for animations
bool apply_json_procedural(
    const json& js, yocto_animation& value, const yocto_scene& scene) {
    if (!serialize_json_objbegin((json&)js, false)) return false;
    if (js.count("make_rotation_axisangle")) {
        for (auto& j : js.at("make_rotation_axisangle")) {
            value.rotation_keyframes.push_back(
                make_rotation_quat(j.get<vec4f>()));
        }
    }
    return true;
}

// Procedural commands for scenes
bool apply_json_procedural(
    const json& js, yocto_scene& value, const yocto_scene& scene) {
    if (!serialize_json_objbegin((json&)js, false)) return false;
    if (js.count("random_instances")) {
        auto& jjs          = js.at("random_instances");
        auto  num          = jjs.value("num", 100);
        auto  seed         = jjs.value("seed", 13);
        auto  shape_offset = (int)scene.shapes.size();
        auto  num_shapes   = 0;
        auto  base         = yocto_shape();
        serialize_json_object((json&)jjs.at("base"), base, scene, false);
        for (auto& j : jjs.at("shapes")) {
            value.shapes.push_back({});
            serialize_json_object((json&)j, value.shapes.back(), scene, false);
            num_shapes++;
        }

        auto [pos, norm, texcoord] = sample_triangles_points(base.triangles,
            base.positions, base.normals, base.texturecoords, num, seed);

        auto rng = make_rng(seed, 17);
        for (auto i = 0; i < num; i++) {
            auto shape = get_random_int(rng, num_shapes - 1);
            value.instances.push_back({});
            value.instances.back().name  = "random_" + std::to_string(i);
            value.instances.back().frame = make_translation_frame(pos[i]);
            value.instances.back().shape = shape + shape_offset;
        }
    }
    return true;
}

// Serialize struct
bool serialize_json_object(
    json& js, yocto_scene& value, const yocto_scene& scene, bool save) {
    static const auto def = yocto_scene();
    if (!serialize_json_objbegin(js, save)) return false;
    if (!serialize_json_value(js, value.name, "name", def.name, save))
        return false;
    if (!serialize_json_objarray(js, value.cameras, "cameras", scene, save))
        return false;
    if (!serialize_json_objarray(js, value.textures, "textures", scene, save))
        return false;
    if (!serialize_json_objarray(
            js, value.voltextures, "voltextures", scene, save))
        return false;
    if (!serialize_json_objarray(js, value.materials, "materials", scene, save))
        return false;
    if (!serialize_json_objarray(js, value.shapes, "shapes", scene, save))
        return false;
    if (!serialize_json_objarray(js, value.surfaces, "surfaces", scene, save))
        return false;
    if (!serialize_json_objarray(js, value.instances, "instances", scene, save))
        return false;
    if (!serialize_json_objarray(
            js, value.environments, "environments", scene, save))
        return false;
    if (!serialize_json_objarray(js, value.nodes, "nodes", scene, save))
        return false;
    if (!serialize_json_objarray(js, value.animations, "animations", scene, save))
        return false;
    if (!serialize_json_procedural(js, value, "!!proc", scene, save))
        return false;
    return true;
}

bool serialize_json_object(json& js, yocto_scene& value, bool save) {
    return serialize_json_object(js, value, value, save);
}

// Load json meshes
bool load_json_meshes(yocto_scene& scene, const string& dirname,
    const load_scene_options& options) {
    if (options.skip_meshes) return true;

    // load shapes
    atomic<bool> exit_error(false);
    parallel_foreach(scene.shapes,
        [&exit_error, &options, &dirname](yocto_shape& shape) {
            if (exit_error) return;
            if (shape.filename == "" || !empty(shape.positions)) return;
            auto filename = normalize_path(dirname + shape.filename);
            if (!load_mesh(filename, shape.points, shape.lines, shape.triangles,
                    shape.quads, shape.positions, shape.normals,
                    shape.texturecoords, shape.colors, shape.radius, false)) {
                if (options.exit_on_error) {
                    exit_error = true;
                    return;
                }
            }
        },
        options.cancel_flag, options.run_serially);
    if (exit_error) return false;

    // load surfaces
    parallel_foreach(scene.surfaces,
        [&exit_error, &options, &dirname](yocto_surface& surface) {
            if (exit_error) return;
            if (surface.filename == "" || !empty(surface.positions)) return;
            auto filename = normalize_path(dirname + surface.filename);
            if (!load_facevarying_mesh(filename, surface.quads_positions,
                    surface.quads_normals, surface.quads_texturecoords,
                    surface.positions, surface.normals, surface.texturecoords,
                    surface.quads_materials)) {
                if (options.exit_on_error) {
                    exit_error = true;
                    return;
                }
            }
        },
        options.cancel_flag, options.run_serially);
    if (exit_error) return false;

    // done
    return true;
}

// Load a scene in the builtin JSON format.
bool load_json_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options) {
    auto scope = log_trace_scoped("loading scene {}", filename);
    // initialize
    scene = {};

    // load jsonz
    auto js = json();
    if (!load_json(filename, js)) return false;

    // deserialize json
    try {
        if (!serialize_json_object(js, scene, false)) {
            log_io_error("could not deserialize json {}", filename);
            return false;
        }
    } catch (...) {
        log_io_error("could not deserialize json {}", filename);
        return false;
    }

    // load meshes and textures
    auto dirname = get_dirname(filename);
    if (!load_json_meshes(scene, dirname, options)) return false;
    if (!load_scene_textures(scene, dirname, options)) return false;

    // fix scene
    if (scene.name == "") scene.name = get_filename(filename);
    add_missing_cameras(scene);
    add_missing_materials(scene);
    add_missing_names(scene);
    update_transforms(scene);

    // done
    return true;
}

// Save json meshes
bool save_json_meshes(const yocto_scene& scene, const string& dirname,
    const save_scene_options& options) {
    if (options.skip_meshes) return true;

    // save shapes
    atomic<bool> exit_error(false);
    parallel_foreach(scene.shapes,
        [&exit_error, &options, &dirname](const yocto_shape& shape) {
            if (exit_error) return;
            if (shape.filename == "") return;
            auto filename = normalize_path(dirname + shape.filename);
            if (!save_mesh(filename, shape.points, shape.lines, shape.triangles,
                    shape.quads, shape.positions, shape.normals,
                    shape.texturecoords, shape.colors, shape.radius)) {
                if (options.exit_on_error) {
                    exit_error = true;
                    return;
                }
            }
        },
        options.cancel_flag, options.run_serially);
    if (exit_error) return false;

    // save surfaces
    parallel_foreach(scene.surfaces,
        [&exit_error, &options, &dirname](const yocto_surface& surface) {
            if (exit_error) return;
            if (surface.filename == "") return;
            auto filename = normalize_path(dirname + surface.filename);
            if (!save_facevarying_mesh(filename, surface.quads_positions,
                    surface.quads_normals, surface.quads_texturecoords,
                    surface.positions, surface.normals, surface.texturecoords,
                    surface.quads_materials)) {
                if (options.exit_on_error) {
                    exit_error = true;
                    return;
                }
            }
        },
        options.cancel_flag, options.run_serially);
    if (exit_error) return false;

    // done
    return true;
}

// Save a scene in the builtin JSON format.
bool save_json_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options) {
    auto scope = log_trace_scoped("saving scene {}", filename);
    // save json
    auto js = json();
    try {
        if (!serialize_json_object(js, (yocto_scene&)scene, true)) {
            log_io_error("could not serialize json {}", filename);
            return false;
        }
    } catch (...) {
        log_io_error("could not serialize json {}", filename);
        return false;
    }
    if (!save_json(filename, js)) return false;

    // save meshes and textures
    auto dirname = get_dirname(filename);
    if (!save_json_meshes(scene, dirname, options)) return false;
    if (!save_scene_textures(scene, dirname, options)) return false;

    // done
    return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// OBJ CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

inline bool operator==(const obj_vertex& a, const obj_vertex& b) {
    return a.position == b.position && a.texturecoord == b.texturecoord &&
           a.normal == b.normal;
}

struct obj_vertex_hash {
    size_t operator()(const obj_vertex& v) const {
        auto vh = std::hash<int>();
        auto h  = (size_t)0;
        for (auto i = 0; i < 3; i++)
            h ^= vh((&v.position)[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

inline bool parse_value(string_view& view, obj_vertex& value) {
    value = obj_vertex{0, 0, 0};
    if (!parse_value(view, value.position)) return false;
    if (view.front() == '/') {
        view.remove_prefix(1);
        if (view.front() == '/') {
            view.remove_prefix(1);
            if (!parse_value(view, value.normal)) return false;
        } else {
            if (!parse_value(view, value.texturecoord)) return false;
            if (view.front() == '/') {
                view.remove_prefix(1);
                if (!parse_value(view, value.normal)) return false;
            }
        }
    }
    return true;
}

// Input for OBJ textures
inline bool parse_value(string_view& view, obj_texture_info& info) {
    // initialize
    info = obj_texture_info();

    // get tokens
    auto tokens = vector<string>();
    while (true) {
        auto token = ""s;
        parse_value(view, token);
        if (token == "") break;
        tokens.push_back(token);
    }
    if (empty(tokens)) return false;

    // texture name
    info.path = normalize_path(tokens.back());

    // texture options
    auto last = string();
    for (auto i = 0; i < tokens.size() - 1; i++) {
        if (tokens[i] == "-bm") info.scale = atof(tokens[i + 1].c_str());
        if (tokens[i] == "-clamp") info.clamp = true;
    }

    return true;
}

// Load obj materials
bool load_mtl(const string& filename, const obj_callbacks& cb,
    const load_obj_options& options) {
    // open file
    auto fs = open(filename, "rt");
    if (!fs) return false;

    // currently parsed material
    auto material = obj_material();
    auto first    = true;

    // read the file line by line
    auto line = ""s;
    while (read_line(fs, line)) {
        // line
        if (line.find('#') != line.npos) line = line.substr(0, line.find('#'));
        auto view = string_view{line};

        // get command
        auto cmd = ""s;
        parse_value(view, cmd);
        if (cmd == "") continue;

        // possible token values
        if (cmd == "newmtl") {
            if (!first && cb.material) cb.material(material);
            first    = false;
            material = obj_material();
            parse_value(view, material.name);
        } else if (cmd == "illum") {
            parse_value(view, material.illum);
        } else if (cmd == "Ke") {
            parse_value(view, material.ke);
        } else if (cmd == "Kd") {
            parse_value(view, material.kd);
        } else if (cmd == "Ks") {
            parse_value(view, material.ks);
        } else if (cmd == "Kt") {
            parse_value(view, material.kt);
        } else if (cmd == "Tf") {
            material.kt = {-1, -1, -1};
            parse_value(view, material.kt);
            if (material.kt.y < 0)
                material.kt = {material.kt.x, material.kt.x, material.kt.x};
            if (options.flip_tr) material.kt = vec3f{1, 1, 1} - material.kt;
        } else if (cmd == "Tr") {
            auto tr = vec3f{-1, -1, -1};
            parse_value(view, tr);
            if (tr.y < 0) tr = {tr.x, tr.x, tr.x};
            material.op = (tr.x + tr.y + tr.z) / 3;
            if (options.flip_tr) material.op = 1 - material.op;
        } else if (cmd == "Ns") {
            parse_value(view, material.ns);
            material.rs = pow(2 / (material.ns + 2), 1 / 4.0f);
            if (material.rs < 0.01f) material.rs = 0;
            if (material.rs > 0.99f) material.rs = 1;
        } else if (cmd == "d") {
            parse_value(view, material.op);
        } else if (cmd == "Pr" || cmd == "rs") {
            parse_value(view, material.rs);
        } else if (cmd == "map_Ke") {
            parse_value(view, material.ke_txt);
        } else if (cmd == "map_Kd") {
            parse_value(view, material.kd_txt);
        } else if (cmd == "map_Ks") {
            parse_value(view, material.ks_txt);
        } else if (cmd == "map_Tr") {
            parse_value(view, material.kt_txt);
        } else if (cmd == "map_d" || cmd == "map_Tr") {
            parse_value(view, material.op_txt);
        } else if (cmd == "map_Pr" || cmd == "map_rs") {
            parse_value(view, material.rs_txt);
        } else if (cmd == "map_occ" || cmd == "occ") {
            parse_value(view, material.occ_txt);
        } else if (cmd == "map_bump" || cmd == "bump") {
            parse_value(view, material.bump_txt);
        } else if (cmd == "map_disp" || cmd == "disp") {
            parse_value(view, material.disp_txt);
        } else if (cmd == "map_norm" || cmd == "norm") {
            parse_value(view, material.norm_txt);
        } else if (cmd == "Ve") {
            parse_value(view, material.ve);
        } else if (cmd == "Va") {
            parse_value(view, material.va);
        } else if (cmd == "Vd") {
            parse_value(view, material.vd);
        } else if (cmd == "Vg") {
            parse_value(view, material.vg);
        } else if (cmd == "map_Vd") {
            parse_value(view, material.vd_txt);
        }
    }

    // issue current material
    if (!first && cb.material) cb.material(material);

    // done
    return true;
}

// Load obj extensions
bool load_objx(const string& filename, const obj_callbacks& cb,
    const load_obj_options& options) {
    // open file
    auto fs = open(filename, "rt");
    if (!fs) return false;

    // read the file line by line
    auto line = ""s;
    while (read_line(fs, line)) {
        // line
        if (line.find('#') != line.npos) line = line.substr(0, line.find('#'));
        auto view = string_view{line.c_str()};

        // get command
        auto cmd = ""s;
        parse_value(view, cmd);
        if (cmd == "") continue;

        // possible token values
        if (cmd == "c") {
            auto camera = obj_camera();
            parse_value(view, camera.name);
            parse_value(view, camera.ortho);
            parse_value(view, camera.width);
            parse_value(view, camera.height);
            parse_value(view, camera.focal);
            parse_value(view, camera.focus);
            parse_value(view, camera.aperture);
            parse_value(view, camera.frame);
            if (cb.camera) cb.camera(camera);
        } else if (cmd == "e") {
            auto environment = obj_environment();
            parse_value(view, environment.name);
            parse_value(view, environment.ke);
            parse_value(view, environment.ke_txt.path);
            parse_value(view, environment.frame);
            if (environment.ke_txt.path == "\"\"") environment.ke_txt.path = "";
            if (cb.environmnet) cb.environmnet(environment);
        } else {
            // unused
        }
    }

    // done
    return true;
}

// Load obj scene
bool load_obj(const string& filename, const obj_callbacks& cb,
    const load_obj_options& options) {
    // open file
    auto fs = open(filename, "rt");
    if (!fs) return false;

    // track vertex size
    auto vert_size = obj_vertex();
    auto verts     = vector<obj_vertex>();  // buffer to avoid reallocation

    // read the file line by line
    auto line = ""s;
    while (read_line(fs, line)) {
        // line
        if (line.find('#') != line.npos) line = line.substr(0, line.find('#'));
        auto view = string_view{line.c_str()};

        // get command
        auto cmd = ""s;
        parse_value(view, cmd);
        if (cmd == "") continue;

        // possible token values
        if (cmd == "v") {
            auto vert = zero3f;
            parse_value(view, vert);
            if (cb.vert) cb.vert(vert);
            vert_size.position += 1;
        } else if (cmd == "vn") {
            auto vert = zero3f;
            parse_value(view, vert);
            if (cb.norm) cb.norm(vert);
            vert_size.normal += 1;
        } else if (cmd == "vt") {
            auto vert = zero2f;
            parse_value(view, vert);
            if (options.flip_texcoord) vert.y = 1 - vert.y;
            if (cb.texcoord) cb.texcoord(vert);
            vert_size.texturecoord += 1;
        } else if (cmd == "f" || cmd == "l" || cmd == "p") {
            verts.clear();
            while (true) {
                auto vert = obj_vertex{};
                parse_value(view, vert);
                if (!vert.position) break;
                if (vert.position < 0)
                    vert.position = vert_size.position + vert.position + 1;
                if (vert.texturecoord < 0)
                    vert.texturecoord = vert_size.texturecoord +
                                        vert.texturecoord + 1;
                if (vert.normal < 0)
                    vert.normal = vert_size.normal + vert.normal + 1;
                verts.push_back(vert);
            }
            if (cmd == "f" && cb.face) cb.face(verts);
            if (cmd == "l" && cb.line) cb.line(verts);
            if (cmd == "p" && cb.point) cb.point(verts);
        } else if (cmd == "o") {
            auto name = ""s;
            parse_value(view, name);
            if (cb.object) cb.object(name);
        } else if (cmd == "usemtl") {
            auto name = ""s;
            parse_value(view, name);
            if (cb.usemtl) cb.usemtl(name);
        } else if (cmd == "g") {
            auto name = ""s;
            parse_value(view, name);
            if (cb.group) cb.group(name);
        } else if (cmd == "s") {
            auto name = ""s;
            parse_value(view, name);
            if (cb.smoothing) cb.smoothing(name);
        } else if (cmd == "mtllib") {
            if (options.geometry_only) continue;
            auto mtlname = ""s;
            parse_value(view, mtlname);
            if (cb.mtllib) cb.mtllib(mtlname);
            auto mtlpath = get_dirname(filename) + mtlname;
            if (!load_mtl(mtlpath, cb, options)) {
                if (options.exit_on_error) return false;
            }
        } else {
            // unused
        }
    }

    // parse extensions if presents
    if (!options.geometry_only) {
        auto extname    = replace_extension(filename, "objx");
        auto ext_exists = exists_file(extname);
        if (ext_exists) {
            if (!load_objx(extname, cb, options)) return false;
        }
    }

    // done
    return true;
}

// Loads an OBJ
bool load_obj_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options) {
    auto scope = log_trace_scoped("loading scene {}", filename);
    scene      = {};

    // splitting policy
    auto split_material  = options.obj_split_shapes;
    auto split_group     = options.obj_split_shapes;
    auto split_smoothing = options.obj_split_shapes;

    // current parsing values
    auto matname   = string();
    auto oname     = string();
    auto gname     = string();
    auto smoothing = true;

    // vertices
    auto opos      = deque<vec3f>();
    auto onorm     = deque<vec3f>();
    auto otexcoord = deque<vec2f>();

    // object maps
    auto tmap = unordered_map<string, int>();
    auto vmap = unordered_map<string, int>();
    auto mmap = unordered_map<string, int>();

    // vertex maps
    auto vertex_map   = unordered_map<obj_vertex, int, obj_vertex_hash>();
    auto pos_map      = unordered_map<int, int>();
    auto norm_map     = unordered_map<int, int>();
    auto texcoord_map = unordered_map<int, int>();

    // add object if needed
    auto is_instance_empty = [](const yocto_scene&     scene,
                                 const yocto_instance& instance) {
        if (instance.shape >= 0)
            return empty(scene.shapes[instance.shape].positions);
        if (instance.surface >= 0)
            return empty(scene.surfaces[instance.surface].positions);
        return true;
    };
    auto add_instance = [&](yocto_scene& scene, const string& objname,
                            const string& groupname) {
        auto instance = yocto_instance();
        instance.name = !empty(objname) ? objname : groupname;
        scene.instances.push_back(instance);
        vertex_map.clear();
        pos_map.clear();
        norm_map.clear();
        texcoord_map.clear();
    };
    // Get current material
    auto get_material_id = [&mmap](const string& matname) {
        if (matname != "") {
            auto it = mmap.find(matname);
            if (it == mmap.end()) {
                log_error("missing material {}", matname);
                return -1;
            }
            return it->second;
        }
        return -1;
    };
    // Parse texture options and name
    auto add_texture = [&scene, &tmap](
                           const obj_texture_info& info, bool force_linear) {
        if (info.path == "") return -1;
        if (tmap.find(info.path) != tmap.end()) {
            return tmap.at(info.path);
        }

        // create texture
        auto texture          = yocto_texture{};
        texture.name          = info.path;
        texture.filename      = info.path;
        texture.clamp_to_edge = info.clamp;
        texture.height_scale  = info.scale;
        texture.ldr_as_linear = force_linear || is_hdr_filename(info.path);
        scene.textures.push_back(texture);
        auto index      = (int)scene.textures.size() - 1;
        tmap[info.path] = index;

        return index;
    };
    // Parse texture options and name
    auto add_voltexture = [&scene, &vmap](
                              const obj_texture_info& info, bool srgb) {
        if (info.path == "") return -1;
        if (vmap.find(info.path) != vmap.end()) {
            return vmap.at(info.path);
        }

        // create texture
        auto texture     = yocto_voltexture{};
        texture.name     = info.path;
        texture.filename = info.path;
        scene.voltextures.push_back(texture);
        auto index      = (int)scene.voltextures.size() - 1;
        vmap[info.path] = index;

        return index;
    };
    // Add  vertices to the current shape
    auto add_verts = [&](const vector<obj_vertex>& verts, yocto_shape& shape) {
        for (auto& vert : verts) {
            auto it = vertex_map.find(vert);
            if (it != vertex_map.end()) continue;
            auto& shape  = scene.shapes.back();
            auto  nverts = (int)shape.positions.size();
            vertex_map.insert(it, {vert, nverts});
            if (vert.position)
                shape.positions.push_back(opos.at(vert.position - 1));
            if (vert.texturecoord)
                shape.texturecoords.push_back(
                    otexcoord.at(vert.texturecoord - 1));
            if (vert.normal) shape.normals.push_back(onorm.at(vert.normal - 1));
        }
    };
    // add vertex
    auto add_fvverts = [&](const vector<obj_vertex>& verts,
                           yocto_surface&            surface) {
        for (auto& vert : verts) {
            if (!vert.position) continue;
            auto pos_it = pos_map.find(vert.position);
            if (pos_it != pos_map.end()) continue;
            auto nverts = (int)surface.positions.size();
            pos_map.insert(pos_it, {vert.position, nverts});
            surface.positions.push_back(opos.at(vert.position - 1));
        }
        for (auto& vert : verts) {
            if (!vert.texturecoord) continue;
            auto texcoord_it = texcoord_map.find(vert.texturecoord);
            if (texcoord_it != texcoord_map.end()) continue;
            auto nverts = (int)surface.texturecoords.size();
            texcoord_map.insert(texcoord_it, {vert.texturecoord, nverts});
            surface.texturecoords.push_back(otexcoord.at(vert.texturecoord - 1));
        }
        for (auto& vert : verts) {
            if (!vert.normal) continue;
            auto norm_it = norm_map.find(vert.normal);
            if (norm_it != norm_map.end()) continue;
            auto nverts = (int)surface.normals.size();
            norm_map.insert(norm_it, {vert.normal, nverts});
            surface.normals.push_back(onorm.at(vert.normal - 1));
        }
    };

    // callbacks
    auto cb     = obj_callbacks();
    cb.vert     = [&](vec3f v) { opos.push_back(v); };
    cb.norm     = [&](vec3f v) { onorm.push_back(v); };
    cb.texcoord = [&](vec2f v) { otexcoord.push_back(v); };
    cb.face     = [&](const vector<obj_vertex>& verts) {
        if (empty(scene.instances)) add_instance(scene, oname, gname);
        if (scene.instances.back().shape < 0 &&
            scene.instances.back().surface < 0) {
            if (options.obj_preserve_face_varying ||
                scene.instances.back().name.find("[yocto::facevarying]") !=
                    string::npos) {
                scene.surfaces.push_back({});
                scene.surfaces.back().name = scene.instances.back().name;
                scene.surfaces.back().materials.push_back(
                    get_material_id(matname));
                scene.instances.back().surface = (int)scene.surfaces.size() - 1;
            } else {
                scene.shapes.push_back({});
                scene.shapes.back().name     = scene.instances.back().name;
                scene.shapes.back().material = get_material_id(matname);
                scene.instances.back().shape = (int)scene.shapes.size() - 1;
            }
        }
        if (scene.instances.back().surface >= 0) {
            auto& surface = scene.surfaces.back();
            add_fvverts(verts, surface);
            if (verts.size() == 4) {
                if (verts[0].position) {
                    surface.quads_positions.push_back(
                        {pos_map.at(verts[0].position),
                            pos_map.at(verts[1].position),
                            pos_map.at(verts[2].position),
                            pos_map.at(verts[3].position)});
                }
                if (verts[0].texturecoord) {
                    surface.quads_texturecoords.push_back(
                        {texcoord_map.at(verts[0].texturecoord),
                            texcoord_map.at(verts[1].texturecoord),
                            texcoord_map.at(verts[2].texturecoord),
                            texcoord_map.at(verts[3].texturecoord)});
                }
                if (verts[0].normal) {
                    surface.quads_normals.push_back({norm_map.at(verts[0].normal),
                        norm_map.at(verts[1].normal), norm_map.at(verts[2].normal),
                        norm_map.at(verts[3].normal)});
                }
            } else {
                if (verts[0].position) {
                    for (auto i = 2; i < verts.size(); i++)
                        surface.quads_positions.push_back(
                            {pos_map.at(verts[0].position),
                                pos_map.at(verts[i - 1].position),
                                pos_map.at(verts[i].position),
                                pos_map.at(verts[i].position)});
                }
                if (verts[0].texturecoord) {
                    for (auto i = 2; i < verts.size(); i++)
                        surface.quads_texturecoords.push_back(
                            {texcoord_map.at(verts[0].texturecoord),
                                texcoord_map.at(verts[i - 1].texturecoord),
                                texcoord_map.at(verts[i].texturecoord),
                                texcoord_map.at(verts[i].texturecoord)});
                }
                if (verts[0].normal) {
                    for (auto i = 2; i < verts.size(); i++)
                        surface.quads_normals.push_back(
                            {norm_map.at(verts[0].normal),
                                norm_map.at(verts[i - 1].normal),
                                norm_map.at(verts[i].normal),
                                norm_map.at(verts[i].normal)});
                }
            }
        } else {
            auto& shape = scene.shapes.back();
            add_verts(verts, shape);
            if (verts.size() == 4) {
                shape.quads.push_back(
                    {vertex_map.at(verts[0]), vertex_map.at(verts[1]),
                        vertex_map.at(verts[2]), vertex_map.at(verts[3])});
            } else {
                for (auto i = 2; i < verts.size(); i++)
                    shape.triangles.push_back({vertex_map.at(verts[0]),
                        vertex_map.at(verts[i - 1]), vertex_map.at(verts[i])});
            }
        }
    };
    cb.line = [&](const vector<obj_vertex>& verts) {
        if (empty(scene.instances)) add_instance(scene, oname, gname);
        if (scene.instances.back().surface >= 0)
            add_instance(scene, oname, gname);
        if (scene.instances.back().shape < 0) {
            scene.shapes.push_back({});
            scene.shapes.back().name     = scene.instances.back().name;
            scene.shapes.back().material = get_material_id(matname);
            scene.instances.back().shape = (int)scene.shapes.size() - 1;
        }
        auto& shape = scene.shapes.back();
        add_verts(verts, shape);
        for (auto i = 1; i < verts.size(); i++)
            shape.lines.push_back(
                {vertex_map.at(verts[i - 1]), vertex_map.at(verts[i])});
    };
    cb.point = [&](const vector<obj_vertex>& verts) {
        if (empty(scene.instances)) add_instance(scene, oname, gname);
        if (scene.instances.back().surface >= 0)
            add_instance(scene, oname, gname);
        if (scene.instances.back().shape < 0) {
            scene.shapes.push_back({});
            scene.shapes.back().name     = scene.instances.back().name;
            scene.shapes.back().material = get_material_id(matname);
            scene.instances.back().shape = (int)scene.shapes.size() - 1;
        }
        auto& shape = scene.shapes.back();
        add_verts(verts, shape);
        for (auto i = 0; i < verts.size(); i++)
            shape.points.push_back(vertex_map.at(verts[i]));
    };
    cb.object = [&](const string& name) {
        oname     = name;
        gname     = "";
        matname   = "";
        smoothing = true;
        add_instance(scene, oname, gname);
    };
    cb.group = [&](const string& name) {
        gname = name;
        if (split_group) {
            add_instance(scene, oname, gname);
        }
    };
    cb.smoothing = [&](const string& name) {
        smoothing = (name == "on");
        if (split_smoothing) {
            add_instance(scene, oname, gname);
        }
    };
    cb.usemtl = [&](const string& name) {
        matname = name;
        if (split_material) {
            add_instance(scene, oname, gname);
        }
    };
    cb.material = [&](const obj_material& omat) {
        auto material                   = yocto_material();
        material.name                   = omat.name;
        material.emission               = omat.ke;
        material.diffuse                = omat.kd;
        material.specular               = omat.ks;
        material.transmission           = omat.kt;
        material.roughness              = omat.rs;
        material.opacity                = omat.op;
        material.emission_texture       = add_texture(omat.ke_txt, false);
        material.diffuse_texture        = add_texture(omat.kd_txt, false);
        material.specular_texture       = add_texture(omat.ks_txt, false);
        material.transmission_texture   = add_texture(omat.kt_txt, false);
        material.opacity_texture        = add_texture(omat.op_txt, true);
        material.roughness_texture      = add_texture(omat.rs_txt, true);
        material.occlusion_texture      = add_texture(omat.occ_txt, true);
        material.bump_texture           = add_texture(omat.bump_txt, true);
        material.displacement_texture   = add_texture(omat.disp_txt, true);
        material.normal_texture         = add_texture(omat.norm_txt, true);
        material.volume_emission        = omat.ve;
        material.volume_albedo          = omat.va;
        material.volume_density         = omat.vd;
        material.volume_phaseg          = omat.vg;
        material.volume_density_texture = add_voltexture(omat.vd_txt, false);
        scene.materials.push_back(material);
        mmap[material.name] = (int)scene.materials.size() - 1;
    };
    cb.camera = [&](const obj_camera& ocam) {
        auto camera           = yocto_camera();
        camera.name           = ocam.name;
        camera.frame          = ocam.frame;
        camera.orthographic   = ocam.ortho;
        camera.film_width     = ocam.width;
        camera.film_height    = ocam.height;
        camera.focal_length   = ocam.focal;
        camera.focus_distance = ocam.focus;
        camera.lens_aperture  = ocam.aperture;
        scene.cameras.push_back(camera);
    };
    cb.environmnet = [&](const obj_environment& oenv) {
        auto environment             = yocto_environment();
        environment.name             = oenv.name;
        environment.frame            = oenv.frame;
        environment.emission         = oenv.ke;
        environment.emission_texture = add_texture(oenv.ke_txt, true);
        scene.environments.push_back(environment);
    };

    // Parse obj
    auto obj_options          = load_obj_options();
    obj_options.exit_on_error = options.exit_on_error;
    obj_options.geometry_only = false;
    if (!load_obj(filename, cb, obj_options)) return false;

    // cleanup empty
    for (auto idx = 0; idx < scene.instances.size(); idx++) {
        if (!is_instance_empty(scene, scene.instances[idx])) continue;
        scene.instances.erase(scene.instances.begin() + idx);
        idx--;
    }

    // merging quads and triangles
    for (auto& shape : scene.shapes) {
        if (empty(shape.triangles) || empty(shape.quads)) continue;
        merge_triangles_and_quads(shape.triangles, shape.quads, false);
    }

    // load textures
    auto dirname = get_dirname(filename);
    if (!load_scene_textures(scene, dirname, options)) return false;

    // fix scene
    scene.name = get_filename(filename);
    add_missing_cameras(scene);
    add_missing_materials(scene);
    add_missing_names(scene);
    update_transforms(scene);

    // done
    return true;
}

bool save_mtl(
    const string& filename, const yocto_scene& scene, bool flip_tr = true) {
    // open
    auto fs = open(filename, "wt");
    if (!fs) return false;

    // for each material, dump all the values
    for (auto& material : scene.materials) {
        print(fs, "newmtl {}\n", material.name);
        print(fs, "  illum 2\n");
        if (material.emission != zero3f)
            print(fs, "  Ke {}\n", material.emission);
        if (material.diffuse != zero3f)
            print(fs, "  Kd {}\n", material.diffuse);
        if (material.specular != zero3f)
            print(fs, "  Ks {}\n", material.specular);
        if (material.transmission != zero3f)
            print(fs, "  Kt {}\n", material.transmission);
        if (material.roughness != 1.0f)
            print(fs, "  Ns {}\n",
                (int)clamp(
                    2 / pow(clamp(material.roughness, 0.0f, 0.99f) + 1e-10f, 4.0f) -
                        2,
                    0.0f, 1.0e9f));
        if (material.opacity != 1.0f) print(fs, "  d {}\n", material.opacity);
        if (material.emission_texture >= 0)
            print(fs, "  map_Ke {}\n",
                scene.textures[material.emission_texture].filename);
        if (material.diffuse_texture >= 0)
            print(fs, "  map_Kd {}\n",
                scene.textures[material.diffuse_texture].filename);
        if (material.specular_texture >= 0)
            print(fs, "  map_Ks {}\n",
                scene.textures[material.specular_texture].filename);
        if (material.transmission_texture >= 0)
            print(fs, "  map_Kt {}\n",
                scene.textures[material.transmission_texture].filename);
        if (material.opacity_texture >= 0 &&
            material.opacity_texture != material.diffuse_texture)
            print(fs, "  map_d  {}\n",
                scene.textures[material.opacity_texture].filename);
        if (material.roughness_texture >= 0)
            print(fs, "  map_Pr {}\n",
                scene.textures[material.roughness_texture].filename);
        if (material.occlusion_texture >= 0)
            print(fs, "  map_occ {}\n",
                scene.textures[material.occlusion_texture].filename);
        if (material.bump_texture >= 0)
            print(fs, "  map_bump {}\n",
                scene.textures[material.bump_texture].filename);
        if (material.displacement_texture >= 0)
            print(fs, "  map_disp {}\n",
                scene.textures[material.displacement_texture].filename);
        if (material.normal_texture >= 0)
            print(fs, "  map_norm {}\n",
                scene.textures[material.normal_texture].filename);
        if (material.volume_emission != zero3f)
            print(fs, "  Ve {}\n", material.volume_emission);
        if (material.volume_density != zero3f)
            print(fs, "  Vd {}\n", material.volume_density);
        if (material.volume_albedo != zero3f)
            print(fs, "  Va {}\n", material.volume_albedo);
        if (material.volume_phaseg != 0)
            print(fs, "  Vg {}\n", material.volume_phaseg);
        if (material.volume_density_texture >= 0)
            print(fs, "  map_Vd {}\n",
                scene.voltextures[material.volume_density_texture].filename);
        print(fs, "\n");
    }

    // done
    return true;
}

bool save_objx(const string& filename, const yocto_scene& scene) {
    // scene
    auto fs = open(filename, "wt");
    if (!fs) return false;

    // cameras
    for (auto& camera : scene.cameras) {
        print(fs, "c {} {} {} {} {} {} {} {}\n", camera.name,
            (int)camera.orthographic, camera.film_width, camera.film_height,
            camera.focal_length, camera.focus_distance, camera.lens_aperture,
            camera.frame);
    }

    // environments
    for (auto& environment : scene.environments) {
        if (environment.emission_texture >= 0) {
            print(fs, "e {} {} {} {}\n", environment.name, environment.emission,
                scene.textures[environment.emission_texture].filename,
                environment.frame);
        } else {
            print(fs, "e {} {} \"\" {}\n", environment.name,
                environment.emission, environment.frame);
        }
    }

    // done
    return true;
}

string to_string(const obj_vertex& v) {
    auto s = std::to_string(v.position);
    if (v.texturecoord) {
        s += "/" + std::to_string(v.texturecoord);
        if (v.normal) s += "/" + std::to_string(v.normal);
    } else {
        if (v.normal) s += "//" + std::to_string(v.normal);
    }
    return s;
}

bool save_obj(const string& filename, const yocto_scene& scene,
    bool flip_texcoord = true) {
    // scene
    auto fs = open(filename, "wt");
    if (!fs) return false;

    // material library
    if (!empty(scene.materials)) {
        auto mtlname = replace_extension(get_filename(filename), "mtl");
        print(fs, "mtllib {}\n", mtlname);
    }

    // shapes
    auto offset = obj_vertex{0, 0, 0};
    for (auto& instance : scene.instances) {
        if (instance.surface >= 0) {
            auto& surface = scene.surfaces[instance.surface];
            print(fs, "o {}\n", instance.name);
            if (!empty(surface.materials) && empty(surface.quads_materials))
                print(fs, "usemtl {}\n",
                    scene.materials[surface.materials.front()].name);
            if (instance.frame == identity_frame3f) {
                for (auto& p : surface.positions) print(fs, "v {}\n", p);
                for (auto& n : surface.normals) print(fs, "vn {}\n", n);
                for (auto& t : surface.texturecoords)
                    print(fs, "vt {}\n",
                        vec2f{t.x, (flip_texcoord) ? 1 - t.y : t.y});
            } else {
                for (auto& pp : surface.positions) {
                    print(fs, "v {}\n", transform_point(instance.frame, pp));
                }
                for (auto& nn : surface.normals) {
                    print(fs, "vn {}\n", transform_direction(instance.frame, nn));
                }
                for (auto& t : surface.texturecoords)
                    print(fs, "vt {}\n",
                        vec2f{t.x, (flip_texcoord) ? 1 - t.y : t.y});
            }
            auto last_material_id = -1;
            for (auto i = 0; i < surface.quads_positions.size(); i++) {
                if (!empty(surface.quads_materials) &&
                    surface.quads_materials[i] != last_material_id) {
                    last_material_id = surface.quads_materials[i];
                    print(fs, "usemtl {}\n",
                        scene.materials[surface.materials[last_material_id]].name);
                }
                if (!empty(surface.texturecoords) && empty(surface.normals)) {
                    auto vert = [offset](int ip, int it) {
                        return obj_vertex{ip + offset.position + 1,
                            it + offset.texturecoord + 1, 0};
                    };
                    auto qp = surface.quads_positions[i];
                    auto qt = surface.quads_texturecoords[i];
                    if (qp.z == qp.w) {
                        print(fs, "f {} {} {}\n", to_string(vert(qp.x, qt.x)),
                            to_string(vert(qp.y, qt.y)),
                            to_string(vert(qp.z, qt.z)));
                    } else {
                        print(fs, "f {} {} {} {}\n", to_string(vert(qp.x, qt.x)),
                            to_string(vert(qp.y, qt.y)),
                            to_string(vert(qp.z, qt.z)),
                            to_string(vert(qp.w, qt.w)));
                    }
                } else if (!empty(surface.texturecoords) &&
                           !empty(surface.normals)) {
                    auto vert = [offset](int ip, int it, int in) {
                        return obj_vertex{ip + offset.position + 1,
                            it + offset.texturecoord + 1, in + offset.normal + 1};
                    };
                    auto qp = surface.quads_positions[i];
                    auto qt = surface.quads_texturecoords[i];
                    auto qn = surface.quads_normals[i];
                    if (qp.z == qp.w) {
                        print(fs, "f {} {} {}\n",
                            to_string(vert(qp.x, qt.x, qn.x)),
                            to_string(vert(qp.y, qt.y, qn.y)),
                            to_string(vert(qp.z, qt.z, qn.z)));
                    } else {
                        print(fs, "f {} {} {} {}\n",
                            to_string(vert(qp.x, qt.x, qn.x)),
                            to_string(vert(qp.y, qt.y, qn.y)),
                            to_string(vert(qp.z, qt.z, qn.z)),
                            to_string(vert(qp.w, qt.w, qn.w)));
                    }
                } else if (!empty(surface.normals)) {
                    auto vert = [offset](int ip, int in) {
                        return obj_vertex{ip + offset.position + 1, 0,
                            in + offset.normal + 1};
                    };
                    auto qp = surface.quads_positions[i];
                    auto qn = surface.quads_normals[i];
                    if (qp.z == qp.w) {
                        print(fs, "f {} {} {}\n", to_string(vert(qp.x, qn.x)),
                            to_string(vert(qp.y, qn.y)),
                            to_string(vert(qp.z, qn.z)));
                    } else {
                        print(fs, "f {} {} {} {}\n", to_string(vert(qp.x, qn.x)),
                            to_string(vert(qp.y, qn.y)),
                            to_string(vert(qp.z, qn.z)),
                            to_string(vert(qp.w, qn.w)));
                    }
                } else {
                    auto vert = [offset](int ip) {
                        return obj_vertex{ip + offset.position + 1, 0, 0};
                    };
                    auto q = surface.quads_positions[i];
                    if (q.z == q.w) {
                        print(fs, "f {} {} {}\n", to_string(vert(q.x)),
                            to_string(vert(q.y)), to_string(vert(q.z)));
                    } else {
                        print(fs, "f {} {} {} {}\n", to_string(vert(q.x)),
                            to_string(vert(q.y)), to_string(vert(q.z)),
                            to_string(vert(q.w)));
                    }
                }
            }
            offset.position += surface.positions.size();
            offset.texturecoord += surface.texturecoords.size();
            offset.normal += surface.normals.size();
        } else if (instance.shape >= 0) {
            auto& shape = scene.shapes[instance.shape];
            print(fs, "o {}\n", instance.name);
            if (shape.material >= 0)
                print(fs, "usemtl {}\n", scene.materials[shape.material].name);
            if (instance.frame == identity_frame3f) {
                for (auto& p : shape.positions) print(fs, "v {}\n", p);
                for (auto& n : shape.normals) print(fs, "vn {}\n", n);
                for (auto& t : shape.texturecoords)
                    print(fs, "vt {}\n",
                        vec2f{t.x, (flip_texcoord) ? 1 - t.y : t.y});
            } else {
                for (auto& pp : shape.positions) {
                    print(fs, "v {}\n", transform_point(instance.frame, pp));
                }
                for (auto& nn : shape.normals) {
                    print(fs, "vn {}\n", transform_direction(instance.frame, nn));
                }
                for (auto& t : shape.texturecoords)
                    print(fs, "vt {}\n",
                        vec2f{t.x, (flip_texcoord) ? 1 - t.y : t.y});
            }
            auto mask = obj_vertex{1, empty(shape.texturecoords) ? 0 : 1,
                empty(shape.normals) ? 0 : 1};
            auto vert = [mask, offset](int i) {
                return obj_vertex{(i + offset.position + 1) * mask.position,
                    (i + offset.texturecoord + 1) * mask.texturecoord,
                    (i + offset.normal + 1) * mask.normal};
            };
            for (auto& p : shape.points) {
                print(fs, "p {}\n", to_string(vert(p)));
            }
            for (auto& l : shape.lines) {
                print(fs, "l {} {}\n", to_string(vert(l.x)),
                    to_string(vert(l.y)));
            }
            for (auto& t : shape.triangles) {
                print(fs, "f {} {} {}\n", to_string(vert(t.x)),
                    to_string(vert(t.y)), to_string(vert(t.z)));
            }
            for (auto& q : shape.quads) {
                if (q.z == q.w) {
                    print(fs, "f {} {} {}\n", to_string(vert(q.x)),
                        to_string(vert(q.y)), to_string(vert(q.z)));
                } else {
                    print(fs, "f {} {} {} {}\n", to_string(vert(q.x)),
                        to_string(vert(q.y)), to_string(vert(q.z)),
                        to_string(vert(q.w)));
                }
            }
            offset.position += shape.positions.size();
            offset.texturecoord += shape.texturecoords.size();
            offset.normal += shape.normals.size();
        } else {
        }
    }

    return true;
}

bool save_obj_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options) {
    auto scope = log_trace_scoped("saving scene {}", filename);
    if (!save_obj(filename, scene, true)) return false;
    if (!empty(scene.materials)) {
        if (!save_mtl(replace_extension(filename, ".mtl"), scene, true))
            return false;
    }
    if (!empty(scene.cameras) || !empty(scene.environments)) {
        if (!save_objx(replace_extension(filename, ".objx"), scene))
            return false;
    }

    // skip textures if needed
    auto dirname = get_dirname(filename);
    if (!save_scene_textures(scene, dirname, options)) return false;

    // done
    return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GLTF CONVESION
// -----------------------------------------------------------------------------
namespace yocto {

static bool startswith(const string& str, const string& substr) {
    if (str.length() < substr.length()) return false;
    for (auto i = 0; i < substr.length(); i++)
        if (str[i] != substr[i]) return false;
    return true;
}

// convert gltf to scene
bool gltf_to_scene(yocto_scene& scene, const json& gltf, const string& dirname) {
    // convert textures
    if (gltf.count("images")) {
        for (auto instance_id = 0; instance_id < gltf.at("images").size();
             instance_id++) {
            auto& gimg       = gltf.at("images").at(instance_id);
            auto  texture    = yocto_texture{};
            texture.name     = gimg.value("name", ""s);
            texture.filename = (startswith(gimg.value("uri", ""s), "data:")) ?
                                   string("[glTF-inline].png") :
                                   gimg.value("uri", ""s);
            scene.textures.push_back(texture);
        }
    }

    // load buffers
    auto bmap = vector<vector<byte>>();
    if (gltf.count("buffers")) {
        bmap.resize(gltf.at("buffers").size());
        for (auto bid = 0; bid < gltf.at("buffers").size(); bid++) {
            auto& gbuf = gltf.at("buffers").at(bid);
            auto& data = bmap.at(bid);
            auto  uri  = gbuf.value("uri", ""s);
            if (uri == "") continue;
            if (startswith(uri, "data:")) {
                // assume it is base64 and find ','
                auto pos = uri.find(',');
                if (pos == uri.npos) {
                    return false;
                }
                // decode
                auto data_char = base64_decode(uri.substr(pos + 1));
                data = vector<unsigned char>((unsigned char*)data_char.c_str(),
                    (unsigned char*)data_char.c_str() + data_char.length());
            } else {
                auto filename = normalize_path(dirname + uri);
                if (!load_binary(filename, data)) return false;
            }
            if (gbuf.value("byteLength", -1) != data.size()) {
                return false;
            }
        }
    }

    // add a texture
    auto add_texture = [&scene, &gltf](const json& ginfo, bool force_linear) {
        if (!gltf.count("images") || !gltf.count("textures")) return -1;
        if (ginfo.is_null() || empty(ginfo)) return -1;
        if (ginfo.value("index", -1) < 0) return -1;
        auto& gtxt = gltf.at("textures").at(ginfo.value("index", -1));
        if (empty(gtxt) || gtxt.value("source", -1) < 0) return -1;
        auto texture_id = gtxt.value("source", -1);
        if (!gltf.count("samplers") || gtxt.value("sampler", -1) < 0)
            return texture_id;
        auto& gsmp = gltf.at("samplers").at(gtxt.value("sampler", -1));
        scene.textures[texture_id].clamp_to_edge = gsmp.value("wrapS", ""s) ==
                                                       "ClampToEdge" ||
                                                   gsmp.value("wrapT", ""s) ==
                                                       "ClampToEdge";
        scene.textures[texture_id].height_scale = gsmp.value("scale", 1.0f) *
                                                  gsmp.value("strength", 1.0f);
        scene.textures[texture_id]
            .ldr_as_linear = force_linear ||
                             is_hdr_filename(scene.textures[texture_id].filename);
        return texture_id;
    };

    // convert materials
    if (gltf.count("materials")) {
        for (auto mid = 0; mid < gltf.at("materials").size(); mid++) {
            auto& gmat        = gltf.at("materials").at(mid);
            auto  material    = yocto_material();
            material.name     = gmat.value("name", ""s);
            material.emission = gmat.value("emissiveFactor", zero3f);
            if (gmat.count("emissiveTexture"))
                material.emission_texture = add_texture(
                    gmat.at("emissiveTexture"), false);
            if (gmat.count("extensions") &&
                gmat.at("extensions").count("KHR_materials_pbrSpecularGlossiness")) {
                material.base_metallic = false;
                material.gltf_textures = true;
                auto& gsg              = gmat.at("extensions")
                                .at("KHR_materials_pbrSpecularGlossiness");
                auto kb = gsg.value("diffuseFactor", vec4f{1, 1, 1, 1});
                material.diffuse  = {kb.x, kb.y, kb.z};
                material.opacity  = kb.w;
                material.specular = gsg.value("specularFactor", vec3f{1, 1, 1});
                material.roughness = 1 - gsg.value("glossinessFactor", 1.0f);
                if (gsg.count("diffuseTexture"))
                    material.diffuse_texture = add_texture(
                        gsg.at("diffuseTexture"), false);
                if (gsg.count("specularGlossinessTexture"))
                    material.specular_texture = add_texture(
                        gsg.at("specularGlossinessTexture"), false);
                material.roughness_texture = material.specular_texture;
            } else if (gmat.count("pbrMetallicRoughness")) {
                material.base_metallic = true;
                material.gltf_textures = true;
                auto& gmr              = gmat.at("pbrMetallicRoughness");
                auto  kb = gmr.value("baseColorFactor", vec4f{1, 1, 1, 1});
                material.diffuse   = {kb.x, kb.y, kb.z};
                material.opacity   = kb.w;
                auto km            = gmr.value("metallicFactor", 1.0f);
                material.specular  = {km, km, km};
                material.roughness = gmr.value("roughnessFactor", 1.0f);
                if (gmr.count("baseColorTexture"))
                    material.diffuse_texture = add_texture(
                        gmr.at("baseColorTexture"), false);
                if (gmr.count("metallicRoughnessTexture"))
                    material.specular_texture = add_texture(
                        gmr.at("metallicRoughnessTexture"), true);
                material.roughness_texture = material.specular_texture;
            }
            if (gmat.count("occlusionTexture"))
                material.occlusion_texture = add_texture(
                    gmat.at("occlusionTexture"), true);
            if (gmat.count("normalTexture"))
                material.normal_texture = add_texture(
                    gmat.at("normalTexture"), true);
            material.double_sided = gmat.value("doubleSided", true);
            scene.materials.push_back(material);
        }
    }

    // get values from accessors
    auto accessor_values =
        [&gltf, &bmap](const json& gacc,
            bool normalize = false) -> vector<std::array<double, 4>> {
        auto gview  = gltf.at("bufferViews").at(gacc.value("bufferView", -1));
        auto data   = bmap.at(gview.value("buffer", -1)).data();
        auto offset = gacc.value("byteOffset", 0) + gview.value("byteOffset", 0);
        auto stride      = gview.value("byteStride", 0);
        auto compTypeNum = gacc.value("componentType", 5123);
        auto count       = gacc.value("count", -1);
        auto type        = gacc.value("type", ""s);
        auto ncomp       = 0;
        if (type == "SCALAR") ncomp = 1;
        if (type == "VEC2") ncomp = 2;
        if (type == "VEC3") ncomp = 3;
        if (type == "VEC4") ncomp = 4;
        auto compSize = 1;
        if (compTypeNum == 5122 || compTypeNum == 5123) {
            compSize = 2;
        }
        if (compTypeNum == 5124 || compTypeNum == 5125 || compTypeNum == 5126) {
            compSize = 4;
        }
        if (!stride) stride = compSize * ncomp;
        auto vals = vector<std::array<double, 4>>(count, {{0.0, 0.0, 0.0, 1.0}});
        for (auto i = 0; i < count; i++) {
            auto d = data + offset + i * stride;
            for (auto c = 0; c < ncomp; c++) {
                if (compTypeNum == 5120) {  // char
                    vals[i][c] = (double)(*(char*)d);
                    if (normalize) vals[i][c] /= SCHAR_MAX;
                } else if (compTypeNum == 5121) {  // byte
                    vals[i][c] = (double)(*(byte*)d);
                    if (normalize) vals[i][c] /= UCHAR_MAX;
                } else if (compTypeNum == 5122) {  // short
                    vals[i][c] = (double)(*(short*)d);
                    if (normalize) vals[i][c] /= SHRT_MAX;
                } else if (compTypeNum == 5123) {  // unsigned short
                    vals[i][c] = (double)(*(unsigned short*)d);
                    if (normalize) vals[i][c] /= USHRT_MAX;
                } else if (compTypeNum == 5124) {  // int
                    vals[i][c] = (double)(*(int*)d);
                    if (normalize) vals[i][c] /= INT_MAX;
                } else if (compTypeNum == 5125) {  // unsigned int
                    vals[i][c] = (double)(*(unsigned int*)d);
                    if (normalize) vals[i][c] /= UINT_MAX;
                } else if (compTypeNum == 5126) {  // float
                    vals[i][c] = (*(float*)d);
                }
                d += compSize;
            }
        }
        return vals;
    };

    // convert meshes
    auto meshes = vector<vector<int>>();
    if (gltf.count("meshes")) {
        for (auto mid = 0; mid < gltf.at("meshes").size(); mid++) {
            auto& gmesh = gltf.at("meshes").at(mid);
            meshes.push_back({});
            auto sid = 0;
            for (auto& gprim : gmesh.value("primitives", json::array())) {
                if (!gprim.count("attributes")) continue;
                auto shape = yocto_shape();
                shape.name = gmesh.value("name", ""s) +
                             ((sid) ? std::to_string(sid) : string());
                sid++;
                for (json::iterator gattr_it = gprim.at("attributes").begin();
                     gattr_it != gprim.at("attributes").end(); ++gattr_it) {
                    auto  semantic = gattr_it.key();
                    auto& gacc     = gltf.at("accessors")
                                     .at(gattr_it.value().get<int>());
                    auto vals = accessor_values(gacc);
                    if (semantic == "POSITION") {
                        shape.positions.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shape.positions.push_back({(float)vals[i][0],
                                (float)vals[i][1], (float)vals[i][2]});
                    } else if (semantic == "NORMAL") {
                        shape.normals.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shape.normals.push_back({(float)vals[i][0],
                                (float)vals[i][1], (float)vals[i][2]});
                    } else if (semantic == "TEXCOORD" ||
                               semantic == "TEXCOORD_0") {
                        shape.texturecoords.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shape.texturecoords.push_back(
                                {(float)vals[i][0], (float)vals[i][1]});
                    } else if (semantic == "COLOR" || semantic == "COLOR_0") {
                        shape.colors.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shape.colors.push_back(
                                {(float)vals[i][0], (float)vals[i][1],
                                    (float)vals[i][2], (float)vals[i][3]});
                    } else if (semantic == "TANGENT") {
                        shape.tangentspaces.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shape.tangentspaces.push_back(
                                {(float)vals[i][0], (float)vals[i][1],
                                    (float)vals[i][2], (float)vals[i][3]});
                        for (auto& t : shape.tangentspaces) t.w = -t.w;
                    } else if (semantic == "RADIUS") {
                        shape.radius.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shape.radius.push_back((float)vals[i][0]);
                    } else {
                        // ignore
                    }
                }
                // indices
                auto mode = gprim.value("mode", 4);
                if (!gprim.count("indices")) {
                    if (mode == 4) {
                        // triangles
                        shape.triangles.reserve(shape.positions.size() / 3);
                        for (auto i = 0; i < shape.positions.size() / 3; i++)
                            shape.triangles.push_back(
                                {i * 3 + 0, i * 3 + 1, i * 3 + 2});
                    } else if (mode == 6) {
                        // triangle fan
                        shape.triangles.reserve(shape.positions.size() - 2);
                        for (auto i = 2; i < shape.positions.size(); i++)
                            shape.triangles.push_back({0, i - 1, i});
                    } else if (mode == 5) {
                        // triangle strip
                        shape.triangles.reserve(shape.positions.size() - 2);
                        for (auto i = 2; i < shape.positions.size(); i++)
                            shape.triangles.push_back({i - 2, i - 1, i});
                    } else if (mode == 1) {
                        // lines
                        shape.lines.reserve(shape.positions.size() / 2);
                        for (auto i = 0; i < shape.positions.size() / 2; i++)
                            shape.lines.push_back({i * 2 + 0, i * 2 + 1});
                    } else if (mode == 2) {
                        // line loop
                        shape.lines.reserve(shape.positions.size());
                        for (auto i = 1; i < shape.positions.size(); i++)
                            shape.lines.push_back({i - 1, i});
                        shape.lines.back() = {(int)shape.positions.size() - 1, 0};
                    } else if (mode == 3) {
                        // line strip
                        shape.lines.reserve(shape.positions.size() - 1);
                        for (auto i = 1; i < shape.positions.size(); i++)
                            shape.lines.push_back({i - 1, i});
                    } else if (mode == -1 || mode == 0) {
                        // points
                        log_error("points not supported");
                    } else {
                        log_error("unknown primitive type");
                    }
                } else {
                    auto indices = accessor_values(
                        gltf.at("accessors").at(gprim.value("indices", -1)),
                        false);
                    if (mode == 4) {
                        // triangles
                        shape.triangles.reserve(indices.size());
                        for (auto i = 0; i < indices.size() / 3; i++)
                            shape.triangles.push_back({(int)indices[i * 3 + 0][0],
                                (int)indices[i * 3 + 1][0],
                                (int)indices[i * 3 + 2][0]});
                    } else if (mode == 6) {
                        // triangle fan
                        shape.triangles.reserve(indices.size() - 2);
                        for (auto i = 2; i < indices.size(); i++)
                            shape.triangles.push_back({(int)indices[0][0],
                                (int)indices[i - 1][0], (int)indices[i][0]});
                    } else if (mode == 5) {
                        // triangle strip
                        shape.triangles.reserve(indices.size() - 2);
                        for (auto i = 2; i < indices.size(); i++)
                            shape.triangles.push_back({(int)indices[i - 2][0],
                                (int)indices[i - 1][0], (int)indices[i][0]});
                    } else if (mode == 1) {
                        // lines
                        shape.lines.reserve(indices.size() / 2);
                        for (auto i = 0; i < indices.size() / 2; i++)
                            shape.lines.push_back({(int)indices[i * 2 + 0][0],
                                (int)indices[i * 2 + 1][0]});
                    } else if (mode == 2) {
                        // line loop
                        shape.lines.reserve(indices.size());
                        for (auto i = 1; i < indices.size(); i++)
                            shape.lines.push_back(
                                {(int)indices[i - 1][0], (int)indices[i][0]});
                        shape.lines.back() = {(int)indices[indices.size() - 1][0],
                            (int)indices[0][0]};
                    } else if (mode == 3) {
                        // line strip
                        shape.lines.reserve(indices.size() - 1);
                        for (auto i = 1; i < indices.size(); i++)
                            shape.lines.push_back(
                                {(int)indices[i - 1][0], (int)indices[i][0]});
                    } else if (mode == -1 || mode == 0) {
                        // points
                        log_error("points not supported");
                    } else {
                        log_error("unknown primitive type");
                    }
                }
                shape.material = gprim.count("material") ?
                                     gprim.value("material", -1) :
                                     -1;
                scene.shapes.push_back(shape);
                meshes.back().push_back((int)scene.shapes.size() - 1);
            }
        }
    }

    // convert cameras
    if (gltf.count("cameras")) {
        for (auto cid = 0; cid < gltf.at("cameras").size(); cid++) {
            auto& gcam          = gltf.at("cameras").at(cid);
            auto  camera        = yocto_camera{};
            camera.name         = gcam.value("name", ""s);
            camera.orthographic = gcam.value("type", ""s) == "orthographic";
            if (camera.orthographic) {
                log_error("orthographic not supported well");
                auto ortho = gcam.value("orthographic", json::object());
                camera.lens_aperture = 0;
                set_camera_perspective(camera, ortho.value("ymag", 0.0f),
                    ortho.value("xmag", 0.0f) / ortho.value("ymag", 0.0f),
                    float_max);
            } else {
                auto persp = gcam.value("perspective", json::object());
                camera.lens_aperture = 0;
                set_camera_perspective(camera, persp.value("yfov", 1.0f),
                    persp.value("aspectRatio", 1.0f), float_max);
            }
            scene.cameras.push_back(camera);
        }
    }

    // convert nodes
    if (gltf.count("nodes")) {
        for (auto nid = 0; nid < gltf.at("nodes").size(); nid++) {
            auto& gnde = gltf.at("nodes").at(nid);
            auto  node = yocto_scene_node{};
            node.name  = gnde.value("name", ""s);
            if (gnde.count("camera")) node.camera = gnde.value("camera", 0);
            node.translation = gnde.value("translation", zero3f);
            node.rotation    = gnde.value("rotation", vec4f{0, 0, 0, 1});
            node.scale       = gnde.value("scale", vec3f{1, 1, 1});
            node.local = mat_to_frame(gnde.value("matrix", identity_mat4f));
            scene.nodes.push_back(node);
        }

        // set up parent pointers
        for (auto nid = 0; nid < gltf.at("nodes").size(); nid++) {
            auto& gnde = gltf.at("nodes").at(nid);
            if (!gnde.count("children")) continue;
            for (auto& cid : gnde.at("children"))
                scene.nodes[cid.get<int>()].parent = nid;
        }

        // set up instances
        for (auto nid = 0; nid < gltf.at("nodes").size(); nid++) {
            auto& gnde = gltf.at("nodes").at(nid);
            if (!gnde.count("mesh")) continue;
            auto  node = scene.nodes[nid];
            auto& shps = meshes.at(gnde.value("mesh", 0));
            if (empty(shps)) continue;
            if (shps.size() == 1) {
                auto instance  = yocto_instance();
                instance.name  = node.name;
                instance.shape = shps[0];
                scene.instances.push_back(instance);
                node.instance = (int)scene.instances.size() - 1;
            } else {
                for (auto shp : shps) {
                    auto& shape    = scene.shapes[shp];
                    auto  instance = yocto_instance();
                    instance.name  = node.name + "_" + shape.name;
                    instance.shape = shp;
                    scene.instances.push_back(instance);
                    auto child     = yocto_scene_node{};
                    child.name     = node.name + "_" + shape.name;
                    child.parent   = nid;
                    child.instance = (int)scene.instances.size() - 1;
                    scene.nodes.push_back(child);
                }
            }
        }
    }

    // convert animations
    if (gltf.count("animations")) {
        for (auto& ganm : gltf.at("animations")) {
            auto aid         = 0;
            auto sampler_map = unordered_map<vec2i, int>();
            for (auto& gchannel : ganm.at("channels")) {
                auto path_ = gchannel.at("target").at("path").get<string>();
                auto path  = -1;
                if (path_ == "translation") path = 0;
                if (path_ == "rotation") path = 1;
                if (path_ == "scale") path = 2;
                if (path_ == "weights") path = 3;
                if (sampler_map.find({gchannel.at("sampler").get<int>(), path}) ==
                    sampler_map.end()) {
                    auto& gsampler = ganm.at("samplers")
                                         .at(gchannel.at("sampler").get<int>());
                    auto animation = yocto_animation{};
                    animation.name = (ganm.count("name") ?
                                             ganm.value("name", ""s) :
                                             "anim") +
                                     std::to_string(aid++);
                    animation.animation_group = ganm.value("name", ""s);
                    auto input_view           = accessor_values(
                        gltf.at("accessors").at(gsampler.value("input", -1)));
                    animation.keyframes_times.resize(input_view.size());
                    for (auto i = 0; i < input_view.size(); i++)
                        animation.keyframes_times[i] = input_view[i][0];
                    auto type = gsampler.value("interpolation", "LINEAR");
                    if (type == "LINEAR")
                        animation.interpolation_type = yocto_interpolation_type::linear;
                    if (type == "STEP")
                        animation.interpolation_type = yocto_interpolation_type::step;
                    if (type == "CUBICSPLINE")
                        animation.interpolation_type = yocto_interpolation_type::bezier;
                    auto output_view = accessor_values(
                        gltf.at("accessors").at(gsampler.value("output", -1)));
                    switch (path) {
                        case 0: {  // translation
                            animation.translation_keyframes.reserve(
                                output_view.size());
                            for (auto i = 0; i < output_view.size(); i++)
                                animation.translation_keyframes.push_back(
                                    {(float)output_view[i][0],
                                        (float)output_view[i][1],
                                        (float)output_view[i][2]});
                        } break;
                        case 1: {  // rotation
                            animation.rotation_keyframes.reserve(
                                output_view.size());
                            for (auto i = 0; i < output_view.size(); i++)
                                animation.rotation_keyframes.push_back(
                                    {(float)output_view[i][0],
                                        (float)output_view[i][1],
                                        (float)output_view[i][2],
                                        (float)output_view[i][3]});
                        } break;
                        case 2: {  // scale
                            animation.scale_keyframes.reserve(output_view.size());
                            for (auto i = 0; i < output_view.size(); i++)
                                animation.scale_keyframes.push_back(
                                    {(float)output_view[i][0],
                                        (float)output_view[i][1],
                                        (float)output_view[i][2]});
                        } break;
                        case 3: {  // weights
                            log_error("weights not supported for now");
#if 0
                        // get a node that it refers to
                        auto ncomp = 0;
                        auto gnode = gltf->get(gchannel->target->node);
                        auto gmesh = gltf->get(gnode->mesh);
                        if (gmesh) {
                            for (auto gshp : gmesh->primitives) {
                                ncomp = max((int)gshp->targets.size(), ncomp);
                            }
                        }
                        if (ncomp) {
                            auto values = vector<float>();
                            values.reserve(output_view.size());
                            for (auto i = 0; i < output_view.size(); i++)
                                values.push_back(output_view.get(i));
                            animation.weights.resize(values.size() / ncomp);
                            for (auto i = 0; i < animation.weights.size(); i++) {
                                animation.weights[i].resize(ncomp);
                                for (auto j = 0; j < ncomp; j++)
                                    animation.weights[i][j] = values[i * ncomp + j];
                            }
                        }
#endif
                        } break;
                        default: { return false; }
                    }
                    sampler_map[{gchannel.at("sampler").get<int>(),
                        path}] = (int)scene.animations.size();
                    scene.animations.push_back(animation);
                }
                scene
                    .animations[sampler_map.at(
                        {gchannel.at("sampler").get<int>(), path})]
                    .node_targets.push_back(
                        (int)gchannel.at("target").at("node").get<int>());
            }
        }
    }

    return true;
}

// Load a scene
bool load_gltf_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options) {
    // initialization
    scene = {};

    // convert json
    auto js = json();
    if (!load_json(filename, js)) return false;
    try {
        if (!gltf_to_scene(scene, js, get_dirname(filename))) return false;
    } catch (...) {
        return false;
    }

    // load textures
    auto dirname = get_dirname(filename);
    if (!load_scene_textures(scene, dirname, options)) return false;

    // fix scene
    scene.name = get_filename(filename);
    add_missing_cameras(scene);
    add_missing_materials(scene);
    add_missing_names(scene);
    update_transforms(scene);

    // fix cameras
    auto bbox = compute_scene_bounds(scene);
    for (auto& camera : scene.cameras) {
        auto center   = (bbox.min + bbox.max) / 2;
        auto distance = dot(-camera.frame.z, center - camera.frame.o);
        if (distance > 0) camera.focus_distance = distance;
    }

    // done
    return true;
}

// convert gltf scene to json
bool scene_to_gltf(const yocto_scene& scene, json& js) {
    // init to emprt object
    js = json::object();

    // start creating json
    js["asset"]["version"] = "2.0";

    // prepare top level nodes
    if (!empty(scene.cameras)) js["cameras"] = json::array();
    if (!empty(scene.textures)) {
        js["textures"] = json::array();
        js["images"]   = json::array();
    }
    if (!empty(scene.materials)) js["materials"] = json::array();
    if (!empty(scene.shapes)) {
        js["meshes"]      = json::array();
        js["buffers"]     = json::array();
        js["bufferViews"] = json::array();
        js["accessors"]   = json::array();
    }
    if (!empty(scene.instances)) js["nodes"] = json::array();
    if (!empty(scene.nodes)) js["nodes"] = json::array();

    // convert cameras
    for (auto& camera : scene.cameras) {
        auto cjs    = json();
        cjs["name"] = camera.name;
        if (!camera.orthographic) {
            cjs["type"]                       = "perspective";
            cjs["perspective"]["aspectRatio"] = camera.film_width /
                                                camera.film_height;
            cjs["perspective"]["znear"] = 0.01f;
        } else {
            cjs["type"]                  = "orthographic";
            cjs["orthographic"]["xmag"]  = camera.film_width / 2;
            cjs["orthographic"]["ymag"]  = camera.film_height / 2;
            cjs["orthographic"]["znear"] = 0.01f;
        }
        js["cameras"].push_back(cjs);
    }

    // textures
    for (auto& texture : scene.textures) {
        auto tjs = json(), ijs = json();
        tjs["source"] = (int)js["images"].size();
        ijs["uri"]    = texture.filename;
        js["images"].push_back(ijs);
        js["textures"].push_back(tjs);
    }

    // material
    for (auto& material : scene.materials) {
        auto mjs           = json();
        mjs["name"]        = material.name;
        mjs["doubleSided"] = material.double_sided;
        if (material.emission != zero3f)
            mjs["emissiveFactor"] = material.emission;
        if (material.emission_texture >= 0)
            mjs["emissiveTexture"]["index"] = material.emission_texture;
        auto kd = vec4f{material.diffuse.x, material.diffuse.y,
            material.diffuse.z, material.opacity};
        if (material.base_metallic) {
            auto mmjs               = json();
            mmjs["baseColorFactor"] = kd;
            mmjs["metallicFactor"]  = material.specular.x;
            mmjs["roughnessFactor"] = material.roughness;
            if (material.diffuse_texture >= 0)
                mmjs["baseColorTexture"]["index"] = material.diffuse_texture;
            if (material.specular_texture >= 0)
                mmjs["metallicRoughnessTexture"]["index"] = material.specular_texture;
            mjs["pbrMetallicRoughness"] = mmjs;
        } else {
            auto mmjs                = json();
            mmjs["diffuseFactor"]    = kd;
            mmjs["specularFactor"]   = material.specular;
            mmjs["glossinessFactor"] = 1 - material.roughness;
            if (material.diffuse_texture >= 0)
                mmjs["diffuseTexture"]["index"] = material.diffuse_texture;
            if (material.specular_texture >= 0)
                mmjs["specularGlossinessTexture"]["index"] = material.specular_texture;
            mjs["extensions"]["KHR_materials_pbrSpecularGlossiness"] = mmjs;
        }
        if (material.normal_texture >= 0)
            mjs["normalTexture"]["index"] = material.normal_texture;
        if (material.occlusion_texture >= 0)
            mjs["occlusionTexture"]["index"] = material.occlusion_texture;
        js["materials"].push_back(mjs);
    }

    // shapes
    for (auto& shape : scene.shapes) {
        auto mjs = json(), bjs = json(), pjs = json();
        auto bid          = js["buffers"].size();
        mjs["name"]       = shape.name;
        mjs["primitives"] = json::array();
        bjs["name"]       = shape.name;
        bjs["byteLength"] = 0;
        bjs["uri"]        = replace_extension(shape.filename, ".bin");
        pjs["material"]   = shape.material;
        auto add_accessor = [&js, &bjs, bid](
                                int count, string type, bool indices = false) {
            auto bytes = count * 4;
            if (type == "VEC2") bytes *= 2;
            if (type == "VEC3") bytes *= 3;
            if (type == "VEC4") bytes *= 4;
            auto ajs = json(), vjs = json();
            vjs["buffer"]        = bid;
            vjs["byteLength"]    = bytes;
            vjs["byteOffset"]    = bjs["byteLength"].get<int>();
            vjs["target"]        = (!indices) ? 34962 : 34963;
            bjs["byteLength"]    = bjs["byteLength"].get<int>() + bytes;
            ajs["bufferView"]    = (int)js["bufferViews"].size();
            ajs["byteOffset"]    = 0;
            ajs["componentType"] = (!indices) ? 5126 : 5125;
            ajs["count"]         = count;
            ajs["type"]          = type;
            js["accessors"].push_back(ajs);
            js["bufferViews"].push_back(vjs);
            return (int)js["accessors"].size() - 1;
        };
        auto nverts = (int)shape.positions.size();
        if (!empty(shape.positions))
            pjs["attributes"]["POSITION"] = add_accessor(nverts, "VEC3");
        if (!empty(shape.normals))
            pjs["attributes"]["NORMAL"] = add_accessor(nverts, "VEC3");
        if (!empty(shape.texturecoords))
            pjs["attributes"]["TEXCOORD_0"] = add_accessor(nverts, "VEC2");
        if (!empty(shape.colors))
            pjs["attributes"]["COLOR_0"] = add_accessor(nverts, "VEC4");
        if (!empty(shape.radius))
            pjs["attributes"]["RADIUS"] = add_accessor(nverts, "SCALAR");
        if (!empty(shape.points)) {
            pjs["indices"] = add_accessor(
                (int)shape.points.size(), "SCALAR", true);
            pjs["mode"] = 1;
        }
        if (!empty(shape.lines)) {
            pjs["indices"] = add_accessor(
                (int)shape.lines.size() * 2, "SCALAR", true);
            pjs["mode"] = 1;
        }
        if (!empty(shape.triangles)) {
            pjs["indices"] = add_accessor(
                (int)shape.triangles.size() * 3, "SCALAR", true);
            pjs["mode"] = 4;
        }
        if (!empty(shape.quads)) {
            auto triangles = convert_quads_to_triangles(shape.quads);
            pjs["indices"] = add_accessor(
                (int)triangles.size() * 3, "SCALAR", true);
            pjs["mode"] = 4;
        }
        mjs["primitives"].push_back(pjs);
        js["meshes"].push_back(mjs);
        js["buffers"].push_back(bjs);
    }

    // nodes
    for (auto& node : scene.nodes) {
        auto njs           = json();
        njs["name"]        = node.name;
        njs["matrix"]      = frame_to_mat(node.local);
        njs["translation"] = node.translation;
        njs["rotation"]    = node.rotation;
        njs["scale"]       = node.scale;
        if (node.camera >= 0) njs["camera"] = node.camera;
        if (node.instance >= 0)
            njs["mesh"] = scene.instances[node.instance].shape;
        if (!empty(node.children)) {
            njs["children"] = json::array();
            for (auto& c : node.children) njs["children"].push_back(c);
        }
        js["nodes"].push_back(njs);
    }

    // animations not supported yet
    if (!empty(scene.animations)) log_error("animation not supported yet");

    // nodes from instances
    if (empty(scene.nodes)) {
        auto camera_id = 0;
        for (auto& camera : scene.cameras) {
            auto njs      = json();
            njs["name"]   = camera.name;
            njs["camera"] = camera_id++;
            njs["matrix"] = frame_to_mat(camera.frame);
            js["nodes"].push_back(njs);
        }
        for (auto& instance : scene.instances) {
            auto njs      = json();
            njs["name"]   = instance.name;
            njs["mesh"]   = instance.shape;
            njs["matrix"] = frame_to_mat(instance.frame);
            js["nodes"].push_back(njs);
        }
    }

    // done
    return true;
}

// save gltf mesh
bool save_gltf_mesh(const string& filename, const yocto_shape& shape) {
    auto scope = log_trace_scoped("saving scene {}", filename);
    auto fs    = open(filename, "wb");
    if (!fs) return false;

    if (!write_values(fs, shape.positions)) return false;
    if (!write_values(fs, shape.normals)) return false;
    if (!write_values(fs, shape.texturecoords)) return false;
    if (!write_values(fs, shape.colors)) return false;
    if (!write_values(fs, shape.radius)) return false;
    if (!write_values(fs, shape.points)) return false;
    if (!write_values(fs, shape.lines)) return false;
    if (!write_values(fs, shape.triangles)) return false;
    auto qtriangles = convert_quads_to_triangles(shape.quads);
    if (!write_values(fs, qtriangles)) return false;

    return true;
}

// Save gltf json
bool save_gltf_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options) {
    // save json
    auto js = json();
    try {
        if (!scene_to_gltf(scene, js)) return false;
    } catch (...) {
        return false;
    }
    if (!save_json(filename, js)) return false;

    // meshes
    auto dirname = get_dirname(filename);
    for (auto& shape : scene.shapes) {
        if (shape.filename == "") continue;
        auto filename = normalize_path(dirname + shape.filename);
        filename      = replace_extension(filename, ".bin");
        if (!save_gltf_mesh(filename, shape)) {
            if (options.exit_on_error) return false;
        }
    }

    // save textures
    if (!save_scene_textures(scene, dirname, options)) return false;

    // done
    return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF PBRT
// -----------------------------------------------------------------------------
namespace yocto {

// convert pbrt to json
bool pbrt_to_json(const string& filename, json& js) {
    auto split = [](const string& str) {
        auto ret = vector<string>();
        if (empty(str)) return ret;
        auto lpos = (size_t)0;
        while (lpos != str.npos) {
            auto pos = str.find_first_of(" \t\n\r", lpos);
            if (pos != str.npos) {
                if (pos > lpos) ret.push_back(str.substr(lpos, pos - lpos));
                lpos = pos + 1;
            } else {
                if (lpos < str.size()) ret.push_back(str.substr(lpos));
                lpos = pos;
            }
        }
        return ret;
    };

    auto is_cmd = [](const vector<string>& tokens, int i) -> bool {
        auto& tok = tokens.at(i);
        return !(tok[0] == '[' || tok[0] == ']' || tok[0] == '\"' ||
                 tok[0] == '-' || tok[0] == '+' || tok[0] == '.' ||
                 std::isdigit(tok[0]));
    };
    auto is_number = [](const vector<string>& tokens, int i) -> bool {
        auto& tok = tokens.at(i);
        return tok[0] == '-' || tok[0] == '+' || tok[0] == '.' ||
               std::isdigit(tok[0]);
    };
    auto parse_string = [](const vector<string>& tokens, int& i) -> string {
        if (tokens[i][0] != '"') {
            log_error("string expected");
            return "";
        }
        auto tok = tokens[i++];
        tok      = tok.substr(1, tok.size() - 2);
        if (tok.find('|') != tok.npos) tok = tok.substr(tok.find('|') + 1);
        return tok;
    };
    auto parse_param = [&](const vector<string>& tokens, int& i,
                           json& js) -> void {
        auto list = false, first = true;
        while (i < tokens.size()) {
            if (is_cmd(tokens, i)) {
                break;
            } else if (tokens[i][0] == '[') {
                list = true;
                i++;
            } else if (tokens[i][0] == ']') {
                list = false;
                i++;
                break;
            } else if (tokens[i][0] == '"') {
                if (!first && !list) break;
                js.push_back(tokens[i].substr(1, tokens[i].size() - 2));
                i++;
                if (!list) break;
            } else {
                if (!first && !list) {
                    log_error("bad options");
                    break;
                }
                js.push_back(atof(tokens[i].c_str()));
                i++;
                if (!list) break;
            }
        }
    };
    auto parse_param_list = [&](const vector<string>& tokens, int& i,
                                json& js) -> void {
        while (i < tokens.size()) {
            if (is_cmd(tokens, i)) break;
            auto name = parse_string(tokens, i);
            js[name]  = json::array();
            parse_param(tokens, i, js.at(name));
            if (js.at(name).size() == 1) {
                js.at(name) = js.at(name).at(0);
            }
        }
    };
    auto parse_param_numbers = [&](const vector<string>& tokens, int& i,
                                   json& js) -> void {
        js["values"] = json::array();
        if (tokens[i][0] == '[') i++;
        while (is_number(tokens, i)) {
            js.at("values").push_back((float)atof(tokens[i++].c_str()));
        }
        if (tokens[i][0] == ']') i++;
    };

    auto fs = open(filename, "rt");
    if (!fs) return false;

    auto pbrt = ""s;
    auto line = ""s;
    while (read_line(fs, line)) {
        if (line.find('#') == line.npos)
            pbrt += line + "\n";
        else
            pbrt += line.substr(0, line.find('#')) + "\n";
    }

    auto re = std::regex("\"(\\w+)\\s+(\\w+)\"");
    pbrt    = std::regex_replace(pbrt, re, "\"$1|$2\"");
    pbrt    = std::regex_replace(pbrt, std::regex("\\["), " [ ");
    pbrt    = std::regex_replace(pbrt, std::regex("\\]"), " ] ");
    js      = json::array();

    auto tokens = split(pbrt);
    auto i      = 0;
    while (i < tokens.size()) {
        if (!is_cmd(tokens, i)) {
            log_error("command expected");
            break;
        }
        auto& tok   = tokens[i++];
        auto  jcmd  = json::object();
        jcmd["cmd"] = tok;
        if (tok == "Transform" || tok == "LookAt" || tok == "Scale" ||
            tok == "Rotate" || tok == "Translate" || tok == "ConcatTransform") {
            parse_param_numbers(tokens, i, jcmd);
        } else if (tok == "Integrator" || tok == "Sampler" ||
                   tok == "PixelFilter" || tok == "Film" || tok == "Camera" ||
                   tok == "Shape" || tok == "AreaLightSource" ||
                   tok == "LightSource") {
            jcmd["type"] = parse_string(tokens, i);
            parse_param_list(tokens, i, jcmd);
        } else if (tok == "Texture") {
            jcmd["name"]       = parse_string(tokens, i);
            jcmd["value_type"] = parse_string(tokens, i);
            jcmd["type"]       = parse_string(tokens, i);
            parse_param_list(tokens, i, jcmd);
        } else if (tok == "MakeNamedMaterial") {
            jcmd["name"] = parse_string(tokens, i);
            parse_param_list(tokens, i, jcmd);
        } else if (tok == "Material") {
            jcmd["type"] = parse_string(tokens, i);
            parse_param_list(tokens, i, jcmd);
        } else if (tok == "NamedMaterial" || tok == "ObjectBegin" ||
                   tok == "ObjectInstance") {
            jcmd["name"] = parse_string(tokens, i);
        } else if (tok == "WorldBegin" || tok == "AttributeBegin" ||
                   tok == "TransformBegin" || tok == "WorldEnd" ||
                   tok == "AttributeEnd" || tok == "TransformEnd" ||
                   tok == "ObjectEnd" || tok == "ReverseOrientation") {
        } else {
            log_error("unsupported command {}", tok);
        }
        js.push_back(jcmd);
    }
    // auto fstr = std::fstream(filename + ".json");
    // fstr << js;
    return true;
}

// Compute the fresnel term for dielectrics. Implementation from
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
vec3f pbrt_fresnel_dielectric(float cosw, const vec3f& eta_) {
    auto eta = eta_;
    if (cosw < 0) {
        eta  = vec3f{1, 1, 1} / eta;
        cosw = -cosw;
    }

    auto sin2 = 1 - cosw * cosw;
    auto eta2 = eta * eta;

    auto cos2t = vec3f{1, 1, 1} - vec3f{sin2, sin2, sin2} / eta2;
    if (cos2t.x < 0 || cos2t.y < 0 || cos2t.z < 0)
        return vec3f{1, 1, 1};  // tir

    auto t0 = vec3f{sqrt(cos2t.x), sqrt(cos2t.y), sqrt(cos2t.z)};
    auto t1 = eta * t0;
    auto t2 = eta * cosw;

    auto rs = (vec3f{cosw, cosw, cosw} - t1) / (vec3f{cosw, cosw, cosw} + t1);
    auto rp = (t0 - t2) / (t0 + t2);

    return (rs * rs + rp * rp) / 2.0f;
}

// Compute the fresnel term for metals. Implementation from
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
vec3f pbrt_fresnel_metal(float cosw, const vec3f& eta, const vec3f& etak) {
    if (etak == zero3f) return pbrt_fresnel_dielectric(cosw, eta);

    cosw       = clamp(cosw, (float)-1, (float)1);
    auto cos2  = cosw * cosw;
    auto sin2  = clamp(1 - cos2, (float)0, (float)1);
    auto eta2  = eta * eta;
    auto etak2 = etak * etak;

    auto t0         = eta2 - etak2 - vec3f{sin2, sin2, sin2};
    auto a2plusb2_2 = t0 * t0 + 4.0f * eta2 * etak2;
    auto a2plusb2   = vec3f{
        sqrt(a2plusb2_2.x), sqrt(a2plusb2_2.y), sqrt(a2plusb2_2.z)};
    auto t1  = a2plusb2 + vec3f{cos2, cos2, cos2};
    auto a_2 = (a2plusb2 + t0) / 2.0f;
    auto a   = vec3f{sqrt(a_2.x), sqrt(a_2.y), sqrt(a_2.z)};
    auto t2  = 2.0f * a * cosw;
    auto rs  = (t1 - t2) / (t1 + t2);

    auto t3 = vec3f{cos2, cos2, cos2} * a2plusb2 +
              vec3f{sin2, sin2, sin2} * vec3f{sin2, sin2, sin2};
    auto t4 = t2 * sin2;
    auto rp = rs * (t3 - t4) / (t3 + t4);

    return (rp + rs) / 2.0f;
}

// load pbrt scenes
bool load_pbrt_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options) {
    auto scope = log_trace_scoped("loading scene {}", filename);
    scene      = yocto_scene{};
    // convert to json
    auto js = json();
    try {
        if (!pbrt_to_json(filename, js)) return false;
    } catch (...) {
        return false;
    }

    auto dirname_ = get_dirname(filename);

    struct stack_item {
        frame3f        frame     = identity_frame3f;
        int            material  = -1;
        yocto_material light_mat = {};
        float          focus = 1, aspect = 1;
        bool           reverse = false;
    };

    // parse
    auto stack = vector<stack_item>();
    stack.push_back(stack_item());

    auto txt_map = unordered_map<string, int>();
    auto mat_map = unordered_map<string, int>();
    auto mid     = 0;

    auto get_vec3f = [](const json& js) -> vec3f {
        if (js.is_number())
            return {js.get<float>(), js.get<float>(), js.get<float>()};
        if (js.is_array() && js.size() == 1)
            return {js.at(0).get<float>(), js.at(0).get<float>(),
                js.at(0).get<float>()};
        if (js.is_array() && js.size() == 3)
            return {js.at(0).get<float>(), js.at(1).get<float>(),
                js.at(2).get<float>()};
        log_error("cannot handle vec3f");
        return zero3f;
    };

    auto get_emission_vec3f = [&get_vec3f](const json& js) -> vec3f {
        if (js.is_array() && js.size() == 2)
            return blackbody_to_rgb(js.at(0).get<float>()) *
                   js.at(1).get<float>();
        return get_vec3f(js);
    };

    auto get_vec4f = [](const json& js) -> vec4f {
        if (js.is_number())
            return {js.get<float>(), js.get<float>(), js.get<float>(),
                js.get<float>()};
        if (js.is_array() && js.size() == 4)
            return {js.at(0).get<float>(), js.at(1).get<float>(),
                js.at(2).get<float>(), js.at(3).get<float>()};
        log_error("cannot handle vec4f");
        return zero4f;
    };

    auto get_mat4f = [](const json& js) -> frame3f {
        if (!js.is_array() || js.size() != 16) {
            log_error("cannot handle vec4f");
            return identity_frame3f;
        }
        float m[16] = {0};
        for (auto i = 0; i < 16; i++) m[i] = js.at(i).get<float>();
        return {{m[0], m[1], m[2]}, {m[4], m[5], m[6]}, {m[8], m[9], m[10]},
            {m[12], m[13], m[14]}};
    };

    auto get_mat3f = [](const json& js) -> frame3f {
        if (!js.is_array() || js.size() != 9) {
            log_error("cannot handle mat3f");
            return identity_frame3f;
        }
        auto m = identity_frame3f;
        for (auto i = 0; i < 9; i++) (&m.x.x)[i] = js.at(i).get<float>();
        return m;
    };

    auto get_vector_vec3i = [](const json& js) -> vector<vec3i> {
        if (!js.is_array() || js.size() % 3) {
            log_error("cannot handle vector<vec3f");
            return {};
        }
        auto vals = vector<vec3i>(js.size() / 3);
        for (auto i = 0; i < vals.size(); i++) {
            vals[i].x = (int)std::round(js.at(i * 3 + 0).get<float>());
            vals[i].y = (int)std::round(js.at(i * 3 + 1).get<float>());
            vals[i].z = (int)std::round(js.at(i * 3 + 2).get<float>());
        }
        return vals;
    };

    auto get_vector_vec3f = [](const json& js) -> vector<vec3f> {
        if (!js.is_array() || js.size() % 3) {
            log_error("cannot handle vector<vec3f>");
            return {};
        }
        auto vals = vector<vec3f>(js.size() / 3);
        for (auto i = 0; i < vals.size(); i++) {
            vals[i].x = js.at(i * 3 + 0).get<float>();
            vals[i].y = js.at(i * 3 + 1).get<float>();
            vals[i].z = js.at(i * 3 + 2).get<float>();
        }
        return vals;
    };

    auto get_vector_vec2f = [](const json& js) -> vector<vec2f> {
        if (!js.is_array() || js.size() % 2) {
            log_error("cannot handle vector<vec3f>");
            return {};
        }
        auto vals = vector<vec2f>(js.size() / 2);
        for (auto i = 0; i < vals.size(); i++) {
            vals[i].x = js.at(i * 2 + 0).get<float>();
            vals[i].y = js.at(i * 2 + 1).get<float>();
        }
        return vals;
    };

    auto get_scaled_texture = [&txt_map, &get_vec3f](
                                  const json& js, vec3f& col, int& txt) {
        if (js.is_string()) {
            col = {1, 1, 1};
            txt = txt_map.at(js.get<string>());
        } else {
            col = get_vec3f(js);
            txt = -1;
        }
    };

    auto use_hierarchy = false;

    unordered_map<string, vector<yocto_instance>> objects;
    for (auto& jcmd : js) {
        auto cmd = jcmd.at("cmd").get<string>();
        if (cmd == "ObjectInstance") {
            use_hierarchy = true;
            break;
        }
    }

    auto lid = 0, sid = 0, cid = 0;
    auto cur_object = ""s;
    for (auto& jcmd : js) {
        auto cmd = jcmd.at("cmd").get<string>();
        if (cmd == "Integrator" || cmd == "Sampler" || cmd == "PixelFilter") {
        } else if (cmd == "Transform") {
            stack.back().frame = get_mat4f(jcmd.at("values"));
        } else if (cmd == "ConcatTransform") {
            stack.back().frame = stack.back().frame *
                                 get_mat4f(jcmd.at("values"));
        } else if (cmd == "Scale") {
            auto v             = get_vec3f(jcmd.at("values"));
            stack.back().frame = stack.back().frame * make_scaling_frame(v);
        } else if (cmd == "Translate") {
            auto v             = get_vec3f(jcmd.at("values"));
            stack.back().frame = stack.back().frame * make_translation_frame(v);
        } else if (cmd == "Rotate") {
            auto v             = get_vec4f(jcmd.at("values"));
            stack.back().frame = stack.back().frame *
                                 make_rotation_frame(
                                     vec3f{v.y, v.z, v.w}, v.x * pif / 180);
        } else if (cmd == "LookAt") {
            auto m             = get_mat3f(jcmd.at("values"));
            stack.back().frame = stack.back().frame *
                                 inverse(make_lookat_frame(m.x, m.y, m.z, true));
            stack.back().focus = length(m.x - m.y);
        } else if (cmd == "ReverseOrientation") {
            stack.back().reverse = !stack.back().reverse;
        } else if (cmd == "Film") {
            stack.back().aspect = jcmd.at("xresolution").get<float>() /
                                  jcmd.at("yresolution").get<float>();
        } else if (cmd == "Camera") {
            auto camera    = yocto_camera{};
            camera.name    = "camera" + std::to_string(cid++);
            camera.frame   = inverse(stack.back().frame);
            camera.frame.z = -camera.frame.z;
            auto focus     = stack.back().focus;
            auto aspect    = stack.back().aspect;
            auto fovy      = 1.0f;
            auto type      = jcmd.at("type").get<string>();
            if (type == "perspective") {
                fovy = jcmd.at("fov").get<float>() * pif / 180;
                if (aspect < 1) fovy = atan(tan(fovy) / aspect);
            } else {
                log_error("{} camera not supported", type);
            }
            set_camera_perspective(camera, fovy, aspect, focus);
            scene.cameras.push_back(camera);
        } else if (cmd == "Texture") {
            auto found = false;
            auto name  = jcmd.at("name").get<string>();
            for (auto& texture : scene.textures) {
                if (texture.name == name) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                scene.textures.push_back({});
                auto& texture         = scene.textures.back();
                texture.name          = jcmd.at("name").get<string>();
                txt_map[texture.name] = (int)scene.textures.size() - 1;
                auto type             = jcmd.at("type").get<string>();
                if (type == "imagemap") {
                    texture.filename = jcmd.at("filename").get<string>();
                    if (get_extension(texture.filename) == "pfm")
                        texture.filename = replace_extension(texture.filename,
                            ".hd"
                            "r");
                } else {
                    log_error("{} texture not supported", type);
                }
            }
        } else if (cmd == "MakeNamedMaterial" || cmd == "Material") {
            auto found = false;
            if (cmd == "MakeNamedMaterial") {
                auto name = jcmd.at("name").get<string>();
                for (auto& material : scene.materials) {
                    if (material.name == name) {
                        found = true;
                        break;
                    }
                }
            }
            if (!found) {
                scene.materials.push_back({});
                auto& material = scene.materials.back();
                if (cmd == "Material") {
                    material.name = "unnamed_mat" + std::to_string(mid++);
                    stack.back().material = (int)scene.materials.size() - 1;
                } else {
                    material.name          = jcmd.at("name").get<string>();
                    mat_map[material.name] = (int)scene.materials.size() - 1;
                }
                auto type = "uber"s;
                if (jcmd.count("type")) type = jcmd.at("type").get<string>();
                if (type == "uber") {
                    if (jcmd.count("Kd"))
                        get_scaled_texture(jcmd.at("Kd"), material.diffuse,
                            material.diffuse_texture);
                    if (jcmd.count("Ks"))
                        get_scaled_texture(jcmd.at("Ks"), material.specular,
                            material.specular_texture);
                    if (jcmd.count("Kt"))
                        get_scaled_texture(jcmd.at("Kt"), material.transmission,
                            material.transmission_texture);
                    if (jcmd.count("opacity")) {
                        auto op     = vec3f{0, 0, 0};
                        auto op_txt = -1;
                        get_scaled_texture(jcmd.at("opacity"), op, op_txt);
                        material.opacity         = (op.x + op.y + op.z) / 3;
                        material.opacity_texture = op_txt;
                    }
                    material.roughness = 0;
                } else if (type == "plastic") {
                    if (jcmd.count("Kd"))
                        get_scaled_texture(jcmd.at("Kd"), material.diffuse,
                            material.diffuse_texture);
                    if (jcmd.count("Ks"))
                        get_scaled_texture(jcmd.at("Ks"), material.specular,
                            material.specular_texture);
                    material.roughness = 0;
                } else if (type == "translucent") {
                    if (jcmd.count("Kd"))
                        get_scaled_texture(jcmd.at("Kd"), material.diffuse,
                            material.diffuse_texture);
                    if (jcmd.count("Ks"))
                        get_scaled_texture(jcmd.at("Ks"), material.specular,
                            material.specular_texture);
                    material.roughness = 0;
                } else if (type == "mix") {
                    auto matname1 = (jcmd.count("namedmaterial1")) ?
                                        jcmd.at("namedmaterial1").get<string>() :
                                        ""s;
                    auto matname2 = (jcmd.count("namedmaterial2")) ?
                                        jcmd.at("namedmaterial2").get<string>() :
                                        ""s;
                    // auto amount = get_vec3f(jcmd.at("amount"));
                    auto matname = (!matname1.empty()) ? matname1 : matname2;
                    for (auto& mat : scene.materials) {
                        if (mat.name == matname) {
                            material = mat;
                            break;
                        }
                    }
                } else if (type == "matte") {
                    material.diffuse = {1, 1, 1};
                    if (jcmd.count("Kd"))
                        get_scaled_texture(jcmd.at("Kd"), material.diffuse,
                            material.diffuse_texture);
                    material.roughness = 1;
                } else if (type == "mirror") {
                    material.diffuse   = {0, 0, 0};
                    material.specular  = {1, 1, 1};
                    material.roughness = 0;
                } else if (type == "metal") {
                    auto eta           = get_vec3f(jcmd.at("eta"));
                    auto k             = get_vec3f(jcmd.at("k"));
                    material.specular  = pbrt_fresnel_metal(1, eta, k);
                    material.roughness = 0;
                } else if (type == "substrate") {
                    if (jcmd.count("Kd"))
                        get_scaled_texture(jcmd.at("Kd"), material.diffuse,
                            material.diffuse_texture);
                    material.specular = {0.04f, 0.04f, 0.04f};
                    if (jcmd.count("Ks"))
                        get_scaled_texture(jcmd.at("Ks"), material.specular,
                            material.specular_texture);
                    material.roughness = 0;
                } else if (type == "glass") {
                    material.specular     = {0.04f, 0.04f, 0.04f};
                    material.transmission = {1, 1, 1};
                    if (jcmd.count("Ks"))
                        get_scaled_texture(jcmd.at("Ks"), material.specular,
                            material.specular_texture);
                    if (jcmd.count("Kt"))
                        get_scaled_texture(jcmd.at("Kt"), material.transmission,
                            material.transmission_texture);
                    material.roughness = 0;
                } else if (type == "mix") {
                    log_warning("mix material not properly supported");
                    if (jcmd.count("namedmaterial1")) {
                        auto mat1 = jcmd.at("namedmaterial1").get<string>();
                        auto saved_name = material.name;
                        material        = scene.materials[mat_map.at(mat1)];
                        material.name   = saved_name;
                    } else {
                        log_error("mix material missing front material");
                    }
                } else {
                    material.diffuse = {1, 0, 0};
                    log_error("{} material not supported", type);
                }
                if (jcmd.count("uroughness")) {
                    auto remap = js.count("remaproughness") &&
                                 js.at("remaproughness").get<bool>();
                    if (jcmd.count("uroughness"))
                        material.roughness = jcmd.at("uroughness").get<float>();
                    // if (!remap) material.rs = material.rs * material.rs;
                    if (remap) log_error("remap roughness not supported");
                }
                if (jcmd.count("roughness")) {
                    auto remap = js.count("remaproughness") &&
                                 js.at("remaproughness").get<bool>();
                    if (jcmd.count("roughness"))
                        material.roughness = jcmd.at("roughness").get<float>();
                    // if (!remap) material.rs = material.rs * material.rs;
                    if (remap) log_error("remap roughness not supported");
                }
                if (stack.back().light_mat.emission != zero3f) {
                    material.emission         = stack.back().light_mat.emission;
                    material.emission_texture = stack.back().light_mat.emission_texture;
                }
            }
        } else if (cmd == "NamedMaterial") {
            stack.back().material = mat_map.at(jcmd.at("name").get<string>());
            if (stack.back().light_mat.emission != zero3f) {
                auto material = yocto_material(
                    scene.materials[stack.back().material]);
                material.name += "_" + std::to_string(lid++);
                material.emission         = stack.back().light_mat.emission;
                material.emission_texture = stack.back().light_mat.emission_texture;
                scene.materials.push_back(material);
                stack.back().material = (int)scene.materials.size() - 1;
            }
        } else if (cmd == "Shape") {
            auto shape = yocto_shape();
            auto type  = jcmd.at("type").get<string>();
            if (type == "plymesh") {
                auto filename  = jcmd.at("filename").get<string>();
                shape.name     = get_filename(filename);
                shape.filename = filename;
                if (!options.skip_meshes) {
                    if (!load_ply_mesh(dirname_ + filename, shape.points,
                            shape.lines, shape.triangles, shape.quads,
                            shape.positions, shape.normals, shape.texturecoords,
                            shape.colors, shape.radius, false))
                        return false;
                }
            } else if (type == "trianglemesh") {
                shape.name     = "mesh" + std::to_string(sid++);
                shape.filename = "models/" + shape.name + ".ply";
                if (jcmd.count("indices"))
                    shape.triangles = get_vector_vec3i(jcmd.at("indices"));
                if (jcmd.count("P"))
                    shape.positions = get_vector_vec3f(jcmd.at("P"));
                if (jcmd.count("N"))
                    shape.normals = get_vector_vec3f(jcmd.at("N"));
                if (jcmd.count("uv")) {
                    shape.texturecoords = get_vector_vec2f(jcmd.at("uv"));
                    for (auto& uv : shape.texturecoords) uv.y = 1 - uv.y;
                }
            } else if (type == "sphere") {
                shape.name     = "sphere" + std::to_string(sid++);
                shape.filename = "models/" + shape.name + ".ply";
                auto radius    = 1.0f;
                if (jcmd.count("radius"))
                    radius = jcmd.at("radius").get<float>();
                tie(shape.quads, shape.positions, shape.normals,
                    shape.texturecoords) = make_sphere_shape({64, 32},
                    2 * radius, {1, 1});
            } else if (type == "disk") {
                shape.name     = "disk" + std::to_string(sid++);
                shape.filename = "models/" + shape.name + ".ply";
                auto radius    = 1.0f;
                if (jcmd.count("radius"))
                    radius = jcmd.at("radius").get<float>();
                tie(shape.quads, shape.positions, shape.normals,
                    shape.texturecoords) = make_disk_shape({32, 16}, 2 * radius,
                    {1, 1});
            } else {
                log_error("{} shape not supported", type);
            }
            auto frame = stack.back().frame;
            auto scl = vec3f{length(frame.x), length(frame.y), length(frame.z)};
            for (auto& p : shape.positions) p *= scl;
            frame = {normalize(frame.x), normalize(frame.y), normalize(frame.z),
                frame.o};
            if (stack.back().reverse) {
                for (auto& t : shape.triangles) swap(t.y, t.z);
            }
            shape.material = stack.back().material;
            scene.shapes.push_back(shape);
            auto instance  = yocto_instance();
            instance.name  = shape.name;
            instance.frame = frame;
            instance.shape = (int)scene.shapes.size() - 1;
            if (cur_object != "") {
                objects[cur_object].push_back(instance);
            } else {
                scene.instances.push_back(instance);
            }
        } else if (cmd == "ObjectInstance") {
            static auto instances = unordered_map<string, int>();
            auto        name      = jcmd.at("name").get<string>();
            auto&       object    = objects.at(name);
            for (auto shape : object) {
                instances[shape.name] += 1;
                auto instance = yocto_instance();
                instance.name = shape.name + "_ist" +
                                std::to_string(instances[shape.name]);
                instance.frame = stack.back().frame * shape.frame;
                instance.shape = shape.shape;
                scene.instances.push_back(instance);
            }
        } else if (cmd == "AreaLightSource") {
            auto type = jcmd.at("type").get<string>();
            if (type == "diffuse") {
                auto ligth_mat         = yocto_material();
                ligth_mat.emission     = get_emission_vec3f(jcmd.at("L"));
                stack.back().light_mat = ligth_mat;
            } else {
                log_error("{} area light not supported", type);
            }
        } else if (cmd == "LightSource") {
            auto type = jcmd.at("type").get<string>();
            if (type == "infinite") {
                auto environment = yocto_environment();
                environment.name = "environment" + std::to_string(lid++);
                // environment.frame =
                // frame3f{{1,0,0},{0,0,-1},{0,-1,0},{0,0,0}}
                // * stack.back().frame;
                environment.frame = stack.back().frame *
                                    frame3f{{1, 0, 0}, {0, 0, 1}, {0, 1, 0},
                                        {0, 0, 0}};
                environment.emission = {1, 1, 1};
                // log_info("stack frame: {}", stack.back().frame);
                // log_info("env   frame: {}", environment.frame);
                if (jcmd.count("scale"))
                    environment.emission *= get_vec3f(jcmd.at("scale"));
                if (jcmd.count("mapname")) {
                    auto texture     = yocto_texture{};
                    texture.filename = jcmd.at("mapname").get<string>();
                    texture.name     = environment.name;
                    scene.textures.push_back(texture);
                    environment.emission_texture = (int)scene.textures.size() - 1;
                }
                scene.environments.push_back(environment);
            } else if (type == "distant") {
                auto distant_dist = 100;
                scene.shapes.push_back({});
                auto& shape = scene.shapes.back();
                shape.name  = "distant" + std::to_string(lid++);
                auto from = vec3f{0, 0, 0}, to = vec3f{0, 0, 0};
                if (jcmd.count("from")) from = get_vec3f(jcmd.at("from"));
                if (jcmd.count("to")) to = get_vec3f(jcmd.at("to"));
                auto dir                 = normalize(from - to);
                auto size                = distant_dist * sin(5 * pif / 180);
                tie(shape.quads, shape.positions, shape.normals,
                    shape.texturecoords) = make_quad_shape({1, 1}, {size, size},
                    {1, 1});
                scene.materials.push_back({});
                auto& material    = scene.materials.back();
                shape.material    = scene.materials.size() - 1;
                material.name     = shape.name;
                material.emission = {1, 1, 1};
                if (jcmd.count("L"))
                    material.emission *= get_emission_vec3f(jcmd.at("L"));
                if (jcmd.count("scale"))
                    material.emission *= get_vec3f(jcmd.at("scale"));
                material.emission *= (distant_dist * distant_dist) /
                                     (size * size);
                auto instance  = yocto_instance();
                instance.name  = shape.name;
                instance.shape = (int)scene.shapes.size() - 1;
                instance.frame = stack.back().frame *
                                 make_lookat_frame(dir * distant_dist, zero3f,
                                     {0, 1, 0}, true);
                scene.instances.push_back(instance);
                log_warning("{} light not properly supported", type);
            } else {
                log_error("{} light not supported", type);
            }
        } else if (cmd == "WorldBegin") {
            stack.push_back(stack_item());
        } else if (cmd == "AttributeBegin") {
            stack.push_back(stack.back());
        } else if (cmd == "ObjectBegin") {
            auto name     = jcmd.at("name").get<string>();
            cur_object    = name;
            objects[name] = {};
        } else if (cmd == "ObjectEnd") {
            cur_object = "";
        } else if (cmd == "TransformBegin") {
            stack.push_back(stack.back());
        } else if (cmd == "WorldEnd" || cmd == "AttributeEnd" ||
                   cmd == "TransformEnd") {
            stack.pop_back();
        } else {
            log_error("{} command not supported", cmd.c_str());
        }
    }
    if (use_hierarchy) {
        auto camera_id = 0;
        for (auto& camera : scene.cameras) {
            auto node   = yocto_scene_node{};
            node.name   = camera.name;
            node.local  = camera.frame;
            node.camera = camera_id++;
            scene.nodes.insert(scene.nodes.begin(), node);
        }
        auto environment_id = 0;
        for (auto& environment : scene.environments) {
            auto node        = yocto_scene_node{};
            node.name        = environment.name;
            node.local       = environment.frame;
            node.environment = environment_id++;
            scene.nodes.push_back(node);
        }
    }

    // load textures
    auto dirname = get_dirname(filename);
    if (!load_scene_textures(scene, dirname, options)) return false;

    // fix scene
    scene.name = get_filename(filename);
    add_missing_cameras(scene);
    add_missing_materials(scene);
    add_missing_names(scene);
    update_transforms(scene);

    return true;
}

// Convert a scene to pbrt format
bool save_pbrt(const string& filename, const yocto_scene& scene) {
    auto scope = log_trace_scoped("saving scene {}", filename);
    auto fs    = open(filename, "wt");
    if (!fs) return false;

#if 0
WorldBegin

#uniform blue - ish illumination from all directions
LightSource "infinite" "rgb L" [.4 .45 .5]

#approximate the sun
LightSource "distant"  "point from" [ -30 40  100 ]
   "blackbody L" [3000 1.5]

AttributeBegin
  Material "glass"
  Shape "sphere" "float radius" 1
AttributeEnd

AttributeBegin
  Texture "checks" "spectrum" "checkerboard"
          "float uscale" [8] "float vscale" [8]
          "rgb tex1" [.1 .1 .1] "rgb tex2" [.8 .8 .8]
  Material "matte" "texture Kd" "checks"
  Translate 0 0 -1
  Shape "trianglemesh"
      "integer indices" [0 1 2 0 2 3]
      "point P" [ -20 -20 0   20 -20 0   20 20 0   -20 20 0 ]
      "float st" [ 0 0   1 0    1 1   0 1 ]
AttributeEnd

WorldEnd
#endif

    // convert camera and settings
    auto& camera         = scene.cameras.front();
    auto  from           = camera.frame.o;
    auto  to             = camera.frame.o - camera.frame.z;
    auto  up             = camera.frame.y;
    auto [width, height] = get_camera_image_size(camera, 0, 720);
    print(fs, "LookAt {} {} {}\n", from, to, up);
    print(fs, "Camera \"perspective\" \"float fov\" {}\n",
        get_camera_fovy(camera) * 180 / pif);

    // save renderer
    print(fs, "Sampler \"random\" \"integer pixelsamples\" [64]\n");
    // fprintf(f, "Sampler \"sobol\" \"interger pixelsamples\" [64]\n");
    print(fs, "Integrator \"path\"\n");
    print(fs,
        "Film \"image\" \"string filename\" [\"{}\"] "
        "\"integer xresolution\" [{}] \"integer yresolution\" [{}]\n",
        replace_extension(filename, "exr"), width, height);

    // start world
    print(fs, "WorldBegin\n");

    // convert textures
    for (auto& texture : scene.textures) {
        print(fs,
            "Texture \"{}\" \"spectrum\" \"imagemap\" "
            "\"string filename\" [\"{}\"]\n",
            texture.name, texture.filename);
    }

    // convert materials
    for (auto& material : scene.materials) {
        print(fs, "MakeNamedMaterial \"{}\" ", material.name);
        print(fs, "\"string type\" \"{}\" ", "uber");
        if (material.diffuse_texture >= 0)
            print(fs, "\"texture Kd\" [\"{}\"] ",
                scene.textures[material.diffuse_texture].name);
        else
            print(fs, "\"rgb Kd\" [{}] ", material.diffuse);
        if (material.specular_texture >= 0)
            print(fs, "\"texture Ks\" [\"{}\"] ",
                scene.textures[material.specular_texture].name);
        else
            print(fs, "\"rgb Ks\" [{}] ", material.specular);
        print(fs, "\"float roughness\" [{}] ", material.roughness);
        print(fs, "\n");
    }

    // convert instances
    for (auto& instance : scene.instances) {
        auto& shape    = scene.shapes[instance.shape];
        auto& material = scene.materials[shape.material];
        print(fs, "AttributeBegin\n");
        print(fs, "TransformBegin\n");
        print(fs, "ConcatTransform [{}]\n", frame_to_mat(instance.frame));
        if (material.emission != zero3f)
            print(fs, "AreaLightSource \"diffuse\" \"rgb L\" [ {} ]\n",
                material.emission);
        print(fs, "NamedMaterial \"{}\"\n", material.name);
        print(fs, "Shape \"plymesh\" \"string filename\" [\"{}\"]\n",
            shape.filename);
        print(fs, "TransformEnd\n");
        print(fs, "AttributeEnd\n");
    }

    // end world
    print(fs, "WorldEnd\n");

    // done
    return true;
}

// Save a pbrt scene
bool save_pbrt_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options) {
    // save json
    if (!save_pbrt(filename, scene)) return false;

    // save meshes
    auto dirname = get_dirname(filename);
    for (auto& shape : scene.shapes) {
        if (shape.filename == "") continue;
        auto filename = normalize_path(dirname + shape.filename);
        if (!save_mesh(filename, shape.points, shape.lines, shape.triangles,
                shape.quads, shape.positions, shape.normals,
                shape.texturecoords, shape.colors, shape.radius)) {
            if (options.exit_on_error) return false;
        }
    }

    // skip textures
    if (!save_scene_textures(scene, dirname, options)) return false;

    // done
    return true;
}

// Attempt to fix pbrt z-up.
void pbrt_flipyz_scene(yocto_scene& scene) {
    // flip meshes
    for (auto& shape : scene.shapes) {
        for (auto& p : shape.positions) swap(p.y, p.z);
        for (auto& n : shape.normals) swap(n.y, n.z);
    }
    for (auto& instance : scene.instances) {
        instance.frame = instance.frame *
                         frame3f{{1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 0, 0}};
    }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF BINARY SCENE FORMAT
// -----------------------------------------------------------------------------
namespace yocto {

// serialize_bin( ) can both save/load data to/from a binary file. The behaviour
// is set by the boolean 'save'. serialize_bin(name, file, true) : writes name
// as binary into file serialize_bin(name, file, false): read file as binary and
// set name

// Serialize type or struct with no allocated resource
template <typename T>
bool serialize_bin_value(T& value, file_stream& fs, bool save) {
    if (save) {
        return write_value(fs, value);
    } else {
        return read_value(fs, value);
    }
}

// Serialize vector
template <typename T>
bool serialize_bin_value(vector<T>& vec, file_stream& fs, bool save) {
    if (save) {
        auto count = (size_t)vec.size();
        if (!write_value(fs, count)) return false;
        if (!write_values(fs, vec)) return false;
        return true;
    } else {
        auto count = (size_t)0;
        if (!read_value(fs, count)) return false;
        vec = vector<T>(count);
        if (!read_values(fs, vec)) return false;
        return true;
    }
}

// Serialize string
bool serialize_bin_value(string& str, file_stream& fs, bool save) {
    if (save) {
        auto count = (size_t)str.size();
        if (!write_value(fs, count)) return false;
        auto vec = vector<char>(str.begin(), str.end());
        if (!write_values(fs, vec)) return false;
        return true;
    } else {
        auto count = (size_t)0;
        if (!read_value(fs, count)) return false;
        auto vec = vector<char>(count);
        if (!read_values(fs, vec)) return false;
        str = {vec.begin(), vec.end()};
        return true;
    }
}

// Serialize image
bool serialize_bin_value(image4f& img, file_stream& fs, bool save) {
    if (save) {
        if (!write_value(fs, img.width)) return false;
        if (!write_value(fs, img.height)) return false;
        if (!write_values(fs, img.width * img.height, data(img))) return false;
        return true;
    } else {
        if (!read_value(fs, img.width)) return false;
        if (!read_value(fs, img.height)) return false;
        img.pixels.resize(img.width * img.height);
        if (!read_values(fs, img.width * img.height, data(img))) return false;
        return true;
    }
}
bool serialize_bin_value(image4b& img, file_stream& fs, bool save) {
    if (save) {
        if (!write_value(fs, img.width)) return false;
        if (!write_value(fs, img.height)) return false;
        if (!write_values(fs, img.width * img.height, data(img))) return false;
        return true;
    } else {
        if (!read_value(fs, img.width)) return false;
        if (!read_value(fs, img.height)) return false;
        img.pixels.resize(img.width * img.height);
        if (!read_values(fs, img.width * img.height, data(img))) return false;
        return true;
    }
}

// Serialize image
bool serialize_bin_value(volume1f& vol, file_stream& fs, bool save) {
    if (save) {
        if (!write_value(fs, vol.width)) return false;
        if (!write_value(fs, vol.height)) return false;
        if (!write_value(fs, vol.depth)) return false;
        if (!write_values(fs, vol.width * vol.height * vol.depth, data(vol)))
            return false;
        return true;
    } else {
        if (!read_value(fs, vol.width)) return false;
        if (!read_value(fs, vol.height)) return false;
        if (!read_value(fs, vol.depth)) return false;
        vol.voxels.resize(vol.width * vol.height * vol.depth);
        if (!read_values(fs, vol.width * vol.height * vol.depth, data(vol)))
            return false;
        return true;
    }
}

// Serialize vector of pointers
template <typename T>
bool serialize_bin_object(vector<T*>& vec, file_stream& fs, bool save) {
    if (save) {
        auto count = (size_t)vec.size();
        if (!serialize_bin_value(count, fs, true)) return false;
        for (auto i = 0; i < vec.size(); ++i) {
            if (!serialize_bin_object(vec[i], fs, true)) return false;
        }
        return true;
    } else {
        auto count = (size_t)0;
        if (!serialize_bin_value(count, fs, false)) return false;
        vec = vector<T*>(count);
        for (auto i = 0; i < vec.size(); ++i) {
            vec[i] = new T();
            if (!serialize_bin_object(vec[i], fs, false)) return false;
        }
        return true;
    }
}

// Serialize vector of pointers
template <typename T>
bool serialize_bin_object(vector<T>& vec, file_stream& fs, bool save) {
    if (save) {
        auto count = (size_t)vec.size();
        if (!serialize_bin_value(count, fs, true)) return false;
        for (auto i = 0; i < vec.size(); ++i) {
            if (!serialize_bin_object(vec[i], fs, true)) return false;
        }
        return true;
    } else {
        auto count = (size_t)0;
        if (!serialize_bin_value(count, fs, false)) return false;
        vec = vector<T>(count);
        for (auto i = 0; i < vec.size(); ++i) {
            vec[i] = T{};
            if (!serialize_bin_object(vec[i], fs, false)) return false;
        }
        return true;
    }
}

// Serialize vector of objects
template <typename T>
bool serialize_bin_object(
    vector<T>& vec, const yocto_scene& scene, file_stream& fs, bool save) {
    if (save) {
        auto count = (size_t)vec.size();
        if (!serialize_bin_value(count, fs, true)) return false;
        for (auto i = 0; i < vec.size(); ++i) {
            if (!serialize_bin_object(vec[i], scene, fs, true)) return false;
        }
        return true;
    } else {
        auto count = (size_t)0;
        if (!serialize_bin_value(count, fs, false)) return false;
        vec = vector<T>(count);
        for (auto i = 0; i < vec.size(); ++i) {
            vec[i] = T{};
            if (!serialize_bin_object(vec[i], scene, fs, false)) return false;
        }
        return true;
    }
}

// Serialize yocto types. This is mostly boiler plate code.
bool serialize_bin_object(yocto_camera& camera, file_stream& fs, bool save) {
    if (!serialize_bin_value(camera.name, fs, save)) return false;
    if (!serialize_bin_value(camera.frame, fs, save)) return false;
    if (!serialize_bin_value(camera.orthographic, fs, save)) return false;
    if (!serialize_bin_value(camera.film_width, fs, save)) return false;
    if (!serialize_bin_value(camera.film_height, fs, save)) return false;
    if (!serialize_bin_value(camera.focal_length, fs, save)) return false;
    if (!serialize_bin_value(camera.focus_distance, fs, save)) return false;
    if (!serialize_bin_value(camera.lens_aperture, fs, save)) return false;
    return true;
}

bool serialize_bin_object(bvh_shape& bvh, file_stream& fs, bool save) {
    if (!serialize_bin_value(bvh.positions, fs, save)) return false;
    if (!serialize_bin_value(bvh.radius, fs, save)) return false;
    if (!serialize_bin_value(bvh.points, fs, save)) return false;
    if (!serialize_bin_value(bvh.lines, fs, save)) return false;
    if (!serialize_bin_value(bvh.triangles, fs, save)) return false;
    if (!serialize_bin_value(bvh.quads, fs, save)) return false;
    if (!serialize_bin_value(bvh.nodes, fs, save)) return false;
    if (!serialize_bin_value(bvh.nodes, fs, save)) return false;
    return true;
}

bool serialize_bin_object(bvh_scene& bvh, file_stream& fs, bool save) {
    if (!serialize_bin_value(bvh.nodes, fs, save)) return false;
    if (!serialize_bin_value(bvh.instances, fs, save)) return false;
    if (!serialize_bin_object(bvh.shape_bvhs, fs, save)) return false;
    if (!serialize_bin_value(bvh.nodes, fs, save)) return false;
    return true;
}

bool serialize_bin_object(
    yocto_shape& shape, const yocto_scene& scene, file_stream& fs, bool save) {
    if (!serialize_bin_value(shape.name, fs, save)) return false;
    if (!serialize_bin_value(shape.filename, fs, save)) return false;
    if (!serialize_bin_value(shape.material, fs, save)) return false;
    if (!serialize_bin_value(shape.subdivision_level, fs, save)) return false;
    if (!serialize_bin_value(shape.catmull_clark, fs, save)) return false;
    if (!serialize_bin_value(shape.compute_vertex_normals, fs, save))
        return false;
    if (!serialize_bin_value(shape.points, fs, save)) return false;
    if (!serialize_bin_value(shape.lines, fs, save)) return false;
    if (!serialize_bin_value(shape.triangles, fs, save)) return false;
    if (!serialize_bin_value(shape.quads, fs, save)) return false;
    if (!serialize_bin_value(shape.positions, fs, save)) return false;
    if (!serialize_bin_value(shape.normals, fs, save)) return false;
    if (!serialize_bin_value(shape.texturecoords, fs, save)) return false;
    if (!serialize_bin_value(shape.colors, fs, save)) return false;
    if (!serialize_bin_value(shape.radius, fs, save)) return false;
    if (!serialize_bin_value(shape.tangentspaces, fs, save)) return false;
    return true;
}

bool serialize_bin_object(yocto_surface& surface, const yocto_scene& scene,
    file_stream& fs, bool save) {
    if (!serialize_bin_value(surface.name, fs, save)) return false;
    if (!serialize_bin_value(surface.filename, fs, save)) return false;
    if (!serialize_bin_value(surface.materials, fs, save)) return false;
    if (!serialize_bin_value(surface.subdivision_level, fs, save)) return false;
    if (!serialize_bin_value(surface.catmull_clark, fs, save)) return false;
    if (!serialize_bin_value(surface.compute_vertex_normals, fs, save))
        return false;
    if (!serialize_bin_value(surface.quads_positions, fs, save)) return false;
    if (!serialize_bin_value(surface.quads_normals, fs, save)) return false;
    if (!serialize_bin_value(surface.quads_texturecoords, fs, save))
        return false;
    if (!serialize_bin_value(surface.quads_materials, fs, save)) return false;
    if (!serialize_bin_value(surface.positions, fs, save)) return false;
    if (!serialize_bin_value(surface.normals, fs, save)) return false;
    if (!serialize_bin_value(surface.texturecoords, fs, save)) return false;
    return true;
}

bool serialize_bin_object(yocto_texture& texture, file_stream& fs, bool save) {
    if (!serialize_bin_value(texture.name, fs, save)) return false;
    if (!serialize_bin_value(texture.filename, fs, save)) return false;
    if (!serialize_bin_value(texture.hdr_image, fs, save)) return false;
    if (!serialize_bin_value(texture.ldr_image, fs, save)) return false;
    if (!serialize_bin_value(texture.clamp_to_edge, fs, save)) return false;
    if (!serialize_bin_value(texture.height_scale, fs, save)) return false;
    if (!serialize_bin_value(texture.no_interpolation, fs, save)) return false;
    if (!serialize_bin_value(texture.ldr_as_linear, fs, save)) return false;
    if (!serialize_bin_value(texture.has_opacity, fs, save)) return false;
    return true;
}

bool serialize_bin_object(yocto_voltexture& texture, file_stream& fs, bool save) {
    if (!serialize_bin_value(texture.name, fs, save)) return false;
    if (!serialize_bin_value(texture.filename, fs, save)) return false;
    if (!serialize_bin_value(texture.volume_data, fs, save)) return false;
    if (!serialize_bin_value(texture.clamp_to_edge, fs, save)) return false;
    return true;
}

bool serialize_bin_object(yocto_environment& environment,
    const yocto_scene& scene, file_stream& fs, bool save) {
    if (!serialize_bin_value(environment.name, fs, save)) return false;
    if (!serialize_bin_value(environment.frame, fs, save)) return false;
    if (!serialize_bin_value(environment.emission, fs, save)) return false;
    if (!serialize_bin_value(environment.emission_texture, fs, save))
        return false;
    return true;
}

bool serialize_bin_object(yocto_material& material, const yocto_scene& scene,
    file_stream& fs, bool save) {
    if (!serialize_bin_value(material.name, fs, save)) return false;
    if (!serialize_bin_value(material.base_metallic, fs, save)) return false;
    if (!serialize_bin_value(material.gltf_textures, fs, save)) return false;
    if (!serialize_bin_value(material.double_sided, fs, save)) return false;
    if (!serialize_bin_value(material.emission, fs, save)) return false;
    if (!serialize_bin_value(material.diffuse, fs, save)) return false;
    if (!serialize_bin_value(material.specular, fs, save)) return false;
    if (!serialize_bin_value(material.transmission, fs, save)) return false;
    if (!serialize_bin_value(material.roughness, fs, save)) return false;
    if (!serialize_bin_value(material.opacity, fs, save)) return false;
    if (!serialize_bin_value(material.fresnel, fs, save)) return false;
    if (!serialize_bin_value(material.refract, fs, save)) return false;
    if (!serialize_bin_value(material.emission_texture, fs, save)) return false;
    if (!serialize_bin_value(material.diffuse_texture, fs, save)) return false;
    if (!serialize_bin_value(material.specular_texture, fs, save)) return false;
    if (!serialize_bin_value(material.transmission_texture, fs, save))
        return false;
    if (!serialize_bin_value(material.roughness_texture, fs, save))
        return false;
    if (!serialize_bin_value(material.opacity_texture, fs, save)) return false;
    if (!serialize_bin_value(material.occlusion_texture, fs, save))
        return false;
    if (!serialize_bin_value(material.bump_texture, fs, save)) return false;
    if (!serialize_bin_value(material.displacement_texture, fs, save))
        return false;
    if (!serialize_bin_value(material.normal_texture, fs, save)) return false;
    if (!serialize_bin_value(material.volume_emission, fs, save)) return false;
    if (!serialize_bin_value(material.volume_albedo, fs, save)) return false;
    if (!serialize_bin_value(material.volume_density, fs, save)) return false;
    if (!serialize_bin_value(material.volume_phaseg, fs, save)) return false;
    if (!serialize_bin_value(material.volume_density_texture, fs, save))
        return false;
    return true;
};

bool serialize_bin_object(yocto_instance& instance, const yocto_scene& scene,
    file_stream& fs, bool save) {
    if (!serialize_bin_value(instance.name, fs, save)) return false;
    if (!serialize_bin_value(instance.frame, fs, save)) return false;
    if (!serialize_bin_value(instance.shape, fs, save)) return false;
    if (!serialize_bin_value(instance.surface, fs, save)) return false;
    return true;
};

bool serialize_scene(yocto_scene& scene, file_stream& fs, bool save) {
    if (!serialize_bin_value(scene.name, fs, save)) return false;
    if (!serialize_bin_object(scene.cameras, fs, save)) return false;
    if (!serialize_bin_object(scene.shapes, scene, fs, save)) return false;
    if (!serialize_bin_object(scene.surfaces, scene, fs, save)) return false;
    if (!serialize_bin_object(scene.textures, fs, save)) return false;
    if (!serialize_bin_object(scene.voltextures, fs, save)) return false;
    if (!serialize_bin_object(scene.materials, scene, fs, save)) return false;
    if (!serialize_bin_object(scene.instances, scene, fs, save)) return false;
    if (!serialize_bin_object(scene.environments, scene, fs, save))
        return false;
    return true;
}

// Load/save a binary dump useful for very fast scene IO.
bool load_ybin_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options) {
    auto scope = log_trace_scoped("loading scene {}", filename);
    auto fs    = open(filename, "rb");
    if (!fs) return false;
    scene = {};
    if (!serialize_scene(scene, fs, false)) return false;
    return true;
}

// Load/save a binary dump useful for very fast scene IO.
bool save_ybin_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options) {
    auto fs = open(filename, "wb");
    if (!fs) return false;
    if (!serialize_scene((yocto_scene&)scene, fs, true)) return false;
    return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Reset mesh data
void reset_mesh_data(vector<int>& points, vector<vec2i>& lines,
    vector<vec3i>& triangles, vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texturecoords, vector<vec4f>& colors,
    vector<float>& radius) {
    points        = {};
    lines         = {};
    triangles     = {};
    quads         = {};
    positions     = {};
    normals       = {};
    texturecoords = {};
    colors        = {};
    radius        = {};
}

// Load ply mesh
bool load_mesh(const string& filename, vector<int>& points, vector<vec2i>& lines,
    vector<vec3i>& triangles, vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texturecoords, vector<vec4f>& colors,
    vector<float>& radius, bool force_triangles) {
    auto ext = get_extension(filename);
    if (ext == "ply" || ext == "PLY") {
        return load_ply_mesh(filename, points, lines, triangles, quads,
            positions, normals, texturecoords, colors, radius, force_triangles);
    } else if (ext == "obj" || ext == "OBJ") {
        return load_obj_mesh(filename, points, lines, triangles, quads,
            positions, normals, texturecoords, force_triangles);
    } else {
        reset_mesh_data(points, lines, triangles, quads, positions, normals,
            texturecoords, colors, radius);
        return false;
    }
}

// Save ply mesh
bool save_mesh(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    const vector<vec4f>& colors, const vector<float>& radius, bool ascii) {
    auto ext = get_extension(filename);
    if (ext == "ply" || ext == "PLY") {
        return save_ply_mesh(filename, points, lines, triangles, quads,
            positions, normals, texturecoords, colors, radius, ascii);
    } else if (ext == "obj" || ext == "OBJ") {
        return save_obj_mesh(filename, points, lines, triangles, quads,
            positions, normals, texturecoords);
    } else {
        return false;
    }
}

bool load_ply_mesh(const string& filename, vector<int>& points,
    vector<vec2i>& lines, vector<vec3i>& triangles, vector<vec4i>& quads,
    vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texturecoords, vector<vec4f>& color, vector<float>& radius,
    bool force_triangles, bool flip_texcoord) {
    // clear
    reset_mesh_data(points, lines, triangles, quads, positions, normals,
        texturecoords, color, radius);

    // load ply
    auto ply = ply_data{};
    if (!load_ply(filename, ply)) {
        log_io_error("empty ply file {}", filename);
        return false;
    }

    // copy vertex data
    for (auto& elem : ply.elements) {
        if (elem.name != "vertex") continue;
        auto count = elem.count;
        for (auto& prop : elem.properties) {
            auto vals        = data(prop.scalars);
            auto copy_floats = [vals, count](auto& vert, const auto& def,
                                   int stride, int offset) {
                if (vert.size() != count) vert.resize(count, def);
                auto dst = (float*)data(vert);
                for (auto i = 0; i < count; i++)
                    dst[i * stride + offset] = vals[i];
            };
            if (prop.name == "x") copy_floats(positions, zero3f, 3, 0);
            if (prop.name == "y") copy_floats(positions, zero3f, 3, 1);
            if (prop.name == "z") copy_floats(positions, zero3f, 3, 2);
            if (prop.name == "nx") copy_floats(normals, zero3f, 3, 0);
            if (prop.name == "ny") copy_floats(normals, zero3f, 3, 1);
            if (prop.name == "nz") copy_floats(normals, zero3f, 3, 2);
            if (prop.name == "u" || prop.name == "s")
                copy_floats(texturecoords, zero2f, 2, 0);
            if (prop.name == "v" || prop.name == "t")
                copy_floats(texturecoords, zero2f, 2, 1);
            if (prop.name == "red") copy_floats(color, vec4f{0, 0, 0, 1}, 4, 0);
            if (prop.name == "green")
                copy_floats(color, vec4f{0, 0, 0, 1}, 4, 1);
            if (prop.name == "blue")
                copy_floats(color, vec4f{0, 0, 0, 1}, 4, 2);
            if (prop.name == "alpha")
                copy_floats(color, vec4f{0, 0, 0, 1}, 4, 3);
            if (prop.name == "radius") copy_floats(radius, 0.0f, 1, 0);
        }
    }

    // fix texture coordinated
    if (flip_texcoord && !empty(texturecoords)) {
        for (auto& uv : texturecoords) uv.y = 1 - uv.y;
    }

    // copy face data
    for (auto& elem : ply.elements) {
        if (elem.name != "face") continue;
        auto count = elem.count;
        for (auto& prop : elem.properties) {
            if (prop.name == "vertex_indices") {
                for (auto fid = 0; fid < count; fid++) {
                    auto& list = prop.lists[fid];
                    auto  num  = (int)prop.scalars[fid];
                    if (num == 4) {
                        quads.push_back({list[0], list[1], list[2], list[3]});
                    } else {
                        for (auto i = 2; i < num; i++)
                            triangles.push_back({list[0], list[i - 1], list[i]});
                    }
                }
            }
        }
    }

    // copy face data
    for (auto& elem : ply.elements) {
        if (elem.name != "line") continue;
        auto count = elem.count;
        for (auto& prop : elem.properties) {
            if (prop.name == "vertex_indices") {
                for (auto fid = 0; fid < count; fid++) {
                    auto& list = prop.lists[fid];
                    auto  num  = (int)prop.scalars[fid];
                    for (auto i = 1; i < num; i++)
                        lines.push_back({list[i], list[i - 1]});
                }
            }
        }
    }

    merge_triangles_and_quads(triangles, quads, force_triangles);

    // done
    return true;
}

// Save ply mesh
bool save_ply_mesh(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    const vector<vec4f>& colors, const vector<float>& radius, bool ascii,
    bool flip_texcoord) {
    auto fs = open(filename, "wb");
    if (!fs) return false;

    // header
    print(fs, "ply\n");
    if (ascii)
        print(fs, "format ascii 1.0\n");
    else
        print(fs, "format binary_little_endian 1.0\n");
    print(fs, "element vertex {}\n", (int)positions.size());
    if (!empty(positions))
        print(fs, "property float x\nproperty float y\nproperty float z\n");
    if (!empty(normals))
        print(fs, "property float nx\nproperty float ny\nproperty float nz\n");
    if (!empty(texturecoords))
        print(fs, "property float u\nproperty float v\n");
    if (!empty(colors))
        print(fs,
            "property float red\nproperty float green\nproperty float blue\nproperty float alpha\n");
    if (!empty(radius)) print(fs, "property float radius\n");
    if (!empty(triangles) || !empty(quads)) {
        print(fs, "element face {}\n", (int)triangles.size() + (int)quads.size());
        print(fs, "property list uchar int vertex_indices\n");
    }
    if (!empty(lines)) {
        print(fs, "element line {}\n", (int)lines.size());
        print(fs, "property list uchar int vertex_indices\n");
    }
    print(fs, "end_header\n");

    // body
    if (ascii) {
        // write vertex data
        for (auto i = 0; i < positions.size(); i++) {
            if (!empty(positions)) print(fs, "{} ", positions[i]);
            if (!empty(normals)) print(fs, "{} ", normals[i]);
            if (!empty(texturecoords))
                print(fs, "{} ",
                    (!flip_texcoord) ?
                        texturecoords[i] :
                        vec2f{texturecoords[i].x, 1 - texturecoords[i].y});
            if (!empty(colors)) print(fs, "{} ", colors[i]);
            if (!empty(radius)) print(fs, "{} ", radius[i]);
            print(fs, "\n");
        }

        // write face data
        for (auto& t : triangles) print(fs, "3 {}\n", t);
        for (auto& q : quads) {
            if (q.z == q.w)
                print(fs, "3 {}\n", vec3i{q.x, q.y, q.z});
            else
                print(fs, "4 {}\n", q);
        }
        for (auto& l : lines) print(fs, "2 {}\n", l);
    } else {
        // write vertex data
        for (auto i = 0; i < positions.size(); i++) {
            if (!empty(positions)) write_value(fs, positions[i]);
            if (!empty(normals)) write_value(fs, normals[i]);
            if (!empty(texturecoords))
                write_value(fs, (!flip_texcoord) ? texturecoords[i] :
                                                   vec2f{texturecoords[i].x,
                                                       1 - texturecoords[i].y});
            if (!empty(colors)) write_value(fs, colors[i]);
            if (!empty(radius)) write_value(fs, radius[i]);
        }

        // write face data
        for (auto& t : triangles) {
            auto n = (byte)3;
            write_value(fs, n);
            write_value(fs, t);
        }
        for (auto& q : quads) {
            if (q.z == q.w) {
                auto n = (byte)3;
                write_value(fs, n);
                write_value(fs, vec3i{q.x, q.y, q.z});
            } else {
                auto n = (byte)4;
                write_value(fs, n);
                write_value(fs, q);
            }
        }
        for (auto& l : lines) {
            auto n = (byte)2;
            write_value(fs, n);
            write_value(fs, l);
        }
    }

    // done
    return true;
}

// Load ply mesh
bool load_obj_mesh(const string& filename, vector<int>& points,
    vector<vec2i>& lines, vector<vec3i>& triangles, vector<vec4i>& quads,
    vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texturecoords, bool force_triangles, bool flip_texcoord) {
    // clear
    auto colors = vector<vec4f>{};
    auto radius = vector<float>{};
    reset_mesh_data(points, lines, triangles, quads, positions, normals,
        texturecoords, colors, radius);

    // obj vertices
    auto opos      = std::deque<vec3f>();
    auto onorm     = std::deque<vec3f>();
    auto otexcoord = std::deque<vec2f>();

    // vertex maps
    auto vertex_map = unordered_map<obj_vertex, int, obj_vertex_hash>();

    // Add  vertices to the current shape
    auto add_verts = [&](const vector<obj_vertex>& verts) {
        for (auto& vert : verts) {
            auto it = vertex_map.find(vert);
            if (it != vertex_map.end()) continue;
            auto nverts = (int)positions.size();
            vertex_map.insert(it, {vert, nverts});
            if (vert.position) positions.push_back(opos.at(vert.position - 1));
            if (vert.texturecoord)
                texturecoords.push_back(otexcoord.at(vert.texturecoord - 1));
            if (vert.normal) normals.push_back(onorm.at(vert.normal - 1));
        }
    };

    auto cb     = obj_callbacks();
    cb.vert     = [&](vec3f v) { opos.push_back(v); };
    cb.norm     = [&](vec3f v) { onorm.push_back(v); };
    cb.texcoord = [&](vec2f v) { otexcoord.push_back(v); };
    cb.face     = [&](const vector<obj_vertex>& verts) {
        add_verts(verts);
        if (verts.size() == 4) {
            quads.push_back({vertex_map.at(verts[0]), vertex_map.at(verts[1]),
                vertex_map.at(verts[2]), vertex_map.at(verts[3])});
        } else {
            for (auto i = 2; i < verts.size(); i++)
                triangles.push_back({vertex_map.at(verts[0]),
                    vertex_map.at(verts[i - 1]), vertex_map.at(verts[i])});
        }
    };
    cb.line = [&](const vector<obj_vertex>& verts) {
        add_verts(verts);
        for (auto i = 1; i < verts.size(); i++)
            lines.push_back(
                {vertex_map.at(verts[i - 1]), vertex_map.at(verts[i])});
    };
    cb.point = [&](const vector<obj_vertex>& verts) {
        add_verts(verts);
        for (auto i = 0; i < verts.size(); i++)
            points.push_back(vertex_map.at(verts[i]));
    };

    // load obj
    auto obj_options          = load_obj_options();
    obj_options.exit_on_error = false;
    obj_options.geometry_only = true;
    obj_options.flip_texcoord = flip_texcoord;
    if (!load_obj(filename, cb, obj_options)) return false;

    // merging quads and triangles
    merge_triangles_and_quads(triangles, quads, force_triangles);

    // done
    return true;
}

// Load ply mesh
bool save_obj_mesh(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    bool flip_texcoord) {
    auto fs = open(filename, "wt");
    if (!fs) return false;

    for (auto& p : positions) print(fs, "v {}\n", p);
    for (auto& n : normals) print(fs, "vn {}\n", n);
    for (auto& t : texturecoords)
        print(fs, "vt {}\n", vec2f{t.x, (flip_texcoord) ? 1 - t.y : t.y});

    auto mask = obj_vertex{
        1, empty(texturecoords) ? 0 : 1, empty(normals) ? 0 : 1};
    auto vert = [mask](int i) {
        return obj_vertex{(i + 1) * mask.position, (i + 1) * mask.texturecoord,
            (i + 1) * mask.normal};
    };

    for (auto& p : points) {
        print(fs, "p {}\n", to_string(vert(p)));
    }
    for (auto& l : lines) {
        print(fs, "l {} {}\n", to_string(vert(l.x)), to_string(vert(l.y)));
    }
    for (auto& t : triangles) {
        print(fs, "f {} {} {}\n", to_string(vert(t.x)), to_string(vert(t.y)),
            to_string(vert(t.z)));
    }
    for (auto& q : quads) {
        if (q.z == q.w) {
            print(fs, "f {} {} {}\n", to_string(vert(q.x)),
                to_string(vert(q.y)), to_string(vert(q.z)));
        } else {
            print(fs, "f {} {} {} {}\n", to_string(vert(q.x)),
                to_string(vert(q.y)), to_string(vert(q.z)), to_string(vert(q.w)));
        }
    }

    return true;
}

// Reset mesh data
void reset_facevarying_mesh_data(vector<vec4i>& quads_positions,
    vector<vec4i>& quads_normals, vector<vec4i>& quads_textuercoords,
    vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texturecoords, vector<int>& quads_materials) {
    quads_positions     = {};
    quads_normals       = {};
    quads_textuercoords = {};
    positions           = {};
    normals             = {};
    texturecoords       = {};
    quads_materials     = {};
}

// Load ply mesh
bool load_facevarying_mesh(const string& filename, vector<vec4i>& quads_positions,
    vector<vec4i>& quads_normals, vector<vec4i>& quads_texturecoords,
    vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texturecoords, vector<int>& quads_materials) {
    auto ext = get_extension(filename);
    if (ext == "obj" || ext == "OBJ") {
        return load_obj_facevarying_mesh(filename, quads_positions,
            quads_normals, quads_texturecoords, positions, normals,
            texturecoords, quads_materials);
    } else {
        reset_facevarying_mesh_data(quads_positions, quads_normals,
            quads_texturecoords, positions, normals, texturecoords,
            quads_materials);
        return false;
    }
}

// Save ply mesh
bool save_facevarying_mesh(const string& filename,
    const vector<vec4i>& quads_positions, const vector<vec4i>& quads_normals,
    const vector<vec4i>& quads_texturecoords, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    const vector<int>& quads_materials, bool ascii) {
    auto ext = get_extension(filename);
    if (ext == "obj" || ext == "OBJ") {
        return save_obj_facevarying_mesh(filename, quads_positions,
            quads_normals, quads_texturecoords, positions, normals,
            texturecoords, quads_materials);
    } else {
        return false;
    }
}

// Load ply mesh
bool load_obj_facevarying_mesh(const string& filename,
    vector<vec4i>& quads_positions, vector<vec4i>& quads_normals,
    vector<vec4i>& quads_texturecoords, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texturecoords,
    vector<int>& quads_materials, bool flip_texcoord) {
    // clear
    reset_facevarying_mesh_data(quads_positions, quads_normals,
        quads_texturecoords, positions, normals, texturecoords, quads_materials);

    // obj vertices
    auto opos      = std::deque<vec3f>();
    auto onorm     = std::deque<vec3f>();
    auto otexcoord = std::deque<vec2f>();

    // vertex maps
    auto pos_map      = unordered_map<int, int>();
    auto texcoord_map = unordered_map<int, int>();
    auto norm_map     = unordered_map<int, int>();

    // material group
    auto material_group      = vector<string>();
    auto current_material_id = -1;

    // add vertex
    auto add_fvverts = [&](const vector<obj_vertex>& verts) {
        for (auto& vert : verts) {
            if (!vert.position) continue;
            auto pos_it = pos_map.find(vert.position);
            if (pos_it != pos_map.end()) continue;
            auto nverts = (int)positions.size();
            pos_map.insert(pos_it, {vert.position, nverts});
            positions.push_back(opos.at(vert.position - 1));
        }
        for (auto& vert : verts) {
            if (!vert.texturecoord) continue;
            auto texcoord_it = texcoord_map.find(vert.texturecoord);
            if (texcoord_it != texcoord_map.end()) continue;
            auto nverts = (int)texturecoords.size();
            texcoord_map.insert(texcoord_it, {vert.texturecoord, nverts});
            texturecoords.push_back(otexcoord.at(vert.texturecoord - 1));
        }
        for (auto& vert : verts) {
            if (!vert.normal) continue;
            auto norm_it = norm_map.find(vert.normal);
            if (norm_it != norm_map.end()) continue;
            auto nverts = (int)normals.size();
            norm_map.insert(norm_it, {vert.normal, nverts});
            normals.push_back(onorm.at(vert.normal - 1));
        }
    };

    auto cb     = obj_callbacks();
    cb.vert     = [&](vec3f v) { opos.push_back(v); };
    cb.norm     = [&](vec3f v) { onorm.push_back(v); };
    cb.texcoord = [&](vec2f v) { otexcoord.push_back(v); };
    cb.face     = [&](const vector<obj_vertex>& verts) {
        add_fvverts(verts);
        if (verts.size() == 4) {
            if (verts[0].position) {
                quads_positions.push_back({pos_map.at(verts[0].position),
                    pos_map.at(verts[1].position), pos_map.at(verts[2].position),
                    pos_map.at(verts[3].position)});
            }
            if (verts[0].texturecoord) {
                quads_texturecoords.push_back(
                    {texcoord_map.at(verts[0].texturecoord),
                        texcoord_map.at(verts[1].texturecoord),
                        texcoord_map.at(verts[2].texturecoord),
                        texcoord_map.at(verts[3].texturecoord)});
            }
            if (verts[0].normal) {
                quads_normals.push_back({norm_map.at(verts[0].normal),
                    norm_map.at(verts[1].normal), norm_map.at(verts[2].normal),
                    norm_map.at(verts[3].normal)});
            }
            quads_materials.push_back(current_material_id);
        } else {
            if (verts[0].position) {
                for (auto i = 2; i < verts.size(); i++)
                    quads_positions.push_back({pos_map.at(verts[0].position),
                        pos_map.at(verts[1].position),
                        pos_map.at(verts[i].position),
                        pos_map.at(verts[i].position)});
            }
            if (verts[0].texturecoord) {
                for (auto i = 2; i < verts.size(); i++)
                    quads_texturecoords.push_back(
                        {texcoord_map.at(verts[0].texturecoord),
                            texcoord_map.at(verts[1].texturecoord),
                            texcoord_map.at(verts[i].texturecoord),
                            texcoord_map.at(verts[i].texturecoord)});
            }
            if (verts[0].normal) {
                for (auto i = 2; i < verts.size(); i++)
                    quads_normals.push_back({norm_map.at(verts[0].normal),
                        norm_map.at(verts[1].normal), norm_map.at(verts[i].normal),
                        norm_map.at(verts[i].normal)});
            }
            for (auto i = 2; i < verts.size(); i++)
                quads_materials.push_back(current_material_id);
        }
    };
    cb.line = [&](const vector<obj_vertex>& verts) {
        log_error("lines not supported!");
    };
    cb.point = [&](const vector<obj_vertex>& verts) {
        log_error("points not supported!");
    };
    cb.usemtl = [&](const string& name) {
        auto pos = std::find(material_group.begin(), material_group.end(), name);
        if (pos == material_group.end()) {
            material_group.push_back(name);
            current_material_id = (int)material_group.size() - 1;
        } else {
            current_material_id = (int)(pos - material_group.begin());
        }
    };

    // load obj
    auto obj_options          = load_obj_options();
    obj_options.exit_on_error = false;
    obj_options.geometry_only = true;
    obj_options.flip_texcoord = flip_texcoord;
    if (!load_obj(filename, cb, obj_options)) return false;

    // cleanup materials ids
    if (std::all_of(quads_materials.begin(), quads_materials.end(),
            [b = quads_materials.front()](auto a) { return a == b; })) {
        quads_materials.clear();
    }

    // done
    return true;
}

// Load ply mesh
bool save_obj_facevarying_mesh(const string& filename,
    const vector<vec4i>& quads_positions, const vector<vec4i>& quads_normals,
    const vector<vec4i>& quads_texturecoords, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    const vector<int>& quads_materials, bool flip_texcoord) {
    auto fs = open(filename, "wt");
    if (!fs) return false;

    for (auto& p : positions) print(fs, "v {}\n", p);
    for (auto& n : normals) print(fs, "vn {}\n", n);
    for (auto& t : texturecoords)
        print(fs, "vt {}\n", vec2f{t.x, (flip_texcoord) ? 1 - t.y : t.y});

    auto fvmask = obj_vertex{
        1, empty(texturecoords) ? 0 : 1, empty(normals) ? 0 : 1};
    auto fvvert = [fvmask](int pi, int ti, int ni) {
        return obj_vertex{(pi + 1) * fvmask.position,
            (ti + 1) * fvmask.texturecoord, (ni + 1) * fvmask.normal};
    };
    auto last_material_id = -1;
    for (auto i = 0; i < quads_positions.size(); i++) {
        if (!empty(quads_materials) && quads_materials[i] != last_material_id) {
            last_material_id = quads_materials[i];
            print(fs, "usemtl material_{}\n", last_material_id);
        }
        auto qp = quads_positions.at(i);
        auto qt = !empty(quads_texturecoords) ? quads_texturecoords.at(i) :
                                                vec4i{-1, -1, -1, -1};
        auto qn = !empty(quads_normals) ? quads_normals.at(i) :
                                          vec4i{-1, -1, -1, -1};
        if (qp.z != qp.w) {
            print(fs, "f {} {} {} {}\n", to_string(fvvert(qp.x, qt.x, qn.x)),
                to_string(fvvert(qp.y, qt.y, qn.y)),
                to_string(fvvert(qp.z, qt.z, qn.z)),
                to_string(fvvert(qp.w, qt.w, qn.w)));
        } else {
            print(fs, "f {} {} {}\n", to_string(fvvert(qp.x, qt.x, qn.x)),
                to_string(fvvert(qp.y, qt.y, qn.y)),
                to_string(fvvert(qp.z, qt.z, qn.z)));
        }
    }

    return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF OBJ IO
// -----------------------------------------------------------------------------
namespace yocto {}

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF PLY IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load ply mesh
bool load_ply(const string& filename, ply_data& ply) {
    // open file
    ply     = {};
    auto fs = open(filename, "rb");
    if (!fs) return false;

    // parse header
    ply        = ply_data{};
    auto ascii = false;
    auto line  = ""s;
    while (read_line(fs, line)) {
        auto view = string_view{line};
        auto cmd  = ""s;
        parse_value(view, cmd);
        if (cmd == "") continue;
        if (cmd == "ply") {
        } else if (cmd == "comment") {
        } else if (cmd == "format") {
            auto fmt = ""s;
            parse_value(view, fmt);
            if (fmt != "ascii" && fmt != "binary_little_endian") return false;
            ascii = fmt == "ascii";
        } else if (cmd == "element") {
            auto elem = ply_element();
            parse_value(view, elem.name);
            parse_value(view, elem.count);
            ply.elements.push_back(elem);
        } else if (cmd == "property") {
            auto prop = ply_property();
            auto type = ""s;
            parse_value(view, type);
            if (type == "list") {
                auto count_type = ""s;
                parse_value(view, count_type);
                auto elem_type = ""s;
                parse_value(view, elem_type);
                if (count_type != "uchar" && count_type != "uint8")
                    log_error("unsupported ply list type");
                if (elem_type != "int" && elem_type != "uint")
                    log_error("unsupported ply list type");
                prop.type = ply_type::ply_int_list;
            } else if (type == "float") {
                prop.type = ply_type::ply_float;
            } else if (type == "uchar" || type == "uint8") {
                prop.type = ply_type::ply_uchar;
            } else if (type == "int") {
                prop.type = ply_type::ply_int;
            } else {
                return false;
            }
            parse_value(view, prop.name);
            prop.scalars.resize(ply.elements.back().count);
            if (prop.type == ply_type::ply_int_list)
                prop.lists.resize(ply.elements.back().count);
            ply.elements.back().properties.push_back(prop);
        } else if (cmd == "end_header") {
            break;
        } else {
            return false;
        }
    }

    // parse content
    if (ascii) {
        for (auto& elem : ply.elements) {
            for (auto vid = 0; vid < elem.count; vid++) {
                if (!read_line(fs, line)) return false;
                auto view = string_view{line};
                for (auto pid = 0; pid < elem.properties.size(); pid++) {
                    auto& prop = elem.properties[pid];
                    if (prop.type == ply_type::ply_float) {
                        auto v = 0.0f;
                        parse_value(view, v);
                        prop.scalars[vid] = v;
                    } else if (prop.type == ply_type::ply_int) {
                        auto v = 0;
                        parse_value(view, v);
                        prop.scalars[vid] = v;
                    } else if (prop.type == ply_type::ply_uchar) {
                        auto vc = (unsigned char)0;
                        auto v  = 0;
                        parse_value(view, v);
                        vc                = (unsigned char)v;
                        prop.scalars[vid] = vc / 255.0f;
                    } else if (prop.type == ply_type::ply_int_list) {
                        auto vc = (unsigned char)0;
                        auto v  = 0;
                        parse_value(view, v);
                        vc                = (unsigned char)v;
                        prop.scalars[vid] = vc;
                        for (auto i = 0; i < (int)prop.scalars[vid]; i++) {
                            auto v = 0;
                            parse_value(view, v);
                            prop.lists[vid][i] = v;
                        }
                    } else {
                        return false;
                    }
                }
            }
        }
    } else {
        for (auto& elem : ply.elements) {
            for (auto vid = 0; vid < elem.count; vid++) {
                for (auto pid = 0; pid < elem.properties.size(); pid++) {
                    auto& prop = elem.properties[pid];
                    if (prop.type == ply_type::ply_float) {
                        auto v = 0.0f;
                        if (!read_value(fs, v)) return false;
                        prop.scalars[vid] = v;
                    } else if (prop.type == ply_type::ply_int) {
                        auto v = 0;
                        if (!read_value(fs, v)) return false;
                        prop.scalars[vid] = v;
                    } else if (prop.type == ply_type::ply_uchar) {
                        auto vc = (unsigned char)0;
                        if (!read_value(fs, vc)) return false;
                        prop.scalars[vid] = vc / 255.0f;
                    } else if (prop.type == ply_type::ply_int_list) {
                        auto vc = (unsigned char)0;
                        if (!read_value(fs, vc)) return false;
                        prop.scalars[vid] = vc;
                        for (auto i = 0; i < (int)prop.scalars[vid]; i++) {
                            if (!read_value(fs, prop.lists[vid][i]))
                                return false;
                        }
                    } else {
                        return false;
                    }
                }
            }
        }
    }

    return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE IO
// -----------------------------------------------------------------------------
namespace yocto {

#if 0

// Load obj mesh
bool load_obj_fvmesh(const string& filename, vector<vec4i>& quads_positions,
    vector<vec3f>& pos, vector<vec4i>& quads_normals, vector<vec3f>& norm,
    vector<vec4i>& quads_texturecoords, vector<vec2f>& texcoord,
    bool flip_texcoord) {
    // clear
    vector<vec4i> quads_colors;
    vector<vec4f> color;
    reset_fvmesh_data(quads_positions, pos, quads_normals, norm,
        quads_texturecoords, texcoord, quads_colors, color);

    // obj vertex
    auto opos      = std::deque<vec3f>();
    auto onorm     = std::deque<vec3f>();
    auto otexcoord = std::deque<vec2f>();

    // vertex maps
    auto pos_map      = unordered_map<int, int>();
    auto texcoord_map = unordered_map<int, int>();
    auto norm_map     = unordered_map<int, int>();

    // add vertex
    auto add_verts = [&](const vector<obj_vertex>& verts) {
        for (auto& vert : verts) {
            if (!vert.position) continue;
            auto pos_it = pos_map.find(vert.position);
            if (pos_it != pos_map.end()) continue;
            auto nverts = (int)pos.size();
            pos_map.insert(pos_it, {vert.position, nverts});
            pos.push_back(opos.at(vert.position - 1));
        }
        for (auto& vert : verts) {
            if (!vert.texturecoord) continue;
            auto texcoord_it = texcoord_map.find(vert.texturecoord);
            if (texcoord_it != texcoord_map.end()) continue;
            auto nverts = (int)texcoord.size();
            texcoord_map.insert(texcoord_it, {vert.texturecoord, nverts});
            texcoord.push_back(otexcoord.at(vert.texturecoord - 1));
        }
        for (auto& vert : verts) {
            if (!vert.normal) continue;
            auto norm_it = norm_map.find(vert.normal);
            if (norm_it != norm_map.end()) continue;
            auto nverts = (int)norm.size();
            norm_map.insert(norm_it, {vert.normal, nverts});
            norm.push_back(onorm.at(vert.normal - 1));
        }
    };

    auto cb     = obj_callbacks();
    cb.vert     = [&](vec3f v) { opos.push_back(v); };
    cb.norm     = [&](vec3f v) { onorm.push_back(v); };
    cb.texcoord = [&](vec2f v) { otexcoord.push_back(v); };
    cb.face     = [&](const vector<obj_vertex>& verts) {
        add_verts(verts);
        if (verts.size() == 4) {
            if (verts[0].position) {
                quads_positions.push_back({pos_map.at(verts[0].position),
                    pos_map.at(verts[1].position), pos_map.at(verts[2].position),
                    pos_map.at(verts[3].position)});
            }
            if (verts[0].texturecoord) {
                quads_texturecoords.push_back(
                    {texcoord_map.at(verts[0].texturecoord),
                        texcoord_map.at(verts[1].texturecoord),
                        texcoord_map.at(verts[2].texturecoord),
                        texcoord_map.at(verts[3].texturecoord)});
            }
            if (verts[0].normal) {
                quads_normals.push_back({norm_map.at(verts[0].normal),
                    norm_map.at(verts[1].normal), norm_map.at(verts[2].normal),
                    norm_map.at(verts[3].normal)});
            }
        } else {
            if (verts[0].position) {
                for (auto i = 2; i < verts.size(); i++)
                    quads_positions.push_back({pos_map.at(verts[0].position),
                        pos_map.at(verts[1].position),
                        pos_map.at(verts[i].position),
                        pos_map.at(verts[i].position)});
            }
            if (verts[0].texturecoord) {
                for (auto i = 2; i < verts.size(); i++)
                    quads_texturecoords.push_back(
                        {texcoord_map.at(verts[0].texturecoord),
                            texcoord_map.at(verts[1].texturecoord),
                            texcoord_map.at(verts[i].texturecoord),
                            texcoord_map.at(verts[i].texturecoord)});
            }
            if (verts[0].normal) {
                for (auto i = 2; i < verts.size(); i++)
                    quads_normals.push_back({norm_map.at(verts[0].normal),
                        norm_map.at(verts[1].normal), norm_map.at(verts[i].normal),
                        norm_map.at(verts[i].normal)});
            }
        }
    };

    // load obj
    return load_obj(filename, cb, true, true, flip_texcoord);
}

// Load ply mesh
bool save_obj_fvmesh(const string& filename, const vector<vec4i>& quads_positions,
    const vector<vec3f>& positions, const vector<vec4i>& quads_normals,
    const vector<vec3f>& normals, const vector<vec4i>& quads_texturecoords,
    const vector<vec2f>& texturecoords, bool flip_texcoord) {
    auto fs = open(filename, "wt");
    if (!fs) return false;

    for (auto& p : positions) print(fs, "v {}\n", p);
    for (auto& n : normals) print(fs, "vn {}\n", n);
    for (auto& t : texturecoords)
        print(fs, "vt {}\n", vec2f{t.x, (flip_texcoord) ? 1 - t.y : t.y});
    auto mask = obj_vertex{
        1, empty(texturecoords) ? 0 : 1, empty(normals) ? 0 : 1};
    auto vert = [mask](int pif, int ti, int ni) {
        return obj_vertex{(pif + 1) * mask.position,
            (ti + 1) * mask.texturecoord, (ni + 1) * mask.normal};
    };
    for (auto i = 0; i < quads_positions.size(); i++) {
        auto qp = quads_positions.at(i);
        auto qt = !empty(quads_texturecoords) ? quads_texturecoords.at(i) :
                                                 vec4i{-1, -1, -1, -1};
        auto qn = !empty(quads_normals) ? quads_normals.at(i) :
                                           vec4i{-1, -1, -1, -1};
        if (qp.z != qp.w)
            print(fs, "f {} {} {} {}\n",
                to_string(vert(qp.x, qt.x, qn.x)).c_str(),
                to_string(vert(qp.y, qt.y, qn.y)).c_str(),
                to_string(vert(qp.z, qt.z, qn.z)).c_str(),
                to_string(vert(qp.w, qt.w, qn.w)).c_str());
        else
            print(fs, "f {} {} {}\n", to_string(vert(qp.x, qt.x, qn.x)).c_str(),
                to_string(vert(qp.y, qt.y, qn.y)).c_str(),
                to_string(vert(qp.z, qt.z, qn.z)).c_str());
    }

    return true;
}

#endif

}  // namespace yocto
