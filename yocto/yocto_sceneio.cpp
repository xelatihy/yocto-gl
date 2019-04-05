//
// Implementation for Yocto/GL Input and Output functions.
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

#include "yocto_sceneio.h"
#include "yocto_obj.h"
#include "yocto_pbrt.h"
#include "yocto_random.h"
#include "yocto_shape.h"
#include "yocto_utils.h"

#include "ext/happly.h"
#define CGLTF_IMPLEMENTATION
#include "ext/cgltf.h"
#include "ext/json.hpp"

#include <array>
#include <climits>
#include <cstdlib>
#include <memory>
#include <regex>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

using nlohmann::json;
using std::unique_ptr;

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON SUPPORT
// -----------------------------------------------------------------------------
namespace yocto {

// Load a JSON object
inline void load_json(const string& filename, json& js) {
    auto text = ""s;
    load_text(filename, text);
    js = json::parse(text);
}

// Save a JSON object
inline void save_json(const string& filename, const json& js) {
    // we have to use streams here since the json library is faster with them
    save_text(filename, js.dump(4));
}

template <typename T, int N>
inline void to_json(json& js, const vec<T, N>& val) {
    nlohmann::to_json(js, (const std::array<T, N>&)val);
}
template <typename T, int N>
inline void from_json(const json& js, vec<T, N>& val) {
    nlohmann::from_json(js, (std::array<T, N>&)val);
}

template <typename T, int N>
inline void to_json(json& js, const frame<T, N>& val) {
    nlohmann::to_json(js, (const std::array<T, N*(N + 1)>&)val);
}
template <typename T, int N>
inline void from_json(const json& js, frame<T, N>& val) {
    nlohmann::from_json(js, (std::array<T, N*(N + 1)>&)val);
}

template <typename T, int N, int M>
inline void to_json(json& js, const mat<T, N, M>& val) {
    nlohmann::to_json(js, (const std::array<T, N * M>&)val);
}
template <typename T, int N, int M>
inline void from_json(const json& js, mat<T, N, M>& val) {
    nlohmann::from_json(js, (std::array<T, N * M>&)val);
}

template <typename T, int N>
inline void to_json(json& js, const bbox<T, N>& val) {
    nlohmann::from_json(js, (std::array<T, N * 2>&)val);
}
template <typename T, int N>
inline void from_json(const json& js, bbox<T, N>& val) {
    nlohmann::to_json(js, (const std::array<T, N * 2>&)val);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF CONVERSION TO/FROM JSON
// -----------------------------------------------------------------------------
namespace yocto {

template <typename T>
inline void to_json(json& js, const image<T>& value) {
    js           = json::object();
    js["width"]  = value.size().x;
    js["height"] = value.size().y;
    js["pixels"] = value._pixels;
}
template <typename T>
inline void from_json(const json& js, image<T>& value) {
    auto width  = js.at("width").get<int>();
    auto height = js.at("height").get<int>();
    auto pixels = js.at("pixels").get<vector<T>>();
    value       = image{{width, height}, (const T*)pixels.data()};
}
template <typename T>
inline void to_json(json& js, const volume<T>& value) {
    js           = json::object();
    js["width"]  = value.size().x;
    js["height"] = value.size().y;
    js["depth"]  = value.size().z;
    js["voxels"] = value._voxels;
}
template <typename T>
inline void from_json(const json& js, volume<T>& value) {
    auto width  = js.at("width").get<int>();
    auto height = js.at("height").get<int>();
    auto depth  = js.at("depth").get<int>();
    auto voxels = js.at("voxels").get<vector<T>>();
    value       = volume{{width, height, depth}, (const T*)voxels.data()};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GENERIC SCENE LOADING
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a scene in the builtin YAML format.
void load_yaml_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options);
void save_yaml_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options);

// Load/save a scene from/to OBJ.
void load_obj_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options);
void save_obj_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options);

// Load/save a scene from/to PLY. Loads/saves only one mesh with no other data.
void load_ply_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options);
void save_ply_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options);

// Load/save a scene from/to glTF.
void load_gltf_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options);
void save_gltf_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options);

// Load/save a scene from/to pbrt. This is not robust at all and only
// works on scene that have been previously adapted since the two renderers
// are too different to match.
void load_pbrt_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options);
void save_pbrt_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options);

// Load/save a binary dump useful for very fast scene IO. This format is not
// an archival format and should only be used as an intermediate format.
void load_ybin_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options);
void save_ybin_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options);

// Load a scene
void load_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options) {
    auto ext = get_extension(filename);
    if (ext == "yaml" || ext == "YAML") {
        load_yaml_scene(filename, scene, options);
    } else if (ext == "obj" || ext == "OBJ") {
        load_obj_scene(filename, scene, options);
    } else if (ext == "gltf" || ext == "GLTF") {
        load_gltf_scene(filename, scene, options);
    } else if (ext == "pbrt" || ext == "PBRT") {
        load_pbrt_scene(filename, scene, options);
    } else if (ext == "ybin" || ext == "YBIN") {
        load_ybin_scene(filename, scene, options);
    } else if (ext == "ply" || ext == "PLY") {
        load_ply_scene(filename, scene, options);
    } else {
        scene = {};
        throw io_error("unsupported scene format " + ext);
    }
}

// Save a scene
void save_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options) {
    auto ext = get_extension(filename);
    if (ext == "yaml" || ext == "YAML") {
        save_yaml_scene(filename, scene, options);
    } else if (ext == "obj" || ext == "OBJ") {
        save_obj_scene(filename, scene, options);
    } else if (ext == "yaml" || ext == "YAML") {
        save_yaml_scene(filename, scene, options);
    } else if (ext == "gltf" || ext == "GLTF") {
        save_gltf_scene(filename, scene, options);
    } else if (ext == "pbrt" || ext == "PBRT") {
        save_pbrt_scene(filename, scene, options);
    } else if (ext == "ybin" || ext == "YBIN") {
        save_ybin_scene(filename, scene, options);
    } else {
        throw io_error("unsupported scene format " + ext);
    }
}

void load_scene_texture(yocto_texture& texture, const string& dirname) {
    if (is_image_preset_filename(texture.uri)) {
        auto [type, nfilename] = get_image_preset_type(texture.uri);
        make_image_preset(texture.hdr_image, texture.ldr_image, type);
        texture.uri = nfilename;
    } else {
        load_image(dirname + texture.uri, texture.hdr_image, texture.ldr_image);
    }
}

void load_scene_voltexture(yocto_voltexture& texture, const string& dirname) {
    if (is_volume_preset_filename(texture.uri)) {
        make_volume_preset(texture.volume_data, get_basename(texture.uri));
        texture.uri = get_noextension(texture.uri) + ".yvol";
    } else {
        load_volume(dirname + texture.uri, texture.volume_data);
    }
}

void load_scene_textures(yocto_scene& scene, const string& dirname,
    const load_scene_options& options) {
    if (options.skip_textures) return;

    // load images
    parallel_foreach(
        scene.textures,
        [&dirname](yocto_texture& texture) {
            if (!texture.hdr_image.empty() || !texture.ldr_image.empty())
                return;
            load_scene_texture(texture, dirname);
        },
        options.cancel_flag, options.run_serially);

    // load volumes
    parallel_foreach(
        scene.voltextures,
        [&dirname](yocto_voltexture& texture) {
            if (!texture.volume_data.empty()) return;
            load_scene_voltexture(texture, dirname);
        },
        options.cancel_flag, options.run_serially);
}

void save_scene_texture(const yocto_texture& texture, const string& dirname) {
    save_image(dirname + texture.uri, texture.hdr_image, texture.ldr_image);
}

void save_scene_voltexture(
    const yocto_voltexture& texture, const string& dirname) {
    save_volume(dirname + texture.uri, texture.volume_data);
}

// helper to save textures
void save_scene_textures(const yocto_scene& scene, const string& dirname,
    const save_scene_options& options) {
    if (options.skip_textures) return;

    // save images
    parallel_foreach(
        scene.textures,
        [&dirname](const yocto_texture& texture) {
            save_scene_texture(texture, dirname);
        },
        options.cancel_flag, options.run_serially);

    // save volumes
    parallel_foreach(
        scene.voltextures,
        [&dirname](const yocto_voltexture& texture) {
            save_scene_voltexture(texture, dirname);
        },
        options.cancel_flag, options.run_serially);
}

void load_scene_shape(yocto_shape& shape, const string& dirname) {
    if (is_shape_preset_filename(shape.uri)) {
        auto [type, nfilename] = get_shape_preset_type(shape.uri);
        make_shape_preset(shape.points, shape.lines, shape.triangles,
            shape.quads, shape.quads_positions, shape.quads_normals,
            shape.quads_texturecoords, shape.positions, shape.normals,
            shape.texturecoords, shape.colors, shape.radius, type);
        shape.uri = nfilename;
    } else {
        load_shape(dirname + shape.uri, shape.points, shape.lines,
            shape.triangles, shape.quads, shape.quads_positions,
            shape.quads_normals, shape.quads_texturecoords, shape.positions,
            shape.normals, shape.texturecoords, shape.colors, shape.radius,
            shape.preserve_facevarying);
    }
}

void save_scene_shape(const yocto_shape& shape, const string& dirname) {
    save_shape(dirname + shape.uri, shape.points, shape.lines, shape.triangles,
        shape.quads, shape.quads_positions, shape.quads_normals,
        shape.quads_texturecoords, shape.positions, shape.normals,
        shape.texturecoords, shape.colors, shape.radius);
}

// Load json meshes
void load_scene_shapes(yocto_scene& scene, const string& dirname,
    const load_scene_options& options) {
    if (options.skip_meshes) return;

    // load shapes
    parallel_foreach(
        scene.shapes,
        [&dirname](yocto_shape& shape) { load_scene_shape(shape, dirname); },
        options.cancel_flag, options.run_serially);
}

// Save json meshes
void save_scene_shapes(const yocto_scene& scene, const string& dirname,
    const save_scene_options& options) {
    if (options.skip_meshes) return;

    // save shapes
    parallel_foreach(
        scene.shapes,
        [&dirname](
            const yocto_shape& shape) { save_scene_shape(shape, dirname); },
        options.cancel_flag, options.run_serially);
}

// check if it is really face varying
bool is_face_varying(const vector<vec4i>& quads_positions,
    const vector<vec4i>& quads_normals, const vector<vec4i>& quads_texcoords) {
    if (quads_positions.empty()) return false;
    if (!quads_normals.empty()) {
        for (auto i = 0; i < quads_positions.size(); i++)
            if (quads_positions[i] != quads_normals[i]) return true;
    }
    if (!quads_texcoords.empty()) {
        for (auto i = 0; i < quads_positions.size(); i++)
            if (quads_positions[i] != quads_texcoords[i]) return true;
    }
    return false;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// YAML SUPPORT
// -----------------------------------------------------------------------------
namespace yocto {

struct yaml_value {
    enum struct type_t { none_t, string_t, boolean_t, number_t, array_t };
    type_t            type = type_t::none_t;
    string            _string;
    bool              _boolean = false;
    array<double, 16> _numbers;
    int               _size = 0;
};

inline void get_value(const yaml_value& yml, int& value) {
    if (yml.type != yaml_value::type_t::number_t)
        throw io_error("int expected");
    value = (int)yml._numbers[0];
}

inline void get_value(const yaml_value& yml, float& value) {
    if (yml.type != yaml_value::type_t::number_t)
        throw io_error("number expected");
    value = (float)yml._numbers[0];
}

inline void get_value(const yaml_value& yml, bool& value) {
    if (yml.type != yaml_value::type_t::boolean_t)
        throw io_error("string expected");
    value = yml._boolean;
}

inline void get_value(const yaml_value& yml, string& value) {
    if (yml.type != yaml_value::type_t::string_t)
        throw io_error("string expected");
    value = yml._string;
}

template <typename T, size_t N>
inline void get_value(const yaml_value& yml, array<T, N>& value) {
    if (yml.type != yaml_value::type_t::array_t || yml._size != N)
        throw io_error("array expected");
    for (auto i = 0; i < N; i++) value[i] = (T)yml._numbers[i];
}

template <typename T, int N>
inline void get_value(const yaml_value& yml, vec<T, N>& value) {
    get_value(yml, (array<T, N>&)value);
}
template <typename T, int N>
inline void get_value(const yaml_value& yml, frame<T, N>& value) {
    get_value(yml, (array<T, N*(N + 1)>&)value);
}
template <typename T, int N, int M>
inline void get_value(const yaml_value& yml, mat<T, N, M>& value) {
    get_value(yml, (array<T, N * M>&)value);
}

inline bool is_string_value(const yaml_value& yml) {
    return yml.type == yaml_value::type_t::string_t;
}

struct yaml_callbacks {
    void object_group(const string& key);
    void object_begin();
    void key_value(const string& key, const yaml_value& value);
};

inline void parse_value(string_view& str, yaml_value& value) {
    skip_whitespace(str);
    if (str.empty()) throw io_error("cannot parse value");
    value = {};
    if (str.front() == '[') {
        str.remove_prefix(1);
        value._size = 0;
        while (true) {
            skip_whitespace(str);
            parse_value(str, value._numbers[value._size++]);
            skip_whitespace(str);
            if (str.empty()) throw io_error("cannot parse value");
            if (str.front() == ']') {
                str.remove_prefix(1);
                break;
            } else if (str.front() == ',') {
                str.remove_prefix(1);
            } else {
                throw io_error("cannot parse value");
            }
            if (value._size >= 16) throw io_error("cannot parse value");
        }
        value.type = yaml_value::type_t::array_t;
    } else if (str.front() == '-' || str.front() == '+' ||
               is_digit(str.front())) {
        parse_value(str, value._numbers[0]);
        value._size = 1;
        value.type  = yaml_value::type_t::number_t;
    } else if (str.front() == '"') {
        parse_value(str, value._string, true);
        value.type = yaml_value::type_t::string_t;
    } else {
        parse_value(str, value._string, false);
        value.type = yaml_value::type_t::string_t;
        if (value._string == "true" || value._string == "True") {
            value._string  = {};
            value._boolean = true;
            value.type     = yaml_value::type_t::boolean_t;
        } else if (value._string == "false" || value._string == "False") {
            value._string  = {};
            value._boolean = false;
            value.type     = yaml_value::type_t::boolean_t;
        }
    }
}

struct load_yaml_options {};

template <typename Callbacks>
inline void load_yaml(const string& filename, Callbacks& callbacks,
    const load_yaml_options& options) {
    // open file
    auto fs = input_file(filename);

    // parsing state
    auto in_objects = false;
    auto in_object  = false;

    // read the file line by line
    char buffer[4096];
    auto value = yaml_value{};
    while (read_line(fs, buffer, sizeof(buffer))) {
        // line
        auto line = string_view{buffer};
        remove_comment_and_newline(line);
        if (line.empty()) continue;
        if (is_whitespace(line)) continue;

        // peek commands
        if (is_space(line.front())) {
            if (!in_objects) throw io_error("bad yaml");
            // indented property
            skip_whitespace(line);
            if (line.empty()) throw io_error("bad yaml");
            if (line.front() == '-') {
                callbacks.object_begin();
                line.remove_prefix(1);
                skip_whitespace(line);
                in_object = true;
            }
            auto key = ""s;
            parse_varname(line, key);
            skip_whitespace(line);
            if (line.empty() || line.front() != ':') throw io_error("bad yaml");
            line.remove_prefix(1);
            skip_whitespace(line);
            parse_value(line, value);
            callbacks.key_value(key, value);
        } else if (is_alpha(line.front())) {
            // new group
            auto key = ""s;
            parse_varname(line, key);
            skip_whitespace(line);
            if (line.empty() || line.front() != ':') throw io_error("bad yaml");
            line.remove_prefix(1);
            if (!line.empty() && !is_whitespace(line))
                throw io_error("bad yaml");
            callbacks.object_group(key);
            in_objects = true;
            in_object  = false;
        } else {
            throw io_error("bad yaml");
        }
    }
}

// Save a scene in the builtin YAML format.
void load_yaml_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options) {
    scene = {};

    struct parse_callbacks : yaml_callbacks {
        yocto_scene&              scene;
        const load_scene_options& options;

        enum struct parsing_type {
            none,
            camera,
            texture,
            voltexture,
            material,
            shape,
            instance,
            environment
        };
        parsing_type type = parsing_type::none;

        unordered_map<string, int> tmap = {{"", -1}};
        unordered_map<string, int> vmap = {{"", -1}};
        unordered_map<string, int> mmap = {{"", -1}};
        unordered_map<string, int> smap = {{"", -1}};

        parse_callbacks(yocto_scene& scene, const load_scene_options& options)
            : scene{scene}, options{options} {}

        void get_ref(const yaml_value& yml, int& value,
            const unordered_map<string, int>& refs) const {
            if (is_string_value(yml)) {
                auto name = ""s;
                get_value(yml, name);
                if (name == "") return;
                try {
                    value = refs.at(name);
                } catch (...) {
                    throw io_error("reference not found " + name);
                }
            } else {
                get_value(yml, value);
                if (value < 0) return;
                if (value >= refs.size()) {
                    throw io_error("reference not found " + to_string(value));
                }
            }
        }

        void object_group(const string& key) {
            if (key == "cameras") {
                type = parsing_type::camera;
            } else if (key == "textures") {
                type = parsing_type::texture;
            } else if (key == "voltextures") {
                type = parsing_type::voltexture;
            } else if (key == "materials") {
                type = parsing_type::material;
            } else if (key == "shapes") {
                type = parsing_type::shape;
            } else if (key == "instances") {
                type = parsing_type::instance;
            } else if (key == "environments") {
                type = parsing_type::environment;
            } else {
                type = parsing_type::none;
                throw io_error("unknown object type");
            }
        }
        void object_begin() {
            switch (type) {
                case parsing_type::camera: scene.cameras.push_back({}); break;
                case parsing_type::texture: scene.textures.push_back({}); break;
                case parsing_type::voltexture:
                    scene.voltextures.push_back({});
                    break;
                case parsing_type::material:
                    scene.materials.push_back({});
                    break;
                case parsing_type::shape: scene.shapes.push_back({}); break;
                case parsing_type::instance:
                    scene.instances.push_back({});
                    break;
                case parsing_type::environment:
                    scene.environments.push_back({});
                    break;
                default: throw io_error("unknown object type");
            }
        }
        void key_value(const string& key, const yaml_value& value) {
            switch (type) {
                case parsing_type::camera: {
                    auto& camera = scene.cameras.back();
                    if (key == "uri") {
                        get_value(value, camera.uri);
                    } else if (key == "frame") {
                        get_value(value, camera.frame);
                    } else if (key == "orthographic") {
                        get_value(value, camera.orthographic);
                    } else if (key == "film_width") {
                        get_value(value, camera.film_width);
                    } else if (key == "film_height") {
                        get_value(value, camera.film_height);
                    } else if (key == "focal_length") {
                        get_value(value, camera.focal_length);
                    } else if (key == "focus_distance") {
                        get_value(value, camera.focus_distance);
                    } else if (key == "lens_aperture") {
                        get_value(value, camera.lens_aperture);
                    } else {
                        throw io_error("unknown property");
                    }
                } break;
                case parsing_type::texture: {
                    auto& texture = scene.textures.back();
                    if (key == "uri") {
                        get_value(value, texture.uri);
                        auto refname = texture.uri;
                        if (is_image_preset_filename(refname)) {
                            auto [_, nname] = get_image_preset_type(refname);
                            refname         = nname;
                        }
                        tmap[refname] = (int)scene.textures.size() - 1;
                    } else if (key == "filename") {
                        get_value(value, texture.uri);
                    } else {
                        throw io_error("unknown property");
                    }
                } break;
                case parsing_type::voltexture: {
                    auto& texture = scene.voltextures.back();
                    if (key == "uri") {
                        get_value(value, texture.uri);
                        auto refname = texture.uri;
                        if (is_image_preset_filename(refname)) {
                            auto [_, nname] = get_image_preset_type(refname);
                            refname         = nname;
                        }
                        vmap[refname] = (int)scene.voltextures.size() - 1;
                    } else {
                        throw io_error("unknown property");
                    }
                } break;
                case parsing_type::material: {
                    auto& material = scene.materials.back();
                    if (key == "uri") {
                        get_value(value, material.uri);
                        mmap[material.uri] = (int)scene.materials.size() - 1;
                    } else if (key == "base_metallic") {
                        get_value(value, material.base_metallic);
                    } else if (key == "gltf_textures") {
                        get_value(value, material.gltf_textures);
                    } else if (key == "emission") {
                        get_value(value, material.emission);
                    } else if (key == "diffuse") {
                        get_value(value, material.diffuse);
                    } else if (key == "specular") {
                        get_value(value, material.specular);
                    } else if (key == "transmission") {
                        get_value(value, material.transmission);
                    } else if (key == "roughness") {
                        get_value(value, material.roughness);
                    } else if (key == "opacity") {
                        get_value(value, material.opacity);
                    } else if (key == "fresnel") {
                        get_value(value, material.fresnel);
                    } else if (key == "refract") {
                        get_value(value, material.refract);
                    } else if (key == "volume_density") {
                        get_value(value, material.volume_density);
                    } else if (key == "volume_albedo") {
                        get_value(value, material.volume_albedo);
                    } else if (key == "volume_phaseg") {
                        get_value(value, material.volume_phaseg);
                    } else if (key == "emission_texture") {
                        get_ref(value, material.emission_texture, tmap);
                    } else if (key == "diffuse_texture") {
                        get_ref(value, material.diffuse_texture, tmap);
                    } else if (key == "specular_texture") {
                        get_ref(value, material.specular_texture, tmap);
                    } else if (key == "transmission_texture") {
                        get_ref(value, material.transmission_texture, tmap);
                    } else if (key == "roughness_texture") {
                        get_ref(value, material.roughness_texture, tmap);
                    } else if (key == "displacement_texture") {
                        get_ref(value, material.displacement_texture, tmap);
                    } else if (key == "normal_texture") {
                        get_ref(value, material.normal_texture, tmap);
                    } else if (key == "volume_density_texture") {
                        get_ref(value, material.volume_density_texture, vmap);
                    } else if (key == "displacement_scale") {
                        get_value(value, material.displacement_scale);
                    } else {
                        throw io_error("unknown property");
                    }
                } break;
                case parsing_type::shape: {
                    auto& shape = scene.shapes.back();
                    if (key == "uri") {
                        get_value(value, shape.uri);
                        auto refname = shape.uri;
                        if (is_shape_preset_filename(refname)) {
                            auto [_, nname] = get_shape_preset_type(refname);
                            refname         = nname;
                        }
                        smap[refname] = (int)scene.shapes.size() - 1;
                    } else if (key == "subdivision_level") {
                        get_value(value, shape.subdivision_level);
                    } else if (key == "catmull_clark") {
                        get_value(value, shape.catmull_clark);
                    } else if (key == "compute_normals") {
                        get_value(value, shape.compute_normals);
                    } else if (key == "preserve_facevarying") {
                        get_value(value, shape.preserve_facevarying);
                    } else {
                        throw io_error("unknown property");
                    }
                } break;
                case parsing_type::instance: {
                    auto& instance = scene.instances.back();
                    if (key == "uri") {
                        get_value(value, instance.uri);
                    } else if (key == "frame") {
                        get_value(value, instance.frame);
                    } else if (key == "shape") {
                        get_ref(value, instance.shape, smap);
                    } else if (key == "material") {
                        get_ref(value, instance.material, mmap);
                    } else {
                        throw io_error("unknown property");
                    }
                } break;
                case parsing_type::environment: {
                    auto& environment = scene.environments.back();
                    if (key == "uri") {
                        get_value(value, environment.uri);
                    } else if (key == "frame") {
                        get_value(value, environment.frame);
                    } else if (key == "emission") {
                        get_value(value, environment.emission);
                    } else if (key == "emission_texture") {
                        get_ref(value, environment.emission_texture, tmap);
                    } else {
                        throw io_error("unknown property");
                    }
                } break;
                default: throw io_error("unknown object type");
            }
        }
    };

    try {
        // Parse obj
        auto yaml_options = load_yaml_options();
        auto cb           = parse_callbacks{scene, options};
        load_yaml(filename, cb, yaml_options);

        // load shape and textures
        auto dirname = get_dirname(filename);
        load_scene_shapes(scene, dirname, options);
        load_scene_textures(scene, dirname, options);
    } catch (const std::exception& e) {
        throw io_error("cannot load scene " + filename + "\n" + e.what());
    }

    // fix scene
    scene.uri = get_filename(filename);
    add_missing_cameras(scene);
    add_missing_materials(scene);
    normalize_uris(scene);
    trim_memory(scene);
    update_transforms(scene);
}

inline void print_yaml_value(output_file& fs, int value) {
    print_value(fs, value);
}
inline void print_yaml_value(output_file& fs, float value) {
    print_value(fs, value);
}
inline void print_yaml_value(output_file& fs, const string& value) {
    print_value(fs, value);
}
inline void print_yaml_value(output_file& fs, const char* value) {
    print_value(fs, value);
}
inline void print_yaml_value(output_file& fs, bool value) {
    print_value(fs, value, true);
}
template <typename T, int N>
inline void print_yaml_value(output_file& fs, const vec<T, N>& value) {
    print_value(fs, value, true);
}
template <typename T, int N>
inline void print_yaml_value(output_file& fs, const frame<T, N>& value) {
    print_value(fs, value, true);
}

// Save yaml
void save_yaml(const string& filename, const yocto_scene& scene) {
    // open file
    auto fs = output_file(filename);

    print_value(
        fs, "# Written by Yocto/GL\n# https://github.com/xelatihy/yocto-gl\n");

    static const auto def_camera      = yocto_camera{};
    static const auto def_texture     = yocto_texture{};
    static const auto def_voltexture  = yocto_voltexture{};
    static const auto def_material    = yocto_material{};
    static const auto def_shape       = yocto_shape{};
    static const auto def_instance    = yocto_instance{};
    static const auto def_environment = yocto_environment{};

    auto print_first = [](output_file& fs, const char* name, auto& value) {
        print_value(fs, "  - ");
        print_value(fs, name);
        print_value(fs, ": ");
        print_yaml_value(fs, value);
        print_value(fs, "\n");
    };
    auto print_optional = [](output_file& fs, const char* name, auto& value,
                              auto& def) {
        if (value == def) return;
        print_value(fs, "    ");
        print_value(fs, name);
        print_value(fs, ": ");
        print_yaml_value(fs, value);
        print_value(fs, "\n");
    };
    auto print_ref = [as_int = false](output_file& fs, const char* name,
                         int value, auto& refs) {
        if (value < 0) return;
        print_value(fs, "    ");
        print_value(fs, name);
        print_value(fs, ": ");
        if (as_int) {
            print_yaml_value(fs, value);
        } else {
            print_yaml_value(fs, refs[value].uri);
        }
        print_value(fs, "\n");
    };

    if (!scene.cameras.empty()) print_value(fs, "\n\ncameras:\n");
    for (auto& camera : scene.cameras) {
        print_first(fs, "uri", camera.uri);
        print_optional(fs, "frame", camera.frame, def_camera.frame);
        print_optional(
            fs, "orthographic", camera.orthographic, def_camera.orthographic);
        print_optional(
            fs, "film_width", camera.film_width, def_camera.film_width);
        print_optional(
            fs, "film_height", camera.film_height, def_camera.film_height);
        print_optional(
            fs, "focal_length", camera.focal_length, def_camera.focal_length);
        print_optional(fs, "focus_distance", camera.focus_distance,
            def_camera.focus_distance);
        print_optional(fs, "lens_aperture", camera.lens_aperture,
            def_camera.lens_aperture);
    }

    if (!scene.textures.empty()) print_value(fs, "\n\ntextures:\n");
    for (auto& texture : scene.textures) {
        print_first(fs, "uri", texture.uri);
    }

    if (!scene.voltextures.empty()) print_value(fs, "\n\nvoltextures:\n");
    for (auto& texture : scene.voltextures) {
        print_first(fs, "uri", texture.uri);
    }

    if (!scene.materials.empty()) print_value(fs, "\n\nmaterials:\n");
    for (auto& material : scene.materials) {
        print_first(fs, "uri", material.uri);
        print_optional(fs, "base_metallic", material.base_metallic,
            def_material.base_metallic);
        print_optional(fs, "gltf_textures", material.gltf_textures,
            def_material.gltf_textures);
        print_optional(
            fs, "emission", material.emission, def_material.emission);
        print_optional(fs, "diffuse", material.diffuse, def_material.diffuse);
        print_optional(
            fs, "specular", material.specular, def_material.specular);
        print_optional(fs, "transmission", material.transmission,
            def_material.transmission);
        print_optional(
            fs, "roughness", material.roughness, def_material.roughness);
        print_optional(fs, "opacity", material.opacity, def_material.opacity);
        print_optional(fs, "fresnel", material.fresnel, def_material.fresnel);
        print_optional(fs, "refract", material.refract, def_material.refract);
        print_ref(
            fs, "emission_texture", material.emission_texture, scene.textures);
        print_ref(
            fs, "diffuse_texture", material.diffuse_texture, scene.textures);
        print_ref(
            fs, "specular_texture", material.specular_texture, scene.textures);
        print_ref(fs, "transmission_texture", material.transmission_texture,
            scene.textures);
        print_ref(fs, "roughness_texture", material.roughness_texture,
            scene.textures);
        print_ref(fs, "displacement_texture", material.displacement_texture,
            scene.textures);
        print_ref(
            fs, "normal_texture", material.normal_texture, scene.textures);
        print_optional(fs, "volume_emission", material.volume_emission,
            def_material.volume_emission);
        print_optional(fs, "volume_albedo", material.volume_albedo,
            def_material.volume_albedo);
        print_optional(fs, "volume_density", material.volume_density,
            def_material.volume_density);
        print_optional(fs, "volume_phaseg", material.volume_phaseg,
            def_material.volume_phaseg);
        print_ref(fs, "volume_density_texture", material.volume_density_texture,
            scene.voltextures);
        // TODO: add displacement scale
    }

    if (!scene.shapes.empty()) print_value(fs, "\n\nshapes:\n");
    for (auto& shape : scene.shapes) {
        print_first(fs, "uri", shape.uri);
        print_optional(fs, "subdivision_level", shape.subdivision_level,
            def_shape.subdivision_level);
        print_optional(
            fs, "catmull_clark", shape.catmull_clark, def_shape.catmull_clark);
        print_optional(fs, "compute_normals", shape.compute_normals,
            def_shape.compute_normals);
        print_optional(fs, "preserve_facevarying", shape.preserve_facevarying,
            def_shape.preserve_facevarying);
    }

    if (!scene.instances.empty()) print_value(fs, "\n\ninstances:\n");
    for (auto& instance : scene.instances) {
        print_first(fs, "uri", instance.uri);
        print_optional(fs, "frame", instance.frame, def_instance.frame);
        print_ref(fs, "shape", instance.shape, scene.shapes);
        print_ref(fs, "material", instance.material, scene.materials);
    }

    if (!scene.environments.empty()) print_value(fs, "\n\nenvironments:\n");
    for (auto& environment : scene.environments) {
        print_first(fs, "uri", environment.uri);
        print_optional(fs, "frame", environment.frame, def_environment.frame);
        print_optional(
            fs, "emission", environment.emission, def_environment.emission);
        print_ref(fs, "emission_texture", environment.emission_texture,
            scene.textures);
    }
}

// Save a scene in the builtin YAML format.
void save_yaml_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options) {
    try {
        // save yaml file
        save_yaml(filename, scene);

        // save meshes and textures
        auto dirname = get_dirname(filename);
        save_scene_shapes(scene, dirname, options);
        save_scene_textures(scene, dirname, options);
    } catch (const std::exception& e) {
        throw io_error("cannot load scene " + filename + "\n" + e.what());
    }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// OBJ CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

// Get material index
void set_obj_material(
    const string& name, int& material_id, const yocto_scene& scene) {
    for (material_id = 0; material_id < scene.materials.size(); material_id++) {
        if (scene.materials[material_id].uri == name) return;
    }
    throw io_error("unknown material " + name);
}

// Loads an OBJ
void load_obj_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options) {
    scene = {};

    struct parse_callbacks : obj_callbacks {
        yocto_scene&              scene;
        const load_scene_options& options;

        // current parsing values
        string mname = ""s;
        string oname = ""s;
        string gname = ""s;

        // vertices
        deque<vec3f> opos      = deque<vec3f>();
        deque<vec3f> onorm     = deque<vec3f>();
        deque<vec2f> otexcoord = deque<vec2f>();

        // object maps
        unordered_map<string, int> tmap = unordered_map<string, int>{{"", -1}};
        unordered_map<string, int> vmap = unordered_map<string, int>{{"", -1}};
        unordered_map<string, int> mmap = unordered_map<string, int>{{"", -1}};

        // vertex maps
        unordered_map<obj_vertex, int> vertex_map =
            unordered_map<obj_vertex, int>();
        unordered_map<int, int> pos_map      = unordered_map<int, int>();
        unordered_map<int, int> norm_map     = unordered_map<int, int>();
        unordered_map<int, int> texcoord_map = unordered_map<int, int>();

        parse_callbacks(yocto_scene& scene, const load_scene_options& options)
            : scene{scene}, options{options} {}

        // add object if needed
        void add_shape() {
            auto shape = yocto_shape{};
            shape.uri  = oname + gname;
            shape.preserve_facevarying =
                options.obj_preserve_face_varying ||
                shape.uri.find("[yocto::facevarying]") != string::npos;
            scene.shapes.push_back(shape);
            auto instance     = yocto_instance{};
            instance.uri      = shape.uri;
            instance.shape    = (int)scene.shapes.size() - 1;
            instance.material = mmap.at(mname);
            scene.instances.push_back(instance);
            vertex_map.clear();
            pos_map.clear();
            norm_map.clear();
            texcoord_map.clear();
        }
        // Parse texture options and name
        int add_texture(const obj_texture_info& info, bool force_linear) {
            if (info.path == "") return -1;
            if (tmap.find(info.path) != tmap.end()) {
                return tmap.at(info.path);
            }

            // create texture
            auto texture = yocto_texture{};
            texture.uri  = info.path;
            texture.uri  = info.path;
            scene.textures.push_back(texture);
            auto index      = (int)scene.textures.size() - 1;
            tmap[info.path] = index;

            return index;
        }
        // Parse texture options and name
        int add_voltexture(const obj_texture_info& info, bool srgb) {
            if (info.path == "") return -1;
            if (vmap.find(info.path) != vmap.end()) {
                return vmap.at(info.path);
            }

            // create texture
            auto texture = yocto_voltexture{};
            texture.uri  = info.path;
            scene.voltextures.push_back(texture);
            auto index      = (int)scene.voltextures.size() - 1;
            vmap[info.path] = index;

            return index;
        }
        // Add  vertices to the current shape
        void add_verts(const vector<obj_vertex>& verts, yocto_shape& shape) {
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
                if (vert.normal)
                    shape.normals.push_back(onorm.at(vert.normal - 1));
                if (shape.normals.size() != 0 &&
                    shape.normals.size() != shape.positions.size()) {
                    while (shape.normals.size() != shape.positions.size())
                        shape.normals.push_back({0, 0, 1});
                }
                if (shape.texturecoords.size() != 0 &&
                    shape.texturecoords.size() != shape.positions.size()) {
                    while (shape.texturecoords.size() != shape.positions.size())
                        shape.texturecoords.push_back({0, 0});
                }
            }
        }
        // add vertex
        void add_fvverts(const vector<obj_vertex>& verts, yocto_shape& shape) {
            for (auto& vert : verts) {
                if (!vert.position) continue;
                auto pos_it = pos_map.find(vert.position);
                if (pos_it != pos_map.end()) continue;
                auto nverts = (int)shape.positions.size();
                pos_map.insert(pos_it, {vert.position, nverts});
                shape.positions.push_back(opos.at(vert.position - 1));
            }
            for (auto& vert : verts) {
                if (!vert.texturecoord) continue;
                auto texcoord_it = texcoord_map.find(vert.texturecoord);
                if (texcoord_it != texcoord_map.end()) continue;
                auto nverts = (int)shape.texturecoords.size();
                texcoord_map.insert(texcoord_it, {vert.texturecoord, nverts});
                shape.texturecoords.push_back(
                    otexcoord.at(vert.texturecoord - 1));
            }
            for (auto& vert : verts) {
                if (!vert.normal) continue;
                auto norm_it = norm_map.find(vert.normal);
                if (norm_it != norm_map.end()) continue;
                auto nverts = (int)shape.normals.size();
                norm_map.insert(norm_it, {vert.normal, nverts});
                shape.normals.push_back(onorm.at(vert.normal - 1));
            }
        }

        // callbacks
        void vert(const vec3f& v) { opos.push_back(v); }
        void norm(const vec3f& v) { onorm.push_back(v); }
        void texcoord(const vec2f& v) { otexcoord.push_back(v); }
        void face(const vector<obj_vertex>& verts) {
            if (scene.shapes.empty()) add_shape();
            if (!scene.shapes.back().positions.empty() &&
                (!scene.shapes.back().lines.empty() ||
                    !scene.shapes.back().points.empty())) {
                add_shape();
            }
            auto& shape = scene.shapes.back();
            if (!shape.preserve_facevarying) {
                add_verts(verts, shape);
                if (verts.size() == 4) {
                    shape.quads.push_back(
                        {vertex_map.at(verts[0]), vertex_map.at(verts[1]),
                            vertex_map.at(verts[2]), vertex_map.at(verts[3])});
                } else {
                    for (auto i = 2; i < verts.size(); i++)
                        shape.triangles.push_back({vertex_map.at(verts[0]),
                            vertex_map.at(verts[i - 1]),
                            vertex_map.at(verts[i])});
                }
            } else {
                add_fvverts(verts, shape);
                if (verts.size() == 4) {
                    if (verts[0].position) {
                        shape.quads_positions.push_back(
                            {pos_map.at(verts[0].position),
                                pos_map.at(verts[1].position),
                                pos_map.at(verts[2].position),
                                pos_map.at(verts[3].position)});
                    }
                    if (verts[0].texturecoord) {
                        shape.quads_texturecoords.push_back(
                            {texcoord_map.at(verts[0].texturecoord),
                                texcoord_map.at(verts[1].texturecoord),
                                texcoord_map.at(verts[2].texturecoord),
                                texcoord_map.at(verts[3].texturecoord)});
                    }
                    if (verts[0].normal) {
                        shape.quads_normals.push_back(
                            {norm_map.at(verts[0].normal),
                                norm_map.at(verts[1].normal),
                                norm_map.at(verts[2].normal),
                                norm_map.at(verts[3].normal)});
                    }
                } else {
                    if (verts[0].position) {
                        for (auto i = 2; i < verts.size(); i++)
                            shape.quads_positions.push_back(
                                {pos_map.at(verts[0].position),
                                    pos_map.at(verts[i - 1].position),
                                    pos_map.at(verts[i].position),
                                    pos_map.at(verts[i].position)});
                    }
                    if (verts[0].texturecoord) {
                        for (auto i = 2; i < verts.size(); i++)
                            shape.quads_texturecoords.push_back(
                                {texcoord_map.at(verts[0].texturecoord),
                                    texcoord_map.at(verts[i - 1].texturecoord),
                                    texcoord_map.at(verts[i].texturecoord),
                                    texcoord_map.at(verts[i].texturecoord)});
                    }
                    if (verts[0].normal) {
                        for (auto i = 2; i < verts.size(); i++)
                            shape.quads_normals.push_back(
                                {norm_map.at(verts[0].normal),
                                    norm_map.at(verts[i - 1].normal),
                                    norm_map.at(verts[i].normal),
                                    norm_map.at(verts[i].normal)});
                    }
                }
            }
        }
        void line(const vector<obj_vertex>& verts) {
            if (scene.shapes.empty()) add_shape();
            if (!scene.shapes.back().positions.empty() &&
                scene.shapes.back().lines.empty()) {
                add_shape();
            }
            auto& shape                = scene.shapes.back();
            shape.preserve_facevarying = false;
            add_verts(verts, shape);
            for (auto i = 1; i < verts.size(); i++)
                shape.lines.push_back(
                    {vertex_map.at(verts[i - 1]), vertex_map.at(verts[i])});
        }
        void point(const vector<obj_vertex>& verts) {
            if (scene.shapes.empty()) add_shape();
            if (!scene.shapes.back().positions.empty() &&
                scene.shapes.back().points.empty()) {
                add_shape();
            }
            auto& shape                = scene.shapes.back();
            shape.preserve_facevarying = false;
            add_verts(verts, shape);
            for (auto i = 0; i < verts.size(); i++)
                shape.points.push_back(vertex_map.at(verts[i]));
        }
        void object(const string& name) {
            oname = name;
            gname = "";
            mname = "";
            add_shape();
        }
        void group(const string& name) {
            gname = name;
            add_shape();
        }
        void usemtl(const string& name) {
            mname = name;
            add_shape();
        }
        void material(const obj_material& omat) {
            auto material                   = yocto_material();
            material.uri                    = omat.name;
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
            material.roughness_texture      = add_texture(omat.rs_txt, true);
            material.displacement_texture   = add_texture(omat.disp_txt, true);
            material.normal_texture         = add_texture(omat.norm_txt, true);
            material.volume_emission        = omat.ve;
            material.volume_albedo          = omat.va;
            material.volume_density         = omat.vd;
            material.volume_phaseg          = omat.vg;
            material.volume_density_texture = add_voltexture(
                omat.vd_txt, false);
            scene.materials.push_back(material);
            mmap[material.uri] = (int)scene.materials.size() - 1;
        }
        void camera(const obj_camera& ocam) {
            auto camera           = yocto_camera();
            camera.uri            = ocam.name;
            camera.frame          = ocam.frame;
            camera.orthographic   = ocam.ortho;
            camera.film_width     = ocam.width;
            camera.film_height    = ocam.height;
            camera.focal_length   = ocam.focal;
            camera.focus_distance = ocam.focus;
            camera.lens_aperture  = ocam.aperture;
            scene.cameras.push_back(camera);
        }
        void environmnet(const obj_environment& oenv) {
            auto environment             = yocto_environment();
            environment.uri              = oenv.name;
            environment.frame            = oenv.frame;
            environment.emission         = oenv.ke;
            environment.emission_texture = add_texture(oenv.ke_txt, true);
            scene.environments.push_back(environment);
        }
        void procedural(const obj_procedural& oproc) {
            auto shape = yocto_shape();
            shape.uri  = oproc.name;
            if (oproc.type == "floor") {
                make_floor_shape(shape.quads, shape.positions, shape.normals,
                    shape.texturecoords,
                    {oproc.level < 0 ? 1 : pow2(oproc.level),
                        oproc.level < 0 ? 20 : pow2(oproc.level)},
                    {oproc.size, oproc.size}, {oproc.size / 2, oproc.size / 2},
                    identity_frame3f);
            } else {
                throw io_error("unknown obj procedural");
            }
            scene.shapes.push_back(shape);
            auto instance  = yocto_instance{};
            instance.uri   = shape.uri;
            instance.shape = (int)scene.shapes.size() - 1;
            if (mmap.find(oproc.material) == mmap.end()) {
                throw io_error("missing material " + oproc.material);
            } else {
                instance.material = mmap.find(oproc.material)->second;
            }
            scene.instances.push_back(instance);
        }
    };

    try {
        // Parse obj
        auto obj_options          = load_obj_options();
        obj_options.geometry_only = false;
        auto cb                   = parse_callbacks{scene, options};
        load_obj(filename, cb, obj_options);

        // cleanup empty
        for (auto idx = 0; idx < scene.shapes.size(); idx++) {
            scene.instances[idx].shape = idx;
            if (!scene.shapes[idx].positions.empty()) continue;
            scene.shapes.erase(scene.shapes.begin() + idx);
            scene.instances.erase(scene.instances.begin() + idx);
            idx--;
        }

        // check if any empty shape is left
        for (auto& shape : scene.shapes) {
            if (shape.positions.empty())
                throw io_error("empty shapes not supported");
        }

        // merging quads and triangles
        for (auto& shape : scene.shapes) {
            if (shape.triangles.empty() || shape.quads.empty()) continue;
            merge_triangles_and_quads(shape.triangles, shape.quads, false);
        }

        // load textures
        auto dirname = get_dirname(filename);
        load_scene_textures(scene, dirname, options);
    } catch (const std::exception& e) {
        throw io_error("cannot load scene " + filename + "\n" + e.what());
    }

    // fix scene
    scene.uri = get_filename(filename);
    add_missing_cameras(scene);
    add_missing_materials(scene);
    normalize_uris(scene);
    trim_memory(scene);
    update_transforms(scene);
}

void save_mtl(
    const string& filename, const yocto_scene& scene, bool flip_tr = true) {
    // open file
    auto fs = output_file(filename);

    // for each material, dump all the values
    for (auto& material : scene.materials) {
        println_values(fs, "newmtl", get_basename(material.uri));
        println_values(fs, "  illum 2");
        println_values(fs, "  Ke", material.emission);
        println_values(fs, "  Kd", material.diffuse);
        println_values(fs, "  Ks", material.specular);
        println_values(fs, "  Kt", material.transmission);
        println_values(fs, "  Ns",
            (int)clamp(
                2 / pow(clamp(material.roughness, 0.0f, 0.99f) + 1e-10f, 4.0f) -
                    2,
                0.0f, 1.0e9f));
        println_values(fs, "  d", material.opacity);
        if (material.emission_texture >= 0)
            println_values(
                fs, "  map_Ke", scene.textures[material.emission_texture].uri);
        if (material.diffuse_texture >= 0)
            println_values(
                fs, "  map_Kd", scene.textures[material.diffuse_texture].uri);
        if (material.specular_texture >= 0)
            println_values(
                fs, "  map_Ks", scene.textures[material.specular_texture].uri);
        if (material.transmission_texture >= 0)
            println_values(fs, "  map_Kt",
                scene.textures[material.transmission_texture].uri);
        if (material.roughness_texture >= 0)
            println_values(
                fs, "  map_Pr", scene.textures[material.roughness_texture].uri);
        if (material.displacement_texture >= 0)
            println_values(fs, "  map_disp",
                scene.textures[material.displacement_texture].uri);
        if (material.normal_texture >= 0)
            println_values(
                fs, "  map_norm", scene.textures[material.normal_texture].uri);
        if (material.volume_emission != zero3f)
            println_values(fs, "  Ve", material.volume_emission);
        if (material.volume_density != zero3f)
            println_values(fs, "  Vd", material.volume_density);
        if (material.volume_albedo != zero3f)
            println_values(fs, "  Va", material.volume_albedo);
        if (material.volume_phaseg != 0)
            println_values(fs, "  Vg", material.volume_phaseg);
        if (material.volume_density_texture >= 0)
            println_values(fs, "  map_Vd",
                scene.voltextures[material.volume_density_texture].uri);
        println_values(fs, "\n");
    }
}

void save_objx(const string& filename, const yocto_scene& scene) {
    // open file
    auto fs = output_file(filename);

    // cameras
    for (auto& camera : scene.cameras) {
        println_values(fs, "c", get_basename(camera.uri),
            (int)camera.orthographic, camera.film_width, camera.film_height,
            camera.focal_length, camera.focus_distance, camera.lens_aperture,
            camera.frame);
    }

    // environments
    for (auto& environment : scene.environments) {
        if (environment.emission_texture >= 0) {
            println_values(fs, "e", get_basename(environment.uri),
                environment.emission,
                scene.textures[environment.emission_texture].uri,
                environment.frame);
        } else {
            println_values(fs, "e", environment.uri, environment.emission,
                "\"\"", environment.frame);
        }
    }
}

inline void print_value(const output_file& fs, const obj_vertex& value) {
    print_value(fs, value.position);
    if (value.texturecoord) {
        print_value(fs, '/');
        print_value(fs, value.texturecoord);
        if (value.normal) {
            print_value(fs, '/');
            print_value(fs, value.normal);
        }
    } else {
        if (value.normal) {
            print_value(fs, '/');
            print_value(fs, '/');
            print_value(fs, value.normal);
        }
    }
}

void save_obj(const string& filename, const yocto_scene& scene,
    bool flip_texcoord = true) {
    // open file
    auto fs = output_file(filename);

    println_values(
        fs, "# Saved by Yocto/GL - https://github.com/xelatihy/yocto-gl\n\n");

    // material library
    if (!scene.materials.empty()) {
        auto mtlname = get_noextension(get_filename(filename)) + ".mtl";
        println_values(fs, "mtllib", mtlname);
    }

    // shapes
    auto offset = obj_vertex{0, 0, 0};
    for (auto& instance : scene.instances) {
        auto& shape = scene.shapes[instance.shape];
        println_values(fs, "o", instance.uri);
        if (instance.material >= 0)
            println_values(fs, "usemtl",
                get_basename(scene.materials[instance.material].uri));
        if (instance.frame == identity_frame3f) {
            for (auto& p : shape.positions) println_values(fs, "v", p);
            for (auto& n : shape.normals) println_values(fs, "vn", n);
            for (auto& t : shape.texturecoords)
                println_values(
                    fs, "vt", vec2f{t.x, (flip_texcoord) ? 1 - t.y : t.y});
        } else {
            for (auto& pp : shape.positions) {
                println_values(fs, "v", transform_point(instance.frame, pp));
            }
            for (auto& nn : shape.normals) {
                println_values(
                    fs, "vn", transform_direction(instance.frame, nn));
            }
            for (auto& t : shape.texturecoords)
                println_values(
                    fs, "vt", vec2f{t.x, (flip_texcoord) ? 1 - t.y : t.y});
        }
        auto mask = obj_vertex{1, shape.texturecoords.empty() ? 0 : 1,
            shape.normals.empty() ? 0 : 1};
        auto vert = [mask, offset](int i) {
            return obj_vertex{(i + offset.position + 1) * mask.position,
                (i + offset.texturecoord + 1) * mask.texturecoord,
                (i + offset.normal + 1) * mask.normal};
        };
        for (auto& p : shape.points) {
            println_values(fs, "p", vert(p));
        }
        for (auto& l : shape.lines) {
            println_values(fs, "l", vert(l.x), vert(l.y));
        }
        for (auto& t : shape.triangles) {
            println_values(fs, "f", vert(t.x), vert(t.y), vert(t.z));
        }
        for (auto& q : shape.quads) {
            if (q.z == q.w) {
                println_values(fs, "f", vert(q.x), vert(q.y), vert(q.z));
            } else {
                println_values(
                    fs, "f", vert(q.x), vert(q.y), vert(q.z), vert(q.w));
            }
        }
        for (auto i = 0; i < shape.quads_positions.size(); i++) {
            if (!shape.texturecoords.empty() && shape.normals.empty()) {
                auto vert = [offset](int ip, int it) {
                    return obj_vertex{ip + offset.position + 1,
                        it + offset.texturecoord + 1, 0};
                };
                auto qp = shape.quads_positions[i];
                auto qt = shape.quads_texturecoords[i];
                if (qp.z == qp.w) {
                    println_values(fs, "f", vert(qp.x, qt.x), vert(qp.y, qt.y),
                        vert(qp.z, qt.z));
                } else {
                    println_values(fs, "f", vert(qp.x, qt.x), vert(qp.y, qt.y),
                        vert(qp.z, qt.z), vert(qp.w, qt.w));
                }
            } else if (!shape.texturecoords.empty() && !shape.normals.empty()) {
                auto vert = [offset](int ip, int it, int in) {
                    return obj_vertex{ip + offset.position + 1,
                        it + offset.texturecoord + 1, in + offset.normal + 1};
                };
                auto qp = shape.quads_positions[i];
                auto qt = shape.quads_texturecoords[i];
                auto qn = shape.quads_normals[i];
                if (qp.z == qp.w) {
                    println_values(fs, "f", vert(qp.x, qt.x, qn.x),
                        vert(qp.y, qt.y, qn.y), vert(qp.z, qt.z, qn.z));
                } else {
                    println_values(fs, "f", vert(qp.x, qt.x, qn.x),
                        vert(qp.y, qt.y, qn.y), vert(qp.z, qt.z, qn.z),
                        vert(qp.w, qt.w, qn.w));
                }
            } else if (!shape.normals.empty()) {
                auto vert = [offset](int ip, int in) {
                    return obj_vertex{
                        ip + offset.position + 1, 0, in + offset.normal + 1};
                };
                auto qp = shape.quads_positions[i];
                auto qn = shape.quads_normals[i];
                if (qp.z == qp.w) {
                    println_values(fs, "f", vert(qp.x, qn.x), vert(qp.y, qn.y),
                        vert(qp.z, qn.z));
                } else {
                    println_values(fs, "f", vert(qp.x, qn.x), vert(qp.y, qn.y),
                        vert(qp.z, qn.z), vert(qp.w, qn.w));
                }
            } else {
                auto vert = [offset](int ip) {
                    return obj_vertex{ip + offset.position + 1, 0, 0};
                };
                auto q = shape.quads_positions[i];
                if (q.z == q.w) {
                    println_values(fs, "f", vert(q.x), vert(q.y), vert(q.z));
                } else {
                    println_values(
                        fs, "f", vert(q.x), vert(q.y), vert(q.z), vert(q.w));
                }
            }
        }
        offset.position += shape.positions.size();
        offset.texturecoord += shape.texturecoords.size();
        offset.normal += shape.normals.size();
    }
}

void save_obj_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options) {
    try {
        save_obj(filename, scene, true);
        if (!scene.materials.empty()) {
            save_mtl(get_noextension(filename) + ".mtl", scene, true);
        }
        if (!scene.cameras.empty() || !scene.environments.empty()) {
            save_objx(get_noextension(filename) + ".objx", scene);
        }

        // skip textures if needed
        auto dirname = get_dirname(filename);
        save_scene_textures(scene, dirname, options);
    } catch (const std::exception& e) {
        throw io_error("cannot save scene " + filename + "\n" + e.what());
    }
}

void print_obj_camera(const yocto_camera& camera) {
    println_values(stdout, "c", camera.uri, (int)camera.orthographic,
        camera.film_width, camera.film_height, camera.focal_length,
        camera.focus_distance, camera.lens_aperture, camera.frame);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PLY CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

void load_ply_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options) {
    scene = {};

    try {
        // load ply mesh
        scene.shapes.push_back({});
        auto& shape = scene.shapes.back();
        load_shape(filename, shape.points, shape.lines, shape.triangles,
            shape.quads, shape.quads_positions, shape.quads_normals,
            shape.quads_texturecoords, shape.positions, shape.normals,
            shape.texturecoords, shape.colors, shape.radius, false);

        // add instance
        auto instance  = yocto_instance{};
        instance.uri   = shape.uri;
        instance.shape = 0;
        scene.instances.push_back(instance);

    } catch (const std::exception& e) {
        throw io_error("cannot load scene " + filename + "\n" + e.what());
    }

    // fix scene
    scene.uri = get_filename(filename);
    add_missing_cameras(scene);
    add_missing_materials(scene);
    normalize_uris(scene);
    trim_memory(scene);
    update_transforms(scene);
}

void save_ply_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options) {
    if (scene.shapes.empty()) {
        throw io_error("cannot save empty scene " + filename);
    }
    try {
        auto& shape = scene.shapes.front();
        save_shape(filename, shape.points, shape.lines, shape.triangles,
            shape.quads, shape.quads_positions, shape.quads_normals,
            shape.quads_texturecoords, shape.positions, shape.normals,
            shape.texturecoords, shape.colors, shape.radius);
    } catch (const std::exception& e) {
        throw io_error("cannot save scene " + filename + "\n" + e.what());
    }
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
void gltf_to_scene(const string& filename, yocto_scene& scene) {
    // load gltf
    auto options = cgltf_options{};
    memset(&options, 0, sizeof(options));
    auto data   = (cgltf_data*)nullptr;
    auto result = cgltf_parse_file(&options, filename.c_str(), &data);
    if (result != cgltf_result_success) {
        throw io_error("could not load gltf " + filename);
    }
    auto gltf = unique_ptr<cgltf_data, void (*)(cgltf_data*)>{data, cgltf_free};
    if (cgltf_load_buffers(&options, data, get_dirname(filename).c_str()) !=
        cgltf_result_success) {
        throw io_error("could not load gltf buffers " + filename);
    }

    // convert textures
    auto imap = unordered_map<cgltf_image*, int>{};
    for (auto tid = 0; tid < gltf->images_count; tid++) {
        auto gimg    = &gltf->images[tid];
        auto texture = yocto_texture{};
        texture.uri  = (startswith(gimg->uri, "data:"))
                          ? string("[glTF-inline].png")
                          : gimg->uri;
        scene.textures.push_back(texture);
        imap[gimg] = tid;
    }

    // add a texture
    auto add_texture = [&imap](
                           const cgltf_texture_view& ginfo, bool force_linear) {
        if (!ginfo.texture || !ginfo.texture->image) return -1;
        auto gtxt = ginfo.texture;
        return imap.at(gtxt->image);
    };

    // convert materials
    auto mmap = unordered_map<cgltf_material*, int>{{nullptr, -1}};
    for (auto mid = 0; mid < gltf->materials_count; mid++) {
        auto gmat         = &gltf->materials[mid];
        auto material     = yocto_material();
        material.uri      = gmat->name ? gmat->name : "";
        material.emission = {gmat->emissive_factor[0], gmat->emissive_factor[1],
            gmat->emissive_factor[2]};
        material.emission_texture = add_texture(gmat->emissive_texture, false);
        if (gmat->has_pbr_specular_glossiness) {
            material.base_metallic = false;
            material.gltf_textures = true;
            auto gsg               = &gmat->pbr_specular_glossiness;
            auto kb = vec4f{gsg->diffuse_factor[0], gsg->diffuse_factor[1],
                gsg->diffuse_factor[2], gsg->diffuse_factor[3]};
            material.diffuse         = {kb.x, kb.y, kb.z};
            material.opacity         = kb.w;
            material.specular        = {gsg->specular_factor[0],
                gsg->specular_factor[1], gsg->specular_factor[2]};
            material.roughness       = 1 - gsg->glossiness_factor;
            material.diffuse_texture = add_texture(gsg->diffuse_texture, false);
            material.specular_texture = add_texture(
                gsg->specular_glossiness_texture, false);
            material.roughness_texture = material.specular_texture;
        } else if (gmat->has_pbr_metallic_roughness) {
            material.base_metallic   = true;
            material.gltf_textures   = true;
            auto gmr                 = &gmat->pbr_metallic_roughness;
            auto kb                  = vec4f{gmr->base_color_factor[0],
                gmr->base_color_factor[1], gmr->base_color_factor[2],
                gmr->base_color_factor[3]};
            material.diffuse         = {kb.x, kb.y, kb.z};
            material.opacity         = kb.w;
            auto km                  = gmr->metallic_factor;
            material.specular        = {km, km, km};
            material.roughness       = gmr->roughness_factor;
            material.diffuse_texture = add_texture(
                gmr->base_color_texture, false);
            material.specular_texture = add_texture(
                gmr->metallic_roughness_texture, true);
            material.roughness_texture = material.specular_texture;
        }
        material.normal_texture = add_texture(gmat->normal_texture, true);
        scene.materials.push_back(material);
        mmap[gmat] = (int)scene.materials.size() - 1;
    }

    // get values from accessors
    auto accessor_values =
        [](const cgltf_accessor* gacc,
            bool normalize = false) -> vector<std::array<double, 4>> {
        auto gview       = gacc->buffer_view;
        auto data        = (byte*)gview->buffer->data;
        auto offset      = gacc->offset + gview->offset;
        auto stride      = gview->stride;
        auto compTypeNum = gacc->component_type;
        auto count       = gacc->count;
        auto type        = gacc->type;
        auto ncomp       = 0;
        if (type == cgltf_type_scalar) ncomp = 1;
        if (type == cgltf_type_vec2) ncomp = 2;
        if (type == cgltf_type_vec3) ncomp = 3;
        if (type == cgltf_type_vec4) ncomp = 4;
        auto compSize = 1;
        if (compTypeNum == cgltf_component_type_r_16 ||
            compTypeNum == cgltf_component_type_r_16u) {
            compSize = 2;
        }
        if (compTypeNum == cgltf_component_type_r_32u ||
            compTypeNum == cgltf_component_type_r_32f) {
            compSize = 4;
        }
        if (!stride) stride = compSize * ncomp;
        auto vals = vector<std::array<double, 4>>(
            count, {{0.0, 0.0, 0.0, 1.0}});
        for (auto i = 0; i < count; i++) {
            auto d = data + offset + i * stride;
            for (auto c = 0; c < ncomp; c++) {
                if (compTypeNum == cgltf_component_type_r_8) {  // char
                    vals[i][c] = (double)(*(char*)d);
                    if (normalize) vals[i][c] /= SCHAR_MAX;
                } else if (compTypeNum == cgltf_component_type_r_8u) {  // byte
                    vals[i][c] = (double)(*(byte*)d);
                    if (normalize) vals[i][c] /= UCHAR_MAX;
                } else if (compTypeNum == cgltf_component_type_r_16) {  // short
                    vals[i][c] = (double)(*(short*)d);
                    if (normalize) vals[i][c] /= SHRT_MAX;
                } else if (compTypeNum ==
                           cgltf_component_type_r_16u) {  // unsigned short
                    vals[i][c] = (double)(*(unsigned short*)d);
                    if (normalize) vals[i][c] /= USHRT_MAX;
                } else if (compTypeNum ==
                           cgltf_component_type_r_32u) {  // unsigned int
                    vals[i][c] = (double)(*(unsigned int*)d);
                    if (normalize) vals[i][c] /= UINT_MAX;
                } else if (compTypeNum ==
                           cgltf_component_type_r_32f) {  // float
                    vals[i][c] = (*(float*)d);
                }
                d += compSize;
            }
        }
        return vals;
    };

    // convert meshes
    auto meshes = unordered_map<cgltf_mesh*, vector<vec2i>>{{nullptr, {}}};
    for (auto mid = 0; mid < gltf->meshes_count; mid++) {
        auto gmesh    = &gltf->meshes[mid];
        meshes[gmesh] = {};
        for (auto sid = 0; sid < gmesh->primitives_count; sid++) {
            auto gprim = &gmesh->primitives[sid];
            if (!gprim->attributes_count) continue;
            auto shape = yocto_shape();
            shape.uri  = (gmesh->name ? gmesh->name : "") +
                        ((sid) ? to_string(sid) : string());
            for (auto aid = 0; aid < gprim->attributes_count; aid++) {
                auto gattr    = &gprim->attributes[aid];
                auto semantic = string(gattr->name ? gattr->name : "");
                auto gacc     = gattr->data;
                auto vals     = accessor_values(gacc);
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
                } else if (semantic == "TEXCOORD" || semantic == "TEXCOORD_0") {
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
            if (!gprim->indices) {
                if (gprim->type == cgltf_primitive_type_triangles) {
                    shape.triangles.reserve(shape.positions.size() / 3);
                    for (auto i = 0; i < shape.positions.size() / 3; i++)
                        shape.triangles.push_back(
                            {i * 3 + 0, i * 3 + 1, i * 3 + 2});
                } else if (gprim->type == cgltf_primitive_type_triangle_fan) {
                    shape.triangles.reserve(shape.positions.size() - 2);
                    for (auto i = 2; i < shape.positions.size(); i++)
                        shape.triangles.push_back({0, i - 1, i});
                } else if (gprim->type == cgltf_primitive_type_triangle_strip) {
                    shape.triangles.reserve(shape.positions.size() - 2);
                    for (auto i = 2; i < shape.positions.size(); i++)
                        shape.triangles.push_back({i - 2, i - 1, i});
                } else if (gprim->type == cgltf_primitive_type_lines) {
                    shape.lines.reserve(shape.positions.size() / 2);
                    for (auto i = 0; i < shape.positions.size() / 2; i++)
                        shape.lines.push_back({i * 2 + 0, i * 2 + 1});
                } else if (gprim->type == cgltf_primitive_type_line_loop) {
                    shape.lines.reserve(shape.positions.size());
                    for (auto i = 1; i < shape.positions.size(); i++)
                        shape.lines.push_back({i - 1, i});
                    shape.lines.back() = {(int)shape.positions.size() - 1, 0};
                } else if (gprim->type == cgltf_primitive_type_line_strip) {
                    shape.lines.reserve(shape.positions.size() - 1);
                    for (auto i = 1; i < shape.positions.size(); i++)
                        shape.lines.push_back({i - 1, i});
                } else if (gprim->type == cgltf_primitive_type_points) {
                    // points
                    throw io_error("points not supported");
                } else {
                    throw io_error("unknown primitive type");
                }
            } else {
                auto indices = accessor_values(gprim->indices);
                if (gprim->type == cgltf_primitive_type_triangles) {
                    shape.triangles.reserve(indices.size() / 3);
                    for (auto i = 0; i < indices.size() / 3; i++)
                        shape.triangles.push_back({(int)indices[i * 3 + 0][0],
                            (int)indices[i * 3 + 1][0],
                            (int)indices[i * 3 + 2][0]});
                } else if (gprim->type == cgltf_primitive_type_triangle_fan) {
                    shape.triangles.reserve(indices.size() - 2);
                    for (auto i = 2; i < indices.size(); i++)
                        shape.triangles.push_back({(int)indices[0][0],
                            (int)indices[i - 1][0], (int)indices[i][0]});
                } else if (gprim->type == cgltf_primitive_type_triangle_strip) {
                    shape.triangles.reserve(indices.size() - 2);
                    for (auto i = 2; i < indices.size(); i++)
                        shape.triangles.push_back({(int)indices[i - 2][0],
                            (int)indices[i - 1][0], (int)indices[i][0]});
                } else if (gprim->type == cgltf_primitive_type_lines) {
                    shape.lines.reserve(indices.size() / 2);
                    for (auto i = 0; i < indices.size() / 2; i++)
                        shape.lines.push_back({(int)indices[i * 2 + 0][0],
                            (int)indices[i * 2 + 1][0]});
                } else if (gprim->type == cgltf_primitive_type_line_loop) {
                    shape.lines.reserve(indices.size());
                    for (auto i = 1; i < indices.size(); i++)
                        shape.lines.push_back(
                            {(int)indices[i - 1][0], (int)indices[i][0]});
                    shape.lines.back() = {(int)indices[indices.size() - 1][0],
                        (int)indices[0][0]};
                } else if (gprim->type == cgltf_primitive_type_line_strip) {
                    shape.lines.reserve(indices.size() - 1);
                    for (auto i = 1; i < indices.size(); i++)
                        shape.lines.push_back(
                            {(int)indices[i - 1][0], (int)indices[i][0]});
                } else if (gprim->type == cgltf_primitive_type_points) {
                    throw io_error("points not supported");
                } else {
                    throw io_error("unknown primitive type");
                }
            }
            scene.shapes.push_back(shape);
            meshes[gmesh].push_back(
                {(int)scene.shapes.size() - 1, mmap.at(gprim->material)});
        }
    }

    // convert cameras
    auto cmap = unordered_map<cgltf_camera*, int>{{nullptr, -1}};
    for (auto cid = 0; cid < gltf->cameras_count; cid++) {
        auto gcam           = &gltf->cameras[cid];
        auto camera         = yocto_camera{};
        camera.uri          = gcam->name ? gcam->name : "";
        camera.orthographic = gcam->type == cgltf_camera_type_orthographic;
        if (camera.orthographic) {
            // throw io_error("orthographic not supported well");
            auto ortho           = &gcam->orthographic;
            camera.lens_aperture = 0;
            camera.orthographic  = true;
            camera.film_width    = ortho->xmag;
            camera.film_height   = ortho->ymag;
        } else {
            auto persp           = &gcam->perspective;
            camera.lens_aperture = 0;
            set_camera_perspectivey(
                camera, persp->yfov, persp->aspect_ratio, float_max);
        }
        scene.cameras.push_back(camera);
        cmap[gcam] = (int)scene.cameras.size() - 1;
    }

    // convert nodes
    auto nmap = unordered_map<cgltf_node*, int>{{nullptr, -1}};
    for (auto nid = 0; nid < gltf->nodes_count; nid++) {
        auto gnde = &gltf->nodes[nid];
        auto node = yocto_scene_node{};
        node.uri  = gnde->name ? gnde->name : "";
        if (gnde->camera) node.camera = cmap.at(gnde->camera);
        if (gnde->has_translation) {
            node.translation = {gnde->translation[0], gnde->translation[1],
                gnde->translation[2]};
        }
        if (gnde->has_rotation) {
            node.rotation = {gnde->rotation[0], gnde->rotation[1],
                gnde->rotation[2], gnde->rotation[3]};
        }
        if (gnde->has_scale) {
            node.scale = {gnde->scale[0], gnde->scale[1], gnde->scale[2]};
        }
        if (gnde->has_matrix) {
            auto m     = gnde->matrix;
            node.local = frame3f(
                mat4f{{m[0], m[1], m[2], m[3]}, {m[4], m[5], m[6], m[7]},
                    {m[8], m[9], m[10], m[11]}, {m[12], m[13], m[14], m[15]}});
        }
        scene.nodes.push_back(node);
        nmap[gnde] = (int)scene.nodes.size();
    }

    // set up parent pointers
    for (auto nid = 0; nid < gltf->nodes_count; nid++) {
        auto gnde = &gltf->nodes[nid];
        if (!gnde->children_count) continue;
        for (auto cid = 0; cid < gnde->children_count; cid++) {
            scene.nodes[nmap.at(gnde->children[cid])].parent = nid;
        }
    }

    // set up instances
    for (auto nid = 0; nid < gltf->nodes_count; nid++) {
        auto gnde = &gltf->nodes[nid];
        if (!gnde->mesh) continue;
        auto& node = scene.nodes[nid];
        auto& shps = meshes.at(gnde->mesh);
        if (shps.empty()) continue;
        if (shps.size() == 1) {
            auto instance     = yocto_instance();
            instance.uri      = node.uri;
            instance.shape    = shps[0].x;
            instance.material = shps[0].y;
            scene.instances.push_back(instance);
            node.instance = (int)scene.instances.size() - 1;
        } else {
            for (auto shp : shps) {
                auto& shape       = scene.shapes[shp.x];
                auto  instance    = yocto_instance();
                instance.uri      = node.uri + "_" + shape.uri;
                instance.shape    = shp.x;
                instance.material = shp.y;
                scene.instances.push_back(instance);
                auto child     = yocto_scene_node{};
                child.uri      = node.uri + "_" + shape.uri;
                child.parent   = nid;
                child.instance = (int)scene.instances.size() - 1;
                scene.nodes.push_back(child);
            }
        }
    }

    // hasher for later
    struct sampler_map_hash {
        size_t operator()(
            const pair<cgltf_animation_sampler*, cgltf_animation_path_type>&
                value) const {
            auto hasher1 = std::hash<cgltf_animation_sampler*>();
            auto hasher2 = std::hash<int>();
            auto h       = (size_t)0;
            h ^= hasher1(value.first) + 0x9e3779b9 + (h << 6) + (h >> 2);
            h ^= hasher2(value.second) + 0x9e3779b9 + (h << 6) + (h >> 2);
            return h;
        }
    };

    // convert animations
    for (auto gid = 0; gid < gltf->animations_count; gid++) {
        auto ganm        = &gltf->animations[gid];
        auto aid         = 0;
        auto sampler_map = unordered_map<
            pair<cgltf_animation_sampler*, cgltf_animation_path_type>, int,
            sampler_map_hash>();
        for (auto cid = 0; cid < ganm->channels_count; cid++) {
            auto gchannel = &ganm->channels[cid];
            auto path     = gchannel->target_path;
            if (sampler_map.find({gchannel->sampler, path}) ==
                sampler_map.end()) {
                auto gsampler  = gchannel->sampler;
                auto animation = yocto_animation{};
                animation.uri  = (ganm->name ? ganm->name : "anim") +
                                to_string(aid++);
                animation.animation_group = ganm->name ? ganm->name : "";
                auto input_view           = accessor_values(gsampler->input);
                animation.keyframes_times.resize(input_view.size());
                for (auto i = 0; i < input_view.size(); i++)
                    animation.keyframes_times[i] = input_view[i][0];
                switch (gsampler->interpolation) {
                    case cgltf_interpolation_type_linear:
                        animation.interpolation_type =
                            yocto_interpolation_type::linear;
                        break;
                    case cgltf_interpolation_type_step:
                        animation.interpolation_type =
                            yocto_interpolation_type::step;
                        break;
                    case cgltf_interpolation_type_cubic_spline:
                        animation.interpolation_type =
                            yocto_interpolation_type::bezier;
                        break;
                }
                auto output_view = accessor_values(gsampler->output);
                switch (path) {
                    case cgltf_animation_path_type_translation: {
                        animation.translation_keyframes.reserve(
                            output_view.size());
                        for (auto i = 0; i < output_view.size(); i++)
                            animation.translation_keyframes.push_back(
                                {(float)output_view[i][0],
                                    (float)output_view[i][1],
                                    (float)output_view[i][2]});
                    } break;
                    case cgltf_animation_path_type_rotation: {
                        animation.rotation_keyframes.reserve(
                            output_view.size());
                        for (auto i = 0; i < output_view.size(); i++)
                            animation.rotation_keyframes.push_back(
                                {(float)output_view[i][0],
                                    (float)output_view[i][1],
                                    (float)output_view[i][2],
                                    (float)output_view[i][3]});
                    } break;
                    case cgltf_animation_path_type_scale: {
                        animation.scale_keyframes.reserve(output_view.size());
                        for (auto i = 0; i < output_view.size(); i++)
                            animation.scale_keyframes.push_back(
                                {(float)output_view[i][0],
                                    (float)output_view[i][1],
                                    (float)output_view[i][2]});
                    } break;
                    case cgltf_animation_path_type_weights: {
                        throw io_error("weights not supported for now");
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
                    default: {
                        throw io_error("bad gltf animation");
                    }
                }
                sampler_map[{gchannel->sampler, path}] =
                    (int)scene.animations.size();
                scene.animations.push_back(animation);
            }
            scene.animations[sampler_map.at({gchannel->sampler, path})]
                .node_targets.push_back(nmap.at(gchannel->target_node));
        }
    }
}

// Load a scene
void load_gltf_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options) {
    // initialization
    scene = {};

    try {
        // load gltf
        gltf_to_scene(filename, scene);

        // load textures
        auto dirname = get_dirname(filename);
        load_scene_textures(scene, dirname, options);

    } catch (const std::exception& e) {
        throw io_error("cannot load scene " + filename + "\n" + e.what());
    }

    // fix scene
    scene.uri = get_filename(filename);
    add_missing_cameras(scene);
    add_missing_materials(scene);
    normalize_uris(scene);
    trim_memory(scene);
    update_transforms(scene);

    // fix cameras
    auto bbox = compute_scene_bounds(scene);
    for (auto& camera : scene.cameras) {
        auto center   = (bbox.min + bbox.max) / 2;
        auto distance = dot(-camera.frame.z, center - camera.frame.o);
        if (distance > 0) camera.focus_distance = distance;
    }
}

// convert gltf scene to json
void scene_to_gltf(const yocto_scene& scene, json& js) {
    // init to emprt object
    js = json::object();

    // start creating json
    js["asset"]["version"] = "2.0";

    // prepare top level nodes
    if (!scene.cameras.empty()) js["cameras"] = json::array();
    if (!scene.textures.empty()) {
        js["textures"] = json::array();
        js["images"]   = json::array();
    }
    if (!scene.materials.empty()) js["materials"] = json::array();
    if (!scene.shapes.empty()) {
        js["meshes"]      = json::array();
        js["buffers"]     = json::array();
        js["bufferViews"] = json::array();
        js["accessors"]   = json::array();
    }
    if (!scene.instances.empty()) js["nodes"] = json::array();
    if (!scene.nodes.empty()) js["nodes"] = json::array();

    // convert cameras
    for (auto& camera : scene.cameras) {
        auto cjs    = json();
        cjs["name"] = camera.uri;
        if (!camera.orthographic) {
            cjs["type"]         = "perspective";
            auto& pcjs          = cjs["perspective"];
            pcjs["yfov"]        = get_camera_fovy(camera);
            pcjs["aspectRatio"] = camera.film_width / camera.film_height;
            pcjs["znear"]       = 0.01f;
        } else {
            cjs["type"]   = "orthographic";
            auto& ocjs    = cjs["orthographic"];
            ocjs["xmag"]  = camera.film_width / 2;
            ocjs["ymag"]  = camera.film_height / 2;
            ocjs["znear"] = 0.01f;
        }
        js["cameras"].push_back(cjs);
    }

    // textures
    for (auto& texture : scene.textures) {
        auto tjs = json(), ijs = json();
        tjs["source"] = (int)js["images"].size();
        ijs["uri"]    = texture.uri;
        js["images"].push_back(ijs);
        js["textures"].push_back(tjs);
    }

    // material
    for (auto& material : scene.materials) {
        auto mjs    = json();
        mjs["name"] = material.uri;
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
                mmjs["metallicRoughnessTexture"]["index"] =
                    material.specular_texture;
            mjs["pbrMetallicRoughness"] = mmjs;
        } else {
            auto mmjs                = json();
            mmjs["diffuseFactor"]    = kd;
            mmjs["specularFactor"]   = material.specular;
            mmjs["glossinessFactor"] = 1 - material.roughness;
            if (material.diffuse_texture >= 0)
                mmjs["diffuseTexture"]["index"] = material.diffuse_texture;
            if (material.specular_texture >= 0)
                mmjs["specularGlossinessTexture"]["index"] =
                    material.specular_texture;
            mjs["extensions"]["KHR_materials_pbrSpecularGlossiness"] = mmjs;
        }
        if (material.normal_texture >= 0)
            mjs["normalTexture"]["index"] = material.normal_texture;
        js["materials"].push_back(mjs);
    }
    // shapes
    auto smap = unordered_map<vec2i, int>{};
    for (auto& instance : scene.instances) {
        if (smap.find({instance.shape, instance.material}) != smap.end())
            continue;
        auto& shape = scene.shapes[instance.shape];
        auto  mjs = json(), bjs = json(), pjs = json();
        auto  bid         = js["buffers"].size();
        mjs["name"]       = shape.uri;
        mjs["primitives"] = json::array();
        bjs["name"]       = shape.uri;
        bjs["byteLength"] = 0;
        bjs["uri"]        = get_noextension(shape.uri) + ".bin";
        pjs["material"]   = instance.material;
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
        if (shape.quads_positions.empty()) {
            auto nverts = (int)shape.positions.size();
            if (!shape.positions.empty())
                pjs["attributes"]["POSITION"] = add_accessor(nverts, "VEC3");
            if (!shape.normals.empty())
                pjs["attributes"]["NORMAL"] = add_accessor(nverts, "VEC3");
            if (!shape.texturecoords.empty())
                pjs["attributes"]["TEXCOORD_0"] = add_accessor(nverts, "VEC2");
            if (!shape.colors.empty())
                pjs["attributes"]["COLOR_0"] = add_accessor(nverts, "VEC4");
            if (!shape.radius.empty())
                pjs["attributes"]["RADIUS"] = add_accessor(nverts, "SCALAR");
            if (!shape.points.empty()) {
                pjs["indices"] = add_accessor(
                    (int)shape.points.size(), "SCALAR", true);
                pjs["mode"] = 1;
            }
            if (!shape.lines.empty()) {
                pjs["indices"] = add_accessor(
                    (int)shape.lines.size() * 2, "SCALAR", true);
                pjs["mode"] = 1;
            }
            if (!shape.triangles.empty()) {
                pjs["indices"] = add_accessor(
                    (int)shape.triangles.size() * 3, "SCALAR", true);
                pjs["mode"] = 4;
            }
            if (!shape.quads.empty()) {
                auto triangles = vector<vec3i>{};
                convert_quads_to_triangles(triangles, shape.quads);
                pjs["indices"] = add_accessor(
                    (int)triangles.size() * 3, "SCALAR", true);
                pjs["mode"] = 4;
            }
        } else {
            auto positions     = vector<vec3f>{};
            auto normals       = vector<vec3f>{};
            auto texturecoords = vector<vec2f>{};
            auto quads         = vector<vec4i>{};
            auto triangles     = vector<vec3i>{};
            convert_facevarying(quads, positions, normals, texturecoords,
                shape.quads_positions, shape.quads_normals,
                shape.quads_texturecoords, shape.positions, shape.normals,
                shape.texturecoords);
            convert_quads_to_triangles(triangles, quads);
            auto nverts = (int)positions.size();
            if (!positions.empty())
                pjs["attributes"]["POSITION"] = add_accessor(nverts, "VEC3");
            if (!normals.empty())
                pjs["attributes"]["NORMAL"] = add_accessor(nverts, "VEC3");
            if (!texturecoords.empty())
                pjs["attributes"]["TEXCOORD_0"] = add_accessor(nverts, "VEC2");
            if (!triangles.empty()) {
                pjs["indices"] = add_accessor(
                    (int)triangles.size() * 3, "SCALAR", true);
                pjs["mode"] = 4;
            }
        }
        mjs["primitives"].push_back(pjs);
        js["meshes"].push_back(mjs);
        js["buffers"].push_back(bjs);
        smap[{instance.shape, instance.material}] = (int)js["meshes"].size() -
                                                    1;
    }

    // nodes
    for (auto& node : scene.nodes) {
        auto njs           = json();
        njs["name"]        = node.uri;
        njs["matrix"]      = mat4f(node.local);
        njs["translation"] = node.translation;
        njs["rotation"]    = node.rotation;
        njs["scale"]       = node.scale;
        if (node.camera >= 0) njs["camera"] = node.camera;
        if (node.instance >= 0) {
            auto& instance = scene.instances[node.instance];
            njs["mesh"]    = smap.at({instance.shape, instance.material});
        }
        if (!node.children.empty()) {
            njs["children"] = json::array();
            for (auto& c : node.children) njs["children"].push_back(c);
        }
        js["nodes"].push_back(njs);
    }

    // animations not supported yet
    if (!scene.animations.empty())
        throw io_error("animation not supported yet");

    // nodes from instances
    if (scene.nodes.empty()) {
        auto camera_id = 0;
        for (auto& camera : scene.cameras) {
            auto njs      = json();
            njs["name"]   = camera.uri;
            njs["camera"] = camera_id++;
            njs["matrix"] = mat4f(camera.frame);
            js["nodes"].push_back(njs);
        }
        for (auto& instance : scene.instances) {
            auto njs      = json();
            njs["name"]   = instance.uri;
            njs["mesh"]   = smap.at({instance.shape, instance.material});
            njs["matrix"] = mat4f(instance.frame);
            js["nodes"].push_back(njs);
        }
    }
}

// save gltf mesh
void save_gltf_mesh(const string& filename, const yocto_shape& shape) {
    // open file
    auto fs = output_file(filename, std::ios::binary);

    if (shape.quads_positions.empty()) {
        write_values(fs, shape.positions);
        write_values(fs, shape.normals);
        write_values(fs, shape.texturecoords);
        write_values(fs, shape.colors);
        write_values(fs, shape.radius);
        write_values(fs, shape.points);
        write_values(fs, shape.lines);
        write_values(fs, shape.triangles);
        auto qtriangles = vector<vec3i>{};
        convert_quads_to_triangles(qtriangles, shape.quads);
        write_values(fs, qtriangles);
    } else {
        auto positions     = vector<vec3f>{};
        auto normals       = vector<vec3f>{};
        auto texturecoords = vector<vec2f>{};
        auto quads         = vector<vec4i>{};
        auto triangles     = vector<vec3i>{};
        convert_facevarying(quads, positions, normals, texturecoords,
            shape.quads_positions, shape.quads_normals,
            shape.quads_texturecoords, shape.positions, shape.normals,
            shape.texturecoords);
        convert_quads_to_triangles(triangles, quads);
        write_values(fs, positions);
        write_values(fs, normals);
        write_values(fs, texturecoords);
        write_values(fs, triangles);
    }
}

// Save gltf json
void save_gltf_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options) {
    try {
        // save json
        auto js               = json::object();
        js["asset"]           = json::object();
        js["asset"]["format"] = "Yocto/Scene";
        js["asset"]["generator"] =
            "Yocto/GL - https://github.com/xelatihy/yocto-gl";
        scene_to_gltf(scene, js);
        save_json(filename, js);

        // meshes
        auto dirname = get_dirname(filename);
        for (auto& shape : scene.shapes) {
            if (shape.uri == "") continue;
            save_gltf_mesh(
                get_noextension(dirname + shape.uri) + ".bin", shape);
        }

        // save textures
        save_scene_textures(scene, dirname, options);
    } catch (const std::exception& e) {
        throw io_error("cannot save scene " + filename + "\n" + e.what());
    }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF PBRT
// -----------------------------------------------------------------------------
namespace yocto {

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
void load_pbrt_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options) {
    scene = yocto_scene{};

    struct parse_callbacks : pbrt_callbacks {
        yocto_scene&              scene;
        const load_scene_options& options;
        const string&             filename;

        parse_callbacks(yocto_scene& scene, const load_scene_options& options,
            const string& filename)
            : scene{scene}, options{options}, filename{filename} {}

        bool verbose                 = false;
        bool remove_contant_textures = true;

        unordered_map<string, yocto_material> mmap =
            unordered_map<string, yocto_material>{{"", {}}};
        unordered_map<string, vec3f> amap = unordered_map<string, vec3f>{
            {"", zero3f}};
        unordered_map<string, int> ammap = unordered_map<string, int>{};
        unordered_map<string, int> tmap  = unordered_map<string, int>{{"", -1}};
        unordered_map<string, vec3f> ctmap = unordered_map<string, vec3f>{
            {"", zero3f}};
        unordered_map<string, bool> timap = unordered_map<string, bool>{
            {"", false}};
        unordered_map<string, vector<yocto_instance>> omap =
            unordered_map<string, vector<yocto_instance>>{};
        string cur_object = ""s;

        float last_film_aspect = -1.0f;

        bool is_constant_texture(const string& name) {
            return ctmap.find(name) != ctmap.end();
        }
        vec3f get_constant_texture_color(const string& name) {
            return ctmap.at(name);
        }

        int get_material(const pbrt_context& ctx) {
            static auto light_id    = 0;
            auto        lookup_name = ctx.material + "_______" + ctx.arealight;
            if (ammap.find(lookup_name) != ammap.end())
                return ammap.at(lookup_name);
            auto material     = mmap.at(ctx.material);
            material.emission = amap.at(ctx.arealight);
            if (material.emission != zero3f) {
                material.uri += "_arealight_" + to_string(light_id++);
            }
            scene.materials.push_back(material);
            ammap[lookup_name] = (int)scene.materials.size() - 1;
            return (int)scene.materials.size() - 1;
        }

        void get_scaled_texture3f(const pbrt_textured<spectrum3f>& textured,
            vec3f& value, int& texture) {
            if (textured.texture == "") {
                value = {textured.value.x, textured.value.y, textured.value.z};
                texture = -1;
            } else if (is_constant_texture(textured.texture)) {
                value   = get_constant_texture_color(textured.texture);
                texture = -1;
            } else {
                value   = {1, 1, 1};
                texture = tmap.at(textured.texture);
            }
        }

        float pbrt_remap_roughness(float roughness) {
            // from pbrt code
            roughness = max(roughness, 1e-3f);
            float x   = log(roughness);
            return 1.62142f + 0.819955f * x + 0.1734f * x * x +
                   0.0171201f * x * x * x + 0.000640711f * x * x * x * x;
        }

        void camera(const pbrt_camera& pcamera, const pbrt_context& ctx) {
            auto camera    = yocto_camera{};
            camera.frame   = inverse((frame3f)ctx.transform_start);
            camera.frame.z = -camera.frame.z;
            if (holds_alternative<pbrt_perspective_camera>(pcamera)) {
                auto& perspective = get<pbrt_perspective_camera>(pcamera);
                auto  aspect      = perspective.frameaspectratio;
                if (aspect < 0) aspect = last_film_aspect;
                if (aspect < 0) aspect = 1;
                if (aspect >= 1) {
                    set_camera_perspectivey(camera, radians(perspective.fov),
                        aspect,
                        clamp(perspective.focaldistance, 1.0e-2f, 1.0e4f));
                } else {
                    set_camera_perspectivex(camera, radians(perspective.fov),
                        aspect,
                        clamp(perspective.focaldistance, 1.0e-2f, 1.0e4f));
                }
            } else if (holds_alternative<pbrt_realistic_camera>(pcamera)) {
                auto& realistic     = get<pbrt_realistic_camera>(pcamera);
                camera.focal_length = max(realistic.approx_focallength, 35.0f) *
                                      0.001f;
                auto aspect = 1.0f;
                if (aspect < 0) aspect = last_film_aspect;
                if (aspect < 0) aspect = 1;
                if (aspect >= 1) {
                    camera.film_height = camera.film_width / aspect;
                } else {
                    camera.film_width = camera.film_height * aspect;
                }
                camera.focus_distance = realistic.focusdistance;
                camera.lens_aperture  = realistic.aperturediameter / 2;
            } else {
                throw io_error(
                    "unsupported Camera type " + to_string(pcamera.index()));
            }
            scene.cameras.push_back(camera);
        }
        void film(const pbrt_film& pfilm, const pbrt_context& ctx) {
            if (holds_alternative<pbrt_film_image>(pfilm)) {
                auto& perspective = get<pbrt_film_image>(pfilm);
                last_film_aspect  = (float)perspective.xresolution /
                                   (float)perspective.yresolution;
                for (auto& camera : scene.cameras) {
                    camera.film_width = camera.film_height * last_film_aspect;
                }
            } else {
                throw io_error("unsupported pbrt type");
            }
        }
        void shape(const pbrt_shape& pshape, const pbrt_context& ctx) {
            static auto shape_id = 0;
            auto        shape    = yocto_shape{};
            shape.uri = "shapes/shape__" + to_string(shape_id++) + ".ply";
            if (holds_alternative<pbrt_trianglemesh_shape>(pshape)) {
                auto& mesh          = get<pbrt_trianglemesh_shape>(pshape);
                shape.positions     = mesh.P;
                shape.normals       = mesh.N;
                shape.texturecoords = mesh.uv;
                for (auto& uv : shape.texturecoords) uv.y = (1 - uv.y);
                shape.triangles = mesh.indices;
            } else if (holds_alternative<pbrt_loopsubdiv_shape>(pshape)) {
                auto& mesh      = get<pbrt_loopsubdiv_shape>(pshape);
                shape.positions = mesh.P;
                shape.triangles = mesh.indices;
                shape.normals.resize(shape.positions.size());
                compute_vertex_normals(
                    shape.normals, shape.triangles, shape.positions);
            } else if (holds_alternative<pbrt_plymesh_shape>(pshape)) {
                auto& mesh = get<pbrt_plymesh_shape>(pshape);
                shape.uri  = mesh.filename;
                if (!options.skip_meshes) {
                    load_shape(get_dirname(filename) + mesh.filename,
                        shape.points, shape.lines, shape.triangles, shape.quads,
                        shape.quads_positions, shape.quads_normals,
                        shape.quads_texturecoords, shape.positions,
                        shape.normals, shape.texturecoords, shape.colors,
                        shape.radius, false);
                }
            } else if (holds_alternative<pbrt_sphere_shape>(pshape)) {
                auto& sphere = get<pbrt_sphere_shape>(pshape);
                make_uvsphere_shape(shape.quads, shape.positions, shape.normals,
                    shape.texturecoords, {64, 32}, 2 * sphere.radius, {1, 1},
                    identity_frame3f);
            } else if (holds_alternative<pbrt_disk_shape>(pshape)) {
                auto& disk = get<pbrt_disk_shape>(pshape);
                make_uvdisk_shape(shape.quads, shape.positions, shape.normals,
                    shape.texturecoords, {32, 16}, 2 * disk.radius, {1, 1},
                    identity_frame3f);
            } else {
                throw io_error(
                    "unsupported shape type " + to_string(pshape.index()));
            }
            scene.shapes.push_back(shape);
            auto instance     = yocto_instance{};
            instance.frame    = (frame3f)ctx.transform_start;
            instance.shape    = (int)scene.shapes.size() - 1;
            instance.material = get_material(ctx);
            if (cur_object == "") {
                scene.instances.push_back(instance);
            } else {
                omap[cur_object].push_back(instance);
            }
        }
        void texture(const pbrt_texture& ptexture, const string& name,
            const pbrt_context& ctx) {
            if (remove_contant_textures &&
                holds_alternative<pbrt_constant_texture>(ptexture)) {
                auto& constant = get<pbrt_constant_texture>(ptexture);
                ctmap[name]    = (vec3f)constant.value.value;
                timap[name]    = false;
                return;
            }
            auto texture = yocto_texture{};
            texture.uri  = "textures/" + name + ".png";
            if (holds_alternative<pbrt_imagemap_texture>(ptexture)) {
                auto& imagemap = get<pbrt_imagemap_texture>(ptexture);
                texture.uri    = imagemap.filename;
            } else if (holds_alternative<pbrt_constant_texture>(ptexture)) {
                auto& constant = get<pbrt_constant_texture>(ptexture);
                texture.ldr_image.resize({1, 1});
                texture.ldr_image[{0, 0}] = float_to_byte(
                    vec4f{(vec3f)constant.value.value, 1});
            } else if (holds_alternative<pbrt_bilerp_texture>(ptexture)) {
                // auto& bilerp   = get<pbrt_bilerp_texture>(ptexture);
                texture.ldr_image.resize({1, 1});
                texture.ldr_image[{0, 0}] = {255, 0, 0, 255};
                if (verbose) printf("texture bilerp not supported well");
            } else if (holds_alternative<pbrt_checkerboard_texture>(ptexture)) {
                auto& checkerboard = get<pbrt_checkerboard_texture>(ptexture);
                auto  rgb1         = checkerboard.tex1.texture == ""
                                ? checkerboard.tex1.value
                                : spectrum3f{0.4f, 0.4f, 0.4f};
                auto rgb2 = checkerboard.tex1.texture == ""
                                ? checkerboard.tex2.value
                                : spectrum3f{0.6f, 0.6f, 0.6f};
                make_checker_image(texture.ldr_image, {1024, 1024}, 16,
                    {float_to_byte(vec3f{rgb1.x, rgb1.y, rgb1.z}), 255},
                    {float_to_byte(vec3f{rgb2.x, rgb2.y, rgb2.z}), 255});
                if (verbose) printf("texture checkerboard not supported well");
            } else if (holds_alternative<pbrt_dots_texture>(ptexture)) {
                // auto& dots   = get<pbrt_dots_texture>(ptexture);
                texture.ldr_image.resize({1, 1});
                texture.ldr_image[{0, 0}] = {255, 0, 0, 255};
                if (verbose) printf("texture dots not supported well");
            } else if (holds_alternative<pbrt_fbm_texture>(ptexture)) {
                auto& fbm = get<pbrt_fbm_texture>(ptexture);
                make_fbm_image(texture.ldr_image, {1024, 1024}, {0, 0, 0, 255},
                    {255, 255, 255, 255}, (float)1, (float)2, (float)0.5f,
                    fbm.octaves);
                if (verbose) printf("texture fbm not supported well");
            } else if (holds_alternative<pbrt_marble_texture>(ptexture)) {
                auto& marble = get<pbrt_marble_texture>(ptexture);
                make_fbm_image(texture.ldr_image, {1024, 1024}, {0, 0, 0, 255},
                    {255, 255, 255, 255}, (float)marble.scale, (float)2,
                    (float)0.5f, marble.octaves);
                if (verbose) printf("texture marble not supported well");
            } else if (holds_alternative<pbrt_mix_texture>(ptexture)) {
                auto& mix = get<pbrt_mix_texture>(ptexture);
                if (timap.at(mix.tex1.texture)) {
                    texture.uri =
                        scene.textures.at(tmap.at(mix.tex1.texture)).uri;
                } else if (timap.at(mix.tex2.texture)) {
                    texture.uri =
                        scene.textures.at(tmap.at(mix.tex2.texture)).uri;
                } else {
                    texture.ldr_image.resize({1, 1});
                    texture.ldr_image[{0, 0}] = {255, 0, 0, 255};
                }
                if (verbose) printf("texture mix not supported well");
            } else if (holds_alternative<pbrt_scale_texture>(ptexture)) {
                auto& scale = get<pbrt_scale_texture>(ptexture);
                if (timap.at(scale.tex1.texture)) {
                    texture.uri =
                        scene.textures.at(tmap.at(scale.tex1.texture)).uri;
                } else if (timap.at(scale.tex2.texture)) {
                    texture.uri =
                        scene.textures.at(tmap.at(scale.tex2.texture)).uri;
                } else {
                    texture.ldr_image.resize({1, 1});
                    texture.ldr_image[{0, 0}] = {255, 0, 0, 255};
                }
                if (verbose) printf("texture scale not supported well");
            } else if (holds_alternative<pbrt_uv_texture>(ptexture)) {
                // auto& uv   = get<pbrt_uv_texture>(ptexture);
                texture.ldr_image.resize({1, 1});
                texture.ldr_image[{0, 0}] = {255, 0, 0, 255};
                if (verbose) printf("texture uv not supported well");
            } else if (holds_alternative<pbrt_windy_texture>(ptexture)) {
                // auto& windy   = get<pbrt_uv_texture>(ptexture);
                texture.ldr_image.resize({1, 1});
                texture.ldr_image[{0, 0}] = {255, 0, 0, 255};
                if (verbose) printf("texture windy not supported well");
            } else if (holds_alternative<pbrt_wrinkled_texture>(ptexture)) {
                // auto& uv   = get<pbrt_wrinkled_texture>(ptexture);
                texture.ldr_image.resize({1, 1});
                texture.ldr_image[{0, 0}] = {255, 0, 0, 255};
                if (verbose) printf("texture wrinkled not supported well");
            } else {
                throw io_error(
                    "texture not supported" + to_string(ptexture.index()));
            }
            scene.textures.push_back(texture);
            tmap[name]  = (int)scene.textures.size() - 1;
            timap[name] = holds_alternative<pbrt_imagemap_texture>(ptexture);
        }
        void material(const pbrt_material& pmaterial, const string& name,
            const pbrt_context& ctx) {
            auto material = yocto_material{};
            material.uri  = name;
            if (holds_alternative<pbrt_uber_material>(pmaterial)) {
                auto& uber = get<pbrt_uber_material>(pmaterial);
                get_scaled_texture3f(
                    uber.Kd, material.diffuse, material.diffuse_texture);
                get_scaled_texture3f(
                    uber.Ks, material.specular, material.specular_texture);
                get_scaled_texture3f(uber.Kt, material.transmission,
                    material.transmission_texture);
                auto op     = vec3f{0, 0, 0};
                auto op_txt = -1;
                get_scaled_texture3f(uber.opacity, op, op_txt);
                material.opacity = (op.x + op.y + op.z) / 3;
                material.roughness =
                    (uber.uroughness.value + uber.vroughness.value) / 2;
                if (uber.remaproughness)
                    material.roughness = pbrt_remap_roughness(
                        material.roughness);
            } else if (holds_alternative<pbrt_plastic_material>(pmaterial)) {
                auto& plastic = get<pbrt_plastic_material>(pmaterial);
                get_scaled_texture3f(
                    plastic.Kd, material.diffuse, material.diffuse_texture);
                // get_scaled_texture3f(
                //     plastic.Ks, material.specular,
                //     material.specular_texture);
                material.specular = {0.04f, 0.04f, 0.04f};
                material.roughness =
                    (plastic.uroughness.value + plastic.vroughness.value) / 2;
                if (plastic.remaproughness)
                    material.roughness = pbrt_remap_roughness(
                        material.roughness);
            } else if (holds_alternative<pbrt_translucent_material>(
                           pmaterial)) {
                auto& translucent = get<pbrt_translucent_material>(pmaterial);
                get_scaled_texture3f(
                    translucent.Kd, material.diffuse, material.diffuse_texture);
                // get_scaled_texture3f(translucent.Ks, material.specular,
                //     material.specular_texture);
                material.specular  = {0.04f, 0.04f, 0.04f};
                material.roughness = (translucent.uroughness.value +
                                         translucent.vroughness.value) /
                                     2;
                if (translucent.remaproughness)
                    material.roughness = pbrt_remap_roughness(
                        material.roughness);
            } else if (holds_alternative<pbrt_matte_material>(pmaterial)) {
                auto& matte = get<pbrt_matte_material>(pmaterial);
                get_scaled_texture3f(
                    matte.Kd, material.diffuse, material.diffuse_texture);
                material.roughness = 1;
            } else if (holds_alternative<pbrt_mirror_material>(pmaterial)) {
                // auto& mirror          =
                // get<pbrt_mirror_material>(pmaterial);
                material.diffuse   = {0, 0, 0};
                material.specular  = {1, 1, 1};
                material.roughness = 0;
            } else if (holds_alternative<pbrt_metal_material>(pmaterial)) {
                auto& metal = get<pbrt_metal_material>(pmaterial);
                auto  eta = zero3f, k = zero3f;
                auto  eta_texture = -1, k_texture = -1;
                get_scaled_texture3f(metal.eta, eta, eta_texture);
                get_scaled_texture3f(metal.k, k, k_texture);
                material.specular = pbrt_fresnel_metal(1, eta, k);
                material.roughness =
                    (metal.uroughness.value + metal.vroughness.value) / 2;
                if (metal.remaproughness)
                    material.roughness = pbrt_remap_roughness(
                        material.roughness);
            } else if (holds_alternative<pbrt_substrate_material>(pmaterial)) {
                auto& substrate = get<pbrt_substrate_material>(pmaterial);
                get_scaled_texture3f(
                    substrate.Kd, material.diffuse, material.diffuse_texture);
                // get_scaled_texture3f(
                //     substrate.Ks, material.specular,
                //     material.specular_texture);
                material.specular = {0.04f, 0.04f, 0.04f};
                material.roughness =
                    (substrate.uroughness.value + substrate.vroughness.value) /
                    2;
                if (substrate.remaproughness)
                    material.roughness = pbrt_remap_roughness(
                        material.roughness);
            } else if (holds_alternative<pbrt_glass_material>(pmaterial)) {
                // auto& glass       = get<pbrt_glass_material>(pmaterial);
                material.specular     = {0.04f, 0.04f, 0.04f};
                material.transmission = {1, 1, 1};
                // get_scaled_texture3f(
                //     glass.Kr, material.specular,
                //     material.specular_texture);
                // get_scaled_texture3f(
                //     glass.Kt, material.transmission,
                //     material.transmission_texture);
                material.roughness = 0;
            } else if (holds_alternative<pbrt_hair_material>(pmaterial)) {
                auto& hair = get<pbrt_hair_material>(pmaterial);
                get_scaled_texture3f(
                    hair.color, material.diffuse, material.diffuse_texture);
                material.roughness = 1;
                if (verbose) printf("hair material not properly supported\n");
            } else if (holds_alternative<pbrt_disney_material>(pmaterial)) {
                auto& disney = get<pbrt_disney_material>(pmaterial);
                get_scaled_texture3f(
                    disney.color, material.diffuse, material.diffuse_texture);
                material.roughness = 1;
                if (verbose) printf("disney material not properly supported\n");
            } else if (holds_alternative<pbrt_kdsubsurface_material>(
                           pmaterial)) {
                auto& kdsubdurface = get<pbrt_kdsubsurface_material>(pmaterial);
                get_scaled_texture3f(kdsubdurface.Kd, material.diffuse,
                    material.diffuse_texture);
                // get_scaled_texture3f(kdsubdurface.Kr, material.specular,
                //     material.specular_texture);
                material.specular  = {0.04f, 0.04f, 0.04f};
                material.roughness = (kdsubdurface.uroughness.value +
                                         kdsubdurface.vroughness.value) /
                                     2;
                if (kdsubdurface.remaproughness)
                    material.roughness = pbrt_remap_roughness(
                        material.roughness);
                if (verbose)
                    printf("kdsubsurface material not properly supported\n");
            } else if (holds_alternative<pbrt_subsurface_material>(pmaterial)) {
                // auto& subdurface           =
                // get<pbrt_subsurface_material>(pmaterial);
                material.diffuse   = {1, 0, 0};
                material.roughness = 1;
                if (verbose)
                    printf("subsurface material not properly supported\n");
            } else if (holds_alternative<pbrt_mix_material>(pmaterial)) {
                auto& mix     = get<pbrt_mix_material>(pmaterial);
                auto  matname = (!mix.namedmaterial1.empty())
                                   ? mix.namedmaterial1
                                   : mix.namedmaterial2;
                material = mmap.at(matname);
                if (verbose) printf("mix material not properly supported\n");
            } else if (holds_alternative<pbrt_fourier_material>(pmaterial)) {
                auto& fourier = get<pbrt_fourier_material>(pmaterial);
                if (holds_alternative<pbrt_plastic_material>(fourier.approx)) {
                    auto& plastic = get<pbrt_plastic_material>(fourier.approx);
                    get_scaled_texture3f(
                        plastic.Kd, material.diffuse, material.diffuse_texture);
                    // get_scaled_texture3f(plastic.Ks, material.specular,
                    //     material.specular_texture);
                    material.specular = {0.04f, 0.04f, 0.04f};
                    material.roughness =
                        (plastic.uroughness.value + plastic.vroughness.value) /
                        2;
                    if (plastic.remaproughness)
                        material.roughness = pbrt_remap_roughness(
                            material.roughness);
                } else if (holds_alternative<pbrt_metal_material>(
                               fourier.approx)) {
                    auto& metal = get<pbrt_metal_material>(fourier.approx);
                    auto  eta = zero3f, k = zero3f;
                    auto  eta_texture = -1, k_texture = -1;
                    get_scaled_texture3f(metal.eta, eta, eta_texture);
                    get_scaled_texture3f(metal.k, k, k_texture);
                    material.specular = pbrt_fresnel_metal(1, eta, k);
                    material.roughness =
                        (metal.uroughness.value + metal.vroughness.value) / 2;
                    if (metal.remaproughness)
                        material.roughness = pbrt_remap_roughness(
                            material.roughness);
                } else if (holds_alternative<pbrt_glass_material>(
                               fourier.approx)) {
                    // auto& glass = get<pbrt_metal_glass>(fourier.approx);
                    material.diffuse      = {0, 0, 0};
                    material.specular     = {0.04f, 0.04f, 0.04f};
                    material.transmission = {1, 1, 1};
                    material.roughness    = 0;
                } else {
                    throw io_error("material type not supported " +
                                   to_string(fourier.approx.index()));
                }
            } else {
                throw io_error("material type not supported " +
                               to_string(pmaterial.index()));
            }
            mmap[name] = material;
        }
        void arealight(const pbrt_arealight& plight, const string& name,
            const pbrt_context& ctx) {
            auto emission = zero3f;
            if (holds_alternative<pbrt_diffuse_arealight>(plight)) {
                auto& diffuse = get<pbrt_diffuse_arealight>(plight);
                emission      = (vec3f)diffuse.L * (vec3f)diffuse.scale;
            } else {
                throw io_error("area light type not supported " +
                               to_string(plight.index()));
            }
            amap[name] = emission;
        }
        void light(const pbrt_light& plight, const pbrt_context& ctx) {
            static auto light_id = 0;
            auto        name     = "light_" + to_string(light_id++);
            if (holds_alternative<pbrt_infinite_light>(plight)) {
                auto& infinite    = get<pbrt_infinite_light>(plight);
                auto  environment = yocto_environment();
                environment.uri   = name;
                // environment.frame =
                // frame3f{{1,0,0},{0,0,-1},{0,-1,0},{0,0,0}}
                // * stack.back().frame;
                environment.frame =
                    (frame3f)ctx.transform_start *
                    frame3f{{1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 0, 0}};
                environment.emission = (vec3f)infinite.scale;
                if (infinite.mapname != "") {
                    auto texture = yocto_texture{};
                    texture.uri  = infinite.mapname;
                    scene.textures.push_back(texture);
                    environment.emission_texture = (int)scene.textures.size() -
                                                   1;
                }
                scene.environments.push_back(environment);
            } else if (holds_alternative<pbrt_distant_light>(plight)) {
                auto& distant      = get<pbrt_distant_light>(plight);
                auto  distant_dist = 100;
                scene.shapes.push_back({});
                auto& shape = scene.shapes.back();
                shape.uri   = name;
                auto dir    = normalize(distant.from - distant.to);
                auto size   = distant_dist * sin(5 * pif / 180);
                make_quad_shape(shape.quads, shape.positions, shape.normals,
                    shape.texturecoords, {1, 1}, {size, size}, {1, 1},
                    identity_frame3f);
                scene.materials.push_back({});
                auto& material    = scene.materials.back();
                material.uri      = shape.uri;
                material.emission = (vec3f)distant.L * (vec3f)distant.scale;
                material.emission *= (distant_dist * distant_dist) /
                                     (size * size);
                auto instance     = yocto_instance();
                instance.uri      = shape.uri;
                instance.shape    = (int)scene.shapes.size() - 1;
                instance.material = (int)scene.materials.size() - 1;
                instance.frame    = (frame3f)ctx.transform_start *
                                 make_lookat_frame(dir * distant_dist, zero3f,
                                     {0, 1, 0}, true);
                scene.instances.push_back(instance);
            } else if (holds_alternative<pbrt_point_light>(plight)) {
                auto& point = get<pbrt_point_light>(plight);
                scene.shapes.push_back({});
                auto& shape = scene.shapes.back();
                shape.uri   = name;
                auto size   = 0.01f;
                make_sphere_shape(shape.quads, shape.positions, shape.normals,
                    shape.texturecoords, 4.0f, size, 1.0f, identity_frame3f);
                scene.materials.push_back({});
                auto& material    = scene.materials.back();
                material.uri      = shape.uri;
                material.emission = (vec3f)point.I * (vec3f)point.scale;
                // TODO: fix emission
                auto instance     = yocto_instance();
                instance.uri      = shape.uri;
                instance.shape    = (int)scene.shapes.size() - 1;
                instance.material = (int)scene.materials.size() - 1;
                instance.frame    = (frame3f)ctx.transform_start *
                                 make_translation_frame(point.from);
                scene.instances.push_back(instance);
            } else if (holds_alternative<pbrt_goniometric_light>(plight)) {
                auto& goniometric = get<pbrt_goniometric_light>(plight);
                scene.shapes.push_back({});
                auto& shape = scene.shapes.back();
                shape.uri   = name;
                auto size   = 0.01f;
                make_sphere_shape(shape.quads, shape.positions, shape.normals,
                    shape.texturecoords, 4.0f, size, 1.0f, identity_frame3f);
                scene.materials.push_back({});
                auto& material    = scene.materials.back();
                material.uri      = shape.uri;
                material.emission = (vec3f)goniometric.I *
                                    (vec3f)goniometric.scale;
                // TODO: fix emission
                auto instance     = yocto_instance();
                instance.uri      = shape.uri;
                instance.shape    = (int)scene.shapes.size() - 1;
                instance.material = (int)scene.materials.size() - 1;
                instance.frame    = (frame3f)ctx.transform_start;
                scene.instances.push_back(instance);
            } else if (holds_alternative<pbrt_spot_light>(plight)) {
                auto& spot = get<pbrt_spot_light>(plight);
                scene.shapes.push_back({});
                auto& shape = scene.shapes.back();
                shape.uri   = name;
                auto size   = 0.01f;
                make_sphere_shape(shape.quads, shape.positions, shape.normals,
                    shape.texturecoords, 4.0f, size, 1.0f, identity_frame3f);
                scene.materials.push_back({});
                auto& material    = scene.materials.back();
                material.uri      = shape.uri;
                material.emission = (vec3f)spot.I * (vec3f)spot.scale;
                // TODO: fix emission
                auto instance     = yocto_instance();
                instance.uri      = shape.uri;
                instance.shape    = (int)scene.shapes.size() - 1;
                instance.material = (int)scene.materials.size() - 1;
                instance.frame    = (frame3f)ctx.transform_start;
                scene.instances.push_back(instance);
            } else {
                throw io_error(
                    "light type not supported " + to_string(plight.index()));
            }
        }
        void begin_object(const pbrt_object& pobject, const pbrt_context& ctx) {
            cur_object       = pobject.name;
            omap[cur_object] = {};
        }
        void end_object(const pbrt_object& pobject, const pbrt_context& ctx) {
            cur_object = "";
        }
        void object_instance(
            const pbrt_object& pobject, const pbrt_context& ctx) {
            auto& pinstances = omap.at(pobject.name);
            for (auto& pinstance : pinstances) {
                auto instance  = yocto_instance();
                instance.frame = (frame3f)ctx.transform_start * pinstance.frame;
                instance.shape = pinstance.shape;
                instance.material = pinstance.material;
                scene.instances.push_back(instance);
            }
        }
    };

    try {
        // Parse pbrt
        auto pbrt_options = load_pbrt_options();
        auto cb           = parse_callbacks{scene, options, filename};
        load_pbrt(filename, cb, pbrt_options);

        // load textures
        auto dirname = get_dirname(filename);
        load_scene_textures(scene, dirname, options);
    } catch (const std::exception& e) {
        throw io_error("cannot load scene " + filename + "\n" + e.what());
    }

    // fix scene
    scene.uri = get_filename(filename);
    add_missing_cameras(scene);
    add_missing_materials(scene);
    normalize_uris(scene);
    trim_memory(scene);
    update_transforms(scene);
}

// Convert a scene to pbrt format
void save_pbrt(const string& filename, const yocto_scene& scene) {
    auto fs = output_file(filename);

    // convert camera and settings
    auto& camera     = scene.cameras.front();
    auto  from       = camera.frame.o;
    auto  to         = camera.frame.o - camera.frame.z;
    auto  up         = camera.frame.y;
    auto  image_size = get_camera_image_size(camera, {0, 720});
    println_values(fs, "LookAt", from, to, up);
    println_values(fs, "Camera \"perspective\" \"float fov\"",
        get_camera_fovy(camera) * 180 / pif);

    // save renderer
    println_values(fs, "Sampler \"random\" \"integer pixelsamples\" [64]");
    // fprintf(f, "Sampler \"sobol\" \"interger pixelsamples\" [64]\n");
    println_values(fs, "Integrator \"path\"");
    println_values(fs, "Film \"image\" \"string filename\" [\"",
        get_noextension(filename) + ".exr", "\"] \"integer xresolution\" [",
        image_size.x, "] \"integer yresolution\" [", image_size.y, "]");

    // start world
    println_values(fs, "WorldBegin");

    // convert textures
    for (auto& texture : scene.textures) {
        println_values(fs, "Texture \"" + get_basename(texture.uri) +
                               "\" \"spectrum\" \"imagemap\" "
                               "\"string filename\" [\"" +
                               texture.uri + "\"]");
    }

    // convert materials
    for (auto& material : scene.materials) {
        println_values(
            fs, "MakeNamedMaterial \"" + get_basename(material.uri) + "\" ");
        println_values(fs, "    \"string type\" \"uber\" ");
        if (material.diffuse_texture >= 0)
            println_values(fs,
                "    \"texture Kd\" [\"" +
                    get_basename(scene.textures[material.diffuse_texture].uri) +
                    "\"] ");
        else
            println_values(fs, "    \"rgb Kd\" [", material.diffuse, "] ");
        if (material.specular_texture >= 0)
            println_values(
                fs, "    \"texture Ks\" [\"" +
                        get_basename(
                            scene.textures[material.specular_texture].uri) +
                        "\"] ");
        else
            println_values(fs, "    \"rgb Ks\" [", material.specular, "] ");
        println_values(
            fs, "    \"float roughness\" [", material.roughness, "] ");
    }

    // convert instances
    for (auto& instance : scene.instances) {
        auto& shape    = scene.shapes[instance.shape];
        auto& material = scene.materials[instance.material];
        println_values(fs, "AttributeBegin");
        println_values(fs, "  TransformBegin");
        println_values(fs, "    ConcatTransform [", mat4f(instance.frame), "]");
        println_values(
            fs, "    NamedMaterial \"" + get_basename(material.uri) + "\"");
        if (material.emission != zero3f)
            println_values(fs, "    AreaLightSource \"diffuse\" \"rgb L\" [ ",
                material.emission, " ]");
        println_values(fs, "    Shape \"plymesh\" \"string filename\" [\"" +
                               get_noextension(shape.uri) + ".ply" + "\"]");
        println_values(fs, "  TransformEnd");
        println_values(fs, "AttributeEnd");
    }

    // end world
    println_values(fs, "WorldEnd");
}

// Save a pbrt scene
void save_pbrt_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options) {
    try {
        // save json
        save_pbrt(filename, scene);

        // save meshes
        auto dirname = get_dirname(filename);
        for (auto& shape : scene.shapes) {
            save_shape(get_noextension(dirname + shape.uri) + ".ply",
                shape.points, shape.lines, shape.triangles, shape.quads,
                shape.quads_positions, shape.quads_normals,
                shape.quads_texturecoords, shape.positions, shape.normals,
                shape.texturecoords, shape.colors, shape.radius);
        }

        // skip textures
        save_scene_textures(scene, dirname, options);

    } catch (const std::exception& e) {
        throw io_error("cannot save scene " + filename + "\n" + e.what());
    }
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

// Serialize vector
template <typename T>
void write_value(output_file& fs, const vector<T>& vec) {
    auto count = (size_t)vec.size();
    write_value(fs, count);
    write_values(fs, vec);
}
template <typename T>
void read_value(input_file& fs, vector<T>& vec) {
    auto count = (size_t)0;
    read_value(fs, count);
    vec = vector<T>(count);
    read_values(fs, vec);
}

// Serialize string
void write_value(output_file& fs, const string& str) {
    auto count = (size_t)str.size();
    write_value(fs, count);
    auto vec = vector<char>(str.begin(), str.end());
    write_values(fs, vec);
}
void read_value(input_file& fs, string& str) {
    auto count = (size_t)0;
    read_value(fs, count);
    auto vec = vector<char>(count);
    read_values(fs, vec);
    str = {vec.begin(), vec.end()};
}

// Serialize image
template <typename T>
void write_value(output_file& fs, const image<T>& img) {
    write_value(fs, img.size());
    write_values(fs, img.data(), (size_t)img.size().x * (size_t)img.size().y);
}
template <typename T>
void read_value(input_file& fs, image<T>& img) {
    auto size = zero2i;
    read_value(fs, size);
    img = {size};
    read_values(fs, img.data(), size.x * size.y);
}

// Serialize image
template <typename T>
void write_value(output_file& fs, const volume<T>& vol) {
    write_value(fs, vol.size());
    write_values(fs, vol.data(),
        (size_t)vol.size().x * (size_t)vol.size().y * (size_t)vol.size().z);
}
template <typename T>
void read_value(input_file& fs, volume<T>& vol) {
    auto size = zero3i;
    read_value(fs, size.x);
    read_value(fs, size.y);
    read_value(fs, size.z);
    vol = {size};
    read_values(fs, vol.data(), size.x * size.y * size.z);
}

// Serialize vector of pointers
template <typename T>
void write_objects(output_file& fs, const vector<T>& vec) {
    auto count = (size_t)vec.size();
    write_value(fs, count);
    for (auto i = 0; i < vec.size(); ++i) {
        write_object(fs, vec[i]);
    }
}
template <typename T>
void read_objects(input_file& fs, vector<T>& vec) {
    auto count = (size_t)0;
    read_value(fs, count);
    vec = vector<T>(count);
    for (auto i = 0; i < vec.size(); ++i) {
        vec[i] = T{};
        read_object(fs, vec[i]);
    }
}

// Serialize yocto types. This is mostly boiler plate code.
void write_object(output_file& fs, const yocto_camera& camera) {
    write_value(fs, camera.uri);
    write_value(fs, camera.frame);
    write_value(fs, camera.orthographic);
    write_value(fs, camera.film_width);
    write_value(fs, camera.film_height);
    write_value(fs, camera.focal_length);
    write_value(fs, camera.focus_distance);
    write_value(fs, camera.lens_aperture);
}
void read_object(input_file& fs, yocto_camera& camera) {
    read_value(fs, camera.uri);
    read_value(fs, camera.frame);
    read_value(fs, camera.orthographic);
    read_value(fs, camera.film_width);
    read_value(fs, camera.film_height);
    read_value(fs, camera.focal_length);
    read_value(fs, camera.focus_distance);
    read_value(fs, camera.lens_aperture);
}

void write_object(output_file& fs, const yocto_shape& shape) {
    write_value(fs, shape.uri);
    write_value(fs, shape.subdivision_level);
    write_value(fs, shape.catmull_clark);
    write_value(fs, shape.compute_normals);
    write_value(fs, shape.preserve_facevarying);
    write_value(fs, shape.points);
    write_value(fs, shape.lines);
    write_value(fs, shape.triangles);
    write_value(fs, shape.quads);
    write_value(fs, shape.quads_positions);
    write_value(fs, shape.quads_normals);
    write_value(fs, shape.quads_texturecoords);
    write_value(fs, shape.positions);
    write_value(fs, shape.normals);
    write_value(fs, shape.texturecoords);
    write_value(fs, shape.colors);
    write_value(fs, shape.radius);
    write_value(fs, shape.tangentspaces);
}
void read_object(input_file& fs, yocto_shape& shape) {
    read_value(fs, shape.uri);
    read_value(fs, shape.subdivision_level);
    read_value(fs, shape.catmull_clark);
    read_value(fs, shape.compute_normals);
    read_value(fs, shape.preserve_facevarying);
    read_value(fs, shape.points);
    read_value(fs, shape.lines);
    read_value(fs, shape.triangles);
    read_value(fs, shape.quads);
    read_value(fs, shape.quads_positions);
    read_value(fs, shape.quads_normals);
    read_value(fs, shape.quads_texturecoords);
    read_value(fs, shape.positions);
    read_value(fs, shape.normals);
    read_value(fs, shape.texturecoords);
    read_value(fs, shape.colors);
    read_value(fs, shape.radius);
    read_value(fs, shape.tangentspaces);
}

void write_object(output_file& fs, const yocto_texture& texture) {
    write_value(fs, texture.uri);
    write_value(fs, texture.hdr_image);
    write_value(fs, texture.ldr_image);
}
void read_object(input_file& fs, yocto_texture& texture) {
    read_value(fs, texture.uri);
    read_value(fs, texture.hdr_image);
    read_value(fs, texture.ldr_image);
}

void write_object(output_file& fs, const yocto_voltexture& texture) {
    write_value(fs, texture.uri);
    write_value(fs, texture.volume_data);
}
void read_object(input_file& fs, yocto_voltexture& texture) {
    read_value(fs, texture.uri);
    read_value(fs, texture.volume_data);
}

void write_object(output_file& fs, const yocto_environment& environment) {
    write_value(fs, environment.uri);
    write_value(fs, environment.frame);
    write_value(fs, environment.emission);
    write_value(fs, environment.emission_texture);
}
void read_object(input_file& fs, yocto_environment& environment) {
    read_value(fs, environment.uri);
    read_value(fs, environment.frame);
    read_value(fs, environment.emission);
    read_value(fs, environment.emission_texture);
}

void write_object(output_file& fs, const yocto_material& material) {
    write_value(fs, material.uri);
    write_value(fs, material.base_metallic);
    write_value(fs, material.gltf_textures);
    write_value(fs, material.emission);
    write_value(fs, material.diffuse);
    write_value(fs, material.specular);
    write_value(fs, material.transmission);
    write_value(fs, material.roughness);
    write_value(fs, material.opacity);
    write_value(fs, material.fresnel);
    write_value(fs, material.refract);
    write_value(fs, material.emission_texture);
    write_value(fs, material.diffuse_texture);
    write_value(fs, material.specular_texture);
    write_value(fs, material.transmission_texture);
    write_value(fs, material.roughness_texture);
    write_value(fs, material.displacement_texture);
    write_value(fs, material.normal_texture);
    write_value(fs, material.volume_emission);
    write_value(fs, material.volume_albedo);
    write_value(fs, material.volume_density);
    write_value(fs, material.volume_phaseg);
    write_value(fs, material.volume_density_texture);
    write_value(fs, material.displacement_scale);
};
void read_object(input_file& fs, yocto_material& material) {
    read_value(fs, material.uri);
    read_value(fs, material.base_metallic);
    read_value(fs, material.gltf_textures);
    read_value(fs, material.emission);
    read_value(fs, material.diffuse);
    read_value(fs, material.specular);
    read_value(fs, material.transmission);
    read_value(fs, material.roughness);
    read_value(fs, material.opacity);
    read_value(fs, material.fresnel);
    read_value(fs, material.refract);
    read_value(fs, material.emission_texture);
    read_value(fs, material.diffuse_texture);
    read_value(fs, material.specular_texture);
    read_value(fs, material.transmission_texture);
    read_value(fs, material.roughness_texture);
    read_value(fs, material.displacement_texture);
    read_value(fs, material.normal_texture);
    read_value(fs, material.volume_emission);
    read_value(fs, material.volume_albedo);
    read_value(fs, material.volume_density);
    read_value(fs, material.volume_phaseg);
    read_value(fs, material.volume_density_texture);
    read_value(fs, material.displacement_scale);
};

void write_object(output_file& fs, const yocto_instance& instance) {
    write_value(fs, instance.uri);
    write_value(fs, instance.frame);
    write_value(fs, instance.shape);
    write_value(fs, instance.material);
};
void read_object(input_file& fs, yocto_instance& instance) {
    read_value(fs, instance.uri);
    read_value(fs, instance.frame);
    read_value(fs, instance.shape);
    read_value(fs, instance.material);
};

void write_object(output_file& fs, const yocto_scene& scene) {
    write_value(fs, scene.uri);
    write_objects(fs, scene.cameras);
    write_objects(fs, scene.shapes);
    write_objects(fs, scene.textures);
    write_objects(fs, scene.voltextures);
    write_objects(fs, scene.materials);
    write_objects(fs, scene.instances);
    write_objects(fs, scene.environments);
}
void read_object(input_file& fs, yocto_scene& scene) {
    read_value(fs, scene.uri);
    read_objects(fs, scene.cameras);
    read_objects(fs, scene.shapes);
    read_objects(fs, scene.textures);
    read_objects(fs, scene.voltextures);
    read_objects(fs, scene.materials);
    read_objects(fs, scene.instances);
    read_objects(fs, scene.environments);
}

// Load/save a binary dump useful for very fast scene IO.
void load_ybin_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options) {
    scene = {};
    try {
        auto fs = input_file(filename, true);
        read_object(fs, scene);
    } catch (const std::exception& e) {
        throw io_error("cannot load scene " + filename + "\n" + e.what());
    }
}

// Load/save a binary dump useful for very fast scene IO.
void save_ybin_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options) {
    try {
        auto fs = output_file(filename, true);
        write_object(fs, scene);
    } catch (const std::exception& e) {
        throw io_error("cannot save scene " + filename + "\n" + e.what());
    }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace yocto {

void make_cornellbox_scene(yocto_scene& scene) {
    scene.uri              = "cornellbox";
    auto& camera           = scene.cameras.emplace_back();
    camera.uri             = "cam";
    camera.frame           = frame3f{{0, 1, 3.9}};
    camera.focal_length    = 0.035;
    camera.lens_aperture   = 0.0;
    camera.film_width      = 0.024;
    camera.film_height     = 0.024;
    auto& floor_mat        = scene.materials.emplace_back();
    floor_mat.uri          = "floor";
    floor_mat.diffuse      = {0.725, 0.71, 0.68};
    auto& ceiling_mat      = scene.materials.emplace_back();
    ceiling_mat.uri        = "ceiling";
    ceiling_mat.diffuse    = {0.725, 0.71, 0.68};
    auto& backwall_mat     = scene.materials.emplace_back();
    backwall_mat.uri       = "backwall";
    backwall_mat.diffuse   = {0.725, 0.71, 0.68};
    auto& rightwall_mat    = scene.materials.emplace_back();
    rightwall_mat.uri      = "rightwall";
    rightwall_mat.diffuse  = {0.14, 0.45, 0.091};
    auto& leftwall_mat     = scene.materials.emplace_back();
    leftwall_mat.uri       = "leftwall";
    leftwall_mat.diffuse   = {0.63, 0.065, 0.05};
    auto& shortbox_mat     = scene.materials.emplace_back();
    shortbox_mat.uri       = "shortbox";
    shortbox_mat.diffuse   = {0.725, 0.71, 0.68};
    auto& tallbox_mat      = scene.materials.emplace_back();
    tallbox_mat.uri        = "tallbox";
    tallbox_mat.diffuse    = {0.725, 0.71, 0.68};
    auto& light_mat        = scene.materials.emplace_back();
    light_mat.uri          = "light";
    light_mat.emission     = {17, 12, 4};
    auto& floor_shp        = scene.shapes.emplace_back();
    floor_shp.uri          = "floor";
    floor_shp.positions    = {{-1, 0, 1}, {1, 0, 1}, {1, 0, -1}, {-1, 0, -1}};
    floor_shp.triangles    = {{0, 1, 2}, {2, 3, 0}};
    auto& ceiling_shp      = scene.shapes.emplace_back();
    ceiling_shp.uri        = "ceiling";
    ceiling_shp.positions  = {{-1, 2, 1}, {-1, 2, -1}, {1, 2, -1}, {1, 2, 1}};
    ceiling_shp.triangles  = {{0, 1, 2}, {2, 3, 0}};
    auto& backwall_shp     = scene.shapes.emplace_back();
    backwall_shp.uri       = "backwall";
    backwall_shp.positions = {{-1, 0, -1}, {1, 0, -1}, {1, 2, -1}, {-1, 2, -1}};
    backwall_shp.triangles = {{0, 1, 2}, {2, 3, 0}};
    auto& rightwall_shp    = scene.shapes.emplace_back();
    rightwall_shp.uri      = "rightwall";
    rightwall_shp.positions = {{1, 0, -1}, {1, 0, 1}, {1, 2, 1}, {1, 2, -1}};
    rightwall_shp.triangles = {{0, 1, 2}, {2, 3, 0}};
    auto& leftwall_shp      = scene.shapes.emplace_back();
    leftwall_shp.uri        = "leftwall";
    leftwall_shp.positions = {{-1, 0, 1}, {-1, 0, -1}, {-1, 2, -1}, {-1, 2, 1}};
    leftwall_shp.triangles = {{0, 1, 2}, {2, 3, 0}};
    auto& shortbox_shp     = scene.shapes.emplace_back();
    shortbox_shp.uri       = "shortbox";
    shortbox_shp.positions = {{0.53, 0.6, 0.75}, {0.7, 0.6, 0.17},
        {0.13, 0.6, 0.0}, {-0.05, 0.6, 0.57}, {-0.05, 0.0, 0.57},
        {-0.05, 0.6, 0.57}, {0.13, 0.6, 0.0}, {0.13, 0.0, 0.0},
        {0.53, 0.0, 0.75}, {0.53, 0.6, 0.75}, {-0.05, 0.6, 0.57},
        {-0.05, 0.0, 0.57}, {0.7, 0.0, 0.17}, {0.7, 0.6, 0.17},
        {0.53, 0.6, 0.75}, {0.53, 0.0, 0.75}, {0.13, 0.0, 0.0},
        {0.13, 0.6, 0.0}, {0.7, 0.6, 0.17}, {0.7, 0.0, 0.17}, {0.53, 0.0, 0.75},
        {0.7, 0.0, 0.17}, {0.13, 0.0, 0.0}, {-0.05, 0.0, 0.57}};
    shortbox_shp.triangles = {{0, 1, 2}, {2, 3, 0}, {4, 5, 6}, {6, 7, 4},
        {8, 9, 10}, {10, 11, 8}, {12, 13, 14}, {14, 15, 12}, {16, 17, 18},
        {18, 19, 16}, {20, 21, 22}, {22, 23, 20}};
    auto& tallbox_shp      = scene.shapes.emplace_back();
    tallbox_shp.uri        = "tallbox";
    tallbox_shp.positions  = {{-0.53, 1.2, 0.09}, {0.04, 1.2, -0.09},
        {-0.14, 1.2, -0.67}, {-0.71, 1.2, -0.49}, {-0.53, 0.0, 0.09},
        {-0.53, 1.2, 0.09}, {-0.71, 1.2, -0.49}, {-0.71, 0.0, -0.49},
        {-0.71, 0.0, -0.49}, {-0.71, 1.2, -0.49}, {-0.14, 1.2, -0.67},
        {-0.14, 0.0, -0.67}, {-0.14, 0.0, -0.67}, {-0.14, 1.2, -0.67},
        {0.04, 1.2, -0.09}, {0.04, 0.0, -0.09}, {0.04, 0.0, -0.09},
        {0.04, 1.2, -0.09}, {-0.53, 1.2, 0.09}, {-0.53, 0.0, 0.09},
        {-0.53, 0.0, 0.09}, {0.04, 0.0, -0.09}, {-0.14, 0.0, -0.67},
        {-0.71, 0.0, -0.49}};
    tallbox_shp.triangles  = {{0, 1, 2}, {2, 3, 0}, {4, 5, 6}, {6, 7, 4},
        {8, 9, 10}, {10, 11, 8}, {12, 13, 14}, {14, 15, 12}, {16, 17, 18},
        {18, 19, 16}, {20, 21, 22}, {22, 23, 20}};
    auto& light_shp        = scene.shapes.emplace_back();
    light_shp.uri          = "light";
    light_shp.positions    = {{-0.25, 1.99, 0.25}, {-0.25, 1.99, -0.25},
        {0.25, 1.99, -0.25}, {0.25, 1.99, 0.25}};
    light_shp.triangles    = {{0, 1, 2}, {2, 3, 0}};
    scene.instances.push_back({"floor", identity_frame3f, 0, 0});
    scene.instances.push_back({"ceiling", identity_frame3f, 1, 1});
    scene.instances.push_back({"backwall", identity_frame3f, 2, 2});
    scene.instances.push_back({"rightwall", identity_frame3f, 3, 3});
    scene.instances.push_back({"leftwall", identity_frame3f, 4, 4});
    scene.instances.push_back({"shortbox", identity_frame3f, 5, 5});
    scene.instances.push_back({"tallbox", identity_frame3f, 6, 6});
    scene.instances.push_back({"light", identity_frame3f, 7, 7});
}

}  // namespace yocto
