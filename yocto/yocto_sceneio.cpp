//
// Implementation for Yocto/Scene.
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

#include "yocto_sceneio.h"
#include "yocto_image.h"
#include "yocto_shape.h"
#include "yocto_utils.h"

#include <climits>
#include <array>
#include <deque>
#include <unordered_map>
using namespace std::string_literals;
#include "ext/json.hpp"

// -----------------------------------------------------------------------------
// GENERIC METHOD
// -----------------------------------------------------------------------------
namespace ygl {

// Load a scene
scene* load_scene(
    const std::string& filename, bool load_textures, bool skip_missing) {
    auto ext = path_extension(filename);
    auto scn = (scene*)nullptr;
    if (ext == ".json" || ext == ".JSON") {
        scn = load_json_scene(filename, load_textures, skip_missing);
    } else if (ext == ".obj" || ext == ".OBJ") {
        scn = load_obj_scene(filename, load_textures, skip_missing);
    } else if (ext == ".gltf" || ext == ".GLTF") {
        scn = load_gltf_scene(filename, load_textures, skip_missing);
    } else {
        throw std::runtime_error("unsupported extension " + ext);
    }
    auto mat = (material*)nullptr;
    for (auto ist : scn->instances) {
        if (ist->mat) continue;
        if (!mat) {
            mat = make_matte_material("<default>", {0.2f, 0.2f, 0.2f});
            scn->materials.push_back(mat);
        }
        ist->mat = mat;
    }
    return scn;
}

// Save a scene
void save_scene(const std::string& filename, const scene* scn,
    bool save_textures, bool skip_missing) {
    auto ext = path_extension(filename);
    if (ext == ".json" || ext == ".JSON") {
        save_json_scene(filename, scn, save_textures, skip_missing);
    } else if (ext == ".obj" || ext == ".OBJ") {
        save_obj_scene(filename, scn, save_textures, skip_missing);
    } else if (ext == ".gltf" || ext == ".GLTF") {
        save_gltf_scene(filename, scn, save_textures, skip_missing);
    } else {
        throw std::runtime_error("unsupported extension " + ext);
    }
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IO UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Fix path separators
std::string fix_path(const std::string& path_) {
    auto path = path_;
    for (auto& c : path)
        if (c == '\\') c = '/';
    return path;
};

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

// Load a text file
std::string load_text(const std::string& filename) {
    // https://stackoverflow.com/questions/2912520/read-file-contents-into-a-string-in-c
    auto f = fopen(filename.c_str(), "r");
    if (!f) throw std::runtime_error("could not load " + filename);
    // Determine file size
    fseek(f, 0, SEEK_END);
    auto size = ftell(f);
    rewind(f);
    // read
    auto buf = std::vector<char>(size);
    if (fread(buf.data(), sizeof(char), size, f) != size)
        throw std::runtime_error("could not load " + filename);
    return std::string(buf.data(), buf.size());
}

// Save a text file
void save_text(const std::string& filename, const std::string& str) {
    auto f = fopen(filename.c_str(), "wb");
    if (!f) throw std::runtime_error("cannot write file " + filename);
    auto num = fwrite(str.c_str(), 1, str.size(), f);
    if (num != str.size())
        throw std::runtime_error("cannot write file " + filename);
    fclose(f);
}

// Load a binary file
std::vector<byte> load_binary(const std::string& filename) {
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
}

// Save a binary file
void save_binary(const std::string& filename, const std::vector<byte>& data) {
    auto f = fopen(filename.c_str(), "wb");
    if (!f) throw std::runtime_error("cannot write file " + filename);
    auto num = fwrite(data.data(), 1, data.size(), f);
    if (num != data.size())
        throw std::runtime_error("cannot write file " + filename);
    fclose(f);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// JSON UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Json alias
using json = nlohmann::json;

// Load a JSON object
json load_json(const std::string& filename) {
    auto txt = load_text(filename);
    return json::parse(txt.begin(), txt.end());
}

// Save a JSON object
void save_json(const std::string& filename, const json& js) {
    save_text(filename, js.dump(2));
}

// json conversions
inline void to_json(json& js, const vec2i& val) {
    js = json((const std::array<int, 2>&)val);
}
inline void from_json(const json& js, vec2i& val) {
    (std::array<int, 2>&)val = js.get<std::array<int, 2>>();
}
inline void to_json(json& js, const vec3i& val) {
    js = json((const std::array<int, 3>&)val);
}
inline void from_json(const json& js, vec3i& val) {
    (std::array<int, 3>&)val = js.get<std::array<int, 3>>();
}
inline void to_json(json& js, const vec4i& val) {
    js = json((const std::array<int, 4>&)val);
}
inline void from_json(const json& js, vec4i& val) {
    (std::array<int, 4>&)val = js.get<std::array<int, 4>>();
}

inline void to_json(json& js, const vec2f& val) {
    js = json((const std::array<float, 2>&)val);
}
inline void from_json(const json& js, vec2f& val) {
    (std::array<float, 2>&)val = js.get<std::array<float, 2>>();
}
inline void to_json(json& js, const vec3f& val) {
    js = json((const std::array<float, 3>&)val);
}
inline void from_json(const json& js, vec3f& val) {
    (std::array<float, 3>&)val = js.get<std::array<float, 3>>();
}
inline void to_json(json& js, const vec4f& val) {
    js = json((const std::array<float, 4>&)val);
}
inline void from_json(const json& js, vec4f& val) {
    (std::array<float, 4>&)val = js.get<std::array<float, 4>>();
}
inline void to_json(json& js, const frame3f& val) {
    js = json((const std::array<float, 12>&)val);
}
inline void from_json(const json& js, frame3f& val) {
    (std::array<float, 12>&)val = js.get<std::array<float, 12>>();
}
inline void to_json(json& js, const mat4f& val) {
    js = json((const std::array<float, 16>&)val);
}
inline void from_json(const json& js, mat4f& val) {
    (std::array<float, 16>&)val = js.get<std::array<float, 16>>();
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// BUILTIN JSON FORMAT
// -----------------------------------------------------------------------------
namespace ygl {

// Serialize struct
void to_json(json& js, const camera& val) {
    js["name"] = val.name;
    js["frame"] = val.frame;
    js["ortho"] = val.ortho;
    js["width"] = val.width;
    js["height"] = val.height;
    js["focal"] = val.focal;
    js["focus"] = val.focus;
    js["aperture"] = val.aperture;
}

// Serialize struct
void from_json(const json& js, camera& val) {
    val.name = js.value("name", ""s);
    val.frame = js.value("frame", identity_frame3f);
    val.ortho = js.value("ortho", false);
    val.width = js.value("width", 0.036f);
    val.height = js.value("height", 0.024f);
    val.focal = js.value("focal", 0.050f);
    val.focus = js.value("focus", flt_max);
    val.aperture = js.value("aperture", 0.0f);
}

// Serialize struct
void to_json(json& js, const texture& val) {
    js["name"] = val.name;
    js["path"] = val.path;
    js["clamp"] = val.clamp;
    js["scale"] = val.scale;
    if (val.path == "") {
        js["width"] = val.width;
        js["height"] = val.height;
        js["img"] = val.img;
    }
}

// Serialize struct
void from_json(const json& js, texture& val) {
    val.name = js.value("name", ""s);
    val.path = js.value("path", ""s);
    val.clamp = js.value("clamp", false);
    val.scale = js.value("scale", 1.0f);
    val.width = js.value("width", 0);
    val.height = js.value("height", 0);
    val.img = js.value("img", std::vector<vec4f>());
}

// Serialize struct
void to_json(json& js, const material& val) {
    js["name"] = val.name;
    js["base_metallic"] = val.base_metallic;
    js["gltf_textures"] = val.gltf_textures;
    js["double_sided"] = val.double_sided;
    js["ke"] = val.ke;
    js["kd"] = val.kd;
    js["ks"] = val.ks;
    js["kt"] = val.kt;
    js["rs"] = val.rs;
    js["op"] = val.op;
    js["fresnel"] = val.fresnel;
    js["refract"] = val.refract;
    if (val.ke_txt) js["ke_txt"] = val.ke_txt->name;
    if (val.kd_txt) js["kd_txt"] = val.kd_txt->name;
    if (val.ks_txt) js["ks_txt"] = val.ks_txt->name;
    if (val.kt_txt) js["kt_txt"] = val.kt_txt->name;
    if (val.rs_txt) js["rs_txt"] = val.rs_txt->name;
    if (val.occ_txt) js["occ_txt"] = val.occ_txt->name;
    if (val.bump_txt) js["bump_txt"] = val.bump_txt->name;
    if (val.disp_txt) js["disp_txt"] = val.disp_txt->name;
    if (val.norm_txt) js["norm_txt"] = val.norm_txt->name;
}

// Serialize struct
void from_json(const json& js, material& val) {
    val.name = js.value("name", ""s);
    val.base_metallic = js.value("base_metallic", false);
    val.gltf_textures = js.value("gltf_textures", false);
    val.double_sided = js.value("double_sided", false);
    val.ke = js.value("ke", zero3f);
    val.kd = js.value("kd", zero3f);
    val.ks = js.value("ks", zero3f);
    val.kt = js.value("kt", zero3f);
    val.rs = js.value("rs", 1.0f);
    val.op = js.value("op", 1.0f);
    val.fresnel = js.value("fresnel", true);
    val.refract = js.value("refract", false);
    if (js.value("ke_txt", ""s) != "")
        val.ke_txt = new texture{js.value("ke_txt", ""s)};
    if (js.value("kd_txt", ""s) != "")
        val.kd_txt = new texture{js.value("kd_txt", ""s)};
    if (js.value("ks_txt", ""s) != "")
        val.ks_txt = new texture{js.value("ks_txt", ""s)};
    if (js.value("kt_txt", ""s) != "")
        val.kt_txt = new texture{js.value("kt_txt", ""s)};
    if (js.value("rs_txt", ""s) != "")
        val.rs_txt = new texture{js.value("rs_txt", ""s)};
    if (js.value("occ_txt", ""s) != "")
        val.occ_txt = new texture{js.value("occ_txt", ""s)};
    if (js.value("bump_txt", ""s) != "")
        val.bump_txt = new texture{js.value("bump_txt", ""s)};
    if (js.value("disp_txt", ""s) != "")
        val.disp_txt = new texture{js.value("disp_txt", ""s)};
    if (js.value("norm_txt", ""s) != "")
        val.norm_txt = new texture{js.value("norm_txt", ""s)};
}

// Serialize struct
void to_json(json& js, const shape& val) {
    js["name"] = val.name;
    js["path"] = val.path;
    if (val.path == "") {
        js["lines"] = val.lines;
        js["triangles"] = val.triangles;
        js["pos"] = val.pos;
        js["norm"] = val.norm;
        js["texcoord"] = val.texcoord;
        js["color"] = val.color;
        js["radius"] = val.radius;
        js["tangsp"] = val.tangsp;
    }
}

// Serialize struct
void from_json(const json& js, shape& val) {
    val.name = js.value("name", ""s);
    val.path = js.value("path", ""s);
    val.lines = js.value("lines", std::vector<vec2i>());
    val.triangles = js.value("triangles", std::vector<vec3i>());
    val.pos = js.value("pos", std::vector<vec3f>());
    val.norm = js.value("norm", std::vector<vec3f>());
    val.texcoord = js.value("texcoord", std::vector<vec2f>());
    val.color = js.value("color", std::vector<vec4f>());
    val.radius = js.value("radius", std::vector<float>());
    val.tangsp = js.value("tangsp", std::vector<vec4f>());
}

// Serialize struct
void to_json(json& js, const subdiv& val) {
    js["name"] = val.name;
    js["path"] = val.path;
    js["level"] = val.level;
    js["catmull_clark"] = val.catmull_clark;
    js["compute_normals"] = val.compute_normals;
    if (val.path == "") {
        js["quads_pos"] = val.quads_pos;
        js["quads_texcoord"] = val.quads_texcoord;
        js["quads_color"] = val.quads_color;
        js["pos"] = val.pos;
        js["texcoord"] = val.texcoord;
        js["color"] = val.color;
    }
}

// Serialize struct
void from_json(const json& js, subdiv& val) {
    val.name = js.value("name", ""s);
    val.path = js.value("path", ""s);
    val.level = js.value("level", 0);
    val.catmull_clark = js.value("catmull_clark", true);
    val.compute_normals = js.value("compute_normals", true);
    val.quads_pos = js.value("quads_pos", std::vector<vec4i>());
    val.quads_texcoord = js.value("quads_texcoord", std::vector<vec4i>());
    val.quads_color = js.value("quads_color", std::vector<vec4i>());
    val.pos = js.value("pos", std::vector<vec3f>());
    val.texcoord = js.value("texcoord", std::vector<vec2f>());
    val.color = js.value("color", std::vector<vec4f>());
}

// Serialize struct
void to_json(json& js, const instance& val) {
    js["name"] = val.name;
    js["frame"] = val.frame;
    if (val.shp) js["shp"] = val.shp->name;
    if (val.mat) js["mat"] = val.mat->name;
    if (val.sbd) js["sbd"] = val.sbd->name;
}

// Serialize struct
void from_json(const json& js, instance& val) {
    val.name = js.value("name", ""s);
    val.frame = js.value("frame", identity_frame3f);
    if (js.value("shp", ""s) != "") val.shp = new shape{js.value("shp", ""s)};
    if (js.value("mat", ""s) != "")
        val.mat = new material{js.value("mat", ""s)};
    if (js.value("sbd", ""s) != "") val.sbd = new subdiv{js.value("sbd", ""s)};
}

// Serialize struct
void to_json(json& js, const environment& val) {
    js["name"] = val.name;
    js["frame"] = val.frame;
    js["ke"] = val.ke;
    if (val.ke_txt) js["ke_txt"] = val.ke_txt->name;
}

// Serialize struct
void from_json(const json& js, environment& val) {
    val.name = js.value("name", ""s);
    val.frame = js.value("frame", identity_frame3f);
    val.ke = js.value("ke", zero3f);
    if (js.value("ke_txt", ""s) != "")
        val.ke_txt = new texture{js.value("ke_txt", ""s)};
}

// Serialize struct
void to_json(json& js, const node& val) {
    js["name"] = val.name;
    js["frame"] = val.frame;
    js["translation"] = val.translation;
    js["rotation"] = val.rotation;
    js["scale"] = val.scale;
    js["weights"] = val.weights;
    if (val.parent) js["parent"] = val.parent->name;
    if (val.cam) js["cam"] = val.cam->name;
    if (val.ist) js["ist"] = val.ist->name;
    if (val.env) js["env"] = val.env->name;
}

// Serialize struct
void from_json(const json& js, node& val) {
    val.name = js.value("name", ""s);
    val.frame = js.value("frame", identity_frame3f);
    val.translation = js.value("translation", zero3f);
    val.rotation = js.value("rotation", zero4f);
    val.scale = js.value("scale", zero3f);
    val.weights = js.value("weights", std::vector<float>());
    if (js.value("parent", ""s) != "")
        val.parent = new node{js.value("parent", ""s)};
    if (js.value("cam", ""s) != "") val.cam = new camera{js.value("cam", ""s)};
    if (js.value("ist", ""s) != "")
        val.ist = new instance{js.value("ist", ""s)};
    if (js.value("env", ""s) != "")
        val.env = new environment{js.value("env", ""s)};
}

// Serialize enum
void to_json(json& js, const animation_type& val) {
    static auto names = std::map<animation_type, std::string>{
        {animation_type::linear, "linear"},
        {animation_type::step, "step"},
        {animation_type::bezier, "bezier"},
    };
    js = names.at(val);
}

// Serialize enum
void from_json(const json& js, animation_type& val) {
    static auto names = std::map<std::string, animation_type>{
        {"linear", animation_type::linear},
        {"step", animation_type::step},
        {"bezier", animation_type::bezier},
    };
    val = names.at(js.get<std::string>());
}

// Serialize struct
void to_json(json& js, const animation& val) {
    js["name"] = val.name;
    js["path"] = val.path;
    js["group"] = val.group;
    js["type"] = val.type;
    if (val.path == "") {
        js["times"] = val.times;
        js["translation"] = val.translation;
        js["rotation"] = val.rotation;
        js["scale"] = val.scale;
    }
    if (!val.targets.empty()) {
        js["targets"] = json::array();
        for (auto v : val.targets) js["targets"].push_back(v->name);
    }
}

// Serialize struct
void from_json(const json& js, animation& val) {
    val.name = js.value("name", ""s);
    val.path = js.value("path", ""s);
    val.group = js.value("group", ""s);
    val.type = js.value("type", animation_type::linear);
    val.times = js.value("times", std::vector<float>());
    val.translation = js.value("translation", std::vector<vec3f>());
    val.rotation = js.value("rotation", std::vector<vec4f>());
    val.scale = js.value("scale", std::vector<vec3f>());
    for (auto& j : js.value("targets", json::array()))
        val.targets.push_back(new node{j.get<std::string>()});
}

// Serialize struct
void to_json(json& js, const scene& val) {
    js["name"] = val.name;
    if (!val.cameras.empty()) {
        js["cameras"] = json::array();
        for (auto v : val.cameras) js["cameras"].push_back(json(*v));
    }
    if (!val.textures.empty()) {
        js["textures"] = json::array();
        for (auto v : val.textures) js["textures"].push_back(json(*v));
    }
    if (!val.materials.empty()) {
        js["materials"] = json::array();
        for (auto v : val.materials) js["materials"].push_back(json(*v));
    }
    if (!val.shapes.empty()) {
        js["shapes"] = json::array();
        for (auto v : val.shapes) js["shapes"].push_back(json(*v));
    }
    if (!val.subdivs.empty()) {
        js["subdivs"] = json::array();
        for (auto v : val.subdivs) js["subdivs"].push_back(json(*v));
    }
    if (!val.instances.empty()) {
        js["instances"] = json::array();
        for (auto v : val.instances) js["instances"].push_back(json(*v));
    }
    if (!val.environments.empty()) {
        js["environments"] = json::array();
        for (auto v : val.environments) js["environments"].push_back(json(*v));
    }
    if (!val.nodes.empty()) {
        js["nodes"] = json::array();
        for (auto v : val.nodes) js["nodes"].push_back(json(*v));
    }
    if (!val.animations.empty()) {
        js["animations"] = json::array();
        for (auto v : val.animations) js["animations"].push_back(*v);
    }
}
void to_json(json& js, const scene*& val) {
    if (!val) {
        js = json();
        return;
    }
    to_json(js, *val);
}

template <typename T>
static std::unordered_map<std::string, T*> make_named_map(
    const std::vector<T*>& elems) {
    auto map = std::unordered_map<std::string, T*>();
    for (auto elem : elems) map[elem->name] = elem;
    return map;
};

// Serialize struct
void from_json(const json& js, scene& val) {
    val.name = js.value("name", ""s);
    for (auto& j : js.value("cameras", json::array()))
        val.cameras.push_back(new camera(j.get<camera>()));
    for (auto& j : js.value("textures", json::array()))
        val.textures.push_back(new texture(j.get<texture>()));
    for (auto& j : js.value("materials", json::array()))
        val.materials.push_back(new material(j.get<material>()));
    for (auto& j : js.value("shapes", json::array()))
        val.shapes.push_back(new shape(j.get<shape>()));
    for (auto& j : js.value("subdivs", json::array()))
        val.subdivs.push_back(new subdiv(j.get<subdiv>()));
    for (auto& j : js.value("instances", json::array()))
        val.instances.push_back(new instance(j.get<instance>()));
    for (auto& j : js.value("environments", json::array()))
        val.environments.push_back(new environment(j.get<environment>()));
    for (auto& j : js.value("nodes", json::array()))
        val.nodes.push_back(new node(j.get<node>()));
    for (auto& j : js.value("animations", json::array()))
        val.animations.push_back(new animation(j.get<animation>()));

    // fix references
    auto cmap = make_named_map(val.cameras);
    auto tmap = make_named_map(val.textures);
    auto mmap = make_named_map(val.materials);
    auto smap = make_named_map(val.shapes);
    auto rmap = make_named_map(val.subdivs);
    auto imap = make_named_map(val.instances);
    auto emap = make_named_map(val.environments);
    auto nmap = make_named_map(val.nodes);
    auto fix_ref = [](auto& map, auto*& ref) {
        if (!ref) return;
        auto name = ref->name;
        delete ref;
        ref = map.at(name);
    };
    for (auto mat : val.materials) {
        fix_ref(tmap, mat->ke_txt);
        fix_ref(tmap, mat->kd_txt);
        fix_ref(tmap, mat->ks_txt);
        fix_ref(tmap, mat->kt_txt);
        fix_ref(tmap, mat->op_txt);
        fix_ref(tmap, mat->rs_txt);
        fix_ref(tmap, mat->occ_txt);
        fix_ref(tmap, mat->norm_txt);
        fix_ref(tmap, mat->bump_txt);
        fix_ref(tmap, mat->disp_txt);
    }
    for (auto ist : val.instances) {
        fix_ref(mmap, ist->mat);
        fix_ref(smap, ist->shp);
        fix_ref(rmap, ist->sbd);
    }
    for (auto env : val.environments) { fix_ref(tmap, env->ke_txt); }
    for (auto nde : val.nodes) {
        fix_ref(nmap, nde->parent);
        fix_ref(cmap, nde->cam);
        fix_ref(imap, nde->ist);
        fix_ref(emap, nde->env);
    }
    for (auto anm : val.animations) {
        for (auto& nde : anm->targets) fix_ref(nmap, nde);
    }
}

// Load a scene in the builtin JSON format.
scene* load_json_scene(
    const std::string& filename, bool load_textures, bool skip_missing) {
    // load json
    auto scn = new scene();
    auto js = load_json(filename);
    from_json(js, *scn);

    // load meshes
    auto dirname = path_dirname(filename);
    for (auto shp : scn->shapes) {
        if (shp->path == "") continue;
        auto filename = dirname + shp->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        try {
            load_mesh(filename, shp->lines, shp->triangles, shp->pos, shp->norm,
                shp->texcoord, shp->color, shp->radius);
        } catch (std::exception&) {
            if (skip_missing) continue;
            throw;
        }
    }

    // skip textures
    if (!load_textures) return scn;

    // set gamma
    auto ldr_gamma = std::unordered_map<texture*, float>{{nullptr, 1.0f}};
    for (auto txt : scn->textures) ldr_gamma[txt] = 2.2f;
    for (auto mat : scn->materials) {
        ldr_gamma[mat->ke_txt] = 2.2f;
        ldr_gamma[mat->kd_txt] = 2.2f;
        ldr_gamma[mat->ks_txt] = 2.2f;
        ldr_gamma[mat->kt_txt] = 2.2f;
        ldr_gamma[mat->rs_txt] = 1;
        ldr_gamma[mat->op_txt] = 1;
        ldr_gamma[mat->norm_txt] = 1;
        ldr_gamma[mat->disp_txt] = 1;
        ldr_gamma[mat->bump_txt] = 1;
    }
    for (auto env : scn->environments) { ldr_gamma[env->ke_txt] = 2.2f; }

    // load images
    for (auto txt : scn->textures) {
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        txt->img =
            load_image(filename, txt->width, txt->height, ldr_gamma.at(txt));
        if (txt->img.empty()) {
            if (skip_missing) continue;
            throw std::runtime_error("cannot laod image " + filename);
        }
    }

    return scn;
}

// Save a scene in the builtin JSON format.
void save_json_scene(const std::string& filename, const scene* scn,
    bool save_textures, bool skip_missing) {
    // save json
    auto js = json(scn);
    save_json(filename, js);

    // save meshes
    auto dirname = path_dirname(filename);
    for (auto shp : scn->shapes) {
        if (shp->path == "") continue;
        auto filename = dirname + shp->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        try {
            save_mesh(filename, shp->lines, shp->triangles, shp->pos, shp->norm,
                shp->texcoord, shp->color, shp->radius);
        } catch (std::exception&) {
            if (skip_missing) continue;
            throw;
        }
    }

    // skip textures
    if (!save_textures) return;

    // set gamma
    auto ldr_gamma = std::unordered_map<texture*, float>{{nullptr, 1.0f}};
    for (auto txt : scn->textures) ldr_gamma[txt] = 2.2f;
    for (auto mat : scn->materials) {
        ldr_gamma[mat->ke_txt] = 2.2f;
        ldr_gamma[mat->kd_txt] = 2.2f;
        ldr_gamma[mat->ks_txt] = 2.2f;
        ldr_gamma[mat->kt_txt] = 2.2f;
        ldr_gamma[mat->rs_txt] = 1;
        ldr_gamma[mat->op_txt] = 1;
        ldr_gamma[mat->norm_txt] = 1;
        ldr_gamma[mat->disp_txt] = 1;
        ldr_gamma[mat->bump_txt] = 1;
    }
    for (auto env : scn->environments) { ldr_gamma[env->ke_txt] = 2.2f; }

    // save images
    for (auto txt : scn->textures) {
        if (txt->img.empty()) continue;
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        try {
            save_image(
                filename, txt->width, txt->height, txt->img, ldr_gamma.at(txt));
        } catch (std::exception&) {
            if (skip_missing) continue;
            throw;
        }
    }
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// OBJ CONVESION
// -----------------------------------------------------------------------------
namespace ygl {

// skip whitespace
inline void parse_skipws(char*& s) {
    while (*s == ' ' || *s == '\t' || *s == '\r' || *s == '\n') s++;
}
// check if whitespace
inline bool parse_isws(char*& s) {
    return *s == ' ' || *s == '\t' || *s == '\r' || *s == '\n';
}
// check if whitespace
inline void parse_remove_comment(char*& s, char cmd_char) {
    auto ss = s;
    while (*ss) {
        if (*ss == cmd_char) {
            *ss = 0;
            break;
        }
        ss++;
    }
}

#if YGL_FASTPARSE
// parse base value
inline int parse_int_base(char*& s) {
    parse_skipws(s);
    auto val = 0;
    auto sn = (*s == '-') ? -1 : 1;
    if (*s == '-' || *s == '+') s++;
    while (*s >= '0' && *s <= '9') val = val * 10 + (*s++ - '0');
    val *= sn;
    return val;
}

// parse base value
inline float parse_float_base(char*& s) {
    parse_skipws(s);
    //    auto ss = s; auto sss = ss;
    auto mantissa = 0, fractional = 0, fractional_length = 0, exponent = 0;
    auto sn = (*s == '-') ? -1 : 1;
    if (*s == '-' || *s == '+') s++;
    while (*s >= '0' && *s <= '9') mantissa = mantissa * 10 + (*s++ - '0');
    if (*s == '.') {
        s++;
        while (*s >= '0' && *s <= '9') {
            fractional = fractional * 10 + (*s++ - '0');
            fractional_length++;
        }
    }
    mantissa *= sn;
    fractional *= sn;
    if (*s == 'e' || *s == 'E') {
        s++;
        auto en = (*s == '-') ? -1 : 1;
        if (*s == '-' || *s == '+') s++;
        while (*s >= '0' && *s <= '9') exponent = exponent * 10 + (*s++ - '0');
        exponent *= en;
    }
    auto dval = (double)mantissa;
    if (fractional)
        dval += fractional * std::pow(10.0, -(double)fractional_length);
    if (exponent) dval *= std::pow(10.0, (double)exponent);
    auto val = (float)dval;
    return val;
#if 0
    auto cval = val;
    sscanf(ss, "%f", &cval);
    if(abs(val - cval) > 0.01f) {
        printf("- %g %g %s\n", val, cval, sss);
    }
    auto len = 0;
    sscanf(s, "%f%n", &val, &len);
    s += len;
#endif
}

// parse base value
inline std::string parse_string_base(char*& s) {
    parse_skipws(s);
    char buf[4096];
    auto val = buf;
    while (*s && !parse_isws(s)) *val++ = *s++;
    *val = 0;
    return buf;
}

#else

// parse base value
inline int parse_int_base(char*& s) {
    auto len = 0;
    auto val = 0;
    sscanf(s, "%d%n", &val, &len);
    s += len;
    return val;
}

// parse base value
inline float parse_float_base(char*& s) {
    auto len = 0;
    auto val = 0.0f;
    sscanf(s, "%f%n", &val, &len);
    s += len;
    return val;
}

// parse base value
inline std::string parse_string_base(char*& s) {
    char buf[4096];
    auto val = buf;
    auto len = 0;
    sscanf(s, "%s%n", buf, &len);
    if (len) {
        s += len;
        val = buf;
    } else {
        val = "";
    }
    return val;
}

#endif

// parse value
inline int parse_int(char*& s) { return parse_int_base(s); }
inline float parse_float(char*& s) { return parse_float_base(s); }
inline std::string parse_string(char*& s) { return parse_string_base(s); }
inline byte parse_uchar(char*& s) { return (byte)parse_int_base(s); }
inline std::vector<std::string> parse_strings(char*& s) {
    auto vals = std::vector<std::string>();
    parse_skipws(s);
    while (*s) {
        auto val = parse_string(s);
        if (val == "") break;
        vals.push_back(val);
        parse_skipws(s);
    }
    return vals;
}
inline float parse_bool(char*& s) { return parse_int(s); }
inline vec2i parse_vec2i(char*& s) {
    auto x = parse_int(s);
    auto y = parse_int(s);
    return {x, y};
}
inline vec2f parse_vec2f(char*& s) {
    auto x = parse_float(s);
    auto y = parse_float(s);
    return {x, y};
}
inline vec3f parse_vec3f(char*& s) {
    auto x = parse_float(s);
    auto y = parse_float(s);
    auto z = parse_float(s);
    return {x, y, z};
}
inline vec4f parse_vec4f(char*& s) {
    auto x = parse_float(s);
    auto y = parse_float(s);
    auto z = parse_float(s);
    auto w = parse_float(s);
    return {x, y, z, w};
}
inline frame3f parse_frame3f(char*& s) {
    auto val = identity_frame3f;
    for (auto i = 0; i < 12; i++) (&val.x.x)[i] = parse_float(s);
    return val;
}
inline bool parse_getline(FILE* f, char* buf, int max_length) {
    if (fgets(buf, max_length, f)) {
        auto ss = buf;
        while (*ss) {
            if (*ss == '\r' || *ss == '\n') {
                *ss = 0;
                break;
            }
            ss++;
        }
        return true;
    } else {
        return false;
    }
}

inline vec3i parse_objvert(char*& s, const vec3i& vert_size) {
    auto token = parse_string(s);
    auto val = vec3i{-1, -1, -1};
    auto i = 0;
    auto sb = (char*)token.c_str();
    while (i < 3 && *sb) {
        (&val.x)[i] = parse_int(sb);
        (&val.x)[i] = ((&val.x)[i] < 0) ? (&vert_size.x)[i] + (&val.x)[i] :
                                          (&val.x)[i] - 1;
        if (*sb != '/') break;
        while (*sb == '/') {
            sb++;
            i++;
        }
    }
    return val;
}

template <typename T>
inline T read_value(FILE* fs) {
    T val;
    if (fread(&val, sizeof(val), 1, fs) != 1)
        throw std::runtime_error("error reading value");
    return val;
}

inline int read_int(FILE* fs) { return read_value<int>(fs); }
inline float read_float(FILE* fs) { return read_value<float>(fs); }
inline byte read_uchar(FILE* fs) { return read_value<byte>(fs); }

// Loads an OBJ
scene* load_obj_scene(const std::string& filename, bool load_textures,
    bool skip_missing, bool split_shapes) {
    auto scn = new scene();

    // parsing policies
    auto flip_texcoord = true;
    auto flip_tr = true;

    // splitting policy
    auto split_material = split_shapes;
    auto split_group = split_shapes;
    auto split_smoothing = split_shapes;

    // open file
    auto fs = fopen(filename.c_str(), "rt");
    if (!fs) throw std::runtime_error("cannot open filename " + filename);

    // current parsing values
    auto matname = std::string();
    auto oname = std::string();
    auto gname = std::string();
    auto smoothing = true;
    auto ist = (instance*)nullptr;

    // vertices
    auto pos = std::deque<vec3f>();
    auto norm = std::deque<vec3f>();
    auto texcoord = std::deque<vec2f>();

    // object maps
    auto tmap = std::unordered_map<std::string, texture*>();
    auto mmap = std::unordered_map<std::string, material*>();

    // vertex maps
    auto name_map = std::unordered_map<std::string, int>();
    auto vert_map = std::unordered_map<vec3i, int>();
    auto pos_map = std::unordered_map<int, int>();
    auto norm_map = std::unordered_map<int, int>();
    auto texcoord_map = std::unordered_map<int, int>();

    // add object if needed
    auto is_instance_empty = [](instance* ist) {
        if (ist->sbd) {
            return ist->sbd->pos.empty();
        } else if (ist->shp) {
            return ist->shp->pos.empty();
        } else {
            return true;
        }
    };
    auto add_instance = [&](scene* scn, const std::string& objname,
                            const std::string& matname,
                            const std::string& groupname, bool smoothing) {
        if (scn->instances.empty() || objname != scn->instances.back()->name ||
            !is_instance_empty(scn->instances.back())) {
            auto ist = new instance();
            scn->instances.push_back(ist);
            ist->shp = new shape();
            scn->shapes.push_back(ist->shp);
        }
        name_map[objname] += 1;
        auto name = (name_map[objname] == 1) ?
                        objname :
                        objname + "_" + std::to_string(name_map[objname] - 1);
        if (objname == "") name = "object" + name;
        auto ist = scn->instances.back();
        ist->name = name;
        if (ist->shp) ist->shp->name = ist->name;
        if (ist->sbd) ist->sbd->name = ist->name;
        if (matname != "") ist->mat = mmap.at(matname);
        vert_map.clear();
        pos_map.clear();
        norm_map.clear();
        texcoord_map.clear();
        return ist;
    };

    // current objet
    ist = add_instance(scn, "", "", "", true);

    // read the file line by line
    char line[4096];
    while (parse_getline(fs, line, sizeof(line))) {
        // prepare to parse
        auto ss = line;
        parse_skipws(ss);
        parse_remove_comment(ss, '#');

        // get command
        auto cmd = parse_string(ss);
        if (cmd == "") continue;

        // possible token values
        if (cmd == "v") {
            pos.push_back(parse_vec3f(ss));
        } else if (cmd == "vn") {
            norm.push_back(parse_vec3f(ss));
        } else if (cmd == "vt") {
            texcoord.push_back(parse_vec2f(ss));
            if (flip_texcoord) texcoord.back().y = 1 - texcoord.back().y;
        } else if (cmd == "f" || cmd == "l") {
            auto num = 0;
            vec3i verts[128];
            int vids[128];
            auto vert_size =
                vec3i{(int)pos.size(), (int)texcoord.size(), (int)norm.size()};
            // elem.material = (int)oobj->materials.size() - 1;
            parse_skipws(ss);
            while (*ss) {
                verts[num++] = parse_objvert(ss, vert_size);
                parse_skipws(ss);
            }
            if (ist->sbd) {
                // TODO: subdivs
            } else if (ist->shp) {
                for (auto i = 0; i < num; i++) {
                    auto it = vert_map.find(verts[i]);
                    if (it == vert_map.end()) {
                        auto nverts = (int)ist->shp->pos.size();
                        vert_map.insert(it, {verts[i], nverts});
                        vids[i] = nverts;
                        if (verts[i].x >= 0)
                            ist->shp->pos.push_back(pos.at(verts[i].x));
                        if (verts[i].y >= 0)
                            ist->shp->texcoord.push_back(
                                texcoord.at(verts[i].y));
                        if (verts[i].z >= 0)
                            ist->shp->norm.push_back(norm.at(verts[i].z));
                    } else {
                        vids[i] = it->second;
                    }
                }
                if (cmd == "f") {
                    for (auto i = 2; i < num; i++)
                        ist->shp->triangles.push_back(
                            {vids[0], vids[i - 1], vids[i]});
                }
                if (cmd == "l") {
                    for (auto i = 1; i < num; i++)
                        ist->shp->lines.push_back({vids[i - 1], vids[i]});
                }
            }
        } else if (cmd == "p") {
            // TODO: support points
        } else if (cmd == "o") {
            oname = parse_string(ss);
            gname = "";
            matname = "";
            smoothing = true;
            ist = add_instance(scn, oname, matname, gname, smoothing);
        } else if (cmd == "usemtl") {
            matname = parse_string(ss);
            if (split_material) {
                ist = add_instance(scn, oname, matname, gname, smoothing);
            } else {
                if (matname != "") ist->mat = mmap.at(matname);
            }
        } else if (cmd == "g") {
            gname = parse_string(ss);
            if (split_group) {
                ist = add_instance(scn, oname, matname, gname, smoothing);
            }
        } else if (cmd == "s") {
            auto name = parse_string(ss);
            smoothing = (name == "on");
            if (split_smoothing) {
                ist = add_instance(scn, oname, matname, gname, smoothing);
            }
        } else if (cmd == "mtllib") {
            auto mtlname = parse_string(ss);
            auto mtlpath = path_dirname(filename) + mtlname;
            // open file
            auto fs = fopen(mtlpath.c_str(), "rt");
            if (!fs)
                throw std::runtime_error("cannot open filename " + mtlpath);

            // Parse texture options and name
            auto add_texture = [scn, &tmap](char*& s) {
                // get tokens
                auto tokens = parse_strings(s);
                if (tokens.empty()) return (texture*)nullptr;

                // texture name
                auto path = tokens.back();
                for (auto& c : path)
                    if (c == '\\') c = '/';
                if (tmap.find(path) != tmap.end()) { return tmap.at(path); }

                // create texture
                auto txt = new texture();
                txt->name = path;
                txt->path = path;
                scn->textures.push_back(txt);
                tmap[path] = txt;

                // texture options
                for (auto i = 0; i < tokens.size(); i++) {
                    if (tokens[i] == "-clamp") txt->clamp = true;
                    // TODO: bump scale
                }

                return txt;
            };

            // add a material preemptively to avoid crashes
            scn->materials.push_back(new material());
            auto mat = scn->materials.back();

            // read the file line by line
            char line[4096];
            while (parse_getline(fs, line, sizeof(line))) {
                // prepare to parse
                auto ss = line;
                parse_skipws(ss);
                parse_remove_comment(ss, '#');

                // get command
                auto cmd = parse_string(ss);
                if (cmd == "") continue;

                // possible token values
                if (cmd == "newmtl") {
                    mat = new material();
                    mat->name = parse_string(ss);
                    scn->materials.push_back(mat);
                    mmap[mat->name] = mat;
                } else if (cmd == "illum") {
                    // auto illum = parse_int(ss);
                    // TODO: something with illum
                } else if (cmd == "Ke") {
                    mat->ke = parse_vec3f(ss);
                } else if (cmd == "Kd") {
                    mat->kd = parse_vec3f(ss);
                } else if (cmd == "Ks") {
                    mat->ks = parse_vec3f(ss);
                } else if (cmd == "Kt") {
                    mat->kt = parse_vec3f(ss);
                } else if (cmd == "Tf") {
                    auto nchan = 0;
                    while (*ss && nchan < 3) {
                        (&mat->kt.x)[nchan++] = parse_float(ss);
                        parse_skipws(ss);
                    }
                    if (nchan < 3) mat->kt = {mat->kt.x, mat->kt.x, mat->kt.x};
                    if (flip_tr) mat->kt = vec3f{1, 1, 1} - mat->kt;
                } else if (cmd == "Tr") {
                    auto nchan = 0;
                    auto tr = zero3f;
                    while (*ss && nchan < 3) {
                        (&tr.x)[nchan++] = parse_float(ss);
                        parse_skipws(ss);
                    }
                    if (nchan < 3) tr = {tr.x, tr.x, tr.x};
                    mat->op = (tr.x + tr.y + tr.z) / 3;
                    if (flip_tr) mat->op = 1 - mat->op;
                } else if (cmd == "Ns") {
                    auto ns = parse_float(ss);
                    mat->rs = pow(2 / (ns + 2), 1 / 4.0f);
                    if (mat->rs < 0.01f) mat->rs = 0;
                    if (mat->rs > 0.99f) mat->rs = 1;
                } else if (cmd == "d") {
                    mat->op = parse_float(ss);
                } else if (cmd == "Pr" || cmd == "rs") {
                    mat->rs = parse_float(ss);
                } else if (cmd == "map_Ke") {
                    mat->ke_txt = add_texture(ss);
                } else if (cmd == "map_Kd") {
                    mat->kd_txt = add_texture(ss);
                } else if (cmd == "map_Ks") {
                    mat->ks_txt = add_texture(ss);
                } else if (cmd == "map_Tr") {
                    mat->kt_txt = add_texture(ss);
                } else if (cmd == "map_d" || cmd == "map_Tr") {
                    mat->op_txt = add_texture(ss);
                } else if (cmd == "map_Pr" || cmd == "map_rs") {
                    mat->rs_txt = add_texture(ss);
                } else if (cmd == "map_occ" || cmd == "occ") {
                    mat->occ_txt = add_texture(ss);
                } else if (cmd == "map_bump" || cmd == "bump") {
                    mat->bump_txt = add_texture(ss);
                } else if (cmd == "map_disp" || cmd == "disp") {
                    mat->disp_txt = add_texture(ss);
                } else if (cmd == "map_norm" || cmd == "norm") {
                    mat->norm_txt = add_texture(ss);
                }
            }

            // remove first fake material
            if (scn->materials.front()->name == "")
                scn->materials.erase(scn->materials.begin());

            // clone
            fclose(fs);
        } else if (cmd == "c") {
            auto cam = new camera();
            cam->name = parse_string(ss);
            cam->ortho = parse_bool(ss);
            cam->width = parse_float(ss);
            cam->height = parse_float(ss);
            cam->focal = parse_float(ss);
            cam->focus = parse_float(ss);
            cam->aperture = parse_float(ss);
            cam->frame = parse_frame3f(ss);
            scn->cameras.push_back(cam);
        } else if (cmd == "e") {
            auto env = new environment();
            env->name = parse_string(ss);
            env->ke = parse_vec3f(ss);
            auto ke_txt = parse_string(ss);
            if (ke_txt != "\"\"") {
                if (tmap.find(ke_txt) == tmap.end()) {
                    auto txt = new texture();
                    txt->name = ke_txt;
                    txt->path = ke_txt;
                    tmap[ke_txt] = txt;
                    scn->textures.push_back(txt);
                }
                env->ke_txt = tmap.at(ke_txt);
            }
            env->frame = parse_frame3f(ss);
            scn->environments.push_back(env);
        } else {
            // unused
        }
    }

    // cleanup empty
    for (auto idx = 0; idx < scn->instances.size(); idx++) {
        if (!is_instance_empty(scn->instances[idx])) continue;
        auto ist = scn->instances[idx];
        if (ist->shp) {
            scn->shapes.erase(
                std::find(scn->shapes.begin(), scn->shapes.end(), ist->shp));
            delete ist->shp;
        }
        if (ist->sbd) {
            scn->subdivs.erase(
                std::find(scn->subdivs.begin(), scn->subdivs.end(), ist->sbd));
            delete ist->sbd;
        }
        delete ist;
        scn->instances.erase(scn->instances.begin() + idx);
        idx--;
    }

    // close file
    fclose(fs);

    // updates
    update_bbox(scn);

    // fix scene
    scn->name = path_filename(filename);
    auto mat = (material*)nullptr;
    for (auto ist : scn->instances) {
        if (ist->mat) continue;
        if (!mat) {
            mat = make_matte_material("<default>", {0.2f, 0.2f, 0.2f});
            scn->materials.push_back(mat);
        }
        ist->mat = mat;
    }

    // skip if needed
    if (!load_textures) return scn;

    // set gamma
    auto ldr_gamma = std::unordered_map<texture*, float>{{nullptr, 1.0f}};
    for (auto txt : scn->textures) ldr_gamma[txt] = 2.2f;
    for (auto mat : scn->materials) {
        ldr_gamma[mat->ke_txt] = 2.2f;
        ldr_gamma[mat->kd_txt] = 2.2f;
        ldr_gamma[mat->ks_txt] = 2.2f;
        ldr_gamma[mat->kt_txt] = 2.2f;
        ldr_gamma[mat->rs_txt] = 1;
        ldr_gamma[mat->op_txt] = 1;
        ldr_gamma[mat->norm_txt] = 1;
        ldr_gamma[mat->disp_txt] = 1;
        ldr_gamma[mat->bump_txt] = 1;
    }
    for (auto env : scn->environments) { ldr_gamma[env->ke_txt] = 2.2f; }

    // load images
    auto dirname = path_dirname(filename);
    for (auto txt : scn->textures) {
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        txt->img =
            load_image(filename, txt->width, txt->height, ldr_gamma.at(txt));
        if (txt->img.empty()) {
            if (skip_missing) continue;
            throw std::runtime_error("cannot laod image " + filename);
        }
    }

    // done
    return scn;
}

void save_obj_scene(const std::string& filename, const scene* scn,
    bool save_textures, bool skip_missing) {
    // scene
    auto f = fopen(filename.c_str(), "wt");
    if (!f) throw std::runtime_error("cannot save file " + filename);

    // material library
    if (!scn->materials.empty()) {
        auto mtlname = path_basename(filename) + ".mtl";
        fprintf(f, "mtllib %s\n", mtlname.c_str());
    }

    // cameras
    for (auto cam : scn->cameras) {
        fprintf(f, "c");
        fprintf(f, " %s", cam->name.c_str());
        fprintf(f, " %d", (int)cam->ortho);
        fprintf(f, " %g", cam->width);
        fprintf(f, " %g", cam->height);
        fprintf(f, " %g", cam->focal);
        fprintf(f, " %g", cam->focus);
        fprintf(f, " %g", cam->aperture);
        for (auto i = 0; i < 12; i++) fprintf(f, " %g", (&cam->frame.x.x)[i]);
        fprintf(f, "\n");
    }

    // environments
    for (auto env : scn->environments) {
        fprintf(f, "e");
        fprintf(f, " %s", env->name.c_str());
        fprintf(f, " %g %g %g", env->ke.x, env->ke.y, env->ke.z);
        fprintf(f, " %s", (env->ke_txt) ? env->ke_txt->path.c_str() : "\"\"");
        for (auto i = 0; i < 12; i++) fprintf(f, " %g", (&env->frame.x.x)[i]);
        fprintf(f, "\n");
    }

    // shapes
    auto offset = zero3i;
    for (auto ist : scn->instances) {
        fprintf(f, "o %s\n", ist->name.c_str());
        if (ist->mat) fprintf(f, "usemtl %s\n", ist->mat->name.c_str());
        if (!ist->sbd) {
            if (ist->frame == identity_frame3f) {
                for (auto& p : ist->shp->pos)
                    fprintf(f, "v %g %g %g\n", p.x, p.y, p.z);
                for (auto& n : ist->shp->norm)
                    fprintf(f, "vn %g %g %g\n", n.x, n.y, n.z);
                for (auto& t : ist->shp->texcoord)
                    fprintf(f, "vt %g %g\n", t.x, t.y);
            } else {
                for (auto& p : ist->shp->pos) {
                    auto tp = transform_point(ist->frame, p);
                    fprintf(f, "v %g %g %g\n", tp.x, tp.y, tp.z);
                }
                for (auto& n : ist->shp->norm) {
                    auto tn = transform_direction(ist->frame, n);
                    fprintf(f, "vn %g %g %g\n", tn.x, tn.y, tn.z);
                }
                for (auto& t : ist->shp->texcoord)
                    fprintf(f, "vt %g %g\n", t.x, t.y);
            }
            if (!ist->shp->texcoord.empty() && !ist->shp->norm.empty()) {
                for (auto& t : ist->shp->triangles)
                    fprintf(f, "f %d/%d/%d %d/%d/%d %d/%d/%d\n",
                        t.x + offset.x + 1, t.x + offset.y + 1,
                        t.x + offset.z + 1, t.y + offset.x + 1,
                        t.y + offset.y + 1, t.y + offset.z + 1,
                        t.z + offset.x + 1, t.z + offset.y + 1,
                        t.z + offset.z + 1);
                for (auto& l : ist->shp->lines)
                    fprintf(f, "l %d/%d/%d %d/%d/%d\n", l.x + offset.x + 1,
                        l.x + offset.y + 1, l.x + offset.z + 1,
                        l.y + offset.x + 1, l.y + offset.y + 1,
                        l.y + offset.z + 1);
            } else if (!ist->shp->texcoord.empty() && ist->shp->norm.empty()) {
                for (auto& t : ist->shp->triangles)
                    fprintf(f, "f %d/%d %d/%d %d/%d\n", t.x + offset.x + 1,
                        t.x + offset.y + 1, t.y + offset.x + 1,
                        t.y + offset.y + 1, t.z + offset.x + 1,
                        t.z + offset.y + 1);
                for (auto& l : ist->shp->lines)
                    fprintf(f, "l %d/%d %d/%d\n", l.x + offset.x + 1,
                        l.x + offset.y + 1, l.y + offset.x + 1,
                        l.y + offset.y + 1);
            } else if (ist->shp->texcoord.empty() && !ist->shp->norm.empty()) {
                for (auto& t : ist->shp->triangles)
                    fprintf(f, "f %d//%d %d//%d %d//%d\n", t.x + offset.x + 1,
                        t.x + offset.z + 1, t.y + offset.x + 1,
                        t.y + offset.z + 1, t.z + offset.x + 1,
                        t.z + offset.z + 1);
                for (auto& l : ist->shp->lines)
                    fprintf(f, "l %d//%d %d//%d\n", l.x + offset.x + 1,
                        l.x + offset.z + 1, l.y + offset.x + 1,
                        l.y + offset.z + 1);
            } else {
                for (auto& t : ist->shp->triangles)
                    fprintf(f, "f %d %d %d\n", t.x + offset.x + 1,
                        t.y + offset.x + 1, t.z + offset.x + 1);
                for (auto& l : ist->shp->lines)
                    fprintf(
                        f, "l %d %d\n", l.x + offset.x + 1, l.y + offset.x + 1);
            }
            offset.x += ist->shp->pos.size();
            offset.y += ist->shp->texcoord.size();
            offset.z += ist->shp->norm.size();
        } else {
            if (ist->frame == identity_frame3f) {
                for (auto& p : ist->sbd->pos)
                    fprintf(f, "v %g %g %g\n", p.x, p.y, p.z);
                for (auto& t : ist->sbd->texcoord)
                    fprintf(f, "vt %g %g\n", t.x, t.y);
            } else {
                for (auto& p : ist->sbd->pos) {
                    auto tp = transform_point(ist->frame, p);
                    fprintf(f, "v %g %g %g\n", tp.x, tp.y, tp.z);
                }
                for (auto& t : ist->sbd->texcoord)
                    fprintf(f, "vt %g %g\n", t.x, t.y);
            }
            if (!ist->sbd->texcoord.empty()) {
                for (auto i = 0; i < ist->sbd->quads_pos.size(); i++) {
                    auto qp = ist->sbd->quads_pos[i];
                    auto qt = ist->sbd->quads_texcoord[i];
                    if (qp.z == qp.w) {
                        fprintf(f, "f %d/%d %d/%d %d/%d\n", qp.x + offset.x + 1,
                            qt.x + offset.y + 1, qp.y + offset.x + 1,
                            qt.y + offset.y + 1, qp.z + offset.x + 1,
                            qt.z + offset.y + 1);
                    } else {
                        fprintf(f, "f %d/%d %d/%d %d/%d %d/%d\n",
                            qp.x + offset.x + 1, qt.x + offset.y + 1,
                            qp.y + offset.x + 1, qt.y + offset.y + 1,
                            qp.z + offset.x + 1, qt.z + offset.y + 1,
                            qp.w + offset.x + 1, qt.w + offset.y + 1);
                    }
                }
            } else {
                for (auto& q : ist->sbd->quads_pos) {
                    if (q.z == q.w) {
                        fprintf(f, "f %d %d %d\n", q.x + offset.x + 1,
                            q.y + offset.x + 1, q.z + offset.x + 1);
                    } else {
                        fprintf(f, "f %d %d %d %d\n", q.x + offset.x + 1,
                            q.y + offset.x + 1, q.z + offset.x + 1,
                            q.w + offset.x + 1);
                    }
                }
            }
            offset.x += ist->sbd->pos.size();
            offset.y += ist->sbd->texcoord.size();
        }
    }

    fclose(f);

    // save materials
    if (scn->materials.empty()) return;

    auto mtlname = replace_path_extension(filename, ".mtl");
    f = fopen(mtlname.c_str(), "wt");
    if (!f) throw std::runtime_error("cannot open filename " + mtlname);

    // for each material, dump all the values
    for (auto mat : scn->materials) {
        fprintf(f, "newmtl %s\n", mat->name.c_str());
        fprintf(f, "  illum %d\n", 2);
        if (mat->ke != zero3f)
            fprintf(f, "  Ke %g %g %g\n", mat->ke.x, mat->ke.y, mat->ke.z);
        if (mat->kd != zero3f)
            fprintf(f, "  Kd %g %g %g\n", mat->kd.x, mat->kd.y, mat->kd.z);
        if (mat->ks != zero3f)
            fprintf(f, "  Ks %g %g %g\n", mat->ks.x, mat->ks.y, mat->ks.z);
        if (mat->kt != zero3f)
            fprintf(f, "  Kt %g %g %g\n", mat->kt.x, mat->kt.y, mat->kt.z);
        if (mat->rs != 1.0f)
            fprintf(f, "  Ns %d\n",
                (int)clamp(2 / pow(mat->rs + 1e-10f, 4.0f) - 2, 0.0f, 1.0e12f));
        if (mat->op != 1.0f) fprintf(f, "  d %g\n", mat->op);
        if (mat->rs != -1.0f) fprintf(f, "  Pr %g\n", mat->rs);
        if (mat->ke_txt) fprintf(f, "  map_Ke %s\n", mat->ke_txt->path.c_str());
        if (mat->kd_txt) fprintf(f, "  map_Kd %s\n", mat->kd_txt->path.c_str());
        if (mat->ks_txt) fprintf(f, "  map_Ks %s\n", mat->ks_txt->path.c_str());
        if (mat->kt_txt) fprintf(f, "  map_Kt %s\n", mat->kt_txt->path.c_str());
        if (mat->op_txt) fprintf(f, "  map_d %s\n", mat->op_txt->path.c_str());
        if (mat->rs_txt) fprintf(f, "  map_Pr %s\n", mat->rs_txt->path.c_str());
        if (mat->occ_txt)
            fprintf(f, "  map_occ %s\n", mat->occ_txt->path.c_str());
        if (mat->bump_txt)
            fprintf(f, "  map_bump %s\n", mat->bump_txt->path.c_str());
        if (mat->disp_txt)
            fprintf(f, "  map_disp %s\n", mat->disp_txt->path.c_str());
        if (mat->norm_txt)
            fprintf(f, "  map_norm %s\n", mat->norm_txt->path.c_str());
        fprintf(f, "\n");
    }

    fclose(f);

    // skip textures if needed
    if (!save_textures) return;

    // set gamma
    auto ldr_gamma = std::unordered_map<texture*, float>{{nullptr, 1.0f}};
    for (auto txt : scn->textures) ldr_gamma[txt] = 2.2f;
    for (auto mat : scn->materials) {
        ldr_gamma[mat->ke_txt] = 2.2f;
        ldr_gamma[mat->kd_txt] = 2.2f;
        ldr_gamma[mat->ks_txt] = 2.2f;
        ldr_gamma[mat->kt_txt] = 2.2f;
        ldr_gamma[mat->rs_txt] = 1;
        ldr_gamma[mat->op_txt] = 1;
        ldr_gamma[mat->norm_txt] = 1;
        ldr_gamma[mat->disp_txt] = 1;
        ldr_gamma[mat->bump_txt] = 1;
    }
    for (auto env : scn->environments) { ldr_gamma[env->ke_txt] = 2.2f; }

    // save images
    auto dirname = path_dirname(filename);
    for (auto txt : scn->textures) {
        if (txt->img.empty()) continue;
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        try {
            save_image(
                filename, txt->width, txt->height, txt->img, ldr_gamma.at(txt));
        } catch (std::exception&) {
            if (skip_missing) continue;
            throw;
        }
    }
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// GLTF CONVESION
// -----------------------------------------------------------------------------
namespace ygl {

static bool startswith(const std::string& str, const std::string& substr) {
    if (str.length() < substr.length()) return false;
    for (auto i = 0; i < substr.length(); i++)
        if (str[i] != substr[i]) return false;
    return true;
}

// Load a scene
scene* load_gltf_scene(
    const std::string& filename, bool load_textures, bool skip_missing) {
    auto gltf = load_json(filename);
    auto scn = new scene();

    // prepare parsing
    auto dirname = path_dirname(filename);

    // convert textures
    if (gltf.count("images")) {
        for (auto iid = 0; iid < gltf.at("images").size(); iid++) {
            auto& gimg = gltf.at("images").at(iid);
            auto txt = new texture();
            txt->name = gimg.value("name", ""s);
            txt->path = (startswith(gimg.value("uri", ""s), "data:")) ?
                            std::string("[glTF-inline].png") :
                            gimg.value("uri", ""s);
            scn->textures.push_back(txt);
        }
    }

    // load buffers
    auto bmap = std::vector<std::vector<byte>>();
    if (gltf.count("buffers")) {
        bmap.resize(gltf.at("buffers").size());
        for (auto bid = 0; bid < gltf.at("buffers").size(); bid++) {
            auto& gbuf = gltf.at("buffers").at(bid);
            auto& data = bmap.at(bid);
            auto uri = gbuf.value("uri", ""s);
            if (uri == "") continue;
            if (startswith(uri, "data:")) {
                // assume it is base64 and find ','
                auto pos = uri.find(',');
                if (pos == uri.npos) {
                    throw std::runtime_error("could not decode base64 data");
                }
                // decode
                auto data_char = base64_decode(uri.substr(pos + 1));
                data = std::vector<unsigned char>(
                    (unsigned char*)data_char.c_str(),
                    (unsigned char*)data_char.c_str() + data_char.length());
            } else {
                data = load_binary(fix_path(dirname + uri));
                if (data.empty()) {
                    throw std::runtime_error("could not load binary file " +
                                             fix_path(dirname + uri));
                }
            }
            if (gbuf.value("byteLength", -1) != data.size()) {
                throw std::runtime_error("mismatched buffer size");
            }
        }
    }

    // add a texture
    auto add_texture = [scn, &gltf](const json& ginfo) {
        if (!gltf.count("images") || !gltf.count("textures"))
            return (texture*)nullptr;
        if (ginfo.is_null() || ginfo.empty()) return (texture*)nullptr;
        if (ginfo.value("index", -1) < 0) return (texture*)nullptr;
        auto& gtxt = gltf.at("textures").at(ginfo.value("index", -1));
        if (gtxt.empty() || gtxt.value("source", -1) < 0)
            return (texture*)nullptr;
        auto txt = scn->textures.at(gtxt.value("source", -1));
        if (!gltf.count("samplers") || gtxt.value("sampler", -1) < 0)
            return txt;
        auto& gsmp = gltf.at("samplers").at(gtxt.value("sampler", -1));
        txt->clamp = gsmp.value("wrapS", ""s) == "ClampToEdge" ||
                     gsmp.value("wrapT", ""s) == "ClampToEdge";
        txt->scale = gsmp.value("scale", 1.0f) * gsmp.value("strength", 1.0f);
        return txt;
    };

    // convert materials
    if (gltf.count("materials")) {
        for (auto mid = 0; mid < gltf.at("materials").size(); mid++) {
            auto& gmat = gltf.at("materials").at(mid);
            auto mat = new material();
            mat->name = gmat.value("name", ""s);
            mat->ke = gmat.value("emissiveFactor", zero3f);
            if (gmat.count("emissiveTexture"))
                mat->ke_txt = add_texture(gmat.at("emissiveTexture"));
            if (gmat.count("extensions") &&
                gmat.at("extensions")
                    .count("KHR_materials_pbrSpecularGlossiness")) {
                mat->base_metallic = false;
                mat->gltf_textures = true;
                auto& gsg = gmat.at("extensions")
                                .at("KHR_materials_pbrSpecularGlossiness");
                auto kb = gsg.value("diffuseFactor", vec4f{1, 1, 1, 1});
                mat->kd = {kb.x, kb.y, kb.z};
                mat->op = kb.w;
                mat->ks = gsg.value("specularFactor", vec3f{1, 1, 1});
                mat->rs = 1 - gsg.value("glossinessFactor", 1.0f);
                if (gsg.count("diffuseTexture"))
                    mat->kd_txt = add_texture(gsg.at("diffuseTexture"));
                if (gsg.count("specularGlossinessTexture"))
                    mat->ks_txt =
                        add_texture(gsg.at("specularGlossinessTexture"));
                mat->rs_txt = mat->ks_txt;
            } else if (gmat.count("pbrMetallicRoughness")) {
                mat->base_metallic = true;
                mat->gltf_textures = true;
                auto& gmr = gmat.at("pbrMetallicRoughness");
                auto kb = gmr.value("baseColorFactor", vec4f{1, 1, 1, 1});
                mat->kd = {kb.x, kb.y, kb.z};
                mat->op = kb.w;
                auto km = gmr.value("metallicFactor", 1.0f);
                mat->ks = {km, km, km};
                mat->rs = gmr.value("roughnessFactor", 1.0f);
                if (gmr.count("baseColorTexture"))
                    mat->kd_txt = add_texture(gmr.at("baseColorTexture"));
                if (gmr.count("metallicRoughnessTexture"))
                    mat->ks_txt =
                        add_texture(gmr.at("metallicRoughnessTexture"));
                mat->rs_txt = mat->ks_txt;
            }
            if (gmat.count("occlusionTexture"))
                mat->occ_txt = add_texture(gmat.at("occlusionTexture"));
            if (gmat.count("normalTexture"))
                mat->norm_txt = add_texture(gmat.at("normalTexture"));
            mat->double_sided = gmat.value("doubleSided", false);
            scn->materials.push_back(mat);
        }
    }

    // get values from accessors
    auto accessor_values =
        [&gltf, &bmap](const json& gacc,
            bool normalize = false) -> std::vector<std::array<double, 4>> {
        auto gview = gltf.at("bufferViews").at(gacc.value("bufferView", -1));
        auto data = bmap.at(gview.value("buffer", -1)).data();
        auto offset =
            gacc.value("byteOffset", 0) + gview.value("byteOffset", 0);
        auto stride = gview.value("byteStride", 0);
        auto compTypeNum = gacc.value("componentType", 5123);
        auto count = gacc.value("count", -1);
        auto type = gacc.value("type", ""s);
        auto ncomp = 0;
        if (type == "SCALAR") ncomp = 1;
        if (type == "VEC2") ncomp = 2;
        if (type == "VEC3") ncomp = 3;
        if (type == "VEC4") ncomp = 4;
        auto compSize = 1;
        if (compTypeNum == 5122 || compTypeNum == 5123) { compSize = 2; }
        if (compTypeNum == 5124 || compTypeNum == 5125 || compTypeNum == 5126) {
            compSize = 4;
        }
        if (!stride) stride = compSize * ncomp;
        auto vals =
            std::vector<std::array<double, 4>>(count, {{0.0, 0.0, 0.0, 1.0}});
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
    auto meshes = std::vector<std::vector<std::pair<shape*, material*>>>();
    if (gltf.count("meshes")) {
        for (auto mid = 0; mid < gltf.at("meshes").size(); mid++) {
            auto& gmesh = gltf.at("meshes").at(mid);
            meshes.push_back({});
            auto sid = 0;
            for (auto& gprim : gmesh.value("primitives", json::array())) {
                if (!gprim.count("attributes")) continue;
                auto shp = new shape();
                shp->name = gmesh.value("name", ""s) +
                            ((sid) ? std::to_string(sid) : std::string());
                sid++;
                for (json::iterator gattr_it = gprim.at("attributes").begin();
                     gattr_it != gprim.at("attributes").end(); ++gattr_it) {
                    auto semantic = gattr_it.key();
                    auto& gacc =
                        gltf.at("accessors").at(gattr_it.value().get<int>());
                    auto vals = accessor_values(gacc);
                    if (semantic == "POSITION") {
                        shp->pos.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shp->pos.push_back({(float)vals[i][0],
                                (float)vals[i][1], (float)vals[i][2]});
                    } else if (semantic == "NORMAL") {
                        shp->norm.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shp->norm.push_back({(float)vals[i][0],
                                (float)vals[i][1], (float)vals[i][2]});
                    } else if (semantic == "TEXCOORD" ||
                               semantic == "TEXCOORD_0") {
                        shp->texcoord.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shp->texcoord.push_back(
                                {(float)vals[i][0], (float)vals[i][1]});
                    } else if (semantic == "COLOR" || semantic == "COLOR_0") {
                        shp->color.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shp->color.push_back(
                                {(float)vals[i][0], (float)vals[i][1],
                                    (float)vals[i][2], (float)vals[i][3]});
                    } else if (semantic == "TANGENT") {
                        shp->tangsp.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shp->tangsp.push_back(
                                {(float)vals[i][0], (float)vals[i][1],
                                    (float)vals[i][2], (float)vals[i][3]});
                        for (auto& t : shp->tangsp) t.w = -t.w;
                    } else if (semantic == "RADIUS") {
                        shp->radius.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shp->radius.push_back((float)vals[i][0]);
                    } else {
                        // ignore
                    }
                }
                // indices
                auto mode = gprim.value("mode", 4);
                if (!gprim.count("indices")) {
                    if (mode == 4) {
                        // triangles
                        shp->triangles.reserve(shp->pos.size() / 3);
                        for (auto i = 0; i < shp->pos.size() / 3; i++)
                            shp->triangles.push_back(
                                {i * 3 + 0, i * 3 + 1, i * 3 + 2});
                    } else if (mode == 6) {
                        // triangle fan
                        shp->triangles.reserve(shp->pos.size() - 2);
                        for (auto i = 2; i < shp->pos.size(); i++)
                            shp->triangles.push_back({0, i - 1, i});
                    } else if (mode == 5) {
                        // triangle strip
                        shp->triangles.reserve(shp->pos.size() - 2);
                        for (auto i = 2; i < shp->pos.size(); i++)
                            shp->triangles.push_back({i - 2, i - 1, i});
                    } else if (mode == 1) {
                        // lines
                        shp->lines.reserve(shp->pos.size() / 2);
                        for (auto i = 0; i < shp->pos.size() / 2; i++)
                            shp->lines.push_back({i * 2 + 0, i * 2 + 1});
                    } else if (mode == 2) {
                        // line loop
                        shp->lines.reserve(shp->pos.size());
                        for (auto i = 1; i < shp->pos.size(); i++)
                            shp->lines.push_back({i - 1, i});
                        shp->lines.back() = {(int)shp->pos.size() - 1, 0};
                    } else if (mode == 3) {
                        // line strip
                        shp->lines.reserve(shp->pos.size() - 1);
                        for (auto i = 1; i < shp->pos.size(); i++)
                            shp->lines.push_back({i - 1, i});
                    } else if (mode == -1 || mode == 0) {
                        // points
                        log_warning("points not supported");
                    } else {
                        throw std::runtime_error("unknown primitive type");
                    }
                } else {
                    auto indices = accessor_values(
                        gltf.at("accessors").at(gprim.value("indices", -1)),
                        false);
                    if (mode == 4) {
                        // triangles
                        shp->triangles.reserve(indices.size());
                        for (auto i = 0; i < indices.size() / 3; i++)
                            shp->triangles.push_back(
                                {(int)indices[i * 3 + 0][0],
                                    (int)indices[i * 3 + 1][0],
                                    (int)indices[i * 3 + 2][0]});
                    } else if (mode == 6) {
                        // triangle fan
                        shp->triangles.reserve(indices.size() - 2);
                        for (auto i = 2; i < indices.size(); i++)
                            shp->triangles.push_back({(int)indices[0][0],
                                (int)indices[i - 1][0], (int)indices[i][0]});
                    } else if (mode == 5) {
                        // triangle strip
                        shp->triangles.reserve(indices.size() - 2);
                        for (auto i = 2; i < indices.size(); i++)
                            shp->triangles.push_back({(int)indices[i - 2][0],
                                (int)indices[i - 1][0], (int)indices[i][0]});
                    } else if (mode == 1) {
                        // lines
                        shp->lines.reserve(indices.size() / 2);
                        for (auto i = 0; i < indices.size() / 2; i++)
                            shp->lines.push_back({(int)indices[i * 2 + 0][0],
                                (int)indices[i * 2 + 1][0]});
                    } else if (mode == 2) {
                        // line loop
                        shp->lines.reserve(indices.size());
                        for (auto i = 1; i < indices.size(); i++)
                            shp->lines.push_back(
                                {(int)indices[i - 1][0], (int)indices[i][0]});
                        shp->lines.back() = {
                            (int)indices[indices.size() - 1][0],
                            (int)indices[0][0]};
                    } else if (mode == 3) {
                        // line strip
                        shp->lines.reserve(indices.size() - 1);
                        for (auto i = 1; i < indices.size(); i++)
                            shp->lines.push_back(
                                {(int)indices[i - 1][0], (int)indices[i][0]});
                    } else if (mode == -1 || mode == 0) {
                        // points
                        log_warning("points not supported");
                    } else {
                        throw std::runtime_error("unknown primitive type");
                    }
                }
                auto mat = (gprim.count("material")) ?
                               scn->materials.at(gprim.value("material", -1)) :
                               nullptr;
                meshes.back().push_back({shp, mat});
                scn->shapes.push_back(shp);
            }
        }
    }

    // convert cameras
    if (gltf.count("cameras")) {
        for (auto cid = 0; cid < gltf.at("cameras").size(); cid++) {
            auto& gcam = gltf.at("cameras").at(cid);
            auto cam = new camera();
            cam->name = gcam.value("name", ""s);
            cam->ortho = gcam.value("type", ""s) == "orthographic";
            if (cam->ortho) {
                log_warning("orthographic not supported well");
                auto ortho = gcam.value("orthographic", json::object());
                set_camera_fovy(cam, ortho.value("ymag", 0.0f),
                    ortho.value("xmag", 0.0f) / ortho.value("ymag", 0.0f));
                cam->near = ortho.value("znear", 0.0f);
                cam->far = ortho.value("zfar", flt_max);
                cam->focus = flt_max;
                cam->aperture = 0;
            } else {
                auto persp = gcam.value("perspective", json::object());
                set_camera_fovy(cam, persp.value("yfov", 1.0f),
                    persp.value("aspectRatio", 1.0f));
                cam->near = persp.value("znear", 0.1f);
                cam->far =
                    persp.count("zfar") ? persp.value("zfar", 0.0f) : flt_max;
                cam->focus = flt_max;
                cam->aperture = 0;
            }
            scn->cameras.push_back(cam);
        }
    }

    // convert nodes
    if (gltf.count("nodes")) {
        for (auto nid = 0; nid < gltf.at("nodes").size(); nid++) {
            auto& gnde = gltf.at("nodes").at(nid);
            auto nde = new node();
            nde->name = gnde.value("name", ""s);
            if (gnde.count("camera"))
                nde->cam = scn->cameras[gnde.value("camera", 0)];
            nde->translation = gnde.value("translation", zero3f);
            nde->rotation = gnde.value("rotation", vec4f{0, 0, 0, 1});
            nde->scale = gnde.value("scale", vec3f{1, 1, 1});
            nde->frame = mat_to_frame(gnde.value("matrix", identity_mat4f));
            scn->nodes.push_back(nde);
        }

        // set up parent pointers
        for (auto nid = 0; nid < gltf.at("nodes").size(); nid++) {
            auto& gnde = gltf.at("nodes").at(nid);
            if (!gnde.count("children")) continue;
            auto nde = scn->nodes[nid];
            for (auto& cid : gnde.at("children"))
                scn->nodes[cid.get<int>()]->parent = nde;
        }

        // set up instances
        for (auto nid = 0; nid < gltf.at("nodes").size(); nid++) {
            auto& gnde = gltf.at("nodes").at(nid);
            if (!gnde.count("mesh")) continue;
            auto nde = scn->nodes[nid];
            auto& shps = meshes.at(gnde.value("mesh", 0));
            if (shps.empty()) continue;
            if (shps.size() == 1) {
                nde->ist = new instance();
                nde->ist->name = nde->name;
                nde->ist->shp = shps[0].first;
                nde->ist->mat = shps[0].second;
                scn->instances.push_back(nde->ist);
            } else {
                for (auto shp : shps) {
                    auto child = new node();
                    child->name = nde->name + "_" + shp.first->name;
                    child->parent = nde;
                    child->ist = new instance();
                    child->ist->name = child->name;
                    child->ist->shp = shp.first;
                    child->ist->mat = shp.second;
                    scn->instances.push_back(child->ist);
                }
            }
        }
    }

    // convert animations
    if (gltf.count("animations")) {
        for (auto& ganm : gltf.at("animations")) {
            auto aid = 0;
            auto sampler_map = std::unordered_map<vec2i, int>();
            for (auto& gchannel : ganm.at("channels")) {
                auto path_ =
                    gchannel.at("target").at("path").get<std::string>();
                auto path = -1;
                if (path_ == "translation") path = 0;
                if (path_ == "rotation") path = 1;
                if (path_ == "scale") path = 2;
                if (path_ == "weights") path = 3;
                if (sampler_map.find({gchannel.at("sampler").get<int>(),
                        path}) == sampler_map.end()) {
                    auto& gsampler = ganm.at("samplers")
                                         .at(gchannel.at("sampler").get<int>());
                    auto anm = new animation();
                    anm->name = (ganm.count("name") ? ganm.value("name", ""s) :
                                                      "anim") +
                                std::to_string(aid++);
                    anm->group = ganm.value("name", ""s);
                    auto input_view = accessor_values(
                        gltf.at("accessors").at(gsampler.value("input", -1)));
                    anm->times.resize(input_view.size());
                    for (auto i = 0; i < input_view.size(); i++)
                        anm->times[i] = input_view[i][0];
                    auto type = gsampler.value("interpolation", "LINEAR");
                    if (type == "LINEAR") anm->type = animation_type::linear;
                    if (type == "STEP") anm->type = animation_type::step;
                    if (type == "CUBICSPLINE")
                        anm->type = animation_type::bezier;
                    auto output_view = accessor_values(
                        gltf.at("accessors").at(gsampler.value("output", -1)));
                    switch (path) {
                        case 0: {  // translation
                            anm->translation.reserve(output_view.size());
                            for (auto i = 0; i < output_view.size(); i++)
                                anm->translation.push_back(
                                    {(float)output_view[i][0],
                                        (float)output_view[i][1],
                                        (float)output_view[i][2]});
                        } break;
                        case 1: {  // rotation
                            anm->rotation.reserve(output_view.size());
                            for (auto i = 0; i < output_view.size(); i++)
                                anm->rotation.push_back(
                                    {(float)output_view[i][0],
                                        (float)output_view[i][1],
                                        (float)output_view[i][2],
                                        (float)output_view[i][3]});
                        } break;
                        case 2: {  // scale
                            anm->scale.reserve(output_view.size());
                            for (auto i = 0; i < output_view.size(); i++)
                                anm->scale.push_back({(float)output_view[i][0],
                                    (float)output_view[i][1],
                                    (float)output_view[i][2]});
                        } break;
                        case 3: {  // weights
                            log_warning("weights not supported for now");
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
                            auto values = std::vector<float>();
                            values.reserve(output_view.size());
                            for (auto i = 0; i < output_view.size(); i++)
                                values.push_back(output_view.get(i));
                            anm->weights.resize(values.size() / ncomp);
                            for (auto i = 0; i < anm->weights.size(); i++) {
                                anm->weights[i].resize(ncomp);
                                for (auto j = 0; j < ncomp; j++)
                                    anm->weights[i][j] = values[i * ncomp + j];
                            }
                        }
#endif
                        } break;
                        default: {
                            throw std::runtime_error(
                                "should not have gotten here");
                        }
                    }
                    sampler_map[{gchannel.at("sampler").get<int>(), path}] =
                        (int)scn->animations.size();
                    scn->animations.push_back(anm);
                }
                scn->animations[sampler_map.at(
                                    {gchannel.at("sampler").get<int>(), path})]
                    ->targets.push_back(
                        scn->nodes
                            [(int)gchannel.at("target").at("node").get<int>()]);
            }
        }
    }

    // compute transforms and bbox
    update_transforms(scn, 0);
    update_bbox(scn);

    // fix elements
    auto mat = (material*)nullptr;
    for (auto ist : scn->instances) {
        if (ist->mat) continue;
        if (!mat) {
            mat = make_matte_material("<default>", {0.2f, 0.2f, 0.2f});
            scn->materials.push_back(mat);
        }
        ist->mat = mat;
    }
    for (auto cam : scn->cameras) {
        auto center = (scn->bbox.min + scn->bbox.max) / 2;
        auto dist = dot(-cam->frame.z, center - cam->frame.o);
        if (dist > 0) cam->focus = dist;
    }

    // skip textures if needed
    if (!load_textures) return scn;

    // set gamma
    auto ldr_gamma = std::unordered_map<texture*, float>{{nullptr, 1.0f}};
    for (auto txt : scn->textures) ldr_gamma[txt] = 2.2f;
    for (auto mat : scn->materials) {
        ldr_gamma[mat->ke_txt] = 2.2f;
        ldr_gamma[mat->kd_txt] = 2.2f;
        ldr_gamma[mat->ks_txt] = (mat->base_metallic) ? 1.0f : 2.2f;
        ldr_gamma[mat->kt_txt] = 2.2f;
        ldr_gamma[mat->rs_txt] = 1;
        ldr_gamma[mat->op_txt] = 1;
        ldr_gamma[mat->norm_txt] = 1;
        ldr_gamma[mat->disp_txt] = 1;
        ldr_gamma[mat->bump_txt] = 1;
    }
    for (auto env : scn->environments) { ldr_gamma[env->ke_txt] = 2.2f; }

    // load images
    for (auto txt : scn->textures) {
        auto filename = fix_path(dirname + txt->path);
        txt->img =
            load_image(filename, txt->width, txt->height, ldr_gamma.at(txt));
        if (txt->img.empty()) {
            if (skip_missing) continue;
            throw std::runtime_error("cannot laod image " + filename);
        }
    }

    // done
    return scn;
}

// Save gltf json
void save_gltf_scene(const std::string& filename, const scene* scn,
    bool save_textures, bool skip_missing) {
    // start creating json
    auto js = json::object();
    js["asset"]["version"] = "2.0";

    // prepare top level nodes
    if (!scn->cameras.empty()) js["cameras"] = json::array();
    if (!scn->textures.empty()) {
        js["textures"] = json::array();
        js["images"] = json::array();
    }
    if (!scn->materials.empty()) js["materials"] = json::array();
    if (!scn->shapes.empty()) {
        js["meshes"] = json::array();
        js["buffers"] = json::array();
        js["bufferViews"] = json::array();
        js["accessors"] = json::array();
    }
    if (!scn->instances.empty()) js["nodes"] = json::array();
    if (!scn->nodes.empty()) js["nodes"] = json::array();

    // convert cameras
    auto cmap = std::unordered_map<camera*, int>();
    for (auto cam : scn->cameras) {
        auto cjs = json();
        cjs["name"] = cam->name;
        if (!cam->ortho) {
            cjs["type"] = "perspective";
            cjs["perspective"]["aspectRatio"] = cam->width / cam->height;
            cjs["perspective"]["yfov"] = eval_camera_fovy(cam);
            cjs["perspective"]["znear"] = cam->near;
        } else {
            cjs["type"] = "orthographic";
            cjs["orthographic"]["xmag"] = cam->width / 2;
            cjs["orthographic"]["ymag"] = cam->height / 2;
            cjs["orthographic"]["zfar"] = cam->far;
            cjs["orthographic"]["znear"] = cam->near;
        }
        cmap[cam] = (int)js["cameras"].size();
        js["cameras"].push_back(cjs);
    }

    // textures
    auto tmap = std::unordered_map<texture*, int>();
    for (auto txt : scn->textures) {
        auto tjs = json(), ijs = json();
        tjs["source"] = (int)js["images"].size();
        ijs["uri"] = txt->path;
        js["images"].push_back(ijs);
        js["textures"].push_back(tjs);
        tmap[txt] = (int)js["textures"].size() - 1;
    }

    // material
    auto mmap = std::unordered_map<material*, int>();
    for (auto mat : scn->materials) {
        auto mjs = json();
        mjs["name"] = mat->name;
        mjs["doubleSided"] = mat->double_sided;
        if (mat->ke != zero3f) mjs["emissiveFactor"] = mat->ke;
        if (mat->ke_txt) mjs["emissiveTexture"]["index"] = tmap.at(mat->ke_txt);
        auto kd = vec4f{mat->kd.x, mat->kd.y, mat->kd.z, mat->op};
        if (mat->base_metallic) {
            auto mmjs = json();
            mmjs["baseColorFactor"] = kd;
            mmjs["metallicFactor"] = mat->ks.x;
            mmjs["roughnessFactor"] = mat->rs;
            if (mat->kd_txt)
                mmjs["baseColorTexture"]["index"] = tmap.at(mat->kd_txt);
            if (mat->ks_txt)
                mmjs["metallicRoughnessTexture"]["index"] =
                    tmap.at(mat->ks_txt);
            mjs["pbrMetallicRoughness"] = mmjs;
        } else {
            auto mmjs = json();
            mmjs["diffuseFactor"] = kd;
            mmjs["specularFactor"] = mat->ks;
            mmjs["glossinessFactor"] = 1 - mat->rs;
            if (mat->kd_txt)
                mmjs["diffuseTexture"]["index"] = tmap.at(mat->kd_txt);
            if (mat->ks_txt)
                mmjs["specularGlossinessTexture"]["index"] =
                    tmap.at(mat->ks_txt);
            mjs["extensions"]["KHR_materials_pbrSpecularGlossiness"] = mmjs;
        }
        if (mat->norm_txt)
            mjs["normalTexture"]["index"] = tmap.at(mat->norm_txt);
        if (mat->occ_txt)
            mjs["occlusionTexture"]["index"] = tmap.at(mat->occ_txt);
        js["materials"].push_back(mjs);
        mmap[mat] = (int)js["materials"].size() - 1;
    }

    // determine shape materials
    auto shape_mats = std::unordered_map<shape*, int>();
    for (auto ist : scn->instances)
        if (ist->mat) shape_mats[ist->shp] = mmap.at(ist->mat);

    // shapes
    auto smap = std::unordered_map<shape*, int>();
    for (auto shp : scn->shapes) {
        auto mjs = json(), bjs = json(), pjs = json();
        auto bid = js["buffers"].size();
        mjs["name"] = shp->name;
        mjs["primitives"] = json::array();
        bjs["name"] = shp->name;
        bjs["byteLength"] = 0;
        bjs["uri"] = replace_path_extension(shp->path, ".bin");
        auto mat_it = shape_mats.find(shp);
        if (mat_it != shape_mats.end()) pjs["material"] = mat_it->second;
        auto add_accessor = [&js, &bjs, bid](int count, std::string type,
                                bool indices = false) {
            auto bytes = count * 4;
            if (type == "VEC2") bytes *= 2;
            if (type == "VEC3") bytes *= 3;
            if (type == "VEC4") bytes *= 4;
            auto ajs = json(), vjs = json();
            vjs["buffer"] = bid;
            vjs["byteLength"] = bytes;
            vjs["byteOffset"] = bjs["byteLength"].get<int>();
            vjs["target"] = (!indices) ? 34962 : 34963;
            bjs["byteLength"] = bjs["byteLength"].get<int>() + bytes;
            ajs["bufferView"] = (int)js["bufferViews"].size();
            ajs["byteOffset"] = 0;
            ajs["componentType"] = (!indices) ? 5126 : 5125;
            ajs["count"] = count;
            ajs["type"] = type;
            js["accessors"].push_back(ajs);
            js["bufferViews"].push_back(vjs);
            return (int)js["accessors"].size() - 1;
        };
        auto nverts = (int)shp->pos.size();
        if (!shp->pos.empty())
            pjs["attributes"]["POSITION"] = add_accessor(nverts, "VEC3");
        if (!shp->norm.empty())
            pjs["attributes"]["NORMAL"] = add_accessor(nverts, "VEC3");
        if (!shp->texcoord.empty())
            pjs["attributes"]["TEXCOORD_0"] = add_accessor(nverts, "VEC2");
        if (!shp->color.empty())
            pjs["attributes"]["COLOR_0"] = add_accessor(nverts, "VEC4");
        if (!shp->radius.empty())
            pjs["attributes"]["RADIUS"] = add_accessor(nverts, "SCALAR");
        if (!shp->lines.empty()) {
            pjs["indices"] =
                add_accessor((int)shp->lines.size() * 2, "SCALAR", true);
            pjs["mode"] = 1;
        }
        if (!shp->triangles.empty()) {
            pjs["indices"] =
                add_accessor((int)shp->triangles.size() * 3, "SCALAR", true);
            pjs["mode"] = 4;
        }
        mjs["primitives"].push_back(pjs);
        js["meshes"].push_back(mjs);
        js["buffers"].push_back(bjs);
        smap[shp] = (int)js["meshes"].size() - 1;
    }

    // nodes
    auto nmap = std::unordered_map<node*, int>();
    for (auto nde : scn->nodes) {
        auto njs = json();
        njs["name"] = nde->name;
        njs["matrix"] = frame_to_mat(nde->frame);
        njs["translation"] = nde->translation;
        njs["rotation"] = nde->rotation;
        njs["scale"] = nde->scale;
        if (nde->cam) njs["camera"] = cmap.at(nde->cam);
        if (nde->ist) njs["mesh"] = smap.at(nde->ist->shp);
        if (!nde->children.empty()) {
            njs["children"] = json::array();
            for (auto c : nde->children) njs["children"].push_back(nmap.at(c));
        }
        js["nodes"].push_back(njs);
        nmap[nde] = (int)js["nodes"].size() - 1;
    }

    // animations not supported yet
    if (!scn->animations.empty()) log_warning("animation not supported yet");

    // nodes from instances
    if (scn->nodes.empty()) {
        for (auto cam : scn->cameras) {
            auto njs = json();
            njs["name"] = cam->name;
            njs["camera"] = cmap.at(cam);
            njs["matrix"] = frame_to_mat(cam->frame);
            js["nodes"].push_back(njs);
        }
        for (auto ist : scn->instances) {
            auto njs = json();
            njs["name"] = ist->name;
            njs["mesh"] = smap.at(ist->shp);
            njs["matrix"] = frame_to_mat(ist->frame);
            js["nodes"].push_back(njs);
        }
    }

    // save
    save_json(filename, js);

    // meshes
    auto dirname = path_dirname(filename);
    for (auto shp : scn->shapes) {
        if (shp->path == "") continue;
        auto filename = dirname + shp->path;
        filename = replace_path_extension(filename, ".bin");
        for (auto& c : filename)
            if (c == '\\') c = '/';
        try {
            auto fs = fopen(filename.c_str(), "wb");
            if (!fs)
                throw std::runtime_error("could not open file " + filename);
            fwrite(shp->pos.data(), 3 * 4, shp->pos.size(), fs);
            fwrite(shp->norm.data(), 3 * 4, shp->norm.size(), fs);
            fwrite(shp->texcoord.data(), 2 * 4, shp->pos.size(), fs);
            fwrite(shp->color.data(), 4 * 4, shp->color.size(), fs);
            fwrite(shp->radius.data(), 4, shp->radius.size(), fs);
            fwrite(shp->lines.data(), 2 * 4, shp->lines.size(), fs);
            fwrite(shp->triangles.data(), 3 * 4, shp->triangles.size(), fs);
            fclose(fs);
        } catch (std::exception&) {
            if (skip_missing) continue;
            throw;
        }
    }

    // skip textures if necessary
    if (!save_textures) return;

    // set gamma
    auto ldr_gamma = std::unordered_map<texture*, float>{{nullptr, 1.0f}};
    for (auto txt : scn->textures) ldr_gamma[txt] = 2.2f;
    for (auto mat : scn->materials) {
        ldr_gamma[mat->ke_txt] = 2.2f;
        ldr_gamma[mat->kd_txt] = 2.2f;
        ldr_gamma[mat->ks_txt] = (!mat->base_metallic) ? 2.2f : 1.0f;
        ldr_gamma[mat->kt_txt] = 2.2f;
        ldr_gamma[mat->rs_txt] = 1;
        ldr_gamma[mat->op_txt] = 1;
        ldr_gamma[mat->norm_txt] = 1;
        ldr_gamma[mat->disp_txt] = 1;
        ldr_gamma[mat->bump_txt] = 1;
    }
    for (auto env : scn->environments) { ldr_gamma[env->ke_txt] = 2.2f; }

    // save images
    for (auto txt : scn->textures) {
        if (txt->img.empty()) continue;
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        try {
            save_image(
                filename, txt->width, txt->height, txt->img, ldr_gamma.at(txt));
        } catch (std::exception&) {
            if (skip_missing) continue;
            throw;
        }
    }
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE IO
// -----------------------------------------------------------------------------
namespace ygl {

// Load ply mesh
void load_mesh(const std::string& filename, std::vector<vec2i>& lines,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<vec4f>& color, std::vector<float>& radius) {
    auto ext = path_extension(filename);
    if (ext == ".ply" || ext == ".PLY") {
        load_ply_mesh(
            filename, lines, triangles, pos, norm, texcoord, color, radius);
    } else if (ext == ".obj" || ext == ".OBJ") {
        load_obj_mesh(filename, lines, triangles, pos, norm, texcoord);
    } else {
        throw std::runtime_error("unsupported mesh extensions " + ext);
    }
}

// Save ply mesh
void save_mesh(const std::string& filename, const std::vector<vec2i>& lines,
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec2f>& texcoord,
    const std::vector<vec4f>& color, const std::vector<float>& radius,
    bool ascii) {
    auto ext = path_extension(filename);
    if (ext == ".ply" || ext == ".PLY") {
        save_ply_mesh(filename, lines, triangles, pos, norm, texcoord, color,
            radius, ascii);
    } else if (ext == ".obj" || ext == ".OBJ") {
        save_obj_mesh(filename, lines, triangles, pos, norm, texcoord);
    } else {
        throw std::runtime_error("unsupported mesh extensions " + ext);
    }
}

#if 1

// Load ply mesh
void load_ply_mesh(const std::string& filename, std::vector<vec2i>& lines,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<vec4f>& color, std::vector<float>& radius) {
    // clear data
    pos.clear();
    norm.clear();
    texcoord.clear();
    color.clear();
    radius.clear();
    lines.clear();
    triangles.clear();

    // ply header parsed
    enum property_type { uchar_type, int_type, float_type, int_list_type };
    struct property {
        std::string name = "";
        property_type type = float_type;
        std::vector<float> scalars;
        std::vector<std::array<int, 8>> lists;
    };
    struct element {
        std::string name = "";
        int count = 0;
        std::vector<property> props;
    };

    auto fs = fopen(filename.c_str(), "r");
    if (!fs) throw std::runtime_error("could not open file " + filename);

    // parse header
    char line[4096];
    auto ascii = false;
    auto elems = std::vector<element>();
    while (parse_getline(fs, line, sizeof(line))) {
        auto ss = line;
        auto cmd = parse_string(ss);
        if (cmd.empty()) continue;
        if (cmd == "ply") {
        } else if (cmd == "comment") {
        } else if (cmd == "format") {
            auto fmt = parse_string(ss);
            if (fmt != "ascii" && fmt != "binary_little_endian")
                throw std::runtime_error("format not supported");
            ascii = fmt == "ascii";
        } else if (cmd == "element") {
            auto elem = element();
            elem.name = parse_string(ss);
            elem.count = parse_int(ss);
            elems.push_back(elem);
        } else if (cmd == "property") {
            auto prop = property();
            auto type = parse_string(ss);
            if (type == "list") {
                auto count_type = parse_string(ss);
                auto elem_type = parse_string(ss);
                if (count_type != "uchar" && count_type != "uint8")
                    throw std::runtime_error("unsupported ply list type");
                if (elem_type != "int")
                    throw std::runtime_error("unsupported ply list type");
                prop.type = int_list_type;
            } else if (type == "float") {
                prop.type = float_type;
            } else if (type == "uchar" || type == "uint8") {
                prop.type = uchar_type;
            } else if (type == "int") {
                prop.type = int_type;
            } else {
                throw std::runtime_error("unsupported ply type");
            }
            prop.name = parse_string(ss);
            prop.scalars.resize(elems.back().count);
            if (prop.type == int_list_type)
                prop.lists.resize(elems.back().count);
            elems.back().props.push_back(prop);
        } else if (cmd == "end_header") {
            break;
        } else {
            throw std::runtime_error("command not supported " + cmd);
        }
    }

    // parse content
    for (auto& elem : elems) {
        for (auto vid = 0; vid < elem.count; vid++) {
            auto ss = (char*)nullptr;
            if (ascii) {
                if (!parse_getline(fs, line, sizeof(line)))
                    throw std::runtime_error("error reading ply");
                ss = line;
            }
            for (auto pid = 0; pid < elem.props.size(); pid++) {
                auto& prop = elem.props[pid];
                if (prop.type == float_type) {
                    prop.scalars[vid] =
                        (ascii) ? parse_float(ss) : read_float(fs);
                } else if (prop.type == int_type) {
                    prop.scalars[vid] = (ascii) ? parse_int(ss) : read_int(fs);
                } else if (prop.type == uchar_type) {
                    prop.scalars[vid] = (ascii) ? parse_uchar(ss) / 255.0 :
                                                  read_uchar(fs) / 255.0;
                } else if (prop.type == int_list_type) {
                    prop.scalars[vid] =
                        (ascii) ? parse_uchar(ss) : read_uchar(fs);
                    for (auto i = 0; i < (int)prop.scalars[vid]; i++)
                        prop.lists[vid][i] =
                            (ascii) ? parse_int(ss) : read_int(fs);
                } else {
                    throw std::runtime_error("unsupported ply type");
                }
            }
        }
    }

    // copy vertex data
    for (auto& elem : elems) {
        if (elem.name != "vertex") continue;
        auto count = elem.count;
        for (auto& prop : elem.props) {
            auto vals = prop.scalars.data();
            auto copy_floats = [vals, count](auto& vert, const auto& def,
                                   int stride, int offset) {
                if (vert.size() != count) vert.resize(count, def);
                auto dst = (float*)vert.data();
                for (auto i = 0; i < count; i++)
                    dst[i * stride + offset] = vals[i];
            };
            if (prop.name == "x") copy_floats(pos, zero3f, 3, 0);
            if (prop.name == "y") copy_floats(pos, zero3f, 3, 1);
            if (prop.name == "z") copy_floats(pos, zero3f, 3, 2);
            if (prop.name == "nx") copy_floats(norm, zero3f, 3, 0);
            if (prop.name == "ny") copy_floats(norm, zero3f, 3, 1);
            if (prop.name == "nz") copy_floats(norm, zero3f, 3, 2);
            if (prop.name == "u") copy_floats(texcoord, zero2f, 2, 0);
            if (prop.name == "v") copy_floats(texcoord, zero2f, 2, 1);
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

    // copy triangle data
    for (auto& elem : elems) {
        if (elem.name != "face") continue;
        auto count = elem.count;
        for (auto& prop : elem.props) {
            if (prop.name == "vertex_indices") {
                for (auto fid = 0; fid < count; fid++) {
                    auto& list = prop.lists[fid];
                    for (auto i = 2; i < (int)prop.scalars[fid]; i++)
                        triangles.push_back({list[0], list[i - 1], list[i]});
                }
            }
        }
    }

    fclose(fs);
}

#else

// Load ply mesh
void load_ply_mesh(const std::string& filename, std::vector<vec2i>& lines,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<vec4f>& color, std::vector<float>& radius) {
    // clear data
    pos.clear();
    norm.clear();
    texcoord.clear();
    color.clear();
    radius.clear();
    lines.clear();
    triangles.clear();

    // open ply and read its header
    auto error_cb = [](p_ply ply, const char* message) {};
    auto ply = ply_open(filename.c_str(), error_cb, 0, nullptr);
    if (!ply) throw std::runtime_error("could not open ply file " + filename);
    if (!ply_read_header(ply))
        throw std::runtime_error("could not read ply file " + filename);

    // read polygons
    static constexpr int poly_size = 8;
    auto polys = std::vector<std::array<int, poly_size>>();
    auto plines = std::vector<std::array<int, poly_size>>();

    // check contained data
    auto elem = (p_ply_element) nullptr;
    while ((elem = ply_get_next_element(ply, elem)) != nullptr) {
        auto elem_name_buf = (const char*)nullptr;
        auto elem_num = (long)0;
        ply_get_element_info(elem, &elem_name_buf, &elem_num);
        auto elem_name = std::string(elem_name_buf);
        if (elem_name == "vertex") {
            auto prop = (p_ply_property) nullptr;
            while ((prop = ply_get_next_property(elem, prop)) != nullptr) {
                auto prop_name_buf = (const char*)nullptr;
                ply_get_property_info(
                    prop, &prop_name_buf, nullptr, nullptr, nullptr);
                auto prop_name = std::string(prop_name_buf);
                if (prop_name == "x") pos.assign(elem_num, zero3f);
                if (prop_name == "nx") norm.assign(elem_num, zero3f);
                if (prop_name == "u" || elem_name == "s" ||
                    elem_name == "texture_u" || elem_name == "texture_s")
                    texcoord.assign(elem_num, zero2f);
                if (prop_name == "red") color.assign(elem_num, {0, 0, 0, 1});
                if (prop_name == "radius") radius.assign(elem_num, 0);
            }
        } else if (elem_name == "face") {
            auto prop = (p_ply_property) nullptr;
            while ((prop = ply_get_next_property(elem, prop)) != nullptr) {
                auto prop_name_buf = (const char*)nullptr;
                ply_get_property_info(
                    prop, &prop_name_buf, nullptr, nullptr, nullptr);
                auto prop_name = std::string(prop_name_buf);
                if (prop_name == "vertex_indices")
                    polys.assign(elem_num, {{0, 0, 0, 0, 0, 0, 0, 0}});
            }
        } else if (elem_name == "line") {
            auto prop = (p_ply_property) nullptr;
            while ((prop = ply_get_next_property(elem, prop)) != nullptr) {
                auto prop_name_buf = (const char*)nullptr;
                ply_get_property_info(
                    prop, &prop_name_buf, nullptr, nullptr, nullptr);
                auto prop_name = std::string(prop_name_buf);
                if (prop_name == "vertex_indices")
                    plines.assign(elem_num, {{0, 0, 0, 0, 0, 0, 0, 0}});
            }
        }
    }

    // vertex data callback
    auto vert_cb = [](p_ply_argument arg) {
        auto data = (float*)nullptr;
        auto idx = (long)0, flags = (long)0;
        ply_get_argument_user_data(arg, (void**)&data, &flags);
        ply_get_argument_element(arg, nullptr, &idx);
        auto stride = (flags & 0x0F0) >> 4;
        auto offset = flags & 0x00F;
        data[idx * stride + offset] = (float)ply_get_argument_value(arg);
        return 1;
    };

    // face data callback
    auto face_cb = [](p_ply_argument arg) {
        auto data = (int*)nullptr;
        auto idx = (long)0, flags = (long)0, len = (long)0, vidx = (long)0;
        ply_get_argument_user_data(arg, (void**)&data, &flags);
        ply_get_argument_element(arg, nullptr, &idx);
        ply_get_argument_property(arg, nullptr, &len, &vidx);
        if (len >= poly_size) throw std::runtime_error("bad face length");
        data[idx * poly_size + (vidx + 1)] = (int)ply_get_argument_value(arg);
        return 1;
    };

    // set up vertex callbacks
    ply_set_read_cb(ply, "vertex", "x", vert_cb, pos.data(), 0x30);
    ply_set_read_cb(ply, "vertex", "y", vert_cb, pos.data(), 0x31);
    ply_set_read_cb(ply, "vertex", "z", vert_cb, pos.data(), 0x32);
    ply_set_read_cb(ply, "vertex", "nx", vert_cb, norm.data(), 0x30);
    ply_set_read_cb(ply, "vertex", "ny", vert_cb, norm.data(), 0x31);
    ply_set_read_cb(ply, "vertex", "nz", vert_cb, norm.data(), 0x32);
    ply_set_read_cb(ply, "vertex", "u", vert_cb, texcoord.data(), 0x20);
    ply_set_read_cb(ply, "vertex", "v", vert_cb, texcoord.data(), 0x21);
    ply_set_read_cb(ply, "vertex", "s", vert_cb, texcoord.data(), 0x20);
    ply_set_read_cb(ply, "vertex", "t", vert_cb, texcoord.data(), 0x21);
    ply_set_read_cb(ply, "vertex", "texture_u", vert_cb, texcoord.data(), 0x20);
    ply_set_read_cb(ply, "vertex", "texture_v", vert_cb, texcoord.data(), 0x21);
    ply_set_read_cb(ply, "vertex", "texture_s", vert_cb, texcoord.data(), 0x20);
    ply_set_read_cb(ply, "vertex", "texture_t", vert_cb, texcoord.data(), 0x21);
    ply_set_read_cb(ply, "vertex", "red", vert_cb, color.data(), 0x40);
    ply_set_read_cb(ply, "vertex", "green", vert_cb, color.data(), 0x41);
    ply_set_read_cb(ply, "vertex", "blue", vert_cb, color.data(), 0x42);
    ply_set_read_cb(ply, "vertex", "alpha", vert_cb, color.data(), 0x43);
    ply_set_read_cb(ply, "vertex", "radius", vert_cb, radius.data(), 0x10);

    // set up triangle and line callbacks
    ply_set_read_cb(ply, "face", "vertex_indices", face_cb, polys.data(), 0x10);
    ply_set_read_cb(
        ply, "line", "vertex_indices", face_cb, plines.data(), 0x10);

    // read file
    if (!ply_read(ply))
        throw std::runtime_error("error reading ply file " + filename);
    ply_close(ply);

    // convert polygons to triangles
    auto ntriangles = 0;
    for (auto& poly : polys) ntriangles += poly[0] - 2;
    triangles.resize(ntriangles);
    auto ti = 0;
    for (auto& poly : polys) {
        for (auto i = 3; i <= poly[0]; i++)
            triangles[ti++] = {poly[1], poly[i - 1], poly[i]};
    }

    // convert polylines to lines
    auto nlines = 0;
    for (auto& poly : plines) nlines += poly[0] - 1;
    lines.resize(nlines);
    auto li = 0;
    for (auto& pline : plines) {
        for (auto i = 1; i <= pline[0]; i++)
            lines[li++] = {pline[i - 1], pline[i]};
    }
}

#endif

// Save ply mesh
void save_ply_mesh(const std::string& filename, const std::vector<vec2i>& lines,
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec2f>& texcoord,
    const std::vector<vec4f>& color, const std::vector<float>& radius,
    bool ascii) {
    auto f = fopen(filename.c_str(), "w");
    if (!f) throw std::runtime_error("cannot save file " + filename);

    // header
    fprintf(f, "ply\n");
    fprintf(f, "format %s 1.0\n", (ascii) ? "ascii" : "binary_little_endian");
    fprintf(f, "element vertex %d\n", (int)pos.size());
    if (!pos.empty())
        fprintf(f, "property float x\nproperty float y\nproperty float z\n");
    if (!norm.empty())
        fprintf(f, "property float nx\nproperty float ny\nproperty float nz\n");
    if (!texcoord.empty()) fprintf(f, "property float u\nproperty float v\n");
    if (!color.empty())
        fprintf(f,
            "property float red\nproperty float green\nproperty float "
            "blue\nproperty float alpha\n");
    if (!radius.empty()) fprintf(f, "property float radius\n");
    if (!triangles.empty()) {
        fprintf(f, "element face %d\n", (int)triangles.size());
        fprintf(f, "property list uchar int vertex_indices\n");
    }
    if (!lines.empty()) {
        fprintf(f, "element line %d\n", (int)triangles.size());
        fprintf(f, "property list uchar int vertex_indices\n");
    }
    fprintf(f, "end_header\n");

    // body
    if (ascii) {
        // write vertex data
        for (auto i = 0; i < pos.size(); i++) {
            if (!pos.empty())
                fprintf(f, "%g %g %g ", pos[i].x, pos[i].y, pos[i].z);
            if (!norm.empty())
                fprintf(f, "%g %g %g ", norm[i].x, norm[i].y, norm[i].z);
            if (!texcoord.empty())
                fprintf(f, "%g %g ", texcoord[i].x, texcoord[i].y);
            if (!color.empty())
                fprintf(f, "%g %g %g %g ", color[i].x, color[i].y, color[i].z,
                    color[i].w);
            if (!radius.empty()) fprintf(f, "%g ", radius[i]);
            fprintf(f, "\n");
        }

        // write face data
        for (auto i = 0; i < triangles.size(); i++) {
            fprintf(f, "3 %d %d %d\n", triangles[i].x, triangles[i].y,
                triangles[i].z);
        }
        for (auto i = 0; i < lines.size(); i++) {
            fprintf(f, "2 %d %d\n", lines[i].x, lines[i].y);
        }
    } else {
        // write vertex data
        for (auto i = 0; i < pos.size(); i++) {
            if (!pos.empty()) fwrite(&pos[i], 4, 3, f);
            if (!norm.empty()) fwrite(&norm[i], 4, 3, f);
            if (!texcoord.empty()) fwrite(&texcoord[i], 4, 2, f);
            if (!color.empty()) fwrite(&color[i], 4, 4, f);
            if (!radius.empty()) fwrite(&radius[i], 4, 1, f);
        }

        // write face data
        for (auto i = 0; i < triangles.size(); i++) {
            auto n = (byte)3;
            fwrite(&n, 1, 1, f);
            fwrite(&triangles[i], 4, 3, f);
        }
        for (auto i = 0; i < lines.size(); i++) {
            auto n = (byte)3;
            fwrite(&n, 1, 1, f);
            fwrite(&lines[i], 4, 2, f);
        }
    }

    // done
    fclose(f);
}

// Load ply mesh
void load_obj_mesh(const std::string& filename, std::vector<vec2i>& lines,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    bool flip_texcoord) {
    // open file
    auto fs = fopen(filename.c_str(), "rt");
    if (!fs) throw std::runtime_error("cannot open filename " + filename);

    pos.clear();
    norm.clear();
    texcoord.clear();
    lines.clear();
    triangles.clear();

    // vertex maps
    auto vert_map = std::unordered_map<vec3i, int>();

    // read the file line by line
    char line[4096];
    while (parse_getline(fs, line, sizeof(line))) {
        // prepare to parse
        auto ss = line;
        parse_skipws(ss);
        parse_remove_comment(ss, '#');

        // get command
        auto cmd = parse_string(ss);
        if (cmd == "") continue;

        // possible token values
        if (cmd == "v") {
            pos.push_back(parse_vec3f(ss));
        } else if (cmd == "vn") {
            norm.push_back(parse_vec3f(ss));
        } else if (cmd == "vt") {
            texcoord.push_back(parse_vec2f(ss));
            if (flip_texcoord) texcoord.back().y = 1 - texcoord.back().y;
        } else if (cmd == "f" || cmd == "l") {
            auto num = 0;
            vec3i verts[128];
            int vids[128];
            auto vert_size =
                vec3i{(int)pos.size(), (int)texcoord.size(), (int)norm.size()};
            // elem.material = (int)oobj->materials.size() - 1;
            parse_skipws(ss);
            while (*ss) {
                verts[num++] = parse_objvert(ss, vert_size);
                parse_skipws(ss);
            }
            for (auto i = 0; i < num; i++) {
                auto it = vert_map.find(verts[i]);
                if (it == vert_map.end()) {
                    auto nverts = (int)pos.size();
                    vert_map.insert(it, {verts[i], nverts});
                    vids[i] = nverts;
                    if (verts[i].x >= 0) pos.push_back(pos.at(verts[i].x));
                    if (verts[i].y >= 0)
                        texcoord.push_back(texcoord.at(verts[i].y));
                    if (verts[i].z >= 0) norm.push_back(norm.at(verts[i].z));
                } else {
                    vids[i] = it->second;
                }
            }
            if (cmd == "f") {
                for (auto i = 2; i < num; i++)
                    triangles.push_back({vids[0], vids[i - 1], vids[i]});
            }
            if (cmd == "l") {
                for (auto i = 1; i < num; i++)
                    lines.push_back({vids[i - 1], vids[i]});
            }
        }
    }

    // close file
    fclose(fs);
}

// Load ply mesh
void save_obj_mesh(const std::string& filename, const std::vector<vec2i>& lines,
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec2f>& texcoord,
    bool flip_texcoord) {
    auto f = fopen(filename.c_str(), "wt");
    if (!f) throw std::runtime_error("cannot save file " + filename);
    for (auto& p : pos) fprintf(f, "v %g %g %g\n", p.x, p.y, p.z);
    for (auto& n : norm) fprintf(f, "vn %g %g %g\n", n.x, n.y, n.z);
    if (flip_texcoord)
        for (auto& t : texcoord) fprintf(f, "vt %g %g\n", t.x, 1 - t.y);
    else
        for (auto& t : texcoord) fprintf(f, "vt %g %g\n", t.x, t.y);
    if (!texcoord.empty() && !norm.empty()) {
        for (auto& t : triangles)
            fprintf(f, "f %d/%d/%d %d/%d/%d %d/%d/%d\n", t.x + 1, t.x + 1,
                t.x + 1, t.y + 1, t.y + 1, t.y + 1, t.z + 1, t.z + 1, t.z + 1);
        for (auto& l : lines)
            fprintf(f, "l %d/%d/%d %d/%d/%d\n", l.x + 1, l.x + 1, l.x + 1,
                l.y + 1, l.y + 1, l.y + 1);
    } else if (!texcoord.empty() && norm.empty()) {
        for (auto& t : triangles)
            fprintf(f, "f %d/%d %d/%d %d/%d\n", t.x + 1, t.x + 1, t.y + 1,
                t.y + 1, t.z + 1, t.z + 1);
        for (auto& l : lines)
            fprintf(f, "l %d/%d %d/%d\n", l.x + 1, l.x + 1, l.y + 1, l.y + 1);
    } else if (texcoord.empty() && !norm.empty()) {
        for (auto& t : triangles)
            fprintf(f, "f %d//%d %d//%d %d//%d\n", t.x + 1, t.x + 1, t.y + 1,
                t.y + 1, t.z + 1, t.z + 1);
        for (auto& l : lines)
            fprintf(f, "l %d//%d %d//%d\n", l.x + 1, l.x + 1, l.y + 1, l.y + 1);
    } else {
        for (auto& t : triangles)
            fprintf(f, "f %d %d %d\n", t.x + 1, t.y + 1, t.z + 1);
        for (auto& l : lines) fprintf(f, "l %d %d\n", l.x + 1, l.y + 1);
    }
    fclose(f);
}

}  // namespace ygl
