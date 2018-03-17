//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//

//
// Converts scenes from https://benedikt-bitterli.me/resources/ to Yocto/GL.
//

#include "../../yocto/yocto_gl.h"
#include "ext/json.hpp"

#define YTUNGSTEN_MITSUBA_OBJ 0
#define YTUNSTEN_ORIGINAL_XFORM 0

using namespace ygl;
using namespace nlohmann;
using namespace std::literals;

const float env_distance = 100.0f;

vec3f to_vec3f(const json& js) {
    if (js.is_number())
        return {js.get<float>(), js.get<float>(), js.get<float>()};
    if (js.is_array())
        return {js.at(0).get<float>(), js.at(1).get<float>(),
            js.at(2).get<float>()};
    throw std::runtime_error("bad json");
};

#if YTUNSTEN_ORIGINAL_XFORM
// from tungsten source
mat3f rotYXZ(const vec3f& rot) {
    auto r = rot * pif / 180.0f;
    float c[] = {std::cos(r.x), std::cos(r.y), std::cos(r.z)};
    float s[] = {std::sin(r.x), std::sin(r.y), std::sin(r.z)};

    return transpose(mat3f{{c[1] * c[2] - s[1] * s[0] * s[2],
                               -c[1] * s[2] - s[1] * s[0] * c[2], -s[1] * c[0]},
        {c[0] * s[2], c[0] * c[2], -s[0]},
        {s[1] * c[2] + c[1] * s[0] * s[2], -s[1] * s[2] + c[1] * s[0] * c[2],
            c[1] * c[0]}});
}

vec3f randomOrtho(const vec3f& a) {
    vec3f res;
    if (std::abs(a.x) > std::abs(a.y))
        res = vec3f(0.0f, 1.0f, 0.0f);
    else
        res = vec3f(1.0f, 0.0f, 0.0f);
    return normalize(cross(a, res));
}

void gramSchmidt(vec3f& a, vec3f& b, vec3f& c) {
    a = normalize(a);
    b -= a * dot(a, b);
    if (dot(b, b) < 1e-5)
        b = randomOrtho(a);
    else
        b = normalize(b);

    c -= a * dot(a, c);
    c -= b * dot(b, c);
    if (dot(c, c) < 1e-5)
        c = cross(a, b);
    else
        c = normalize(c);
}

frame3f xform_to_frame(const json& js, float* focus = nullptr) {
    if (js.is_array()) {
        auto m = identity_mat4f;
        for (auto i = 0; i < 16; i++) ((float*)&m)[i] = js.at(i).get<float>();
        return mat_to_frame(transpose(m));
    } else if (js.is_object()) {
        auto x = vec3f(1.0f, 0.0f, 0.0f);
        auto y = vec3f(0.0f, 1.0f, 0.0f);
        auto z = vec3f(0.0f, 0.0f, 1.0f);

        auto pos = zero3f;
        if (js.count("position")) pos = to_vec3f(js.at("position"));

        bool explicitX = false, explicitY = false, explicitZ = false;

        if (js.count("look_at")) {
            auto lookAt = to_vec3f(js.at("look_at"));
            z = lookAt - pos;
            if (focus) *focus = length(pos - lookAt);
            explicitZ = true;
        }

        if (js.count("up")) {
            y = to_vec3f(js.at("up"));
            explicitY = true;
        }

        if (js.count("x_axis")) {
            x = to_vec3f(js.at("x_axis"));
            explicitX = true;
        }
        if (js.count("y_axis")) {
            y = to_vec3f(js.at("y_axis"));
            explicitY = true;
        }
        if (js.count("z_axis")) {
            z = to_vec3f(js.at("z_axis"));
            explicitZ = true;
        }

        int id =
            (explicitZ ? 4 : 0) + (explicitY ? 2 : 0) + (explicitX ? 1 : 0);
        switch (id) {
            case 0: gramSchmidt(z, y, x); break;
            case 1: gramSchmidt(x, z, y); break;
            case 2: gramSchmidt(y, z, x); break;
            case 3: gramSchmidt(y, x, z); break;
            case 4: gramSchmidt(z, y, x); break;
            case 5: gramSchmidt(z, x, y); break;
            case 6: gramSchmidt(z, y, x); break;
            case 7: gramSchmidt(z, y, x); break;
        }

        if (dot(cross(x, y), z) < 0.0f) {
            if (!explicitX)
                x = -x;
            else if (!explicitY)
                y = -y;
            else
                z = -z;
        }

        if (js.count("scale")) {
            auto scale = to_vec3f(js.at("scale"));
            x *= scale.x;
            y *= scale.y;
            z *= scale.z;
        }

        if (js.count("rotation")) {
            auto rot = to_vec3f(js.at("rotation"));
            auto xform = rotYXZ(rot);
            x = xform * x;
            y = xform * y;
            z = xform * z;
        }

        return {x, y, z, pos};
    } else {
        throw std::runtime_error(
            "Parameter has wrong type: Expecting a matrix value here");
    }
}

#else

// from tungsten source
frame3f rotation_yxz_frame(const vec3f& r) {
    float c[] = {std::cos(r.x), std::cos(r.y), std::cos(r.z)};
    float s[] = {std::sin(r.x), std::sin(r.y), std::sin(r.z)};

    auto m = mat3f{{c[1] * c[2] - s[1] * s[0] * s[2],
                       -c[1] * s[2] - s[1] * s[0] * c[2], -s[1] * c[0]},
        {c[0] * s[2], c[0] * c[2], -s[0]},
        {s[1] * c[2] + c[1] * s[0] * s[2], -s[1] * s[2] + c[1] * s[0] * c[2],
            c[1] * c[0]}};
    m = transpose(m);
    return {m.x, m.y, m.z, zero3f};
}

std::pair<frame3f, vec3f> xform_to_scaled_frame(
    const json& js, float* focus = nullptr) {
    if (js.is_array()) {
        auto m = identity_mat4f;
        for (auto i = 0; i < 16; i++) ((float*)&m)[i] = js.at(i).get<float>();
        auto f = mat_to_frame(m);
        auto scl = vec3f{length(f.x), length(f.y), length(f.z)};
        return {{normalize(f.x), normalize(f.y), normalize(f.z), f.o}, scl};
    } else if (js.count("look_at")) {
        auto from = to_vec3f(js.at("position"));
        auto to = to_vec3f(js.at("look_at"));
        auto up = to_vec3f(js.at("up"));
        if (focus) *focus = length(from - to);
        auto scl =
            (js.count("scale")) ? to_vec3f(js.at("scale")) : vec3f{1, 1, 1};
        return {lookat_frame(from, to, up), scl};
    } else {
        auto pos =
            (js.count("position")) ? to_vec3f(js.at("position")) : zero3f;
        auto rot =
            (js.count("rotation")) ? to_vec3f(js.at("rotation")) : zero3f;
        auto scl =
            (js.count("scale")) ? to_vec3f(js.at("scale")) : vec3f{1, 1, 1};
        rot *= pif / 180;
        // return translation_frame(pos) * rotation_frame(vec3f{0, 1, 0}, rot.x)
        // *
        //        rotation_frame(vec3f{1, 0, 0}, rot.y) *
        //        rotation_frame(vec3f{0, 0, 1}, rot.z) * scaling_frame(scl);
        // return translation_frame(pos) * rotation_frame(vec3f{0, 0, 1}, rot.x)
        // *
        //        rotation_frame(vec3f{1, 0, 0}, rot.y) *
        //        rotation_frame(vec3f{0, 1, 0}, rot.z) * scaling_frame(scl);
        return {translation_frame(pos) * rotation_yxz_frame(rot), scl};
    }
    log_error("bad xform");
    return {identity_frame3f, vec3f{1, 1, 1}};
}

#endif

shape* make_shape(const std::string& type, const json& js) {
    static auto tcount = std::unordered_map<std::string, int>();
    auto shp = new shape();
    shp->name = type + std::to_string(tcount[type]++);
    if (type == "quad") {
        shp->pos = {{-0.5f, 0.0f, -0.5f}, {0.5f, 0.0f, -0.5f},
            {0.5f, 0.0f, 0.5f}, {-0.5f, 0.0f, 0.5f}};
        shp->norm = {{0, 1, 0}, {0, 1, 0}, {0, 1, 0}, {0, 1, 0}};
        shp->texcoord = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
        shp->triangles = {{0, 1, 2}, {0, 2, 3}};
    } else if (type == "cube") {
        make_cube(shp->quads, shp->pos, shp->norm, shp->texcoord, {1, 1, 1},
            {1, 1, 1}, {1, 1, 1});
    } else if (type == "sphere") {
        make_sphere(shp->quads, shp->pos, shp->norm, shp->texcoord, {64, 32}, 2,
            {1, 1});
    } else if (type == "infinite_sphere") {
        make_sphere(shp->quads, shp->pos, shp->norm, shp->texcoord, {128, 64},
            1000, {1, 1});
        for (auto& n : shp->norm) n = -n;
        for (auto& t : shp->triangles) std::swap(t.y, t.z);
        for (auto& q : shp->quads) std::swap(q.y, q.w);
    } else if (type == "infinite_sphere_cap") {
        // flipped
        shp->pos = {{-0.5f, 0.0f, -0.5f}, {0.5f, 0.0f, -0.5f},
            {0.5f, 0.0f, 0.5f}, {-0.5f, 0.0f, 0.5f}};
        shp->norm = {{0, -1, 0}, {0, -1, 0}, {0, -1, 0}, {0, -1, 0}};
        shp->texcoord = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
        shp->triangles = {{0, 2, 1}, {0, 3, 2}};
        auto angle = js.at("cap_angle").get<float>() * pif / 180;
        for (auto& p : shp->pos) p = vec3f{0, env_distance, 0} + p * 2 * env_distance * std::sin(angle);
    } else
        log_error("unknow type " + type);
    return shp;
}

shape* load_wo3(const std::string& filename) {
    auto fs = fopen(filename.c_str(), "rb");
    if (!fs) throw std::runtime_error("cannot open file " + filename);

    // uint64 numVerts, numTris;
    // FileUtils::streamRead(stream, numVerts);
    // verts.resize(size_t(numVerts));
    // FileUtils::streamRead(stream, verts);
    // FileUtils::streamRead(stream, numTris);
    // tris.resize(size_t(numTris));
    // FileUtils::streamRead(stream, tris);

    struct vertex {
        vec3f pos;
        vec3f norm;
        vec2f texcoord;
    };

    auto nverts = (uint64_t)0;
    fread(&nverts, sizeof(uint64_t), 1, fs);
    auto verts = std::vector<vertex>(nverts);
    fread(verts.data(), sizeof(vertex), nverts, fs);
    auto ntriangles = (uint64_t)0;
    fread(&ntriangles, sizeof(uint64_t), 1, fs);
    auto triangles = std::vector<vec4i>(ntriangles);
    fread(triangles.data(), sizeof(vec4i), ntriangles, fs);

    fclose(fs);

    auto shp = new shape();
    shp->name = path_basename(filename);
    for (auto& vert : verts) {
        shp->pos.push_back(vert.pos);
        shp->norm.push_back(vert.norm);
        shp->texcoord.push_back(vert.texcoord);
    }
    for (auto& t : triangles) { shp->triangles.push_back({t.x, t.y, t.z}); }
    for (auto& uv : shp->texcoord) uv.y = 1 - uv.y;

    return shp;
}

texture* add_texture(const std::string& path_, scene* scn) {
    auto path = path_;
    if (ygl::path_extension(path) == ".pfm")
        path = ygl::replace_path_extension(path, ".hdr");
    if (path == "") return nullptr;
    for (auto t : scn->textures)
        if (t->path == path) return t;
    auto txt = new texture();
    txt->name = path;
    txt->path = path;
    scn->textures.push_back(txt);
    return txt;
}

std::pair<vec3f, texture*> add_scaled_texture(const json& js, scene* scn) {
    if (js.is_number() || js.is_array()) {
        return {to_vec3f(js), nullptr};
    } else if (js.is_string()) {
        return {vec3f{1, 1, 1}, add_texture(js.get<std::string>(), scn)};
    } else if (js.is_object()) {
        if (js.count("type") && js.at("type").get<std::string>() != "bitmap") {
            log_error("canont handle texture type");
            return {vec3f{1, 1, 1}, nullptr};
        }
        auto kx = vec3f{1, 1, 1};
        auto txt = (texture*)nullptr;
        if (js.count("scale")) kx = to_vec3f(js.at("scale"));
        if (js.count("file"))
            txt = add_texture(js.at("file").get<std::string>(), scn);
        return {kx, txt};
    } else {
        log_error("cannot handle texture");
        return {vec3f{1, 1, 1}, nullptr};
    }
}

auto metal_ior = std::unordered_map<std::string, std::pair<vec3f, vec3f>>{
    {"Ag", {vec3f(0.1552646489f, 0.1167232965f, 0.1383806959f),
               vec3f(4.8283433224f, 3.1222459278f, 2.1469504455f)}},
    {"Al", {vec3f(1.6574599595f, 0.8803689579f, 0.5212287346f),
               vec3f(9.2238691996f, 6.2695232477f, 4.8370012281f)}},
    {"AlSb", {vec3f(-0.0485225705f, 4.1427547893f, 4.6697691348f),
                 vec3f(-0.0363741915f, 0.0937665154f, 1.3007390124f)}},
    {"Au", {vec3f(0.1431189557f, 0.3749570432f, 1.4424785571f),
               vec3f(3.9831604247f, 2.3857207478f, 1.6032152899f)}},
    {"Cr", {vec3f(4.3696828663f, 2.9167024892f, 1.6547005413f),
               vec3f(5.2064337956f, 4.2313645277f, 3.7549467933f)}},
    {"Cu", {vec3f(0.2004376970f, 0.9240334304f, 1.1022119527f),
               vec3f(3.9129485033f, 2.4528477015f, 2.1421879552f)}},
    {"Li", {vec3f(0.2657871942f, 0.1956102432f, 0.2209198538f),
               vec3f(3.5401743407f, 2.3111306542f, 1.6685930000f)}},
    {"TiN", {vec3f(1.6484691607f, 1.1504482522f, 1.3797795097f),
                vec3f(3.3684596226f, 1.9434888540f, 1.1020123347f)}},
    {"VN", {vec3f(2.8656011588f, 2.1191817791f, 1.9400767149f),
               vec3f(3.0323264950f, 2.0561075580f, 1.6162930914f)}},
    {"V", {vec3f(4.2775126218f, 3.5131538236f, 2.7611257461f),
              vec3f(3.4911844504f, 2.8893580874f, 3.1116965117f)}},
    {"W", {vec3f(4.3707029924f, 3.3002972445f, 2.9982666528f),
              vec3f(3.5006778591f, 2.6048652781f, 2.2731930614f)}}};

void parse_bsdf(material* mat, const json& js, scene* scn) {
    if (js.count("name")) mat->name = js.at("name");
    mat->double_sided = true;
    auto type = js.at("type").get<std::string>();
    auto albedo = vec3f{1, 1, 1};
    auto albedo_txt = (texture*)nullptr;
    auto bump_scale = vec3f{1, 1, 1};
    auto bump_txt = (texture*)nullptr;
    auto alpha = vec3f{1, 1, 1};
    auto alpha_txt = (texture*)nullptr;
    if (js.count("albedo")) {
        std::tie(albedo, albedo_txt) = add_scaled_texture(js.at("albedo"), scn);
    }
    if (js.count("alpha")) {
        std::tie(alpha, alpha_txt) = add_scaled_texture(js.at("alpha"), scn);
    }
    if (js.count("bump")) {
        std::tie(bump_scale, bump_txt) = add_scaled_texture(js.at("bump"), scn);
    }
    if (type == "null") {
        mat->kd = zero3f;
    } else if (type == "lambert") {
        mat->kd = albedo;
        mat->kd_txt.txt = albedo_txt;
        mat->bump_txt.txt = bump_txt;
    } else if (type == "plastic" || type == "rough_plastic") {
        mat->kd = albedo;
        mat->kd_txt.txt = albedo_txt;
        mat->ks = vec3f{0.04f, 0.04f, 0.04f};
        mat->rs =
            (js.count("roughness")) ? sqrt(js.at("roughness").get<float>()) : 0;
        mat->bump_txt.txt = bump_txt;
    } else if (type == "conductor" || type == "rough_conductor") {
        if (js.count("material")) {
            auto ctype = js.at("material").get<std::string>();
            auto eta = metal_ior.at(ctype).first;
            auto k = metal_ior.at(ctype).second;
            mat->ks = albedo * fresnel_metal(1, eta, k);
        } else {
            mat->ks = albedo;
        }
        mat->rs = (js.count("roughness")) ? sqrt(js.at("roughness").get<float>()) : 0;
        mat->bump_txt.txt = bump_txt;
    } else if (type == "mirror") {
        mat->ks = albedo;
        mat->rs = 0;
        mat->bump_txt.txt = bump_txt;
    } else if (type == "transparency") {
        parse_bsdf(mat, js.at("base"), scn);
        mat->op = alpha.x;
        if (alpha_txt) log_warning("missing alpha texture");
    } else if (type == "smooth_coat") {
        parse_bsdf(mat, js.at("substrate"), scn);
        auto sigma = to_vec3f(js.at("sigma_a")) * to_vec3f(js.at("thickness"));
        auto kx = vec3f{std::exp(-sigma.x * 2), std::exp(-sigma.y * 2),
            std::exp(-sigma.z * 2)};
#if 1
        // mat->kr = {0.04f,0.04f,0.04f};
        mat->ks *= kx;
#else
        mat->ks = {0.04f, 0.04f, 0.04f};
        mat->rs = 0;
        mat->kd = kx;
#endif
        log_warning("bad translation of smooth coat");
    } else if (type == "thinsheet") {
        mat->ks = vec3f{0.04f, 0.04f, 0.04f};
        mat->rs = (js.count("roughness")) ? js.at("roughness").get<float>() : 0;
        mat->kt = albedo;
        mat->kt_txt.txt = albedo_txt;
        mat->bump_txt.txt = bump_txt;
    } else if (type == "dielectric" || type == "rough_dielectric") {
        mat->ks = vec3f{0.04f, 0.04f, 0.04f};
        mat->rs = (js.count("roughness")) ? js.at("roughness").get<float>() : 0;
        mat->kt = albedo;
        mat->kt_txt.txt = albedo_txt;
        mat->bump_txt.txt = bump_txt;
    } else {
        println("Cannot handle {} material", type);
    }
}

material* add_material(const json& js, scene* scn, shape* shp) {
    if (js.count("bsdf") && js.at("bsdf").is_string() &&
        !js.count("emission") && !js.count("power")) {
        auto name = js.at("bsdf").get<std::string>();
        for (auto m : scn->materials)
            if (m->name == name) return m;
        throw std::runtime_error("bad bsdf name");
    }
    static auto mids = std::map<std::string, int>();
    auto mat = new material();
    scn->materials.push_back(mat);
    if (js.count("bsdf")) {
        auto jbsdf = js.at("bsdf");
        if (jbsdf.is_string()) {
            auto name = jbsdf.get<std::string>();
            for (auto m : scn->materials)
                if (m->name == name) *mat = *m;
        } else {
            parse_bsdf(mat, jbsdf, scn);
            if (mat->name == "") mat->name = "unnamed_mat";
        }
    }
    mat->name += "_" + std::to_string(mids[mat->name]++);
    if (js.count("emission")) {
        std::tie(mat->ke, mat->ke_txt.txt) =
            add_scaled_texture(js.at("emission"), scn);
    }
    if (js.count("power")) {
        mat->ke = to_vec3f(js.at("power"));
        auto area = 0.0f;
        for (auto t : shp->triangles)
            area += triangle_area(shp->pos[t.x], shp->pos[t.y], shp->pos[t.z]);
        for (auto q : shp->quads)
            area += quad_area(
                shp->pos[q.x], shp->pos[q.y], shp->pos[q.z], shp->pos[q.w]);
        mat->ke /= (area * pif);
    }
    if (!contains(scn->materials, mat)) scn->materials.push_back(mat);
    return mat;
}

scene* load_tungsten(const std::string& filename, bool facet_non_smooth) {
    auto dirname = path_dirname(filename);
    auto meshdir = dirname;
    if (meshdir.back() != '/') meshdir += '/';

    // load json
    std::ifstream jsstream(filename.c_str());
    if (!jsstream) throw std::runtime_error("could not load json " + filename);
    auto js = json();
    try {
        jsstream >> js;
    } catch (const std::exception& e) {
        throw std::runtime_error(
            std::string("could not load json with error ") + e.what());
    }

    // create scene
    auto scn = new scene();

    // convert cameras
    if (js.count("camera")) {
        auto& jcam = js.at("camera");
        auto cam = new camera();
        cam->name = "cam";
        if (jcam.count("transform")) {
            auto s = zero3f;
            std::tie(cam->frame, s) =
                xform_to_scaled_frame(jcam.at("transform"), &cam->focus);
        }
        cam->aspect = 1;
        if (jcam.count("resolution")) {
            auto& jres = jcam.at("resolution");
            if (jres.is_array()) {
                cam->aspect = (float)jres.at(0).get<float>() /
                              (float)jres.at(1).get<float>();
            }
        }
        if (jcam.count("fov")) {
            auto fov = jcam.at("fov").get<float>() * pif / 180.0f;
            // cam->yfov = jcam.at("fov").get<float>() * pif / 180.0f;
            cam->yfov = fov / cam->aspect;
        }
        scn->cameras.push_back(cam);
    }

    // convert materials
    if (js.count("bsdfs")) {
        auto& jmats = js.at("bsdfs");
        for (auto& jmat : jmats) {
            auto mat = new material();
            scn->materials.push_back(mat);
            parse_bsdf(mat, jmat, scn);
        }
    }

    // primitives
    auto shp_map = std::unordered_map<std::string, shape*>();
    auto used_meshes = std::unordered_set<std::string>();
    if (js.count("primitives")) {
        auto& jprims = js.at("primitives");
        for (auto& jprim : jprims) {
            auto type = jprim.at("type").get<std::string>();
            auto smooth =
                (jprim.count("smooth")) ? jprim.at("smooth").get<bool>() : true;
            auto cull = (jprim.count("backface_culling")) ?
                            jprim.at("backface_culling").get<bool>() :
                            true;
            auto compute_normals =
                (jprim.count("recompute_normals")) ?
                    jprim.at("recompute_normals").get<bool>() :
                    true;
            auto frame = identity_frame3f;
            auto scale = vec3f{1, 1, 1};
            if (!jprim.at("transform").empty()) {
                std::tie(frame, scale) =
                    xform_to_scaled_frame(jprim.at("transform"));
            }
            auto shp = (shape*)nullptr;
            if (type == "mesh") {
                auto path = jprim.at("file").get<std::string>();
                if (path_extension(path) == ".obj") {
                    path = replace_path_extension(path, ".obj");
                    auto oscn = load_scene(dirname + path);
                    if (!oscn || oscn->shapes.size() != 1 ||
                        oscn->shapes.size() != 1)
                        throw std::runtime_error("bad obj " + path);
                    shp = oscn->shapes.at(0);
                } else if (path_extension(path) == ".wo3") {
                    shp = load_wo3(dirname + path);
                } else
                    throw std::runtime_error(
                        "cannot handle mesh format " + path_extension(path));
                if (contains(used_meshes, path))
                    log_error("already used mesh " + path);
                used_meshes.insert(path);
            } else if (type == "quad" || type == "sphere" || type == "cube" ||
                       type == "infinite_sphere" ||
                       type == "infinite_sphere_cap") {
                cull = false;
                smooth = true;
                compute_normals = false;
                shp = make_shape(type, jprim);
            } else if (type == "curves" || type == "disk") {
                // TODO: implement these later
            } else {
                log_error("unknown shape type " + type);
                continue;
            }
            if (!shp) continue;
            shp->frame = frame;
            for (auto& p : shp->pos) p *= scale;
            shp->groups.push_back({"", add_material(jprim, scn, shp), false});
            if (cull) log_warning("culling enabled for type {}", type);
            if (!smooth) {
                shp->groups.at(0).faceted = true;
                if (facet_non_smooth && !shp->triangles.empty()) {
                    facet_triangles(shp->triangles, shp->pos, shp->norm,
                        shp->texcoord, shp->color, shp->radius);
                }
                // log_warning("faceting enabled");
            }
            if (compute_normals)
                log_warning("recomputing normals for type {}", type);
            if (type == "infinite_sphere") {
                auto env = new environment();
                env->name = shp->name;
                env->frame = shp->frame;
                env->frame = shp->frame * rotation_frame<float>({0, 1, 0}, pif);
                env->ke = shp->groups.at(0).mat->ke;
                env->ke_txt = shp->groups.at(0).mat->ke_txt;
                scn->environments.push_back(env);
                //     scn->materials.erase(std::find(scn->materials.begin(),
                //     scn->materials.end(), shp->groups.at(0).mat)); delete
                //     shp->groups.at(0).mat;
            } else if (type == "infinite_sphere_cap") {
                println("f {}", shp->frame);
                println(
                    "y {}", transform_direction(shp->frame, vec3f{0, 1, 0}));
                println(
                    "z {}", transform_direction(shp->frame, vec3f{0, 0, 1}));
                shp->groups.at(0).mat->ke *= env_distance * env_distance;
                scn->shapes.push_back(shp);
            } else {
                scn->shapes.push_back(shp);
            }
        }
    }

    // convert shapes
    return scn;
}

int main(int argc, char** argv) {
    // parse command line
    auto parser = ygl::make_parser(
        argc, argv, "ytungsten", "convert tungsten files to yocto");
    auto facet_non_smooth = ygl::parse_flag(
        parser, "--facet-non-smooth", "-f", "facet non smooth meshes", false);
    auto outfilename = ygl::parse_opt(
        parser, "--output", "-o", "output scene filename", "out.obj"s);
    auto filename =
        ygl::parse_arg(parser, "scene", "input scene filenames", ""s);
    if (ygl::should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // load image
    auto scn = load_tungsten(filename, facet_non_smooth);

    // save scene
    system(("mkdir -p " + path_dirname(outfilename)).c_str());
    auto opts = save_options();
    opts.save_textures = false;
    save_scene(outfilename, scn, opts);

    return 0;
}
