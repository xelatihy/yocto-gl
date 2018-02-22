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

frame3f xform_to_frame(const json& js, float* focus = nullptr) {
    if (js.is_array()) {
        auto m = identity_mat4f;
        for (auto i = 0; i < 16; i++) ((float*)&m)[i] = js.at(i).get<float>();
    } else if (js.count("look_at")) {
        auto from = to_vec3f(js.at("position"));
        auto to = to_vec3f(js.at("look_at"));
        auto up = to_vec3f(js.at("up"));
        if (focus) *focus = length(from - to);
        return lookat_frame(from, to, up);
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
        return translation_frame(pos) * rotation_yxz_frame(rot) *
               scaling_frame(scl);
    }
    log_error("bad xform");
    return identity_frame3f;
    //    throw std::runtime_error("xform not handled");
}

#endif

std::shared_ptr<shape> make_shape(const std::string& type) {
    static auto tcount = std::unordered_map<std::string, int>();
    auto shp = std::make_shared<shape>();
    shp->name = type + std::to_string(tcount[type]++);
    if (type == "quad") {
        shp->pos = {
            {-0.5f, 0.0f, -0.5f},
            {0.5f, 0.0f, -0.5f},
            {0.5f, 0.0f, 0.5f},
            {-0.5f, 0.0f, 0.5f},
        };
        shp->norm = {{0, 1, 0}, {0, 1, 0}, {0, 1, 0}, {0, 1, 0}};
        shp->texcoord = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
        shp->triangles = {{0, 1, 2}, {0, 2, 3}};
    } else if (type == "cube") {
        make_uvcube(shp->quads, shp->pos, shp->norm, shp->texcoord, 0);
        for (auto& p : shp->pos) p *= 0.5f;
    } else if (type == "sphere") {
        make_uvsphere(shp->quads, shp->pos, shp->norm, shp->texcoord, 4);
    } else
        log_error("unknow type " + type);
    return shp;
}

std::shared_ptr<shape> load_wo3(const std::string& filename) {
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

    auto shp = std::make_shared<shape>();
    shp->name = path_basename(filename);
    for (auto& vert : verts) {
        shp->pos.push_back(vert.pos);
        shp->norm.push_back(vert.norm);
        shp->texcoord.push_back(vert.texcoord);
    }
    for (auto& t : triangles) { shp->triangles.push_back({t.x, t.y, t.z}); }

    return shp;
}

std::shared_ptr<texture> add_texture(const std::string& path, const std::shared_ptr<scene>& scn) {
    if (path == "") return nullptr;
    for (auto t : scn->textures)
        if (t->path == path) return t;
    auto txt = std::make_shared<texture>();
    txt->name = path;
    txt->path = path;
    scn->textures.push_back(txt);
    return txt;
}

auto metal_ior = std::unordered_map<std::string, std::pair<vec3f, vec3f>>{
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

std::shared_ptr<material> add_bsdf(const json& js, const std::shared_ptr<scene>& scn) {
    static auto count = 0;
    auto mat = std::make_shared<material>();
    if (js.count("name"))
        mat->name = js.at("name");
    else
        mat->name = "unnamed_mat" + std::to_string(count++);
    mat->double_sided = true;
    auto type = js.at("type").get<std::string>();
    auto albedo = vec3f{1, 1, 1};
    auto albedo_txt = std::shared_ptr<texture>();
    auto bump_txt = std::shared_ptr<texture>();
    auto alpha = vec3f{1, 1, 1};
    auto alpha_txt = std::shared_ptr<texture>();
    if (js.count("albedo")) {
        auto& jalbedo = js.at("albedo");
        if (jalbedo.is_string()) {
            albedo_txt = add_texture(jalbedo.get<std::string>(), scn);
        } else {
            albedo = to_vec3f(jalbedo);
        }
    }
    if (js.count("alpha")) {
        auto& jalpha = js.at("alpha");
        if (jalpha.is_string()) {
            alpha_txt = add_texture(jalpha.get<std::string>(), scn);
        } else {
            alpha = to_vec3f(jalpha);
        }
    }
    if (js.count("bump")) {
        bump_txt = add_texture(js.at("bump").get<std::string>(), scn);
    }
    if (type == "null") {
        mat->kd = zero3f;
    } else if (type == "lambert") {
        mat->kd = albedo;
        mat->kd_txt = albedo_txt;
        mat->bump_txt = bump_txt;
    } else if (type == "plastic" || type == "rough_plastic") {
        mat->kd = albedo;
        mat->kd_txt = albedo_txt;
        mat->ks = vec3f{0.04f, 0.04f, 0.04f};
        mat->rs =
            (js.count("roughness")) ? sqrt(js.at("roughness").get<float>()) : 0;
        mat->bump_txt = bump_txt;
    } else if (type == "conductor" || type == "rough_conductor") {
        if (js.count("material")) {
            auto ctype = js.at("material").get<std::string>();
            auto eta = metal_ior.at(ctype).first;
            auto k = metal_ior.at(ctype).second;
            mat->ks = albedo * fresnel_metal(1, eta, k);
        } else {
            mat->ks = albedo;
        }
        mat->rs = (js.count("roughness")) ? js.at("roughness").get<float>() : 0;
        mat->bump_txt = bump_txt;
    } else if (type == "mirror") {
        mat->ks = albedo;
        mat->rs = 0;
        mat->bump_txt = bump_txt;
    } else if (type == "transparency") {
        auto base_bsdf = add_bsdf(js.at("base"), scn);
        base_bsdf->op = alpha.x;
        if (alpha_txt) log_warning("missing alpha texture");
        mat = base_bsdf;
        // base_bsdf->op_txt = alpha_txt;
    } else if (type == "thinsheet" || type == "dielectric" ||
               type == "rough_dielectric") {
        mat->ks = vec3f{0.04f, 0.04f, 0.04f};
        mat->rs = (js.count("roughness")) ? js.at("roughness").get<float>() : 0;
        mat->kt = albedo;
        mat->kt_txt = albedo_txt;
        mat->bump_txt = bump_txt;
    } else {
        println("Cannot handle {} material", type);
    }
    scn->materials.push_back(mat);
    return mat;
}

std::shared_ptr<material> add_material(const json& js, const std::shared_ptr<scene>& scn) {
    auto mat = std::shared_ptr<material>();
    auto jbsdf = js.at("bsdf");
    if (jbsdf.is_string()) {
        auto name = jbsdf.get<std::string>();
        for (auto m : scn->materials)
            if (m->name == name) mat = m;
    } else
        mat = add_bsdf(jbsdf, scn);
    if (js.count("emission")) {
        // FIXME
        mat->ke = to_vec3f(js.at("emission"));
    }
    return mat;
}

std::shared_ptr<scene> load_tungsten(
    const std::string& filename, const std::string& meshdirname) {
    auto dirname = path_dirname(filename);
    auto meshdir = (meshdirname == "") ? dirname : meshdirname;
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
    auto scn = std::make_shared<scene>();

    // convert cameras
    if (js.count("camera")) {
        auto& jcam = js.at("camera");
        auto cam = std::make_shared<camera>();
        cam->name = "cam";
        if (jcam.count("transform")) {
            cam->frame = xform_to_frame(jcam.at("transform"), &cam->focus);
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
        for (auto& jmat : jmats) add_bsdf(jmat, scn);
    }

    // primitives
    auto shp_map = std::unordered_map<std::string, shape_group*>();
    auto used_meshes = std::unordered_set<std::string>();
    if (js.count("primitives")) {
        auto& jprims = js.at("primitives");
        for (auto& jprim : jprims) {
            auto type = jprim.at("type").get<std::string>();
            auto frame = identity_frame3f;
            if (!jprim.at("transform").empty())
                frame = xform_to_frame(jprim.at("transform"));
            auto shp = std::shared_ptr<shape>{};
            if (type == "mesh") {
                auto path = jprim.at("file").get<std::string>();
                if (path_extension(path) == ".obj") {
                    path = replace_path_extension(path, ".obj");
                    auto oscn = load_scene(dirname + path);
                    if (!oscn || oscn->shapes.size() != 1 ||
                        oscn->shapes.at(0)->shapes.size() != 1)
                        throw std::runtime_error("bad obj " + path);
                    shp = oscn->shapes.at(0)->shapes.at(0);
                } else if (path_extension(path) == ".wo3") {
                    shp = load_wo3(dirname + path);
                } else
                    throw std::runtime_error(
                        "cannot handle mesh format " + path_extension(path));
                if (contains(used_meshes, path))
                    log_error("already used mesh " + path);
                used_meshes.insert(path);
            } else if (type == "quad" || type == "sphere" || type == "cube") {
                shp = make_shape(type);
            } else if (type == "curves" || type == "infinite_sphere" ||
                       type == "disk") {
                // TODO: implement these later
            } else {
                log_error("unknown shape type " + type);
            }
            if (!shp) continue;
            if (frame != identity_frame3f) {
                for (auto& p : shp->pos) p = transform_point(frame, p);
                for (auto& n : shp->norm) n = transform_direction(frame, n);
            }
            shp->mat = add_material(jprim, scn);
            auto sgr = std::make_shared<shape_group>();
            sgr->shapes.push_back(shp);
            sgr->name = sgr->shapes.at(0)->name;
            scn->shapes.push_back(sgr);
        }
    }

    // convert shapes
    return scn;
}

int main(int argc, char** argv) {
    // parse command line
    auto parser = ygl::make_parser(
        argc, argv, "ytungsten", "convert tungsten files to yocto");
    auto meshdirname =
        ygl::parse_opt(parser, "--meshdir", "-m", "mesh input directory", ""s);
    auto outfilename = ygl::parse_opt(
        parser, "--output", "-o", "output scene filename", "out.obj"s);
    auto filename =
        ygl::parse_arg(parser, "scene", "input scene filenames", ""s);
    if (ygl::should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // load image
    auto scn = load_tungsten(filename, meshdirname);

    // save scene
    system(("mkdir -p " + path_dirname(outfilename)).c_str());
    auto opts = save_options();
    opts.save_textures = false;
    save_scene(outfilename, scn, opts);

    return 0;
}
