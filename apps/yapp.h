//
// YSCENE: helper scene object to write demo code for YOCTO/GL library.
//

//
// LICENSE:
//
// Copyright (c) 2016 Fabio Pellacini
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

#ifndef _YAPP_H_
#define _YAPP_H_

// clang-format off
#ifndef __APPLE__
#include <GL/glew.h>
#else
#include <OpenGL/gl.h>
#include <OpenGL/gl3.h>
#endif
#include <GLFW/glfw3.h>
// clang-format on

#include "../yocto/yocto_cmd.h"
#include "../yocto/yocto_gltf.h"
#include "../yocto/yocto_obj.h"
#include "../yocto/yocto_shape.h"
#include "../yocto/yocto_glu.h"

namespace yapp {

//
// Using directives
//
using namespace ym;
using std::string;

//
// Scene geometry
//
struct shape {
    // whole shape data
    string name;     // shape name
    int matid = -1;  // index in the material array (-1 if not found)
    frame3f frame = identity_frame3f;  // frame

    // shape elements
    vector<int> points;       // points
    vector<vec2i> lines;      // lines
    vector<vec3i> triangles;  // triangles

    // vertex data
    vector<vec3f> pos;       // per-vertex position (3 float)
    vector<vec3f> norm;      // per-vertex normals (3 float)
    vector<vec2f> texcoord;  // per-vertex texcoord (2 float)
    vector<vec3f> color;     // [extension] per-vertex color (3 float)
    vector<float> radius;    // [extension] per-vertex radius (1 float)
};

//
// Scene Material
//
struct material {
    // whole material data
    string name;  // material name

    // color information
    vec3f ke = zero3f;  // emission color
    vec3f kd = zero3f;  // diffuse color
    vec3f ks = zero3f;  // specular color
    float rs = 0.0001;  // roughness

    // indices in the texture array (-1 if not found)
    int ke_txt = -1;
    int kd_txt = -1;
    int ks_txt = -1;
    int rs_txt = -1;
};

//
// Scene Texture
//
struct texture {
    string path;       // path
    image<vec4f> hdr;  // if loaded, hdr data
    image<vec4b> ldr;  // if loaded, ldr data
};

//
// Scene Camera
//
struct camera {
    string name;                       // name
    frame3f frame = identity_frame3f;  // frame
    bool ortho = false;                // ortho camera
    float yfov = tan(pif / 3);         // vertical field of view
    float aspect = 16.0f / 9.0f;       // aspect ratio
    float aperture = 0;                // lens aperture
    float focus = 1;                   // focus distance
};

//
// Envinonment map
//
struct environment {
    string name;     // name
    int matid = -1;  // index of material in material array (-1 if not found)
    frame3f frame = identity_frame3f;  // frame
};

//
// Asset
//
struct scene {
    vector<shape> shapes;              // shape array
    vector<material> materials;        // material array
    vector<texture> textures;          // texture array
    vector<camera> cameras;            // camera array
    vector<environment> environments;  // environment array
};

//
// backward compatible calls
//
inline int get_etype(const shape& shape) {
    if (!shape.points.empty()) {
        assert(shape.lines.empty() && shape.triangles.empty());
        return 1;
    } else if (!shape.lines.empty()) {
        assert(shape.points.empty() && shape.triangles.empty());
        return 2;
    } else if (!shape.triangles.empty()) {
        assert(shape.points.empty() && shape.lines.empty());
        return 3;
    } else
        return 0;
}

inline int get_nelems(const shape& shape) {
    auto et = get_etype(shape);
    if (et == 1) return (int)shape.points.size();
    if (et == 2) return (int)shape.lines.size();
    if (et == 3) return (int)shape.triangles.size();
    return 0;
}

inline const int* get_elems(const shape& shape) {
    auto et = get_etype(shape);
    if (et == 1) return shape.points.data();
    if (et == 2) return (int*)shape.lines.data();
    if (et == 3) return (int*)shape.triangles.data();
    return nullptr;
}

inline image<vec4b> make_image4b(int w, int h, int nc, const unsigned char* d) {
    auto img = image<vec4b>({w, h}, vec4b(0, 0, 0, 255));
    for (auto i = 0; i < w * h; i++) {
        for (auto c = 0; c < nc; c++) img.data()[i][c] = d[i * nc + c];
    }
    return img;
}

inline image<vec4f> make_image4f(int w, int h, int nc, const float* d) {
    auto img = image<vec4f>({w, h}, vec4f(0, 0, 0, 1));
    for (auto i = 0; i < w * h; i++) {
        for (auto c = 0; c < nc; c++) img.data()[i][c] = d[i * nc + c];
    }
    return img;
}

//
// Conversions
//
template <typename T, size_t N, unsigned long M>
inline mat<T, N, M> vcast(const array<array<T, N>, M>& v) {
    return *reinterpret_cast<const mat<T, N, M>*>(&v);
}
template <typename T, size_t N, size_t M>
inline array<array<T, N>, M> vcast(const mat<T, N, M>& v) {
    return *reinterpret_cast<const array<array<T, N>, M>*>(&v);
}
template <typename T, size_t N, size_t M>
inline mat<T, N, M> fvcast(const array<T, N * M>& v) {
    return *reinterpret_cast<const mat<T, N, M>*>(&v);
}
template <typename T, size_t N, size_t M>
inline array<T, N * M> fvcast(const mat<T, N, M>& v) {
    return *reinterpret_cast<const array<T, N * M>*>(&v);
}

//
// Loads a scene from obj.
//
inline bool load_obj_scene(const string& filename, scene& scene,
                           string& errmsg) {
    // clear scene
    scene = yapp::scene();

    // load scene
    yobj::obj obj;
    if (!yobj::load_obj(filename, obj, errmsg)) return false;

    // flatten to scene
    auto fl_scene = yobj::flatten_obj(obj);

    // load textures
    if (!yobj::load_textures(fl_scene, ycmd::get_dirname(filename), errmsg)) {
        printf("%s", errmsg.c_str());
        return false;
    }

    // convert cameras
    for (auto& fl_cam : fl_scene.cameras) {
        scene.cameras.emplace_back();
        auto& cam = scene.cameras.back();
        cam.name = fl_cam.name;
        cam.frame = to_frame(vcast(fl_cam.xform));
        cam.ortho = fl_cam.ortho;
        cam.yfov = fl_cam.yfov;
        cam.aspect = fl_cam.aspect;
        cam.aperture = fl_cam.aperture;
        cam.focus = fl_cam.focus;
    }

    // convert shapes
    for (auto& fl_shape : fl_scene.shapes) {
        scene.shapes.emplace_back();
        auto& sh = scene.shapes.back();
        sh.name = fl_shape.name;
        sh.frame = identity_frame3f;
        sh.matid = fl_shape.matid;
        sh.pos = vector<vec3f>(fl_shape.pos.begin(), fl_shape.pos.end());
        sh.norm = vector<vec3f>(fl_shape.norm.begin(), fl_shape.norm.end());
        sh.texcoord =
            vector<vec2f>(fl_shape.texcoord.begin(), fl_shape.texcoord.end());
        sh.color = vector<vec3f>(fl_shape.color.begin(), fl_shape.color.end());
        sh.radius = fl_shape.radius;
        sh.points = fl_shape.points;
        sh.lines = vector<vec2i>(fl_shape.lines.begin(), fl_shape.lines.end());
        sh.triangles =
            vector<vec3i>(fl_shape.triangles.begin(), fl_shape.triangles.end());
    }

    // convert materials
    for (auto& fl_mat : fl_scene.materials) {
        scene.materials.emplace_back();
        auto& mat = scene.materials.back();
        mat.name = fl_mat.name;
        mat.ke = fl_mat.ke;
        mat.kd = fl_mat.kd;
        mat.ks = fl_mat.ks;
        mat.rs = fl_mat.rs;
        mat.ke_txt = fl_mat.ke_txt;
        mat.kd_txt = fl_mat.kd_txt;
        mat.ks_txt = fl_mat.ks_txt;
        mat.rs_txt = fl_mat.rs_txt;
    }

    // convert textures
    for (auto& fl_txt : fl_scene.textures) {
        scene.textures.emplace_back();
        auto& txt = scene.textures.back();
        txt.path = fl_txt.path;
        if (!fl_txt.datab.empty())
            txt.ldr = make_image4b(fl_txt.width, fl_txt.height, fl_txt.ncomp,
                                   fl_txt.datab.data());
        if (!fl_txt.dataf.empty())
            txt.hdr = make_image4f(fl_txt.width, fl_txt.height, fl_txt.ncomp,
                                   fl_txt.dataf.data());
    }

    // convert envs
    for (auto& fl_env : fl_scene.environments) {
        scene.environments.emplace_back();
        auto& env = scene.environments.back();
        env.name = fl_env.name;
        env.frame = to_frame(vcast(fl_env.xform));
        env.matid = fl_env.matid;
    }

    // done
    return true;
}

//
// Saves a scene to obj.
//
inline bool save_obj_scene(const string& filename, const scene& scene,
                           string& errmsg) {
    // flatten to scene
    auto fl_scene = yobj::fl_scene();

    // convert cameras
    for (auto& cam : scene.cameras) {
        fl_scene.cameras.emplace_back();
        auto& fl_cam = fl_scene.cameras.back();
        fl_cam.name = cam.name;
        fl_cam.xform = vcast(to_mat(cam.frame));
        fl_cam.ortho = cam.ortho;
        fl_cam.yfov = cam.yfov;
        fl_cam.aspect = cam.aspect;
        fl_cam.aperture = cam.aperture;
        fl_cam.focus = cam.focus;
    }

    // convert shapes
    for (auto& shape : scene.shapes) {
        fl_scene.shapes.emplace_back();
        auto& fl_shape = fl_scene.shapes.back();
        fl_shape.name = shape.name;
        assert(shape.frame == identity_frame3f);
        fl_shape.matid = shape.matid;
        fl_shape.pos =
            vector<array<float, 3>>(shape.pos.begin(), shape.pos.end());
        fl_shape.norm =
            vector<array<float, 3>>(shape.norm.begin(), shape.norm.end());
        fl_shape.texcoord = vector<array<float, 2>>(shape.texcoord.begin(),
                                                    shape.texcoord.end());
        fl_shape.color =
            vector<array<float, 3>>(shape.color.begin(), shape.color.end());
        fl_shape.radius = shape.radius;
        fl_shape.points = shape.points;
        fl_shape.lines =
            vector<array<int, 2>>(shape.lines.begin(), shape.lines.end());
        fl_shape.triangles = vector<array<int, 3>>(shape.triangles.begin(),
                                                   shape.triangles.end());
    }

    // convert materials
    for (auto& mat : scene.materials) {
        fl_scene.materials.emplace_back();
        auto& fl_mat = fl_scene.materials.back();
        fl_mat.name = mat.name;
        fl_mat.ke = mat.ke;
        fl_mat.kd = mat.kd;
        fl_mat.ks = mat.ks;
        fl_mat.rs = mat.rs;
        fl_mat.ke_txt = mat.ke_txt;
        fl_mat.kd_txt = mat.kd_txt;
        fl_mat.ks_txt = mat.ks_txt;
        fl_mat.rs_txt = mat.rs_txt;
    }

    // convert textures
    for (auto& txt : scene.textures) {
        fl_scene.textures.emplace_back();
        auto& fl_txt = fl_scene.textures.back();
        fl_txt.path = txt.path;
    }

    // convert envs
    for (auto& env : scene.environments) {
        fl_scene.environments.emplace_back();
        auto& fl_env = fl_scene.environments.back();
        fl_env.name = env.name;
        fl_env.xform = vcast(to_mat(env.frame));
        fl_env.matid = env.matid;
    }

    // save obj
    auto obj = yobj::unflatten_obj(fl_scene);
    if (!yobj::save_obj(filename, obj, errmsg)) return false;

    // done
    return true;
}

//
// Saves a scene to gltf.
//
inline bool save_gltf_scene(const string& filename, const scene& scene,
                            string& errmsg) {
    // flatten to scene
    auto fl_scene = ygltf::fl_gltf();

    // convert cameras
    for (auto& cam : scene.cameras) {
        fl_scene.cameras.emplace_back();
        auto& fl_cam = fl_scene.cameras.back();
        fl_cam.name = cam.name;
        fl_cam.xform = fvcast(to_mat(cam.frame));
        fl_cam.ortho = cam.ortho;
        fl_cam.yfov = cam.yfov;
        fl_cam.aspect = cam.aspect;
    }

    // convert shapes
    for (auto& shape : scene.shapes) {
        fl_scene.meshes.emplace_back();
        auto& fl_mesh = fl_scene.meshes.back();
        fl_mesh.name = shape.name;
        fl_mesh.xform = fvcast(to_mat(shape.frame));
        fl_mesh.primitives.push_back((int)fl_scene.primitives.size());
        fl_scene.primitives.emplace_back();
        auto& fl_shape = fl_scene.primitives.back();
        fl_shape.material = shape.matid;
        fl_shape.pos =
            vector<array<float, 3>>(shape.pos.begin(), shape.pos.end());
        fl_shape.norm =
            vector<array<float, 3>>(shape.norm.begin(), shape.norm.end());
        fl_shape.texcoord = vector<array<float, 2>>(shape.texcoord.begin(),
                                                    shape.texcoord.end());
        fl_shape.color =
            vector<array<float, 3>>(shape.color.begin(), shape.color.end());
        fl_shape.radius = shape.radius;
        fl_shape.points = shape.points;
        fl_shape.lines =
            vector<array<int, 2>>(shape.lines.begin(), shape.lines.end());
        fl_shape.triangles = vector<array<int, 3>>(shape.triangles.begin(),
                                                   shape.triangles.end());
    }

    // convert materials
    for (auto& mat : scene.materials) {
        fl_scene.materials.emplace_back();
        auto& fl_mat = fl_scene.materials.back();
        fl_mat.name = mat.name;
        fl_mat.ke = mat.ke;
        fl_mat.kd = mat.kd;
        fl_mat.ks = mat.ks;
        fl_mat.rs = mat.rs;
        fl_mat.ke_txt = mat.ke_txt;
        fl_mat.kd_txt = mat.kd_txt;
        fl_mat.ks_txt = mat.ks_txt;
        fl_mat.rs_txt = mat.rs_txt;
    }

    // convert textures
    for (auto& txt : scene.textures) {
        fl_scene.textures.emplace_back();
        auto& fl_txt = fl_scene.textures.back();
        fl_txt.path = txt.path;
    }

    // save gltf
    auto gltf =
        ygltf::unflatten_gltf(fl_scene, ycmd::get_filename(filename) + ".bin");
    if (!ygltf::save_gltf(filename, gltf, errmsg, true, false, false))
        return false;

    // done
    return true;
}

//
// Load gltf scene
//
inline bool load_gltf_scene(const string& filename, scene& scene, bool binary,
                            string& errmsg) {
    // clear scene
    scene = yapp::scene();

    // load scene
    ygltf::glTF_t gltf;
    if (binary) {
        if (!ygltf::load_binary_gltf(filename, gltf, errmsg, true, false))
            return false;
    } else {
        if (!ygltf::load_gltf(filename, gltf, errmsg, true, false))
            return false;
    }

    // flatten to scene
    auto fl_scene = ygltf::flatten_gltf(gltf, gltf.scene);

    // convert cameras
    for (auto& fl_cam : fl_scene.cameras) {
        scene.cameras.emplace_back();
        auto& cam = scene.cameras.back();
        cam.name = fl_cam.name;
        cam.frame = to_frame(fvcast<float, 4, 4>(fl_cam.xform));
        cam.ortho = fl_cam.ortho;
        cam.yfov = fl_cam.yfov;
        cam.aspect = fl_cam.aspect;
        cam.aperture = 0;
        cam.focus = 1;
    }

    // convert shapes
    for (auto& fl_mesh : fl_scene.meshes) {
        for (auto& fl_shape_id : fl_mesh.primitives) {
            auto& fl_shape = fl_scene.primitives.at(fl_shape_id);
            scene.shapes.emplace_back();
            auto& sh = scene.shapes.back();
            sh.name = fl_mesh.name;
            sh.frame = to_frame(fvcast<float, 4, 4>(fl_mesh.xform));
            sh.matid = fl_shape.material;
            sh.pos = vector<vec3f>(fl_shape.pos.begin(), fl_shape.pos.end());
            sh.norm = vector<vec3f>(fl_shape.norm.begin(), fl_shape.norm.end());
            sh.texcoord = vector<vec2f>(fl_shape.texcoord.begin(),
                                        fl_shape.texcoord.end());
            sh.color =
                vector<vec3f>(fl_shape.color.begin(), fl_shape.color.end());
            sh.radius = fl_shape.radius;
            sh.points = fl_shape.points;
            sh.lines =
                vector<vec2i>(fl_shape.lines.begin(), fl_shape.lines.end());
            sh.triangles = vector<vec3i>(fl_shape.triangles.begin(),
                                         fl_shape.triangles.end());
        }
    }

    // convert materials
    for (auto& fl_mat : fl_scene.materials) {
        scene.materials.emplace_back();
        auto& mat = scene.materials.back();
        mat.name = fl_mat.name;
        mat.ke = fl_mat.ke;
        mat.kd = fl_mat.kd;
        mat.ks = fl_mat.ks;
        mat.rs = fl_mat.rs;
        mat.ke_txt = fl_mat.ke_txt;
        mat.kd_txt = fl_mat.kd_txt;
        mat.ks_txt = fl_mat.ks_txt;
        mat.rs_txt = fl_mat.rs_txt;
    }

    // convert textures
    for (auto& fl_txt : fl_scene.textures) {
        scene.textures.emplace_back();
        auto& txt = scene.textures.back();
        txt.path = fl_txt.path;
        if (!fl_txt.datab.empty())
            txt.ldr = make_image4b(fl_txt.width, fl_txt.height, fl_txt.ncomp,
                                   fl_txt.datab.data());
        if (!fl_txt.dataf.empty())
            txt.hdr = make_image4f(fl_txt.width, fl_txt.height, fl_txt.ncomp,
                                   fl_txt.dataf.data());
    }

    // done
    return true;
}

//
// Load scene
//
inline bool load_scene(const string& filename, scene& scene, string& errmsg) {
    // clear scene
    scene = yapp::scene();

    // get extension
    auto ext = ycmd::get_extension(filename);
    if (ext == ".obj") {
        // load obj
        if (!load_obj_scene(filename, scene, errmsg)) return false;
    } else if (ext == ".gltf") {
        // load obj
        if (!load_gltf_scene(filename, scene, false, errmsg)) return false;
    } else if (ext == ".glb") {
        // load obj
        if (!load_gltf_scene(filename, scene, true, errmsg)) return false;
    } else {
        errmsg = "unknown file type";
        return false;
    }

    // check textures and patch them up if needed
    for (auto& txt : scene.textures) {
        if (txt.hdr.empty() && txt.ldr.empty()) {
            printf("unable to load texture %s\n", txt.path.c_str());
            txt.ldr = image<vec4b>({1, 1}, vec4b(255, 255, 255, 255));
        }
    }

    // ensure normals
    for (auto& shape : scene.shapes) {
        if (!shape.norm.empty()) continue;
        shape.norm.resize(shape.pos.size());
        yshape::compute_normals(shape.points, shape.lines, shape.triangles,
                                shape.pos, shape.norm);
    }

    // ensure radius is necessary
    for (auto& shape : scene.shapes) {
        if (shape.points.empty() && shape.lines.empty()) continue;
        if (!shape.radius.empty()) continue;
        shape.radius.resize(shape.pos.size(), 0.001f);
    }

    // make camera if not there
    if (!scene.cameras.size()) {
        // find scene bounds
        auto bbox = invalid_bbox3f;
        for (auto& shape : scene.shapes) {
            for (auto& p : shape.pos) bbox += transform_point(shape.frame, p);
        }
        auto bbox_center = center(bbox);
        auto bbox_size = diagonal(bbox);
        auto bbox_msize = fmax(bbox_size[0], fmax(bbox_size[1], bbox_size[2]));
        // create camera
        auto cam = yapp::camera();
        // set up camera
        auto camera_dir = vec3f{1, 0.4f, 1};
        auto from = camera_dir * bbox_msize + bbox_center;
        auto to = bbox_center;
        auto up = vec3f{0, 1, 0};
        cam.frame = lookat_frame3(from, to, up);
        cam.ortho = false;
        cam.aspect = 16.0f / 9.0f;
        cam.yfov = 2 * atan(0.5f);
        cam.aperture = 0;
        cam.focus = length(to - from);
        scene.cameras.push_back(cam);
    } else {
        auto bbox = invalid_bbox3f;
        for (auto& shape : scene.shapes) {
            for (auto& p : shape.pos) bbox += transform_point(shape.frame, p);
        }
        for (auto& cam : scene.cameras) {
            if (!cam.focus) cam.focus = length(cam.frame.o() - bbox.center());
        }
    }

    return true;
}

//
// Load scene
//
inline scene load_scene(const string& filename) {
    string errmsg;
    auto sc = scene();
    if (!load_scene(filename, sc, errmsg)) {
        printf("error loading scene: %s\n", errmsg.c_str());
        return {};
    }
    return sc;
}

//
// Init shading
//
inline void init_shade(const yapp::scene& scene, yglu::uint& shade_prog,
                       yglu::uint& shade_vao, vector<yglu::uint>& shade_txt,
                       vector<array<yglu::uint, 7>>& shade_vbo) {
    yglu::stdshader::make_program(&shade_prog, &shade_vao);
    for (auto& txt : scene.textures) {
        if (!txt.hdr.empty()) {
            shade_txt.push_back(yglu::modern::make_texture(
                txt.hdr.size()[0], txt.hdr.size()[1], 4, (float*)txt.hdr.data(),
                true, true, true));
        } else if (!txt.ldr.empty()) {
            shade_txt.push_back(yglu::modern::make_texture(
                txt.ldr.size()[0], txt.ldr.size()[1], 4,
                (unsigned char*)txt.ldr.data(), true, true, true));
        } else
            assert(false);
    }
    for (auto& shape : scene.shapes) {
        shade_vbo.push_back({0, 0, 0, 0});
        auto& vbos = shade_vbo.back();
        if (!shape.pos.empty())
            vbos[0] = yglu::modern::make_buffer((int)shape.pos.size(),
                                                3 * sizeof(float),
                                                shape.pos.data(), false, false);
        if (!shape.norm.empty())
            vbos[1] = yglu::modern::make_buffer(
                (int)shape.norm.size(), 3 * sizeof(float), shape.norm.data(),
                false, false);
        if (!shape.texcoord.empty())
            vbos[2] = yglu::modern::make_buffer(
                (int)shape.texcoord.size(), 2 * sizeof(float),
                shape.texcoord.data(), false, false);
        if (!shape.color.empty())
            vbos[3] = yglu::modern::make_buffer(
                (int)shape.color.size(), 3 * sizeof(float), shape.color.data(),
                false, false);
        if (!shape.points.empty())
            vbos[4] =
                yglu::modern::make_buffer((int)shape.points.size(), sizeof(int),
                                          shape.points.data(), true, false);
        if (!shape.lines.empty())
            vbos[5] = yglu::modern::make_buffer(
                (int)shape.lines.size(), 2 * sizeof(int), shape.lines.data(),
                true, false);
        if (!shape.triangles.empty())
            vbos[6] = yglu::modern::make_buffer(
                (int)shape.triangles.size(), 3 * sizeof(int),
                shape.triangles.data(), true, false);
    }
}

//
// Display a scene
//
inline void shade(const yapp::scene& scene, int cur_camera, yglu::uint prog,
                  yglu::uint vao, const vector<yglu::uint>& txt,
                  vector<array<yglu::uint, 7>>& vbo, const vec4f& background,
                  float exposure, float gamma_, bool wireframe, bool edges,
                  bool camera_lights, const vec3f& amb) {
    // begin frame
    glEnable(GL_DEPTH_TEST);
    glClearDepth(1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glDisable(GL_CULL_FACE);

    if (wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    auto& cam = scene.cameras[cur_camera];
    auto camera_xform = to_mat(cam.frame);
    auto camera_view = to_mat(inverse(cam.frame));
    auto camera_proj = perspective_mat4(cam.yfov, cam.aspect, 0.1f, 10000.0f);

    yglu::stdshader::begin_frame(prog, vao, camera_lights, exposure, gamma_,
                                 vcast(camera_xform), vcast(camera_view),
                                 vcast(camera_proj));

    if (!camera_lights) {
        auto nlights = 0;
        array<vec3f, 16> light_pos, light_ke;
        array<yglu::ltype, 16> light_type;
        for (auto&& shape : scene.shapes) {
            if (shape.matid < 0) continue;
            auto&& mat = scene.materials[shape.matid];
            if (mat.ke == zero3f) continue;
            for (auto p : shape.points) {
                if (nlights >= 16) continue;
                light_pos[nlights] = shape.pos[p];
                light_pos[nlights] =
                    transform_point(shape.frame, light_pos[nlights]);
                light_ke[nlights] = mat.ke;
                light_type[nlights] = yglu::ltype::point;
                nlights++;
            }
        }
        yglu::stdshader::set_lights(
            prog, amb, nlights, (yglu::float3*)light_pos.data(),
            (yglu::float3*)light_ke.data(), light_type.data());
    }

    for (auto sid = 0; sid < scene.shapes.size(); sid++) {
        auto&& shape = scene.shapes[sid];
        yglu::stdshader::begin_shape(prog, vcast(to_mat(shape.frame)));

        if (shape.matid >= 0) {
            auto&& mat = scene.materials[shape.matid];

#define __txt(i) ((i >= 0) ? txt[i] : 0)
            yglu::stdshader::set_material(
                prog, mat.ke, mat.kd, mat.ks, mat.rs, __txt(mat.ke_txt),
                __txt(mat.kd_txt), __txt(mat.ks_txt), __txt(mat.rs_txt), false);
        } else {
            auto kd = vec3f{0.8f, 0.8f, 0.8f};
            yglu::stdshader::set_material(prog, {0, 0, 0}, kd, {0, 0, 0}, 0, 0,
                                          0, 0, 0, false);
        }

        auto bids = vbo[sid];
        yglu::stdshader::set_vert(prog, bids[0], bids[1], bids[2], bids[3]);

        yglu::stdshader::draw_points(prog, (int)shape.points.size(), bids[4]);
        yglu::stdshader::draw_lines(prog, (int)shape.lines.size(), bids[5]);
        yglu::stdshader::draw_triangles(prog, (int)shape.triangles.size(),
                                        bids[6]);

        if (edges && !wireframe) {
            yglu::stdshader::set_material(prog, {0, 0, 0}, {0, 0, 0}, {0, 0, 0},
                                          0, 0, 0, 0, 0, false);

            glLineWidth(2);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glDepthRange(0, 0.999999);
            yglu::stdshader::draw_triangles(prog, (int)shape.triangles.size(),
                                            bids[6]);
            glDepthRange(0, 1);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glLineWidth(1);
        }

        yglu::stdshader::end_shape();
    }

    yglu::stdshader::end_frame();

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

//
// Init shading
//
inline void init_draw(const yapp::scene& scene, vector<yglu::uint>& shade_txt) {
    for (auto& txt : scene.textures) {
        if (!txt.hdr.empty()) {
            shade_txt.push_back(yglu::legacy::make_texture(
                txt.hdr.size()[0], txt.hdr.size()[1], 4, (float*)txt.hdr.data(),
                true, true));
        } else if (!txt.ldr.empty()) {
            shade_txt.push_back(yglu::legacy::make_texture(
                txt.ldr.size()[0], txt.ldr.size()[1], 4,
                (unsigned char*)txt.ldr.data(), true, true));
        } else
            assert(false);
    }
}

//
// Draw a scene
//
inline void draw(const yapp::scene& scene, int cur_camera,
                 const vector<yglu::uint>& txt, const vec4f& background,
                 float exposure, float gamma_, bool wireframe, bool edges,
                 bool camera_lights, const vec3f& amb) {
    // begin frame
    glEnable(GL_DEPTH_TEST);
    glClearDepth(1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glDisable(GL_CULL_FACE);

    if (wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    auto& cam = scene.cameras[cur_camera];
    auto camera_xform = to_mat(cam.frame);
    auto camera_view = to_mat(inverse(cam.frame));
    auto camera_proj = perspective_mat4(cam.yfov, cam.aspect, 0.1f, 10000.0f);

    auto nlights = 0;
    array<vec3f, 16> light_pos, light_ke;
    array<yglu::ltype, 16> light_type;

    yglu::legacy::begin_frame(vcast(camera_xform), vcast(camera_view),
                              vcast(camera_proj), camera_lights, true);

    if (!camera_lights) {
        for (auto&& shape : scene.shapes) {
            if (shape.matid < 0) continue;
            auto&& mat = scene.materials[shape.matid];
            if (mat.ke == zero3f) continue;
            for (auto p : shape.points) {
                if (nlights >= 16) continue;
                light_pos[nlights] = shape.pos[p];
                light_pos[nlights] =
                    transform_point(shape.frame, light_pos[nlights]);
                light_ke[nlights] = mat.ke;
                light_type[nlights] = yglu::ltype::point;
                nlights++;
            }
        }
        yglu::legacy::set_lights(amb, nlights, (yglu::float3*)light_pos.data(),
                                 (yglu::float3*)light_ke.data(),
                                 light_type.data());
    }

    for (auto&& shape : scene.shapes) {
        yglu::legacy::begin_shape(vcast(to_mat(shape.frame)));

        if (shape.matid >= 0) {
            auto&& mat = scene.materials[shape.matid];

#define __txt(i) ((i >= 0) ? txt[i] : 0)
            yglu::legacy::set_material(
                mat.ke, mat.kd, mat.ks,
                yglu::legacy::specular_roughness_to_exponent(mat.rs),
                __txt(mat.kd_txt), true);
        } else {
            auto kd = vec3f{0.8f, 0.8f, 0.8f};
            yglu::legacy::set_material({0, 0, 0}, kd, {0, 0, 0}, 0, 0, true);
        }

        yglu::legacy::draw_points((int)shape.points.size(), shape.points.data(),
                                  (yglu::float3*)shape.pos.data(),
                                  (yglu::float3*)shape.norm.data(),
                                  (yglu::float2*)shape.texcoord.data(),
                                  (yglu::float3*)shape.color.data());
        yglu::legacy::draw_lines(
            (int)shape.lines.size(), (yglu::int2*)shape.lines.data(),
            (yglu::float3*)shape.pos.data(), (yglu::float3*)shape.norm.data(),
            (yglu::float2*)shape.texcoord.data(),
            (yglu::float3*)shape.color.data());
        yglu::legacy::draw_triangles(
            (int)shape.triangles.size(), (yglu::int3*)shape.triangles.data(),
            (yglu::float3*)shape.pos.data(), (yglu::float3*)shape.norm.data(),
            (yglu::float2*)shape.texcoord.data(),
            (yglu::float3*)shape.color.data());

        if (edges && !wireframe) {
            yglu::legacy::set_material({0, 0, 0}, {0, 0, 0}, {0, 0, 0}, 0, 0,
                                       true);

            glLineWidth(2);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glDepthRange(0, 0.999999);
            yglu::legacy::draw_triangles((int)shape.triangles.size(),
                                         (yglu::int3*)shape.triangles.data(),
                                         (yglu::float3*)shape.pos.data(),
                                         (yglu::float3*)shape.norm.data(),
                                         (yglu::float2*)shape.texcoord.data(),
                                         (yglu::float3*)shape.color.data());
            glDepthRange(0, 1);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glLineWidth(1);
        }

        yglu::legacy::end_shape();
    }

    yglu::legacy::end_frame();

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

}  // namespace

#endif
