//
// YAPP: helper scene object to write demo code for YOCTO/GL library.
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

#include "../yocto/yocto_cmd.h"
#include "../yocto/yocto_gltf.h"
#include "../yocto/yocto_obj.h"
#include "../yocto/yocto_shape.h"

namespace yapp {

//
// Typedefs for vec/mat types
//
using float2 = std::array<float, 2>;
using float3 = std::array<float, 3>;
using float3x4 = std::array<std::array<float, 3>, 4>;
using float4x4 = std::array<std::array<float, 4>, 4>;
using int2 = std::array<int, 2>;
using int3 = std::array<int, 3>;

//
// Scene geometry
//
struct shape {
    // whole shape data
    std::string name;  // shape name
    int matid = -1;    // index in the material array (-1 if not found)
    float3x4 frame = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}}};  // frame

    // shape elements
    std::vector<int> points;      // points
    std::vector<int2> lines;      // lines
    std::vector<int3> triangles;  // triangles

    // vertex data
    std::vector<float3> pos;       // per-vertex position (3 float)
    std::vector<float3> norm;      // per-vertex normals (3 float)
    std::vector<float2> texcoord;  // per-vertex texcoord (2 float)
    std::vector<float3> color;     // [extension] per-vertex color (3 float)
    std::vector<float> radius;     // [extension] per-vertex radius (1 float)
};

//
// Scene Material
//
struct material {
    // whole material data
    std::string name;  // material name

    // color information
    float3 ke = {0, 0, 0};  // emission color
    float3 kd = {0, 0, 0};  // diffuse color
    float3 ks = {0, 0, 0};  // specular color
    float rs = 0.0001;      // roughness

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
    std::string path;                // path
    int width = 0;                   // width
    int height = 0;                  // height
    int ncomp = 0;                   // number of components
    std::vector<float> hdr;          // if loaded, hdr data
    std::vector<unsigned char> ldr;  // if loaded, ldr data
};

//
// Scene Camera
//
struct camera {
    std::string name;                                                 // name
    float3x4 frame = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}}};  // frame
    bool ortho = false;                  // ortho camera
    float yfov = std::tan(3.1416f / 3);  // vertical field of view
    float aspect = 16.0f / 9.0f;         // aspect ratio
    float aperture = 0;                  // lens aperture
    float focus = 1;                     // focus distance
};

//
// Envinonment map
//
struct environment {
    std::string name;  // name
    int matid = -1;    // index of material in material array (-1 if not found)
    float3x4 frame = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}}};  // frame
};

//
// Asset
//
struct scene {
    std::vector<shape*> shapes;              // shape array
    std::vector<material*> materials;        // material array
    std::vector<texture*> textures;          // texture array
    std::vector<camera*> cameras;            // camera array
    std::vector<environment*> environments;  // environment array

    ~scene();
};

//
// Destructor
//
inline scene::~scene() {
    for (auto v : shapes) delete v;
    for (auto v : materials) delete v;
    for (auto v : textures) delete v;
    for (auto v : cameras) delete v;
    for (auto v : environments) delete v;
}

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

//
// Loads a scene from obj.
//
inline scene* load_obj_scene(const std::string& filename) {
    // load scene
    auto obj = std::unique_ptr<yobj::obj>(yobj::load_obj(filename));

    // flatten to scene
    auto fl_scene =
        std::unique_ptr<yobj::fl_scene>(yobj::flatten_obj(obj.get()));

    // cleanup
    obj.reset(nullptr);

    // load textures
    yobj::load_textures(fl_scene.get(), ycmd::get_dirname(filename));

    // init scene
    auto sc = std::unique_ptr<scene>(new scene());

    // convert cameras
    for (auto fl_cam : fl_scene->cameras) {
        auto cam = new camera();
        cam->name = fl_cam->name;
        cam->frame = ym::to_frame(ym::mat4f(fl_cam->xform));
        cam->ortho = fl_cam->ortho;
        cam->yfov = fl_cam->yfov;
        cam->aspect = fl_cam->aspect;
        cam->aperture = fl_cam->aperture;
        cam->focus = fl_cam->focus;
        sc->cameras.push_back(cam);
    }

    // convert shapes
    for (auto fl_shape : fl_scene->shapes) {
        auto sh = new shape();
        sh->name = fl_shape->name;
        sh->frame = ym::identity_frame3f;
        sh->matid = fl_shape->matid;
        sh->pos = fl_shape->pos;
        sh->norm = fl_shape->norm;
        sh->texcoord = fl_shape->texcoord;
        sh->color = fl_shape->color;
        sh->radius = fl_shape->radius;
        sh->points = fl_shape->points;
        sh->lines = fl_shape->lines;
        sh->triangles = fl_shape->triangles;
        sc->shapes.push_back(sh);
    }

    // convert materials
    for (auto fl_mat : fl_scene->materials) {
        auto mat = new material();
        mat->name = fl_mat->name;
        mat->ke = fl_mat->ke;
        mat->kd = fl_mat->kd;
        mat->ks = fl_mat->ks;
        mat->rs = fl_mat->rs;
        mat->ke_txt = fl_mat->ke_txt;
        mat->kd_txt = fl_mat->kd_txt;
        mat->ks_txt = fl_mat->ks_txt;
        mat->rs_txt = fl_mat->rs_txt;
        sc->materials.push_back(mat);
    }

    // convert textures
    for (auto fl_txt : fl_scene->textures) {
        auto txt = new texture();
        txt->path = fl_txt->path;
        txt->width = fl_txt->width;
        txt->height = fl_txt->height;
        txt->ncomp = fl_txt->ncomp;
        txt->ldr = fl_txt->datab;
        txt->hdr = fl_txt->dataf;
        sc->textures.push_back(txt);
    }

    // convert envs
    for (auto fl_env : fl_scene->environments) {
        auto env = new environment();
        env->name = fl_env->name;
        env->frame = ym::to_frame(ym::mat4f(fl_env->xform));
        env->matid = fl_env->matid;
        sc->environments.push_back(env);
    }

    // done
    return sc.release();
}

//
// Saves a scene to obj.
//
inline void save_obj_scene(const std::string& filename, const scene* sc) {
    // flatten to scene
    auto fl_scene = std::unique_ptr<yobj::fl_scene>(new yobj::fl_scene());

    // convert cameras
    for (auto cam : sc->cameras) {
        auto fl_cam = new yobj::fl_camera();
        fl_cam->name = cam->name;
        fl_cam->xform = ym::to_mat((ym::frame3f)cam->frame);
        fl_cam->ortho = cam->ortho;
        fl_cam->yfov = cam->yfov;
        fl_cam->aspect = cam->aspect;
        fl_cam->aperture = cam->aperture;
        fl_cam->focus = cam->focus;
        fl_scene->cameras.push_back(fl_cam);
    }

    // convert shapes
    for (auto shape : sc->shapes) {
        auto fl_shape = new yobj::fl_shape();
        fl_shape->name = shape->name;
        fl_shape->matid = shape->matid;
        fl_shape->pos.resize(shape->pos.size());
        for (auto i = 0; i < shape->pos.size(); i++) {
            fl_shape->pos[i] = ym::transform_point((ym::frame3f)shape->frame,
                                                   (ym::vec3f)shape->pos[i]);
        }
        fl_shape->norm.resize(shape->norm.size());
        for (auto i = 0; i < shape->pos.size(); i++) {
            fl_shape->norm[i] = ym::transform_direction(
                (ym::frame3f)shape->frame, (ym::vec3f)shape->norm[i]);
        }
        fl_shape->texcoord = std::vector<std::array<float, 2>>(
            shape->texcoord.begin(), shape->texcoord.end());
        fl_shape->color = std::vector<std::array<float, 3>>(
            shape->color.begin(), shape->color.end());
        fl_shape->radius = shape->radius;
        fl_shape->points = shape->points;
        fl_shape->lines = std::vector<std::array<int, 2>>(shape->lines.begin(),
                                                          shape->lines.end());
        fl_shape->triangles = std::vector<std::array<int, 3>>(
            shape->triangles.begin(), shape->triangles.end());
        fl_scene->shapes.push_back(fl_shape);
    }

    // convert materials
    for (auto mat : sc->materials) {
        auto fl_mat = new yobj::fl_material();
        fl_mat->name = mat->name;
        fl_mat->ke = mat->ke;
        fl_mat->kd = mat->kd;
        fl_mat->ks = mat->ks;
        fl_mat->rs = mat->rs;
        fl_mat->ke_txt = mat->ke_txt;
        fl_mat->kd_txt = mat->kd_txt;
        fl_mat->ks_txt = mat->ks_txt;
        fl_mat->rs_txt = mat->rs_txt;
        fl_scene->materials.push_back(fl_mat);
    }

    // convert textures
    for (auto txt : sc->textures) {
        auto fl_txt = new yobj::fl_texture();
        fl_txt->path = txt->path;
        fl_scene->textures.push_back(fl_txt);
    }

    // convert envs
    for (auto env : sc->environments) {
        auto fl_env = new yobj::fl_environment();
        fl_env->name = env->name;
        fl_env->xform = ym::to_mat((ym::frame3f)env->frame);
        fl_env->matid = env->matid;
        fl_scene->environments.push_back(fl_env);
    }

    // save obj
    auto obj = std::unique_ptr<yobj::obj>(yobj::unflatten_obj(fl_scene.get()));
    yobj::save_obj(filename, obj.get());
}

//
// Saves a scene to gltf.
//
inline void save_gltf_scene(const std::string& filename, const scene* sc) {
    // flatten to scene
    auto fl_scene = std::unique_ptr<ygltf::fl_gltf>(new ygltf::fl_gltf());

    // convert cameras
    for (auto cam : sc->cameras) {
        auto fl_cam = new ygltf::fl_camera();
        fl_cam->name = cam->name;
        fl_cam->xform = ym::to_mat((ym::frame3f)cam->frame);
        fl_cam->ortho = cam->ortho;
        fl_cam->yfov = cam->yfov;
        fl_cam->aspect = cam->aspect;
        fl_scene->cameras.push_back(fl_cam);
    }

    // convert shapes
    for (auto shape : sc->shapes) {
        auto fl_mesh = new ygltf::fl_mesh();
        fl_mesh->name = shape->name;
        fl_mesh->xform = ym::to_mat((ym::frame3f)shape->frame);
        fl_mesh->primitives.push_back((int)fl_scene->primitives.size());
        auto fl_shape = new ygltf::fl_primitives();
        fl_shape->material = shape->matid;
        fl_shape->pos = std::vector<std::array<float, 3>>(shape->pos.begin(),
                                                          shape->pos.end());
        fl_shape->norm = std::vector<std::array<float, 3>>(shape->norm.begin(),
                                                           shape->norm.end());
        fl_shape->texcoord = std::vector<std::array<float, 2>>(
            shape->texcoord.begin(), shape->texcoord.end());
        fl_shape->color = std::vector<std::array<float, 3>>(
            shape->color.begin(), shape->color.end());
        fl_shape->radius = shape->radius;
        fl_shape->points = shape->points;
        fl_shape->lines = std::vector<std::array<int, 2>>(shape->lines.begin(),
                                                          shape->lines.end());
        fl_shape->triangles = std::vector<std::array<int, 3>>(
            shape->triangles.begin(), shape->triangles.end());
        fl_scene->primitives.push_back(fl_shape);
        fl_scene->meshes.push_back(fl_mesh);
    }

    // convert materials
    for (auto mat : sc->materials) {
        auto fl_mat = new ygltf::fl_material();
        fl_mat->name = mat->name;
        fl_mat->ke = mat->ke;
        fl_mat->kd = mat->kd;
        fl_mat->ks = mat->ks;
        fl_mat->rs = mat->rs;
        fl_mat->ke_txt = mat->ke_txt;
        fl_mat->kd_txt = mat->kd_txt;
        fl_mat->ks_txt = mat->ks_txt;
        fl_mat->rs_txt = mat->rs_txt;
        fl_scene->materials.push_back(fl_mat);
    }

    // convert textures
    for (auto txt : sc->textures) {
        auto fl_txt = new ygltf::fl_texture();
        fl_txt->path = txt->path;
        fl_scene->textures.push_back(fl_txt);
    }

    // save gltf
    auto gltf = std::unique_ptr<ygltf::glTF_t>(ygltf::unflatten_gltf(
        fl_scene.get(), ycmd::get_filename(filename) + ".bin"));
    ygltf::save_gltf(filename, gltf.get(), true, false, false);
}

//
// Load gltf scene
//
inline scene* load_gltf_scene(const std::string& filename, bool binary) {
    // load scene
    auto gltf = std::unique_ptr<ygltf::glTF_t>(
        (binary) ? ygltf::load_binary_gltf(filename, true, false)
                 : ygltf::load_gltf(filename, true, false));

    // flatten to scene
    auto fl_scene = std::unique_ptr<ygltf::fl_gltf>(
        ygltf::flatten_gltf(gltf.get(), gltf->scene));

    // init scene
    auto sc = std::unique_ptr<scene>(new scene());

    // convert cameras
    for (auto fl_cam : fl_scene->cameras) {
        auto cam = new camera();
        cam->name = fl_cam->name;
        cam->frame = ym::to_frame(ym::mat4f(fl_cam->xform));
        cam->ortho = fl_cam->ortho;
        cam->yfov = fl_cam->yfov;
        cam->aspect = fl_cam->aspect;
        cam->aperture = 0;
        cam->focus = 1;
        sc->cameras.push_back(cam);
    }

    // convert shapes
    for (auto fl_mesh : fl_scene->meshes) {
        for (auto& fl_shape_id : fl_mesh->primitives) {
            auto fl_shape = fl_scene->primitives.at(fl_shape_id);
            auto sh = new shape();
            sh->name = fl_mesh->name;
            sh->frame = ym::to_frame(ym::mat4f(fl_mesh->xform));
            sh->matid = fl_shape->material;
            sh->pos = fl_shape->pos;
            sh->norm = fl_shape->norm;
            sh->texcoord = fl_shape->texcoord;
            sh->color = fl_shape->color;
            sh->radius = fl_shape->radius;
            sh->points = fl_shape->points;
            sh->lines = fl_shape->lines;
            sh->triangles = fl_shape->triangles;
            sc->shapes.push_back(sh);
        }
    }

    // convert materials
    for (auto fl_mat : fl_scene->materials) {
        auto mat = new material();
        mat->name = fl_mat->name;
        mat->ke = fl_mat->ke;
        mat->kd = fl_mat->kd;
        mat->ks = fl_mat->ks;
        mat->rs = fl_mat->rs;
        mat->ke_txt = fl_mat->ke_txt;
        mat->kd_txt = fl_mat->kd_txt;
        mat->ks_txt = fl_mat->ks_txt;
        mat->rs_txt = fl_mat->rs_txt;
        sc->materials.push_back(mat);
    }

    // convert textures
    for (auto fl_txt : fl_scene->textures) {
        auto txt = new texture();
        txt->path = fl_txt->path;
        txt->width = fl_txt->width;
        txt->height = fl_txt->height;
        txt->width = fl_txt->width;
        txt->height = fl_txt->height;
        txt->ncomp = fl_txt->ncomp;
        txt->ldr = fl_txt->datab;
        txt->hdr = fl_txt->dataf;
        sc->textures.push_back(txt);
    }

    // done
    return sc.release();
}

//
// Save scene
//
inline void save_scene(const std::string& filename, const scene* sc) {
    auto ext = ycmd::get_extension(filename);
    if (ext == ".obj") {
        save_obj_scene(filename, sc);
    } else if (ext == ".gltf") {
        save_gltf_scene(filename, sc);
    } else {
        throw std::invalid_argument("unknown file type");
    }
}

//
// Load scene
//
inline scene* load_scene(const std::string& filename) {
    // declare scene
    auto sc = std::unique_ptr<scene>();

    // get extension
    auto ext = ycmd::get_extension(filename);
    if (ext == ".obj") {
        // load obj
        sc = std::unique_ptr<scene>(load_obj_scene(filename));
    } else if (ext == ".gltf") {
        // load obj
        sc = std::unique_ptr<scene>(load_gltf_scene(filename, false));
    } else if (ext == ".glb") {
        // load obj
        sc = std::unique_ptr<scene>(load_gltf_scene(filename, true));
    } else {
        throw std::invalid_argument("unknown file type");
    }

    // check textures and patch them up if needed
    for (auto txt : sc->textures) {
        if (txt->hdr.empty() && txt->ldr.empty()) {
            printf("unable to load texture %s\n", txt->path.c_str());
            txt->width = 1;
            txt->height = 1;
            txt->ncomp = 4;
            txt->ldr = {255, 255, 255, 255};
        }
    }

    // ensure normals
    for (auto sh : sc->shapes) {
        if (!sh->norm.empty()) continue;
        sh->norm.resize(sh->pos.size());
        yshape::compute_normals(
            (int)sh->points.size(), sh->points.data(), (int)sh->lines.size(),
            sh->lines.data(), (int)sh->triangles.size(), sh->triangles.data(),
            (int)sh->pos.size(), sh->pos.data(), sh->norm.data());
    }

    // ensure radius is necessary
    for (auto sh : sc->shapes) {
        if (sh->points.empty() && sh->lines.empty()) continue;
        if (!sh->radius.empty()) continue;
        sh->radius.resize(sh->pos.size(), 0.001f);
    }

    // make camera if not there
    if (!sc->cameras.size()) {
        // find scene bounds
        auto bbox = ym::invalid_bbox3f;
        for (auto sh : sc->shapes) {
            for (auto p : sh->pos)
                bbox +=
                    ym::transform_point((ym::frame3f)sh->frame, (ym::vec3f)p);
        }
        auto center = ym::center(bbox);
        auto bbox_size = ym::diagonal(bbox);
        auto bbox_msize =
            ym::max(bbox_size[0], ym::max(bbox_size[1], bbox_size[2]));
        // create camera
        auto cam = new yapp::camera();
        // set up camera
        auto camera_dir = ym::vec3f{1, 0.4f, 1};
        auto from = camera_dir * bbox_msize + center;
        auto to = center;
        auto up = ym::vec3f{0, 1, 0};
        cam->frame = ym::lookat_frame3(from, to, up);
        cam->ortho = false;
        cam->aspect = 16.0f / 9.0f;
        cam->yfov = 2 * atan(0.5f);
        cam->aperture = 0;
        cam->focus = ym::length(to - from);
        sc->cameras.push_back(cam);
    } else {
        auto bbox = ym::invalid_bbox3f;
        for (auto sh : sc->shapes) {
            for (auto p : sh->pos)
                bbox +=
                    ym::transform_point((ym::frame3f)sh->frame, (ym::vec3f)p);
        }
        for (auto cam : sc->cameras) {
            if (!cam->focus) {
                auto ddir = ym::dot((ym::vec3f)cam->frame[2],
                                    ym::center(bbox) -
                                        ym::pos((ym::frame3f)cam->frame));
                cam->focus = (ddir > 0) ? 1 : -ddir;
            }
        }
    }

    return sc.release();
}

struct params {
    std::string filename;
    std::string imfilename;

    yapp::scene* scene;

    ~params() {
        if (scene) delete scene;
    }
};

inline void init_params(params* pars, ycmd::parser* parser) {
    // params
    auto imfilename = ycmd::parse_opt<std::string>(parser, "--output", "-o",
                                                   "image filename", "out.hdr");
    auto filename = ycmd::parse_arg<std::string>(parser, "scene",
                                                 "scene filename", "", true);

    // filenames
    pars->filename = filename;
    pars->imfilename = imfilename;

    // loading scene
    if (filename.empty()) return;
    pars->scene = yapp::load_scene(filename);
}

}  // namespace

#endif
