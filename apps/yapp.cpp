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
// YAPP IMPLEMENTATION
//

#include "yapp.h"

#include "../yocto/yocto_math.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../yocto/stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "../yocto/stb_image.h"

#include "tinyply.h"

#include <fstream>

namespace yapp {

//
// Destructor
//
scene::~scene() {
    for (auto v : shapes) delete v;
    for (auto v : materials) delete v;
    for (auto v : textures) delete v;
    for (auto v : cameras) delete v;
    for (auto v : environments) delete v;
}

//
// backward compatible calls
//
int get_etype(const shape& shape) {
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

int get_nelems(const shape& shape) {
    auto et = get_etype(shape);
    if (et == 1) return (int)shape.points.size();
    if (et == 2) return (int)shape.lines.size();
    if (et == 3) return (int)shape.triangles.size();
    return 0;
}

const int* get_elems(const shape& shape) {
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
    auto fl_scene = std::unique_ptr<yobj::fl_obj>(yobj::flatten_obj(obj.get()));

    // cleanup
    obj.reset(nullptr);

    // load textures
    yobj::load_textures(fl_scene.get(), ycmd::get_dirname(filename), true);

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
    for (auto fl_mesh : fl_scene->meshes) {
        for (auto prim_id : fl_mesh->primitives) {
            auto fl_prim = fl_scene->primitives[prim_id];
            auto sh = new shape();
            sh->name = fl_mesh->name;
            sh->frame = ym::identity_frame3f;
            sh->matid = fl_prim->material;
            sh->pos = fl_prim->pos;
            sh->norm = fl_prim->norm;
            sh->texcoord = fl_prim->texcoord;
            sh->color = fl_prim->color;
            sh->radius = fl_prim->radius;
            sh->points = fl_prim->points;
            sh->lines = fl_prim->lines;
            sh->triangles = fl_prim->triangles;
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
    auto fl_scene = std::unique_ptr<yobj::fl_obj>(new yobj::fl_obj);

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
        auto fl_mesh = new yobj::fl_mesh();
        fl_mesh->name = shape->name;
        fl_mesh->primitives.push_back((int)fl_scene->primitives.size());
        auto fl_prim = new yobj::fl_primitives();
        fl_prim->material = shape->matid;
        fl_prim->pos.resize(shape->pos.size());
        for (auto i = 0; i < shape->pos.size(); i++) {
            fl_prim->pos[i] = ym::transform_point((ym::frame3f)shape->frame,
                                                  (ym::vec3f)shape->pos[i]);
        }
        fl_prim->norm.resize(shape->norm.size());
        for (auto i = 0; i < shape->pos.size(); i++) {
            fl_prim->norm[i] = ym::transform_direction(
                (ym::frame3f)shape->frame, (ym::vec3f)shape->norm[i]);
        }
        fl_prim->texcoord = std::vector<std::array<float, 2>>(
            shape->texcoord.begin(), shape->texcoord.end());
        fl_prim->color = std::vector<std::array<float, 3>>(shape->color.begin(),
                                                           shape->color.end());
        fl_prim->radius = shape->radius;
        fl_prim->points = shape->points;
        fl_prim->lines = std::vector<std::array<int, 2>>(shape->lines.begin(),
                                                         shape->lines.end());
        fl_prim->triangles = std::vector<std::array<int, 3>>(
            shape->triangles.begin(), shape->triangles.end());
        fl_scene->primitives.push_back(fl_prim);
        fl_scene->meshes.push_back(fl_mesh);
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
        auto fl_prim = new ygltf::fl_primitives();
        fl_prim->material = shape->matid;
        fl_prim->pos = std::vector<std::array<float, 3>>(shape->pos.begin(),
                                                         shape->pos.end());
        fl_prim->norm = std::vector<std::array<float, 3>>(shape->norm.begin(),
                                                          shape->norm.end());
        fl_prim->texcoord = std::vector<std::array<float, 2>>(
            shape->texcoord.begin(), shape->texcoord.end());
        fl_prim->color = std::vector<std::array<float, 3>>(shape->color.begin(),
                                                           shape->color.end());
        fl_prim->radius = shape->radius;
        fl_prim->points = shape->points;
        fl_prim->lines = std::vector<std::array<int, 2>>(shape->lines.begin(),
                                                         shape->lines.end());
        fl_prim->triangles = std::vector<std::array<int, 3>>(
            shape->triangles.begin(), shape->triangles.end());
        fl_scene->primitives.push_back(fl_prim);
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
        fl_scene.get(), ycmd::get_basename(filename) + ".bin"));
    ygltf::save_gltf(filename, gltf.get(), true, false, false);
}

//
// Load gltf scene
//
inline scene* load_gltf_scene(const std::string& filename, bool binary) {
    // load scene
    auto gltf = std::unique_ptr<ygltf::glTF_t>(
        (binary) ? ygltf::load_binary_gltf(filename, true, false, true, true)
                 : ygltf::load_gltf(filename, true, false, true, true));

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
        for (auto fl_prim_id : fl_mesh->primitives) {
            auto fl_prim = fl_scene->primitives.at(fl_prim_id);
            auto sh = new shape();
            sh->name = fl_mesh->name;
            sh->frame = ym::to_frame(ym::mat4f(fl_mesh->xform));
            sh->matid = fl_prim->material;
            sh->pos = fl_prim->pos;
            sh->norm = fl_prim->norm;
            sh->texcoord = fl_prim->texcoord;
            sh->color = fl_prim->color;
            sh->radius = fl_prim->radius;
            sh->points = fl_prim->points;
            sh->lines = fl_prim->lines;
            sh->triangles = fl_prim->triangles;
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
// Loads a scene from ply.
//
inline scene* load_ply_scene(const std::string& filename) {
    // preallocate vertex and element data
    auto sh = std::unique_ptr<shape>(new shape());
    sh->matid = 0;

    // Tinyply can and will throw exceptions at you!
    try {
        // Read the file and create a std::istringstream suitable
        std::ifstream ss(filename, std::ios::binary);

        // Parse the ASCII header fields
        tinyply::PlyFile file(ss);

        // vertex data
        auto pos = std::vector<float>();
        auto npos = file.request_properties_from_element("vertex",
                                                         {"x", "y", "z"}, pos);
        auto norm = std::vector<float>();
        auto nnorm = file.request_properties_from_element(
            "vertex", {"nx", "ny", "nz"}, norm);

        auto kd = std::vector<float>();
        auto nkd = file.request_properties_from_element(
            "vertex", {"kdr", "kdg", "kdb"}, kd);

        auto ks = std::vector<float>();
        auto nks = file.request_properties_from_element(
            "vertex", {"ksr", "ksg", "ksb"}, ks);

        auto rs = std::vector<float>();
        auto nrs = file.request_properties_from_element("vertex", {"rs"}, rs);

        // load triangle data
        std::vector<uint32_t> faces;
        auto ntriangle = file.request_properties_from_element(
            "face", {"vertex_indices"}, faces, 3);

        // Now populate the vectors...
        file.read(ss);

        // Set vertex data
        if (npos)
            sh->pos.assign((float3*)pos.data(), (float3*)pos.data() + npos);
        if (nnorm)
            sh->norm.assign((float3*)norm.data(), (float3*)norm.data() + nnorm);
        if (nkd) sh->kd.assign((float3*)kd.data(), (float3*)kd.data() + nkd);
        if (nks) sh->ks.assign((float3*)ks.data(), (float3*)ks.data() + nks);
        if (nrs) sh->rs.assign(rs.data(), rs.data() + nrs);
        if (ntriangle)
            sh->triangles.assign((int3*)faces.data(),
                                 (int3*)faces.data() + ntriangle);
    } catch (const std::exception& e) {
        throw;
    }

    // init scene
    auto sc = std::unique_ptr<scene>(new scene());

    // set shape
    sc->shapes.push_back(sh.release());

    // create material
    auto mat = new material();
    mat->name = "default";
    mat->ke = {0, 0, 0};
    mat->kd = {1, 1, 1};
    mat->ks = {1, 1, 1};
    mat->rs = 0.1f;
    mat->ke_txt = -1;
    mat->kd_txt = -1;
    mat->ks_txt = -1;
    mat->rs_txt = -1;
    sc->materials.push_back(mat);

    // done
    return sc.release();
}

//
// Save scene
//
void save_scene(const std::string& filename, const scene* sc) {
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
scene* load_scene(const std::string& filename, float scale) {
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
    } else if (ext == ".ply") {
        // load ply
        sc = std::unique_ptr<scene>(load_ply_scene(filename));
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

    // scale if necessary
    if (scale != 1.0f) {
        if (scale != 1.0f) {
            for (auto sh : sc->shapes) {
                for (auto& p : sh->pos) (ym::vec3f&)p *= scale;
            }
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

//
// Load env map
//
void load_envmap(scene* scn, const std::string& filename, float scale) {
    if (filename.empty()) return;
    // texture
    auto txt = new texture();
    txt->path = filename;
    auto pixels =
        stbi_loadf(filename.c_str(), &txt->width, &txt->height, &txt->ncomp, 0);
    txt->hdr.assign(pixels, pixels + txt->width * txt->height * txt->ncomp);
    free(pixels);
    scn->textures.push_back(txt);
    // material
    auto mat = new material();
    mat->name = "env_mat";
    mat->ke = {scale, scale, scale};
    mat->ke_txt = (int)scn->textures.size() - 1;
    scn->materials.push_back(mat);
    // environment
    auto env = new environment();
    env->name = "env";
    env->matid = (int)scn->materials.size() - 1;
    env->frame = ym::lookat_frame3(ym::vec3f{0, 0, 1}, ym::vec3f{0, 0, 0},
                                   ym::vec3f{0, 1, 0});
    scn->environments.push_back(env);
}

std::vector<int4> make_trace_blocks(int w, int h, int bs) {
    std::vector<int4> blocks;
    for (int j = 0; j < h; j += bs) {
        for (int i = 0; i < w; i += bs) {
            blocks.push_back({i, j, ym::min(bs, w - i), ym::min(bs, h - j)});
        }
    }
    return blocks;
}

void save_image(const std::string& filename, int width, int height,
                const float4* hdr, float exposure, float gamma,
                bool srgb_output) {
    auto ext = ycmd::get_extension(filename);
    auto tone_mapped = (const float4*)nullptr;
    auto tone_mapped_buffer = std::vector<ym::vec4f>();
    if (exposure != 0 || gamma != 1) {
        tone_mapped_buffer =
            std::vector<ym::vec4f>(width * height, {0, 0, 0, 0});
        ym::exposure_gamma(width, height, 4, (const float*)hdr,
                           (float*)tone_mapped_buffer.data(), exposure, gamma,
                           false);
        tone_mapped = tone_mapped_buffer.data();
    } else {
        tone_mapped = hdr;
    }
    if (ext == ".hdr") {
        stbi_write_hdr(filename.c_str(), width, height, 4,
                       (float*)tone_mapped->data());
    } else if (ext == ".png") {
        auto ldr = std::vector<ym::vec4b>(width * height, {0, 0, 0, 0});
        if (srgb_output) {
            ym::linear_to_srgb(width, height, 4, (const float*)tone_mapped,
                               (unsigned char*)ldr.data());
        } else {
            ym::linear_to_byte(width, height, 4, (const float*)tone_mapped,
                               (unsigned char*)ldr.data());
        }
        stbi_write_png(filename.c_str(), width, height, 4, ldr.data(),
                       width * 4);
    } else {
        printf("supports only hdr and png for image writing\n");
        return;
    }
}

ybvh::scene* make_bvh(const yapp::scene* scene) {
    auto scene_bvh = ybvh::make_scene((int)scene->shapes.size());
    auto sid = 0;
    for (auto shape : scene->shapes) {
        if (!shape->points.empty()) {
            ybvh::set_point_shape(scene_bvh, sid++, shape->frame,
                                  (int)shape->points.size(),
                                  shape->points.data(), (int)shape->pos.size(),
                                  shape->pos.data(), shape->radius.data());
        } else if (!shape->lines.empty()) {
            ybvh::set_line_shape(scene_bvh, sid++, shape->frame,
                                 (int)shape->lines.size(), shape->lines.data(),
                                 (int)shape->pos.size(), shape->pos.data(),
                                 shape->radius.data());

        } else if (!shape->triangles.empty()) {
            ybvh::set_triangle_shape(
                scene_bvh, sid++, shape->frame, (int)shape->triangles.size(),
                shape->triangles.data(), (int)shape->pos.size(),
                shape->pos.data(), shape->radius.data());

        } else {
            ybvh::set_point_shape(scene_bvh, sid++, shape->frame,
                                  (int)shape->pos.size(), shape->pos.data(),
                                  shape->radius.data());
        }
    }
    ybvh::build_bvh(scene_bvh);
    return scene_bvh;
}

ytrace::scene* make_trace_scene(const yapp::scene* scene,
                                const ybvh::scene* scene_bvh, int camera) {
    auto trace_scene = ytrace::make_scene(
        (int)scene->cameras.size(), (int)scene->shapes.size(),
        (int)scene->materials.size(), (int)scene->textures.size(),
        (int)scene->environments.size());

    auto cid = 0;
    for (auto cam : scene->cameras) {
        ytrace::set_camera(trace_scene, cid++, cam->frame, cam->yfov,
                           cam->aspect, cam->aperture, cam->focus);
    }

    auto tid = 0;
    for (auto txt : scene->textures) {
        if (!txt->hdr.empty()) {
            ytrace::set_texture(trace_scene, tid++, txt->width, txt->height,
                                txt->ncomp, txt->hdr.data());
        } else if (!txt->ldr.empty()) {
            ytrace::set_texture(trace_scene, tid++, txt->width, txt->height,
                                txt->ncomp, txt->ldr.data());
        } else
            assert(false);
    }

    auto eid = 0;
    for (auto env : scene->environments) {
        auto mat = scene->materials[env->matid];
        ytrace::set_environment(trace_scene, eid++, env->frame, mat->ke,
                                mat->ke_txt);
    }

    auto mid = 0;
    for (auto mat : scene->materials) {
        ytrace::set_material(trace_scene, mid++, mat->ke, mat->kd, mat->ks,
                             mat->rs, mat->ke_txt, mat->kd_txt, mat->ks_txt,
                             mat->ks_txt);
    }

    auto sid = 0;
    for (auto shape : scene->shapes) {
        if (!shape->points.empty()) {
            ytrace::set_point_shape(trace_scene, sid++, shape->frame,
                                    shape->matid, (int)shape->points.size(),
                                    shape->points.data(),
                                    (int)shape->pos.size(), shape->pos.data(),
                                    shape->norm.data(), shape->texcoord.data(),
                                    shape->color.data(), shape->radius.data());
        } else if (!shape->lines.empty()) {
            ytrace::set_line_shape(trace_scene, sid++, shape->frame,
                                   shape->matid, (int)shape->lines.size(),
                                   shape->lines.data(), (int)shape->pos.size(),
                                   shape->pos.data(), shape->norm.data(),
                                   shape->texcoord.data(), shape->color.data(),
                                   shape->radius.data());

        } else if (!shape->triangles.empty()) {
            ytrace::set_triangle_shape(
                trace_scene, sid++, shape->frame, shape->matid,
                (int)shape->triangles.size(), shape->triangles.data(),
                (int)shape->pos.size(), shape->pos.data(), shape->norm.data(),
                shape->texcoord.data(), shape->color.data());

        } else {
        }
        if (!shape->ke.empty() || !shape->kd.empty() || !shape->ks.empty() ||
            !shape->rs.empty()) {
            ytrace::set_vert_material(trace_scene, sid - 1, shape->ke.data(),
                                      shape->kd.data(), shape->ks.data(),
                                      shape->rs.data());
        }
    }

    ytrace::set_intersection_callbacks(
        trace_scene, (void*)scene_bvh,
        [](auto ctx, auto o, auto d, auto tmin, auto tmax) {
            auto scene_bvh = (ybvh::scene*)ctx;
            auto isec = ybvh::intersect_ray(scene_bvh, o, d, tmin, tmax, false);
            auto ipt = ytrace::intersect_point();
            ipt.dist = isec.dist;
            ipt.sid = isec.sid;
            ipt.eid = isec.eid;
            ipt.euv = {isec.euv[0], isec.euv[1], isec.euv[2]};
            return ipt;
        },
        [](auto ctx, auto o, auto d, auto tmin, auto tmax) {
            auto scene_bvh = (ybvh::scene*)ctx;
            return (bool)ybvh::intersect_ray(scene_bvh, o, d, tmin, tmax, true);
        });

    ytrace::set_logging_callbacks(trace_scene, nullptr, ycmd::log_msgfv);

    ytrace::init_lights(trace_scene);

    return trace_scene;
}

#if 0
ybvh::scene* make_bvh(const yapp::scene* scene) {
    auto scene_bvh = ybvh::make_scene((int)scene->shapes.size());
    auto sid = 0;
    for (auto shape : scene->shapes) {
        if (!shape->points.empty()) {
            ybvh::set_point_shape(scene_bvh, sid++, shape->frame,
                                  (int)shape->points.size(),
                                  shape->points.data(), (int)shape->pos.size(),
                                  shape->pos.data(), shape->radius.data());
        } else if (!shape->lines.empty()) {
            ybvh::set_line_shape(scene_bvh, sid++, shape->frame,
                                 (int)shape->lines.size(), shape->lines.data(),
                                 (int)shape->pos.size(), shape->pos.data(),
                                 shape->radius.data());

        } else if (!shape->triangles.empty()) {
            ybvh::set_triangle_shape(
                scene_bvh, sid++, shape->frame, (int)shape->triangles.size(),
                shape->triangles.data(), (int)shape->pos.size(),
                shape->pos.data(), shape->radius.data());

        } else {
            ybvh::set_point_shape(scene_bvh, sid++, shape->frame,
                                  (int)shape->pos.size(), shape->pos.data(),
                                  shape->radius.data());
        }
    }
    ybvh::build_bvh(scene_bvh);
    return scene_bvh;
}
#endif

ysym::scene* make_rigid_scene(const yapp::scene* scene,
                              ybvh::scene* scene_bvh) {
    // allocate scene
    auto rigid_scene = ysym::make_scene((int)scene->shapes.size());

    // add each shape
    auto sid = 0;
    for (auto shape : scene->shapes) {
        auto mat = scene->materials[shape->matid];
        auto density =
            (shape->name != "floor" && ym::length((ym::vec3f)mat->ke) == 0 &&
             !shape->triangles.empty())
                ? 1.0f
                : 0.0f;
        ysym::set_body(rigid_scene, sid++, shape->frame, {0, 0, 0}, {0, 0, 0},
                       density, (int)shape->triangles.size(),
                       shape->triangles.data(), (int)shape->pos.size(),
                       shape->pos.data());
    }

    // set up final bvh
    scene_bvh = make_bvh(scene);

    // setup collisions
    ysym::set_overlap_callbacks(
        rigid_scene, scene_bvh,
        [](auto ctx, std::vector<ysym::int2>* overlaps) {
            auto scene_bvh = (ybvh::scene*)ctx;
            ybvh::overlap_shape_bounds(scene_bvh, scene_bvh, false, true, true,
                                       overlaps);
        },
        [](auto ctx, int sid, const ysym::float3& pt, float max_dist) {
            auto scene_bvh = (ybvh::scene*)ctx;
            auto overlap =
                ybvh::overlap_point(scene_bvh, sid, pt, max_dist, false);
            return *(ysym::overlap_point*)&overlap;
        },
        [](auto ctx, int sid1, int sid2, float max_dist,
           std::vector<std::pair<ysym::overlap_point, ysym::int2>>* overlaps_) {
            auto scene_bvh = (ybvh::scene*)ctx;
            auto overlaps =
                (std::vector<std::pair<ybvh::point, ybvh::int2>>*)overlaps_;
            ybvh::overlap_verts(scene_bvh, scene_bvh, sid1, sid2, true,
                                max_dist, true, overlaps);
        },
        [](auto ctx, auto rigid_scene, int nshapes) {
            auto scene_bvh = (ybvh::scene*)ctx;
            for (auto sid = 0; sid < nshapes; sid++) {
                ybvh::set_shape_frame(scene_bvh, sid,
                                      ysym::get_body_frame(rigid_scene, sid));
            }
            ybvh::refit_bvh(scene_bvh);
        });

    // initialize
    ysym::init_simulation(rigid_scene);

    return rigid_scene;
}

void simulate_step(yapp::scene* scene, ysym::scene* rigid_scene, float dt) {
    ysym::advance_simulation(rigid_scene, dt);
    for (auto sid = 0; sid < scene->shapes.size(); sid++) {
        scene->shapes[sid]->frame = ysym::get_body_frame(rigid_scene, sid);
    }
}

params* init_params(const std::string& help, int argc, char** argv,
                    bool trace_params, bool sym_params, bool shade_params,
                    bool ui_params) {
    static auto rtype_names = std::vector<std::pair<std::string, int>>{
        {"default", (int)ytrace::rng_type::def},
        {"uniform", (int)ytrace::rng_type::uniform},
        {"stratified", (int)ytrace::rng_type::stratified},
        {"cmjs", (int)ytrace::rng_type::cmjs}};
    static auto stype_names = std::vector<std::pair<std::string, int>>{
        {"default", (int)ytrace::shader_type::def},
        {"eye", (int)ytrace::shader_type::eyelight},
        {"direct", (int)ytrace::shader_type::direct},
        {"path", (int)ytrace::shader_type::pathtrace}};

    // parser
    auto parser = ycmd::make_parser(argc, argv, help.c_str());

    // parameters
    auto pars = new params();

    // render
    if (trace_params || shade_params) {
        pars->exposure = ycmd::parse_optf(parser, "--exposure", "-e",
                                          "hdr image exposure", 0);
        pars->gamma =
            ycmd::parse_optf(parser, "--gamma", "-g", "hdr image gamma", 1);
        pars->srgb =
            ycmd::parse_optb(parser, "--srgb", "", "hdr srgb output", true);
        auto aspect = ycmd::parse_optf(parser, "--aspect", "-a", "image aspect",
                                       16.0f / 9.0f);
        auto res = ycmd::parse_opti(parser, "--resolution", "-r",
                                    "image resolution", 720);
        pars->render_params.camera_id =
            ycmd::parse_opti(parser, "--camera", "-C", "camera", 0);
        pars->envmap_filename = ycmd::parse_opts(parser, "--envmap_filename",
                                                 "", "environment map", "");
        pars->envmap_scale = ycmd::parse_optf(parser, "--envmap_scale", "",
                                              "environment map scale", 1);
        pars->save_progressive =
            ycmd::parse_opti(parser, "--save_progressive", "",
                             "frames to save progressive images", 0);
        auto amb =
            ycmd::parse_optf(parser, "--ambient", "", "ambient factor", 0);

        pars->width = (int)std::round(aspect * res);
        pars->height = res;
        pars->render_params.amb = {amb, amb, amb};
    }

    if (shade_params) {
        auto camera_lights = ycmd::parse_flag(parser, "--camera_lights", "-c",
                                              "enable camera lights", false);
        pars->render_params.stype = (camera_lights)
                                        ? ytrace::shader_type::eyelight
                                        : ytrace::shader_type::direct;
    }

    // params
    if (trace_params) {
        pars->render_params.rtype = (ytrace::rng_type)ycmd::parse_opte(
            parser, "--random", "", "random type", (int)ytrace::rng_type::def,
            rtype_names);
        pars->render_params.stype = (ytrace::shader_type)ycmd::parse_opte(
            parser, "--integrator", "-i", "integrator type",
            (int)ytrace::shader_type::def, stype_names);
        auto camera_lights = ycmd::parse_flag(parser, "--camera_lights", "-c",
                                              "enable camera lights", false);
        pars->nthreads = ycmd::parse_opti(
            parser, "--threads", "-t", "number of threads [0 for default]", 0);
        pars->block_size =
            ycmd::parse_opt<int>(parser, "--block_size", "", "block size", 32);
        pars->render_params.nsamples =
            ycmd::parse_opti(parser, "--samples", "-s", "image samples", 256);

        if (camera_lights) {
            pars->render_params.stype = ytrace::shader_type::eyelight;
        }
    }

    if (sym_params) {
        pars->dt = ycmd::parse_optf(parser, "--delta_time", "-dt", "delta time",
                                    1 / 60.0f);
        pars->nframes = ycmd::parse_opti(parser, "--nframes", "-n",
                                         "number of frames", 1000);
        pars->outfilename = ycmd::parse_opts(parser, "--output", "-o",
                                             "output filename", "out.%04d.obj");
    }

    if (ui_params) {
        pars->no_ui =
            ycmd::parse_flag(parser, "--no-ui", "", "run without ui", false);
        pars->legacy_gl = ycmd::parse_flag(parser, "--legacy_opengl", "-L",
                                           "uses legacy OpenGL", false);
    }

    // params
    pars->scene_scale =
        ycmd::parse_optf(parser, "--scale", "", "scale scene", 1.0f);
    pars->imfilename =
        ycmd::parse_opts(parser, "--output", "-o", "image filename", "out.hdr");
    pars->filename =
        ycmd::parse_args(parser, "scene", "scene filename", "", true);

    // check parsing
    ycmd::check_parser(parser);

    // done
    return pars;
}

//
// Logging
//
void set_default_loggers() {
    auto loggers = ycmd::get_default_loggers();
    loggers->push_back(ycmd::make_stdout_logger());
    loggers->push_back(
        ycmd::make_file_logger("yocto.log", true, ycmd::log_level_verbose));
}

}  // namespace
