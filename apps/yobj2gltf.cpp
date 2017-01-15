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

#include "../yocto/yocto_cmd.h"
#include "../yocto/yocto_gltf.h"
#include "../yocto/yocto_obj.h"

#include <algorithm>
#include <array>
#include <memory>

#define YOBJ2GLTF_VERBOSE

using float3 = std::array<float, 3>;
using float3x2 = std::array<std::array<float, 3>, 2>;

ygltf::fl_gltf* convert(const yobj::fl_obj* obj) {
    auto gltf = new ygltf::fl_gltf();

    // convert primitives
    for (auto oprim : obj->primitives) {
        auto gprim = new ygltf::fl_primitives();
        gprim->material = oprim->material;
        gprim->pos = oprim->pos;
        gprim->norm = oprim->norm;
        gprim->texcoord = oprim->texcoord;
        gprim->color = oprim->color;
        gprim->points = oprim->points;
        gprim->lines = oprim->lines;
        gprim->triangles = oprim->triangles;
        gltf->primitives.push_back(gprim);
    }

    // convert meshes
    for (auto omesh : obj->meshes) {
        auto gmesh = new ygltf::fl_mesh();
        gmesh->name = omesh->name;
        gmesh->primitives = omesh->primitives;
        gltf->meshes.push_back(gmesh);
    }

    // convert textures
    for (auto otxt : obj->textures) {
        auto gtxt = new ygltf::fl_texture();
        gtxt->path = otxt->path;
        gtxt->width = otxt->width;
        gtxt->height = otxt->height;
        gtxt->ncomp = otxt->ncomp;
        gtxt->dataf = otxt->dataf;
        gtxt->datab = otxt->datab;
        gltf->textures.push_back(gtxt);
    }

    // convert materials
    for (auto omat : obj->materials) {
        auto gmat = new ygltf::fl_material();
        gmat->name = omat->name;
        gmat->ke = omat->ke;
        gmat->kd = omat->kd;
        gmat->ks = omat->ks;
        gmat->rs = omat->rs;
        gmat->ke_txt = omat->ke_txt;
        gmat->kd_txt = omat->kd_txt;
        gmat->ks_txt = omat->ks_txt;
        gmat->rs_txt = omat->rs_txt;
        gltf->materials.push_back(gmat);
    }

    // convert cameras
    for (auto ocam : obj->cameras) {
        auto gcam = new ygltf::fl_camera();
        gcam->name = ocam->name;
        gcam->xform = ocam->xform;
        gcam->ortho = ocam->ortho;
        gcam->yfov = ocam->yfov;
        gcam->aspect = ocam->aspect;
        gltf->cameras.push_back(gcam);
    }

    // done
    return gltf;
}

template <typename T>
float3x2 compute_bounds(const T* obj) {
    auto bbox = float3x2();
    if (obj->primitives.empty()) return bbox;
    bbox[0] = obj->primitives[0]->pos[0];
    bbox[1] = obj->primitives[0]->pos[0];
    for (auto prim : obj->primitives) {
        for (auto p : prim->pos) {
            for (auto i = 0; i < 3; i++) {
                bbox[0][i] = std::min(bbox[0][i], p[i]);
                bbox[1][i] = std::max(bbox[1][i], p[i]);
            }
        }
    }
    return bbox;
}

template <typename T>
void print_xxx_info(const T* obj) {
    auto nobjs = (int)obj->meshes.size();
    auto ngroups = (int)obj->primitives.size();
    auto nverts = 0, nnorms = 0, ntexcoords = 0, npoints = 0, nlines = 0,
         ntriangles = 0;
    auto unique_prim = true;
    for (auto prim : obj->primitives) {
        nverts += prim->pos.size();
        nnorms += prim->norm.size();
        ntexcoords += prim->texcoord.size();
        npoints += prim->points.size();
        nlines += prim->lines.size();
        ntriangles += prim->triangles.size();
        if (((prim->points.empty() ? 0 : 1) + (prim->lines.empty() ? 0 : 1) +
             (prim->triangles.empty() ? 0 : 1)) <= 1)
            unique_prim = false;
    }

    auto bbox = compute_bounds(obj);
    auto bboxc =
        float3{(bbox[1][0] + bbox[0][0]) / 2, (bbox[1][1] + bbox[0][1]) / 2,
               (bbox[1][2] + bbox[0][2]) / 2};
    auto bboxs = float3{bbox[1][0] - bbox[0][0], bbox[1][1] - bbox[0][1],
                        bbox[1][2] - bbox[0][2]};

    printf("number of objects:    %d\n", nobjs);
    printf("number of groups:     %d\n", ngroups);
    printf("unique primitive:     %d\n", unique_prim);
    printf("number of vertices:   %d\n", nverts);
    printf("number of normals:    %d\n", nnorms);
    printf("number of texcoords:  %d\n", ntexcoords);
    printf("number of points:     %d\n", npoints);
    printf("number of lines:      %d\n", nlines);
    printf("number of triangles:  %d\n", ntriangles);
    printf("\n");
    printf("bbox min:    %g %g %g\n", bbox[0][0], bbox[0][1], bbox[0][2]);
    printf("bbox max:    %g %g %g\n", bbox[1][0], bbox[1][1], bbox[1][2]);
    printf("bbox center: %g %g %g\n", bboxc[0], bboxc[1], bboxc[2]);
    printf("bbox size:   %g %g %g\n", bboxs[0], bboxs[1], bboxs[2]);
    printf("\n");

#ifdef YOBJ2GLTF_VERBOSE
    auto pid = 0;
    for (auto prim : obj->primitives) {
        printf("primitive:            %d\n", pid++);
        printf("number of vertices:   %d\n", (int)prim->pos.size());
        printf("number of normals:    %d\n", (int)prim->norm.size());
        printf("number of texcoords:  %d\n", (int)prim->texcoord.size());
        printf("number of points:     %d\n", (int)prim->points.size());
        printf("number of lines:      %d\n", (int)prim->lines.size());
        printf("number of triangles:  %d\n", (int)prim->triangles.size());
        printf("\n");
    }
#endif
}

void print_obj_info(const yobj::fl_obj* obj) { print_xxx_info(obj); }
void print_gltf_info(const ygltf::fl_gltf* gltf) { print_xxx_info(gltf); }

void scale_obj(yobj::fl_obj* obj, float scale) {
    for (auto prim : obj->primitives) {
        for (auto& p : prim->pos) {
            p[0] *= scale;
            p[1] *= scale;
            p[2] *= scale;
        }
    }
}

void flipy_texcoord_obj(yobj::fl_obj* obj) {
    for (auto prim : obj->primitives) {
        if (prim->texcoord.empty()) continue;
        auto tmin = prim->texcoord[0][1], tmax = prim->texcoord[0][1];
        for (auto t : prim->texcoord) {
            tmin = std::min(tmin, t[1]);
            tmax = std::min(tmax, t[1]);
        }
        if (tmin >= 0 && tmax <= 1) {
            for (auto& t : prim->texcoord) t[1] = 1 - t[1];
        } else {
            for (auto& t : prim->texcoord) t[1] = 1 - t[1];
        }
    }
}

int main(int argc, char** argv) {
    // command line params
    auto parser = ycmd::make_parser(argc, argv, "converts obj to gltf");
    auto no_flipy_texcoord =
        ycmd::parse_flag(parser, "--no_flipy_texcoord", "",
                         "disable texcoord vertical flipping");
    auto scale =
        ycmd::parse_optf(parser, "--scale", "", "scale the model", 1.0f);
    auto print_info = ycmd::parse_flag(parser, "--print_info", "-i",
                                       "print information", false);
    auto validate = ycmd::parse_flag(parser, "--validate", "",
                                     "validate after saving", false);
    auto no_save = ycmd::parse_flag(parser, "--no_save", "-e",
                                    "exit without saving", false);
    auto filename_in =
        ycmd::parse_args(parser, "filename_in", "input filename", "", true);
    auto filename_out =
        ycmd::parse_args(parser, "filename_out", "output filename", "", false);
    ycmd::check_parser(parser);

    // set output filename if not present
    if (filename_out.empty()) {
        filename_out = ycmd::get_dirname(filename_in) +
                       ycmd::get_basename(filename_in) + ".gltf";
    }

    // load obj
    auto obj_ = yobj::load_obj(filename_in);
    auto obj = yobj::flatten_obj(obj_);
    delete obj_;

    // print information
    if (print_info) {
        printf("OBJ information ------------------------\n");
        print_obj_info(obj);
    }

    // scale
    if (scale != 1.0f) scale_obj(obj, scale);
    if (!no_flipy_texcoord) flipy_texcoord_obj(obj);

    // print infomation again if needed
    if (print_info and (scale != 1.0f || !no_flipy_texcoord)) {
        printf("OBJ post-correction information -------\n");
        print_obj_info(obj);
    }

    // convert
    auto gltf = convert(obj);

    // print information
    if (print_info) {
        printf("glTF information ---------------------\n");
        print_gltf_info(gltf);
    }

    // exit without saving
    if (no_save) return 0;

    // save gltf
    auto gltf_ =
        ygltf::unflatten_gltf(gltf, ycmd::get_basename(filename_out) + ".bin");
    ygltf::save_gltf(filename_out, gltf_, true, false, false);
    delete gltf_;

    // validate
    if (validate) {
        auto vgltf_ = ygltf::load_gltf(filename_out, true, false, false);
        auto vgltf = ygltf::flatten_gltf(vgltf_);
        if (print_info) {
            printf("glTF validate information -----------------\n");
            print_gltf_info(vgltf);
        }
    }

    // cleanup
    delete gltf;
    delete obj;
    return 0;
}
