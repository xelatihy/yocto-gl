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

#include "../yocto/yocto_gltf.h"
#include "../yocto/yocto_obj.h"
#include "../yocto/yocto_utils.h"

#include <algorithm>
#include <array>
#include <memory>

#define YOBJ2GLTF_VERBOSE

template <typename T>
static inline int index(const std::vector<T*>& vec, T* val) {
    auto pos = std::find(vec.begin(), vec.end(), val);
    if (pos != vec.end()) return (int)(pos - vec.begin());
    return -1;
}

ygltf::scene_group* obj2gltf(const yobj::scene* obj, bool add_scene) {
    auto gltf = new ygltf::scene_group();

    // convert textures
    for (auto otxt : obj->textures) {
        auto gtxt = new ygltf::texture();
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
        auto gmat = new ygltf::material();
        gmat->name = omat->name;
        gmat->emission = omat->ke;
        gmat->emission_txt =
            (index(obj->textures, omat->ke_txt) < 0) ?
                nullptr :
                gltf->textures[index(obj->textures, omat->ke_txt)];
        gmat->metallic_roughness = new ygltf::material_metallic_rooughness();
        auto gmr = gmat->metallic_roughness;
        gmr->base = omat->kd;
        gmr->opacity = omat->opacity;
        gmr->metallic = 0;
        gmr->roughness = omat->rs;
        gmr->base_txt = (index(obj->textures, omat->kd_txt) < 0) ?
                            nullptr :
                            gltf->textures[index(obj->textures, omat->kd_txt)];
        gltf->materials.push_back(gmat);
    }

    // convert meshes
    for (auto omesh : obj->meshes) {
        auto gmesh = new ygltf::mesh();
        gmesh->name = omesh->name;
        for (auto oprim : omesh->shapes) {
            auto gprim = new ygltf::shape();
            gprim->material =
                (index(obj->materials, oprim->mat) < 0) ?
                    nullptr :
                    gltf->materials[index(obj->materials, oprim->mat)];
            gprim->pos = oprim->pos;
            gprim->norm = oprim->norm;
            gprim->texcoord = oprim->texcoord;
            gprim->color = oprim->color;
            gprim->points = oprim->points;
            gprim->lines = oprim->lines;
            gprim->triangles = oprim->triangles;
            gmesh->shapes.push_back(gprim);
        }
        gltf->meshes.push_back(gmesh);
    }

    if (add_scene) {
        // init nodes
        auto scn = new ygltf::scene();
        scn->name = "scene";
        gltf->default_scene = scn;
        gltf->scenes.push_back(scn);

        // convert instances
        if (obj->instances.empty()) {
            for (auto msh : gltf->meshes) {
                auto gnode = new ygltf::node();
                gnode->name = msh->name;
                gnode->mesh = msh;
                scn->nodes.push_back(gnode);
                gltf->root_nodes.push_back(gnode);
                gltf->nodes.push_back(gnode);
            }
        } else {
            for (auto oist : obj->instances) {
                auto gnode = new ygltf::node();
                gnode->name = oist->name;
                gnode->matrix = oist->xform;
                gnode->mesh = gltf->meshes[index(obj->meshes, oist->mesh)];
                gltf->root_nodes.push_back(gnode);
                scn->nodes.push_back(gnode);
                gltf->nodes.push_back(gnode);
            }
        }

        // convert cameras
        if (obj->cameras.empty()) {
            ygltf::add_default_cameras(gltf);
        } else {
            // TODO: convert cameras
            for (auto ocam : obj->cameras) {
                auto gcam = new ygltf::camera();
                gcam->name = ocam->name;
                gcam->ortho = ocam->ortho;
                gcam->yfov = ocam->yfov;
                gcam->aspect = ocam->aspect;
                gcam->focus = ocam->focus;
                gcam->aperture = ocam->aperture;
                gltf->cameras.push_back(gcam);
                auto gnode = new ygltf::node();
                gnode->name = ocam->name;
                gnode->matrix = ocam->xform;
                gnode->camera = gcam;
                gltf->root_nodes.push_back(gnode);
                scn->nodes.push_back(gnode);
                gltf->nodes.push_back(gnode);
            }
        }
    }

    // done
    return gltf;
}

void print_obj_info(const yobj::scene* obj) {
    auto nobjs = (int)obj->meshes.size();
    auto nverts = 0, nnorms = 0, ntexcoords = 0, npoints = 0, nlines = 0,
         ntriangles = 0, ngroups = 0;
    auto unique_prim = true;
    for (auto msh : obj->meshes) {
        for (auto prim : msh->shapes) {
            ngroups += 1;
            nverts += prim->pos.size();
            nnorms += prim->norm.size();
            ntexcoords += prim->texcoord.size();
            npoints += prim->points.size();
            nlines += prim->lines.size();
            ntriangles += prim->triangles.size();
            if (((prim->points.empty() ? 0 : 1) +
                    (prim->lines.empty() ? 0 : 1) +
                    (prim->triangles.empty() ? 0 : 1)) <= 1)
                unique_prim = false;
        }
    }

    auto bbox = compute_scene_bounds(obj);
    auto bboxc = ym::vec3f{(bbox[1][0] + bbox[0][0]) / 2,
        (bbox[1][1] + bbox[0][1]) / 2, (bbox[1][2] + bbox[0][2]) / 2};
    auto bboxs = ym::vec3f{bbox[1][0] - bbox[0][0], bbox[1][1] - bbox[0][1],
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
    auto pid = 0, mid = 0;
    for (auto mesh : obj->meshes) {
        printf("mesh:            %d\n", mid++);
        for (auto prim : mesh->shapes) {
            printf("primitive:            %d\n", pid++);
            printf("number of vertices:   %d\n", (int)prim->pos.size());
            printf("number of normals:    %d\n", (int)prim->norm.size());
            printf("number of texcoords:  %d\n", (int)prim->texcoord.size());
            printf("number of points:     %d\n", (int)prim->points.size());
            printf("number of lines:      %d\n", (int)prim->lines.size());
            printf("number of triangles:  %d\n", (int)prim->triangles.size());
            printf("\n");
        }
    }
#endif
}

void print_gltf_info(const ygltf::scene_group* gltf) {
    auto nobjs = (int)gltf->meshes.size();
    auto nverts = 0, nnorms = 0, ntexcoords = 0, npoints = 0, nlines = 0,
         ntriangles = 0, ngroups = 0;
    auto unique_prim = true;
    for (auto msh : gltf->meshes) {
        for (auto prim : msh->shapes) {
            nverts += prim->pos.size();
            nnorms += prim->norm.size();
            ntexcoords += prim->texcoord.size();
            npoints += prim->points.size();
            nlines += prim->lines.size();
            ntriangles += prim->triangles.size();
            if (((prim->points.empty() ? 0 : 1) +
                    (prim->lines.empty() ? 0 : 1) +
                    (prim->triangles.empty() ? 0 : 1)) <= 1)
                unique_prim = false;
        }
    }

    auto bbox = compute_scene_bounds(gltf);
    auto bboxc = ym::vec3f{(bbox[1][0] + bbox[0][0]) / 2,
        (bbox[1][1] + bbox[0][1]) / 2, (bbox[1][2] + bbox[0][2]) / 2};
    auto bboxs = ym::vec3f{bbox[1][0] - bbox[0][0], bbox[1][1] - bbox[0][1],
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
    auto pid = 0, mid = 0;
    for (auto msh : gltf->meshes) {
        printf("mesh:            %d\n", mid++);
        for (auto prim : msh->shapes) {
            printf("primitive:            %d\n", pid++);
            printf("number of vertices:   %d\n", (int)prim->pos.size());
            printf("number of normals:    %d\n", (int)prim->norm.size());
            printf("number of texcoords:  %d\n", (int)prim->texcoord.size());
            printf("number of points:     %d\n", (int)prim->points.size());
            printf("number of lines:      %d\n", (int)prim->lines.size());
            printf("number of triangles:  %d\n", (int)prim->triangles.size());
            printf("\n");
        }
    }
#endif
}

void scale_obj(yobj::scene* obj, float scale) {
    for (auto msh : obj->meshes) {
        for (auto prim : msh->shapes) {
            for (auto& p : prim->pos) {
                p[0] *= scale;
                p[1] *= scale;
                p[2] *= scale;
            }
        }
    }
}

void flipy_texcoord_obj(yobj::scene* obj) {
    for (auto msh : obj->meshes) {
        for (auto prim : msh->shapes) {
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
}

int main(int argc, char** argv) {
    // command line params
    auto parser = yu::cmdline::make_parser(argc, argv, "converts obj to gltf");
    auto no_flipy_texcoord = parse_flag(parser, "--no_flipy_texcoord", "",
        "disable texcoord vertical flipping");
    auto scale = parse_optf(parser, "--scale", "", "scale the model", 1.0f);
    auto print_info =
        parse_flag(parser, "--print_info", "-i", "print information", false);
    auto validate =
        parse_flag(parser, "--validate", "", "validate after saving", false);
    auto no_save =
        parse_flag(parser, "--no_save", "-e", "exit without saving", false);
    auto filename_in =
        parse_args(parser, "filename_in", "input filename", "", true);
    auto filename_out =
        parse_args(parser, "filename_out", "output filename", "", false);
    check_parser(parser);

    // set output filename if not present
    if (filename_out.empty()) {
        filename_out = yu::path::get_dirname(filename_in) +
                       yu::path::get_basename(filename_in) + ".gltf";
    }

    // load obj
    auto obj = yobj::load_scene(filename_in, false);

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
    auto gltf = obj2gltf(obj, false);

    // print information
    if (print_info) {
        printf("glTF information ---------------------\n");
        print_gltf_info(gltf);
    }

    // exit without saving
    if (no_save) return 0;

    // save gltf
    ygltf::save_scenes(filename_out, gltf, false);

    // validate
    if (validate) {
        auto vgltf_ = ygltf::load_gltf(filename_out, true, false);
        auto vgltf = ygltf::gltf_to_scenes(vgltf_);
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
