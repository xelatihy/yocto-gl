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

#include <memory>

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

int main(int argc, char** argv) {
    // command line params
    auto parser = ycmd::make_parser(argc, argv, "converts obj to gltf");
    auto filename_in = ycmd::parse_arg<std::string>(parser, "filename_in",
                                                    "input filename", "", true);
    auto filename_out = ycmd::parse_arg<std::string>(
        parser, "filename_out", "output filename", "", false);
    ycmd::check_parser(parser);

    // set output filename if not present
    if (filename_out.empty()) {
        filename_out = ycmd::get_basename(filename_in) + ".gltf";
    }

    // load obj
    auto obj_ = yobj::load_obj(filename_in);
    auto obj = yobj::flatten_obj(obj_);
    delete obj_;

    // convert
    auto gltf = convert(obj);

    // save gltf
    auto gltf_ =
        ygltf::unflatten_gltf(gltf, ycmd::get_basename(filename_out) + ".bin");
    ygltf::save_gltf(filename_out, gltf_, true, false, false);
    delete gltf_;

    // cleanup
    delete gltf;
    delete obj;
}
