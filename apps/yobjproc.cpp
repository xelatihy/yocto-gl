//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
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

#include "../yocto/ygl.h"
#include "../yocto/yglio.h"

std::string to_string(const ygl::obj_vertex& v) {
    auto s = std::to_string(v.pos);
    if (v.texcoord) {
        s += "/" + std::to_string(v.texcoord);
        if (v.norm) s += "/" + std::to_string(v.norm);
    } else {
        if (v.norm) s += "//" + std::to_string(v.norm);
    }
    return s;
}

int main(int argc, char* argv[]) {
    // parse command line
    auto parser = ygl::make_cmdline_parser(
        argc, argv, "Process obj files directly", "yobjproc");
    auto translation = ygl::parse_vec3f(
        parser, "--translation,-t", ygl::zero3f, "translation");
    auto scale =
        ygl::parse_vec3f(parser, "--scale,-s", ygl::vec3f{1, 1, 1}, "scale");
    auto print_info =
        ygl::parse_flag(parser, "--print-info,-i", false, "print obj info");
    auto output = ygl::parse_string(
        parser, "--output,-o", "out.obj", "output obj scene", true);
    auto filename = ygl::parse_string(
        parser, "filename", "img.obj", "input obj filename", true);
    ygl::check_cmdline(parser);

    // prepare stats
    auto npos = 0, nnorm = 0, ntexcoord = 0;
    auto nobjects = 0, ngroups = 0, nusemtl = 0, nmaterials = 0;
    auto nfaces = 0, nplines = 0, nlines = 0, nppoints = 0, npoints = 0;
    auto ntriangles = 0, nquads = 0, npolys = 0;
    auto bbox = ygl::invalid_bbox3f;
    auto tbox = ygl::invalid_bbox3f;

    // prepare file to output
    auto fs = fopen(output.c_str(), "wt");
    if (!fs) {
        printf("cannot open file %s\n", output.c_str());
        exit(1);
    }

    // obj callbacks
    auto cb = ygl::obj_callbacks();

    // vertex data
    cb.vert = [&](const ygl::vec3f& v) {
        auto tv = translation + scale * v;
        fprintf(fs, "v  %6g %6g %6g\n", tv.x, tv.y, tv.z);
        bbox += v;
        tbox += tv;
        npos += 1;
    };
    cb.norm = [&](const ygl::vec3f& v) {
        fprintf(fs, "vn %6g %6g %6g\n", v.x, v.y, v.z);
        nnorm += 1;
    };
    cb.texcoord = [&](const ygl::vec2f& v) {
        fprintf(fs, "vt %6g %6g\n", v.x, v.y);
        ntexcoord += 1;
    };
    cb.face = [&](const std::vector<ygl::obj_vertex>& verts) {
        fprintf(fs, "f");
        for (auto v : verts) fprintf(fs, " %s", to_string(v).c_str());
        fprintf(fs, "\n");
        nfaces += 1;
        if (verts.size() == 3) ntriangles += 1;
        if (verts.size() == 4) nquads += 1;
        if (verts.size() > 4) npolys += 1;
    };
    cb.line = [&](const std::vector<ygl::obj_vertex>& verts) {
        fprintf(fs, "l");
        for (auto v : verts) fprintf(fs, " %s", to_string(v).c_str());
        fprintf(fs, "\n");
        nplines += 1;
        nlines += (int)verts.size() - 1;
    };
    cb.point = [&](const std::vector<ygl::obj_vertex>& verts) {
        fprintf(fs, "p");
        for (auto v : verts) fprintf(fs, " %s", to_string(v).c_str());
        fprintf(fs, "\n");
        nppoints += 1;
        npoints += (int)verts.size();
    };
    cb.object = [&](const std::string& name) {
        fprintf(fs, "o %s\n", name.c_str());
        nobjects += 1;
    };
    cb.group = [&](const std::string& name) {
        fprintf(fs, "g %s\n", name.c_str());
        ngroups += 1;
    };
    cb.usemtl = [&](const std::string& name) {
        fprintf(fs, "usemtl %s\n", name.c_str());
        nusemtl += 1;
    };
    cb.mtllib = [&](const std::string& name) {
        fprintf(fs, "mtllib %s\n", name.c_str());
    };
    cb.material = [&](auto&) { nmaterials += 1; };

    // parse and write
    try {
        ygl::load_obj(filename, cb, false, false);
    } catch (std::exception& e) {
        printf(
            "could not load obj %s\nerror: %s\n", filename.c_str(), e.what());
        exit(1);
    }

    // print info
    if (print_info) {
        auto center = (bbox.max + bbox.min) / 2;
        auto size = bbox.max - bbox.min;
        printf(
            "bbox min: % 6g % 6g % 6g\n", bbox.min.x, bbox.min.y, bbox.min.z);
        printf(
            "bbox max: % 6g % 6g % 6g\n", bbox.max.x, bbox.max.y, bbox.max.z);
        printf("bbox cen: % 6g % 6g % 6g\n", center.x, center.y, center.z);
        printf("bbox siz: % 6g % 6g % 6g\n", size.x, size.y, size.z);
        if (translation != ygl::zero3f || scale != ygl::vec3f{1, 1, 1}) {
            auto center = (tbox.max + tbox.min) / 2;
            auto size = tbox.max - tbox.min;
            printf("tbox min: % 6g % 6g % 6g\n", tbox.min.x, tbox.min.y,
                tbox.min.z);
            printf("tbox max: % 6g % 6g % 6g\n", tbox.max.x, tbox.max.y,
                tbox.max.z);
            printf("tbox cen: % 6g % 6g % 6g\n", center.x, center.y, center.z);
            printf("tbox siz: % 6g % 6g % 6g\n", size.x, size.y, size.z);
        }
        printf("faces:    % 6d % 6d % 6d % 6d\n", nfaces, ntriangles, nquads,
            npolys);
        printf("lines:    % 6d % 6d\n", nplines, nlines);
        printf("point:    % 6d % 6d\n", nppoints, npoints);
        printf("objs:     % 6d % 6d\n", nobjects, ngroups);
        printf("mats:     % 6d % 6d\n", nmaterials, nusemtl);
    }

    // done
    fclose(fs);
}
