//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
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

#include <yocto/yocto_commonio.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_shape.h>
using namespace yocto;

#include "ext/filesystem.hpp"
namespace sfs = ghc::filesystem;

using namespace std::string_literals;

// Shape presets used ofr testing.
bool make_shape_preset(vector<int>& points, vector<vec2i>& lines,
    vector<vec3i>& triangles, vector<vec4i>& quads, vector<vec4i>& quadspos,
    vector<vec4i>& quadsnorm, vector<vec4i>& quadstexcoord,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<vec3f>& colors, vector<float>& radius, const string& type,
    string& error) {
  auto set_quads = [&](quads_shape&& shape) {
    quads     = shape.quads;
    positions = shape.positions;
    normals   = shape.normals;
    texcoords = shape.texcoords;
  };
  auto set_triangles = [&](triangles_shape&& shape) {
    triangles = shape.triangles;
    positions = shape.positions;
    normals   = shape.normals;
    texcoords = shape.texcoords;
  };
  auto set_lines = [&](lines_shape&& shape) {
    lines     = shape.lines;
    positions = shape.positions;
    normals   = shape.normals;
    texcoords = shape.texcoords;
    radius    = shape.radius;
  };
  auto set_points = [&](points_shape&& shape) {
    points    = shape.points;
    positions = shape.positions;
    normals   = shape.normals;
    texcoords = shape.texcoords;
    radius    = shape.radius;
  };
  auto set_fvquads = [&](quads_fvshape&& shape) {
    quadspos      = shape.quadspos;
    quadsnorm     = shape.quadsnorm;
    quadstexcoord = shape.quadstexcoord;
    positions     = shape.positions;
    normals       = shape.normals;
    texcoords     = shape.texcoords;
  };

  if (type == "default-quad") {
    set_quads(make_rect());
  } else if (type == "default-quady") {
    set_quads(make_recty());
  } else if (type == "default-cube") {
    set_quads(make_box());
  } else if (type == "default-cube-rounded") {
    set_quads(make_rounded_box());
  } else if (type == "default-sphere") {
    set_quads(make_sphere());
  } else if (type == "default-disk") {
    set_quads(make_disk());
  } else if (type == "default-disk-bulged") {
    set_quads(make_bulged_disk());
  } else if (type == "default-quad-bulged") {
    set_quads(make_bulged_rect());
  } else if (type == "default-uvsphere") {
    set_quads(make_uvsphere());
  } else if (type == "default-uvsphere-flipcap") {
    set_quads(make_capped_uvsphere());
  } else if (type == "default-uvdisk") {
    set_quads(make_uvdisk());
  } else if (type == "default-uvcylinder") {
    set_quads(make_uvcylinder());
  } else if (type == "default-uvcylinder-rounded") {
    set_quads(make_rounded_uvcylinder({32, 32, 32}));
  } else if (type == "default-geosphere") {
    set_triangles(make_geosphere());
  } else if (type == "default-floor") {
    set_quads(make_floor());
  } else if (type == "default-floor-bent") {
    set_quads(make_bent_floor());
  } else if (type == "default-matball") {
    set_quads(make_sphere());
  } else if (type == "default-hairball") {
    auto base = make_sphere(pow2(5), 0.8);
    set_lines(make_hair(base, {4, 65536}, {0.2, 0.2}, {0.002, 0.001}));
  } else if (type == "default-hairball-interior") {
    set_quads(make_sphere(pow2(5), 0.8));
  } else if (type == "default-suzanne") {
    set_quads(make_monkey());
  } else if (type == "default-cube-facevarying") {
    set_fvquads(make_fvbox());
  } else if (type == "default-sphere-facevarying") {
    set_fvquads(make_fvsphere());
  } else if (type == "default-quady-displaced") {
    set_quads(make_recty({256, 256}));
  } else if (type == "default-sphere-displaced") {
    set_quads(make_sphere(128));
  } else if (type == "test-cube") {
    set_quads(make_rounded_box(
        {32, 32, 32}, {0.075f, 0.075f, 0.075f}, {1, 1, 1}, 0.3 * 0.075f));
    for (auto& p : positions) p += {0, 0.075, 0};
  } else if (type == "test-uvsphere") {
    set_quads(make_uvsphere({32, 32}, 0.075));
    for (auto& p : positions) p += {0, 0.075, 0};
  } else if (type == "test-uvsphere-flipcap") {
    set_quads(make_capped_uvsphere({32, 32}, 0.075, {1, 1}, 0.3 * 0.075));
    for (auto& p : positions) p += {0, 0.075, 0};
  } else if (type == "test-sphere") {
    set_quads(make_sphere(32, 0.075f, 1));
    for (auto& p : positions) p += {0, 0.075, 0};
  } else if (type == "test-sphere-displaced") {
    set_quads(make_sphere(128, 0.075f, 1));
    for (auto& p : positions) p += {0, 0.075, 0};
  } else if (type == "test-disk") {
    set_quads(make_disk(32, 0.075f, 1));
    for (auto& p : positions) p += {0, 0.075, 0};
  } else if (type == "test-uvcylinder") {
    set_quads(make_rounded_uvcylinder(
        {32, 32, 32}, {0.075, 0.075}, {1, 1, 1}, 0.3 * 0.075));
    for (auto& p : positions) p += {0, 0.075, 0};
  } else if (type == "test-floor") {
    set_quads(make_floor({1, 1}, {2, 2}, {20, 20}));
  } else if (type == "test-quad") {
    set_quads(make_rect({1, 1}, {0.075, 0.075}, {1, 1}));
  } else if (type == "test-quady") {
    set_quads(make_recty({1, 1}, {0.075, 0.075}, {1, 1}));
  } else if (type == "test-quad-displaced") {
    set_quads(make_rect({256, 256}, {0.075, 0.075}, {1, 1}));
  } else if (type == "test-quady-displaced") {
    set_quads(make_recty({256, 256}, {0.075, 0.075}, {1, 1}));
  } else if (type == "test-matball") {
    set_quads(make_sphere(32, 0.075));
    for (auto& p : positions) p += {0, 0.075, 0};
  } else if (type == "test-hairball1") {
    auto base = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : base.positions) p += {0, 0.075, 0};
    set_lines(make_hair(base, {4, 65536}, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f}, {0.03, 100}));
  } else if (type == "test-hairball2") {
    auto base = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : base.positions) p += {0, 0.075, 0};
    set_lines(make_hair(base, {4, 65536}, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f}));
  } else if (type == "test-hairball3") {
    auto base = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : base.positions) p += {0, 0.075, 0};
    set_lines(make_hair(base, {4, 65536}, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f}, {0, 0}, {0.5, 128}));
  } else if (type == "test-hairball-interior") {
    set_quads(make_sphere(32, 0.075f * 0.8f, 1));
    for (auto& p : positions) p += {0, 0.075, 0};
  } else if (type == "test-suzanne-subdiv") {
    set_quads(make_monkey(0.075f * 0.8f));
    for (auto& p : positions) p += {0, 0.075, 0};
  } else if (type == "test-cube-subdiv") {
    // set_quads(make_cube( 0.075f);
    set_fvquads(make_fvcube(0.075f));
    // make_fvbox(quadspos, quadsnorm, quadstexcoord, positions, normals,
    //      texcoords, {1, 1, 1}, {0.075f, 0.075f, 0.075f});
    for (auto& p : positions) p += {0, 0.075, 0};
  } else if (type == "test-arealight1") {
    set_quads(make_rect({1, 1}, {0.2, 0.2}));
  } else if (type == "test-arealight2") {
    set_quads(make_rect({1, 1}, {0.2, 0.2}));
  } else if (type == "test-largearealight1") {
    set_quads(make_rect({1, 1}, {0.4, 0.4}));
  } else if (type == "test-largearealight2") {
    set_quads(make_rect({1, 1}, {0.4, 0.4}));
  } else if (type == "test-pointlight1") {
    set_points(make_point(0));
  } else if (type == "test-pointlight2") {
    set_points(make_point(0));
  } else if (type == "test-point") {
    set_points(make_points(1));
  } else if (type == "test-points") {
    set_points(make_points(4096));
  } else if (type == "test-points-random") {
    set_points(make_random_points(4096, {0.2, 0.2, 0.2}));
  } else if (type == "test-particles") {
    set_points(make_points(4096));
  } else if (type == "test-cloth") {
    set_quads(make_rect({64, 64}, {0.2, 0.2}));
  } else if (type == "test-clothy") {
    set_quads(make_recty({64, 64}, {0.2, 0.2}));
  } else {
    error = "unknown preset";
    return false;
  }
  return true;
}

// Shape presets used ofr testing.
bool make_shape_preset(
    generic_shape& shape, const string& type, string& error) {
  return make_shape_preset(shape.points, shape.lines, shape.triangles,
      shape.quads, shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
      shape.positions, shape.normals, shape.texcoords, shape.colors,
      shape.radius, type, error);
}

int main(int argc, const char* argv[]) {
  // command line parameters
  auto facevarying   = false;
  auto positiononly  = false;
  auto trianglesonly = false;
  auto smooth        = false;
  auto faceted       = false;
  auto rotate        = zero3f;
  auto scale         = vec3f{1, 1, 1};
  auto uscale        = 1.0f;
  auto translate     = zero3f;
  auto info          = false;
  auto output        = "out.ply"s;
  auto filename      = "mesh.ply"s;

  // parse command line
  auto cli = make_cli("ymshproc", "Applies operations on a triangle mesh");
  add_option(cli, "--facevarying", facevarying, "Preserve facevarying");
  add_option(cli, "--positiononly", positiononly, "Remove all but positions");
  add_option(cli, "--trianglesonly", trianglesonly, "Remove all but triangles");
  add_option(cli, "--smooth", smooth, "Compute smooth normals");
  add_option(cli, "--faceted", faceted, "Remove normals");
  add_option(cli, "--rotatey,-ry", rotate.y, "Rotate around y axis");
  add_option(cli, "--rotatex,-rx", rotate.x, "Rotate around x axis");
  add_option(cli, "--rotatez,-rz", rotate.z, "Rotate around z axis");
  add_option(cli, "--translatey,-ty", translate.y, "Translate along y axis");
  add_option(cli, "--translatex,-tx", translate.x, "Translate along x axis");
  add_option(cli, "--translatez,-tz", translate.z, "Translate along z axis");
  add_option(cli, "--scale,-s", uscale, "Scale along xyz axes");
  add_option(cli, "--scaley,-sy", scale.y, "Scale along y axis");
  add_option(cli, "--scalex,-sx", scale.x, "Scale along x axis");
  add_option(cli, "--scalez,-sz", scale.z, "Scale along z axis");
  add_option(cli, "--info,-i", info, "print mesh info");
  add_option(cli, "--output,-o", output, "output mesh");
  add_option(cli, "mesh", filename, "input mesh", true);
  parse_cli(cli, argc, argv);

  // mesh data
  auto shape = generic_shape{};

  // load mesh
  auto ioerror = ""s;
  print_progress("load shape", 0, 1);
  auto ext      = sfs::path(filename).extension().string();
  auto basename = sfs::path(filename).stem().string();
  if (ext == ".ypreset") {
    if (!make_shape_preset(shape, basename, ioerror)) print_fatal(ioerror);
  } else {
    if (!load_shape(filename, shape, ioerror, facevarying))
      print_fatal(ioerror);
  }
  print_progress("load shape", 1, 1);

  // remove data
  if (positiononly) {
    shape.normals       = {};
    shape.texcoords     = {};
    shape.colors        = {};
    shape.radius        = {};
    shape.quadsnorm     = {};
    shape.quadstexcoord = {};
    if (!shape.quadspos.empty()) swap(shape.quads, shape.quadspos);
  }

  // convert data
  if (trianglesonly) {
    if (!shape.quadspos.empty())
      throw std::runtime_error("cannot convert facevarying data to triangles");
    if (!shape.quads.empty()) {
      shape.triangles = quads_to_triangles(shape.quads);
      shape.quads     = {};
    }
  }

  // print info
  if (info) {
    print_info("shape stats ------------");
    auto stats = shape_stats(shape);
    for (auto& stat : stats) print_info(stat);
  }

  // transform
  if (uscale != 1) scale *= uscale;
  if (translate != zero3f || rotate != zero3f || scale != vec3f{1, 1, 1}) {
    print_progress("transform shape", 0, 1);
    auto xform = translation_frame(translate) * scaling_frame(scale) *
                 rotation_frame({1, 0, 0}, radians(rotate.x)) *
                 rotation_frame({0, 0, 1}, radians(rotate.z)) *
                 rotation_frame({0, 1, 0}, radians(rotate.y));
    for (auto& p : shape.positions) p = transform_point(xform, p);
    for (auto& n : shape.normals)
      n = transform_normal(xform, n, max(scale) != min(scale));
    print_progress("transform shape", 1, 1);
  }

  // compute normals
  if (smooth) {
    print_progress("smooth shape", 0, 1);
    if (!shape.points.empty()) {
      shape.normals = vector<vec3f>{shape.positions.size(), {0, 0, 1}};
    } else if (!shape.lines.empty()) {
      shape.normals = compute_tangents(shape.lines, shape.positions);
    } else if (!shape.triangles.empty()) {
      shape.normals = compute_normals(shape.triangles, shape.positions);
    } else if (!shape.quads.empty()) {
      shape.normals = compute_normals(shape.quads, shape.positions);
    } else if (!shape.quadspos.empty()) {
      shape.normals = compute_normals(shape.quadspos, shape.positions);
      if (!shape.quadspos.empty()) shape.quadsnorm = shape.quadspos;
    }
    print_progress("smooth shape", 1, 1);
  }

  // remove normals
  if (faceted) {
    print_progress("facet shape", 0, 1);
    shape.normals   = {};
    shape.quadsnorm = {};
    print_progress("facet shape", 1, 1);
  }

  if (info) {
    print_info("shape stats ------------");
    auto stats = shape_stats(shape);
    for (auto& stat : stats) print_info(stat);
  }

  // save mesh
  print_progress("save shape", 0, 1);
  if (!save_shape(output, shape, ioerror, facevarying)) print_fatal(ioerror);
  print_progress("save shape", 1, 1);

  // done
  return 0;
}
