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
using namespace yocto::math;
namespace shp = yocto::shape;
namespace cli = yocto::commonio;

#include "ext/filesystem.hpp"
namespace sfs = ghc::filesystem;

using std::string;
using std::vector;
using namespace std::string_literals;

// Shape presets used ofr testing.
bool make_shape_preset(vector<int>& points, std::vector<vec2i>& lines,
    std::vector<vec3i>& triangles, std::vector<vec4i>& quads,
    std::vector<vec4i>& quadspos, std::vector<vec4i>& quadsnorm,
    std::vector<vec4i>& quadstexcoord, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    std::vector<vec3f>& colors, std::vector<float>& radius,
    const std::string& type, std::string& error) {
  if (type == "default-quad") {
    shp::make_recty(quads, positions, normals, texcoords);
  } else if (type == "default-quady") {
    shp::make_rect(quads, positions, normals, texcoords);
  } else if (type == "default-cube") {
    shp::make_box(quads, positions, normals, texcoords);
  } else if (type == "default-cube-rounded") {
    shp::make_rounded_box(quads, positions, normals, texcoords);
  } else if (type == "default-sphere") {
    shp::make_sphere(quads, positions, normals, texcoords);
  } else if (type == "default-disk") {
    shp::make_disk(quads, positions, normals, texcoords);
  } else if (type == "default-disk-bulged") {
    shp::make_bulged_disk(quads, positions, normals, texcoords);
  } else if (type == "default-quad-bulged") {
    shp::make_bulged_rect(quads, positions, normals, texcoords);
  } else if (type == "default-uvsphere") {
    shp::make_uvsphere(quads, positions, normals, texcoords);
  } else if (type == "default-uvsphere-flipcap") {
    shp::make_capped_uvsphere(quads, positions, normals, texcoords);
  } else if (type == "default-uvdisk") {
    shp::make_uvdisk(quads, positions, normals, texcoords);
  } else if (type == "default-uvcylinder") {
    shp::make_uvcylinder(quads, positions, normals, texcoords);
  } else if (type == "default-uvcylinder-rounded") {
    shp::make_rounded_uvcylinder(
        quads, positions, normals, texcoords, {32, 32, 32});
  } else if (type == "default-geosphere") {
    shp::make_geosphere(triangles, positions);
  } else if (type == "default-floor") {
    shp::make_floor(quads, positions, normals, texcoords);
  } else if (type == "default-floor-bent") {
    shp::make_bent_floor(quads, positions, normals, texcoords);
  } else if (type == "default-matball") {
    shp::make_sphere(quads, positions, normals, texcoords);
  } else if (type == "default-hairball") {
    auto base_triangles = std::vector<vec3i>{};
    auto base_quads     = std::vector<vec4i>{};
    auto base_positions = std::vector<vec3f>{};
    auto base_normals   = std::vector<vec3f>{};
    auto base_texcoords = std::vector<vec2f>{};
    shp::make_sphere(
        base_quads, base_positions, base_normals, base_texcoords, pow2(5), 0.8);
    shp::make_hair(lines, positions, normals, texcoords, radius, base_triangles,
        base_quads, base_positions, base_normals, base_texcoords, {4, 65536},
        {0.2, 0.2}, {0.002, 0.001});
  } else if (type == "default-hairball-interior") {
    shp::make_sphere(quads, positions, normals, texcoords, pow2(5), 0.8);
  } else if (type == "default-suzanne") {
    shp::make_monkey(quads, positions);
  } else if (type == "default-cube-facevarying") {
    shp::make_fvbox(
        quadspos, quadsnorm, quadstexcoord, positions, normals, texcoords);
  } else if (type == "default-sphere-facevarying") {
    shp::make_fvsphere(
        quadspos, quadsnorm, quadstexcoord, positions, normals, texcoords);
  } else if (type == "default-quady-displaced") {
    shp::make_recty(quads, positions, normals, texcoords, {256, 256});
  } else if (type == "default-sphere-displaced") {
    shp::make_sphere(quads, positions, normals, texcoords, 128);
  } else if (type == "test-cube") {
    shp::make_rounded_box(quads, positions, normals, texcoords, {32, 32, 32},
        {0.075f, 0.075f, 0.075f}, {1, 1, 1}, 0.3 * 0.075f);
    for (auto& p : positions) p += {0, 0.075, 0};
  } else if (type == "test-uvsphere") {
    shp::make_uvsphere(quads, positions, normals, texcoords, {32, 32}, 0.075);
    for (auto& p : positions) p += {0, 0.075, 0};
  } else if (type == "test-uvsphere-flipcap") {
    shp::make_capped_uvsphere(quads, positions, normals, texcoords, {32, 32},
        0.075, {1, 1}, 0.3 * 0.075);
    for (auto& p : positions) p += {0, 0.075, 0};
  } else if (type == "test-sphere") {
    shp::make_sphere(quads, positions, normals, texcoords, 32, 0.075f, 1);
    for (auto& p : positions) p += {0, 0.075, 0};
  } else if (type == "test-sphere-displaced") {
    shp::make_sphere(quads, positions, normals, texcoords, 128, 0.075f, 1);
    for (auto& p : positions) p += {0, 0.075, 0};
  } else if (type == "test-disk") {
    shp::make_disk(quads, positions, normals, texcoords, 32, 0.075f, 1);
    for (auto& p : positions) p += {0, 0.075, 0};
  } else if (type == "test-uvcylinder") {
    shp::make_rounded_uvcylinder(quads, positions, normals, texcoords,
        {32, 32, 32}, {0.075, 0.075}, {1, 1, 1}, 0.3 * 0.075);
    for (auto& p : positions) p += {0, 0.075, 0};
  } else if (type == "test-floor") {
    shp::make_floor(
        quads, positions, normals, texcoords, {1, 1}, {2, 2}, {20, 20});
  } else if (type == "test-quad") {
    shp::make_rect(
        quads, positions, normals, texcoords, {1, 1}, {0.075, 0.075}, {1, 1});
  } else if (type == "test-quady") {
    shp::make_recty(
        quads, positions, normals, texcoords, {1, 1}, {0.075, 0.075}, {1, 1});
  } else if (type == "test-quad-displaced") {
    shp::make_rect(quads, positions, normals, texcoords, {256, 256},
        {0.075, 0.075}, {1, 1});
  } else if (type == "test-quady-displaced") {
    shp::make_recty(quads, positions, normals, texcoords, {256, 256},
        {0.075, 0.075}, {1, 1});
  } else if (type == "test-matball") {
    shp::make_sphere(quads, positions, normals, texcoords, 32, 0.075);
    for (auto& p : positions) p += {0, 0.075, 0};
  } else if (type == "test-hairball1") {
    auto base_triangles = std::vector<vec3i>{};
    auto base_quads     = std::vector<vec4i>{};
    auto base_positions = std::vector<vec3f>{};
    auto base_normals   = std::vector<vec3f>{};
    auto base_texcoords = std::vector<vec2f>{};
    shp::make_sphere(base_quads, base_positions, base_normals, base_texcoords,
        32, 0.075f * 0.8f, 1);
    for (auto& p : base_positions) p += {0, 0.075, 0};
    shp::make_hair(lines, positions, normals, texcoords, radius, base_triangles,
        base_quads, base_positions, base_normals, base_texcoords, {4, 65536},
        {0.1f * 0.15f, 0.1f * 0.15f}, {0.001f * 0.15f, 0.0005f * 0.15f},
        {0.03, 100});
  } else if (type == "test-hairball2") {
    auto base_triangles = std::vector<vec3i>{};
    auto base_quads     = std::vector<vec4i>{};
    auto base_positions = std::vector<vec3f>{};
    auto base_normals   = std::vector<vec3f>{};
    auto base_texcoords = std::vector<vec2f>{};
    shp::make_sphere(base_quads, base_positions, base_normals, base_texcoords,
        32, 0.075f * 0.8f, 1);
    for (auto& p : base_positions) p += {0, 0.075, 0};
    shp::make_hair(lines, positions, normals, texcoords, radius, base_triangles,
        base_quads, base_positions, base_normals, base_texcoords, {4, 65536},
        {0.1f * 0.15f, 0.1f * 0.15f}, {0.001f * 0.15f, 0.0005f * 0.15f});
  } else if (type == "test-hairball3") {
    auto base_triangles = std::vector<vec3i>{};
    auto base_quads     = std::vector<vec4i>{};
    auto base_positions = std::vector<vec3f>{};
    auto base_normals   = std::vector<vec3f>{};
    auto base_texcoords = std::vector<vec2f>{};
    shp::make_sphere(base_quads, base_positions, base_normals, base_texcoords,
        32, 0.075f * 0.8f, 1);
    for (auto& p : base_positions) p += {0, 0.075, 0};
    shp::make_hair(lines, positions, normals, texcoords, radius, base_triangles,
        base_quads, base_positions, base_normals, base_texcoords, {4, 65536},
        {0.1f * 0.15f, 0.1f * 0.15f}, {0.001f * 0.15f, 0.0005f * 0.15f}, {0, 0},
        {0.5, 128});
  } else if (type == "test-hairball-interior") {
    shp::make_sphere(
        quads, positions, normals, texcoords, 32, 0.075f * 0.8f, 1);
    for (auto& p : positions) p += {0, 0.075, 0};
  } else if (type == "test-suzanne-subdiv") {
    shp::make_monkey(quads, positions, 0.075f * 0.8f);
    for (auto& p : positions) p += {0, 0.075, 0};
  } else if (type == "test-cube-subdiv") {
    // make_cube(quads, positions, normals, texcoords, 0.075f);
    shp::make_fvcube(quadspos, quadsnorm, quadstexcoord, positions, normals,
        texcoords, 0.075f);
    // make_fvbox(quadspos, quadsnorm, quadstexcoord, positions, normals,
    //      texcoords, {1, 1, 1}, {0.075f, 0.075f, 0.075f});
    for (auto& p : positions) p += {0, 0.075, 0};
  } else if (type == "test-arealight1") {
    shp::make_rect(quads, positions, normals, texcoords, {1, 1}, {0.2, 0.2});
  } else if (type == "test-arealight2") {
    shp::make_rect(quads, positions, normals, texcoords, {1, 1}, {0.2, 0.2});
  } else if (type == "test-largearealight1") {
    shp::make_rect(quads, positions, normals, texcoords, {1, 1}, {0.4, 0.4});
  } else if (type == "test-largearealight2") {
    shp::make_rect(quads, positions, normals, texcoords, {1, 1}, {0.4, 0.4});
  } else if (type == "test-pointlight1") {
    shp::make_point(points, positions, normals, texcoords, radius, 0);
  } else if (type == "test-pointlight2") {
    shp::make_point(points, positions, normals, texcoords, radius, 0);
  } else if (type == "test-points") {
    shp::make_points(points, positions, normals, texcoords, radius, 4096);
  } else if (type == "test-points-random") {
    shp::make_random_points(
        points, positions, normals, texcoords, radius, 4096, {0.2, 0.2, 0.2});
  } else if (type == "test-particles") {
    shp::make_points(points, positions, normals, texcoords, radius, 4096);
  } else if (type == "test-cloth") {
    shp::make_rect(quads, positions, normals, texcoords, {32, 32}, {0.2, 0.2});
  } else if (type == "test-clothy") {
    shp::make_recty(quads, positions, normals, texcoords, {32, 32}, {0.2, 0.2});
  } else {
    error = "unknown preset";
    return false;
  }
  return true;
}

// Shape presets used ofr testing.
bool make_shape_preset(vector<int>& points, std::vector<vec2i>& lines,
    std::vector<vec3i>& triangles, std::vector<vec4i>& quads,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, std::vector<vec3f>& colors,
    std::vector<float>& radius, const std::string& type, std::string& error) {
  auto quadspos      = std::vector<vec4i>{};
  auto quadsnorm     = std::vector<vec4i>{};
  auto quadstexcoord = std::vector<vec4i>{};
  if (!make_shape_preset(points, lines, triangles, quads, quadspos, quadsnorm,
          quadstexcoord, positions, normals, texcoords, colors, radius, type,
          error))
    return false;
  if (!quadspos.empty()) throw std::runtime_error("bad preset type");
  return true;
}

// Shape presets used ofr testing.
bool make_shape_preset(vector<vec4i>& quadspos, std::vector<vec4i>& quadsnorm,
    std::vector<vec4i>& quadstexcoord, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    const std::string& type, std::string& error) {
  auto points    = std::vector<int>{};
  auto lines     = std::vector<vec2i>{};
  auto triangles = std::vector<vec3i>{};
  auto quads     = std::vector<vec4i>{};
  auto colors    = std::vector<vec3f>{};
  auto radius    = std::vector<float>{};
  if (!make_shape_preset(points, lines, triangles, quads, quadspos, quadsnorm,
          quadstexcoord, positions, normals, texcoords, colors, radius, type,
          error))
    return false;
  if (quadspos.empty()) throw std::runtime_error("bad preset type");
  return true;
}

int main(int argc, const char* argv[]) {
  // command line parameters
  auto facevarying          = false;
  auto positiononly         = false;
  auto trianglesonly        = false;
  auto smooth               = false;
  auto faceted              = false;
  auto rotate               = zero3f;
  auto scale                = vec3f{1};
  auto uscale               = 1.0f;
  auto translate            = zero3f;
  auto info                 = false;
  auto geodesic_source      = -1;
  int  p0                   = -1;
  int  p1                   = -1;
  int  p2                   = -1;
  auto num_geodesic_samples = 0;
  auto geodesic_scale       = 30.0f;
  auto slice                = false;
  auto output               = "out.ply"s;
  auto filename             = "mesh.ply"s;

  // parse command line
  auto cli = cli::make_cli("ymshproc", "Applies operations on a triangle mesh");
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
  add_option(cli, "--geodesic-source,-g", geodesic_source, "Geodesic source");
  add_option(cli, "--path-vertex0,-p0", p0, "Path vertex 0");
  add_option(cli, "--path-vertex1,-p1", p1, "Path vertex 1");
  add_option(cli, "--path-vertex2,-p2", p2, "Path vertex 2");
  add_option(cli, "--num-geodesic-samples", num_geodesic_samples,
      "Number of sampled geodesic sources");
  add_option(cli, "--geodesic-scale", geodesic_scale, "Geodesic scale");
  add_option(cli, "--slice", slice, "Slice mesh along field isolines");
  add_option(cli, "--output,-o", output, "output mesh");
  add_option(cli, "mesh", filename, "input mesh", true);
  parse_cli(cli, argc, argv);

  // mesh data
  auto positions     = std::vector<vec3f>{};
  auto normals       = std::vector<vec3f>{};
  auto texcoords     = std::vector<vec2f>{};
  auto colors        = std::vector<vec3f>{};
  auto radius        = std::vector<float>{};
  auto points        = std::vector<int>{};
  auto lines         = std::vector<vec2i>{};
  auto triangles     = std::vector<vec3i>{};
  auto quads         = std::vector<vec4i>{};
  auto quadspos      = std::vector<vec4i>{};
  auto quadsnorm     = std::vector<vec4i>{};
  auto quadstexcoord = std::vector<vec4i>{};

  // load mesh
  auto ioerror = ""s;
  cli::print_progress("load shape", 0, 1);
  if (!facevarying) {
    auto ext      = sfs::path(filename).extension().string();
    auto basename = sfs::path(filename).stem().string();
    if (ext == ".ypreset") {
      if (!make_shape_preset(points, lines, triangles, quads, positions,
              normals, texcoords, colors, radius, basename, ioerror))
        cli::print_fatal(ioerror);
    } else {
      if (!shp::load_shape(filename, points, lines, triangles, quads, positions,
              normals, texcoords, colors, radius, ioerror))
        cli::print_fatal(ioerror);
    }
  } else {
    auto ext      = sfs::path(filename).extension().string();
    auto basename = sfs::path(filename).stem().string();
    if (ext == ".ypreset") {
      if (!make_shape_preset(quadspos, quadsnorm, quadstexcoord, positions,
              normals, texcoords, basename, ioerror))
        cli::print_fatal(ioerror);
    } else {
      if (!shp::load_fvshape(filename, quadspos, quadsnorm, quadstexcoord,
              positions, normals, texcoords, ioerror))
        cli::print_fatal(ioerror);
    }
  }
  cli::print_progress("load shape", 1, 1);

  // remove data
  if (positiononly) {
    normals       = {};
    texcoords     = {};
    colors        = {};
    radius        = {};
    quadsnorm     = {};
    quadstexcoord = {};
    if (!quadspos.empty()) swap(quads, quadspos);
  }

  // convert data
  if (trianglesonly) {
    if (!quadspos.empty())
      throw std::runtime_error("cannot convert facevarying data to triangles");
    if (!quads.empty()) {
      triangles = shp::quads_to_triangles(quads);
      quads     = {};
    }
  }

  // print info
  if (info) {
    cli::print_info("shape stats ------------");
    auto stats = shp::shape_stats(points, lines, triangles, quads, quadspos,
        quadsnorm, quadstexcoord, positions, normals, texcoords, colors,
        radius);
    for (auto& stat : stats) cli::print_info(stat);
  }

  // transform
  if (uscale != 1) scale *= uscale;
  if (translate != zero3f || rotate != zero3f || scale != vec3f{1}) {
    cli::print_progress("transform shape", 0, 1);
    auto xform = translation_frame(translate) * scaling_frame(scale) *
                 rotation_frame({1, 0, 0}, radians(rotate.x)) *
                 rotation_frame({0, 0, 1}, radians(rotate.z)) *
                 rotation_frame({0, 1, 0}, radians(rotate.y));
    for (auto& p : positions) p = transform_point(xform, p);
    for (auto& n : normals)
      n = transform_normal(xform, n, max(scale) != min(scale));
    cli::print_progress("transform shape", 1, 1);
  }

  // compute normals
  if (smooth) {
    cli::print_progress("smooth shape", 0, 1);
    if (!points.empty()) {
      normals = std::vector<vec3f>{positions.size(), {0, 0, 1}};
    } else if (!lines.empty()) {
      normals = shp::compute_tangents(lines, positions);
    } else if (!triangles.empty()) {
      normals = shp::compute_normals(triangles, positions);
    } else if (!quads.empty()) {
      normals = shp::compute_normals(quads, positions);
    } else if (!quadspos.empty()) {
      normals = shp::compute_normals(quadspos, positions);
      if (!quadspos.empty()) quadsnorm = quadspos;
    }
    cli::print_progress("smooth shape", 1, 1);
  }

  // remove normals
  if (faceted) {
    cli::print_progress("facet shape", 0, 1);
    normals   = {};
    quadsnorm = {};
    cli::print_progress("facet shape", 1, 1);
  }

  // compute geodesics and store them as colors
  if (geodesic_source >= 0 || num_geodesic_samples > 0) {
    cli::print_progress("compute geodesic", 0, 1);
    auto adjacencies = shp::face_adjacencies(triangles);
    auto solver  = shp::make_geodesic_solver(triangles, adjacencies, positions);
    auto sources = std::vector<int>();
    if (geodesic_source >= 0) {
      sources = {geodesic_source};
    } else {
      sources = shp::sample_vertices_poisson(solver, num_geodesic_samples);
    }
    auto field = shp::compute_geodesic_distances(solver, sources);

    if (slice) {
      auto tags = std::vector<int>(triangles.size(), 0);
      shp::meandering_triangles(
          field, geodesic_scale, 0, 1, 2, triangles, tags, positions, normals);
      for (int i = 0; i < triangles.size(); i++) {
        if (tags[i] == 1) triangles[i] = {-1, -1, -1};
      }
    } else {
      colors = std::vector<vec3f>(positions.size());
      for (int i = 0; i < colors.size(); ++i) {
        colors[i] = vec3f(sinf(geodesic_scale * field[i]));
      }
      // distance_to_color(shape.colors, field, geodesic_scale);
    }
    cli::print_progress("compute geodesic", 1, 1);
  }

  if (p0 != -1) {
    cli::print_progress("cut mesh", 0, 1);
    auto tags        = std::vector<int>(triangles.size(), 0);
    auto adjacencies = shp::face_adjacencies(triangles);
    auto solver = shp::make_geodesic_solver(triangles, adjacencies, positions);

    auto               paths = std::vector<shp::surface_path>();
    std::vector<float> fields[3];
    fields[0] = shp::compute_geodesic_distances(solver, {p0});
    fields[1] = shp::compute_geodesic_distances(solver, {p1});
    fields[2] = shp::compute_geodesic_distances(solver, {p2});
    for (int i = 0; i < 3; ++i) {
      for (auto& f : fields[i]) f = -f;
    }

    paths.push_back(shp::integrate_field(
        triangles, positions, adjacencies, tags, 0, fields[1], p0, p1));

    paths.push_back(shp::integrate_field(
        triangles, positions, adjacencies, tags, 0, fields[2], p1, p2));

    paths.push_back(shp::integrate_field(
        triangles, positions, adjacencies, tags, 0, fields[0], p2, p0));

    auto plines     = std::vector<vec2i>{};
    auto ppositions = std::vector<vec3f>{};
    for (int i = 0; i < 3; i++) {
      auto pos  = make_positions_from_path(paths[i], positions);
      auto line = std::vector<vec2i>(pos.size() - 1);
      for (int k = 0; k < line.size(); k++) {
        line[k] = {k, k + 1};
        line[k] += (int)lines.size();
      }
      plines.insert(plines.end(), line.begin(), line.end());
      ppositions.insert(ppositions.end(), pos.begin(), pos.end());
    }
    points    = {};
    lines     = plines;
    triangles = {};
    quads     = {};
    positions = ppositions;
    normals   = {};
    texcoords = {};
    colors    = {};
    radius    = {};
    cli::print_progress("cut mesh", 1, 1);
  }

  if (info) {
    cli::print_info("shape stats ------------");
    auto stats = shp::shape_stats(points, lines, triangles, quads, quadspos,
        quadsnorm, quadstexcoord, positions, normals, texcoords, colors,
        radius);
    for (auto& stat : stats) cli::print_info(stat);
  }

  // save mesh
  cli::print_progress("save shape", 0, 1);
  if (!quadspos.empty()) {
    if (!shp::save_fvshape(output, quadspos, quadsnorm, quadstexcoord,
            positions, normals, texcoords, ioerror))
      cli::print_fatal(ioerror);
  } else {
    if (!shp::save_shape(output, points, lines, triangles, quads, positions,
            normals, texcoords, colors, radius, ioerror))
      cli::print_fatal(ioerror);
  }
  cli::print_progress("save shape", 1, 1);

  // done
  return 0;
}
