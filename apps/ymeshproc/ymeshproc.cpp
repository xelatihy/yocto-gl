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
#include <yocto/yocto_mesh.h>
#include <yocto/yocto_shape.h>
using namespace yocto;

// Shape presets used ofr testing.
bool make_mesh_preset(vector<vec3i>& triangles, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, vector<vec3f>& colors,
    const string& type, string& error) {
  auto set_quads = [&](quads_shape&& shape) {
    triangles = quads_to_triangles(shape.quads);
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
  auto set_fvquads = [&](quads_fvshape&& shape) {
    triangles = quads_to_triangles(shape.quadspos);
    positions = shape.positions;
    normals   = (shape.quadspos == shape.quadsnorm) ? shape.normals
                                                  : vector<vec3f>{};
    texcoords = (shape.quadspos == shape.quadstexcoord) ? shape.texcoords
                                                        : vector<vec2f>{};
  };

  if (type == "mesh-quad") {
    set_quads(make_rect());
  } else if (type == "mesh-quady") {
    set_quads(make_recty());
  } else if (type == "mesh-cube") {
    set_quads(make_box());
  } else if (type == "mesh-cube-rounded") {
    set_quads(make_rounded_box());
  } else if (type == "mesh-sphere") {
    set_quads(make_sphere());
  } else if (type == "mesh-disk") {
    set_quads(make_disk());
  } else if (type == "mesh-disk-bulged") {
    set_quads(make_bulged_disk());
  } else if (type == "mesh-quad-bulged") {
    set_quads(make_bulged_rect());
  } else if (type == "mesh-uvsphere") {
    set_quads(make_uvsphere());
  } else if (type == "mesh-uvsphere-flipcap") {
    set_quads(make_capped_uvsphere());
  } else if (type == "mesh-uvdisk") {
    set_quads(make_uvdisk());
  } else if (type == "mesh-uvcylinder") {
    set_quads(make_uvcylinder());
  } else if (type == "mesh-uvcylinder-rounded") {
    set_quads(make_rounded_uvcylinder({32, 32, 32}));
  } else if (type == "mesh-geosphere") {
    set_triangles(make_geosphere());
  } else if (type == "mesh-floor") {
    set_quads(make_floor());
  } else if (type == "mesh-floor-bent") {
    set_quads(make_bent_floor());
  } else if (type == "mesh-matball") {
    set_quads(make_sphere());
  } else if (type == "mesh-hairball-interior") {
    set_quads(make_sphere(pow2(5), 0.8));
  } else if (type == "mesh-suzanne") {
    set_quads(make_monkey());
  } else if (type == "mesh-cube-facevarying") {
    set_fvquads(make_fvbox());
  } else if (type == "mesh-sphere-facevarying") {
    set_fvquads(make_fvsphere());
  } else {
    error = "unknown preset";
    return false;
  }

  return true;
}

int main(int argc, const char* argv[]) {
  // command line parameters
  auto smooth               = false;
  auto faceted              = false;
  auto rotate               = zero3f;
  auto scale                = vec3f{1, 1, 1};
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
  auto cli = make_cli("ymshproc", "Applies operations on a triangle mesh");
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
  auto positions = vector<vec3f>{};
  auto normals   = vector<vec3f>{};
  auto texcoords = vector<vec2f>{};
  auto colors    = vector<vec3f>{};
  auto triangles = vector<vec3i>{};
  auto lines     = vector<vec2i>{};  // for line output

  // load mesh
  auto ioerror = ""s;
  print_progress("load mesh", 0, 1);
  if (path{filename}.extension().string() == ".ypreset") {
    if (!make_mesh_preset(triangles, positions, normals, texcoords, colors,
            path{filename}.stem().string(), ioerror))
      print_fatal(ioerror);
  } else {
    if (!load_mesh(filename, triangles, positions, normals, texcoords, colors,
            ioerror))
      print_fatal(ioerror);
  }
  print_progress("load mesh", 1, 1);

  // print info
  if (info) {
    print_info("shape stats ------------");
    auto stats = mesh_stats(triangles, positions, normals, texcoords, colors);
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
    for (auto& p : positions) p = transform_point(xform, p);
    for (auto& n : normals)
      n = transform_normal(xform, n, max(scale) != min(scale));
    print_progress("transform shape", 1, 1);
  }

  // compute normals
  if (smooth) {
    print_progress("smooth shape", 0, 1);
    normals = compute_normals(triangles, positions);
    print_progress("smooth shape", 1, 1);
  }

  // remove normals
  if (faceted) {
    print_progress("facet shape", 0, 1);
    normals = {};
    print_progress("facet shape", 1, 1);
  }

  // compute geodesics and store them as colors
  if (geodesic_source >= 0 || num_geodesic_samples > 0) {
    print_progress("compute geodesic", 0, 1);
    auto adjacencies = face_adjacencies(triangles);
    auto solver      = make_geodesic_solver(triangles, adjacencies, positions);
    auto sources     = vector<int>();
    if (geodesic_source >= 0) {
      sources = {geodesic_source};
    } else {
      sources = sample_vertices_poisson(solver, num_geodesic_samples);
    }
    auto field = compute_geodesic_distances(solver, sources);

    if (slice) {
      auto tags = vector<int>(triangles.size(), 0);
      meandering_triangles(
          field, geodesic_scale, 0, 1, 2, triangles, tags, positions, normals);
      for (int i = 0; i < triangles.size(); i++) {
        if (tags[i] == 1) triangles[i] = {-1, -1, -1};
      }
    } else {
      colors = vector<vec3f>(positions.size());
      for (int i = 0; i < colors.size(); ++i) {
        auto c    = sinf(geodesic_scale * field[i]);
        colors[i] = vec3f{c, c, c};
      }
      // distance_to_color(shape.colors, field, geodesic_scale);
    }
    print_progress("compute geodesic", 1, 1);
  }

  if (p0 != -1) {
    print_progress("cut mesh", 0, 1);
    auto tags        = vector<int>(triangles.size(), 0);
    auto adjacencies = face_adjacencies(triangles);
    auto solver      = make_geodesic_solver(triangles, adjacencies, positions);

    auto          paths = vector<surface_path>();
    vector<float> fields[3];
    fields[0] = compute_geodesic_distances(solver, {p0});
    fields[1] = compute_geodesic_distances(solver, {p1});
    fields[2] = compute_geodesic_distances(solver, {p2});
    for (int i = 0; i < 3; ++i) {
      for (auto& f : fields[i]) f = -f;
    }

    paths.push_back(integrate_field(
        triangles, positions, adjacencies, tags, 0, fields[1], p0, p1));

    paths.push_back(integrate_field(
        triangles, positions, adjacencies, tags, 0, fields[2], p1, p2));

    paths.push_back(integrate_field(
        triangles, positions, adjacencies, tags, 0, fields[0], p2, p0));

    auto plines     = vector<vec2i>{};
    auto ppositions = vector<vec3f>{};
    for (int i = 0; i < 3; i++) {
      auto pos  = make_positions_from_path(paths[i], positions);
      auto line = vector<vec2i>(pos.size() - 1);
      for (int k = 0; k < line.size(); k++) {
        line[k] = {k, k + 1};
        line[k] += (int)lines.size();
      }
      plines.insert(plines.end(), line.begin(), line.end());
      ppositions.insert(ppositions.end(), pos.begin(), pos.end());
    }
    lines     = lines;
    triangles = {};
    positions = positions;
    normals   = {};
    texcoords = {};
    colors    = {};
    print_progress("cut mesh", 1, 1);
  }

  if (info) {
    print_info("mesh stats ------------");
    auto stats = mesh_stats(triangles, positions, normals, texcoords, colors);
    for (auto& stat : stats) print_info(stat);
  }

  // save mesh
  if (!triangles.empty()) {
    print_progress("save mesh", 0, 1);
    if (!save_mesh(
            output, triangles, positions, normals, texcoords, colors, ioerror))
      print_fatal(ioerror);
    print_progress("save mesh", 1, 1);
  }
  if (!lines.empty()) {
    print_progress("save lines", 0, 1);
    if (!save_lines(
            output, lines, positions, normals, texcoords, colors, ioerror))
      print_fatal(ioerror);
    print_progress("save lines", 1, 1);
  }

  // done
  return 0;
}
