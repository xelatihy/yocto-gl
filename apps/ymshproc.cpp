//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_shape.h"
using namespace yocto;

#include "ext/CLI11.hpp"

int main(int argc, char** argv) {
  // command line parameters
  auto geodesic_source = -1;
  auto facevarying     = false;
  auto normals         = false;
  auto rotate          = zero3f;
  auto scale           = vec3f{1};
  auto uscale          = 1.0f;
  auto translate       = zero3f;
  auto output          = "out.ply"s;
  auto filename        = "mesh.ply"s;

  // parse command line
  auto parser = CLI::App{"Applies operations on a triangle mesh"};
  parser.add_option("--geodesic-source,-g", geodesic_source, "Geodesic source");
  parser.add_flag("--facevarying", facevarying, "Preserve facevarying");
  parser.add_flag("--normals", normals, "Compute smooth normals");
  parser.add_option("--rotatey", rotate.y, "Rotate around y axis");
  parser.add_option("--rotatex", rotate.x, "Rotate around x axis");
  parser.add_option("--rotatez", rotate.z, "Rotate around z axis");
  parser.add_option("--translatey", translate.y, "Translate along y axis");
  parser.add_option("--translatex", translate.x, "Translate along x axis");
  parser.add_option("--translatez", translate.z, "Translate along z axis");
  parser.add_option("--scale", uscale, "Scale along xyz axes");
  parser.add_option("--scaley", scale.y, "Scale along y axis");
  parser.add_option("--scalex", scale.x, "Scale along x axis");
  parser.add_option("--scalez", scale.z, "Scale along z axis");
  parser.add_option("--output,-o", output, "output mesh")->required(true);
  parser.add_option("mesh", filename, "input mesh")->required(true);
  try {
    parser.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return parser.exit(e);
  }

  // load mesh
  auto shape = yocto_shape{};
  try {
    printf("loading shape");
    auto load_timer = timer();
    load_shape(filename, shape.points, shape.lines, shape.triangles,
        shape.quads, shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
        shape.positions, shape.normals, shape.texcoords, shape.colors,
        shape.radius, facevarying);
    printf(" in %s\n", load_timer.elapsedf().c_str());
  } catch (const std::exception& e) {
    printf("%s\n", e.what());
    exit(1);
  }

  // transform
  if (uscale != 1) scale *= uscale;
  if (translate != zero3f || rotate != zero3f || scale != vec3f{1}) {
    printf("transforming shape");
    auto transform_timer = timer();
    auto xform = translation_frame(translate) * scaling_frame(scale) *
                 rotation_frame({1, 0, 0}, radians(rotate.x)) *
                 rotation_frame({0, 0, 1}, radians(rotate.z)) *
                 rotation_frame({0, 1, 0}, radians(rotate.y));
    for (auto& p : shape.positions) p = transform_point(xform, p);
    for (auto& n : shape.normals)
      n = transform_normal(xform, n, max(scale) != min(scale));
    printf(" in %s\n", transform_timer.elapsedf().c_str());
  }

  // compute normals
  if (normals) {
    printf("computing normals");
    auto transform_timer = timer();
    shape.normals = compute_normals(shape);
    if (!shape.quadspos.empty()) shape.quadsnorm = shape.quadspos;
    printf(" in %s\n", transform_timer.elapsedf().c_str());
  }

  // compute geodesics and store them as colors
  if (geodesic_source >= 0) {
    printf("computing geodesics");
    auto transform_timer = timer();
    auto solver = geodesic_solver{};
    init_geodesic_solver(solver, shape.triangles, shape.positions);
    auto distances = vector<float>{};
    compute_geodesic_distances(solver, distances, {geodesic_source});
    shape.colors = vector<vec4f>{};
    convert_distance_to_color(shape.colors, distances);
    printf(" in %s\n", transform_timer.elapsedf().c_str());
  }

  // save mesh
  try {
    printf("saving shape");
    auto save_timer = timer();
    save_shape(output, shape.points, shape.lines, shape.triangles, shape.quads,
        shape.quadspos, shape.quadsnorm, shape.quadstexcoord, shape.positions,
        shape.normals, shape.texcoords, shape.colors, shape.radius);
    printf(" in %s\n", save_timer.elapsedf().c_str());
  } catch (const std::exception& e) {
    printf("%s\n", e.what());
    exit(1);
  }

  // done
  return 0;
}
