//
// LICENSE:
//
// Copyright (c) 2016 -- 2022 Fabio Pellacini
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

#include <yocto/yocto_cli.h>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_gui.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_scene.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>

using namespace yocto;
using namespace std::string_literals;

// main function
void run(const vector<string>& args) {
  // parameters
  auto shapename    = "shape.ply"s;
  auto outname      = "out.ply"s;
  auto facevarying  = false;
  auto info         = false;
  auto smooth       = false;
  auto facet        = false;
  auto aspositions  = false;
  auto astriangles  = false;
  auto translate    = vec3f{0, 0, 0};
  auto rotate       = vec3f{0, 0, 0};
  auto scale        = vec3f{1, 1, 1};
  auto subdivisions = 0;
  auto catmullclark = false;
  auto toedges      = false;
  auto tovertices   = false;

  // parse command line
  auto cli = make_cli("yconverts", "convert shapes");
  add_option(cli, "shape", shapename, "input shape");
  add_option(cli, "output", outname, "output shape");
  add_option(cli, "facevarying", facevarying, "process as facevarying");
  add_option(cli, "smooth", smooth, "smooth normals");
  add_option(cli, "facet", facet, "facet normals");
  add_option(cli, "aspositions", aspositions, "remove all but positions");
  add_option(cli, "astriangles", astriangles, "convert to triangles");
  add_option(cli, "translate", (array<float, 3>&)translate, "translate shape");
  add_option(cli, "scale", (array<float, 3>&)scale, "scale shape");
  add_option(cli, "rotate", (array<float, 3>&)rotate, "rotate shape");
  add_option(cli, "subdivisions", subdivisions, "apply subdivision");
  add_option(cli, "catmullclark", catmullclark, "subdivide as Catmull-Clark");
  add_option(cli, "toedges", toedges, "convert shape to edges");
  add_option(cli, "tovertices", tovertices, "convert shape to vertices");
  parse_cli(cli, args);

  // start converting

  // switch between facevarying and not
  if (!facevarying) {
    // load mesh
    auto shape = load_shape(shapename, true);

    // remove data
    if (aspositions) {
      shape.normals   = {};
      shape.texcoords = {};
      shape.colors    = {};
      shape.radius    = {};
    }

    // convert data
    if (astriangles) {
      if (!shape.quads.empty()) {
        shape.triangles = quads_to_triangles(shape.quads);
        shape.quads     = {};
      }
    }

    // print stats
    if (info) {
      print_info("shape stats ------------");
      auto stats = shape_stats(shape);
      for (auto& stat : stats) print_info(stat);
    }

    // subdivision
    if (subdivisions > 0) {
      shape = subdivide_shape(shape, subdivisions, catmullclark);
    }

    // transform
    if (translate != vec3f{0, 0, 0} || rotate != vec3f{0, 0, 0} ||
        scale != vec3f{1, 1, 1}) {
      auto translation = translation_frame(translate);
      auto scaling     = scaling_frame(scale);
      auto rotation    = rotation_frame({1, 0, 0}, radians(rotate.x)) *
                      rotation_frame({0, 0, 1}, radians(rotate.z)) *
                      rotation_frame({0, 1, 0}, radians(rotate.y));
      auto xform = translation * scaling * rotation;
      for (auto& p : shape.positions) p = transform_point(xform, p);
      auto nonuniform_scaling = min(scale) != max(scale);
      for (auto& n : shape.normals)
        n = transform_normal(xform, n, nonuniform_scaling);
    }

    // convert to edges
    if (toedges) {
      // check faces
      if (shape.triangles.empty() && shape.quads.empty()) {
        throw io_error{"empty shape " + shapename};
      }

      // convert to edges
      auto edges = !shape.triangles.empty() ? get_edges(shape.triangles)
                                            : get_edges(shape.quads);
      shape      = lines_to_cylinders(edges, shape.positions, 4, 0.001f);
    }

    // convert to vertices
    if (tovertices) {
      // convert to spheres
      shape = points_to_spheres(shape.positions);
    }

    // compute normals
    if (smooth) {
      if (!shape.points.empty()) {
        shape.normals = vector<vec3f>{shape.positions.size(), {0, 0, 1}};
      } else if (!shape.lines.empty()) {
        shape.normals = lines_tangents(shape.lines, shape.positions);
      } else if (!shape.triangles.empty()) {
        shape.normals = triangles_normals(shape.triangles, shape.positions);
      } else if (!shape.quads.empty()) {
        shape.normals = quads_normals(shape.quads, shape.positions);
      }
    }

    // remove normals
    if (facet) {
      shape.normals = {};
    }

    if (info) {
      print_info("shape stats ------------");
      auto stats = shape_stats(shape);
      for (auto& stat : stats) print_info(stat);
    }

    // save mesh
    save_shape(outname, shape, true);
  } else {
    // load mesh
    auto shape = load_fvshape(shapename, true);

    // remove data
    if (aspositions) {
      shape.normals       = {};
      shape.texcoords     = {};
      shape.quadsnorm     = {};
      shape.quadstexcoord = {};
    }

    // print info
    if (info) {
      print_info("shape stats ------------");
      auto stats = fvshape_stats(shape);
      for (auto& stat : stats) print_info(stat);
    }

    // subdivision
    if (subdivisions > 0) {
      shape = subdivide_fvshape(shape, subdivisions, catmullclark);
    }

    // transform
    if (translate != vec3f{0, 0, 0} || rotate != vec3f{0, 0, 0} ||
        scale != vec3f{1, 1, 1}) {
      auto translation = translation_frame(translate);
      auto scaling     = scaling_frame(scale);
      auto rotation    = rotation_frame({1, 0, 0}, radians(rotate.x)) *
                      rotation_frame({0, 0, 1}, radians(rotate.z)) *
                      rotation_frame({0, 1, 0}, radians(rotate.y));
      auto xform = translation * scaling * rotation;
      for (auto& p : shape.positions) p = transform_point(xform, p);
      auto nonuniform_scaling = min(scale) != max(scale);
      for (auto& n : shape.normals)
        n = transform_normal(xform, n, nonuniform_scaling);
    }

    // compute normals
    if (smooth) {
      if (!shape.quadspos.empty()) {
        shape.normals = quads_normals(shape.quadspos, shape.positions);
        if (!shape.quadspos.empty()) shape.quadsnorm = shape.quadspos;
      }
    }

    // remove normals
    if (facet) {
      shape.normals   = {};
      shape.quadsnorm = {};
    }

    if (info) {
      print_info("shape stats ------------");
      auto stats = fvshape_stats(shape);
      for (auto& stat : stats) print_info(stat);
    }

    // save mesh
    save_fvshape(outname, shape, true);
  }
}

// Main
int main(int argc, const char* argv[]) {
  try {
    run({argv, argv + argc});
    return 0;
  } catch (const std::exception& error) {
    print_error(error.what());
    return 1;
  }
}
