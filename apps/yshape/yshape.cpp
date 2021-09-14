//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
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

#include <fmt/core.h>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_gui.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_scene.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>

#include "ext/CLI11.hpp"

using namespace yocto;

// convert params
struct convert_params {
  string shape        = "shape.ply";
  string output       = "out.ply";
  bool   info         = false;
  bool   smooth       = false;
  bool   facet        = false;
  bool   aspositions  = false;
  bool   astriangles  = false;
  vec3f  translate    = {0, 0, 0};
  vec3f  rotate       = {0, 0, 0};
  vec3f  scale        = {1, 1, 1};
  int    subdivisions = 0;
  bool   catmullclark = false;
  bool   toedges      = false;
  bool   tovertices   = false;
};

void add_options(CLI::App& cli, convert_params& params) {
  cli.add_option("shape", params.shape, "Input shape.");
  cli.add_option("--output", params.output, "Output shape.");
  cli.add_flag("--smooth", params.smooth, "Smooth normals.");
  cli.add_flag("--facet", params.facet, "Facet normals.");
  cli.add_flag(
      "--aspositions", params.aspositions, "Remove all but positions.");
  cli.add_flag("--astriangles", params.astriangles, "Convert to triangles.");
  cli.add_option(
      "--translate", (array<float, 3>&)params.translate, "Translate shape.");
  cli.add_option("--scale", (array<float, 3>&)params.scale, "Scale shape.");
  cli.add_option("--rotate", (array<float, 3>&)params.rotate, "Rotate shape.");
  cli.add_option("--subdivisions", params.subdivisions, "Apply subdivision.");
  cli.add_flag(
      "--catmullclark", params.catmullclark, "Catmull-Clark subdivision.");
  cli.add_flag("--toedges", params.toedges, "Convert shape to edges.");
  cli.add_flag("--tovertices", params.tovertices, "Convert shape to vertices.");
}

// convert images
int run_convert(const convert_params& params) {
  fmt::print("converting {}\n", params.shape);

  // load mesh
  auto error = string{};
  auto shape = shape_data{};
  if (!load_shape(params.shape, shape, error, true)) {
    fmt::print("error: cannot load {}\n", params.shape);
    return 1;
  }

  // remove data
  if (params.aspositions) {
    shape.normals   = {};
    shape.texcoords = {};
    shape.colors    = {};
    shape.radius    = {};
  }

  // convert data
  if (params.astriangles) {
    if (!shape.quads.empty()) {
      shape.triangles = quads_to_triangles(shape.quads);
      shape.quads     = {};
    }
  }

  // print stats
  if (params.info) {
    fmt::print("shape stats ------------\n");
    auto stats = shape_stats(shape);
    for (auto& stat : stats) fmt::print("{}\n", stat);
  }

  // subdivision
  if (params.subdivisions > 0) {
    shape = subdivide_shape(shape, params.subdivisions, params.catmullclark);
  }

  // transform
  if (params.translate != vec3f{0, 0, 0} || params.rotate != vec3f{0, 0, 0} ||
      params.scale != vec3f{1, 1, 1}) {
    auto translation = translation_frame(params.translate);
    auto scaling     = scaling_frame(params.scale);
    auto rotation    = rotation_frame({1, 0, 0}, radians(params.rotate.x)) *
                    rotation_frame({0, 0, 1}, radians(params.rotate.z)) *
                    rotation_frame({0, 1, 0}, radians(params.rotate.y));
    auto xform = translation * scaling * rotation;
    for (auto& p : shape.positions) p = transform_point(xform, p);
    auto nonuniform_scaling = min(params.scale) != max(params.scale);
    for (auto& n : shape.normals)
      n = transform_normal(xform, n, nonuniform_scaling);
  }

  // convert to edges
  if (params.toedges) {
    // check faces
    if (shape.triangles.empty() && shape.quads.empty()) {
      fmt::print("error: empty shape {}\n", params.shape);
      return 1;
    }

    // convert to edges
    auto edges = !shape.triangles.empty() ? get_edges(shape.triangles)
                                          : get_edges(shape.quads);
    shape      = lines_to_cylinders(edges, shape.positions, 4, 0.001f);
  }

  // convert to vertices
  if (params.tovertices) {
    // convert to spheres
    shape = points_to_spheres(shape.positions);
  }

  // compute normals
  if (params.smooth) {
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
  if (params.facet) {
    shape.normals = {};
  }

  if (params.info) {
    fmt::print("shape stats ------------\n");
    auto stats = shape_stats(shape);
    for (auto& stat : stats) fmt::print("{}\n", stat);
  }

  // save mesh
  if (!save_shape(params.output, shape, error, true)) {
    fmt::print("error: cannot save {}\n", params.output);
    return 1;
  }

  // done
  return 0;
}

// fvconvert params
struct fvconvert_params {
  string shape        = "shape.obj";
  string output       = "out.obj";
  bool   info         = false;
  bool   smooth       = false;
  bool   facet        = false;
  bool   aspositions  = false;
  vec3f  translate    = {0, 0, 0};
  vec3f  rotate       = {0, 0, 0};
  vec3f  scale        = {1, 1, 1};
  int    subdivisions = 0;
  bool   catmullclark = false;
};

void add_options(CLI::App& cli, fvconvert_params& params) {
  cli.add_option("shape", params.shape, "Input shape.");
  cli.add_option("--output", params.output, "Output shape.");
  cli.add_flag("--smooth", params.smooth, "Smooth normals.");
  cli.add_flag("--facet", params.facet, "Facet normals.");
  cli.add_flag(
      "--aspositions", params.aspositions, "Remove all but positions.");
  cli.add_option(
      "--translate", (array<float, 3>&)params.translate, "Translate shape.");
  cli.add_option("--scale", (array<float, 3>&)params.scale, "Scale shape.");
  cli.add_option("--rotate", (array<float, 3>&)params.rotate, "Rotate shape.");
  cli.add_option("--subdivisions", params.subdivisions, "Apply subdivision.");
  cli.add_flag(
      "--catmullclark", params.catmullclark, "Catmull-Clark subdivision.");
}

// convert images
int run_fvconvert(const fvconvert_params& params) {
  fmt::print("converting {}\n", params.shape);

  // load mesh
  auto error = string{};
  auto shape = fvshape_data{};
  if (!load_fvshape(params.shape, shape, error, true)) {
    fmt::print("error: cannot load {}\n", params.shape);
    return 1;
  }

  // remove data
  if (params.aspositions) {
    shape.normals       = {};
    shape.texcoords     = {};
    shape.quadsnorm     = {};
    shape.quadstexcoord = {};
  }

  // print info
  if (params.info) {
    fmt::print("shape stats ------------\n");
    auto stats = fvshape_stats(shape);
    for (auto& stat : stats) fmt::print("{}\n", stat);
  }

  // subdivision
  if (params.subdivisions > 0) {
    shape = subdivide_fvshape(shape, params.subdivisions, params.catmullclark);
  }

  // transform
  if (params.translate != vec3f{0, 0, 0} || params.rotate != vec3f{0, 0, 0} ||
      params.scale != vec3f{1, 1, 1}) {
    auto translation = translation_frame(params.translate);
    auto scaling     = scaling_frame(params.scale);
    auto rotation    = rotation_frame({1, 0, 0}, radians(params.rotate.x)) *
                    rotation_frame({0, 0, 1}, radians(params.rotate.z)) *
                    rotation_frame({0, 1, 0}, radians(params.rotate.y));
    auto xform = translation * scaling * rotation;
    for (auto& p : shape.positions) p = transform_point(xform, p);
    auto nonuniform_scaling = min(params.scale) != max(params.scale);
    for (auto& n : shape.normals)
      n = transform_normal(xform, n, nonuniform_scaling);
  }

  // compute normals
  if (params.smooth) {
    if (!shape.quadspos.empty()) {
      shape.normals = quads_normals(shape.quadspos, shape.positions);
      if (!shape.quadspos.empty()) shape.quadsnorm = shape.quadspos;
    }
  }

  // remove normals
  if (params.facet) {
    shape.normals   = {};
    shape.quadsnorm = {};
  }

  if (params.info) {
    fmt::print("shape stats ------------\n");
    auto stats = fvshape_stats(shape);
    for (auto& stat : stats) fmt::print("{}\n", stat);
  }

  // save mesh
  if (!save_fvshape(params.output, shape, error, true)) {
    fmt::print("error: cannot save {}\n", params.output);
    return 1;
  }

  // done
  return 0;
}

// view params
struct view_params {
  string shape  = "shape.ply";
  string output = "out.ply";
  bool   addsky = false;
};

void add_options(CLI::App& cli, view_params& params) {
  cli.add_option("shape", params.shape, "Input shape.");
  cli.add_option("--output", params.output, "Output shape.");
  cli.add_flag("--addsky", params.addsky, "Add sky.");
}

// view shapes
int run_view(const view_params& params) {
  fmt::print("viewing {}\n", params.shape);

  // load shape
  auto error = string{};
  auto shape = shape_data{};
  if (!load_shape(params.shape, shape, error, true)) {
    fmt::print("error: cannot load {}\n", params.shape);
    return 1;
  }

  // make scene
  auto scene = make_shape_scene(shape, params.addsky);

  // run view
  show_trace_gui("yshape", params.shape, scene);

  // done
  return 0;
}

struct heightfield_params {
  string image     = "heightfield.png";
  string output    = "out.ply";
  bool   smooth    = false;
  float  height    = 1.0f;
  bool   info      = false;
  vec3f  translate = {0, 0, 0};
  vec3f  rotate    = {0, 0, 0};
  vec3f  scale     = {1, 1, 1};
};

void add_options(CLI::App& cli, heightfield_params& params) {
  cli.add_option("image", params.image, "Input image.");
  cli.add_option("--output", params.output, "Output shape.");
  cli.add_flag("--smooth", params.smooth, "Smoooth normals.");
  cli.add_option("--height", params.height, "Shape height.");
  cli.add_flag("--info", params.info, "Print info.");
  cli.add_option(
      "--translate", (array<float, 3>&)params.translate, "Translate shape.");
  cli.add_option("--scale", (array<float, 3>&)params.scale, "Scale shape.");
  cli.add_option("--rotate", (array<float, 3>&)params.rotate, "Rotate shape.");
}

int run_heightfield(const heightfield_params& params) {
  // load image
  auto error = string{};
  auto image = image_data{};
  if (!load_image(params.image, image, error)) {
    fmt::print("error: cannot load {}\n", params.image);
    return 1;
  }

  // adjust height
  if (params.height != 1) {
    for (auto& pixel : image.pixels) pixel *= params.height;
  }

  // create heightfield
  auto shape = make_heightfield({image.width, image.height}, image.pixels);
  if (!params.smooth) shape.normals.clear();

  // print info
  if (params.info) {
    fmt::print("shape stats ------------\n");
    auto stats = shape_stats(shape);
    for (auto& stat : stats) fmt::print("{}\n", stat);
  }

  // transform
  if (params.translate != vec3f{0, 0, 0} || params.rotate != vec3f{0, 0, 0} ||
      params.scale != vec3f{1, 1, 1}) {
    auto translation = translation_frame(params.translate);
    auto scaling     = scaling_frame(params.scale);
    auto rotation    = rotation_frame({1, 0, 0}, radians(params.rotate.x)) *
                    rotation_frame({0, 0, 1}, radians(params.rotate.z)) *
                    rotation_frame({0, 1, 0}, radians(params.rotate.y));
    auto xform = translation * scaling * rotation;
    for (auto& p : shape.positions) p = transform_point(xform, p);
    auto nonuniform_scaling = min(params.scale) != max(params.scale);
    for (auto& n : shape.normals)
      n = transform_normal(xform, n, nonuniform_scaling);
  }

  // save mesh
  if (!save_shape(params.output, shape, error, true)) {
    fmt::print("error: cannot save {}\n", params.output);
    return 1;
  }

  // done
  return 0;
}

struct hair_params {
  string shape   = "shape.ply";
  string output  = "out.ply";
  int    hairs   = 65536;
  int    steps   = 8;
  float  length  = 0.02f;
  float  noise   = 0.001f;
  float  gravity = 0.0005f;
  float  radius  = 0.0001f;
};

void add_options(CLI::App& cli, hair_params& params) {
  cli.add_option("shape", params.shape, "Input shape.");
  cli.add_option("--output", params.output, "Output shape.");
  cli.add_option("--hairs", params.hairs, "Number of hairs.");
  cli.add_option("--steps", params.steps, "Hair steps.");
  cli.add_option("--length", params.length, "Hair length.");
  cli.add_option("--noise", params.noise, "Noise weight.");
  cli.add_option("--gravity", params.gravity, "Gravity scale.");
  cli.add_option("--radius", params.radius, "Hair radius.");
}

int run_hair(const hair_params& params) {
  // load mesh
  auto error = string{};
  auto shape = shape_data{};
  if (!load_shape(params.shape, shape, error)) {
    fmt::print("error: cannot load {}\n", params.shape);
    return 1;
  }

  // generate hair
  auto hair = make_hair2(shape, {params.steps, params.hairs},
      {params.length, params.length}, {params.radius, params.radius},
      params.noise, params.gravity);

  // save mesh
  if (!save_shape(params.output, hair, error, true)) {
    fmt::print("error: cannot save {}\n", params.output);
    return 1;
  }

  // done
  return 0;
}

struct sample_params {
  string shape   = "shape.ply";
  string output  = "out.ply";
  int    samples = 4096;
};

void add_options(CLI::App& cli, sample_params& params) {
  cli.add_option("shape", params.shape, "Input shape.");
  cli.add_option("--output", params.output, "Output shape.");
  cli.add_option("--samples", params.samples, "Number of samples.");
}

int run_sample(const sample_params& params) {
  // load mesh
  auto error = string{};
  auto shape = shape_data{};
  if (!load_shape(params.shape, shape, error)) {
    fmt::print("error: cannot load {}\n", params.shape);
    return 1;
  }

  // generate samples
  auto samples = sample_shape(shape, params.samples);

  // sample shape
  auto sshape = shape_data{};
  for (auto& sample : samples) {
    sshape.points.push_back((int)sshape.points.size());
    sshape.positions.push_back(eval_position(shape, sample.element, sample.uv));
    sshape.radius.push_back(0.001f);
  }

  // save mesh
  if (!save_shape(params.output, sshape, error)) {
    fmt::print("error: cannot save {}\n", params.output);
    return 1;
  }

  // done
  return 0;
}

struct glview_params {
  string shape  = "shape.ply";
  bool   addsky = false;
};

// Cli
void add_options(CLI::App& cli, glview_params& params) {
  cli.add_option("shape", params.shape, "Input shape.");
  cli.add_flag("--addsky", params.addsky, "Add sky.");
}

int run_glview(const glview_params& params) {
  // loading shape
  auto error = string{};
  auto shape = shape_data{};
  if (!load_shape(params.shape, shape, error)) {
    fmt::print("error: cannot load {}\n", params.shape);
    return 1;
  }

  // make scene
  auto scene = make_shape_scene(shape, params.addsky);

  // run viewer
  show_shade_gui("yshape", params.shape, scene, {});

  // done
  return 0;
}

struct app_params {
  string             command     = "convert";
  convert_params     convert     = {};
  fvconvert_params   fvconvert   = {};
  view_params        view        = {};
  heightfield_params heightfield = {};
  hair_params        hair        = {};
  sample_params      sample      = {};
  glview_params      glview      = {};
};

// Run
int main(int argc, const char** argv) {
  // command line parameters
  auto params = app_params{};
  auto cli    = CLI::App("Process and view shapes");
  add_options(
      *cli.add_subcommand("convert", "Convert shapes."), params.convert);
  add_options(*cli.add_subcommand("fvconvert", "Convert face-varying shapes."),
      params.fvconvert);
  add_options(*cli.add_subcommand("view", "View shapes."), params.view);
  add_options(*cli.add_subcommand("heightfield", "Create an heightfield."),
      params.heightfield);
  add_options(
      *cli.add_subcommand("hair", "Grow hairs on a shape."), params.hair);
  add_options(*cli.add_subcommand("sample", "Sample shapepoints on a shape."),
      params.sample);
  add_options(
      *cli.add_subcommand("glview", "View shapes with OpenGL."), params.glview);
  cli.require_subcommand(1);
  try {
    cli.parse(argc, argv);
    params.command = cli.get_subcommands().front()->get_name();
  } catch (const CLI::ParseError& e) {
    return cli.exit(e);
  }

  // dispatch commands
  if (params.command == "convert") {
    return run_convert(params.convert);
  } else if (params.command == "fvconvert") {
    return run_fvconvert(params.fvconvert);
  } else if (params.command == "view") {
    return run_view(params.view);
  } else if (params.command == "heightfield") {
    return run_heightfield(params.heightfield);
  } else if (params.command == "hair") {
    return run_hair(params.hair);
  } else if (params.command == "sample") {
    return run_sample(params.sample);
  } else if (params.command == "glview") {
    return run_glview(params.glview);
  } else {
    fmt::print("error: unknown command\n");
    return 1;
  }
}
