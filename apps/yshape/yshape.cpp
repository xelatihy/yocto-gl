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

void add_options(cli_command& cli, convert_params& params) {
  add_option(cli, "shape", params.shape, "input shape");
  add_option(cli, "output", params.output, "output shape");
  add_option(cli, "smooth", params.smooth, "smooth normals");
  add_option(cli, "facet", params.facet, "facet normals");
  add_option(
      cli, "aspositions", params.aspositions, "remove all but positions");
  add_option(cli, "astriangles", params.astriangles, "convert to triangles");
  add_option(
      cli, "translate", (array<float, 3>&)params.translate, "translate shape");
  add_option(cli, "scale", (array<float, 3>&)params.scale, "scale shape");
  add_option(cli, "rotate", (array<float, 3>&)params.rotate, "rotate shape");
  add_option(cli, "subdivisions", params.subdivisions, "apply subdivision");
  add_option(
      cli, "catmullclark", params.catmullclark, "subdivide as Catmull-Clark");
  add_option(cli, "toedges", params.toedges, "convert shape to edges");
  add_option(cli, "tovertices", params.tovertices, "convert shape to vertices");
}

// convert images
void run_convert(const convert_params& params) {
  // load mesh
  auto shape = load_shape(params.shape, true);

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
    print_info("shape stats ------------");
    auto stats = shape_stats(shape);
    for (auto& stat : stats) print_info(stat);
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
      throw io_error{"empty shape " + params.shape};
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
    print_info("shape stats ------------");
    auto stats = shape_stats(shape);
    for (auto& stat : stats) print_info(stat);
  }

  // save mesh
  save_shape(params.output, shape, true);
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

void add_options(cli_command& cli, fvconvert_params& params) {
  add_option(cli, "shape", params.shape, "input shape");
  add_option(cli, "output", params.output, "output shape");
  add_option(cli, "smooth", params.smooth, "smooth normals");
  add_option(cli, "facet", params.facet, "facet normals");
  add_option(
      cli, "aspositions", params.aspositions, "remove all but positions");
  add_option(
      cli, "translate", (array<float, 3>&)params.translate, "translate shape");
  add_option(cli, "scale", (array<float, 3>&)params.scale, "scale shape");
  add_option(cli, "rotate", (array<float, 3>&)params.rotate, "rotate shape");
  add_option(cli, "subdivisions", params.subdivisions, "apply subdivision");
  add_option(
      cli, "catmullclark", params.catmullclark, "subdivide as Catmull-Clark");
}

// convert images
void run_fvconvert(const fvconvert_params& params) {
  // load mesh
  auto shape = load_fvshape(params.shape, true);

  // remove data
  if (params.aspositions) {
    shape.normals       = {};
    shape.texcoords     = {};
    shape.quadsnorm     = {};
    shape.quadstexcoord = {};
  }

  // print info
  if (params.info) {
    print_info("shape stats ------------");
    auto stats = fvshape_stats(shape);
    for (auto& stat : stats) print_info(stat);
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
    print_info("shape stats ------------");
    auto stats = fvshape_stats(shape);
    for (auto& stat : stats) print_info(stat);
  }

  // save mesh
  save_fvshape(params.output, shape, true);
}

// view params
struct view_params {
  string shape  = "shape.ply";
  string output = "out.ply";
  bool   addsky = false;
};

void add_options(cli_command& cli, view_params& params) {
  add_option(cli, "shape", params.shape, "input shape");
  add_option(cli, "output", params.output, "output shape");
  add_option(cli, "addsky", params.addsky, "add sky");
}

// view shapes
void run_view(const view_params& params) {
  // load shape
  auto shape = load_shape(params.shape, true);

  // make scene
  auto scene = make_shape_scene(shape, params.addsky);

  // run view
  show_trace_gui("yshape", params.shape, scene);
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

void add_options(cli_command& cli, heightfield_params& params) {
  add_option(cli, "image", params.image, "input image");
  add_option(cli, "output", params.output, "output shape");
  add_option(cli, "smooth", params.smooth, "smoooth normals");
  add_option(cli, "height", params.height, "shape height");
  add_option(cli, "info", params.info, "print info");
  add_option(
      cli, "translate", (array<float, 3>&)params.translate, "translate shape");
  add_option(cli, "scale", (array<float, 3>&)params.scale, "scale shape");
  add_option(cli, "rotate", (array<float, 3>&)params.rotate, "rotate shape");
}

void run_heightfield(const heightfield_params& params) {
  // load image
  auto image = load_image(params.image);

  // adjust height
  if (params.height != 1) {
    for (auto& pixel : image.pixels) pixel *= params.height;
  }

  // create heightfield
  auto shape = make_heightfield({image.width, image.height}, image.pixels);
  if (!params.smooth) shape.normals.clear();

  // print info
  if (params.info) {
    print_info("shape stats ------------");
    auto stats = shape_stats(shape);
    for (auto& stat : stats) print_info(stat);
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
  save_shape(params.output, shape, true);
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

void add_options(cli_command& cli, hair_params& params) {
  add_option(cli, "shape", params.shape, "input shape");
  add_option(cli, "output", params.output, "output shape");
  add_option(cli, "hairs", params.hairs, "number of hairs");
  add_option(cli, "steps", params.steps, "hair steps");
  add_option(cli, "length", params.length, "hair length");
  add_option(cli, "noise", params.noise, "noise weight");
  add_option(cli, "gravity", params.gravity, "gravity scale");
  add_option(cli, "radius", params.radius, "hair radius");
}

void run_hair(const hair_params& params) {
  // load mesh
  auto shape = load_shape(params.shape);

  // generate hair
  auto hair = make_hair2(shape, {params.steps, params.hairs},
      {params.length, params.length}, {params.radius, params.radius},
      params.noise, params.gravity);

  // save mesh
  save_shape(params.output, hair, true);
}

struct sample_params {
  string shape   = "shape.ply";
  string output  = "out.ply";
  int    samples = 4096;
};

void add_options(cli_command& cli, sample_params& params) {
  add_option(cli, "shape", params.shape, "input shape");
  add_option(cli, "output", params.output, "output shape");
  add_option(cli, "samples", params.samples, "number of samples");
}

void run_sample(const sample_params& params) {
  // load mesh
  auto shape = load_shape(params.shape);

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
  save_shape(params.output, sshape);
}

struct glview_params {
  string shape  = "shape.ply";
  bool   addsky = false;
};

// Cli
void add_options(cli_command& cli, glview_params& params) {
  add_option(cli, "shape", params.shape, "input shape");
  add_option(cli, "addsky", params.addsky, "add sky");
}

void run_glview(const glview_params& params) {
  // loading shape
  auto shape = load_shape(params.shape);

  // make scene
  auto scene = make_shape_scene(shape, params.addsky);

  // run viewer
  show_shade_gui("yshape", params.shape, scene, {});
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
  try {
    // command line parameters
    auto params = app_params{};
    auto cli    = make_cli("yshape", "Process and view shapes");
    add_command_var(cli, params.command);
    add_command(cli, "convert", params.convert, "convert shapes");
    add_command(
        cli, "fvconvert", params.fvconvert, "convert face-varying shapes");
    add_command(cli, "view", params.view, "view shapes");
    add_command(
        cli, "heightfield", params.heightfield, "create an heightfield");
    add_command(cli, "hair", params.hair, "grow hairs on a shape");
    add_command(cli, "sample", params.sample, "sample points on a shape");
    add_command(cli, "glview", params.glview, "view shapes with OpenGL");
    parse_cli(cli, argc, argv);

    // dispatch commands
    if (params.command == "convert") {
      run_convert(params.convert);
    } else if (params.command == "fvconvert") {
      run_fvconvert(params.fvconvert);
    } else if (params.command == "view") {
      run_view(params.view);
    } else if (params.command == "heightfield") {
      run_heightfield(params.heightfield);
    } else if (params.command == "hair") {
      run_hair(params.hair);
    } else if (params.command == "sample") {
      run_sample(params.sample);
    } else if (params.command == "glview") {
      run_glview(params.glview);
    } else {
      throw io_error{"unknown command"};
    }
  } catch (const std::exception& error) {
    print_error(error.what());
    return 1;
  }

  // done
  return 0;
}
