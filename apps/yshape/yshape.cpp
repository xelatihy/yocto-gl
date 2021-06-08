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

#include <yocto/yocto_cli.h>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_scene.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
#if YOCTO_OPENGL == 1
#include <yocto_gui/yocto_glview.h>
#endif
using namespace yocto;

// convert params
struct convert_params {
  string shape       = "shape.ply";
  string output      = "out.ply";
  bool   info        = false;
  bool   smooth      = false;
  bool   facet       = false;
  bool   aspositions = false;
  bool   astriangles = false;
  vec3f  translate   = {0, 0, 0};
  vec3f  rotate      = {0, 0, 0};
  vec3f  scale       = {1, 1, 1};
  float  scaleu      = 1;
  bool   toedges     = false;
  bool   tovertices  = false;
};

void add_command(const cli_command& cli, const string& name,
    convert_params& params, const string& usage) {
  auto cmd = add_command(cli, name, usage);
  add_argument(cmd, "shape", params.shape, "Input shape.");
  add_option(cmd, "output", params.output, "Output shape.");
  add_option(cmd, "smooth", params.smooth, "Smooth normals.");
  add_option(cmd, "facet", params.facet, "Facet normals.");
  add_option(
      cmd, "aspositions", params.aspositions, "Remove all but positions.");
  add_option(cmd, "astriangles", params.astriangles, "Convert to triangles.");
  add_option(cmd, "translatex", params.translate.x, "Translate shape.");
  add_option(cmd, "translatey", params.translate.y, "Translate shape.");
  add_option(cmd, "translatez", params.translate.z, "Translate shape.");
  add_option(cmd, "scalex", params.scale.x, "Scale shape.");
  add_option(cmd, "scaley", params.scale.y, "Scale shape.");
  add_option(cmd, "scalez", params.scale.z, "Scale shape.");
  add_option(cmd, "scaleu", params.scaleu, "Scale shape.");
  add_option(cmd, "rotatex", params.rotate.x, "Rotate shape.");
  add_option(cmd, "rotatey", params.rotate.y, "Rotate shape.");
  add_option(cmd, "rotatez", params.rotate.z, "Rotate shape.");
  add_option(cmd, "toedges", params.toedges, "Convert shape to edges.");
  add_option(
      cmd, "tovertices", params.tovertices, "Convert shape to vertices.");
}

// convert images
int run_convert(const convert_params& params) {
  // shape data
  auto shape = scene_shape{};

  // load mesh
  auto ioerror = ""s;
  if (!load_shape(params.shape, shape, ioerror, true)) print_fatal(ioerror);

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

  // transform
  if (params.translate != vec3f{0, 0, 0} || params.rotate != vec3f{0, 0, 0} ||
      params.scale != vec3f{1, 1, 1} || params.scaleu != 1) {
    auto translation = translation_frame(params.translate);
    auto scaling     = scaling_frame(params.scale * params.scaleu);
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
    if (shape.triangles.empty() && shape.quads.empty())
      print_fatal("empty faces");

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
  if (!save_shape(params.output, shape, ioerror, true)) print_fatal(ioerror);

  // done
  return 0;
}

// fvconvert params
struct fvconvert_params {
  string shape       = "shape.obj";
  string output      = "out.obj";
  bool   info        = false;
  bool   smooth      = false;
  bool   facet       = false;
  bool   aspositions = false;
  vec3f  translate   = {0, 0, 0};
  vec3f  rotate      = {0, 0, 0};
  vec3f  scale       = {1, 1, 1};
  float  scaleu      = 1;
};

void add_command(const cli_command& cli, const string& name,
    fvconvert_params& params, const string& usage) {
  auto cmd = add_command(cli, name, usage);
  add_argument(cmd, "shape", params.shape, "Input shape.");
  add_option(cmd, "output", params.output, "Output shape.");
  add_option(cmd, "smooth", params.smooth, "Smooth normals.");
  add_option(cmd, "facet", params.facet, "Facet normals.");
  add_option(
      cmd, "aspositions", params.aspositions, "Remove all but positions.");
  add_option(cmd, "translatex", params.translate.x, "Translate shape.");
  add_option(cmd, "translatey", params.translate.y, "Translate shape.");
  add_option(cmd, "translatez", params.translate.z, "Translate shape.");
  add_option(cmd, "scalex", params.scale.x, "Scale shape.");
  add_option(cmd, "scaley", params.scale.y, "Scale shape.");
  add_option(cmd, "scalez", params.scale.z, "Scale shape.");
  add_option(cmd, "scaleu", params.scaleu, "Scale shape.");
  add_option(cmd, "rotatex", params.rotate.x, "Rotate shape.");
  add_option(cmd, "rotatey", params.rotate.y, "Rotate shape.");
  add_option(cmd, "rotatez", params.rotate.z, "Rotate shape.");
}

// convert images
int run_fvconvert(const fvconvert_params& params) {
  // mesh data
  auto shape = scene_fvshape{};

  // load mesh
  auto ioerror = ""s;
  if (path_filename(params.shape) == ".ypreset") {
    if (!make_fvshape_preset(shape, path_basename(params.shape), ioerror))
      print_fatal(ioerror);
  } else {
    if (!load_fvshape(params.shape, shape, ioerror)) print_fatal(ioerror);
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
    print_info("shape stats ------------");
    auto stats = fvshape_stats(shape);
    for (auto& stat : stats) print_info(stat);
  }

  // transform
  if (params.translate != vec3f{0, 0, 0} || params.rotate != vec3f{0, 0, 0} ||
      params.scale != vec3f{1, 1, 1} || params.scaleu != 1) {
    auto translation = translation_frame(params.translate);
    auto scaling     = scaling_frame(params.scale * params.scaleu);
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
  if (!save_fvshape(params.output, shape, ioerror, true)) print_fatal(ioerror);

  // done
  return 0;
}

// view params
struct view_params {
  string shape  = "shape.ply";
  string output = "out.ply";
  bool   addsky = false;
};

void add_command(const cli_command& cli, const string& name,
    view_params& params, const string& usage) {
  auto cmd = add_command(cli, name, usage);
  add_argument(cmd, "shape", params.shape, "Input shape.");
  add_option(cmd, "output", params.output, "Output shape.");
  add_option(cmd, "addsky", params.addsky, "Add sky.");
}

#ifndef YOCTO_OPENGL

// view shapes
int run_view(const view_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

// view shapes
int run_view(const view_params& params) {
  // load shape
  auto shape   = scene_shape{};
  auto ioerror = ""s;
  if (path_filename(params.shape) == ".ypreset") {
    if (!make_shape_preset(shape, path_basename(params.shape), ioerror))
      print_fatal(ioerror);
  } else {
    if (!load_shape(params.shape, shape, ioerror, true)) print_fatal(ioerror);
  }

  // make scene
  auto scene = make_shape_scene(shape, params.addsky);

  // run view
  view_scene("yshape", params.shape, scene);

  // done
  return 0;
}

#endif

struct heightfield_params {
  string image     = "heightfield.png"s;
  string output    = "out.ply"s;
  bool   smooth    = false;
  float  height    = 1.0f;
  bool   info      = false;
  vec3f  translate = {0, 0, 0};
  vec3f  rotate    = {0, 0, 0};
  vec3f  scale     = {1, 1, 1};
  float  scaleu    = 1;
};

void add_command(const cli_command& cli, const string& name,
    heightfield_params& params, const string& usage) {
  auto cmd = add_command(cli, name, usage);
  add_argument(cmd, "image", params.image, "Input image.");
  add_option(cmd, "output", params.output, "Output shape.");
  add_option(cmd, "smooth", params.smooth, "Smoooth normals.");
  add_option(cmd, "height", params.height, "Shape height.");
  add_option(cmd, "info", params.info, "Print info.");
  add_option(cmd, "translatex", params.translate.x, "Translate shape.");
  add_option(cmd, "translatey", params.translate.y, "Translate shape.");
  add_option(cmd, "translatez", params.translate.z, "Translate shape.");
  add_option(cmd, "scalex", params.scale.x, "Scale shape.");
  add_option(cmd, "scaley", params.scale.y, "Scale shape.");
  add_option(cmd, "scalez", params.scale.z, "Scale shape.");
  add_option(cmd, "scaleu", params.scaleu, "Scale shape.");
  add_option(cmd, "rotatex", params.rotate.x, "Rotate shape.");
  add_option(cmd, "rotatey", params.rotate.y, "Rotate shape.");
  add_option(cmd, "rotatez", params.rotate.z, "Rotate shape.");
}

int run_heightfield(const heightfield_params& params) {
  // load mesh
  auto image   = color_image{};
  auto ioerror = ""s;
  if (!load_image(params.image, image, ioerror)) print_fatal(ioerror);

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
  if (!save_shape(params.output, shape, ioerror, true)) print_fatal(ioerror);

  // done
  return 0;
}

struct hair_params {
  string shape   = "shape.ply"s;
  string output  = "out.ply"s;
  int    hairs   = 65536;
  int    steps   = 8;
  float  length  = 0.02f;
  float  noise   = 0.001f;
  float  gravity = 0.0005f;
  float  radius  = 0.0001f;
};

void add_command(const cli_command& cli, const string& name,
    hair_params& params, const string& usage) {
  auto cmd = add_command(cli, name, usage);
  add_argument(cmd, "shape", params.shape, "Input shape.");
  add_option(cmd, "output", params.output, "Output shape.");
  add_option(cmd, "hairs", params.hairs, "Number of hairs.");
  add_option(cmd, "steps", params.steps, "Hair steps.");
  add_option(cmd, "length", params.length, "Hair length.");
  add_option(cmd, "noise", params.noise, "Noise weight.");
  add_option(cmd, "gravity", params.gravity, "Gravity scale.");
  add_option(cmd, "radius", params.radius, "Hair radius.");
}

int run_hair(const hair_params& params) {
  // load mesh
  auto shape   = scene_shape{};
  auto ioerror = ""s;
  if (!load_shape(params.shape, shape, ioerror)) print_fatal(ioerror);

  // generate hair
  auto hair = make_hair2(shape, {params.steps, params.hairs},
      {params.length, params.length}, {params.radius, params.radius},
      params.noise, params.gravity);

  // save mesh
  if (!save_shape(params.output, hair, ioerror, true)) print_fatal(ioerror);

  // done
  return 0;
}

struct sample_params {
  string shape   = "shape.ply"s;
  string output  = "out.ply"s;
  int    samples = 4096;
};

void add_command(const cli_command& cli, const string& name,
    sample_params& params, const string& usage) {
  auto cmd = add_command(cli, name, usage);
  add_argument(cmd, "shape", params.shape, "Input shape.");
  add_option(cmd, "output", params.output, "Output shape.");
  add_option(cmd, "samples", params.samples, "Number of samples.");
}

int run_sample(const sample_params& params) {
  // load mesh
  auto shape   = scene_shape{};
  auto ioerror = ""s;
  if (!load_shape(params.shape, shape, ioerror)) print_fatal(ioerror);

  // generate samples
  auto samples = sample_shape(shape, params.samples);

  // sample shape
  auto sshape = scene_shape{};
  for (auto& sample : samples) {
    sshape.points.push_back((int)sshape.points.size());
    sshape.positions.push_back(eval_position(shape, sample.element, sample.uv));
    sshape.radius.push_back(0.001f);
  }

  // save mesh
  if (!save_shape(params.output, sshape, ioerror, true)) print_fatal(ioerror);

  // done
  return 0;
}

struct glview_params {
  string shape  = "shape.ply";
  bool   addsky = false;
};

// Cli
void add_command(const cli_command& cli, const string& name,
    glview_params& params, const string& usage) {
  auto cmd = add_command(cli, name, usage);
  add_argument(cmd, "shape", params.shape, "Input shape.");
  add_option(cmd, "addsky", params.addsky, "Add sky.");
}

#ifndef YOCTO_OPENGL

// view shapes
int run_glview(const glview_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

int run_glview(const glview_params& params) {
  // loading shape
  auto ioerror = ""s;
  auto shape   = scene_shape{};
  if (!load_shape(params.shape, shape, ioerror, true)) print_fatal(ioerror);

  // make scene
  auto scene = make_shape_scene(shape, params.addsky);

  // run viewer
  glview_scene("yshape", params.shape, scene, {});

  // done
  return 0;
}

#endif

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

// Cli
cli_state make_commands(
    const string& name, app_params& params, const string& usage) {
  auto cli = make_cli(name, usage);
  set_command_var(cli, params.command);
  add_command(cli, "convert", params.convert, "Convert shapes.");
  add_command(
      cli, "fvconvert", params.fvconvert, "Convert face-varying shapes.");
  add_command(cli, "view", params.view, "View shapes.");
  add_command(cli, "heightfield", params.heightfield, "Create an heightfield.");
  add_command(cli, "hair", params.hair, "Grow hairs on a shape.");
  add_command(cli, "sample", params.sample, "Sample shapepoints on a shape.");
  add_command(cli, "glview", params.glview, "View shapes with OpenGL.");
  return cli;
}

// Parse cli
void parse_cli(app_params& params, int argc, const char** argv) {
  auto cli = make_commands("yshape", params, "Process and view shapes.");
  parse_cli(cli, argc, argv);
}

int main(int argc, const char* argv[]) {
  // command line parameters
  auto params = app_params{};
  parse_cli(params, argc, argv);

  // dispatch commands
  if (params.command == "convert") {
    return run_convert(params.convert);
  } else if (params.command == "fvconvert") {
    return run_view(params.view);
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
    return print_fatal("unknown command " + params.command);
  }
}
