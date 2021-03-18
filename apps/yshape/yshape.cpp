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
};

void add_command(cli_command& cli, const string& name, convert_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_argument(cmd, "shape", value.shape, "Input shape.");
  add_option(cmd, "output", value.output, "Output shape.");
  add_option(cmd, "smooth", value.smooth, "Smooth normals.");
  add_option(cmd, "facet", value.facet, "Facet normals.");
  add_option(
      cmd, "aspositions", value.aspositions, "Remove all but positions.");
  add_option(cmd, "astriangles", value.astriangles, "Convert to triangles.");
  add_option(cmd, "translatex", value.translate.x, "Translate shape.");
  add_option(cmd, "translatey", value.translate.y, "Translate shape.");
  add_option(cmd, "translatez", value.translate.z, "Translate shape.");
  add_option(cmd, "scalex", value.scale.x, "Scale shape.");
  add_option(cmd, "scaley", value.scale.y, "Scale shape.");
  add_option(cmd, "scalez", value.scale.z, "Scale shape.");
  add_option(cmd, "scaleu", value.scaleu, "Scale shape.");
  add_option(cmd, "rotatex", value.rotate.x, "Rotate shape.");
  add_option(cmd, "rotatey", value.rotate.y, "Rotate shape.");
  add_option(cmd, "rotatez", value.rotate.z, "Rotate shape.");
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

void add_command(cli_command& cli, const string& name, fvconvert_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_argument(cmd, "shape", value.shape, "Input shape.");
  add_option(cmd, "output", value.output, "Output shape.");
  add_option(cmd, "smooth", value.smooth, "Smooth normals.");
  add_option(cmd, "facet", value.facet, "Facet normals.");
  add_option(
      cmd, "aspositions", value.aspositions, "Remove all but positions.");
  add_option(cmd, "translatex", value.translate.x, "Translate shape.");
  add_option(cmd, "translatey", value.translate.y, "Translate shape.");
  add_option(cmd, "translatez", value.translate.z, "Translate shape.");
  add_option(cmd, "scalex", value.scale.x, "Scale shape.");
  add_option(cmd, "scaley", value.scale.y, "Scale shape.");
  add_option(cmd, "scalez", value.scale.z, "Scale shape.");
  add_option(cmd, "scaleu", value.scaleu, "Scale shape.");
  add_option(cmd, "rotatex", value.rotate.x, "Rotate shape.");
  add_option(cmd, "rotatey", value.rotate.y, "Rotate shape.");
  add_option(cmd, "rotatez", value.rotate.z, "Rotate shape.");
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

void add_command(cli_command& cli, const string& name, view_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_argument(cmd, "shape", value.shape, "Input shape.");
  add_option(cmd, "output", value.output, "Output shape.");
  add_option(cmd, "addsky", value.addsky, "Add sky.");
}

#ifndef YOCTO_OPENGL

// view shapes
int run_view(const view_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

// view shapes
int run_view(const view_params& params) {
  // shape data
  auto shape = scene_shape{};

  // load mesh
  auto ioerror = ""s;
  if (path_filename(params.shape) == ".ypreset") {
    if (!make_shape_preset(shape, path_basename(params.shape), ioerror))
      print_fatal(ioerror);
  } else {
    if (!load_shape(params.shape, shape, ioerror, true)) print_fatal(ioerror);
  }

  // run view
  view_shape("yshape", params.shape, shape, params.addsky);

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

void add_command(cli_command& cli, const string& name,
    heightfield_params& value, const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_argument(cmd, "image", value.image, "Input image.");
  add_option(cmd, "output", value.output, "Output shape.");
  add_option(cmd, "smooth", value.smooth, "Smoooth normals.");
  add_option(cmd, "height", value.height, "Shape height.");
  add_option(cmd, "info", value.info, "Print info.");
  add_option(cmd, "translatex", value.translate.x, "Translate shape.");
  add_option(cmd, "translatey", value.translate.y, "Translate shape.");
  add_option(cmd, "translatez", value.translate.z, "Translate shape.");
  add_option(cmd, "scalex", value.scale.x, "Scale shape.");
  add_option(cmd, "scaley", value.scale.y, "Scale shape.");
  add_option(cmd, "scalez", value.scale.z, "Scale shape.");
  add_option(cmd, "scaleu", value.scaleu, "Scale shape.");
  add_option(cmd, "rotatex", value.rotate.x, "Rotate shape.");
  add_option(cmd, "rotatey", value.rotate.y, "Rotate shape.");
  add_option(cmd, "rotatez", value.rotate.z, "Rotate shape.");
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

struct glview_params {
  string shape = "shape.ply";
};

// Cli
void add_command(cli_command& cli, const string& name, glview_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_argument(cmd, "shape", value.shape, "Input shape.");
}

#ifndef YOCTO_OPENGL

// view shapes
int run_glview(const glview_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

static scene_model make_shapescene(const scene_shape& ioshape_) {
  // Frame camera
  auto camera_frame = [](float lens, float aspect,
                          float film = 0.036) -> frame3f {
    auto camera_dir  = normalize(vec3f{0, 0.5, 1});
    auto bbox_radius = 2.0f;
    auto camera_dist = bbox_radius * lens / (film / aspect);
    return lookat_frame(camera_dir * camera_dist, {0, 0, 0}, {0, 1, 0});
  };

  // init scene
  auto scene = scene_model{};

  // rescale shape to unit
  auto ioshape = ioshape_;
  auto bbox    = invalidb3f;
  for (auto& pos : ioshape.positions) bbox = merge(bbox, pos);
  for (auto& pos : ioshape.positions) pos -= center(bbox);
  for (auto& pos : ioshape.positions) pos /= max(size(bbox));
  // TODO(fabio): this should be a math function

  // camera
  auto& camera  = scene.cameras.emplace_back();
  camera.frame  = camera_frame(0.050, 16.0f / 9.0f, 0.036);
  camera.lens   = 0.050;
  camera.aspect = 16.0f / 9.0f;
  camera.film   = 0.036;
  camera.focus  = length(camera.frame.o - center(bbox));

  // material
  auto& shape_material     = scene.materials.emplace_back();
  shape_material.type      = scene_material_type::glossy;
  shape_material.color     = {0.5, 1, 0.5};
  shape_material.roughness = 0.2;

  // shapes
  scene.shapes.emplace_back(ioshape);

  // instances
  auto& shape_instance    = scene.instances.emplace_back();
  shape_instance.shape    = 0;
  shape_instance.material = 0;

  // done
  return scene;
}

int run_glview(const glview_params& params) {
  // loading shape
  auto ioerror = ""s;
  auto shape   = scene_shape{};
  if (!load_shape(params.shape, shape, ioerror, true)) print_fatal(ioerror);

  // create scene
  auto scene = make_shapescene(shape);

  // run viewer
  glview_scene(scene, params.shape, "");

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
  glview_params      glview      = {};
};

// Cli
void add_commands(cli_command& cli, const string& name, app_params& value,
    const string& usage) {
  cli = make_cli(name, usage);
  add_command_name(cli, "command", value.command, "Command.");
  add_command(cli, "convert", value.convert, "Convert shapes.");
  add_command(
      cli, "fvconvert", value.fvconvert, "Convert face-varying shapes.");
  add_command(cli, "view", value.view, "View shapes.");
  add_command(cli, "heightfield", value.heightfield, "Create an heightfield.");
  add_command(cli, "glview", value.glview, "View shapes with OpenGL.");
}

// Parse cli
void parse_cli(app_params& params, int argc, const char** argv) {
  auto cli = cli_command{};
  add_commands(cli, "yshape", params, "Process and view shapes.");
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
  } else if (params.command == "glview") {
    return run_glview(params.glview);
  } else {
    return print_fatal("unknown command " + params.command);
  }
}
