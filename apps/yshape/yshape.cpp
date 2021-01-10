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

#include <yocto/yocto_commonio.h>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_math.h>
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
};

// Json IO
void serialize_value(json_mode mode, json_value& json, convert_params& value,
    const string& description) {
  serialize_object(mode, json, value, description);
  serialize_property(mode, json, value.shape, "shape", "Input shape.", true);
  serialize_property(mode, json, value.output, "output", "Output shape.");
  serialize_property(mode, json, value.smooth, "smooth", "Smooth normals.");
  serialize_property(mode, json, value.facet, "facet", "Facet normals.");
  serialize_property(mode, json, value.aspositions, "aspositions",
      "Remove all but positions.");
  serialize_property(
      mode, json, value.astriangles, "astriangles", "Convert to triangles.");
  serialize_property(mode, json, (array<float, 3>&)value.translate, "translate",
      "Translate shape.");
  serialize_property(
      mode, json, (array<float, 3>&)value.scale, "scale", "Scale shape.");
  serialize_property(
      mode, json, (array<float, 3>&)value.rotate, "rotate", "Rotate shape.");
  serialize_clipositionals(mode, json, {"shape"});
  serialize_clialternates(mode, json, {{"output", "o"}});
}

// convert images
int run_convert(const convert_params& params) {
  // shape data
  auto shape = shape_data{};

  // load mesh
  auto ioerror = ""s;
  print_progress("load shape", 0, 1);
  if (!load_shape(params.shape, shape, ioerror, false)) print_fatal(ioerror);
  print_progress("load shape", 1, 1);

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
      params.scale != vec3f{1, 1, 1}) {
    print_progress("transform shape", 0, 1);
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
    print_progress("transform shape", 1, 1);
  }

  // compute normals
  if (params.smooth) {
    print_progress("smooth shape", 0, 1);
    if (!shape.points.empty()) {
      shape.normals = vector<vec3f>{shape.positions.size(), {0, 0, 1}};
    } else if (!shape.lines.empty()) {
      shape.normals = lines_tangents(shape.lines, shape.positions);
    } else if (!shape.triangles.empty()) {
      shape.normals = triangles_normals(shape.triangles, shape.positions);
    } else if (!shape.quads.empty()) {
      shape.normals = quads_normals(shape.quads, shape.positions);
    }
    print_progress("smooth shape", 1, 1);
  }

  // remove normals
  if (params.facet) {
    print_progress("facet shape", 0, 1);
    shape.normals = {};
    print_progress("facet shape", 1, 1);
  }

  if (params.info) {
    print_info("shape stats ------------");
    auto stats = shape_stats(shape);
    for (auto& stat : stats) print_info(stat);
  }

  // save mesh
  print_progress("save shape", 0, 1);
  if (!save_shape(params.output, shape, ioerror, false)) print_fatal(ioerror);
  print_progress("save shape", 1, 1);

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
};

// Json IO
void serialize_value(json_mode mode, json_value& json, fvconvert_params& value,
    const string& description) {
  serialize_object(mode, json, value, description);
  serialize_property(mode, json, value.shape, "shape", "Input shape.", true);
  serialize_property(mode, json, value.output, "output", "Output shape.");
  serialize_property(mode, json, value.smooth, "smooth", "Smooth normals.");
  serialize_property(mode, json, value.facet, "facet", "Facet normals.");
  serialize_property(mode, json, value.aspositions, "aspositions",
      "Remove all but positions.");
  serialize_property(mode, json, (array<float, 3>&)value.translate, "translate",
      "Translate shape.");
  serialize_property(
      mode, json, (array<float, 3>&)value.scale, "scale", "Scale shape.");
  serialize_property(
      mode, json, (array<float, 3>&)value.rotate, "rotate", "Rotate shape.");
  serialize_clipositionals(mode, json, {"shape"});
  serialize_clialternates(mode, json, {{"output", "o"}});
}

// convert images
int run_fvconvert(const fvconvert_params& params) {
  // mesh data
  auto shape = fvshape_data{};

  // load mesh
  auto ioerror = ""s;
  print_progress("load shape", 0, 1);
  if (path_filename(params.shape) == ".ypreset") {
    if (!make_fvshape_preset(shape, path_basename(params.shape), ioerror))
      print_fatal(ioerror);
  } else {
    if (!load_fvshape(params.shape, shape, ioerror)) print_fatal(ioerror);
  }
  print_progress("load shape", 1, 1);

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
      params.scale != vec3f{1, 1, 1}) {
    print_progress("transform shape", 0, 1);
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
    print_progress("transform shape", 1, 1);
  }

  // compute normals
  if (params.smooth) {
    print_progress("smooth shape", 0, 1);
    if (!shape.quadspos.empty()) {
      shape.normals = quads_normals(shape.quadspos, shape.positions);
      if (!shape.quadspos.empty()) shape.quadsnorm = shape.quadspos;
    }
    print_progress("smooth shape", 1, 1);
  }

  // remove normals
  if (params.facet) {
    print_progress("facet shape", 0, 1);
    shape.normals   = {};
    shape.quadsnorm = {};
    print_progress("facet shape", 1, 1);
  }

  if (params.info) {
    print_info("shape stats ------------");
    auto stats = fvshape_stats(shape);
    for (auto& stat : stats) print_info(stat);
  }

  // save mesh
  print_progress("save shape", 0, 1);
  if (!save_fvshape(params.output, shape, ioerror, true)) print_fatal(ioerror);
  print_progress("save shape", 1, 1);

  // done
  return 0;
}

// view params
struct view_params {
  string shape  = "shape.ply";
  string output = "out.ply";
  bool   addsky = false;
};

// Json IO
void serialize_value(json_mode mode, json_value& json, view_params& value,
    const string& description) {
  serialize_object(mode, json, value, description);
  serialize_property(mode, json, value.shape, "shape", "Input shape.", true);
  serialize_property(mode, json, value.output, "output", "Output shape.");
  serialize_property(mode, json, value.addsky, "addsky", "Add sky.");
  serialize_clipositionals(mode, json, {"shape"});
  serialize_clialternates(mode, json, {{"output", "o"}});
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
  auto shape = shape_data{};

  // load mesh
  auto ioerror = ""s;
  print_progress("load shape", 0, 1);
  if (path_filename(params.shape) == ".ypreset") {
    if (!make_shape_preset(shape, path_basename(params.shape), ioerror))
      print_fatal(ioerror);
  } else {
    if (!load_shape(params.shape, shape, ioerror, false)) print_fatal(ioerror);
  }
  print_progress("load shape", 1, 1);

  // run view
  view_shape("yshape", params.shape, shape, params.addsky, print_progress);

  // done
  return 0;
}

#endif

struct heightfield_params {
  string image     = "heightfield.png"s;
  string output    = "out.ply"s;
  bool   smooth    = false;
  float  height    = 1.0f;
  vec3f  scale     = {1, 1, 1};
  vec3f  rotate    = {0, 0, 0};
  vec3f  translate = {0, 0, 0};
  bool   info      = false;
};

// Json IO
void serialize_value(json_mode mode, json_value& json,
    heightfield_params& value, const string& description) {
  serialize_object(mode, json, value, description);
  serialize_property(mode, json, value.image, "image", "Input image.", true);
  serialize_property(mode, json, value.output, "output", "Output shape.");
  serialize_property(mode, json, value.smooth, "smooth", "Smoooth normals.");
  serialize_property(mode, json, value.height, "height", "Shape height.");
  serialize_property(
      mode, json, (array<float, 3>&)value.scale, "scale", "Scale shape.");
  serialize_property(
      mode, json, (array<float, 3>&)value.rotate, "rotate", "Rotate shape.");
  serialize_property(mode, json, (array<float, 3>&)value.scale, "translate",
      "Translate shape.");
  serialize_property(mode, json, value.info, "info", "Print info.");
  serialize_clipositionals(mode, json, {"image"});
  serialize_clialternates(mode, json, {{"output", "o"}});
}

int run_heightfield(const heightfield_params& params) {
  // load mesh
  auto image   = image_data{};
  auto ioerror = ""s;
  print_progress("load image", 0, 1);
  if (!load_image(params.image, image, ioerror)) print_fatal(ioerror);
  print_progress("load image", 1, 1);

  // convert to float
  if (!!image.pixelsf.empty())
    image = convert_image(image, image.linear, false);

  // adjust height
  if (params.height != 1) {
    for (auto& pixel : image.pixelsf) pixel *= params.height;
  }

  // create heightfield
  auto shape = make_heightfield({image.width, image.height}, image.pixelsf);
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
    print_progress("transform shape", 0, 1);
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
    print_progress("transform shape", 1, 1);
  }
  // save mesh
  print_progress("save shape", 0, 1);
  if (!save_shape(params.output, shape, ioerror)) print_fatal(ioerror);
  print_progress("save shape", 1, 1);

  // done
  return 0;
}

struct app_params {
  string           command   = "convert";
  convert_params   convert   = {};
  fvconvert_params fvconvert = {};
  view_params      view      = {};
};

// Json IO
void serialize_value(json_mode mode, json_value& json, app_params& value,
    const string& description) {
  serialize_object(mode, json, value, description);
  serialize_command(mode, json, value.command, "command", "Command.");
  serialize_property(mode, json, value.convert, "convert", "Convert shapes.");
  serialize_property(
      mode, json, value.fvconvert, "fvconvert", "Convert face-varying shapes.");
  serialize_property(mode, json, value.view, "view", "View shapes.");
}

int main(int argc, const char* argv[]) {
  // command line parameters
  auto params = app_params{};
  parse_cli(params, "Process and view shapes", argc, argv);

  // dispatch commands
  if (params.command == "convert") {
    return run_convert(params.convert);
  } else if (params.command == "fvconvert") {
    return run_view(params.view);
  } else if (params.command == "view") {
    return run_view(params.view);
  } else {
    return print_fatal("unknown command " + params.command);
  }
}
