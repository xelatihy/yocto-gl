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
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_shape.h>
#if YOCTO_OPENGL == 1
// include vieweer later
#endif
using namespace yocto;

// Shape presets used ofr testing.
bool make_shape_preset(vector<int>& points, vector<vec2i>& lines,
    vector<vec3i>& triangles, vector<vec4i>& quads, vector<vec4i>& quadspos,
    vector<vec4i>& quadsnorm, vector<vec4i>& quadstexcoord,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<vec4f>& colors, vector<float>& radius, const string& type,
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
  auto shape = generic_shape{};

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
      shape.normals = compute_tangents(shape.lines, shape.positions);
    } else if (!shape.triangles.empty()) {
      shape.normals = compute_normals(shape.triangles, shape.positions);
    } else if (!shape.quads.empty()) {
      shape.normals = compute_normals(shape.quads, shape.positions);
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
  auto shape = generic_shape{};

  // load mesh
  auto ioerror = ""s;
  print_progress("load shape", 0, 1);
  if (path_filename(params.shape) == ".ypreset") {
    if (!make_shape_preset(shape, path_basename(params.shape), ioerror))
      print_fatal(ioerror);
  } else {
    if (!load_shape(params.shape, shape, ioerror, true)) print_fatal(ioerror);
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
    if (!shape.quadspos.empty()) {
      shape.normals = compute_normals(shape.quadspos, shape.positions);
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
    auto stats = shape_stats(shape);
    for (auto& stat : stats) print_info(stat);
  }

  // save mesh
  print_progress("save shape", 0, 1);
  if (!save_shape(params.output, shape, ioerror, true)) print_fatal(ioerror);
  print_progress("save shape", 1, 1);

  // done
  return 0;
}

// view params
struct view_params {
  vector<string> shapes = {"shape.ply"};
  string         output = "out.ply";
};

// Json IO
void serialize_value(json_mode mode, json_value& json, view_params& value,
    const string& description) {
  serialize_object(mode, json, value, description);
  serialize_property(mode, json, value.shapes, "shapes", "Input shapes.", true);
  serialize_property(mode, json, value.output, "output", "Output shape.");
  serialize_clipositionals(mode, json, {"shapes"});
  serialize_clialternates(mode, json, {{"output", "o"}});
}

#ifndef YOCTO_OPENGL

// view shapes
int run_view(const view_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

#if 1

// view shapes
int run_view(const view_params& params) {
  return print_fatal("Not implemented yet");
}

#else

// view shapes
int run_view(const view_params& params) {
  // open viewer
  auto viewer_guard = make_sceneview("yshape");
  auto viewer       = viewer_guard.get();

  // set image
  for (auto& filename : params.shapes) {
    // load
    auto hdr     = image<vec4f>{};
    auto ldr     = image<vec4b>{};
    auto ioerror = string{};
    if (is_preset_filename(filename)) {
      if (is_preset_hdr(filename)) {
        if (!make_image_preset(path_basename(filename), hdr, ioerror))
          return print_fatal(ioerror);
      } else {
        if (!make_image_preset(path_basename(filename), ldr, ioerror))
          return print_fatal(ioerror);
      }
    } else if (is_hdr_filename(filename)) {
      if (!load_image(filename, hdr, ioerror)) return print_fatal(ioerror);
    } else {
      if (!load_image(filename, ldr, ioerror)) return print_fatal(ioerror);
    }

    // push image to the viewer
    if (!hdr.empty()) {
      set_image(viewer, filename, hdr);
    } else {
      set_image(viewer, filename, ldr);
    }
  }

  // run view
  run_view(viewer);

  // done
  return 0;
}

#endif

#endif

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
      mode, json, value.convert, "convert", "Convert face-varying shapes.");
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
