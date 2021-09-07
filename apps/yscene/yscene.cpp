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
#include <yocto/yocto_math.h>
#include <yocto/yocto_scene.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
#include <yocto/yocto_trace.h>
#if YOCTO_OPENGL == 1
#include <yocto_gui/yocto_glview.h>
#endif
using namespace yocto;

// convert params
struct convert_params {
  string scene     = "scene.ply";
  string output    = "out.ply";
  bool   info      = false;
  bool   validate  = false;
  string copyright = "";
};

// Cli
void add_options(const cli_command& cli, convert_params& params) {
  add_argument(cli, "scene", params.scene, "Input scene.");
  add_option(cli, "output", params.output, "Output scene.");
  add_option(cli, "info", params.info, "Print info.");
  add_option(cli, "validate", params.validate, "Validate scene.");
  add_option(cli, "copyright", params.copyright, "Set scene copyright.");
}

// convert images
void run_convert(const convert_params& params) {
  // load scene
  auto error = string{};
  auto scene = scene_data{};
  print_progress_begin("load scene");
  if (!load_scene(params.scene, scene, error)) print_fatal(error);
  print_progress_end();

  // copyright
  if (params.copyright != "") {
    scene.copyright = params.copyright;
  }

  // validate scene
  if (params.validate) {
    for (auto& error : scene_validation(scene)) print_info("error: " + error);
  }

  // print info
  if (params.info) {
    print_info("scene stats ------------");
    for (auto stat : scene_stats(scene)) print_info(stat);
  }

  // tesselate if needed
  if (!scene.subdivs.empty()) {
    print_progress_begin("tesselate subdivs");
    tesselate_subdivs(scene);
    print_progress_end();
  }

  // save scene
  print_progress_begin("save scene");
  make_scene_directories(params.output, scene);
  save_scene(params.output, scene);
  print_progress_end();
}

// info params
struct info_params {
  string scene    = "scene.ply";
  bool   validate = false;
};

// Cli
void add_options(const cli_command& cli, info_params& params) {
  add_argument(cli, "scene", params.scene, "Input scene.");
  add_option(cli, "validate", params.validate, "Validate scene.");
}

// print info for scenes
void run_info(const info_params& params) {
  // load scene
  auto error = string{};
  print_progress_begin("load scene");
  auto scene = scene_data{};
  if (!load_scene(params.scene, scene, error)) print_fatal(error);
  print_progress_end();

  // validate scene
  if (params.validate) {
    for (auto& error : scene_validation(scene)) print_info("error: " + error);
  }

  // print info
  print_info("scene stats ------------");
  for (auto stat : scene_stats(scene)) print_info(stat);
}

// render params
struct render_params : trace_params {
  string scene     = "scene.json";
  string output    = "out.png";
  string camname   = "";
  bool   addsky    = false;
  string envname   = "";
  bool   savebatch = false;
};

// Cli
void add_options(const cli_command& cli, render_params& params) {
  add_argument(cli, "scene", params.scene, "Scene filename.");
  add_option(cli, "output", params.output, "Output filename.");
  add_option(cli, "camera", params.camname, "Camera name.");
  add_option(cli, "addsky", params.addsky, "Add sky.");
  add_option(cli, "envname", params.envname, "Add environment map.");
  add_option(cli, "savebatch", params.savebatch, "Save batch.");
  add_option(
      cli, "resolution", params.resolution, "Image resolution.", {1, 4096});
  add_option(
      cli, "sampler", params.sampler, "Sampler type.", trace_sampler_names);
  add_option(cli, "falsecolor", params.falsecolor, "False color type.",
      trace_falsecolor_names);
  add_option(cli, "samples", params.samples, "Number of samples.", {1, 4096});
  add_option(cli, "bounces", params.bounces, "Number of bounces.", {1, 128});
  add_option(cli, "denoise", params.denoise, "Enable denoiser.");
  add_option(cli, "batch", params.batch, "Sample batch.");
  add_option(cli, "clamp", params.clamp, "Clamp params.", {10, flt_max});
  add_option(cli, "nocaustics", params.nocaustics, "Disable caustics.");
  add_option(cli, "envhidden", params.envhidden, "Hide environment.");
  add_option(cli, "tentfilter", params.tentfilter, "Filter image.");
  add_option(cli, "embreebvh", params.embreebvh, "Use Embree as BVH.");
  add_option(
      cli, "highqualitybvh", params.highqualitybvh, "Use high quality BVH.");
  add_option(cli, "exposure", params.exposure, "Exposure value.");
  add_option(cli, "filmic", params.filmic, "Filmic tone mapping.");
  add_option(cli, "noparallel", params.noparallel, "Disable threading.");
}

// convert images
void run_render(const render_params& params_) {
  // copy params
  auto params = params_;

  // scene loading
  auto error = string{};
  print_progress_begin("load scene");
  auto scene = scene_data{};
  if (!load_scene(params.scene, scene, error)) print_fatal(error);
  print_progress_end();

  // add sky
  if (params.addsky) add_sky(scene);

  // add environment
  if (!params.envname.empty()) {
    print_progress_begin("add environment");
    add_environment(scene, params.envname);
    print_progress_end();
  }

  // camera
  params.camera = find_camera(scene, params.camname);

  // tesselation
  if (!scene.subdivs.empty()) {
    print_progress_begin("tesselate subdivs");
    tesselate_subdivs(scene);
    print_progress_end();
  }

  // build bvh
  print_progress_begin("build bvh");
  auto bvh = make_bvh(scene, params);
  print_progress_end();

  // init renderer
  print_progress_begin("build lights");
  auto lights = make_lights(scene, params);
  print_progress_end();

  // fix renderer type if no lights
  if (lights.lights.empty() && is_sampler_lit(params)) {
    print_info("no lights presents, image will be black");
    params.sampler = trace_sampler_type::eyelight;
  }

  // state
  print_progress_begin("init state");
  auto state = make_state(scene, params);
  print_progress_end();

  // render
  print_progress_begin("render image", params.samples);
  for (auto sample = 0; sample < params.samples; sample++) {
    trace_samples(state, scene, bvh, lights, params);
    if (params.savebatch && state.samples % params.batch == 0) {
      auto image = params.denoise ? get_denoised(state) : get_render(state);
      auto ext = "-s" + std::to_string(sample) + path_extension(params.output);
      auto outfilename = replace_extension(params.output, ext);
      if (!is_hdr_filename(params.output))
        image = tonemap_image(image, params.exposure, params.filmic);
      if (!save_image(outfilename, image, error)) print_fatal(error);
    }
    print_progress_next();
  }

  // save image
  print_progress_begin("save image");
  auto image = params.denoise ? get_denoised(state) : get_render(state);
  if (!is_hdr_filename(params.output))
    image = tonemap_image(image, params.exposure, params.filmic);
  if (!save_image(params.output, image, error)) print_fatal(error);
  print_progress_end();
}

// convert params
struct view_params : trace_params {
  string scene   = "scene.json";
  string output  = "out.png";
  string camname = "";
  bool   addsky  = false;
  string envname = "";
};

// Cli
void add_options(const cli_command& cli, view_params& params) {
  add_argument_with_config(
      cli, "scene", params.scene, "Scene filename.", "yscene.json");
  add_option(cli, "output", params.output, "Output filename.");
  add_option(cli, "camera", params.camname, "Camera name.");
  add_option(cli, "addsky", params.addsky, "Add sky.");
  add_option(cli, "envname", params.envname, "Add environment map.");
  add_option(
      cli, "resolution", params.resolution, "Image resolution.", {1, 4096});
  add_option(
      cli, "sampler", params.sampler, "Sampler type.", trace_sampler_names);
  add_option(cli, "falsecolor", params.falsecolor, "False color type.",
      trace_falsecolor_names);
  add_option(cli, "samples", params.samples, "Number of samples.", {1, 4096});
  add_option(cli, "bounces", params.bounces, "Number of bounces.", {1, 128});
  add_option(cli, "denoise", params.denoise, "Enable denoiser.");
  add_option(cli, "batch", params.batch, "Sample batch.");
  add_option(cli, "clamp", params.clamp, "Clamp params.", {10, flt_max});
  add_option(cli, "nocaustics", params.nocaustics, "Disable caustics.");
  add_option(cli, "envhidden", params.envhidden, "Hide environment.");
  add_option(cli, "tentfilter", params.tentfilter, "Filter image.");
  add_option(cli, "embreebvh", params.embreebvh, "Use Embree as BVH.");
  add_option(
      cli, "highqualitybvh", params.highqualitybvh, "Use high quality BVH.");
  add_option(cli, "exposure", params.exposure, "Exposure value.");
  add_option(cli, "filmic", params.filmic, "Filmic tone mapping.");
  add_option(cli, "noparallel", params.noparallel, "Disable threading.");
}

#ifndef YOCTO_OPENGL

// view scene
void run_view(const view_params& params) { print_fatal("Opengl not compiled"); }

#else

// view scene
void run_view(const view_params& params_) {
  // copy params
  auto params = params_;

  // load scene
  auto error = string{};
  print_progress_begin("load scene");
  auto scene = scene_data{};
  if (!load_scene(params.scene, scene, error)) print_fatal(error);
  print_progress_end();

  // add sky
  if (params.addsky) add_sky(scene);

  // add environment
  if (!params.envname.empty()) {
    print_progress_begin("add environment");
    add_environment(scene, params.envname);
    print_progress_end();
  }

  // tesselation
  if (!scene.subdivs.empty()) {
    print_progress_begin("tesselate subdivs");
    tesselate_subdivs(scene);
    print_progress_end();
  }

  // find camera
  params.camera = find_camera(scene, params.camname);

  // run view
  view_scene("yscene", params.scene, scene, params);
}

#endif

struct glview_params {
  string scene   = "scene.json"s;
  string camname = "";
};

// Cli
void add_options(const cli_command& cli, glview_params& params) {
  add_argument(cli, "scene", params.scene, "Input scene.");
  add_option(cli, "camera", params.camname, "Camera name.");
}

#ifndef YOCTO_OPENGL

// view scene
void run_glview(const glview_params& params) {
  print_fatal("Opengl not compiled");
}

#else

void run_glview(const glview_params& params_) {
  // copy params
  auto params = params_;

  // loading scene
  auto error = string{};
  print_progress_begin("load scene");
  auto scene = scene_data{};
  if (!load_scene(params.scene, scene, error)) print_fatal(error);
  print_progress_end();

  // tesselation
  if (!scene.subdivs.empty()) {
    print_progress_begin("tesselate subdivs");
    tesselate_subdivs(scene);
    print_progress_end();
  }

  // camera
  auto glparams   = glscene_params{};
  glparams.camera = find_camera(scene, params.camname);

  // run viewer
  glview_scene("yscene", params.scene, scene, glparams);
}

#endif

struct app_params {
  string         command = "convert";
  convert_params convert = {};
  info_params    info    = {};
  render_params  render  = {};
  view_params    view    = {};
  glview_params  glview  = {};
};

// Cli
void add_options(const cli_command& cli, app_params& params) {
  set_command_var(cli, params.command);
  add_command(cli, "convert", params.convert, "Convert scenes.");
  add_command(cli, "info", params.info, "Print scenes info.");
  add_command(cli, "render", params.render, "Render scenes.");
  add_command(cli, "view", params.view, "View scenes.");
  add_command(cli, "glview", params.glview, "View scenes with OpenGL.");
}

// Run
void run(const vector<string>& args) {
  // command line parameters
  auto error  = string{};
  auto params = app_params{};
  auto cli    = make_cli("yscene", params, "Process and view scenes.");
  if (!parse_cli(cli, args, error)) print_fatal(error);

  // dispatch commands
  if (params.command == "convert") {
    return run_convert(params.convert);
  } else if (params.command == "info") {
    return run_info(params.info);
  } else if (params.command == "render") {
    return run_render(params.render);
  } else if (params.command == "view") {
    return run_view(params.view);
  } else if (params.command == "glview") {
    return run_glview(params.glview);
  } else {
    print_fatal("yscene; unknown command");
  }
}

// Main
int main(int argc, const char* argv[]) { run(make_cli_args(argc, argv)); }
