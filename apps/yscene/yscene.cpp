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

#include <yocto/yocto_math.h>
#include <yocto/yocto_scene.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
#include <yocto/yocto_trace.h>
#if YOCTO_OPENGL == 1
#include <yocto/yocto_glview.h>
#endif
#include <fmt/chrono.h>
#include <fmt/core.h>

using namespace yocto;

#include <CLI/CLI.hpp>
#include <chrono>
using clock_ = std::chrono::high_resolution_clock;

// convert params
struct convert_params {
  string scene     = "scene.ply";
  string output    = "out.ply";
  bool   info      = false;
  bool   validate  = false;
  string copyright = "";
};

// Cli
void add_options(CLI::App& cli, convert_params& params) {
  cli.add_option("scene", params.scene, "Input scene.");
  cli.add_option("--output", params.output, "Output scene.");
  cli.add_flag("--info", params.info, "Print info.");
  cli.add_flag("--validate", params.validate, "Validate scene.");
  cli.add_option("--copyright", params.copyright, "Set scene copyright.");
}

// convert images
int run_convert(const convert_params& params) {
  fmt::print("converting {}\n", params.scene);
  auto start_time = clock_::now();

  // load scene
  auto error = string{};
  auto scene = scene_data{};
  start_time = clock_::now();
  if (!load_scene(params.scene, scene, error)) {
    fmt::print("error: cannot load {}\n", params.scene);
    return 1;
  }
  fmt::print("load scene: {:%H:%M:%S}\n", clock_::now() - start_time);

  // copyright
  if (params.copyright != "") {
    scene.copyright = params.copyright;
  }

  // validate scene
  if (params.validate) {
    for (auto& error : scene_validation(scene))
      fmt::print("error: {}\n", error);
  }

  // print info
  if (params.info) {
    fmt::print("scene stats ------------\n");
    for (auto stat : scene_stats(scene)) fmt::print("{}\n", stat);
  }

  // tesselate if needed
  if (!scene.subdivs.empty()) {
    auto start_time = clock_::now();
    tesselate_subdivs(scene);
    fmt::print("tesselate subdivs: {}\n", clock_::now() - start_time);
  }

  // save scene
  start_time = clock_::now();
  make_scene_directories(params.output, scene);
  if (!save_scene(params.output, scene, error)) {
    fmt::print("error: cannot save {}\n", params.output);
    return 1;
  }
  fmt::print("save scene: {:%H:%M:%S}\n", clock_::now() - start_time);

  // done
  return 0;
}

// info params
struct info_params {
  string scene    = "scene.ply";
  bool   validate = false;
};

// Cli
void add_options(CLI::App& cli, info_params& params) {
  cli.add_option("scene", params.scene, "Input scene.");
  cli.add_flag("--validate", params.validate, "Validate scene.");
}

// print info for scenes
int run_info(const info_params& params) {
  fmt::print("info for {}\n", params.scene);
  auto start_time = clock_::now();

  // load scene
  auto error = string{};
  start_time = clock_::now();
  auto scene = scene_data{};
  if (!load_scene(params.scene, scene, error)) {
    fmt::print("error: cannot load {}\n", params.scene);
    return 1;
  }
  fmt::print("load scene: {:%H:%M:%S}\n", clock_::now() - start_time);

  // validate scene
  if (params.validate) {
    for (auto& error : scene_validation(scene))
      fmt::print("error: {}\n", error);
  }

  // print info
  fmt::print("scene stats ------------\n");
  for (auto stat : scene_stats(scene)) fmt::print("{}\n", stat);

  // done
  return 0;
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
void add_options(CLI::App& cli, render_params& params) {
  cli.add_option("scene", params.scene, "Scene filename.");
  cli.add_option("--output", params.output, "Output filename.");
  cli.add_option("--camera", params.camname, "Camera name.");
  cli.add_flag("--addsky", params.addsky, "Add sky.");
  cli.add_option("--envname", params.envname, "Add environment map.");
  cli.add_flag("--savebatch", params.savebatch, "Save batch.");
  cli.add_option("--resolution", params.resolution, "Image resolution.");
  cli.add_option("--sampler", params.sampler, "Sampler type.")
      ->transform(CLI::CheckedTransformer(trace_sampler_labels));
  cli.add_option("--falsecolor", params.falsecolor, "False color type.")
      ->transform(CLI::CheckedTransformer(trace_sampler_labels));
  cli.add_option("--samples", params.samples, "Number of samples.");
  cli.add_option("--bounces", params.bounces, "Number of bounces.");
  cli.add_flag("--denoise", params.denoise, "Enable denoiser.");
  cli.add_option("--batch", params.batch, "Sample batch.");
  cli.add_option("--clamp", params.clamp, "Clamp params.");
  cli.add_flag("--nocaustics", params.nocaustics, "Disable caustics.");
  cli.add_flag("--envhidden", params.envhidden, "Hide environment.");
  cli.add_flag("--tentfilter", params.tentfilter, "Filter image.");
  cli.add_flag("--embreebvh", params.embreebvh, "Use Embree as BVH.");
  cli.add_flag(
      "--highqualitybvh", params.highqualitybvh, "Use high quality BVH.");
  cli.add_option("--exposure", params.exposure, "Exposure value.");
  cli.add_flag("--filmic", params.filmic, "Filmic tone mapping.");
  cli.add_flag("--noparallel", params.noparallel, "Disable threading.");
}

// convert images
int run_render(const render_params& params_) {
  fmt::print("rendering {}\n", params_.scene);
  auto start_time = clock_::now();

  // copy params
  auto params = params_;

  // scene loading
  auto error = string{};
  start_time = clock_::now();
  auto scene = scene_data{};
  if (!load_scene(params.scene, scene, error)) {
    fmt::print("error: cannot load {}\n", params.scene);
    return 1;
  }
  fmt::print("load scene: {:%H:%M:%S}\n", clock_::now() - start_time);

  // add sky
  if (params.addsky) add_sky(scene);

  // add environment
  if (!params.envname.empty()) {
    add_environment(scene, params.envname);
  }

  // camera
  params.camera = find_camera(scene, params.camname);

  // tesselation
  if (!scene.subdivs.empty()) {
    tesselate_subdivs(scene);
  }

  // build bvh
  auto bvh = make_bvh(scene, params);

  // init renderer
  auto lights = make_lights(scene, params);

  // fix renderer type if no lights
  if (lights.lights.empty() && is_sampler_lit(params)) {
    fmt::print("no lights presents, image will be black\n");
    params.sampler = trace_sampler_type::eyelight;
  }

  // state
  auto state = make_state(scene, params);

  // render
  start_time = clock_::now();
  for (auto sample = 0; sample < params.samples; sample++) {
    auto start_sample = clock_::now();
    trace_samples(state, scene, bvh, lights, params);
    if (params.savebatch && state.samples % params.batch == 0) {
      auto image = params.denoise ? get_denoised(state) : get_render(state);
      auto ext = "-s" + std::to_string(sample) + path_extension(params.output);
      auto outfilename = replace_extension(params.output, ext);
      if (!is_hdr_filename(params.output))
        image = tonemap_image(image, params.exposure, params.filmic);
      if (!save_image(outfilename, image, error)) {
        fmt::print("error: cannot save {}\n", outfilename);
        return 1;
      }
    }
    fmt::print("render sample {}/{}: {:%H:%M:%S}\n", sample, params.samples,
        clock_::now() - start_sample);
  }
  fmt::print("render image: {:%H:%M:%S}\n", clock_::now() - start_time);

  // save image
  start_time = clock_::now();
  auto image = params.denoise ? get_denoised(state) : get_render(state);
  if (!is_hdr_filename(params.output))
    image = tonemap_image(image, params.exposure, params.filmic);
  if (!save_image(params.output, image, error)) {
    fmt::print("error: cannot save {}\n", params.output);
    return 1;
  }
  fmt::print("save image: {:%H:%M:%S}\n", clock_::now() - start_time);

  // done
  return 0;
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
void add_options(CLI::App& cli, view_params& params) {
  cli.add_option("scene", params.scene, "Scene filename.");
  cli.add_option("--output", params.output, "Output filename.");
  cli.add_option("--camera", params.camname, "Camera name.");
  cli.add_flag("--addsky", params.addsky, "Add sky.");
  cli.add_option("--envname", params.envname, "Add environment map.");
  cli.add_option("--resolution", params.resolution, "Image resolution.");
  cli.add_option("--sampler", params.sampler, "Sampler type.")
      ->transform(CLI::CheckedTransformer(trace_sampler_labels));
  cli.add_option("--falsecolor", params.falsecolor, "False color type.")
      ->transform(CLI::CheckedTransformer(trace_falsecolor_labels));
  cli.add_option("--samples", params.samples, "Number of samples.");
  cli.add_option("--bounces", params.bounces, "Number of bounces.");
  cli.add_flag("--denoise", params.denoise, "Enable denoiser.");
  cli.add_option("--batch", params.batch, "Sample batch.");
  cli.add_option("--clamp", params.clamp, "Clamp params.");
  cli.add_flag("--nocaustics", params.nocaustics, "Disable caustics.");
  cli.add_flag("--envhidden", params.envhidden, "Hide environment.");
  cli.add_flag("--tentfilter", params.tentfilter, "Filter image.");
  cli.add_flag("--embreebvh", params.embreebvh, "Use Embree as BVH.");
  cli.add_flag(
      "--highqualitybvh", params.highqualitybvh, "Use high quality BVH.");
  cli.add_option("--exposure", params.exposure, "Exposure value.");
  cli.add_flag("--filmic", params.filmic, "Filmic tone mapping.");
  cli.add_flag("--noparallel", params.noparallel, "Disable threading.");
}

#ifndef YOCTO_OPENGL

// view scene
int run_view(const view_params& params) {
  fmt::print("error: opengl not compiled\n");
  return 1;
}

#else

// view scene
int run_view(const view_params& params_) {
  fmt::print("viewing {}\n", params_.scene);
  auto start_time = clock_::now();

  // copy params
  auto params = params_;

  // load scene
  auto error = string{};
  start_time = clock_::now();
  auto scene = scene_data{};
  if (!load_scene(params.scene, scene, error)) {
    fmt::print("error: cannot load {}\n", params.scene);
    return 1;
  }
  fmt::print("load scene: {:%H:%M:%S}\n", clock_::now() - start_time);

  // add sky
  if (params.addsky) add_sky(scene);

  // add environment
  if (!params.envname.empty()) {
    add_environment(scene, params.envname);
  }

  // tesselation
  if (!scene.subdivs.empty()) {
    tesselate_subdivs(scene);
  }

  // find camera
  params.camera = find_camera(scene, params.camname);

  // run view
  show_trace_gui("yscene", params.scene, scene, params);

  // done
  return 0;
}

#endif

struct glview_params {
  string scene   = "scene.json";
  string camname = "";
};

// Cli
void add_options(CLI::App& cli, glview_params& params) {
  cli.add_option("scene", params.scene, "Input scene.");
  cli.add_option("--camera", params.camname, "Camera name.");
}

#ifndef YOCTO_OPENGL

// view scene
int run_glview(const glview_params& params) {
  fmt::print("error: opengl not compiled\n");
  return 1;
}

#else

int run_glview(const glview_params& params_) {
  fmt::print("viewing {}\n", params_.scene);

  // copy params
  auto params = params_;

  // loading scene
  auto error = string{};
  auto scene = scene_data{};
  if (!load_scene(params.scene, scene, error)) {
    fmt::print("error: cannot load {}\n", params.scene);
    return 1;
  }

  // tesselation
  if (!scene.subdivs.empty()) {
    tesselate_subdivs(scene);
  }

  // camera
  auto viewparams   = shade_params{};
  viewparams.camera = find_camera(scene, params.camname);

  // run viewer
  show_shade_gui("yscene", params.scene, scene, viewparams);

  // done
  return 0;
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

// Run
int main(int argc, const char* argv[]) {
  // command line parameters
  auto params = app_params{};
  auto cli    = CLI::App("Process and view scenes");
  add_options(
      *cli.add_subcommand("convert", "Convert scenes."), params.convert);
  add_options(*cli.add_subcommand("info", "Print scenes info."), params.info);
  add_options(*cli.add_subcommand("render", "Render scenes."), params.render);
  add_options(*cli.add_subcommand("view", "View scenes."), params.view);
  add_options(
      *cli.add_subcommand("glview", "View scenes with OpenGL."), params.glview);
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
  } else if (params.command == "info") {
    return run_info(params.info);
  } else if (params.command == "render") {
    return run_render(params.render);
  } else if (params.command == "view") {
    return run_view(params.view);
  } else if (params.command == "glview") {
    return run_glview(params.glview);
  } else {
    fmt::print("error: unknown command\n");
    return 1;
  }
}
