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

#include <yocto/yocto_gui.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_scene.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
#include <yocto/yocto_trace.h>

using namespace yocto;

#include <chrono>
#include <filesystem>
namespace fs = std::filesystem;

#include "ext/CLI11.hpp"

int64_t now() {
  return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}
string format_duration(int64_t duration) {
  auto elapsed = duration / 1000000;  // milliseconds
  auto hours   = (int)(elapsed / 3600000);
  elapsed %= 3600000;
  auto mins = (int)(elapsed / 60000);
  elapsed %= 60000;
  auto secs  = (int)(elapsed / 1000);
  auto msecs = (int)(elapsed % 1000);
  char buffer[256];
  snprintf(
      buffer, sizeof(buffer), "%02d:%02d:%02d.%03d", hours, mins, secs, msecs);
  return buffer;
}

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
void run_convert(const convert_params& params) {
  std::cout << "converting " + params.scene + "\n";
  auto start = now();

  // load scene
  start      = now();
  auto scene = load_scene(params.scene);
  std::cout << "load scene: " + format_duration(now() - start) + "\n";

  // copyright
  if (params.copyright != "") {
    scene.copyright = params.copyright;
  }

  // validate scene
  if (params.validate) {
    auto errors = scene_validation(scene);
    for (auto& error : errors) std::cerr << "error: " + error + "\n";
    if (!errors.empty()) throw io_error{"invalid scene"};
  }

  // print info
  if (params.info) {
    std::cout << "scene stats ------------\n";
    for (auto stat : scene_stats(scene)) std::cout << stat + "\n";
  }

  // tesselate if needed
  if (!scene.subdivs.empty()) {
    start = now();
    tesselate_subdivs(scene);
    std::cout << "tesselate subdivs: " + format_duration(now() - start) + "\n";
  }

  // save scene
  start = now();
  make_scene_directories(params.output, scene);
  save_scene(params.output, scene);
  std::cout << "save scene: " + format_duration(now() - start) + "\n";
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
void run_info(const info_params& params) {
  std::cout << "info for " + params.scene + "\n";
  auto start = now();

  // load scene
  start      = now();
  auto scene = load_scene(params.scene);
  std::cout << "load scene: " + format_duration(now() - start) + "\n";

  // validate scene
  if (params.validate) {
    for (auto& error : scene_validation(scene))
      std::cerr << "error: " + error + "\n";
  }

  // print info
  std::cout << "scene stats ------------\n";
  for (auto stat : scene_stats(scene)) std::cout << stat + "\n";
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
void run_render(const render_params& params_) {
  std::cout << "rendering " + params_.scene + "\n";
  auto start = now();

  // copy params
  auto params = params_;

  // scene loading
  start      = now();
  auto scene = load_scene(params.scene);
  std::cout << "load scene: " + format_duration(now() - start) + "\n";

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
    std::cout << "no lights presents, image will be black\n";
    params.sampler = trace_sampler_type::eyelight;
  }

  // state
  auto state = make_state(scene, params);

  // render
  start = now();
  for (auto sample = 0; sample < params.samples; sample++) {
    auto sample_start = now();
    trace_samples(state, scene, bvh, lights, params);
    std::cout << ("render sample " + std::to_string(sample) + "/" +
                  std::to_string(params.samples) + ":" +
                  format_duration(now() - sample_start) + "\n");
    if (params.savebatch && state.samples % params.batch == 0) {
      auto image = params.denoise ? get_denoised(state) : get_render(state);
      auto outfilename = fs::path(params.output)
                             .replace_extension(
                                 "-s" + std::to_string(sample) +
                                 fs::path(params.output).extension().string())
                             .string();
      if (!is_hdr_filename(params.output))
        image = tonemap_image(image, params.exposure, params.filmic);
      save_image(outfilename, image);
    }
  }
  std::cout << "render image: " + format_duration(now() - start) + "\n";

  // save image
  start      = now();
  auto image = params.denoise ? get_denoised(state) : get_render(state);
  if (!is_hdr_filename(params.output))
    image = tonemap_image(image, params.exposure, params.filmic);
  save_image(params.output, image);
  std::cout << "save image: " + format_duration(now() - start) + "\n";
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

// view scene
void run_view(const view_params& params_) {
  std::cout << "viewing " + params_.scene + "\n";
  auto start = now();

  // copy params
  auto params = params_;

  // load scene
  start      = now();
  auto scene = load_scene(params.scene);
  std::cout << "load scene: " + format_duration(now() - start) + "\n";

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
}

struct glview_params {
  string scene   = "scene.json";
  string camname = "";
};

// Cli
void add_options(CLI::App& cli, glview_params& params) {
  cli.add_option("scene", params.scene, "Input scene.");
  cli.add_option("--camera", params.camname, "Camera name.");
}

void run_glview(const glview_params& params_) {
  std::cout << "viewing " + params_.scene + "\n";

  // copy params
  auto params = params_;

  // loading scene
  auto scene = load_scene(params.scene);

  // tesselation
  if (!scene.subdivs.empty()) {
    tesselate_subdivs(scene);
  }

  // camera
  auto viewparams   = shade_params{};
  viewparams.camera = find_camera(scene, params.camname);

  // run viewer
  show_shade_gui("yscene", params.scene, scene, viewparams);
}

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
  try {
    if (params.command == "convert") {
      run_convert(params.convert);
    } else if (params.command == "info") {
      run_info(params.info);
    } else if (params.command == "render") {
      run_render(params.render);
    } else if (params.command == "view") {
      run_view(params.view);
    } else if (params.command == "glview") {
      run_glview(params.glview);
    } else {
      throw io_error{"unknown command"};
    }
  } catch (const io_error& error) {
    std::cerr << "error: " << error.what() << "\n";
    return 1;
  }

  // done
  return 0;
}
