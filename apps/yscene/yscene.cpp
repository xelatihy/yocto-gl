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
#include <yocto/yocto_gui.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_scene.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
#include <yocto/yocto_trace.h>

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
void add_options(cli_command& cli, convert_params& params) {
  add_option(cli, "scene", params.scene, "input scene");
  add_option(cli, "output", params.output, "output scene");
  add_option(cli, "info", params.info, "print info");
  add_option(cli, "validate", params.validate, "validate scene");
  add_option(cli, "copyright", params.copyright, "set scene copyright");
}

// convert images
void run_convert(const convert_params& params) {
  print_info("converting {}", params.scene);
  auto timer = simple_timer{};

  // load scene
  timer      = simple_timer{};
  auto scene = load_scene(params.scene);
  print_info("load scene: {}", elapsed_formatted(timer));

  // copyright
  if (params.copyright != "") {
    scene.copyright = params.copyright;
  }

  // validate scene
  if (params.validate) {
    auto errors = scene_validation(scene);
    for (auto& error : errors) print_error(error);
    if (!errors.empty()) throw io_error{"invalid scene"};
  }

  // print info
  if (params.info) {
    print_info("scene stats ------------");
    for (auto stat : scene_stats(scene)) print_info(stat);
  }

  // tesselate if needed
  if (!scene.subdivs.empty()) {
    timer = simple_timer{};
    tesselate_subdivs(scene);
    print_info("tesselate subdivs: {}", elapsed_formatted(timer));
  }

  // save scene
  timer = simple_timer{};
  make_scene_directories(params.output, scene);
  save_scene(params.output, scene);
  print_info("save scene: {}", elapsed_formatted(timer));
}

// info params
struct info_params {
  string scene    = "scene.ply";
  bool   validate = false;
};

// Cli
void add_options(cli_command& cli, info_params& params) {
  add_option(cli, "scene", params.scene, "input scene");
  add_option(cli, "validate", params.validate, "validate scene");
}

// print info for scenes
void run_info(const info_params& params) {
  print_info("info for {}", params.scene);
  auto timer = simple_timer{};

  // load scene
  timer      = simple_timer{};
  auto scene = load_scene(params.scene);
  print_info("load scene: {}" + elapsed_formatted(timer));

  // validate scene
  if (params.validate) {
    for (auto& error : scene_validation(scene)) print_error(error);
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
void add_options(cli_command& cli, render_params& params) {
  add_option(cli, "scene", params.scene, "scene filename");
  add_option(cli, "output", params.output, "output filename");
  add_option(cli, "camera", params.camname, "camera name");
  add_option(cli, "addsky", params.addsky, "add sky");
  add_option(cli, "envname", params.envname, "add environment");
  add_option(cli, "savebatch", params.savebatch, "save batch");
  add_option(cli, "resolution", params.resolution, "image resolution");
  add_option(
      cli, "sampler", params.sampler, "sampler type", trace_sampler_labels);
  add_option(cli, "falsecolor", params.falsecolor, "false color type",
      trace_falsecolor_labels);
  add_option(cli, "samples", params.samples, "number of samples");
  add_option(cli, "bounces", params.bounces, "number of bounces");
  add_option(cli, "denoise", params.denoise, "enable denoiser");
  add_option(cli, "batch", params.batch, "sample batch");
  add_option(cli, "clamp", params.clamp, "clamp params");
  add_option(cli, "nocaustics", params.nocaustics, "disable caustics");
  add_option(cli, "envhidden", params.envhidden, "hide environment");
  add_option(cli, "tentfilter", params.tentfilter, "filter image");
  add_option(cli, "embreebvh", params.embreebvh, "use Embree bvh");
  add_option(cli, "highqualitybvh", params.highqualitybvh, "high quality bvh");
  add_option(cli, "noparallel", params.noparallel, "disable threading");
}

// convert images
void run_render(const render_params& params_) {
  print_info("rendering {}", params_.scene);
  auto timer = simple_timer{};

  // copy params
  auto params = params_;

  // scene loading
  timer      = simple_timer{};
  auto scene = load_scene(params.scene);
  print_info("load scene: {}", elapsed_formatted(timer));

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
  timer    = simple_timer{};
  auto bvh = make_trace_bvh(scene, params);
  print_info("build bvh: {}", elapsed_formatted(timer));

  // init renderer
  auto lights = make_trace_lights(scene, params);

  // fix renderer type if no lights
  if (lights.lights.empty() && is_sampler_lit(params)) {
    print_info("no lights presents, image will be black");
    params.sampler = trace_sampler_type::eyelight;
  }

  // state
  auto state = make_trace_state(scene, params);

  // render
  timer = simple_timer{};
  for (auto sample : range(0, params.samples, params.batch)) {
    auto sample_timer = simple_timer{};
    trace_samples(state, scene, bvh, lights, params);
    print_info("render sample {}/{}: {}", state.samples, params.samples,
        elapsed_formatted(sample_timer));
    if (params.savebatch && state.samples % params.batch == 0) {
      auto image     = get_image(state);
      auto batchname = replace_extension(params.output,
          "-" + std::to_string(state.samples) + path_extension(params.output));
      save_image(batchname, image);
    }
  }
  print_info("render image: {}", elapsed_formatted(timer));

  // save image
  timer      = simple_timer{};
  auto image = get_image(state);
  save_image(params.output, image);
  print_info("save image: {}", elapsed_formatted(timer));
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
void add_options(cli_command& cli, view_params& params) {
  add_option(cli, "scene", params.scene, "scene filename");
  add_option(cli, "output", params.output, "output filename");
  add_option(cli, "camera", params.camname, "camera name");
  add_option(cli, "addsky", params.addsky, "add sky");
  add_option(cli, "envname", params.envname, "add environment");
  add_option(cli, "resolution", params.resolution, "image resolution");
  add_option(
      cli, "sampler", params.sampler, "sampler type", trace_sampler_labels);
  add_option(cli, "falsecolor", params.falsecolor, "false color type",
      trace_falsecolor_labels);
  add_option(cli, "samples", params.samples, "number of samples");
  add_option(cli, "bounces", params.bounces, "number of bounces");
  add_option(cli, "denoise", params.denoise, "enable denoiser");
  add_option(cli, "batch", params.batch, "sample batch");
  add_option(cli, "clamp", params.clamp, "clamp params");
  add_option(cli, "nocaustics", params.nocaustics, "disable caustics");
  add_option(cli, "envhidden", params.envhidden, "hide environment");
  add_option(cli, "tentfilter", params.tentfilter, "filter image");
  add_option(cli, "embreebvh", params.embreebvh, "use Embree as BVH");
  add_option(
      cli, "--highqualitybvh", params.highqualitybvh, "use high quality BVH");
  add_option(cli, "noparallel", params.noparallel, "disable threading");
}

// view scene
void run_view(const view_params& params_) {
  print_info("viewing {}", params_.scene);
  auto timer = simple_timer{};

  // copy params
  auto params = params_;

  // load scene
  timer      = simple_timer{};
  auto scene = load_scene(params.scene);
  print_info("load scene: {}", elapsed_formatted(timer));

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
void add_options(cli_command& cli, glview_params& params) {
  add_option(cli, "scene", params.scene, "input scene");
  add_option(cli, "camera", params.camname, "camera name");
}

void run_glview(const glview_params& params_) {
  print_info("viewing {}", params_.scene);

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
  try {
    // command line parameters
    auto params = app_params{};
    auto cli    = make_cli("yscene", "process and view scenes");
    add_command_var(cli, params.command);
    add_command(cli, "convert", params.convert, "convert scenes");
    add_command(cli, "info", params.info, "print scenes info");
    add_command(cli, "render", params.render, "render scenes");
    add_command(cli, "view", params.view, "view scenes");
    add_command(cli, "glview", params.glview, "view scenes with OpenGL");
    parse_cli(cli, argc, argv);

    // dispatch commands
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
  } catch (const std::exception& error) {
    print_error(error.what());
    return 1;
  }

  // done
  return 0;
}
