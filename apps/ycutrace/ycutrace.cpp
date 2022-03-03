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
#include <yocto/yocto_cutrace.h>
#include <yocto/yocto_gui.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_scene.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
#include <yocto/yocto_trace.h>

using namespace yocto;

#include <filesystem>
namespace fs = std::filesystem;

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
  add_option(cli, "exposure", params.exposure, "exposure value");
  add_option(cli, "filmic", params.filmic, "filmic tone mapping");
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

  // triangulation
  for (auto& shape : scene.shapes) {
    if (shape.quads.empty()) continue;
    shape.triangles = quads_to_triangles(shape.quads);
    shape.quads     = {};
  }

  // slice params
  auto params__ = (cutrace_params&)params;

  // initialize context
  timer        = simple_timer{};
  auto context = make_cutrace_context(params__);
  print_info("init gpu: {}", elapsed_formatted(timer));

  // upload scene to the gpu
  timer        = simple_timer{};
  auto cuscene = make_cutrace_scene(scene, params__);
  print_info("upload scene: {}", elapsed_formatted(timer));

  // build bvh
  timer    = simple_timer{};
  auto bvh = make_cutrace_bvh(context, cuscene, scene, params__);
  print_info("build bvh: {}", elapsed_formatted(timer));

  // init lights
  auto lights = make_cutrace_lights(scene, params__);

  // state
  auto state = make_cutrace_state(scene, params__);

  // render
  timer = simple_timer{};
  trace_start(context, state, cuscene, bvh, lights, scene, params__);
  for (auto sample : range(0, params__.samples, params__.batch)) {
    auto sample_timer = simple_timer{};
    trace_samples(context, state, cuscene, bvh, lights, scene, params__);
    print_info("render sample {}/{}: {}", state.samples, params.samples,
        elapsed_formatted(sample_timer));
    if (params.savebatch && state.samples % params.batch == 0) {
      auto image       = params.denoise ? get_denoised_image(state)
                                        : get_rendered_image(state);
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
  print_info("render image: {}", elapsed_formatted(timer));

  // save image
  timer      = simple_timer{};
  auto image = params.denoise ? get_denoised_image(state)
                              : get_rendered_image(state);
  if (!is_hdr_filename(params.output))
    image = tonemap_image(image, params.exposure, params.filmic);
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
  add_option(cli, "exposure", params.exposure, "exposure value");
  add_option(cli, "filmic", params.filmic, "filmic tone mapping");
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
  show_trace_gui("ycuda", params.scene, scene, params);
}

struct app_params {
  string        command = "render";
  render_params render  = {};
  view_params   view    = {};
};

// Run
int main(int argc, const char* argv[]) {
  try {
    // command line parameters
    auto params = app_params{};
    auto cli    = make_cli("ycuda", "render and view scenes with cuda");
    add_command_var(cli, params.command);
    add_command(cli, "render", params.render, "render scenes");
    add_command(cli, "view", params.view, "view scenes");
    parse_cli(cli, argc, argv);

    // dispatch commands
    if (params.command == "render") {
      run_render(params.render);
    } else if (params.command == "view") {
      run_view(params.view);
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
