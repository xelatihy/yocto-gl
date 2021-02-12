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
#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_trace.h>
#if YOCTO_OPENGL == 1
#include <yocto_gui/yocto_glview.h>
#endif
using namespace yocto;

// render params
struct render_params : trace_params {
  string scene     = "scene.json";
  string output    = "out.png";
  string camname   = "";
  bool   addsky    = false;
  bool   savebatch = false;
};

// Cli
void add_command(cli_command& cli, const string& name, render_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_positional(cmd, "scene", value.scene, "Scene filename.");
  add_optional(cmd, "output", value.output, "Output filename.", {}, "o");
  add_optional(cmd, "camera", value.camname, "Camera name.", {}, "c");
  add_optional(cmd, "addsky", value.addsky, "Add sky.");
  add_optional(cmd, "savebatch", value.savebatch, "Save batch.");
  add_optional(
      cmd, "resolution", value.resolution, "Image resolution.", {1, 4096}, "r");
  add_optional(
      cmd, "sampler", value.sampler, "Sampler type.", trace_sampler_names, "t");
  add_optional(cmd, "falsecolor", value.falsecolor, "False color type.",
      trace_falsecolor_names, "F");
  add_optional(
      cmd, "samples", value.samples, "Number of samples.", {1, 4096}, "s");
  add_optional(
      cmd, "bounces", value.bounces, "Number of bounces.", {1, 128}, "b");
  add_optional(cmd, "clamp", value.clamp, "Clamp value.", {10, flt_max});
  add_optional(cmd, "nocaustics", value.nocaustics, "Disable caustics.");
  add_optional(cmd, "envhidden", value.envhidden, "Hide environment.");
  add_optional(cmd, "tentfilter", value.tentfilter, "Filter image.");
  add_optional(cmd, "bvh", value.bvh, "Bvh type.", trace_bvh_names);
  add_optional(cmd, "noparallel", value.noparallel, "Disable threading.");
}

// convert images
int run_render(const render_params& params_) {
  // copy params
  auto params = params_;

  // scene loading
  auto scene   = scene_scene{};
  auto ioerror = string{};
  if (!load_scene(params.scene, scene, ioerror, print_progress))
    return print_fatal(ioerror);

  // add sky
  if (params.addsky) add_sky(scene);

  // tesselation
  tesselate_shapes(scene, print_progress);

  // build bvh
  auto bvh = trace_bvh{};
  init_bvh(bvh, scene, params, print_progress);

  // init renderer
  auto lights = trace_lights{};
  init_lights(lights, scene, params, print_progress);

  // fix renderer type if no lights
  if (lights.lights.empty() && is_sampler_lit(params)) {
    print_info("no lights presents, image will be black");
    params.sampler = trace_sampler_type::eyelight;
  }

  // camera
  params.camera = find_camera(scene, params.camname);

  // render
  auto render = trace_image(scene, bvh, lights, params, print_progress,
      [savebatch = params.savebatch, output = params.output](
          const image<vec4f>& render, int sample, int samples) {
        if (!savebatch) return;
        auto ext = "-s" + std::to_string(sample + samples) +
                   path_extension(output);
        auto outfilename = replace_extension(output, ext);
        auto ioerror     = ""s;
        print_progress("save image", sample, samples);
        if (!save_image(outfilename, render, ioerror)) print_fatal(ioerror);
      });

  // save image
  print_progress("save image", 0, 1);
  if (!save_image(params.output, render, ioerror)) return print_fatal(ioerror);
  print_progress("save image", 1, 1);

  // done
  return 0;
}

// convert params
struct view_params : trace_params {
  string scene   = "scene.json";
  string output  = "out.png";
  string camname = "";
  bool   addsky  = false;
};

// Cli
void add_command(cli_command& cli, const string& name, view_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_positional(cmd, "scene", value.scene, "Scene filename.");
  add_optional(cmd, "output", value.output, "Output filename.", {}, "o");
  add_optional(cmd, "camera", value.camname, "Camera name.", {}, "c");
  add_optional(cmd, "addsky", value.addsky, "Add sky.");
  add_optional(
      cmd, "resolution", value.resolution, "Image resolution.", {1, 4096}, "r");
  add_optional(
      cmd, "sampler", value.sampler, "Sampler type.", trace_sampler_names, "t");
  add_optional(cmd, "falsecolor", value.falsecolor, "False color type.",
      trace_falsecolor_names, "F");
  add_optional(
      cmd, "samples", value.samples, "Number of samples.", {1, 4096}, "s");
  add_optional(
      cmd, "bounces", value.bounces, "Number of bounces.", {1, 128}, "b");
  add_optional(cmd, "clamp", value.clamp, "Clamp value.", {10, flt_max});
  add_optional(cmd, "nocaustics", value.nocaustics, "Disable caustics.");
  add_optional(cmd, "envhidden", value.envhidden, "Hide environment.");
  add_optional(cmd, "tentfilter", value.tentfilter, "Filter image.");
  add_optional(cmd, "bvh", value.bvh, "Bvh type.", trace_bvh_names);
  add_optional(cmd, "noparallel", value.noparallel, "Disable threading.");
}

#ifndef YOCTO_OPENGL

// convert images
int run_view(const view_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

// Parameter conversions
void from_params(const gui_params& uiparams, trace_params& params) {
  params.camera     = uiparams.at("camera");
  params.resolution = uiparams.at("resolution");
  params.sampler    = uiparams.at("sampler");
  params.falsecolor = uiparams.at("falsecolor");
  params.samples    = uiparams.at("samples");
  params.bounces    = uiparams.at("bounces");
  params.clamp      = uiparams.at("clamp");
  params.nocaustics = uiparams.at("nocaustics");
  params.envhidden  = uiparams.at("envhidden");
  params.tentfilter = uiparams.at("tentfilter");
}
void to_params(gui_params& uiparams, const trace_params& params,
    const vector<string>& camera_names) {
  uiparams["camera"]     = {params.camera, camera_names};
  uiparams["resolution"] = {params.resolution, {128, 4096}};
  uiparams["sampler"]    = {params.sampler, trace_sampler_names};
  uiparams["falsecolor"] = {params.falsecolor, trace_falsecolor_names};
  uiparams["samples"]    = {params.samples, {1, 4096}};
  uiparams["bounces"]    = {params.bounces, {1, 64}};
  uiparams["clamp"]      = {params.clamp, {10, 1000}};
  uiparams["nocaustics"] = params.nocaustics;
  uiparams["envhidden"]  = params.envhidden;
  uiparams["tentfilter"] = params.tentfilter;
}

// interactive render
int run_view(const view_params& params_) {
  // copy params
  auto params = params_;

  // open viewer
  auto viewer = make_imageviewer("yimage");

  // scene loading
  auto scene   = scene_scene{};
  auto ioerror = string{};
  if (!load_scene(params.scene, scene, ioerror, print_progress))
    return print_fatal(ioerror);

  // create camera names if missing
  if (scene.camera_names.empty()) {
    for (auto& camera : scene.cameras) {
      scene.camera_names.push_back(
          "camera" + std::to_string(&camera - scene.cameras.data()));
    }
  }

  // add sky
  if (params.addsky) add_sky(scene);

  // tesselation
  tesselate_shapes(scene, print_progress);

  // build bvh
  auto bvh = trace_bvh{};
  init_bvh(bvh, scene, params, print_progress);

  // init renderer
  auto lights = trace_lights{};
  init_lights(lights, scene, params, print_progress);

  // fix renderer type if no lights
  if (lights.lights.empty() && is_sampler_lit(params)) {
    print_info("no lights presents, image will be black");
    params.sampler = trace_sampler_type::eyelight;
  }

  // camera
  params.camera = find_camera(scene, params.camname);

  // init state
  auto state = trace_state{};

  // render start
  trace_start(
      state, scene, bvh, lights, params,
      [&](const string& message, int sample, int nsamples) {
        set_param(viewer, "render", "sample", {sample, {0, 4096}, true});
        print_progress(message, sample, nsamples);
      },
      [&](const image<vec4f>& render, int current, int total) {
        set_image(viewer, "render", render);
      });

  // show rendering params
  auto uiparams = gui_params{};
  to_params(uiparams, params, scene.camera_names);
  set_params(viewer, "render", "Render", uiparams);

  // set callback
  set_params_callback(
      viewer, [&](const string& name, const gui_params& uiparams) {
        if (name != "render") return;
        if (uiparams.empty()) return;
        trace_stop(state);
        from_params(uiparams, params);
        trace_start(
            state, scene, bvh, lights, params,
            [&](const string& message, int sample, int nsamples) {
              set_param(viewer, "render", "sample", {sample, {1, 4096}, true});
              print_progress(message, sample, nsamples);
            },
            [&](const image<vec4f>& render, int current, int total) {
              set_image(viewer, "render", render);
            });
      });

  set_input_callback(viewer, [&](const string& name, const gui_input& input) {
    if ((input.mouse_left || input.mouse_right) &&
        input.mouse_pos != input.mouse_last) {
      trace_stop(state);
      auto dolly  = 0.0f;
      auto pan    = zero2f;
      auto rotate = zero2f;
      if (input.mouse_left && !input.modifier_shift)
        rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
      if (input.mouse_right)
        dolly = (input.mouse_pos.x - input.mouse_last.x) / 100.0f;
      auto& camera = scene.cameras[params.camera];
      if (input.mouse_left && input.modifier_shift)
        pan = (input.mouse_pos - input.mouse_last) * camera.focus / 200.0f;
      pan.x                                = -pan.x;
      std::tie(camera.frame, camera.focus) = camera_turntable(
          camera.frame, camera.focus, rotate, dolly, pan);
      trace_start(
          state, scene, bvh, lights, params,
          [&](const string& message, int sample, int nsamples) {
            set_param(viewer, "render", "sample", {sample, {1, 4096}, true});
            print_progress(message, sample, nsamples);
          },
          [&](const image<vec4f>& render, int current, int total) {
            set_image(viewer, "render", render);
          });
    }
  });

  // run view
  run_viewer(viewer);

  // stop
  trace_stop(state);

  // // save image
  // print_progress("save image", 0, 1);
  // if (!save_image(params.output, render, ioerror)) return
  // print_fatal(ioerror); print_progress("save image", 1, 1);

  // done
  return 0;
}

#endif

struct app_params {
  string        command = "render";
  render_params render  = {};
  view_params   view    = {};
};

// Cli
void add_commands(cli_command& cli, const string& name, app_params& value,
    const string& usage) {
  cli = make_cli(name, usage);
  add_command_name(cli, "command", value.command, "Command.");
  add_command(cli, "render", value.render, "Render offline.");
  add_command(cli, "view", value.view, "Render interactively.");
}

// Parse cli
void parse_cli(app_params& params, int argc, const char** argv) {
  auto cli = cli_command{};
  add_commands(cli, "ytrace", params, "Render images from scenes");
  parse_cli(cli, argc, argv);
}

int main(int argc, const char* argv[]) {
  // parse cli
  auto params = app_params{};
  parse_cli(params, argc, argv);

  // dispatch commands
  if (params.command == "render") {
    return run_render(params.render);
  } else if (params.command == "view") {
    return run_view(params.view);
  } else {
    print_fatal("unknown command");
    return 1;
  }
}
