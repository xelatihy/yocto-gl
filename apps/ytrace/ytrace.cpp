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
#include <yocto/yocto_json.h>
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
  string camera    = "";
  bool   addsky    = false;
  bool   savebatch = false;
};

// Cli
void add_command(cli_state& cli, const string& name, render_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_positional(cmd, "scene", value.scene, "Scene filename.");
  add_optional(cmd, "output", value.output, "Output filename.", "o");
  add_optional(cmd, "camera", value.camera, "Camera name.", "c");
  add_optional(cmd, "addsky", value.addsky, "Add sky.");
  add_optional(cmd, "savebatch", value.savebatch, "Save batch.");
  add_optional(cmd, "resolution", value.resolution, "Image resolution.", "r");
  add_optional(cmd, "sampler", value.sampler, "Sampler type.", "t");
  add_optional(cmd, "falsecolor", value.falsecolor, "False color type.", "F");
  add_optional(cmd, "samples", value.samples, "Number of samples.", "s");
  add_optional(cmd, "bounces", value.bounces, "Number of bounces.", "b");
  add_optional(cmd, "clamp", value.clamp, "Clamp value.");
  add_optional(cmd, "nocaustics", value.nocaustics, "Disable caustics.");
  add_optional(cmd, "envhidden", value.envhidden, "Hide environment.");
  add_optional(cmd, "tentfilter", value.tentfilter, "Filter image.");
  add_optional(cmd, "seed", value.seed, "Random seed.");
  add_optional(cmd, "bvh", value.bvh, "Bvh type.");
  add_optional(cmd, "noparallel", value.noparallel, "Disable threading.");
}

// convert images
int run_render(const render_params& params) {
  // scene loading
  auto scene   = scene_scene{};
  auto ioerror = string{};
  if (!load_scene(params.scene, scene, ioerror, print_progress))
    return print_fatal(ioerror);

  // add sky
  if (params.addsky) add_sky(scene);

  // get camera
  auto camera_handle = find_camera(scene, params.camera);

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
  }

  // render
  auto& camera = scene.cameras[camera_handle];
  auto  render = trace_image(scene, camera, bvh, lights, params, print_progress,
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
  string scene  = "scene.json";
  string output = "out.png";
  string camera = "";
  bool   addsky = false;
};

// Cli
void add_command(cli_state& cli, const string& name, view_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_positional(cmd, "scene", value.scene, "Scene filename.");
  add_optional(cmd, "output", value.output, "Output filename.", "o");
  add_optional(cmd, "camera", value.camera, "Camera name.", "c");
  add_optional(cmd, "addsky", value.addsky, "Add sky.");
}

#ifndef YOCTO_OPENGL

// convert images
int run_view(const view_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

// interactive render
int run_view(const view_params& params) {
  // open viewer
  auto viewer = make_imageviewer("yimage");

  // scene loading
  auto scene   = scene_scene{};
  auto ioerror = string{};
  if (!load_scene(params.scene, scene, ioerror, print_progress))
    return print_fatal(ioerror);

  // add sky
  if (params.addsky) add_sky(scene);

  // get camera
  auto camera_handle = find_camera(scene, params.camera);

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
  }

  // init state
  auto state = trace_state{};

  // render start
  auto& camera = scene.cameras[camera_handle];
  trace_start(
      state, scene, camera, bvh, lights, params,
      [&](const string& message, int sample, int nsamples) {
        set_widget(viewer, "render", "sample", to_json(sample),
            to_schema(sample, "Current sample"));
        print_progress(message, sample, nsamples);
      },
      [&](const image<vec4f>& render, int current, int total) {
        set_image(viewer, "render", render);
      });

  // show rendering params
  set_widgets(
      viewer, "render", to_json(params), to_schema(params, "Render params"));

  // set callback
  set_callback(viewer, [&](const string& name, const json_value& uiparams,
                           const gui_input& input) {
    if (name != "render") return;
    if (!uiparams.is_null()) {
      trace_stop(state);
      (view_params&)params = from_json<view_params>(uiparams);
      // show rendering params
      set_widgets(viewer, "render", to_json(params),
          to_schema(params, "Render params"));
      trace_start(
          state, scene, camera, bvh, lights, params,
          [&](const string& message, int sample, int nsamples) {
            set_widget(viewer, "render", "sample", to_json(sample),
                to_schema(sample, "Current sample"));
            print_progress(message, sample, nsamples);
          },
          [&](const image<vec4f>& render, int current, int total) {
            set_image(viewer, "render", render);
          });
    } else if ((input.mouse_left || input.mouse_right) &&
               input.mouse_pos != input.mouse_last) {
      trace_stop(state);
      auto dolly  = 0.0f;
      auto pan    = zero2f;
      auto rotate = zero2f;
      if (input.mouse_left && !input.modifier_shift)
        rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
      if (input.mouse_right)
        dolly = (input.mouse_pos.x - input.mouse_last.x) / 100.0f;
      if (input.mouse_left && input.modifier_shift)
        pan = (input.mouse_pos - input.mouse_last) * camera.focus / 200.0f;
      pan.x                                = -pan.x;
      std::tie(camera.frame, camera.focus) = camera_turntable(
          camera.frame, camera.focus, rotate, dolly, pan);
      trace_start(
          state, scene, camera, bvh, lights, params,
          [&](const string& message, int sample, int nsamples) {
            set_widget(viewer, "render", "sample", to_json(sample),
                to_schema(sample, "Current sample"));
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
void add_commands(cli_state& cli, const string& name, app_params& value,
    const string& usage) {
  cli = make_cli(name, usage);
  add_command_name(cli, "command", value.command, "Command.");
  add_command(cli, "render", value.render, "Render offline.");
  add_command(cli, "view", value.view, "Render interactively.");
}

// Parse cli
void parse_cli(app_params& params, int argc, const char** argv) {
  auto cli = cli_state{};
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
