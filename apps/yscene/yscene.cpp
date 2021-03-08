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
  string scene    = "scene.ply";
  string output   = "out.ply";
  bool   info     = false;
  bool   validate = false;
};

// Cli
void add_command(cli_command& cli, const string& name, convert_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_argument(cmd, "scene", value.scene, "Input scene.");
  add_option(cmd, "output", value.output, "Output scene.", {}, "o");
  add_option(cmd, "info", value.info, "Print info.");
  add_option(cmd, "validate", value.validate, "Validate scene.");
}

// convert images
int run_convert(const convert_params& params) {
  // load scene
  auto scene   = scene_scene{};
  auto ioerror = ""s;
  if (!load_scene(params.scene, scene, ioerror)) print_fatal(ioerror);

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
  tesselate_subdivs(scene);

  // make a directory if needed
  if (!make_scene_directories(params.output, scene, ioerror))
    print_fatal(ioerror);

  // save scene
  if (!save_scene(params.output, scene, ioerror)) print_fatal(ioerror);

  // done
  return 0;
}

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
  add_argument(cmd, "scene", value.scene, "Scene filename.");
  add_option(cmd, "output", value.output, "Output filename.", {}, "o");
  add_option(cmd, "camera", value.camname, "Camera name.", {}, "c");
  add_option(cmd, "addsky", value.addsky, "Add sky.");
  add_option(cmd, "savebatch", value.savebatch, "Save batch.");
  add_option(
      cmd, "resolution", value.resolution, "Image resolution.", {1, 4096}, "r");
  add_option(
      cmd, "sampler", value.sampler, "Sampler type.", trace_sampler_names, "t");
  add_option(cmd, "falsecolor", value.falsecolor, "False color type.",
      trace_falsecolor_names, "F");
  add_option(
      cmd, "samples", value.samples, "Number of samples.", {1, 4096}, "s");
  add_option(
      cmd, "bounces", value.bounces, "Number of bounces.", {1, 128}, "b");
  add_option(cmd, "clamp", value.clamp, "Clamp value.", {10, flt_max});
  add_option(cmd, "nocaustics", value.nocaustics, "Disable caustics.");
  add_option(cmd, "envhidden", value.envhidden, "Hide environment.");
  add_option(cmd, "tentfilter", value.tentfilter, "Filter image.");
  add_option(cmd, "embreebvh", value.embreebvh, "Use Embree as BVH.");
  add_option(
      cmd, "highqualitybvh", value.highqualitybvh, "Use high quality BVH.");
  add_option(cmd, "noparallel", value.noparallel, "Disable threading.");
}

// convert images
int run_render(const render_params& params_) {
  // copy params
  auto params = params_;

  // scene loading
  auto scene   = scene_scene{};
  auto ioerror = string{};
  if (!load_scene(params.scene, scene, ioerror)) return print_fatal(ioerror);

  // add sky
  if (params.addsky) add_sky(scene);

  // camera
  params.camera = find_camera(scene, params.camname);

  // tesselation
  tesselate_subdivs(scene);

  // build bvh
  auto bvh = make_bvh(scene, params);

  // init renderer
  auto lights = make_lights(scene, params);

  // fix renderer type if no lights
  if (lights.lights.empty() && is_sampler_lit(params)) {
    print_info("no lights presents, image will be black");
    params.sampler = trace_sampler_type::eyelight;
  }

  // state
  auto state = make_state(scene, params);
  auto image = make_image(state.width, state.height, true, false);

  // render
  trace_image(
      image, state, scene, bvh, lights, params, [&](int sample, int samples) {
        if (!params.savebatch) return;
        auto ext = "-s" + std::to_string(sample + samples) +
                   path_extension(params.output);
        auto outfilename = replace_extension(params.output, ext);
        auto ioerror     = ""s;
        log_progress("save image", sample, samples);
        if (!save_image(outfilename, image, ioerror)) print_fatal(ioerror);
      });

  // save image
  log_progress("save image", 0, 1);
  if (!save_image(params.output, image, ioerror)) return print_fatal(ioerror);
  log_progress("save image", 1, 1);

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
  add_argument(cmd, "scene", value.scene, "Scene filename.");
  add_option(cmd, "output", value.output, "Output filename.", {}, "o");
  add_option(cmd, "camera", value.camname, "Camera name.");
  add_option(cmd, "addsky", value.addsky, "Add sky.");
  add_option(
      cmd, "resolution", value.resolution, "Image resolution.", {1, 4096}, "r");
  add_option(
      cmd, "sampler", value.sampler, "Sampler type.", trace_sampler_names, "t");
  add_option(cmd, "falsecolor", value.falsecolor, "False color type.",
      trace_falsecolor_names, "F");
  add_option(
      cmd, "samples", value.samples, "Number of samples.", {1, 4096}, "s");
  add_option(
      cmd, "bounces", value.bounces, "Number of bounces.", {1, 128}, "b");
  add_option(cmd, "clamp", value.clamp, "Clamp value.", {10, flt_max});
  add_option(cmd, "nocaustics", value.nocaustics, "Disable caustics.");
  add_option(cmd, "envhidden", value.envhidden, "Hide environment.");
  add_option(cmd, "tentfilter", value.tentfilter, "Filter image.");
  add_option(cmd, "embreebvh", value.embreebvh, "Use Embree as BVH.");
  add_option(
      cmd, "highqualitybvh", value.highqualitybvh, "Use high quality BVH.");
  add_option(cmd, "noparallel", value.noparallel, "Disable threading.");
}

#ifndef YOCTO_OPENGL

// view scene
int run_view(const view_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

// view scene
int run_view(const view_params& params_) {
  // copy params
  auto params = params_;

  // load scene
  auto scene   = scene_scene{};
  auto ioerror = ""s;
  if (!load_scene(params.scene, scene, ioerror)) print_fatal(ioerror);

  // add sky
  if (params.addsky) add_sky(scene);

  // tesselation
  tesselate_subdivs(scene);

  // find camera
  params.camera = find_camera(scene, params.camname);

  // run view
  view_scene("yscene", params.scene, scene, params);

  // done
  return 0;
}

#endif

struct glview_params {
  string scene = "scene.json"s;
};

// Cli
void add_command(cli_command& cli, const string& name, glview_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_argument(cmd, "scene", value.scene, "Input scene.");
}

#ifndef YOCTO_OPENGL

// view scene
int run_glview(const glview_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

int run_glview(const glview_params& params) {
  // loading scene
  auto ioerror = ""s;
  auto scene   = scene_scene{};
  if (!load_scene(params.scene, scene, ioerror)) print_fatal(ioerror);

  // tesselation
  tesselate_subdivs(scene);

  // run viewer
  glview_scene(
      scene, params.scene, "",
      [](gui_window* win, const gui_input& input, scene_scene& scene,
          shade_scene& glscene) {},
      [](gui_window* win, const gui_input& input, scene_scene& scene,
          shade_scene& glscene) {});

  // done
  return 0;
}

#endif

struct app_params {
  string         command = "convert";
  convert_params convert = {};
  render_params  render  = {};
  view_params    view    = {};
  glview_params  glview  = {};
};

// Cli
void add_commands(cli_command& cli, const string& name, app_params& value,
    const string& usage) {
  cli = make_cli(name, usage);
  add_command_name(cli, "command", value.command, "Command.");
  add_command(cli, "convert", value.convert, "Convert scenes.");
  add_command(cli, "render", value.render, "Render scenes.");
  add_command(cli, "view", value.view, "View scenes.");
  add_command(cli, "glview", value.glview, "View scenes with OpenGL.");
}

// Parse cli
void parse_cli(app_params& params, int argc, const char** argv) {
  auto cli = cli_command{};
  add_commands(cli, "yscene", params, "Process and view scenes.");
  parse_cli(cli, argc, argv);
}

int main(int argc, const char* argv[]) {
  // command line parameters
  auto params = app_params{};
  parse_cli(params, argc, argv);
  set_log_level(true);

  // dispatch commands
  if (params.command == "convert") {
    return run_convert(params.convert);
  } else if (params.command == "render") {
    return run_render(params.render);
  } else if (params.command == "view") {
    return run_view(params.view);
  } else if (params.command == "glview") {
    return run_glview(params.glview);
  } else {
    return print_fatal("unknown command " + params.command);
  }
}
