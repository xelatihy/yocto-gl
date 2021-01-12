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
#include <yocto/yocto_math.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
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
void add_command(cli_command& cli, const string& name, convert_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_positional(cmd, "scene", value.scene, "Input scene.");
  add_optional(cmd, "output", value.output, "Output scene.", "o");
  add_optional(cmd, "info", value.info, "Print info.");
  add_optional(cmd, "validate", value.validate, "Validate scene.");
  add_optional(cmd, "copyright", value.copyright, "Set scene copyright.");
}

// convert images
int run_convert(const convert_params& params) {
  // load scene
  auto scene   = scene_scene{};
  auto ioerror = ""s;
  if (!load_scene(params.scene, scene, ioerror, print_progress))
    print_fatal(ioerror);

  // copyright
  if (params.copyright != "") {
    scene.asset.copyright = params.copyright;
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
  if (path_extension(params.output) != ".json") {
    tesselate_shapes(scene, print_progress);
  }

  // make a directory if needed
  if (!make_directory(path_dirname(params.output), ioerror))
    print_fatal(ioerror);
  if (!scene.shapes.empty()) {
    if (!make_directory(
            path_join(path_dirname(params.output), "shapes"), ioerror))
      print_fatal(ioerror);
  }
  if (!scene.textures.empty()) {
    if (!make_directory(
            path_join(path_dirname(params.output), "textures"), ioerror))
      print_fatal(ioerror);
  }

  // save scene
  if (!save_scene(params.output, scene, ioerror, print_progress))
    print_fatal(ioerror);

  // done
  return 0;
}

// convert params
struct view_params {
  string scene  = "scene.json";
  string output = "out.png";
  string camera = "";
  bool   addsky = false;
};

// Cli
void add_command(cli_command& cli, const string& name, view_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_positional(cmd, "scene", value.scene, "Scene filename.");
  add_optional(cmd, "output", value.output, "Output filename.", "o");
  add_optional(cmd, "camera", value.camera, "Camera name.");
  add_optional(cmd, "addsky", value.addsky, "Add sky.");
}

#ifndef YOCTO_OPENGL

// view scene
int run_view(const view_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

// view scene
int run_view(const view_params& params) {
  // load scene
  auto scene   = scene_scene{};
  auto ioerror = ""s;
  if (!load_scene(params.scene, scene, ioerror, print_progress))
    print_fatal(ioerror);

  // add sky
  if (params.addsky) add_sky(scene);

  // tesselation
  tesselate_shapes(scene, print_progress);

  // run view
  view_scene("yscene", params.scene, scene, find_camera(scene, params.camera),
      print_progress);

  // done
  return 0;
}

#endif

struct app_params {
  string         command = "convert";
  convert_params convert = {};
  view_params    view    = {};
};

// Cli
void add_commands(cli_command& cli, const string& name, app_params& value,
    const string& usage) {
  cli = make_cli(name, usage);
  add_command_name(cli, "command", value.command, "Command.");
  add_command(cli, "convert", value.convert, "Convert shapes.");
  add_command(cli, "view", value.view, "View shapes.");
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

  // dispatch commands
  if (params.command == "convert") {
    return run_convert(params.convert);
  } else if (params.command == "view") {
    return run_view(params.view);
  } else {
    return print_fatal("unknown command " + params.command);
  }
}
