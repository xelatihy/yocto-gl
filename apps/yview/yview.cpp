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
#include <yocto/yocto_gui.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_scene.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
#include <yocto/yocto_trace.h>

using namespace yocto;

#include <filesystem>
namespace fs = std::filesystem;

struct app_params {
  string scene   = "scene.json";
  string camname = "";
};

// Cli
app_params parse_params(int argc, const char* argv[]) {
  auto params = app_params{};

  auto cli = make_cli("ytrace", "render with raytracing");
  add_option(cli, "scene", params.scene, "input scene");
  add_option(cli, "camera", params.camname, "camera name");
  parse_cli(cli, argc, argv);

  return params;
}

void view_interactively(const app_params& params_) {
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

// Run
int main(int argc, const char* argv[]) {
  try {
    // command line parameters
    auto params = parse_params(argc, argv);

    // dispatch commands
    if (!params.scene.empty()) {
      view_interactively(params);
    } else {
      throw io_error{"missing scene"};
    }
  } catch (const std::exception& error) {
    print_error(error.what());
    return 1;
  }

  // done
  return 0;
}
