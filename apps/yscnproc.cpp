//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#include <memory>

#include "../yocto/yocto_commonio.h"
#include "../yocto/yocto_image.h"
#include "../yocto/yocto_math.h"
#include "../yocto/yocto_sceneio.h"
using std::string;
using namespace std::string_literals;

#include "ext/filesystem.hpp"
namespace fs = ghc::filesystem;

void make_dir(const std::string& dirname) {
  if (fs::exists(dirname)) return;
  try {
    fs::create_directories(dirname);
  } catch (...) {
    ycl::print_fatal("cannot create directory " + dirname);
  }
}

int main(int argc, const char* argv[]) {
  // command line parameters
  auto validate = false;
  auto info     = false;
  auto output   = "out.json"s;
  auto filename = "scene.json"s;

  // parse command line
  auto cli = ycl::make_cli("yscnproc", "Process scene");
  add_option(cli, "--info,-i", info, "print scene info");
  add_option(cli, "--validate/--no-validate", validate, "Validate scene");
  add_option(cli, "--output,-o", output, "output scene");
  add_option(cli, "scene", filename, "input scene", true);
  parse_cli(cli, argc, argv);

  // load scene
  auto scene_guard = std::make_unique<ysc::model>();
  auto scene       = scene_guard.get();
  auto ioerror     = ""s;
  if (!load_scene(filename, scene, ioerror, ycl::print_progress))
    ycl::print_fatal(ioerror);

  // validate scene
  if (validate) {
    for (auto& error : scene_validation(scene))
      ycl::print_info("error: " + error);
  }

  // print info
  if (info) {
    ycl::print_info("scene stats ------------");
    for (auto stat : scene_stats(scene)) ycl::print_info(stat);
  }

  // tesselate if needed
  if (fs::path(output).extension() != ".json") {
    for (auto iosubdiv : scene->subdivs) {
      tesselate_subdiv(scene, iosubdiv);
    }
  }

  // make a directory if needed
  make_dir(fs::path(output).parent_path());
  if (!scene->shapes.empty())
    make_dir(fs::path(output).parent_path() / "shapes");
  if (!scene->subdivs.empty())
    make_dir(fs::path(output).parent_path() / "subdivs");
  if (!scene->textures.empty())
    make_dir(fs::path(output).parent_path() / "textures");
  if (!scene->instances.empty())
    make_dir(fs::path(output).parent_path() / "instances");

  // save scene
  if (!save_scene(output, scene, ioerror, ycl::print_progress))
    ycl::print_fatal(ioerror);

  // done
  return 0;
}
