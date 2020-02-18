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

#include "../yocto/yocto_commonio.h"
#include "../yocto/yocto_image.h"
#include "../yocto/yocto_math.h"
#include "../yocto/yocto_sceneio.h"
using namespace yocto;

#include <memory>
#include <set>
using std::make_unique;

#include "ext/filesystem.hpp"
namespace fs = ghc::filesystem;

bool mkdir(const string& dir) {
  if (dir == "" || dir == "." || dir == ".." || dir == "./" || dir == "../")
    return true;
#ifndef _MSC_VER
  system(("mkdir -p " + dir).c_str());
  return true;
#else
  system(("mkdir " + dir).c_str());
  return true;
#endif
}

int main(int argc, const char* argv[]) {
  // command line parameters
  auto mesh_filenames     = false;
  auto shape_directory    = "shapes/"s;
  auto subdiv_directory   = "subdivs/"s;
  auto instance_directory = "instances/"s;
  auto validate           = false;
  auto info               = false;
  auto output             = "out.json"s;
  auto filename           = "scene.json"s;

  // parse command line
  auto cli = make_cli("yscnproc", "Process scene");
  add_option(cli, "--mesh-filenames", mesh_filenames, "Add mesh filenames.");
  add_option(cli, "--shape-directory", shape_directory,
      "Shape directory when adding names.");
  add_option(cli, "--subdiv-directory", subdiv_directory,
      "Subdiv directory when adding names.");
  add_option(cli, "--info,-i", info, "print scene info");
  add_option(cli, "--validate/--no-validate", validate, "Validate scene");
  add_option(cli, "--output,-o", output, "output scene");
  add_option(cli, "scene", filename, "input scene", true);
  parse_cli(cli, argc, argv);

  // load scene
  auto scene_guard = make_unique<sceneio_model>();
  auto scene       = scene_guard.get();
  auto ioerror     = ""s;
  if (!load_scene(filename, scene, ioerror, print_progress))
    print_fatal(ioerror);

  // validate scene
  if (validate) {
    for (auto& error : scene_validation(scene)) print_info("error: " + error);
  }

  // print info
  if (info) {
    print_info("scene stats ------------");
    for (auto stat : scene_stats(scene)) print_info(stat);
  }

  // tesselate if needed
  if (fs::path(output).extension() != ".json") {
    for (auto iosubdiv : scene->subdivs) {
      tesselate_subdiv(scene, iosubdiv);
    }
  }

  // add missing mesh names if necessary
  if (!shape_directory.empty() && shape_directory.back() != '/')
    shape_directory += '/';

  // add missing mesh names if necessary
  if (!instance_directory.empty() && instance_directory.back() != '/')
    instance_directory += '/';

  // make a directory if needed
  auto dirname  = fs::path(output).parent_path();
  auto dirnames = std::set<fs::path>{};
  if (!dirname.empty()) dirnames.insert(dirname);
  if (!scene->shapes.empty()) dirnames.insert(dirname / "shapes");
  if (!scene->subdivs.empty()) dirnames.insert(dirname / "subdivs");
  if (!scene->textures.empty()) dirnames.insert(dirname / "textures");
  if (!scene->instances.empty()) dirnames.insert(dirname / "instances");
  for (auto& dir : dirnames) {
    if (!mkdir(dir)) {
      throw std::runtime_error{"cannot create directory " + output};
    }
  }

  // save scene
  if (!save_scene(output, scene, ioerror, print_progress)) print_fatal(ioerror);

  // done
  return 0;
}
