//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
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
using namespace yocto;

#include <memory>

#include "ext/filesystem.hpp"
namespace sfs = ghc::filesystem;

// Shape presets used ofr testing.
bool make_preset(scene_model* scene, const string& type, string& error) {
  if (type == "cornellbox") {
    make_cornellbox(scene);
    return true;
  } else {
    error = "unknown preset";
    return false;
  }
  return true;
}

void make_dir(const string& dirname) {
  if (sfs::exists(dirname)) return;
  try {
    sfs::create_directories(dirname);
  } catch (...) {
    print_fatal("cannot create directory " + dirname);
  }
}

int main(int argc, const char* argv[]) {
  // command line parameters
  auto validate  = false;
  auto info      = false;
  auto copyright = ""s;
  auto output    = "out.json"s;
  auto filename  = "scene.json"s;

  // parse command line
  auto cli = make_cli("yscnproc", "Process scene");
  add_option(cli, "--info,-i", info, "print scene info");
  add_option(cli, "--copyright,-c", copyright, "copyright string");
  add_option(cli, "--validate/--no-validate", validate, "Validate scene");
  add_option(cli, "--output,-o", output, "output scene");
  add_option(cli, "scene", filename, "input scene", true);
  parse_cli(cli, argc, argv);

  // load scene
  auto ext         = sfs::path(filename).extension().string();
  auto basename    = sfs::path(filename).stem().string();
  auto scene_guard = std::make_unique<scene_model>();
  auto scene       = scene_guard.get();
  auto ioerror     = ""s;
  if (ext == ".ypreset") {
    print_progress("make preset", 0, 1);
    if (!make_preset(scene, basename, ioerror)) print_fatal(ioerror);
    print_progress("make preset", 1, 1);
  } else {
    if (!load_scene(filename, scene, ioerror, print_progress))
      print_fatal(ioerror);
  }

  // copyright
  if (copyright != "") {
    scene->copyright = copyright;
  }

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
  if (sfs::path(output).extension() != ".json") {
    tesselate_shapes(scene, print_progress);
  }

  // make a directory if needed
  make_dir(sfs::path(output).parent_path());
  if (!scene->shapes.empty())
    make_dir(sfs::path(output).parent_path() / "shapes");
  if (!scene->textures.empty())
    make_dir(sfs::path(output).parent_path() / "textures");
  if (!scene->instances.empty())
    make_dir(sfs::path(output).parent_path() / "instances");

  // save scene
  if (!save_scene(output, scene, ioerror, print_progress)) print_fatal(ioerror);

  // done
  return 0;
}
