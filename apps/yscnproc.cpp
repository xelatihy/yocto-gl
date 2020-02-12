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

#include "../yocto/yocto_image.h"
#include "../yocto/yocto_math.h"
#include "../yocto/yocto_sceneio.h"
using namespace yocto;

#include <memory>
#include <set>
using std::make_unique;
using std::set;

#include "ext/CLI11.hpp"
#include "ext/Timer.hpp"
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

int run_app(int argc, const char** argv) {
  // command line parameters
  auto mesh_filenames     = false;
  auto shape_directory    = "shapes/"s;
  auto subdiv_directory   = "subdivs/"s;
  auto instance_directory = "instances/"s;
  auto obj_instances      = false;
  auto validate           = false;
  auto info               = false;
  auto output             = "out.json"s;
  auto filename           = "scene->json"s;

  // parse command line
  auto cli = CLI::App{"Process scene"};
  cli.add_option("--mesh-filenames", mesh_filenames, "Add mesh filenames.");
  cli.add_option("--shape-directory", shape_directory,
      "Shape directory when adding names.");
  cli.add_option("--subdiv-directory", subdiv_directory,
      "Subdiv directory when adding names.");
  cli.add_flag("--obj-instances", obj_instances, "preserve instances in obj");
  cli.add_option("--info,-i", info, "print scene info");
  cli.add_flag("--validate", validate, "Validate scene");
  cli.add_option("--output,-o", output, "output scene", true);
  cli.add_option("scene", filename, "input scene", true);
  try {
    cli.parse(argc, argv);
  } catch (CLI::ParseError& e) {
    return cli.exit(e);
  }

  // load scene
  auto scene = shared_ptr<sceneio_model>{};
  {
    auto timer = CLI::AutoTimer("loading scene");
    scene      = load_scene(filename);
  }

  // validate scene
  if (validate) {
    auto timer = CLI::AutoTimer("validate");
    for (auto& error : scene_validation(scene)) std::cout << error << "\n";
  }

  // print info
  if (info) {
    std::cout << "scene stats ------------\n";
    for (auto stat : scene_stats(scene)) std::cout << stat << "\n";
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
  auto dirnames = set<fs::path>{};
  if (!dirname.empty()) dirnames.insert(dirname);
  for (auto& shape : scene->shapes)
    dirnames.insert(dirname / fs::path(shape->name).parent_path());
  for (auto texture : scene->textures)
    dirnames.insert(dirname / fs::path(texture->name).parent_path());
  for (auto instance : scene->instances)
    dirnames.insert(dirname / fs::path(instance->name).parent_path());
  for (auto& dir : dirnames) {
    if (!mkdir(dir)) {
      throw std::runtime_error{"cannot create directory " + output};
    }
  }

  // save scene
  {
    auto timer = CLI::AutoTimer("save");
    save_scene(output, scene, obj_instances);
  }

  // done
  return 0;
}

int main(int argc, const char* argv[]) {
  try {
    return run_app(argc, argv);
  } catch (std::exception& e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
}
