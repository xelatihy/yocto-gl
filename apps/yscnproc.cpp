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
#include <unordered_set>
using std::make_shared;
using std::unordered_set;

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

void run_app(int argc, const char** argv) {
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
  auto cli = make_cli("yscnproc", "Process scene");
  add_cli_option(
      cli, "--mesh-filenames", mesh_filenames, "Add mesh filenames.");
  add_cli_option(cli, "--shape-directory", shape_directory,
      "Shape directory when adding names.");
  add_cli_option(cli, "--subdiv-directory", subdiv_directory,
      "Subdiv directory when adding names.");
  add_cli_option(
      cli, "--obj-instances", obj_instances, "preserve instances in obj");
  add_cli_option(cli, "--info,-i", info, "print scene info");
  add_cli_option(cli, "--validate", validate, "Validate scene");
  add_cli_option(cli, "--output,-o", output, "output scene", true);
  add_cli_option(cli, "scene", filename, "input scene", true);
  parse_cli(cli, argc, argv);

  // load scene
  auto scene      = make_shared<sceneio_model>();
  auto load_timer = print_timed("loading scene");
  load_scene(filename, scene.get());
  print_elapsed(load_timer);

  // validate scene
  if (validate) {
    auto validate_timer = print_timed("validating scene");
    auto errors         = scene_validation(scene.get());
    print_elapsed(validate_timer);
    for (auto& error : errors) print_info(error);
  }

  // print info
  if (info) {
    print_info("scene stats ------------");
    for (auto stat : scene_stats(scene.get())) print_info(stat);
  }

  // tesselate if needed
  if (get_extension(output) != ".yaml" && get_extension(output) != ".json") {
    for (auto& iosubdiv : scene->subdivs) {
      tesselate_subdiv(scene.get(), iosubdiv);
      iosubdiv = {};
    }
  }

  // add missing mesh names if necessary
  if (!shape_directory.empty() && shape_directory.back() != '/')
    shape_directory += '/';

  // add missing mesh names if necessary
  if (!instance_directory.empty() && instance_directory.back() != '/')
    instance_directory += '/';

  // make a directory if needed
  auto dirname  = get_dirname(output);
  auto dirnames = unordered_set<string>{};
  if (!dirname.empty()) dirnames.insert(dirname);
  for (auto& shape : scene->shapes)
    dirnames.insert(dirname + get_dirname(shape->name));
  for (auto texture : scene->textures)
    dirnames.insert(dirname + get_dirname(texture->name));
  for (auto instance : scene->instances)
    dirnames.insert(dirname + get_dirname(instance->name));
  for (auto& dir : dirnames) {
    if (!mkdir(dir)) {
      print_fatal("cannot create directory " + output);
    }
  }

  // save scene
  auto save_timer = print_timed("saving scene");
  save_scene(output, scene.get(), obj_instances);
  print_elapsed(save_timer);
}

int main(int argc, const char* argv[]) {
  try {
    run_app(argc, argv);
    return 0;
  } catch (std::exception& e) {
    print_fatal(e.what());
    return 1;
  }
}
