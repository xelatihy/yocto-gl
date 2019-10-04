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
#include "../yocto/yocto_math.h"
#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_sceneio.h"
using namespace yocto;

#include <unordered_set>
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

int main(int argc, const char** argv) {
  // command line parameters
  auto notextures       = false;
  auto mesh_filenames   = false;
  auto shape_directory  = "shapes/"s;
  auto subdiv_directory = "subdivs/"s;
  auto uniform_txt      = false;
  auto obj_instances    = false;
  auto validate         = false;
  auto info             = false;
  auto output           = "out.json"s;
  auto filename         = "scene.json"s;

  // parse command line
  auto cli = make_cli("yscnproc", "Process scene");
  add_cli_option(cli, "--notextures", notextures, "Disable textures.");
  add_cli_option(
      cli, "--mesh-filenames", mesh_filenames, "Add mesh filenames.");
  add_cli_option(cli, "--shape-directory", shape_directory,
      "Shape directory when adding names.");
  add_cli_option(cli, "--subdiv-directory", subdiv_directory,
      "Subdiv directory when adding names.");
  add_cli_option(
      cli, "--uniform-textures", uniform_txt, "uniform texture formats");
  add_cli_option(
      cli, "--obj-instances", obj_instances, "preserve instances in obj");
  add_cli_option(cli, "--info,-i", info, "print scene info");
  add_cli_option(cli, "--validate", validate, "Validate scene");
  add_cli_option(cli, "--output,-o", output, "output scene", true);
  add_cli_option(cli, "scene", filename, "input scene", true);
  if (!parse_cli(cli, argc, argv)) exit(1);

  // fix options
  auto load_prms         = load_params();
  auto save_prms         = save_params();
  load_prms.notextures   = notextures;
  save_prms.notextures   = notextures;
  save_prms.objinstances = obj_instances;

  // load scene
  auto scene = yocto_scene{};
  try {
    auto timer = print_timed("loading scene");
    load_scene(filename, scene, load_prms);
  } catch (const std::exception& e) {
    print_fatal(e.what());
  }

  // validate scene
  if (validate) {
    auto timer = print_timed("validating scene");
    print_validation(scene);
  }

  // print info
  if (info) {
    for (auto stat : format_stats(scene)) print_info(stat);
  }

  // change texture names
  if (uniform_txt) {
    for (auto& texture : scene.textures) {
      auto ext = get_extension(texture.filename);
      if (is_hdr_filename(texture.filename)) {
        texture.filename = replace_extension(texture.filename, ".hdr");
      } else {
        texture.filename = replace_extension(texture.filename, ".png");
      }
    }
  }

  // tesselating scene
  {
    auto timer = print_timed("tesselating scene");
    tesselate_subdivs(scene);
  }

  // add missing mesh names if necessary
  if (!shape_directory.empty() && shape_directory.back() != '/')
    shape_directory += '/';
  if (mesh_filenames) {
    auto sid = 0;
    for (auto& shape : scene.shapes) {
      if (!shape.quadspos.empty()) {
        shape.filename = shape_directory + "shape_" + std::to_string(sid) +
                         ".obj";
      } else {
        shape.filename = shape_directory + "shape_" + std::to_string(sid) +
                         ".ply";
      }
      sid++;
    }
  }

  // make a directory if needed
  auto dirname  = get_dirname(output);
  auto dirnames = unordered_set<string>{};
  if (!dirname.empty()) dirnames.insert(dirname);
  for (auto& shape : scene.shapes)
    dirnames.insert(dirname + get_dirname(shape.filename));
  for (auto& texture : scene.textures)
    dirnames.insert(dirname + get_dirname(texture.filename));
  for (auto& dir : dirnames) {
    if (!mkdir(dir)) {
      print_fatal("cannot create directory " + output);
    }
  }

  // save scene
  try {
    auto timer = print_timed("saving scene");
    save_scene(output, scene, save_prms);
  } catch (const std::exception& e) {
    print_fatal(e.what());
  }

  // done
  return 0;
}
