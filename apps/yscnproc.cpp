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

#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_sceneio.h"
using namespace yocto;

#include <unordered_set>
using std::unordered_set;

#include "ext/filesystem.hpp"
namespace fs = ghc::filesystem;

#include "ext/CLI11.hpp"

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

int main(int argc, char** argv) {
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
  auto parser = CLI::App{"Process scene"};
  parser.add_flag("--notextures", notextures, "Disable textures.");
  parser.add_flag("--mesh-filenames", mesh_filenames, "Add mesh filenames.");
  parser.add_option("--shape-directory", shape_directory,
      "Shape directory when adding names.");
  parser.add_option("--subdiv-directory", subdiv_directory,
      "Subdiv directory when adding names.");
  parser.add_flag("--uniform-textures", uniform_txt, "uniform texture formats");
  parser.add_flag(
      "--obj-instances", obj_instances, "preserve instances in obj");
  parser.add_flag("--info,-i", info, "print scene info");
  parser.add_flag("--validate", validate, "Validate scene");
  parser.add_option("--output,-o", output, "output scene")->required(true);
  parser.add_option("scene", filename, "input scene")->required(true);
  try {
    parser.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return parser.exit(e);
  }
  setbuf(stdout, nullptr);

  // fix options
  auto load_prms         = load_params();
  auto save_prms         = save_params();
  load_prms.notextures   = notextures;
  save_prms.notextures   = notextures;
  save_prms.objinstances = obj_instances;

  // load scene
  auto scene = yocto_scene{};
  try {
    printf("loading scene");
    auto load_timer = timer();
    load_scene(filename, scene, load_prms);
    printf(" in %s\n", load_timer.elapsedf().c_str());
  } catch (const std::exception& e) {
    printf("%s\n", e.what());
    exit(1);
  }

  // validate scene
  if (validate) {
    printf("validating scene");
    auto validate_timer = timer();
    print_validation(scene);
    printf(" in %s\n", validate_timer.elapsedf().c_str());
  }

  // print info
  if (info) printf("%s\n", format_stats(scene).c_str());

  // change texture names
  if (uniform_txt) {
    for (auto& texture : scene.textures) {
      auto ext = fs::path(texture.uri).extension().string();
      if (is_hdr_filename(texture.uri)) {
        texture.uri = fs::path(texture.uri).replace_extension(".hdr");
      } else {
        texture.uri = fs::path(texture.uri).replace_extension(".png");
      }
    }
  }

  // tesselating scene
  printf("tesselating scene");
  auto tesselate_timer = timer();
  tesselate_subdivs(scene);
  printf(" in %s\n", tesselate_timer.elapsedf().c_str());

  // add missing mesh names if necessary
  if (!shape_directory.empty() && shape_directory.back() != '/')
    shape_directory += '/';
  if (mesh_filenames) {
    auto sid = 0;
    for (auto& shape : scene.shapes) {
      if (!shape.quadspos.empty()) {
        shape.uri = shape_directory + "shape_" + std::to_string(sid) + ".obj";
      } else {
        shape.uri = shape_directory + "shape_" + std::to_string(sid) + ".ply";
      }
      sid++;
    }
  }
  // add missing mesh names if necessary
  if (!subdiv_directory.empty() && subdiv_directory.back() != '/')
    subdiv_directory += '/';
  if (mesh_filenames) {
    auto sid = 0;
    for (auto& subdiv : scene.subdivs) {
      if (!subdiv.quadspos.empty()) {
        subdiv.uri = subdiv_directory + "subdiv_" + std::to_string(sid) +
                     ".obj";
      } else {
        subdiv.uri = subdiv_directory + "subdiv_" + std::to_string(sid) +
                     ".ply";
      }
      sid++;
    }
  }

  // make a directory if needed
  auto dirname  = fs::path(output).parent_path();
  auto dirnames = unordered_set<string>{};
  if (!dirname.empty()) dirnames.insert(dirname);
  for (auto& shape : scene.shapes)
    dirnames.insert(dirname / fs::path(shape.uri).parent_path());
  for (auto& subdiv : scene.subdivs)
    dirnames.insert(dirname / fs::path(subdiv.uri).parent_path());
  for (auto& texture : scene.textures)
    dirnames.insert(dirname / fs::path(texture.uri).parent_path());
  for (auto& dir : dirnames) {
    if (!mkdir(fs::path(dir))) {
      printf("cannot create directory %s\n", fs::path(output).c_str());
    }
  }

  // save scene
  try {
    printf("saving scene");
    auto save_timer = timer();
    save_scene(output, scene, save_prms);
    printf(" in %s\n", save_timer.elapsedf().c_str());
  } catch (const std::exception& e) {
    printf("%s\n", e.what());
    exit(1);
  }

  // done
  return 0;
}
