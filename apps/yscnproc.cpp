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
#include "../yocto/yocto_utils.h"
using namespace yocto;

#include "ext/CLI11.hpp"

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

int main(int argc, char** argv) {
  // command line parameters
  auto notextures       = false;
  auto mesh_filenames   = false;
  auto shape_directory  = "shapes/"s;
  auto subdiv_directory = "subdivs/"s;
  auto uniform_txt      = false;
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
  parser.add_flag("--info,-i", info, "print scene info");
  parser.add_flag("--validate", validate, "Validate scene");
  parser.add_option("--output,-o", output, "output scene")->required(true);
  parser.add_option("scene", filename, "input scene")->required(true);
  try {
    parser.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return parser.exit(e);
  }

  // fix options
  auto load_prms       = load_params();
  auto save_prms       = save_params();
  load_prms.notextures = notextures;
  save_prms.notextures = notextures;

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
  if (info) print_info(format_stats(scene));

  // change texture names
  if (uniform_txt) {
    for (auto& texture : scene.textures) {
      auto ext = get_extension(texture.uri);
      if (is_hdr_filename(texture.uri)) {
        texture.uri = get_noextension(texture.uri) + ".hdr";
      } else {
        texture.uri = get_noextension(texture.uri) + ".png";
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
        shape.uri = shape_directory + "shape_" + to_string(sid) + ".obj";
      } else {
        shape.uri = shape_directory + "shape_" + to_string(sid) + ".ply";
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
        subdiv.uri = subdiv_directory + "subdiv_" + to_string(sid) + ".obj";
      } else {
        subdiv.uri = subdiv_directory + "subdiv_" + to_string(sid) + ".ply";
      }
      sid++;
    }
  }

  // make a directory if needed
  auto dirname  = get_dirname(output);
  auto dirnames = unordered_set<string>{dirname};
  for (auto& shape : scene.shapes)
    dirnames.insert(dirname + get_dirname(shape.uri));
  for (auto& subdiv : scene.subdivs)
    dirnames.insert(dirname + get_dirname(subdiv.uri));
  for (auto& texture : scene.textures)
    dirnames.insert(dirname + get_dirname(texture.uri));
  for (auto& dir : dirnames) {
    if (!mkdir(get_dirname(dir))) {
      print_fatal("cannot create directory " + get_dirname(output));
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
