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
using namespace std;

#include "ext/CLI11.hpp"
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

// progress callback
void print_progress(const string& message, int current, int total) {
  static auto pad = [](const string& str, int n) -> string {
    return string(max(0, n - str.size()), '0') + str;
  };
  static auto pade = [](const string& str, int n) -> string {
    return str + string(max(0, n - str.size()), ' ');
  };
  using clock               = std::chrono::high_resolution_clock;
  static int64_t start_time = 0;
  if (current == 0) start_time = clock::now().time_since_epoch().count();
  auto elapsed = clock::now().time_since_epoch().count() - start_time;
  elapsed /= 1000000;  // millisecs
  auto mins  = pad(to_string(elapsed / 60000), 2);
  auto secs  = pad(to_string((elapsed % 60000) / 1000), 2);
  auto msecs = pad(to_string((elapsed % 60000) % 1000), 3);
  auto n     = (int)(30 * (float)current / (float)total);
  auto bar   = "[" + pade(string(n, '='), 30) + "]";
  auto line  = bar + " " + mins + ":" + secs + "." + msecs + " " +
              pade(message, 30);
  printf("\r%s\r", line.c_str());
  if (current == total) printf("\n");
  fflush(stdout);
}

int run_app(int argc, const char** argv) {
  // command line parameters
  auto mesh_filenames     = false;
  auto shape_directory    = "shapes/"s;
  auto subdiv_directory   = "subdivs/"s;
  auto instance_directory = "instances/"s;
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
  auto scene_guard = make_unique<sceneio_model>();
  auto scene       = scene_guard.get();
  load_scene(filename, scene, print_progress);

  // validate scene
  if (validate) {
    for (auto& error : scene_validation(scene))
      printf("error: %s\n", error.c_str());
  }

  // print info
  if (info) {
    printf("scene stats ------------\n");
    for (auto stat : scene_stats(scene)) printf("%s\n", stat.c_str());
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
  save_scene(output, scene, print_progress);

  // done
  return 0;
}

int main(int argc, const char* argv[]) {
  try {
    return run_app(argc, argv);
  } catch (std::exception& e) {
    fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
