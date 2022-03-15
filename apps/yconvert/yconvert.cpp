//
// LICENSE:
//
// Copyright (c) 2016 -- 2022 Fabio Pellacini
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
using namespace std::string_literals;

// main function
void run(const vector<string>& args) {
  // parameters
  auto scenename = "scene.json"s;
  auto outname   = "out.json"s;
  auto info      = false;
  auto validate  = false;
  auto copyright = ""s;

  // parse command line
  auto cli = make_cli("yconvert", "convert scenes, shapes and images");
  add_option(cli, "scene", scenename, "input scene");
  add_option(cli, "output", outname, "output scenename");
  add_option(cli, "info", info, "print info");
  add_option(cli, "validate", validate, "run validate");
  add_option(cli, "copyright", copyright, "set scene copyright");
  parse_cli(cli, args);

  // start converting
  print_info("converting {}", scenename);
  auto timer = simple_timer{};

  // load scene
  timer      = simple_timer{};
  auto scene = load_scene(scenename);
  print_info("load scene: {}", elapsed_formatted(timer));

  // copyright
  if (copyright != "") {
    scene.copyright = copyright;
  }

  // validate scene
  if (validate) {
    auto errors = scene_validation(scene);
    for (auto& error : errors) print_error(error);
    if (!errors.empty()) throw io_error{"invalid scene"};
  }

  // print info
  if (info) {
    print_info("scene stats ------------");
    for (auto stat : scene_stats(scene)) print_info(stat);
  }

  // tesselate if needed
  if (!scene.subdivs.empty()) {
    timer = simple_timer{};
    tesselate_subdivs(scene);
    print_info("tesselate subdivs: {}", elapsed_formatted(timer));
  }

  // save scene
  timer = simple_timer{};
  make_scene_directories(outname, scene);
  save_scene(outname, scene);
  print_info("save scene: {}", elapsed_formatted(timer));
}

// Main
int main(int argc, const char* argv[]) {
  try {
    run({argv, argv + argc});
    return 0;
  } catch (const std::exception& error) {
    print_error(error.what());
    return 1;
  }
}
