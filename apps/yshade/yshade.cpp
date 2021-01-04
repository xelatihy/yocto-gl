//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
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
#include <yocto/yocto_json.h>

#include "yshade_scene.h"
#include "yshade_sculpt.h"
#include "yshade_shape.h"

using namespace yocto;

struct app_params {
  string         command = "sculpt";
  sculpt__params sculpt  = {};
  scene_params   scene   = {};
  shape_params   shape   = {};
};

// Json IO
void serialize_value(json_mode mode, json_value& json, app_params& value,
    const string& description) {
  serialize_object(mode, json, value, description);
  serialize_command(mode, json, value.command, "command", "Command.");
  serialize_property(mode, json, value.scene, "scene", "View scenes.");
  serialize_property(mode, json, value.shape, "shape", "View shapes.");
  serialize_property(mode, json, value.sculpt, "sculpt", "Sculpt shapes.");
}

int main(int argc, const char* argv[]) {
  // command line parameters
  auto params = app_params{};
  parse_cli(params, "View and edit interactively", argc, argv);

  // dispatch commands
  if (params.command == "scene") {
    return run_scene(params.scene);
  } else if (params.command == "shape") {
    return run_shape(params.shape);
  } else if (params.command == "sculpt") {
    return run_sculpt(params.sculpt);
  } else {
    return print_fatal("unknown command " + params.command);
  }
}
