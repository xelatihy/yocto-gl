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

#include <yocto/yocto_cli.h>
#include <yocto/yocto_color.h>
#include <yocto/yocto_gui.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_scene.h>
#include <yocto/yocto_sceneio.h>

using namespace yocto;

// grade params
struct app_params : colorgrade_params {
  string image       = "image.png";
  string output      = "out.png";
  bool   interactive = true;
};

// Cli
app_params parse_params(int argc, const char* argv[]) {
  auto params = app_params{};

  auto cli = make_cli("ycolorgrade", "adjust image colors");
  add_option(cli, "image", params.image, "Input image.");
  add_option(cli, "output", params.output, "Output image.");
  parse_cli(cli, argc, argv);

  return params;
}

// grade images
void colorgrade_offline(const app_params& params) {
  // load image
  auto image = load_image(params.image);

  // apply color grade
  image = colorgrade_image(image, params);

  // save image
  save_image(params.output, image);
}

// grade images
void colorgrade_interactive(const app_params& params) {
  // load image
  auto image = load_image(params.image);

  // run viewer
  show_colorgrade_gui("ycolorgrade", params.image, image);
}

// Main
int main(int argc, const char* argv[]) {
  try {
    // command line parameters
    auto params = parse_params(argc, argv);

    // dispatch commands
    if (!params.interactive) {
      colorgrade_offline(params);
    } else {
      colorgrade_interactive(params);
    }
  } catch (const std::exception& error) {
    print_error(error.what());
    return 1;
  }

  // done
  return 0;
}
