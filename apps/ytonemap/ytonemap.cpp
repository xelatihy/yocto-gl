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
#include <yocto/yocto_sceneio.h>

using namespace yocto;
using namespace std::string_literals;

// main function
void run(const vector<string>& args) {
  // parameters
  auto filename    = "image.png"s;
  auto output      = "out.png"s;
  auto exposure    = 0.0f;
  auto filmic      = false;
  auto width       = 0;
  auto height      = 0;
  auto interactive = false;

  // parse command line
  auto cli = make_cli("ytonemap", "tonemap image");
  add_option(cli, "image", filename, "Input image.");
  add_option(cli, "output", output, "Output image.");
  add_option(cli, "exposure", exposure, "Tonemap exposure.");
  add_option(cli, "filmic", filmic, "Tonemap filmic.");
  add_option(cli, "width", width, "Resize width.");
  add_option(cli, "height", height, "Resize height.");
  add_option(cli, "interactive", interactive, "Run interactively.");
  parse_cli(cli, args);

  // load
  auto image = load_image(filename);

  // resize if needed
  if (width != 0 || height != 0) {
    image = resize_image(image, width, height);
  }

  // switch between interactive and offline
  if (!interactive) {
    // tonemap if needed
    if (image.linear && is_ldr_filename(output)) {
      image = tonemap_image(image, exposure, filmic);
    }

    // save
    save_image(output, image);
  } else {
    // run viewer
    show_image_gui("ytonemap", filename, image);
  }
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
