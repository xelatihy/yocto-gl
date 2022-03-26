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
#include <yocto/yocto_color.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_sceneio.h>

using namespace yocto;
using namespace std::string_literals;

// main function
void run(const vector<string>& args) {
  // parameters
  auto imagename  = "image.png"s;
  auto alphaname  = "alpha.png"s;
  auto outname    = "out.png"s;
  auto from_color = false;
  auto from_black = false;
  auto to_color   = false;

  // parse command line
  auto cli = make_cli("yimalpha", "set image alpha");
  add_option(cli, "image", imagename, "Input image");
  add_option(cli, "alpha", alphaname, "Alpha image");
  add_option(cli, "output", outname, "Output image");
  add_option(cli, "from-color", from_color, "Alpha from color");
  add_option(cli, "from-black", from_black, "Alpha from black");
  add_option(cli, "to-color", to_color, "color from alpha");
  parse_cli(cli, args);

  // load
  auto image = load_image(imagename);
  auto alpha = load_image(alphaname);

  // check sizes
  if (image.width != alpha.width || image.height != alpha.height)
    throw io_error("different image sizes");

  // check types
  if (image.linear != alpha.linear) throw io_error("different image types");

  // edit alpha
  auto out = make_image(image.width, image.height, image.linear);
  for (auto idx : range(image.pixels.size())) {
    auto calpha = alpha.pixels[idx];
    auto alpha_ = from_color   ? mean(xyz(calpha))
                  : from_black ? (mean(xyz(calpha)) > 0.01 ? 1.0f : 0.0f)
                               : calpha.w;
    if (to_color) {
      out.pixels[idx] = {alpha_, alpha_, alpha_, alpha_};
    } else {
      auto color      = image.pixels[idx];
      color.w         = alpha_;
      out.pixels[idx] = color;
    }
  }

  // save
  save_image(outname, out);
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
