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
  auto filename1 = "image1.png"s;
  auto filename2 = "image2.png"s;
  auto outname   = ""s;
  auto signal    = false;
  auto threshold = 0.0f;

  auto cli = make_cli("yimdiff", "tonemap image");
  add_option(cli, "image1", filename1, "Input image 1.");
  add_option(cli, "image2", filename2, "Input image 2.");
  add_option(cli, "output", outname, "Output image.");
  add_option(cli, "signal", signal, "Error on diff.");
  add_option(cli, "threshold", threshold, "Diff threshold.");
  parse_cli(cli, args);

  // load
  auto image1 = load_image(filename1);
  auto image2 = load_image(filename2);

  // check sizes
  if (image1.width != image2.width || image1.height != image2.height)
    throw io_error("different image sizes");

  // check types
  if (image1.linear != image2.linear) throw io_error("different image types");

  // compute diff
  auto diff = image_difference(image1, image2, true);

  // save
  if (outname != "") save_image(outname, diff);

  // check diff
  if (signal) {
    for (auto& c : diff.pixels) {
      if (max(xyz(c)) > threshold) {
        throw std::runtime_error{"image content differ"};
      }
    }
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
