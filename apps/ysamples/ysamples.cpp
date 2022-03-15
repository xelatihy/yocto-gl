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
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_scene.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>

using namespace yocto;
using namespace std::string_literals;

// main function
void run(const vector<string>& args) {
  // parameters
  auto shapename = "shape.ply"s;
  auto outname   = "out.ply"s;
  auto ashairs   = false;
  auto samples   = 4096;
  auto hairs     = 65536;
  auto steps     = 8;
  auto length    = 0.02f;
  auto noise     = 0.001f;
  auto gravity   = 0.0005f;
  auto radius    = 0.0001f;

  // parse command line
  auto cli = make_cli("ysamples", "sample shapes");
  add_option(cli, "shape", shapename, "input shape");
  add_option(cli, "output", outname, "output shape");
  add_option(cli, "ashairs", ashairs, "as hairs");
  add_option(cli, "samples", samples, "number of samples");
  add_option(cli, "hairs", hairs, "number of hairs");
  add_option(cli, "steps", steps, "hair steps");
  add_option(cli, "length", length, "hair length");
  add_option(cli, "noise", noise, "noise weight");
  add_option(cli, "gravity", gravity, "gravity scale");
  add_option(cli, "radius", radius, "hair radius");

  // load mesh
  auto shape = load_shape(shapename);

  if (!ashairs) {
    // generate samples
    auto points = sample_shape(shape, samples);

    // sample shape
    auto sshape = shape_data{};
    for (auto& point : points) {
      sshape.points.push_back((int)sshape.points.size());
      sshape.positions.push_back(eval_position(shape, point.element, point.uv));
      sshape.radius.push_back(radius * 10);
    }

    // save mesh
    save_shape(outname, sshape);
  } else {
    // generate hair
    auto sshape = make_hair2(shape, {steps, hairs}, {length, length},
        {radius, radius}, noise, gravity);

    // save mesh
    save_shape(outname, sshape, true);
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
