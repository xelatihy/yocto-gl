//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
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
#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_trace.h>
using namespace yocto;
using namespace std::string_literals;

#include <OpenImageDenoise/oidn.hpp>
#include <iostream>

static bool oidn_image_denoise(image<vec4f>& denoised,
    const image<vec4f>& color, const image<vec4f>& albedo,
    const image<vec4f>& normal, bool hdr, std::string& error,
    progress_callback progress_cb) {
  if (!normal.empty() && albedo.empty()) {
    error =
        "cannot use normal feature image without specifying an albedo feature image.";
    return false;
  }

  // All feature images must be the same size as the color image
  if (!albedo.empty() && albedo.imsize() != color.imsize()) {
    error = "albedo image size doesn't match color image size";
    return false;
  }
  if (!normal.empty() && normal.imsize() != color.imsize()) {
    error = "normal image size doesn't match color image size";
    return false;
  }

  auto [width, height] = color.imsize();
  auto device          = oidn::newDevice();
  device.commit();

  // create a new ray tracing denoise filter
  auto filter = device.newFilter("RT");

  // set the color image of the filter
  filter.set("hdr", hdr);
  filter.setImage(
      "color", (void*)color.data(), oidn::Format::Float4, width, height);

  // set albedo if present
  if (!albedo.empty()) {
    filter.setImage(
        "albedo", (void*)albedo.data(), oidn::Format::Float4, width, height);
  }

  // set normal if present
  if (!normal.empty()) {
    filter.setImage(
        "normal", (void*)normal.data(), oidn::Format::Float4, width, height);
  }

  // initialize 'out' image to the correct size and set it as filter output
  denoised.resize(color.imsize());
  filter.setImage(
      "output", (void*)denoised.data(), oidn::Format::Float4, width, height);

  // register the user provided progress callback in the filter
  if (progress_cb) {
    auto prog_monitor = [](void* uptr, double prog) {
      int  current     = (int)(prog * 100);
      auto progress_cb = *(progress_callback*)uptr;
      progress_cb("denoise image", current, 100);
      return true;
    };
    filter.setProgressMonitorFunction(prog_monitor, (void*)&progress_cb);
  }

  // excecute the filter
  filter.commit();
  filter.execute();

  // check and return eventual oidn errors
  const char* errorMessage;
  if (device.getError(errorMessage) != oidn::Error::None) {
    error = errorMessage;
    return false;
  }

  return true;
}

int main(int argc, const char* argv[]) {
  auto outname     = "out.png"s;
  auto filename    = "img.hdr"s;
  auto albedo_name = ""s;
  auto normal_name = ""s;

  // parse cli arguments and assure their correctness
  auto cli = make_cli(
      "yimagedenoise", "Denoise images using Intel Open Image Denoise");
  add_optional(cli, "alb", albedo_name, "albedo feature image filename", "a");
  add_optional(cli, "nrm", normal_name, "normal feature image filename", "n");
  add_optional(cli, "output", outname, "output image filename", "o");
  add_positional(cli, "filename", filename, "input image filename", true);
  parse_cli(cli, argc, argv);

  if (normal_name != "" && albedo_name == "")
    print_fatal("cannot use normal feature image without an albedo one");
  if (normal_name != "" && !is_hdr_filename(normal_name))
    print_fatal("normal feature image must be provided in pfm or exr format");

  // error string for yocto library functions
  auto error = ""s;

  // load all the provided images
  auto hdr   = is_hdr_filename(filename);
  auto color = image<vec4f>{}, albedo = image<vec4f>{}, normal = image<vec4f>{};

  if (!load_image(filename, color, error)) {
    print_fatal(error);
  }
  if (albedo_name != "" && !load_image(albedo_name, albedo, error)) {
    print_fatal(error);
  }
  if (normal_name != "" && !load_image(normal_name, normal, error)) {
    print_fatal(error);
  }

  // call the denoiser passing in the noisy image and the provided feature
  // images
  auto denoised = image<vec4f>{};
  if (!oidn_image_denoise(
          denoised, color, albedo, normal, hdr, error, print_progress))
    print_fatal(error);

  // finally save the result
  print_progress("save image", 0, 1);
  if (!save_image(outname, denoised, error)) print_fatal(error);
  print_progress("save image", 1, 1);
}
