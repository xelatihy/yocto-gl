//
// Library for image denoising using Intel Open Image Denoise
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//

#include "yocto_denoise.h"

#include <OpenImageDenoise/oidn.hpp>

namespace yocto::denoise {

static bool oidn_image_denoise(const image<vec3f> &color, bool hdr,
    const image<vec3f> *albedo, const image<vec3f> *normal, image<vec3f> &out,
    std::string &error, progress_callback progress_cb) {
  if (normal && !albedo)
    throw std::runtime_error{
        "cannot use normal feature image without specifying an albedo feature image."};

  // All feature images must be the same size as the color image
  if (albedo && albedo->imsize() != color.imsize()) {
    error = "albedo image size doesn't match color image size";
    return false;
  }
  if (normal && normal->imsize() != color.imsize()) {
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
      "color", (void *)color.data(), oidn::Format::Float3, width, height);

  // set albedo if present
  if (albedo)
    filter.setImage(
        "albedo", (void *)albedo->data(), oidn::Format::Float3, width, height);

  // set normal if present
  if (normal)
    filter.setImage(
        "normal", (void *)normal->data(), oidn::Format::Float3, width, height);

  // initialize 'out' image to the correct size and set it as filter output
  out.resize(color.imsize());
  filter.setImage(
      "output", (void *)out.data(), oidn::Format::Float3, width, height);

  // register the user provided progress callback in the filter
  if (progress_cb) {
    auto prog_monitor = [](void *uptr, double prog) {
      int  current     = (int)(prog * 100);
      auto progress_cb = *(progress_callback *)uptr;
      progress_cb("denoise image", current, 100);
      return true;
    };
    filter.setProgressMonitorFunction(prog_monitor, (void *)&progress_cb);
  }

  // excecute the filter
  filter.commit();
  filter.execute();

  // check and return eventual oidn errors
  const char *errorMessage;
  if (device.getError(errorMessage) != oidn::Error::None) {
    error = errorMessage;
    return false;
  }

  return true;
}

bool oidn_image_denoise(const image<vec3f> &color, bool hdr, image<vec3f> &out,
    std::string &error, progress_callback progress_cb) {
  return oidn_image_denoise(
      color, hdr, nullptr, nullptr, out, error, progress_cb);
}

bool oidn_image_denoise(const image<vec3f> &color, bool hdr,
    const image<vec3f> &albedo, image<vec3f> &out, std::string &error,
    progress_callback progress_cb) {
  return oidn_image_denoise(
      color, hdr, &albedo, nullptr, out, error, progress_cb);
}

bool oidn_image_denoise(const image<vec3f> &color, bool hdr,
    const image<vec3f> &albedo, const image<vec3f> &normal, image<vec3f> &out,
    std::string &error, progress_callback progress_cb) {
  return oidn_image_denoise(
      color, hdr, &albedo, &normal, out, error, progress_cb);
}

}  // namespace yocto::denoise