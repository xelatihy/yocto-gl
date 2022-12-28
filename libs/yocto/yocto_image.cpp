//
// Implementation for Yocto/Image.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2022 Fabio Pellacini
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

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_image.h"

#include <stb_image/stb_image_resize.h>

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Resize an image.
template <size_t N>
array2d<vec<float, N>> resize_image(
    const array2d<vec<float, N>>& image, const vec2s& extents_) {
  // determine new size
  auto extents = extents_;
  auto aspect  = (double)image.extent(0) / (double)image.extent(1);
  if (extents == vec2s{0, 0})
    throw std::invalid_argument{"bad image size in resize"};
  if (extents[1] == 0) {
    extents = {extents[0], (size_t)round(extents[0] / aspect)};
  } else if (extents[0] == 0) {
    extents = {(size_t)round(extents[1] * aspect), extents[1]};
  }

  // alpha channel index
  auto alpha_index = N == 4 ? 3 : STBIR_ALPHA_CHANNEL_NONE;

  // resize
  auto result = array2d<vec<float, N>>(extents);
  stbir_resize_float_generic((float*)image.data(), (int)image.extent(0),
      (int)image.extent(1), (int)(sizeof(float) * N * image.extent(0)),
      (float*)result.data(), (int)result.extent(0), (int)result.extent(1),
      (int)(sizeof(float) * N * result.extent(0)), N, alpha_index, 0,
      STBIR_EDGE_CLAMP, STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
  return result;
}

// Explicit instantiations
template array2d<vec1f> resize_image(const array2d<vec1f>&, const vec2s&);
template array2d<vec2f> resize_image(const array2d<vec2f>&, const vec2s&);
template array2d<vec3f> resize_image(const array2d<vec3f>&, const vec2s&);
template array2d<vec4f> resize_image(const array2d<vec4f>&, const vec2s&);

// Resize an image.
template <size_t N>
array2d<vec<byte, N>> resize_image(
    const array2d<vec<byte, N>>& image, const vec2s& extents_) {
  // determine new size
  auto extents = extents_;
  auto aspect  = (double)image.extent(0) / (double)image.extent(1);
  if (extents == vec2s{0, 0})
    throw std::invalid_argument{"bad image size in resize"};
  if (extents[1] == 0) {
    extents = {extents[0], (size_t)round(extents[0] / aspect)};
  } else if (extents[0] == 0) {
    extents = {(size_t)round(extents[1] * aspect), extents[1]};
  }

  // alpha channel index
  auto alpha_index = N == 4 ? 3 : STBIR_ALPHA_CHANNEL_NONE;

  // resize
  auto result = array2d<vec<byte, N>>(extents);
  stbir_resize_uint8_generic((byte*)image.data(), (int)image.extent(0),
      (int)image.extent(1), (int)(sizeof(byte) * N * image.extent(0)),
      (byte*)result.data(), (int)result.extent(0), (int)result.extent(1),
      (int)(sizeof(byte) * N * result.extent(0)), N, alpha_index, 0,
      STBIR_EDGE_CLAMP, STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
  return result;
}

// Explicit instantiations
template array2d<vec1b> resize_image(const array2d<vec1b>&, const vec2s&);
template array2d<vec2b> resize_image(const array2d<vec2b>&, const vec2s&);
template array2d<vec3b> resize_image(const array2d<vec3b>&, const vec2s&);
template array2d<vec4b> resize_image(const array2d<vec4b>&, const vec2s&);

}  // namespace yocto
