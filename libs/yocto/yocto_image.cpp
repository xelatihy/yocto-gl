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

#include <future>
#include <memory>
#include <stdexcept>

#include "yocto_color.h"
#include "yocto_noise.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::unique_ptr;

}  // namespace yocto

// -----------------------------------------------------------------------------
// PARALLEL HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename T, typename Func>
inline void parallel_for_batch(T num, T batch, Func&& func) {
  auto              futures  = vector<std::future<void>>{};
  auto              nthreads = std::thread::hardware_concurrency();
  std::atomic<T>    next_idx(0);
  std::atomic<bool> has_error(false);
  for (auto thread_id = 0; thread_id < (int)nthreads; thread_id++) {
    futures.emplace_back(std::async(
        std::launch::async, [&func, &next_idx, &has_error, num, batch]() {
          try {
            while (true) {
              auto start = next_idx.fetch_add(batch);
              if (start >= num) break;
              if (has_error) break;
              auto end = std::min(num, start + batch);
              for (auto i = (T)start; i < end; i++) func(i);
            }
          } catch (...) {
            has_error = true;
            throw;
          }
        }));
  }
  for (auto& f : futures) f.get();
}

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename T, typename Func>
inline void parallel_for_batch(array<T, 2> num, Func&& func) {
  auto              futures  = vector<std::future<void>>{};
  auto              nthreads = std::thread::hardware_concurrency();
  std::atomic<T>    next_idx(0);
  std::atomic<bool> has_error(false);
  for (auto thread_id = 0; thread_id < (int)nthreads; thread_id++) {
    futures.emplace_back(
        std::async(std::launch::async, [&func, &next_idx, &has_error, num]() {
          try {
            while (true) {
              auto j = next_idx.fetch_add(1);
              if (j >= num[1]) break;
              if (has_error) break;
              for (auto i = (T)0; i < num[0]; i++) func({i, j});
            }
          } catch (...) {
            has_error = true;
            throw;
          }
        }));
  }
  for (auto& f : futures) f.get();
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Resize an image.
array2d<vec4f> resize_image(
    const array2d<vec4f>& image, const vec2uz& extents_) {
  // determine new size
  auto extents = extents_;
  auto aspect  = (double)image.extent(0) / (double)image.extent(1);
  if (extents == vec2uz{0, 0})
    throw std::invalid_argument{"bad image size in resize"};
  if (extents[1] == 0) {
    extents = {extents[0], (size_t)round(extents[0] / aspect)};
  } else if (extents[0] == 0) {
    extents = {(size_t)round(extents[1] * aspect), extents[1]};
  }

  // resize
  auto result = array2d<vec4f>(extents);
  stbir_resize_float_generic((float*)image.data(), (int)image.extent(0),
      (int)image.extent(1), (int)(sizeof(vec4f) * image.extent(0)),
      (float*)result.data(), (int)result.extent(0), (int)result.extent(1),
      (int)(sizeof(vec4f) * result.extent(0)), 4, 3, 0, STBIR_EDGE_CLAMP,
      STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
  return result;
}

}  // namespace yocto
