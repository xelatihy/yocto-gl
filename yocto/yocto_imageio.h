//
// # Yocto/ImageIO: Tiny library for image input and output
//
// Yocto/ImageIO provides loading and saving functionality for images.
// We support PNG, JPG, TGA, HDR, EXR formats.
//
// Error reporting is done through exceptions using the `io_error` exception.
//
// ## Image Loading and Saving
//
// 1. load and save images with `load_image()` and `save_image()`
//
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#ifndef _YOCTO_IMAGEIO_H_
#define _YOCTO_IMAGEIO_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_image.h"

// -----------------------------------------------------------------------------
// IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Check if an image is HDR based on filename.
bool is_hdr_filename(const string& filename);

// Loads/saves a 1-4 channels float image in linear color space.
template <int N>
void load_image(const string& filename, image<vec<float, N>>& img);
template <int N>
void save_image(const string& filename, const image<vec<float, N>>& img);
template <int N>
void load_image_from_memory(
    const byte* data, int data_size, image<vec<float, N>>& img);

// Loads/saves a 1-4 byte image in sRGB color space.
template <int N>
void load_image(const string& filename, image<vec<byte, N>>& img);
template <int N>
void save_image(const string& filename, const image<vec<byte, N>>& img);
template <int N>
void load_image_from_memory(
    const byte* data, int data_size, image<vec<byte, N>>& img);
template <int N>
void load_image_from_memory(
    const byte* data, int data_size, image<vec<byte, N>>& img);

// Convenience helper that saves an HDR images as wither a linear HDR file or
// a tonemapped LDR file depending on file name
template <int N>
inline void save_tonemapped_image(const string& filename,
    const image<vec<float, N>>& hdr, float exposure, bool filmic, bool srgb) {
    if (is_hdr_filename(filename)) {
        save_image(filename, hdr);
    } else {
        auto ldr = image<vec<byte, N>>{hdr.size()};
        tonemap_image8(ldr, hdr, exposure, filmic, srgb);
        save_image(filename, ldr);
    }
}

// imageio error
struct imageio_error : runtime_error {
    explicit imageio_error(const char* msg) : runtime_error{msg} {}
    explicit imageio_error(const std::string& msg) : runtime_error{msg} {}
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// VOLUME IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Loads/saves a 1 channel volume.
void load_volume(const string& filename, volume1f& vol);
void save_volume(const string& filename, const volume1f& vol);

}  // namespace yocto

// -----------------------------------------------------------------------------
// BUILTIN IMAGES
// -----------------------------------------------------------------------------
namespace yocto {

// Loads/saves a 1-4 channel builtin image.
template <int N>
void load_builtin_image(const string& name, image<vec<byte, N>>& img);
template <int N>
void load_builtin_image(const string& name, image<vec<float, N>>& img);

// Save with a logo embedded
template <typename T, int N>
inline void save_image_with_logo(
    const string& filename, const image<vec<T, N>>& img) {
    auto logo = image<vec<T, N>>{};
    load_builtin_image("logo-render", logo);
    auto img_copy = img;
    auto offset   = img.size() - logo.size() - 8;
    set_image_region(img_copy, logo, offset);
    save_image(filename, img_copy);
}

// Convenience helper that saves an HDR images as wither a linear HDR file or
// a tonemapped LDR file depending on file name
template <int N>
inline void save_tonemapped_image_with_logo(const string& filename,
    const image<vec<float, N>>& hdr, float exposure, bool filmic, bool srgb) {
    if (is_hdr_filename(filename)) {
        save_image_with_logo(filename, hdr);
    } else {
        auto ldr = image<vec<byte, N>>{hdr.size()};
        tonemap_image8(ldr, hdr, exposure, filmic, srgb);
        save_image_with_logo(filename, ldr);
    }
}

}  // namespace yocto

#endif
