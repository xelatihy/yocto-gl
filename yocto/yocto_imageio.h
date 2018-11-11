//
// # Yocto/ImageIO: Tiny library for image input and output
//
// Yocto/ImageIO provides loading and saving functionality for images.
// We support PNG, JPG, TGA, HDR, EXR formats.
//
// We do not use exception as the API for returning errors, although they might
// be used internally in the implementastion of the methods. In load functions,
// as error is signaled by returning an empty object or a null pointer. In
// save functions, errors are returned with the supplied boolean. In the future,
// we will also provide return types with error codes.
//
// ## Image Loading and Saving
//
// 1. load and save images with `load_image()` and `save_image()`
//
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
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

#include <string>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

using std::string;

}

// -----------------------------------------------------------------------------
// IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Check if an image is HDR based on filename.
bool is_hdr_filename(const string& filename);

// Loads/saves a 4 channel float image in linear color space.
bool load_image(const string& filename, image<vec4f>& img);
bool save_image(const string& filename, const image<vec4f>& img);
bool load_image_from_memory(const byte* data, int data_size, image<vec4f>& img);

// Loads/saves a 4 channel byte image in sRGB color space.
bool load_image(const string& filename, image<vec4b>& img);
bool save_image(const string& filename, const image<vec4b>& img);
bool load_image_from_memory(const byte* data, int data_size, image<vec4b>& img);
bool load_image_from_memory(const byte* data, int data_size, image<vec4b>& img);

// Convenience helper that saves an HDR images as wither a linear HDR file or
// a tonemapped LDR file depending on file name
bool save_tonemapped_image(const string& filename, const image<vec4f>& hdr,
    float exposure = 0, bool filmic = false, bool srgb = true);

}  // namespace yocto

// -----------------------------------------------------------------------------
// VOLUME IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Loads/saves a 1 channel volume.
bool load_volume(const string& filename, volume<float>& vol);
bool save_volume(const string& filename, const volume<float>& vol);

}  // namespace yocto

#endif
