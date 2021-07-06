//
// # Yocto/ImageIO: Image serialization
//
// Yocto/ImageIO supports loading and saving images from Png, Jpg, Tga, Bmp,
// Pnm, Hdr, Exr, Pfm.
// Yocto/ImageIO is implemented in `yocto_imageio.h` and `yocto_imageio.cpp`
// and depends on `stb_image.h`, `stb_image_write.h`, and `tinyexr.h`.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
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

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "yocto_image.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::pair;
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Check if an image is HDR or LDR based on filename.
bool is_hdr_filename(const string& filename);
bool is_ldr_filename(const string& filename);

// Loads/saves a 4 channels float/byte image in linear/srgb color space.
image_data load_image(const string& filename);
void       load_image(const string& filename, image_data& image);
void       save_image(const string& filename, const image_data& image);

// Make presets. Supported mostly in IO.
image_data make_image_preset(const string& type);

// Loads/saves a 4 channels float/byte image in linear/srgb color space.
bool load_image(const string& filename, image_data& img, string& error);
bool save_image(const string& filename, const image_data& img, string& error);

// Make presets. Supported mostly in IO.
bool make_image_preset(image_data& image, const string& type, string& error);

}  // namespace yocto

#endif
