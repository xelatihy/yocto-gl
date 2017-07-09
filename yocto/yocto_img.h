///
/// # Yocto/Img
///
/// Utilities for loading/saving images, mostly to avoid linking
/// problems and code duplication across other yocto libraries. Functions as
/// a wrapper to other libraries.
///
/// This library depends in yocto_math.h, stb_image.h, stb_image_write.h
/// stb_image_resize.h.
///
///
/// ## Usage
///
/// 1. load images with `load_image()`
/// 2. save images with `save_image()`
/// 2. resize images with `resize_image()`
///
///
/// ## History
///
/// - v 0.1: initial release
///
namespace yimg {}

//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
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

#ifndef _YIMG_H_
#define _YIMG_H_

#include <array>
#include <string>
#include <vector>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

///
/// Utilitied for reading and writing images.
///
namespace yimg {

// typedefs
using byte = unsigned char;

///
/// Loads an ldr image.
///
void load_image(
    const std::string& filename, int& w, int& h, int& ncomp, byte*& ldr);

///
/// Loads an hdr image.
///
void load_image(
    const std::string& filename, int& w, int& h, int& ncomp, float*& hdr);

///
/// Loads an image. Uses extension to determine which format to load.
/// Suppoted formats are PNG, JPG, HDR.
///
void load_image(const std::string& filename, int& w, int& h, int& ncomp,
    float*& hdr, byte*& ldr);

///
/// Loads an image.
///
void load_image(const std::string& filename, int& w, int& h, int& ncomp,
    std::vector<float>& hdr, std::vector<byte>& ldr);

///
/// Saves an hdr image. Uses extension to determine which format to load.
/// Suppoted formats are HDR.
///
void save_image(const std::string& filename, int width, int height, int ncomp,
    const float* hdr);

///
/// Saves an ldr image. Uses extension to determine which format to load.
/// Suppoted formats are PNG.
///
void save_image(const std::string& filename, int width, int height, int ncomp,
    const byte* ldr);

///
/// Loads an image from memory.
///
void load_image_from_memory(const std::string& fmt, const byte* data,
    int length, int& w, int& h, int& ncomp, float*& hdr, byte*& ldr);

///
/// Loads an image from memory.
///
void load_image_from_memory(const std::string& fmt, const byte* data,
    int length, int& w, int& h, int& ncomp, std::vector<float>& hdr,
    std::vector<byte>& ldr);

///
/// Resize image.
///
void resize_image(int width, int height, int ncomp, const float* img,
    int res_width, int res_height, float* res_img);

///
/// Resize image.
///
void resize_image(int width, int height, int ncomp, const byte* img,
    int res_width, int res_height, byte* res_img);

///
/// Resize image.
///
void resize_image(int width, int height, int ncomp,
    const std::vector<float>& img, int res_width, int res_height,
    std::vector<float>& res_img);

///
/// Resize image.
///
void resize_image(int width, int height, int ncomp,
    const std::vector<byte>& img, int res_width, int res_height,
    std::vector<byte>& res_img);

}  // namespace yimg

#endif
