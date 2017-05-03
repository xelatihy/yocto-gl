///
/// YOCTO_IMAGE: utilities for loading/saving images, mostly to avoid linking
/// problems and code duplication across other yocto libraries.
///
///
/// USAGE:
///
/// 1. load images with load_image()
/// 2. save images with save_image()
///
/// Notes: Images own their own memory and are cleaned up automatically when
/// destroyed. If you want to get the buffer without copying use release_hdr()
/// and release_ldr().
///
///
/// COMPILATION:
///
/// To use the library include the .h and compile the .cpp. To use this library
/// as a header-only library, define YIMG_INLINE before including this file.
///
/// Image loading and savuing depends on stb libraries. If not desired, disabled
/// by definining YIMG_NO_STBIMAGE before compiling the .cpp file.
///
///
/// HISTORY:
/// - v 0.1: initial release
///
namespace ycmd {}

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

// compilation options
#ifdef YIMG_INLINE
#define YIMG_API inline
#else
#define YIMG_API
#endif

#include <array>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

namespace yimg {

///
/// Typedefs
///
using byte = unsigned char;
using float3 = std::array<float, 3>;
using byte3 = std::array<byte, 3>;
using float4 = std::array<float, 4>;
using byte4 = std::array<byte, 4>;

///
/// Simple image structure to ease the passing of parameters. Memory is not
/// managed.
///
struct simage {
    /// image width
    int width = 0;

    //. image height
    int height = 0;

    /// number of components
    int ncomp = 0;

    /// float data for hdr images if loaded
    float* hdr = nullptr;

    /// char data for ldr images if loaded
    byte* ldr = nullptr;

    /// default constructor
    simage() {}

    /// allocating constructor
    simage(int width, int height, int ncomp, bool ishdr)
        : width(width)
        , height(height)
        , ncomp(ncomp)
        , hdr((ishdr) ? new float[width * height * ncomp] : nullptr)
        , ldr((ishdr) ? nullptr : new byte[width * height * ncomp]) {}

    /// destructor
    ~simage() {
        if (ldr) delete[] ldr;
        if (hdr) delete[] hdr;
    }
};

///
/// Initializes an image
///
YIMG_API simage* make_image(int width, int height, int ncomp, bool hdr);

///
/// Removes the hdr buffer from the image.
///
YIMG_API float* release_hdr(simage* img);

///
/// Removes the ldr buffer from the image.
///
YIMG_API byte* release_ldr(simage* img);

///
/// Loads an image. Uses extension to determine which format to load.
/// Suppoted formats are PNG, JPG, HDR.
///
YIMG_API simage* load_image(const std::string& filename);

///
/// Loads an image. Unwrap params compared to previous one.
///
YIMG_API void load_image(const std::string& filename, int& w, int& h,
    int& ncomp, float*& hdr, byte*& ldr);

///
/// Loads an image flipping Y.
///
YIMG_API simage* load_image_flipy(const std::string& filename);

///
/// Loads an image from memory.
///
YIMG_API simage* load_image_from_memory(
    const std::string& fmt, byte* data, int length);

///
/// Saves an image. Uses extension and image content to determine
/// which format to save to. Suppoted formats are PNG, HDR.
///
YIMG_API void save_image(const std::string& filename, const simage* img);

///
/// Saves an image. Unwrap params compared to previous one.
///
YIMG_API void save_image(const std::string& filename, int width, int height,
    int ncomp, const float* hdr, const byte* ldr);

//
// Resize an image. If width or height are less than 0, they are set
// automatically to maintain aspect ratio.
//
YIMG_API simage* resize_image(const simage* img, int res_width, int res_height);

//
// Resize image. Unwrap params compared to previous one.
//
YIMG_API void resize_image(int width, int height, int ncomp, const float* hdr,
    const byte* ldr, int& res_width, int& res_height, float*& res_hdr,
    byte*& res_ldr);

//
// Tone mapping configurations
//
enum struct tonemap_type { def = 0, linear, srgb, gamma, filmic };

///
/// Apply tone mapping operator to a pixel.
///
YIMG_API byte3 tonemap_pixel(
    const float3& hdr, float exposure, tonemap_type tm, float gamma);

///
/// Apply tone mapping operator to a pixel.
///
YIMG_API byte4 tonemap_pixel(
    const float4& hdr, float exposure, tonemap_type tm, float gamma);

///
/// Tone mapping HDR to LDR images.
///
YIMG_API simage* tonemap_image(
    simage* img, float exposure, tonemap_type tm, float gamma);

///
/// Tone mapping HDR to LDR images.
///
YIMG_API void tonemap_image(int width, int height, int ncomp, const float* hdr,
    byte* ldr, float exposure, tonemap_type tm, float gamma);
}  // namespace yimg

// -----------------------------------------------------------------------------
// INCLUDE FOR HEADER-ONLY MODE
// -----------------------------------------------------------------------------

#ifdef YIMG_INLINE
#include "yocto_img.cpp"
#endif

#endif
