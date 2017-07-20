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

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF YOCTO_IMG
// -----------------------------------------------------------------------------

#include "yocto_img.h"

#include <cmath>

#ifndef YIMG_NO_STBIMAGE

#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wunused-variable"
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#endif

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_STATIC
#include "ext/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_WRITE_STATIC
#include "ext/stb_image_write.h"

#define STB_IMAGE_RESIZE_IMPLEMENTATION
#define STB_IMAGE_RESIZE_STATIC
#include "ext/stb_image_resize.h"

#ifndef _WIN32
#pragma GCC diagnostic pop
#endif

#endif

namespace yimg {
//
// Get extension (including '.').
//
std::string get_extension(const std::string& filename) {
    auto pos = filename.rfind('.');
    if (pos == std::string::npos) return "";
    return filename.substr(pos);
}

//
// Loads an ldr image.
//
ym::image4b load_image4b(const std::string& filename) {
    auto w = 0, h = 0, c = 0;
    auto pixels =
        std::unique_ptr<byte>(stbi_load(filename.c_str(), &w, &h, &c, 4));
    if (!pixels) return {};
    return ym::image4b(w, h, (ym::vec4b*)pixels.get());
}

//
// Loads an hdr image.
//
ym::image4f load_image4f(const std::string& filename) {
    auto w = 0, h = 0, c = 0;
    auto pixels =
        std::unique_ptr<float>(stbi_loadf(filename.c_str(), &w, &h, &c, 4));
    if (!pixels) return {};
    return ym::image4f(w, h, (ym::vec4f*)pixels.get());
}

//
// Saves an ldr image.
//
bool save_image4b(const std::string& filename, const ym::image4b& img) {
    if (get_extension(filename) == ".png") {
        return stbi_write_png(filename.c_str(), img.width(), img.height(), 4,
            (byte*)img.data(), img.width() * 4);
    }
    return false;
}

//
// Saves an hdr image.
//
bool save_image4f(const std::string& filename, const ym::image4f& img) {
    if (get_extension(filename) == ".hdr") {
        return stbi_write_hdr(
            filename.c_str(), img.width(), img.height(), 4, (float*)img.data());
    }
    return false;
}

//
// Loads an image
//
float* load_imagef(
    const std::string& filename, int& width, int& height, int& ncomp) {
    return stbi_loadf(filename.c_str(), &width, &height, &ncomp, 0);
}

//
// Loads an image
//
byte* load_image(
    const std::string& filename, int& width, int& height, int& ncomp) {
    return stbi_load(filename.c_str(), &width, &height, &ncomp, 0);
}

//
// Loads an image from memory.
//
float* load_imagef_from_memory(const std::string& fmt, const byte* data,
    int length, int& width, int& height, int& ncomp) {
    return stbi_loadf_from_memory(data, length, &width, &height, &ncomp, 0);
}

//
// Loads an image from memory.
//
byte* load_image_from_memory(const std::string& fmt, const byte* data,
    int length, int& width, int& height, int& ncomp) {
    return stbi_load_from_memory(data, length, &width, &height, &ncomp, 0);
}

//
// Saves an image
//
bool save_imagef(const std::string& filename, int width, int height, int ncomp,
    const float* hdr) {
    if (get_extension(filename) == ".hdr") {
        return stbi_write_hdr(filename.c_str(), width, height, ncomp, hdr);
    } else
        return false;
}

//
// Saves an image
//
bool save_image(const std::string& filename, int width, int height, int ncomp,
    const byte* ldr) {
    if (get_extension(filename) == ".png") {
        return stbi_write_png(
            filename.c_str(), width, height, ncomp, ldr, width * ncomp);
    } else
        return false;
}

//
// Resize image.
//
void resize_image(int width, int height, int ncomp, const float* hdr,
    int res_width, int res_height, float* res_hdr) {
    auto img_stride = (int)sizeof(float) * width * ncomp;
    auto res_stride = (int)sizeof(float) * res_width * ncomp;
    stbir_resize_float(hdr, width, height, img_stride, res_hdr, res_width,
        res_height, res_stride, ncomp);
}

//
// Resize image.
//
void resize_image(int width, int height, int ncomp, const byte* ldr,
    int res_width, int res_height, byte* res_ldr) {
    auto img_stride = (int)sizeof(byte) * width * ncomp;
    auto res_stride = (int)sizeof(byte) * res_width * ncomp;
    stbir_resize_uint8_srgb(ldr, width, height, img_stride, res_ldr, res_width,
        res_height, res_stride, ncomp,
        (ncomp == 4) ? 3 : STBIR_ALPHA_CHANNEL_NONE, 0);
}

//
// Resize image.
//
void resize_image(int width, int height, int ncomp,
    const std::vector<float>& img, int res_width, int res_height,
    std::vector<float>& res_img) {
    res_img.resize(res_width * res_height);
    resize_image(width, height, ncomp, img.data(), res_width, res_height,
        res_img.data());
}

//
// Resize image.
//
void resize_image(int width, int height, int ncomp,
    const std::vector<byte>& img, int res_width, int res_height,
    std::vector<byte>& res_img) {
    res_img.resize(res_width * res_height);
    resize_image(width, height, ncomp, img.data(), res_width, res_height,
        res_img.data());
}

}  // namespace yimg
