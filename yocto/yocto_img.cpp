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
static inline std::string _get_extension(const std::string& filename) {
    auto pos = filename.rfind('.');
    if (pos == std::string::npos) return "";
    return filename.substr(pos);
}

//
// Loads an image
//
void load_image(const std::string& filename, int& width, int& height,
    int& ncomp, float*& hdr) {
    hdr = stbi_loadf(filename.c_str(), &width, &height, &ncomp, 0);
    if (!hdr) { throw std::runtime_error("cannot load image " + filename); }
}

//
// Loads an image
//
void load_image(const std::string& filename, int& width, int& height,
    int& ncomp, byte*& ldr) {
    ldr = stbi_load(filename.c_str(), &width, &height, &ncomp, 0);
    if (!ldr) { throw std::runtime_error("cannot load image " + filename); }
}

//
// Loads an image from memory.
//
void load_image_from_memory(const std::string& fmt, const byte* data,
    int length, int& width, int& height, int& ncomp, float*& hdr) {
    hdr = stbi_loadf_from_memory(data, length, &width, &height, &ncomp, 0);
    if (!*hdr) {
        throw std::runtime_error(
            "cannot load image from memory with format " + fmt);
    }
}

//
// Loads an image from memory.
//
void load_image_from_memory(const std::string& fmt, const byte* data,
    int length, int& width, int& height, int& ncomp, byte*& ldr) {
    ldr = stbi_load_from_memory(data, length, &width, &height, &ncomp, 0);
    if (!ldr) {
        throw std::runtime_error(
            "cannot load image from memory with format " + fmt);
    }
}

//
// Loads an image from memory.
//
void load_image_from_memory(const std::string& fmt, byte* data, int length,
    int& width, int& height, int& ncomp, float*& hdr, byte*& ldr) {
    if (fmt == "hdr")
        load_image_from_memory(fmt, data, length, width, height, ncomp, hdr);
    else
        load_image_from_memory(fmt, data, length, width, height, ncomp, ldr);
}

//
// Loads an image
//
void load_image(const std::string& filename, int& width, int& height,
    int& ncomp, float*& hdr, byte*& ldr) {
    if (_get_extension(filename) == ".hdr")
        load_image(filename, width, height, ncomp, hdr);
    else
        load_image(filename, width, height, ncomp, ldr);
}

//
// Loads an image
//
void load_image(const std::string& filename, int& width, int& height,
    int& ncomp, std::vector<float>& hdr, std::vector<byte>& ldr) {
    auto hdr_ = (float*)nullptr;
    auto ldr_ = (byte*)nullptr;
    load_image(filename, width, height, ncomp, hdr_, ldr_);
    if (hdr_) {
        hdr = std::vector<float>(hdr_, hdr_ + width * height * ncomp);
        delete hdr_;
    }
    if (ldr_) {
        ldr = std::vector<unsigned char>(ldr_, ldr_ + width * height * ncomp);
        delete ldr_;
    }
}

//
// Loads an image from memory.
//
void load_image_from_memory(const std::string& fmt, const byte* data,
    int length, int& width, int& height, int& ncomp, std::vector<float>& hdr,
    std::vector<byte>& ldr) {
    auto hdr_ = (float*)nullptr;
    auto ldr_ = (byte*)nullptr;
    load_image_from_memory(
        fmt, (byte*)data, length, width, height, ncomp, hdr_, ldr_);
    if (hdr_) {
        hdr = std::vector<float>(hdr_, hdr_ + width * height * ncomp);
        delete hdr_;
    }
    if (ldr_) {
        ldr = std::vector<unsigned char>(ldr_, ldr_ + width * height * ncomp);
        delete ldr_;
    }
}

//
// Saves an image
//
void save_image(const std::string& filename, int width, int height, int ncomp,
    const float* hdr) {
    if (_get_extension(filename) == ".hdr") {
        stbi_write_hdr(filename.c_str(), width, height, ncomp, hdr);
    } else {
        throw std::invalid_argument("unsupported output extension " + filename);
    }
}

//
// Saves an image
//
void save_image(const std::string& filename, int width, int height, int ncomp,
    const byte* ldr) {
    if (_get_extension(filename) == ".png") {
        if (!ldr) throw std::invalid_argument("ldr data required");
        stbi_write_png(
            filename.c_str(), width, height, ncomp, ldr, width * ncomp);
    } else {
        throw std::invalid_argument("unsupported output extension " + filename);
    }
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
