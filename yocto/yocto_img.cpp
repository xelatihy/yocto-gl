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
// Initializes an image
//
YIMG_API simage* make_image(int width, int height, int ncomp, bool hdr) {
    return new simage(width, height, ncomp, hdr);
}

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
YIMG_API simage* load_image(const std::string& filename) {
    auto img = new simage();
    if (_get_extension(filename) == ".hdr") {
        img->hdr = stbi_loadf(
            filename.c_str(), &img->width, &img->height, &img->ncomp, 0);
    } else {
        img->ldr = stbi_load(
            filename.c_str(), &img->width, &img->height, &img->ncomp, 0);
    }
    if (!img->ldr && !img->hdr) {
        delete img;
        throw std::runtime_error("cannot load image " + filename);
    }
    return img;
}

//
// Loads an image
//
YIMG_API void load_image(const std::string& filename, int& width, int& height,
    int& ncomp, float*& hdr, byte*& ldr) {
    if (_get_extension(filename) == ".hdr") {
        hdr = stbi_loadf(filename.c_str(), &width, &height, &ncomp, 0);
    } else {
        ldr = stbi_load(filename.c_str(), &width, &height, &ncomp, 0);
    }
    if (!hdr && !ldr) {
        throw std::runtime_error("cannot load image " + filename);
    }
}

//
// Loads an image flipping Y.
//
YIMG_API simage* load_image_flipy(const std::string& filename) {
    stbi_set_flip_vertically_on_load(true);
    auto img = load_image(filename);
    stbi_set_flip_vertically_on_load(false);
    return img;
}

//
// Loads an image from memory.
//
YIMG_API simage* load_image_from_memory(
    const std::string& fmt, byte* data, int length) {
    auto img = new simage();
    if (fmt == "hdr") {
        img->hdr = stbi_loadf_from_memory(
            data, length, &img->width, &img->height, &img->ncomp, 0);
    } else {
        img->ldr = stbi_load_from_memory(
            data, length, &img->width, &img->height, &img->ncomp, 0);
    }
    if (!img->ldr && !img->hdr) {
        delete img;
        throw std::runtime_error(
            "cannot load image from memory with format " + fmt);
    }
    return img;
}

//
// Saves an image
//
YIMG_API void save_image(const std::string& filename, const simage* img) {
    save_image(
        filename, img->width, img->height, img->ncomp, img->hdr, img->ldr);
}

//
// Saves an image
//
YIMG_API void save_image(const std::string& filename, int width, int height,
    int ncomp, const float* hdr, const byte* ldr) {
    if (_get_extension(filename) == ".hdr") {
        if (!hdr) throw std::invalid_argument("hdr data required");
        stbi_write_hdr(filename.c_str(), width, height, ncomp, hdr);
    } else if (_get_extension(filename) == ".png") {
        if (!ldr) throw std::invalid_argument("ldr data required");
        stbi_write_png(
            filename.c_str(), width, height, ncomp, ldr, width * ncomp);
    } else {
        throw std::invalid_argument("unsupported output extension " + filename);
    }
}

//
// Removes the hdr buffer from the image.
//
YIMG_API float* release_hdr(simage* img) {
    auto ptr = img->hdr;
    img->hdr = nullptr;
    return ptr;
}

//
// Removes the ldr buffer from the image.
//
YIMG_API byte* release_ldr(simage* img) {
    auto ptr = img->ldr;
    img->ldr = nullptr;
    return ptr;
}

//
// Resize image.
//
YIMG_API simage* resize_image(
    const simage* img, int res_width, int res_height) {
    auto res = new simage();
    resize_image(img->width, img->height, img->ncomp, img->hdr, img->ldr,
        res->width, res->height, res->hdr, res->ldr);
    return res;
}

//
// Resize image.
//
YIMG_API void resize_image(int width, int height, int ncomp, const float* hdr,
    const byte* ldr, int& res_width, int& res_height, float*& res_hdr,
    byte*& res_ldr) {
    if (res_width < 0 && res_height < 0)
        throw std::invalid_argument("at least argument should be >0");
    if (res_width < 0)
        res_width = (int)std::round(width * (res_height / (float)height));
    if (res_height < 0)
        res_height = (int)std::round(height * (res_width / (float)width));
    if (hdr) {
        res_hdr = new float[width * height * ncomp];
        auto img_stride = sizeof(float) * width * ncomp;
        auto res_stride = sizeof(float) * res_width * ncomp;
        stbir_resize_float(hdr, width, height, img_stride, res_hdr, res_width,
            res_height, res_stride, ncomp);
    } else {
        res_ldr = new byte[width * height * ncomp];
        auto img_stride = sizeof(byte) * width * ncomp;
        auto res_stride = sizeof(byte) * res_width * ncomp;
        stbir_resize_uint8_srgb(ldr, width, height, img_stride, res_ldr,
            res_width, res_height, res_stride, ncomp,
            (ncomp == 4) ? 3 : STBIR_ALPHA_CHANNEL_NONE, 0);
    }
}

//
// Tone mapping HDR to LDR images.
//
YIMG_API simage* tonemap_image(
    simage* hdr, float exposure, tonemap_type tm, float gamma) {
    if (!hdr->hdr) throw std::invalid_argument("tonemap hdr only");
    auto ldr = make_image(hdr->width, hdr->height, hdr->ncomp, false);
    tonemap_image(hdr->width, hdr->height, hdr->ncomp, hdr->hdr, ldr->ldr,
        exposure, tm, gamma);
    return ldr;
}

//
// float to byte
//
static inline byte _float_to_byte(float x) {
    return std::max(0, std::min(int(x * 256), 255));
}

#if 1
//
// filmic curve
// https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
//
static inline float _filmic_tonemapping(float x) {
    // rescaler ange
    x *= 2.05f;
    float a = 2.51f;
    float b = 0.03f;
    float c = 2.43f;
    float d = 0.59f;
    float e = 0.14f;
    auto y = ((x * (a * x + b)) / (x * (c * x + d) + e));
    return std::pow(std::min(1.0f, std::max(0.0f, y)), 1 / 2.2f);
}
#else
static inline float _filmic_tonemapping(float x) {
    auto y =
        (x * (x * (x * (x * 2708.7142 + 6801.1525) + 1079.5474) + 1.1614649) -
            0.00004139375) /
        (x * (x * (x * (x * 983.38937 + 4132.0662) + 2881.6522) + 128.35911) +
            1.0);
    return (float)std::max(y, 0.0);
}
#endif

//
// Apply tone mapping operator to a pixel.
//
YIMG_API byte3 tonemap_pixel(
    const float3& hdr, float exposure, tonemap_type tm, float gamma) {
    auto s = std::pow(2, exposure);
    switch (tm) {
        case tonemap_type::linear:
            return {_float_to_byte(s * hdr[0]), _float_to_byte(s * hdr[1]),
                _float_to_byte(s * hdr[2])};
        case tonemap_type::def:
        case tonemap_type::srgb: {
            return {_float_to_byte(std::pow(s * hdr[0], 1 / 2.2f)),
                _float_to_byte(std::pow(s * hdr[1], 1 / 2.2f)),
                _float_to_byte(std::pow(s * hdr[2], 1 / 2.2f))};
        } break;
        case tonemap_type::gamma: {
            return {_float_to_byte(std::pow(s * hdr[0], 1 / gamma)),
                _float_to_byte(std::pow(s * hdr[1], 1 / gamma)),
                _float_to_byte(std::pow(s * hdr[2], 1 / gamma))};
        } break;
        case tonemap_type::filmic: {
            return {_float_to_byte(_filmic_tonemapping(s * hdr[0])),
                _float_to_byte(_filmic_tonemapping(s * hdr[1])),
                _float_to_byte(_filmic_tonemapping(s * hdr[2]))};
        } break;
    }
}

//
// Apply tone mapping operator to a pixel.
//
YIMG_API byte4 tonemap_pixel(
    const float4& hdr, float exposure, tonemap_type tm, float gamma) {
    auto ldr = byte4();
    (byte3&)ldr = tonemap_pixel((const float3&)hdr, exposure, tm, gamma);
    ldr[3] = _float_to_byte(hdr[3]);
    return ldr;
}

//
// Tone mapping HDR to LDR images.
//
YIMG_API void tonemap_image(int width, int height, int ncomp, const float* hdr,
    byte* ldr, float exposure, tonemap_type tm, float gamma) {
    if (ncomp < 3 || ncomp > 4)
        throw std::invalid_argument("tonemap supports 3-4 channels only");
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto h = hdr + (j * width + i) * ncomp;
            auto l = ldr + (j * width + i) * ncomp;
            if (ncomp == 3)
                *(byte3*)l = tonemap_pixel(*(float3*)h, exposure, tm, gamma);
            else
                *(byte4*)l = tonemap_pixel(*(float4*)h, exposure, tm, gamma);
        }
    }
}

}  // namespace yimg
