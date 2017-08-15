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
#include <map>
#include <memory>

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

#define TINYEXR_IMPLEMENTATION
#include "ext/tinyexr.h"

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
// Check if an image is HDR based on filename
//
bool is_hdr_filename(const std::string& filename) {
    auto ext = get_extension(filename);
    return ext == ".hdr" || ext == ".exr";
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
    auto ext = get_extension(filename);
    auto w = 0, h = 0, c = 0;
    auto pixels = std::unique_ptr<float>(nullptr);
    if (ext == ".exr") {
        auto pixels_ = (float*)nullptr;
        if (!LoadEXR(&pixels_, &w, &h, filename.c_str(), nullptr))
            pixels = std::unique_ptr<float>(pixels_);
    } else {
        pixels =
            std::unique_ptr<float>(stbi_loadf(filename.c_str(), &w, &h, &c, 4));
    }
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
    } else if (get_extension(filename) == ".jpg") {
        return stbi_write_jpg(filename.c_str(), img.width(), img.height(), 4,
            (byte*)img.data(), 75);
    } else {
        return false;
    }
}

//
// Saves an hdr image.
//
bool save_image4f(const std::string& filename, const ym::image4f& img) {
    if (get_extension(filename) == ".hdr") {
        return stbi_write_hdr(
            filename.c_str(), img.width(), img.height(), 4, (float*)img.data());
    } else if (get_extension(filename) == ".exr") {
        return !SaveEXR(
            (float*)img.data(), img.width(), img.height(), 4, filename.c_str());
    } else {
        return false;
    }
}

//
// Loads an image
//
std::vector<float> load_imagef(
    const std::string& filename, int& width, int& height, int& ncomp) {
    auto pixels = stbi_loadf(filename.c_str(), &width, &height, &ncomp, 0);
    if (!pixels) return {};
    auto ret = std::vector<float>(pixels, pixels + width * height * ncomp);
    free(pixels);
    return ret;
}

//
// Loads an image
//
std::vector<byte> load_image(
    const std::string& filename, int& width, int& height, int& ncomp) {
    auto pixels = stbi_load(filename.c_str(), &width, &height, &ncomp, 0);
    if (!pixels) return {};
    auto ret = std::vector<byte>(pixels, pixels + width * height * ncomp);
    free(pixels);
    return ret;
}

//
// Loads an image from memory.
//
std::vector<float> load_imagef_from_memory(const std::string& filename,
    const byte* data, int length, int& width, int& height, int& ncomp) {
    auto pixels =
        stbi_loadf_from_memory(data, length, &width, &height, &ncomp, 0);
    if (!pixels) return {};
    auto ret = std::vector<float>(pixels, pixels + width * height * ncomp);
    free(pixels);
    return ret;
}

//
// Loads an image from memory.
//
std::vector<byte> load_image_from_memory(const std::string& filename,
    const byte* data, int length, int& width, int& height, int& ncomp) {
    auto pixels =
        stbi_load_from_memory(data, length, &width, &height, &ncomp, 0);
    if (!pixels) return {};
    auto ret = std::vector<byte>(pixels, pixels + width * height * ncomp);
    free(pixels);
    return ret;
}

//
// Saves an image
//
bool save_imagef(const std::string& filename, int width, int height, int ncomp,
    const float* hdr) {
    if (get_extension(filename) == ".hdr") {
        return stbi_write_hdr(filename.c_str(), width, height, ncomp, hdr);
    } else {
        return false;
    }
}

//
// Saves an image
//
bool save_image(const std::string& filename, int width, int height, int ncomp,
    const byte* ldr) {
    if (get_extension(filename) == ".png") {
        return stbi_write_png(
            filename.c_str(), width, height, ncomp, ldr, width * ncomp);
    } else if (get_extension(filename) == ".jpg") {
        return stbi_write_jpg(filename.c_str(), width, height, ncomp, ldr, 75);
    } else {
        return false;
    }
}

static const auto filter_map = std::map<resize_filter, stbir_filter>{
    {resize_filter::def, STBIR_FILTER_DEFAULT},
    {resize_filter::box, STBIR_FILTER_BOX},
    {resize_filter::triangle, STBIR_FILTER_TRIANGLE},
    {resize_filter::cubic_spline, STBIR_FILTER_CUBICBSPLINE},
    {resize_filter::catmull_rom, STBIR_FILTER_CATMULLROM},
    {resize_filter::mitchell, STBIR_FILTER_MITCHELL}};

static const auto edge_map = std::map<resize_edge, stbir_edge>{
    {resize_edge::def, STBIR_EDGE_CLAMP},
    {resize_edge::clamp, STBIR_EDGE_CLAMP},
    {resize_edge::reflect, STBIR_EDGE_REFLECT},
    {resize_edge::wrap, STBIR_EDGE_WRAP}, {resize_edge::zero, STBIR_EDGE_ZERO}};

///
/// Resize image.
///
void resize_image(const ym::image4f& img, ym::image4f& res_img,
    resize_filter filter, resize_edge edge, bool premultiplied_alpha) {
    stbir_resize_float_generic((float*)img.data(), img.width(), img.height(),
        sizeof(ym::vec4f) * img.width(), (float*)res_img.data(),
        res_img.width(), res_img.height(), sizeof(ym::vec4f) * res_img.width(),
        4, 3, (premultiplied_alpha) ? STBIR_FLAG_ALPHA_PREMULTIPLIED : 0,
        edge_map.at(edge), filter_map.at(filter), STBIR_COLORSPACE_LINEAR,
        nullptr);
}

///
/// Resize image.
///
void resize_image(const ym::image4b& img, ym::image4b& res_img,
    resize_filter filter, resize_edge edge, bool premultiplied_alpha) {
    stbir_resize_uint8_generic((unsigned char*)img.data(), img.width(),
        img.height(), sizeof(ym::vec4b) * img.width(),
        (unsigned char*)res_img.data(), res_img.width(), res_img.height(),
        sizeof(ym::vec4b) * res_img.width(), 4, 3,
        (premultiplied_alpha) ? STBIR_FLAG_ALPHA_PREMULTIPLIED : 0,
        edge_map.at(edge), filter_map.at(filter), STBIR_COLORSPACE_LINEAR,
        nullptr);
}

}  // namespace yimg
