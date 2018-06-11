//
// Implementation for Yocto/Image.
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

#include "yocto_image.h"
#include "yocto_utils.h"

#include <cstdlib>

#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wunused-variable"
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#endif

#ifndef __clang_analyzer__

#define STB_IMAGE_IMPLEMENTATION
#include "ext/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "ext/stb_image_write.h"

#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "ext/stb_image_resize.h"

#define TINYEXR_IMPLEMENTATION
#include "ext/tinyexr.h"

#endif

#ifndef _WIN32
#pragma GCC diagnostic pop
#endif

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR COLOR UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Convert between CIE XYZ and xyY
vec3f xyz_to_xyY(vec3f xyz) {
    if (xyz == zero3f) return zero3f;
    return {xyz.x / (xyz.x + xyz.y + xyz.z), xyz.y / (xyz.x + xyz.y + xyz.z),
        xyz.y};
}
// Convert between CIE XYZ and xyY
vec3f xyY_to_xyz(vec3f xyY) {
    if (xyY.y == 0) return zero3f;
    return {xyY.x * xyY.z / xyY.y, xyY.z, (1 - xyY.x - xyY.y) * xyY.z / xyY.y};
}
// Convert between CIE XYZ and RGB
vec3f xyz_to_rgb(vec3f xyz) {
    // from http://www.brucelindbloom.com/index.html?Eqn_RGB_to_XYZ.html
    if (xyz == zero3f) return zero3f;
    return {+3.2404542f * xyz.x - 1.5371385f * xyz.y - 0.4985314f * xyz.z,
        -0.9692660f * xyz.x + 1.8760108f * xyz.y + 0.0415560f * xyz.z,
        +0.0556434f * xyz.x - 0.2040259f * xyz.y + 1.0572252f * xyz.z};
}
// Convert between CIE XYZ and RGB
vec3f rgb_to_xyz(vec3f rgb) {
    // from http://www.brucelindbloom.com/index.html?Eqn_RGB_to_XYZ.html
    if (rgb == zero3f) return zero3f;
    return {0.4124564f * rgb.x + 0.3575761f * rgb.y + 0.1804375f * rgb.z,
        0.2126729f * rgb.x + 0.7151522f * rgb.y + 0.0721750f * rgb.z,
        0.0193339f * rgb.x + 0.1191920f * rgb.y + 0.9503041f * rgb.z};
}

// Convert HSV to RGB
vec3f hsv_to_rgb(vec3f hsv) {
    // from Imgui.cpp
    auto h = hsv.x, s = hsv.y, v = hsv.z;
    if (hsv.y == 0.0f) return {v, v, v};

    h = fmodf(h, 1.0f) / (60.0f / 360.0f);
    int i = (int)h;
    float f = h - (float)i;
    float p = v * (1.0f - s);
    float q = v * (1.0f - s * f);
    float t = v * (1.0f - s * (1.0f - f));

    switch (i) {
        case 0: return {v, t, p};
        case 1: return {q, v, p};
        case 2: return {p, v, t};
        case 3: return {p, q, v};
        case 4: return {t, p, v};
        case 5: return {v, p, q};
        default: return {v, p, q};
    }
}
vec3f rgb_to_hsv(vec3f rgb) {
    // from Imgui.cpp
    auto r = rgb.x, g = rgb.y, b = rgb.z;
    float K = 0.f;
    if (g < b) {
        std::swap(g, b);
        K = -1.f;
    }
    if (r < g) {
        std::swap(r, g);
        K = -2.f / 6.f - K;
    }

    float chroma = r - (g < b ? g : b);
    return {
        fabsf(K + (g - b) / (6.f * chroma + 1e-20f)), chroma / (r + 1e-20f), r};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGE TYPE
// -----------------------------------------------------------------------------
namespace ygl {

image4f make_image4f(int width, int height, vec4f c) {
    return image4f{
        width, height, std::vector<vec4f>{(size_t)width * (size_t)height, c}};
}

image4b make_image4b(int width, int height, vec4b c) {
    return image4b{
        width, height, std::vector<vec4b>{(size_t)width * (size_t)height, c}};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Conversion between linear and gamma-encoded images.
image4f gamma_to_linear(const image4f& srgb, float gamma) {
    auto lin = make_image4f(srgb.width, srgb.height);
    for (auto i = 0; i < srgb.pxl.size(); i++) {
        xyz(lin.pxl[i]) = gamma_to_linear(xyz(srgb.pxl[i]), gamma);
        lin.pxl[i].w = srgb.pxl[i].w;
    }
    return lin;
}
image4f linear_to_gamma(const image4f& lin, float gamma) {
    auto srgb = make_image4f(lin.width, lin.height);
    for (auto i = 0; i < lin.pxl.size(); i++) {
        xyz(srgb.pxl[i]) = linear_to_gamma(xyz(lin.pxl[i]), gamma);
        srgb.pxl[i].w = lin.pxl[i].w;
    }
    return srgb;
}

// Conversion from/to floats.
image4f byte_to_float(const image4b& bt) {
    auto fl = make_image4f(bt.width, bt.height);
    for (auto i = 0; i < bt.pxl.size(); i++)
        fl.pxl[i] = byte_to_float(bt.pxl[i]);
    return fl;
}
image4b float_to_byte(const image4f& fl) {
    auto bt = make_image4b(fl.width, fl.height);
    for (auto i = 0; i < fl.pxl.size(); i++)
        bt.pxl[i] = float_to_byte(fl.pxl[i]);
    return bt;
}

vec3f filmic_tonemap(vec3f hdr) {
    // https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
    auto x = hdr;
    // x *= 0.6; // brings it back to ACES range
    x = (x * (2.51f * x + vec3f{0.03f, 0.03f, 0.03f})) /
        (x * (2.43f * x + vec3f{0.59f, 0.59f, 0.59f}) +
            vec3f{0.14f, 0.14f, 0.14f});
    return x;
}

// Tone mapping HDR to LDR images.
image4f expose_image(const image4f& hdr, float exposure) {
    if (!exposure) return hdr;
    auto ldr = make_image4f(hdr.width, hdr.height);
    auto scale = pow(2.0f, exposure);
    for (auto i = 0; i < hdr.pxl.size(); i++) {
        xyz(ldr.pxl[i]) = xyz(hdr.pxl[i]) * scale;
        ldr.pxl[i].w = hdr.pxl[i].w;
    }
    return ldr;
}

// Tone mapping HDR to LDR images.
image4f filmic_tonemap_image(const image4f& hdr) {
    auto ldr = make_image4f(hdr.width, hdr.height);
    for (auto i = 0; i < hdr.pxl.size(); i++) {
        xyz(ldr.pxl[i]) = filmic_tonemap(xyz(hdr.pxl[i]));
        ldr.pxl[i].w = hdr.pxl[i].w;
    }
    return ldr;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGEIO
// -----------------------------------------------------------------------------
namespace ygl {

// Pfm load
float* load_pfm(const char* filename, int* w, int* h, int* nc, int req) {
    auto split = [](const std::string& str) {
        auto ret = std::vector<std::string>();
        if (str.empty()) return ret;
        auto lpos = (size_t)0;
        while (lpos != str.npos) {
            auto pos = str.find_first_of(" \t\n\r", lpos);
            if (pos != str.npos) {
                if (pos > lpos) ret.push_back(str.substr(lpos, pos - lpos));
                lpos = pos + 1;
            } else {
                if (lpos < str.size()) ret.push_back(str.substr(lpos));
                lpos = pos;
            }
        }
        return ret;
    };

    auto f = fopen(filename, "rb");
    if (!f) return nullptr;

    // buffer
    char buf[256];
    auto toks = std::vector<std::string>();

    // read magic
    if (!fgets(buf, 256, f)) return nullptr;
    toks = split(buf);
    if (toks[0] == "Pf")
        *nc = 1;
    else if (toks[0] == "PF")
        *nc = 3;
    else
        return nullptr;

    // read w, h
    if (!fgets(buf, 256, f)) return nullptr;
    toks = split(buf);
    *w = atoi(toks[0].c_str());
    *h = atoi(toks[1].c_str());

    // read scale
    if (!fgets(buf, 256, f)) return nullptr;
    toks = split(buf);
    auto s = atof(toks[0].c_str());

    // read the data (flip y)
    auto npixels = (*w) * (*h);
    auto nvalues = (*w) * (*h) * (*nc);
    auto nrow = (*w) * (*nc);
    auto pixels = new float[nvalues];
    for (auto j = *h - 1; j >= 0; j--) {
        if (fread(pixels + j * nrow, sizeof(float), nrow, f) != nrow) {
            delete[] pixels;
            return nullptr;
        }
    }

    // done reading
    fclose(f);

    // endian conversion
    if (s > 0) {
        for (auto i = 0; i < nvalues; ++i) {
            auto dta = (uint8_t*)(pixels + i);
            std::swap(dta[0], dta[3]);
            std::swap(dta[1], dta[2]);
        }
    }

    // scale
    auto scl = (s > 0) ? s : -s;
    if (scl != 1) {
        for (auto i = 0; i < nvalues; i++) pixels[i] *= scl;
    }

    // proper number of channels
    if (!req || *nc == req) return pixels;

    // pack into channels
    if (req < 0 || req > 4) {
        delete[] pixels;
        return nullptr;
    }
    auto cpixels = new float[req * npixels];
    for (auto i = 0; i < npixels; i++) {
        auto vp = pixels + i * (*nc);
        auto cp = cpixels + i * req;
        if (*nc == 1) {
            switch (req) {
                case 1: cp[0] = vp[0]; break;
                case 2:
                    cp[0] = vp[0];
                    cp[1] = vp[0];
                    break;
                case 3:
                    cp[0] = vp[0];
                    cp[1] = vp[0];
                    cp[2] = vp[0];
                    break;
                case 4:
                    cp[0] = vp[0];
                    cp[1] = vp[0];
                    cp[2] = vp[0];
                    cp[3] = 1;
                    break;
            }
        } else {
            switch (req) {
                case 1: cp[0] = vp[0]; break;
                case 2:
                    cp[0] = vp[0];
                    cp[1] = vp[1];
                    break;
                case 3:
                    cp[0] = vp[0];
                    cp[1] = vp[1];
                    cp[2] = vp[2];
                    break;
                case 4:
                    cp[0] = vp[0];
                    cp[1] = vp[1];
                    cp[2] = vp[2];
                    cp[3] = 1;
                    break;
            }
        }
    }
    delete[] pixels;
    return cpixels;
}

// save pfm
bool save_pfm(const char* filename, int w, int h, int nc, const float* pixels) {
    auto f = fopen(filename, "wb");
    if (!f) return false;

    fprintf(f, "%s\n", (nc == 1) ? "Pf" : "PF");
    fprintf(f, "%d %d\n", w, h);
    fprintf(f, "-1\n");
    if (nc == 1 || nc == 3) {
        fwrite(pixels, sizeof(float), w * h * nc, f);
    } else {
        for (auto i = 0; i < w * h; i++) {
            auto vz = 0.0f;
            auto v = pixels + i * nc;
            fwrite(v + 0, sizeof(float), 1, f);
            fwrite(v + 1, sizeof(float), 1, f);
            if (nc == 2)
                fwrite(&vz, sizeof(float), 1, f);
            else
                fwrite(v + 2, sizeof(float), 1, f);
        }
    }

    fclose(f);

    return true;
}

// check hdr extensions
bool is_hdr_filename(const std::string& filename) {
    auto ext = path_extension(filename);
    return ext == ".hdr" || ext == ".exr" || ext == ".pfm";
}

// Loads an hdr image.
image4f load_image(const std::string& filename, float ldr_gamma) {
    auto ext = path_extension(filename);
    auto img = image4f();
    if (ext == ".exr") {
        auto pixels = (vec4f*)nullptr;
        if (LoadEXR((float**)&pixels, &img.width, &img.height, filename.c_str(),
                nullptr) < 0)
            throw std::runtime_error("could not load image " + filename);
        if (!pixels)
            throw std::runtime_error("could not load image " + filename);

        img.pxl = std::vector<vec4f>(pixels, pixels + img.width * img.height);
        free(pixels);
    } else if (ext == ".pfm") {
        auto ncomp = 0;
        auto pixels = (vec4f*)load_pfm(
            filename.c_str(), &img.width, &img.height, &ncomp, 4);
        if (!pixels)
            throw std::runtime_error("could not load image " + filename);

        img.pxl = std::vector<vec4f>(pixels, pixels + img.width * img.height);
        free(pixels);
    } else if (ext == ".hdr") {
        auto ncomp = 0;
        auto pixels = (vec4f*)stbi_loadf(
            filename.c_str(), &img.width, &img.height, &ncomp, 4);
        if (!pixels)
            throw std::runtime_error("could not load image " + filename);

        img.pxl = std::vector<vec4f>(pixels, pixels + img.width * img.height);
        free(pixels);
    } else {
        auto img8 = image4b();
        auto ncomp = 0;
        auto pixels = (vec4b*)stbi_load(
            filename.c_str(), &img8.width, &img8.height, &ncomp, 4);
        if (!pixels)
            throw std::runtime_error("could not load image " + filename);
        img8.pxl = std::vector<vec4b>(pixels, pixels + img.width * img.height);
        free(pixels);
        img = gamma_to_linear(byte_to_float(img8), ldr_gamma);
    }
    return img;
}

// Saves an hdr image.
void save_image(
    const std::string& filename, const image4f& img, float ldr_gamma) {
    if (path_extension(filename) == ".png") {
        auto ldr = float_to_byte(linear_to_gamma(img));
        if (!stbi_write_png(filename.c_str(), img.width, img.height, 4,
                (byte*)ldr.pxl.data(), img.width * 4))
            throw std::runtime_error("could not save image " + filename);
    } else if (path_extension(filename) == ".jpg") {
        auto ldr = float_to_byte(linear_to_gamma(img));
        if (!stbi_write_jpg(filename.c_str(), img.width, img.height, 4,
                (byte*)ldr.pxl.data(), 75))
            throw std::runtime_error("could not save image " + filename);
    } else if (path_extension(filename) == ".tga") {
        auto ldr = float_to_byte(linear_to_gamma(img));
        if (!stbi_write_tga(filename.c_str(), img.width, img.height, 4,
                (byte*)ldr.pxl.data()))
            throw std::runtime_error("could not save image " + filename);
    } else if (path_extension(filename) == ".bmp") {
        auto ldr = float_to_byte(linear_to_gamma(img));
        if (!stbi_write_bmp(filename.c_str(), img.width, img.height, 4,
                (byte*)ldr.pxl.data()))
            throw std::runtime_error("could not save image " + filename);
    } else if (path_extension(filename) == ".hdr") {
        if (!stbi_write_hdr(filename.c_str(), img.width, img.height, 4,
                (float*)img.pxl.data()))
            throw std::runtime_error("could not save image " + filename);
    } else if (path_extension(filename) == ".pfm") {
        if (!save_pfm(filename.c_str(), img.width, img.height, 4,
                (float*)img.pxl.data()))
            throw std::runtime_error("could not save image " + filename);
    } else if (path_extension(filename) == ".exr") {
        if (!SaveEXR((float*)img.pxl.data(), img.width, img.height, 4,
                filename.c_str()))
            throw std::runtime_error("could not save image " + filename);
    } else {
        throw std::runtime_error(
            "unsupported image format " + path_extension(filename));
    }
}

// Loads an hdr image.
image4f load_image_from_memory(
    const byte* data, int data_size, float ldr_gamma) {
    auto img = image4f();
    stbi_ldr_to_hdr_gamma(ldr_gamma);
    auto ncomp = 0;
    auto pixels = (vec4f*)stbi_loadf_from_memory(
        data, data_size, &img.width, &img.height, &ncomp, 4);
    stbi_ldr_to_hdr_gamma(2.2f);
    if (!pixels) throw std::runtime_error("could not decode image from memory");
    img.pxl = std::vector<vec4f>(pixels, pixels + img.width * img.height);
    delete pixels;
    return img;
}

// Loads an hdr image.
image4b load_image4b(const std::string& filename, float ldr_gamma) {
    auto img = image4b();
    stbi_hdr_to_ldr_gamma(ldr_gamma);
    auto ncomp = 0;
    auto pixels =
        (vec4b*)stbi_load(filename.c_str(), &img.width, &img.height, &ncomp, 4);
    stbi_hdr_to_ldr_gamma(2.2f);
    if (!pixels) throw std::runtime_error("could not decode image from memory");
    img.pxl = std::vector<vec4b>(pixels, pixels + img.width * img.height);
    delete pixels;
    return img;
}

// Saves an hdr image.
void save_image4b(
    const std::string& filename, const image4b& img, float ldr_gamma) {
    if (path_extension(filename) == ".png") {
        if (!stbi_write_png(filename.c_str(), img.width, img.height, 4,
                (byte*)img.pxl.data(), img.width * 4))
            throw std::runtime_error("could not save image " + filename);
    } else if (path_extension(filename) == ".jpg") {
        if (!stbi_write_jpg(filename.c_str(), img.width, img.height, 4,
                (byte*)img.pxl.data(), 75))
            throw std::runtime_error("could not save image " + filename);
    } else if (path_extension(filename) == ".tga") {
        if (!stbi_write_tga(filename.c_str(), img.width, img.height, 4,
                (byte*)img.pxl.data()))
            throw std::runtime_error("could not save image " + filename);
    } else if (path_extension(filename) == ".bmp") {
        if (!stbi_write_bmp(filename.c_str(), img.width, img.height, 4,
                (byte*)img.pxl.data()))
            throw std::runtime_error("could not save image " + filename);
    } else if (path_extension(filename) == ".hdr") {
        auto hdr = linear_to_gamma(byte_to_float(img), ldr_gamma);
        if (!stbi_write_hdr(filename.c_str(), img.width, img.height, 4,
                (float*)hdr.pxl.data()))
            throw std::runtime_error("could not save image " + filename);
    } else if (path_extension(filename) == ".pfm") {
        auto hdr = linear_to_gamma(byte_to_float(img), ldr_gamma);
        if (!save_pfm(filename.c_str(), img.width, img.height, 4,
                (float*)hdr.pxl.data()))
            throw std::runtime_error("could not save image " + filename);
    } else if (path_extension(filename) == ".exr") {
        auto hdr = linear_to_gamma(byte_to_float(img), ldr_gamma);
        if (!SaveEXR((float*)hdr.pxl.data(), img.width, img.height, 4,
                filename.c_str()))
            throw std::runtime_error("could not save image " + filename);
    } else {
        throw std::runtime_error(
            "unsupported image format " + path_extension(filename));
    }
}

// Loads an hdr image.
image4b load_image4b_from_memory(
    const byte* data, int data_size, int& width, int& height, float ldr_gamma) {
    auto img = image4b();
    stbi_hdr_to_ldr_gamma(ldr_gamma);
    auto ncomp = 0;
    auto pixels = (vec4b*)stbi_load_from_memory(
        data, data_size, &img.width, &img.height, &ncomp, 4);
    stbi_hdr_to_ldr_gamma(2.2f);
    if (!pixels) throw std::runtime_error("could not decode image from memory");
    img.pxl = std::vector<vec4b>(pixels, pixels + width * height);
    delete pixels;
    return img;
}

// Resize image.
std::vector<vec4f> resize_image(int width, int height,
    const std::vector<vec4f>& img, int res_width, int res_height) {
    auto res_img = std::vector<vec4f>(res_width * res_height);
    stbir_resize_float_generic((float*)img.data(), width, height,
        sizeof(vec4f) * width, (float*)res_img.data(), res_width, res_height,
        sizeof(vec4f) * res_width, 4, 3, 0, STBIR_EDGE_CLAMP,
        STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
    return res_img;
}

// Resize image.
image4f resize_image(const image4f& img, int res_width, int res_height) {
    if (!res_width && !res_height) throw std::runtime_error("bad image size");
    if (!res_width)
        res_width = (int)round(img.width * (res_height / (float)img.height));
    if (!res_height)
        res_height = (int)round(img.height * (res_width / (float)img.width));
    auto res_img = make_image4f(res_width, res_height);
    stbir_resize_float_generic((float*)img.pxl.data(), img.width, img.height,
        sizeof(vec4f) * img.width, (float*)res_img.pxl.data(), res_width,
        res_height, sizeof(vec4f) * res_width, 4, 3, 0, STBIR_EDGE_CLAMP,
        STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
    return res_img;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGE EXAMPLES
// -----------------------------------------------------------------------------
namespace ygl {

// Make a grid image
image4f make_grid_image(int width, int height, int tiles, vec4f c0, vec4f c1) {
    auto img = make_image4f(width, height);
    auto tile = width / tiles;
    for (int j = 0; j < width; j++) {
        for (int i = 0; i < height; i++) {
            auto c = i % tile == 0 || i % tile == tile - 1 || j % tile == 0 ||
                     j % tile == tile - 1;
            img.pxl[i + j * width] = (c) ? c0 : c1;
        }
    }
    return img;
}

// Make a checkerboard image
image4f make_checker_image(
    int width, int height, int tiles, vec4f c0, vec4f c1) {
    auto img = make_image4f(width, height);
    auto tile = width / tiles;
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto c = (i / tile + j / tile) % 2 == 0;
            img.pxl[i + j * width] = (c) ? c0 : c1;
        }
    }
    return img;
}

// Make an image with bumps and dimples.
image4f make_bumpdimple_image(int width, int height, int tiles) {
    auto img = make_image4f(width, height);
    auto tile = width / tiles;
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto c = (i / tile + j / tile) % 2 == 0;
            auto ii = i % tile - tile / 2, jj = j % tile - tile / 2;
            auto r =
                sqrt(float(ii * ii + jj * jj)) / sqrt(float(tile * tile) / 4);
            auto h = 0.5f;
            if (r < 0.5f) { h += (c) ? (0.5f - r) : -(0.5f - r); }
            img.pxl[i + j * width] = {h, h, h, 1};
        }
    }
    return img;
}

// Make a uv colored grid
image4f make_ramp_image(int width, int height, vec4f c0, vec4f c1) {
    auto img = make_image4f(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto u = (float)i / (float)width;
            img.pxl[i + j * width] = c0 * (1 - u) + c1 * u;
        }
    }
    return img;
}

// Make a gamma ramp image
image4f make_gammaramp_imagef(int width, int height) {
    auto img = make_image4f(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto u = j / float(height - 1);
            if (i < width / 3) u = pow(u, 2.2f);
            if (i > (width * 2) / 3) u = pow(u, 1 / 2.2f);
            img.pxl[i + j * width] = {u, u, u, 1};
        }
    }
    return img;
}

// Make an image color with red/green in the [0,1] range. Helpful to visualize
// uv texture coordinate application.
image4f make_uv_image(int width, int height) {
    auto img = make_image4f(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            img.pxl[i + j * width] = {
                i / (float)(width - 1), j / (float)(height - 1), 0, 1};
        }
    }
    return img;
}

// Make a uv colored grid
image4f make_uvgrid_image(int width, int height, int tiles, bool colored) {
    auto img = make_image4f(width, height);
    auto tile = width / tiles;
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto ii = i / tile, jj = j / tile;
            auto ww = width / tile, hh = height / tile;
            auto ph =
                (((256 / (ww * hh)) * (ii + jj * ww) - 64 + 256) % 256) / 360.f;
            auto pv = 0.5f;
            auto ps = 0.8f;
            if (i % (tile / 2) && j % (tile / 2)) {
                if ((i / tile + j / tile) % 2)
                    pv += 0.05f;
                else
                    pv -= 0.05f;
            } else {
                pv = 0.8f;
                ps = 0.2f;
            }
            auto rgb = (colored) ? hsv_to_rgb({ph, ps, pv}) : vec3f{pv, pv, pv};
            img.pxl[i + (height - j - 1) * width] = {rgb.x, rgb.y, rgb.z, 1};
        }
    }
    return img;
}

// Comvert a bump map to a normal map.
image4f bump_to_normal_map(const image4f& img, float scale) {
    auto norm = make_image4f(img.width, img.height);
    for (int j = 0; j < img.height; j++) {
        for (int i = 0; i < img.width; i++) {
            auto i1 = (i + 1) % img.width, j1 = (j + 1) % img.height;
            auto p00 = img.pxl[i + j * img.width],
                 p10 = img.pxl[i1 + j * img.width],
                 p01 = img.pxl[i + j1 * img.width];
            auto g00 = (float(p00.x) + float(p00.y) + float(p00.z)) / (3 * 255);
            auto g01 = (float(p01.x) + float(p01.y) + float(p01.z)) / (3 * 255);
            auto g10 = (float(p10.x) + float(p10.y) + float(p10.z)) / (3 * 255);
            auto n = vec3f{scale * (g00 - g10), scale * (g00 - g01), 1.0f};
            n.y = -n.y;  // make green pointing up, even if y axis points down
            n = normalize(n) * 0.5f + vec3f{0.5f, 0.5f, 0.5f};
            norm.pxl[i + j * img.width] = {n.x, n.y, n.z, 1};
        }
    }
    return norm;
}

// Implementation of sunsky modified heavily from pbrt
image4f make_sunsky_image(int width, int height, float thetaSun,
    float turbidity, bool has_sun, vec3f ground_albedo) {
    auto wSun = vec3f{0, cos(thetaSun), sin(thetaSun)};

    // sunSpectralRad =  ComputeAttenuatedSunlight(thetaS, turbidity);
    auto sunAngularRadius = 9.35e-03f / 2;  // Wikipedia
    auto thetaS = thetaSun;

    auto t1 = thetaSun, t2 = thetaSun * thetaSun,
         t3 = thetaSun * thetaSun * thetaSun;
    auto T = turbidity;
    auto T2 = turbidity * turbidity;

    auto zenith_xyY = vec3f{
        (+0.00165f * t3 - 0.00374f * t2 + 0.00208f * t1 + 0) * T2 +
            (-0.02902f * t3 + 0.06377f * t2 - 0.03202f * t1 + 0.00394f) * T +
            (+0.11693f * t3 - 0.21196f * t2 + 0.06052f * t1 + 0.25885f),
        (+0.00275f * t3 - 0.00610f * t2 + 0.00316f * t1 + 0) * T2 +
            (-0.04214f * t3 + 0.08970f * t2 - 0.04153f * t1 + 0.00515f) * T +
            (+0.15346f * t3 - 0.26756f * t2 + 0.06669f * t1 + 0.26688f),
        1000 * (4.0453f * T - 4.9710f) *
                tan((4.0f / 9.0f - T / 120.0f) * (pi - 2 * t1)) -
            .2155f * T + 2.4192f};

    auto perez_A_xyY = vec3f{-0.01925f * T - 0.25922f, -0.01669f * T - 0.26078f,
        +0.17872f * T - 1.46303f};
    auto perez_B_xyY = vec3f{-0.06651f * T + 0.00081f, -0.09495f * T + 0.00921f,
        -0.35540f * T + 0.42749f};
    auto perez_C_xyY = vec3f{-0.00041f * T + 0.21247f, -0.00792f * T + 0.21023f,
        -0.02266f * T + 5.32505f};
    auto perez_D_xyY = vec3f{-0.06409f * T - 0.89887f, -0.04405f * T - 1.65369f,
        +0.12064f * T - 2.57705f};
    auto perez_E_xyY = vec3f{-0.00325f * T + 0.04517f, -0.01092f * T + 0.05291f,
        -0.06696f * T + 0.37027f};

    auto perez_f = [thetaS](float A, float B, float C, float D, float E,
                       float theta, float gamma, float zenith) -> float {
        auto den = ((1 + A * exp(B)) *
                    (1 + C * exp(D * thetaS) + E * cos(thetaS) * cos(thetaS)));
        auto num = ((1 + A * exp(B / cos(theta))) *
                    (1 + C * exp(D * gamma) + E * cos(gamma) * cos(gamma)));
        return zenith * num / den;
    };

    auto sky = [&perez_f, perez_A_xyY, perez_B_xyY, perez_C_xyY, perez_D_xyY,
                   perez_E_xyY, zenith_xyY](auto theta, auto gamma) -> vec3f {
        auto x = perez_f(perez_A_xyY.x, perez_B_xyY.x, perez_C_xyY.x,
            perez_D_xyY.x, perez_E_xyY.x, theta, gamma, zenith_xyY.x);
        auto y = perez_f(perez_A_xyY.y, perez_B_xyY.y, perez_C_xyY.y,
            perez_D_xyY.y, perez_E_xyY.y, theta, gamma, zenith_xyY.y);
        auto Y = perez_f(perez_A_xyY.z, perez_B_xyY.z, perez_C_xyY.z,
            perez_D_xyY.z, perez_E_xyY.z, theta, gamma, zenith_xyY.z);
        return xyz_to_rgb(xyY_to_xyz({x, y, Y})) / 10000.0f;
    };

    // compute sun luminance
    // TODO: how this relates to zenith intensity?
    auto sun_ko = vec3f{0.48f, 0.75f, 0.14f};
    auto sun_kg = vec3f{0.1f, 0.0f, 0.0f};
    auto sun_kwa = vec3f{0.02f, 0.0f, 0.0f};
    auto sun_sol = vec3f{20000.0f, 27000.0f, 30000.0f};
    auto sun_lambda = vec3f{680, 530, 480};
    auto sun_beta = 0.04608365822050f * turbidity - 0.04586025928522f;
    auto sun_m =
        1.0f / (cos(thetaSun) + 0.000940f * pow(1.6386f - thetaSun, -1.253f));

    auto sun_le = zero3f;
    for (auto i = 0; i < 3; i++) {
        auto tauR =
            exp(-sun_m * 0.008735f * pow((&sun_lambda.x)[i] / 1000, -4.08f));
        auto tauA =
            exp(-sun_m * sun_beta * pow((&sun_lambda.x)[i] / 1000, -1.3f));
        auto tauO = exp(-sun_m * (&sun_ko.x)[i] * .35f);
        auto tauG = exp(-1.41f * (&sun_kg.x)[i] * sun_m /
                        pow(1 + 118.93f * (&sun_kg.x)[i] * sun_m, 0.45f));
        auto tauWA =
            exp(-0.2385f * (&sun_kwa.x)[i] * 2.0f * sun_m /
                pow(1 + 20.07f * (&sun_kwa.x)[i] * 2.0f * sun_m, 0.45f));
        (&sun_le.x)[i] = (&sun_sol.x)[i] * tauR * tauA * tauO * tauG * tauWA;
    }

    auto sun = [has_sun, sunAngularRadius, sun_le](auto theta, auto gamma) {
        return (has_sun && gamma < sunAngularRadius) ? sun_le / 10000.0f :
                                                       zero3f;
    };

    auto img = make_image4f(width, height, {0, 0, 0});
    for (auto j = 0; j < height / 2; j++) {
        auto theta = pi * ((j + 0.5f) / height);
        theta = clamp(theta, 0.0f, pi / 2 - flt_eps);
        for (int i = 0; i < width; i++) {
            auto phi = 2 * pi * (float(i + 0.5f) / width);
            auto w =
                vec3f{cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta)};
            auto gamma = acos(clamp(dot(w, wSun), -1.0f, 1.0f));
            auto col = sky(theta, gamma) + sun(theta, gamma);
            img.pxl[i + j * width] = {col.x, col.y, col.z, 1};
        }
    }

    if (ground_albedo != zero3f) {
        auto ground = zero3f;
        for (auto j = 0; j < height / 2; j++) {
            auto theta = pi * ((j + 0.5f) / height);
            for (int i = 0; i < width; i++) {
                auto pxl = img.pxl[i + j * width];
                auto le = vec3f{pxl.x, pxl.y, pxl.z};
                auto angle = sin(theta) * 4 * pi / (width * height);
                ground += le * (ground_albedo / pi) * cos(theta) * angle;
            }
        }
        for (auto j = height / 2; j < height; j++) {
            for (int i = 0; i < width; i++) {
                img.pxl[i + j * width] = {ground.x, ground.y, ground.z, 1};
            }
        }
    }

    return img;
}

// Make an image of multiple lights.
image4f make_lights_image(int width, int height, vec3f le, int nlights,
    float langle, float lwidth, float lheight) {
    auto img = make_image4f(width, height, {0, 0, 0, 1});
    for (auto j = 0; j < height / 2; j++) {
        auto theta = pi * ((j + 0.5f) / height);
        theta = clamp(theta, 0.0f, pi / 2 - flt_eps);
        if (fabs(theta - langle) > lheight / 2) continue;
        for (int i = 0; i < width; i++) {
            auto phi = 2 * pi * (float(i + 0.5f) / width);
            auto inlight = false;
            for (auto l = 0; l < nlights; l++) {
                auto lphi = 2 * pi * (l + 0.5f) / nlights;
                inlight = inlight || fabs(phi - lphi) < lwidth / 2;
            }
            img.pxl[i + j * width] = {le.x, le.y, le.z, 1};
        }
    }
    return img;
}

// Make a noise image. Wrap works only if both width and height are powers of
// two.
image4f make_noise_image(int width, int height, float scale, bool wrap) {
    auto wrap3i = (wrap) ? vec3i{width, height, 2} : zero3i;
    auto img = make_image4f(width, height);
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto p = vec3f{i / (float)width, j / (float)height, 0.5f} * scale;
            auto g = perlin_noise(p, wrap3i);
            g = clamp(0.5f + 0.5f * g, 0.0f, 1.0f);
            img.pxl[i + j * width] = {g, g, g, 1};
        }
    }
    return img;
}

// Make a noise image. Wrap works only if both width and height are powers of
// two.
image4f make_fbm_image(int width, int height, float scale, float lacunarity,
    float gain, int octaves, bool wrap) {
    auto wrap3i = (wrap) ? vec3i{width, height, 2} : zero3i;
    auto img = make_image4f(width, height);
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto p = vec3f{i / (float)width, j / (float)height, 0.5f} * scale;
            auto g = perlin_fbm_noise(p, lacunarity, gain, octaves, wrap3i);
            g = clamp(0.5f + 0.5f * g, 0.0f, 1.0f);
            img.pxl[i + j * width] = {g, g, g, 1};
        }
    }
    return img;
}

// Make a noise image. Wrap works only if both width and height are powers of
// two.
image4f make_ridge_image(int width, int height, float scale, float lacunarity,
    float gain, float offset, int octaves, bool wrap) {
    auto wrap3i = (wrap) ? vec3i{width, height, 2} : zero3i;
    auto img = make_image4f(width, height);
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto p = vec3f{i / (float)width, j / (float)height, 0.5f} * scale;
            auto g = perlin_ridge_noise(
                p, lacunarity, gain, offset, octaves, wrap3i);
            g = clamp(g, 0.0f, 1.0f);
            img.pxl[i + j * width] = {g, g, g, 1};
        }
    }
    return img;
}

// Make a noise image. Wrap works only if both width and height are powers of
// two.
image4f make_turbulence_image(int width, int height, float scale,
    float lacunarity, float gain, int octaves, bool wrap) {
    auto wrap3i = (wrap) ? vec3i{width, height, 2} : zero3i;
    auto img = make_image4f(width, height);
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto p = vec3f{i / (float)width, j / (float)height, 0.5f} * scale;
            auto g =
                perlin_turbulence_noise(p, lacunarity, gain, octaves, wrap3i);
            g = clamp(g, 0.0f, 1.0f);
            img.pxl[i + j * width] = {g, g, g, 1};
        }
    }
    return img;
}

}  // namespace ygl
