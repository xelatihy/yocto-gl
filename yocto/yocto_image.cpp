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

// Element-wise float to byte conversion.
vec4b float_to_byte(const vec4f& a) {
    return {float_to_byte(a.x), float_to_byte(a.y), float_to_byte(a.z),
        float_to_byte(a.w)};
}
vec4f byte_to_float(const vec4b& a) {
    return {byte_to_float(a.x), byte_to_float(a.y), byte_to_float(a.z),
        byte_to_float(a.w)};
}

// Conversion between linear and gamma-encoded images.
vec4f gamma_to_linear(const vec4f& srgb, float gamma) {
    return {gamma_to_linear(srgb.x, gamma), gamma_to_linear(srgb.y, gamma),
        gamma_to_linear(srgb.z, gamma), srgb.w};
}
vec4f linear_to_gamma(const vec4f& lin, float gamma) {
    return {linear_to_gamma(lin.x, gamma), linear_to_gamma(lin.y, gamma),
        linear_to_gamma(lin.z, gamma), lin.w};
}
vec3f gamma_to_linear(const vec3f& srgb, float gamma) {
    return {gamma_to_linear(srgb.x, gamma), gamma_to_linear(srgb.y, gamma),
        gamma_to_linear(srgb.z, gamma)};
}
vec3f linear_to_gamma(const vec3f& lin, float gamma) {
    return {linear_to_gamma(lin.x, gamma), linear_to_gamma(lin.y, gamma),
        linear_to_gamma(lin.z, gamma)};
}

// Approximate luminance estimation
float luminance(const vec3f& a) { return (a.x + a.y + a.z) / 3; }

// Convert between CIE XYZ and xyY
vec3f xyz_to_xyY(const vec3f& xyz) {
    if (xyz == zero3f) return zero3f;
    return {xyz.x / (xyz.x + xyz.y + xyz.z), xyz.y / (xyz.x + xyz.y + xyz.z),
        xyz.y};
}
// Convert between CIE XYZ and xyY
vec3f xyY_to_xyz(const vec3f& xyY) {
    if (xyY.y == 0) return zero3f;
    return {xyY.x * xyY.z / xyY.y, xyY.z, (1 - xyY.x - xyY.y) * xyY.z / xyY.y};
}
// Convert between CIE XYZ and RGB
vec3f xyz_to_rgb(const vec3f& xyz) {
    // from http://www.brucelindbloom.com/index.html?Eqn_RGB_to_XYZ.html
    if (xyz == zero3f) return zero3f;
    return {+3.2404542f * xyz.x - 1.5371385f * xyz.y - 0.4985314f * xyz.z,
        -0.9692660f * xyz.x + 1.8760108f * xyz.y + 0.0415560f * xyz.z,
        +0.0556434f * xyz.x - 0.2040259f * xyz.y + 1.0572252f * xyz.z};
}
// Convert between CIE XYZ and RGB
vec3f rgb_to_xyz(const vec3f& rgb) {
    // from http://www.brucelindbloom.com/index.html?Eqn_RGB_to_XYZ.html
    if (rgb == zero3f) return zero3f;
    return {0.4124564f * rgb.x + 0.3575761f * rgb.y + 0.1804375f * rgb.z,
        0.2126729f * rgb.x + 0.7151522f * rgb.y + 0.0721750f * rgb.z,
        0.0193339f * rgb.x + 0.1191920f * rgb.y + 0.9503041f * rgb.z};
}

// Convert HSV to RGB
vec3f hsv_to_rgb(const vec3f& hsv) {
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
vec3f rgb_to_hsv(const vec3f& rgb) {
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
// IMPLEMENTATION FOR IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Conversion between linear and gamma-encoded images.
std::vector<vec4f> gamma_to_linear(
    const std::vector<vec4f>& srgb, float gamma) {
    auto lin = std::vector<vec4f>(srgb.size());
    for (auto i = 0; i < srgb.size(); i++)
        lin[i] = gamma_to_linear(srgb[i], gamma);
    return lin;
}
std::vector<vec4f> linear_to_gamma(const std::vector<vec4f>& lin, float gamma) {
    auto srgb = std::vector<vec4f>(lin.size());
    for (auto i = 0; i < lin.size(); i++)
        srgb[i] = linear_to_gamma(lin[i], gamma);
    return srgb;
}
std::vector<vec3f> gamma_to_linear(
    const std::vector<vec3f>& srgb, float gamma) {
    auto lin = std::vector<vec3f>(srgb.size());
    for (auto i = 0; i < srgb.size(); i++)
        lin[i] = gamma_to_linear(srgb[i], gamma);
    return lin;
}
std::vector<vec3f> linear_to_gamma(const std::vector<vec3f>& lin, float gamma) {
    auto srgb = std::vector<vec3f>(lin.size());
    for (auto i = 0; i < lin.size(); i++)
        srgb[i] = linear_to_gamma(lin[i], gamma);
    return srgb;
}
// Conversion between linear and gamma-encoded images.
std::vector<float> gamma_to_linear(
    const std::vector<float>& srgb, int ncomp, float gamma) {
    auto lin = std::vector<float>(srgb.size());
    if (ncomp == 1) {
        for (auto i = 0; i < srgb.size(); i++) {
            lin[i] = gamma_to_linear(srgb[i], gamma);
        }
    } else if (ncomp == 2) {
        for (auto i = 0; i < srgb.size() / 2; i++) {
            lin[i * 2 + 0] = gamma_to_linear(srgb[i * 2 + 0], gamma);
            lin[i * 2 + 1] = srgb[i * 2 + 1];
        }
    } else if (ncomp == 3) {
        for (auto i = 0; i < srgb.size() / 3; i++) {
            lin[i * 3 + 0] = gamma_to_linear(srgb[i * 3 + 0], gamma);
            lin[i * 3 + 1] = gamma_to_linear(srgb[i * 3 + 1], gamma);
            lin[i * 3 + 2] = gamma_to_linear(srgb[i * 3 + 2], gamma);
        }
    } else if (ncomp == 4) {
        for (auto i = 0; i < srgb.size() / 4; i++) {
            lin[i * 4 + 0] = gamma_to_linear(srgb[i * 4 + 0], gamma);
            lin[i * 4 + 1] = gamma_to_linear(srgb[i * 4 + 1], gamma);
            lin[i * 4 + 2] = gamma_to_linear(srgb[i * 4 + 2], gamma);
            lin[i * 4 + 3] = srgb[i * 4 + 3];
        }
    } else {
        throw std::runtime_error("number of components not supported");
    }
    return lin;
}
std::vector<float> linear_to_gamma(
    const std::vector<float>& lin, int ncomp, float gamma) {
    auto srgb = std::vector<float>(lin.size());
    if (ncomp == 1) {
        for (auto i = 0; i < lin.size(); i++) {
            srgb[i] = linear_to_gamma(srgb[i], gamma);
        }
    } else if (ncomp == 2) {
        for (auto i = 0; i < lin.size() / 2; i++) {
            srgb[i * 2 + 0] = linear_to_gamma(lin[i * 2 + 0], gamma);
            srgb[i * 2 + 1] = lin[i * 2 + 1];
        }
    } else if (ncomp == 3) {
        for (auto i = 0; i < lin.size() / 3; i++) {
            srgb[i * 3 + 0] = linear_to_gamma(lin[i * 3 + 0], gamma);
            srgb[i * 3 + 1] = linear_to_gamma(lin[i * 3 + 1], gamma);
            srgb[i * 3 + 2] = linear_to_gamma(lin[i * 3 + 2], gamma);
        }
    } else if (ncomp == 4) {
        for (auto i = 0; i < lin.size() / 4; i++) {
            srgb[i * 4 + 0] = linear_to_gamma(lin[i * 4 + 0], gamma);
            srgb[i * 4 + 1] = linear_to_gamma(lin[i * 4 + 1], gamma);
            srgb[i * 4 + 2] = linear_to_gamma(lin[i * 4 + 2], gamma);
            srgb[i * 4 + 3] = lin[i * 4 + 3];
        }
    } else {
        throw std::runtime_error("number of components not supported");
    }
    return srgb;
}

// Conversion from/to floats.
std::vector<vec4f> byte_to_float(const std::vector<vec4b>& bt) {
    auto fl = std::vector<vec4f>(bt.size());
    for (auto i = 0; i < bt.size(); i++) fl[i] = byte_to_float(bt[i]);
    return fl;
}
std::vector<vec4b> float_to_byte(const std::vector<vec4f>& fl) {
    auto bt = std::vector<vec4b>(fl.size());
    for (auto i = 0; i < fl.size(); i++) bt[i] = float_to_byte(fl[i]);
    return bt;
}
// Conversion from/to floats.
std::vector<float> byte_to_float(const std::vector<byte>& bt) {
    auto fl = std::vector<float>(bt.size());
    for (auto i = 0; i < bt.size(); i++) fl[i] = byte_to_float(bt[i]);
    return fl;
}
std::vector<byte> float_to_byte(const std::vector<float>& fl) {
    auto bt = std::vector<byte>(fl.size());
    for (auto i = 0; i < fl.size(); i++) bt[i] = float_to_byte(fl[i]);
    return bt;
}

// Conversion between different number of channels
std::vector<vec4f> imagef_to_image4f(const std::vector<float>& ax, int ncomp) {
    auto npixels = (int)ax.size() / ncomp;
    auto a4 = std::vector<vec4f>(npixels);
    if (ncomp == 1) {
        for (auto i = 0; i < npixels; i++) { a4[i] = {ax[i], ax[i], ax[i], 1}; }
    } else if (ncomp == 2) {
        for (auto i = 0; i < npixels; i++) {
            a4[i] = {
                ax[i * 2 + 0], ax[i * 2 + 0], ax[i * 2 + 0], ax[i * 2 + 1]};
        }
    } else if (ncomp == 3) {
        for (auto i = 0; i < npixels; i++) {
            a4[i] = {ax[i * 3 + 0], ax[i * 3 + 1], ax[i * 3 + 2], 1};
        }
    } else if (ncomp == 4) {
        for (auto i = 0; i < npixels; i++) {
            a4[i] = {
                ax[i * 4 + 0], ax[i * 4 + 1], ax[i * 4 + 2], ax[i * 4 + 3]};
        }
    } else {
        throw std::runtime_error("number of components not supported");
    }
    return a4;
}
std::vector<float> image4f_to_imagef(const std::vector<vec4f>& a4, int ncomp) {
    auto npixels = (int)a4.size();
    auto ax = std::vector<float>(npixels * ncomp);
    if (ncomp == 1) {
        for (auto i = 0; i < npixels; i++) {
            ax[i] = (a4[i].x + a4[i].y + a4[i].z) / 3.0f;
        }
    } else if (ncomp == 2) {
        for (auto i = 0; i < npixels; i++) {
            ax[i * 2 + 0] = (a4[i].x + a4[i].y + a4[i].z) / 3.0f;
            ax[i * 2 + 1] = a4[i].w;
        }
    } else if (ncomp == 3) {
        for (auto i = 0; i < npixels; i++) {
            ax[i * 3 + 0] = a4[i].x;
            ax[i * 3 + 1] = a4[i].y;
            ax[i * 3 + 2] = a4[i].z;
        }
    } else if (ncomp == 4) {
        for (auto i = 0; i < npixels; i++) {
            ax[i * 4 + 0] = a4[i].x;
            ax[i * 4 + 1] = a4[i].y;
            ax[i * 4 + 2] = a4[i].z;
            ax[i * 4 + 3] = a4[i].w;
        }
    } else {
        throw std::runtime_error("number of components not supported");
    }
    return ax;
}
std::vector<vec4b> imageb_to_image4b(const std::vector<byte>& ax, int ncomp) {
    auto npixels = (int)ax.size() / ncomp;
    auto a4 = std::vector<vec4b>(npixels);
    if (ncomp == 1) {
        for (auto i = 0; i < npixels; i++) {
            a4[i] = {ax[i], ax[i], ax[i], (byte)255};
        }
    } else if (ncomp == 2) {
        for (auto i = 0; i < npixels; i++) {
            a4[i] = {
                ax[i * 2 + 0], ax[i * 2 + 0], ax[i * 2 + 0], ax[i * 2 + 1]};
        }
    } else if (ncomp == 3) {
        for (auto i = 0; i < npixels; i++) {
            a4[i] = {ax[i * 3 + 0], ax[i * 3 + 1], ax[i * 3 + 2], (byte)255};
        }
    } else if (ncomp == 4) {
        for (auto i = 0; i < npixels; i++) {
            a4[i] = {
                ax[i * 4 + 0], ax[i * 4 + 1], ax[i * 4 + 2], ax[i * 4 + 3]};
        }
    } else {
        throw std::runtime_error("number of components not supported");
    }
    return a4;
}
std::vector<byte> image4b_to_imageb(const std::vector<vec4b>& a4, int ncomp) {
    auto npixels = (int)a4.size();
    auto ax = std::vector<byte>(npixels * ncomp);
    if (ncomp == 1) {
        for (auto i = 0; i < npixels; i++) {
            ax[i] = (byte)((a4[i].x + a4[i].y + a4[i].z) / 3.0f);
        }
    } else if (ncomp == 2) {
        for (auto i = 0; i < npixels; i++) {
            ax[i * 2 + 0] = (byte)((a4[i].x + a4[i].y + a4[i].z) / 3.0f);
            ax[i * 2 + 1] = a4[i].w;
        }
    } else if (ncomp == 3) {
        for (auto i = 0; i < npixels; i++) {
            ax[i * 3 + 0] = a4[i].x;
            ax[i * 3 + 1] = a4[i].y;
            ax[i * 3 + 2] = a4[i].z;
        }
    } else if (ncomp == 4) {
        for (auto i = 0; i < npixels; i++) {
            ax[i * 4 + 0] = a4[i].x;
            ax[i * 4 + 1] = a4[i].y;
            ax[i * 4 + 2] = a4[i].z;
            ax[i * 4 + 3] = a4[i].w;
        }
    } else {
        throw std::runtime_error("number of components not supported");
    }
    return ax;
}
// Conversion between different number of channels
std::vector<vec3f> imagef_to_image3f(const std::vector<float>& ax, int ncomp) {
    auto npixels = (int)ax.size() / ncomp;
    auto a3 = std::vector<vec3f>(npixels);
    if (ncomp == 1) {
        for (auto i = 0; i < npixels; i++) { a3[i] = {ax[i], ax[i], ax[i]}; }
    } else if (ncomp == 2) {
        for (auto i = 0; i < npixels; i++) {
            a3[i] = {
                ax[i * 2 + 0], ax[i * 2 + 0], ax[i * 2 + 0]};
        }
    } else if (ncomp == 3) {
        for (auto i = 0; i < npixels; i++) {
            a3[i] = {ax[i * 3 + 0], ax[i * 3 + 1], ax[i * 3 + 2]};
        }
    } else if (ncomp == 4) {
        for (auto i = 0; i < npixels; i++) {
            a3[i] = {
                ax[i * 4 + 0], ax[i * 4 + 1], ax[i * 4 + 2]};
        }
    } else {
        throw std::runtime_error("number of components not supported");
    }
    return a3;
}
std::vector<float> image3f_to_imagef(const std::vector<vec3f>& a3, int ncomp) {
    auto npixels = (int)a3.size();
    auto ax = std::vector<float>(npixels * ncomp);
    if (ncomp == 1) {
        for (auto i = 0; i < npixels; i++) {
            ax[i] = (a3[i].x + a3[i].y + a3[i].z) / 3.0f;
        }
    } else if (ncomp == 2) {
        for (auto i = 0; i < npixels; i++) {
            ax[i * 2 + 0] = (a3[i].x + a3[i].y + a3[i].z) / 3.0f;
            ax[i * 2 + 1] = 1;
        }
    } else if (ncomp == 3) {
        for (auto i = 0; i < npixels; i++) {
            ax[i * 3 + 0] = a3[i].x;
            ax[i * 3 + 1] = a3[i].y;
            ax[i * 3 + 2] = a3[i].z;
        }
    } else if (ncomp == 4) {
        for (auto i = 0; i < npixels; i++) {
            ax[i * 4 + 0] = a3[i].x;
            ax[i * 4 + 1] = a3[i].y;
            ax[i * 4 + 2] = a3[i].z;
            ax[i * 4 + 3] = 1;
        }
    } else {
        throw std::runtime_error("number of components not supported");
    }
    return ax;
}

vec3f filmic_tonemap(const vec3f& hdr) {
    // https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
    auto x = hdr;
    // x *= 0.6; // brings it back to ACES range
    x = (x * (2.51f * x + vec3f{0.03f, 0.03f, 0.03f})) /
        (x * (2.43f * x + vec3f{0.59f, 0.59f, 0.59f}) +
            vec3f{0.14f, 0.14f, 0.14f});
    return x;
}

// Tone mapping HDR to LDR images.
std::vector<vec4f> expose_image(const std::vector<vec4f>& hdr, float exposure) {
    if (!exposure) return hdr;
    auto ldr = std::vector<vec4f>(hdr.size());
    auto scale = pow(2.0f, exposure);
    for (auto i = 0; i < hdr.size(); i++)
        ldr[i] = {
            hdr[i].x * scale, hdr[i].y * scale, hdr[i].z * scale, hdr[i].w};
    return ldr;
}

// Tone mapping HDR to LDR images.
std::vector<vec4f> filmic_tonemap_image(const std::vector<vec4f>& hdr) {
    auto ldr = std::vector<vec4f>(hdr.size());
    for (auto i = 0; i < hdr.size(); i++) {
        auto h = filmic_tonemap({hdr[i].x, hdr[i].y, hdr[i].z});
        ldr[i] = {h.x, h.y, h.z, hdr[i].w};
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

// Loads an ldr image.
std::vector<vec4b> load_image4b(
    const std::string& filename, int& width, int& height) {
    auto ncomp = 0;
    auto pixels = load_imageb(filename, width, height, ncomp);
    return imageb_to_image4b(pixels, ncomp);
}

// Loads an hdr image.
std::vector<vec4f> load_image4f(
    const std::string& filename, int& width, int& height, float ldr_gamma) {
    auto ncomp = 0;
    auto pixels = load_imagef(filename, width, height, ncomp);
    return imagef_to_image4f(pixels, ncomp);
}

// Saves an ldr image.
bool save_image4b(const std::string& filename, int width, int height,
    const std::vector<vec4b>& img) {
    auto pxl = image4b_to_imageb(img, 4);
    return save_imageb(filename, width, height, 4, pxl);
}

// Saves an hdr image.
bool save_image4f(const std::string& filename, int width, int height,
    const std::vector<vec4f>& img, float ldr_gamma) {
    auto pxl = image4f_to_imagef(img, 4);
    return save_imagef(filename, width, height, 4, pxl, ldr_gamma);
}

// Loads an ldr image.
std::vector<vec4b> load_image4b_from_memory(
    const byte* data, int data_size, int& width, int& height) {
    auto ncomp = 0;
    auto pixels = load_imageb_from_memory(data, data_size, width, height, ncomp);
    return imageb_to_image4b(pixels, ncomp);
}

// Loads an hdr image.
std::vector<vec4f> load_image4f_from_memory(
    const byte* data, int data_size, int& width, int& height, float ldr_gamma) {
    auto ncomp = 0;
    auto pixels = load_imagef_from_memory(data, data_size, width, height, ncomp);
    return imagef_to_image4f(pixels, ncomp);
}

// Loads an ldr image.
std::vector<byte> load_imageb(
    const std::string& filename, int& width, int& height, int& ncomp) {
    auto pixels = stbi_load(filename.c_str(), &width, &height, &ncomp, 0);
    if (!pixels) return {};
    auto img = std::vector<byte>(pixels, pixels + width * height * ncomp);
    free(pixels);
    return img;
}

// Loads an hdr image.
std::vector<float> load_imagef(const std::string& filename, int& width,
    int& height, int& ncomp, float ldr_gamma) {
    auto ext = path_extension(filename);
    if (ext == ".exr") {
        auto pixels = (float*)nullptr;
        if (LoadEXR((float**)&pixels, &width, &height, filename.c_str(),
                nullptr) < 0)
            return {};
        if (!pixels) return {};
        ncomp = 4;
        auto img = std::vector<float>(pixels, pixels + width * height * ncomp);
        free(pixels);
        return img;
    } else if (ext == ".pfm") {
        auto pixels = load_pfm(filename.c_str(), &width, &height, &ncomp, 0);
        if (!pixels) return {};
        auto img = std::vector<float>(pixels, pixels + width * height * ncomp);
        free(pixels);
        return img;
    } else if (ext == ".hdr") {
        auto pixels = stbi_loadf(filename.c_str(), &width, &height, &ncomp, 0);
        if (!pixels) return {};
        auto img = std::vector<float>(pixels, pixels + width * height * ncomp);
        free(pixels);
        return img;
    } else {
        auto pixels = stbi_load(filename.c_str(), &width, &height, &ncomp, 0);
        if (!pixels) return {};
        auto img8 = std::vector<byte>(pixels, pixels + width * height * ncomp);
        free(pixels);
        auto img = gamma_to_linear(byte_to_float(img8), ncomp, ldr_gamma);
        return img;
    }
}

// Saves an ldr image.
bool save_imageb(const std::string& filename, int width, int height, int ncomp,
    const std::vector<byte>& img) {
    if (path_extension(filename) == ".png") {
        return stbi_write_png(
            filename.c_str(), width, height, ncomp, img.data(), width * ncomp);
    } else if (path_extension(filename) == ".jpg") {
        return stbi_write_jpg(
            filename.c_str(), width, height, ncomp, img.data(), 75);
    } else if (path_extension(filename) == ".tga") {
        return stbi_write_tga(
            filename.c_str(), width, height, ncomp, img.data());
    } else if (path_extension(filename) == ".bmp") {
        return stbi_write_bmp(
            filename.c_str(), width, height, ncomp, img.data());
    } else {
        return false;
    }
}

// Saves an hdr image.
bool save_imagef(const std::string& filename, int width, int height, int ncomp,
    const std::vector<float>& img, float ldr_gamma) {
    if (path_extension(filename) == ".png") {
        auto ldr = float_to_byte(linear_to_gamma(img, ncomp));
        return stbi_write_png(filename.c_str(), width, height, ncomp,
            (byte*)ldr.data(), width * ncomp);
    } else if (path_extension(filename) == ".jpg") {
        auto ldr = float_to_byte(linear_to_gamma(img, ncomp));
        return stbi_write_jpg(
            filename.c_str(), width, height, ncomp, (byte*)ldr.data(), 75);
    } else if (path_extension(filename) == ".tga") {
        auto ldr = float_to_byte(linear_to_gamma(img, ncomp));
        return stbi_write_tga(
            filename.c_str(), width, height, ncomp, (byte*)ldr.data());
    } else if (path_extension(filename) == ".bmp") {
        auto ldr = float_to_byte(linear_to_gamma(img, ncomp));
        return stbi_write_bmp(
            filename.c_str(), width, height, ncomp, (byte*)ldr.data());
    } else if (path_extension(filename) == ".hdr") {
        return stbi_write_hdr(
            filename.c_str(), width, height, ncomp, (float*)img.data());
    } else if (path_extension(filename) == ".pfm") {
        return save_pfm(filename.c_str(), width, height, ncomp, img.data());
    } else if (path_extension(filename) == ".exr") {
        return !SaveEXR(
            (float*)img.data(), width, height, ncomp, filename.c_str());
    } else {
        return false;
    }
}

// Loads an ldr image.
std::vector<byte> load_imageb_from_memory(
    const byte* data, int data_size, int& width, int& height, int& ncomp) {
    auto pixels =
        stbi_load_from_memory(data, data_size, &width, &height, &ncomp, 0);
    if (!pixels) return {};
    auto img = std::vector<byte>(pixels, pixels + width * height * ncomp);
    free(pixels);
    return img;
}

// Loads an hdr image.
std::vector<float> load_imagef_from_memory(const byte* data, int data_size,
    int& width, int& height, int& ncomp, float ldr_gamma) {
    stbi_ldr_to_hdr_gamma(ldr_gamma);
    auto pixels =
        stbi_loadf_from_memory(data, data_size, &width, &height, &ncomp, 0);
    stbi_ldr_to_hdr_gamma(2.2f);
    if (!pixels) return {};
    auto img = std::vector<float>(pixels, pixels + width * height * ncomp);
    delete pixels;
    return img;
}

static const auto resize_filter_map = std::map<resize_filter, stbir_filter>{
    {resize_filter::def, STBIR_FILTER_DEFAULT},
    {resize_filter::box, STBIR_FILTER_BOX},
    {resize_filter::triangle, STBIR_FILTER_TRIANGLE},
    {resize_filter::cubic_spline, STBIR_FILTER_CUBICBSPLINE},
    {resize_filter::catmull_rom, STBIR_FILTER_CATMULLROM},
    {resize_filter::mitchell, STBIR_FILTER_MITCHELL}};

static const auto resize_edge_map = std::map<resize_edge, stbir_edge>{
    {resize_edge::def, STBIR_EDGE_CLAMP},
    {resize_edge::clamp, STBIR_EDGE_CLAMP},
    {resize_edge::reflect, STBIR_EDGE_REFLECT},
    {resize_edge::wrap, STBIR_EDGE_WRAP}, {resize_edge::zero, STBIR_EDGE_ZERO}};

// Resize image.
std::vector<vec4f> resize_image(int width, int height,
    const std::vector<vec4f>& img, int res_width, int res_height,
    resize_filter filter, resize_edge edge, bool premultiplied_alpha) {
    auto res_img = std::vector<vec4f>(res_width * res_height);
    stbir_resize_float_generic((float*)img.data(), width, height,
        sizeof(vec4f) * width, (float*)res_img.data(), res_width, res_height,
        sizeof(vec4f) * res_width, 4, 3,
        (premultiplied_alpha) ? STBIR_FLAG_ALPHA_PREMULTIPLIED : 0,
        resize_edge_map.at(edge), resize_filter_map.at(filter),
        STBIR_COLORSPACE_LINEAR, nullptr);
    return res_img;
}
std::vector<vec4b> resize_image(int width, int height,
    const std::vector<vec4b>& img, int res_width, int res_height,
    resize_filter filter, resize_edge edge, bool premultiplied_alpha) {
    auto res_img = std::vector<vec4b>(res_width * res_height);
    stbir_resize_uint8_generic((unsigned char*)img.data(), width, height,
        sizeof(vec4b) * width, (unsigned char*)res_img.data(), res_width,
        res_height, sizeof(vec4b) * res_width, 4, 3,
        (premultiplied_alpha) ? STBIR_FLAG_ALPHA_PREMULTIPLIED : 0,
        resize_edge_map.at(edge), resize_filter_map.at(filter),
        STBIR_COLORSPACE_LINEAR, nullptr);
    return res_img;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGE EXAMPLES
// -----------------------------------------------------------------------------
namespace ygl {

// Make a grid image
std::vector<vec3f> make_grid_image(
    int width, int height, int tiles, const vec3f& c0, const vec3f& c1) {
    auto img = std::vector<vec3f>(width * height);
    auto tile = width / tiles;
    for (int j = 0; j < width; j++) {
        for (int i = 0; i < height; i++) {
            auto c = i % tile == 0 || i % tile == tile - 1 || j % tile == 0 ||
                     j % tile == tile - 1;
            img[i + j * width] = (c) ? c0 : c1;
        }
    }
    return img;
}

// Make a checkerboard image
std::vector<vec3f> make_checker_image(
    int width, int height, int tiles, const vec3f& c0, const vec3f& c1) {
    auto img = std::vector<vec3f>(width * height);
    auto tile = width / tiles;
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto c = (i / tile + j / tile) % 2 == 0;
            img[i + j * width] = (c) ? c0 : c1;
        }
    }
    return img;
}

// Make an image with bumps and dimples.
std::vector<vec3f> make_bumpdimple_image(int width, int height, int tiles) {
    auto img = std::vector<vec3f>(width * height);
    auto tile = width / tiles;
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto c = (i / tile + j / tile) % 2 == 0;
            auto ii = i % tile - tile / 2, jj = j % tile - tile / 2;
            auto r =
                sqrt(float(ii * ii + jj * jj)) / sqrt(float(tile * tile) / 4);
            auto h = 0.5f;
            if (r < 0.5f) { h += (c) ? (0.5f - r) : -(0.5f - r); }
            img[i + j * width] = {h, h, h};
        }
    }
    return img;
}

// Make a uv colored grid
std::vector<vec3f> make_ramp_image(
    int width, int height, const vec3f& c0, const vec3f& c1) {
    auto img = std::vector<vec3f>(width * height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto u = (float)i / (float)width;
            img[i + j * width] = c0 * (1 - u) + c1 * u;
        }
    }
    return img;
}

// Make a gamma ramp image
std::vector<vec3f> make_gammaramp_imagef(int width, int height) {
    auto img = std::vector<vec3f>(width * height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto u = j / float(height - 1);
            if (i < width / 3) u = pow(u, 2.2f);
            if (i > (width * 2) / 3) u = pow(u, 1 / 2.2f);
            img[i + j * width] = {u, u, u};
        }
    }
    return img;
}

// Make an image color with red/green in the [0,1] range. Helpful to visualize
// uv texture coordinate application.
std::vector<vec3f> make_uv_image(int width, int height) {
    auto img = std::vector<vec3f>(width * height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            img[i + j * width] = {
                i / (float)(width - 1), j / (float)(height - 1), 0};
        }
    }
    return img;
}

// Make a uv colored grid
std::vector<vec3f> make_uvgrid_image(
    int width, int height, int tiles, bool colored) {
    auto img = std::vector<vec3f>(width * height);
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
            img[i + (height - j - 1) * width] = rgb;
        }
    }
    return img;
}

// Comvert a bump map to a normal map.
std::vector<vec3f> bump_to_normal_map(
    int width, int height, const std::vector<vec3f>& img, float scale) {
    auto norm = std::vector<vec3f>(img.size());
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto i1 = (i + 1) % width, j1 = (j + 1) % height;
            auto p00 = img[i + j * width], p10 = img[i1 + j * width],
                 p01 = img[i + j1 * width];
            auto g00 = (float(p00.x) + float(p00.y) + float(p00.z)) / (3 * 255);
            auto g01 = (float(p01.x) + float(p01.y) + float(p01.z)) / (3 * 255);
            auto g10 = (float(p10.x) + float(p10.y) + float(p10.z)) / (3 * 255);
            auto n = vec3f{scale * (g00 - g10), scale * (g00 - g01), 1.0f};
            n.y = -n.y;  // make green pointing up, even if y axis points down
            n = normalize(n) * 0.5f + vec3f{0.5f, 0.5f, 0.5f};
            norm[i + j * width] = n;
        }
    }
    return norm;
}

// Implementation of sunsky modified heavily from pbrt
std::vector<vec3f> make_sunsky_image(int width, int height, float thetaSun,
    float turbidity, bool has_sun, const vec3f& ground_albedo) {
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

    auto img = std::vector<vec3f>(width * height, {0, 0, 0});
    for (auto j = 0; j < height / 2; j++) {
        auto theta = pi * ((j + 0.5f) / height);
        theta = clamp(theta, 0.0f, pi / 2 - flt_eps);
        for (int i = 0; i < width; i++) {
            auto phi = 2 * pi * (float(i + 0.5f) / width);
            auto w =
                vec3f{cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta)};
            auto gamma = acos(clamp(dot(w, wSun), -1.0f, 1.0f));
            auto col = sky(theta, gamma) + sun(theta, gamma);
            img[i + j * width] = {col.x, col.y, col.z};
        }
    }

    if (ground_albedo != zero3f) {
        auto ground = zero3f;
        for (auto j = 0; j < height / 2; j++) {
            auto theta = pi * ((j + 0.5f) / height);
            for (int i = 0; i < width; i++) {
                auto pxl = img[i + j * width];
                auto le = vec3f{pxl.x, pxl.y, pxl.z};
                auto angle = sin(theta) * 4 * pi / (width * height);
                ground += le * (ground_albedo / pi) * cos(theta) * angle;
            }
        }
        for (auto j = height / 2; j < height; j++) {
            for (int i = 0; i < width; i++) {
                img[i + j * width] = ground;
            }
        }
    }

    return img;
}

// Make an image of multiple lights.
std::vector<vec3f> make_lights_image(int width, int height, const vec3f& le,
    int nlights, float langle, float lwidth, float lheight) {
    auto img = std::vector<vec3f>(width * height, {0, 0, 0});
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
            img[i + j * width] = le;
        }
    }
    return img;
}

// Make a noise image. Wrap works only if both width and height are powers of
// two.
std::vector<vec3f> make_noise_image(
    int width, int height, float scale, bool wrap) {
    auto wrap3i = (wrap) ? vec3i{width, height, 2} : zero3i;
    auto img = std::vector<vec3f>(width * height);
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto p = vec3f{i / (float)width, j / (float)height, 0.5f} * scale;
            auto g = perlin_noise(p, wrap3i);
            g = clamp(0.5f + 0.5f * g, 0.0f, 1.0f);
            img[i + j * width] = {g, g, g};
        }
    }
    return img;
}

// Make a noise image. Wrap works only if both width and height are powers of
// two.
std::vector<vec3f> make_fbm_image(int width, int height, float scale,
    float lacunarity, float gain, int octaves, bool wrap) {
    auto wrap3i = (wrap) ? vec3i{width, height, 2} : zero3i;
    auto img = std::vector<vec3f>(width * height);
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto p = vec3f{i / (float)width, j / (float)height, 0.5f} * scale;
            auto g = perlin_fbm_noise(p, lacunarity, gain, octaves, wrap3i);
            g = clamp(0.5f + 0.5f * g, 0.0f, 1.0f);
            img[i + j * width] = {g, g, g};
        }
    }
    return img;
}

// Make a noise image. Wrap works only if both width and height are powers of
// two.
std::vector<vec3f> make_ridge_image(int width, int height, float scale,
    float lacunarity, float gain, float offset, int octaves, bool wrap) {
    auto wrap3i = (wrap) ? vec3i{width, height, 2} : zero3i;
    auto img = std::vector<vec3f>(width * height);
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto p = vec3f{i / (float)width, j / (float)height, 0.5f} * scale;
            auto g = perlin_ridge_noise(
                p, lacunarity, gain, offset, octaves, wrap3i);
            g = clamp(g, 0.0f, 1.0f);
            img[i + j * width] = {g, g, g};
        }
    }
    return img;
}

// Make a noise image. Wrap works only if both width and height are powers of
// two.
std::vector<vec3f> make_turbulence_image(int width, int height, float scale,
    float lacunarity, float gain, int octaves, bool wrap) {
    auto wrap3i = (wrap) ? vec3i{width, height, 2} : zero3i;
    auto img = std::vector<vec3f>(width * height);
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto p = vec3f{i / (float)width, j / (float)height, 0.5f} * scale;
            auto g =
                perlin_turbulence_noise(p, lacunarity, gain, octaves, wrap3i);
            g = clamp(g, 0.0f, 1.0f);
            img[i + j * width] = {g, g, g};
        }
    }
    return img;
}

}  // namespace ygl
