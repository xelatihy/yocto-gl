//
// Implementation for Yocto/ImageIO.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#include "yocto_imageio.h"
#include "yocto_utils.h"

#include <climits>
#include <cstdlib>
#include <memory>

#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wunused-variable"
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#endif

// #ifndef __clang_analyzer__

#include "ext/json.hpp"
#include "ext/stb_image.h"
#include "ext/stb_image_write.h"
#include "ext/tinyexr.h"

// #endif

#ifndef _WIN32
#pragma GCC diagnostic pop
#endif

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

using nlohmann::json;
using std::unique_ptr;

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON SUPPORT
// -----------------------------------------------------------------------------
namespace yocto {

// Load a JSON object
inline void load_json(const string& filename, json& js) {
    auto text = ""s;
    load_text(filename, text);
    js = json::parse(text);
}

// Save a JSON object
inline void save_json(const string& filename, const json& js) {
    // we have to use streams here since the json library is faster with them
    save_text(filename, js.dump(4));
}

template <typename T, int N>
inline void to_json(json& js, const vec<T, N>& val) {
    nlohmann::to_json(js, (const std::array<T, N>&)val);
}
template <typename T, int N>
inline void from_json(const json& js, vec<T, N>& val) {
    nlohmann::from_json(js, (std::array<T, N>&)val);
}

template <typename T, int N>
inline void to_json(json& js, const frame<T, N>& val) {
    nlohmann::to_json(js, (const std::array<T, N*(N + 1)>&)val);
}
template <typename T, int N>
inline void from_json(const json& js, frame<T, N>& val) {
    nlohmann::from_json(js, (std::array<T, N*(N + 1)>&)val);
}

template <typename T, int N, int M>
inline void to_json(json& js, const mat<T, N, M>& val) {
    nlohmann::to_json(js, (const std::array<T, N * M>&)val);
}
template <typename T, int N, int M>
inline void from_json(const json& js, mat<T, N, M>& val) {
    nlohmann::from_json(js, (std::array<T, N * M>&)val);
}

template <typename T, int N>
inline void to_json(json& js, const bbox<T, N>& val) {
    nlohmann::from_json(js, (std::array<T, N * 2>&)val);
}
template <typename T, int N>
inline void from_json(const json& js, bbox<T, N>& val) {
    nlohmann::to_json(js, (const std::array<T, N * 2>&)val);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGEIO
// -----------------------------------------------------------------------------
namespace yocto {

// Split a string
vector<string> split_string(const string& str) {
    auto ret = vector<string>();
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
}

// Pfm load
float* load_pfm(const char* filename, int* w, int* h, int* nc, int req) {
    auto fs = fopen(filename, "rb");
    if (!fs) return nullptr;
    auto fs_guard = unique_ptr<FILE, void (*)(FILE*)>{
        fs, [](FILE* f) { fclose(f); }};

    // buffer
    char buffer[4096];
    auto toks = vector<string>();

    // read magic
    if (!fgets(buffer, sizeof(buffer), fs)) return nullptr;
    toks = split_string(buffer);
    if (toks[0] == "Pf")
        *nc = 1;
    else if (toks[0] == "PF")
        *nc = 3;
    else
        return nullptr;

    // read w, h
    if (!fgets(buffer, sizeof(buffer), fs)) return nullptr;
    toks = split_string(buffer);
    *w   = atoi(toks[0].c_str());
    *h   = atoi(toks[1].c_str());

    // read scale
    if (!fgets(buffer, sizeof(buffer), fs)) return nullptr;
    toks   = split_string(buffer);
    auto s = atof(toks[0].c_str());

    // read the data (flip y)
    auto npixels = (size_t)(*w) * (size_t)(*h);
    auto nvalues = npixels * (size_t)(*nc);
    auto nrow    = (size_t)(*w) * (size_t)(*nc);
    auto pixels  = unique_ptr<float[]>(new float[nvalues]);
    for (auto j = *h - 1; j >= 0; j--) {
        if (fread(pixels.get() + j * nrow, sizeof(float), nrow, fs) != nrow)
            return nullptr;
    }

    // endian conversion
    if (s > 0) {
        for (auto i = 0; i < nvalues; ++i) {
            auto dta = (uint8_t*)(pixels.get() + i);
            swap(dta[0], dta[3]);
            swap(dta[1], dta[2]);
        }
    }

    // scale
    auto scl = (s > 0) ? s : -s;
    if (scl != 1) {
        for (auto i = 0; i < nvalues; i++) pixels[i] *= scl;
    }

    // proper number of channels
    if (!req || *nc == req) return pixels.release();

    // pack into channels
    if (req < 0 || req > 4) {
        return nullptr;
    }
    auto cpixels = unique_ptr<float[]>(new float[req * npixels]);
    for (auto i = 0ull; i < npixels; i++) {
        auto vp = pixels.get() + i * (*nc);
        auto cp = cpixels.get() + i * req;
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
    return cpixels.release();
}

// save pfm
bool save_pfm(const char* filename, int w, int h, int nc, const float* pixels) {
    auto fs = fopen(filename, "wb");
    if (!fs) return false;
    auto fs_guard = unique_ptr<FILE, void (*)(FILE*)>{
        fs, [](FILE* f) { fclose(f); }};

    if (fprintf(fs, "%s\n", (nc == 1) ? "Pf" : "PF") < 0) return false;
    if (fprintf(fs, "%d %d\n", w, h) < 0) return false;
    if (fprintf(fs, "-1\n") < 0) return false;
    if (nc == 1 || nc == 3) {
        if (fwrite(pixels, sizeof(float), w * h * nc, fs) != w * h * nc)
            return false;
    } else {
        for (auto i = 0; i < w * h; i++) {
            auto vz = 0.0f;
            auto v  = pixels + i * nc;
            if (fwrite(v + 0, sizeof(float), 1, fs) != 1) return false;
            if (fwrite(v + 1, sizeof(float), 1, fs) != 1) return false;
            if (nc == 2) {
                if (fwrite(&vz, sizeof(float), 1, fs) != 1) return false;
            } else {
                if (fwrite(v + 2, sizeof(float), 1, fs) != 1) return false;
            }
        }
    }

    return true;
}

// Pnm load
byte* load_pnm(const char* filename, int* w, int* h, int* nc, int req) {
    auto fs = fopen(filename, "rb");
    if (!fs) return nullptr;
    auto fs_guard = unique_ptr<FILE, void (*)(FILE*)>{
        fs, [](FILE* f) { fclose(f); }};

    // read magic
    char magic[2];
    if (fscanf(fs, "%c%c", magic + 0, magic + 1) != 2) return nullptr;
    if (magic[0] == 'P' && magic[1] == '2')
        *nc = 1;
    else if (magic[0] == 'P' && magic[1] == '3')
        *nc = 3;
    else
        return nullptr;

    // read w, h, nc
    if (fscanf(fs, "%d %d", w, h) != 2) return nullptr;

    // read max
    auto max = 0;
    if (fscanf(fs, "%d", &max) != 1) return nullptr;
    if (max > 255) return nullptr;

    // read the data (flip y)
    auto npixels = (size_t)(*w) * (size_t)(*h);
    auto nvalues = npixels * (size_t)(*nc);
    auto pixels  = unique_ptr<byte[]>(new byte[nvalues]);
    for (auto i = 0; i < nvalues; i++) {
        if (fscanf(fs, "%hhu", &pixels[i]) != 1) return nullptr;
    }

    // proper number of channels
    if (!req || *nc == req) return pixels.release();

    // pack into channels
    if (req < 0 || req > 4) {
        return nullptr;
    }
    auto cpixels = unique_ptr<byte[]>(new byte[req * npixels]);
    for (auto i = 0ull; i < npixels; i++) {
        auto vp = pixels.get() + i * (*nc);
        auto cp = cpixels.get() + i * req;
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
                    cp[3] = (byte)255;
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
                    cp[3] = (byte)255;
                    break;
            }
        }
    }
    return cpixels.release();
}

// Pnm load
byte* load_pnm_from_string(const char* data, int* w, int* h, int* nc, int req) {
    // read magic
    auto offset = 0;
    char magic[256];
    if (sscanf(data, "%s%n", magic, &offset) != 1) return nullptr;
    if (magic == "P2"s)
        *nc = 1;
    else if (magic == "P3"s)
        *nc = 3;
    else
        return nullptr;

    // read w, h, nc
    data += offset + 1;
    if (sscanf(data, "%d %d%n", w, h, &offset) != 2) return nullptr;

    // read max
    data += offset + 1;
    auto max = 0;
    if (sscanf(data, "%d%n", &max, &offset) != 1) return nullptr;
    if (max > 255) return nullptr;

    // read the data (flip y)
    auto npixels = (size_t)(*w) * (size_t)(*h);
    auto nvalues = npixels * (size_t)(*nc);
    auto pixels  = unique_ptr<byte[]>(new byte[nvalues]);
    for (auto i = 0; i < nvalues; i++) {
        data += offset + 1;
        if (sscanf(data, "%hhu%n", &pixels[i], &offset) != 1) return nullptr;
    }

    // proper number of channels
    if (!req || *nc == req) return pixels.release();

    // pack into channels
    if (req < 0 || req > 4) {
        return nullptr;
    }
    auto cpixels = unique_ptr<byte[]>(new byte[req * npixels]);
    for (auto i = 0ull; i < npixels; i++) {
        auto vp = pixels.get() + i * (*nc);
        auto cp = cpixels.get() + i * req;
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
                    cp[3] = (byte)255;
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
                    cp[3] = (byte)255;
                    break;
            }
        }
    }
    return cpixels.release();
}

// save pnm
bool save_pnm(const char* filename, int w, int h, int nc, const byte* pixels) {
    auto fs = fopen(filename, "wb");
    if (!fs) return false;
    auto fs_guard = unique_ptr<FILE, void (*)(FILE*)>{
        fs, [](FILE* f) { fclose(f); }};

    if (fprintf(fs, "%s\n", (nc == 1) ? "P2" : "P3") < 0) return false;
    if (fprintf(fs, "%d %d\n", w, h) < 0) return false;
    if (fprintf(fs, "255\n") < 0) return false;
    for (auto j = 0; j < h; j++) {
        for (auto i = 0; i < w; i++) {
            auto v = pixels + (j * w + i) * nc;
            if (nc == 1 || nc == 2) {
                if (fprintf(fs, "%d ", (int)v[0]) < 0) return false;
            } else {
                if (fprintf(fs, "%d %d %d ", (int)v[0], (int)v[1], (int)v[2]) <
                    0)
                    return false;
            }
        }
        if (fprintf(fs, "\n") < 0) return false;
    }

    return true;
}

// load pfm image
template <int N>
void load_pfm_image(const string& filename, image<vec<float, N>>& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = load_pfm(filename.c_str(), &width, &height, &ncomp, N);
    if (!pixels) {
        throw imageio_error("error loading image " + filename);
    }
    img = image{{width, height}, (const vec<float, N>*)pixels};
    delete[] pixels;
}
template <int N>
void save_pfm_image(const string& filename, const image<vec<float, N>>& img) {
    if (!save_pfm(filename.c_str(), img.size().x, img.size().y, N,
            (float*)img.data())) {
        throw imageio_error("error saving image " + filename);
    }
}

// load pfm image
template <int N>
void load_pnm_image(const string& filename, image<vec<byte, N>>& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = load_pnm(filename.c_str(), &width, &height, &ncomp, N);
    if (!pixels) {
        throw imageio_error("error loading image " + filename);
    }
    img = image{{width, height}, (const vec<byte, N>*)pixels};
    delete[] pixels;
}
template <int N>
void save_pnm_image(const string& filename, const image<vec<byte, N>>& img) {
    if (!save_pnm(filename.c_str(), img.size().x, img.size().y, N,
            (byte*)img.data())) {
        throw imageio_error("error saving image " + filename);
    }
}
template <int N>
void load_pnm_image_from_string(const char* data, image<vec<byte, N>>& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = load_pnm_from_string(data, &width, &height, &ncomp, N);
    if (!pixels) {
        throw imageio_error("error loading image from string");
    }
    img = image{{width, height}, (const vec<byte, N>*)pixels};
    delete[] pixels;
}

// load exr image weith tiny exr
static const char* get_tinyexr_error(int error) {
    switch (error) {
        case TINYEXR_ERROR_INVALID_MAGIC_NUMBER: return "INVALID_MAGIC_NUMBER";
        case TINYEXR_ERROR_INVALID_EXR_VERSION: return "INVALID_EXR_VERSION";
        case TINYEXR_ERROR_INVALID_ARGUMENT: return "INVALID_ARGUMENT";
        case TINYEXR_ERROR_INVALID_DATA: return "INVALID_DATA";
        case TINYEXR_ERROR_INVALID_FILE: return "INVALID_FILE";
        // case TINYEXR_ERROR_INVALID_PARAMETER: return "INVALID_PARAMETER";
        case TINYEXR_ERROR_CANT_OPEN_FILE: return "CANT_OPEN_FILE";
        case TINYEXR_ERROR_UNSUPPORTED_FORMAT: return "UNSUPPORTED_FORMAT";
        case TINYEXR_ERROR_INVALID_HEADER: return "INVALID_HEADER";
        default: throw imageio_error("unknown tinyexr error");
    }
}

template <int N>
void load_exr_image(const string& filename, image<vec<float, N>>& img) {
    // TODO
    if (N != 4) throw runtime_error("bad number of channels");
    auto width = 0, height = 0;
    auto pixels = (float*)nullptr;
    if (auto error = LoadEXR(
            &pixels, &width, &height, filename.c_str(), nullptr);
        error < 0) {
        throw imageio_error("error loading image " + filename + "("s +
                            get_tinyexr_error(error) + ")"s);
    }
    if (!pixels) {
        throw imageio_error("error loading image " + filename);
    }
    img = image{{width, height}, (const vec<float, N>*)pixels};
    free(pixels);
}
template <int N>
void save_exr_image(const string& filename, const image<vec<float, N>>& img) {
    // TODO
    if (N != 4) throw runtime_error("bad number of channels");
    if (!SaveEXR((float*)img.data(), img.size().x, img.size().y, N,
            filename.c_str())) {
        throw imageio_error("error saving image " + filename);
    }
}

// load an image using stbi library
template <int N>
void load_stb_image(const string& filename, image<vec<byte, N>>& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_load(filename.c_str(), &width, &height, &ncomp, N);
    if (!pixels) {
        throw imageio_error("error loading image " + filename);
    }
    img = image{{width, height}, (const vec<byte, N>*)pixels};
    free(pixels);
}
template <int N>
void load_stb_image(const string& filename, image<vec<float, N>>& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_loadf(filename.c_str(), &width, &height, &ncomp, N);
    if (!pixels) {
        throw imageio_error("error loading image " + filename);
    }
    img = image{{width, height}, (const vec<float, N>*)pixels};
    free(pixels);
}

// save an image with stbi
template <int N>
void save_png_image(const string& filename, const image<vec<byte, N>>& img) {
    if (!stbi_write_png(filename.c_str(), img.size().x, img.size().y, N,
            img.data(), img.size().x * 4)) {
        throw imageio_error("error saving image " + filename);
    }
}
template <int N>
void save_jpg_image(const string& filename, const image<vec<byte, N>>& img) {
    if (!stbi_write_jpg(
            filename.c_str(), img.size().x, img.size().y, 4, img.data(), 75)) {
        throw imageio_error("error saving image " + filename);
    }
}
template <int N>
void save_tga_image(const string& filename, const image<vec<byte, N>>& img) {
    if (!stbi_write_tga(
            filename.c_str(), img.size().x, img.size().y, 4, img.data())) {
        throw imageio_error("error saving image " + filename);
    }
}
template <int N>
void save_bmp_image(const string& filename, const image<vec<byte, N>>& img) {
    if (!stbi_write_bmp(
            filename.c_str(), img.size().x, img.size().y, 4, img.data())) {
        throw imageio_error("error saving image " + filename);
    }
}
template <int N>
void save_hdr_image(const string& filename, const image<vec<float, N>>& img) {
    if (!stbi_write_hdr(filename.c_str(), img.size().x, img.size().y, 4,
            (float*)img.data())) {
        throw imageio_error("error saving image " + filename);
    }
}

// load an image using stbi library
template <int N>
void load_stb_image_from_memory(
    const byte* data, int data_size, image<vec<byte, N>>& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_load_from_memory(
        data, data_size, &width, &height, &ncomp, 4);
    if (!pixels) {
        throw imageio_error("error loading in-memory image");
    }
    img = image{{width, height}, (const vec<byte, N>*)pixels};
    free(pixels);
}
template <int N>
void load_stb_image_from_memory(
    const byte* data, int data_size, image<vec<float, N>>& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_loadf_from_memory(
        data, data_size, &width, &height, &ncomp, 4);
    if (!pixels) {
        throw imageio_error("error loading in-memory image {}");
    }
    img = image{{width, height}, (const vec<float, N>*)pixels};
    free(pixels);
}

void apply_json_procedural(const json& js, image<vec4f>& img) {
    auto type   = js.value("type", ""s);
    auto width  = js.value("width", 1024);
    auto height = js.value("height", 1024);
    if (type == "sky" && width < height * 2) width = height * 2;
    img.resize({width, height});
    if (type == "") {
        img = image{{width, height}, zero4f};
    } else if (type == "grid") {
        make_grid_image(img, js.value("tile", 8),
            js.value("c0", vec4f{0.2f, 0.2f, 0.2f, 1}),
            js.value("c1", vec4f{0.5f, 0.5f, 0.5f, 1}));
    } else if (type == "checker") {
        make_checker_image(img, js.value("tile", 8),
            js.value("c0", vec4f{0.2f, 0.2f, 0.2f, 1}),
            js.value("c1", vec4f{0.5f, 0.5f, 0.5f, 1}));
    } else if (type == "bump") {
        make_bumpdimple_image(img, js.value("tile", 8),
            js.value("c0", vec4f{0, 0, 0, 1}),
            js.value("c1", vec4f{1, 1, 1, 1}));
    } else if (type == "uvramp") {
        make_uvramp_image(img);
    } else if (type == "gammaramp") {
        make_gammaramp_image(img, js.value("c0", vec4f{0, 0, 0, 1}),
            js.value("c1", vec4f{1, 1, 1, 1}));
    } else if (type == "blackbodyramp") {
        make_blackbodyramp_image(img);
    } else if (type == "uvgrid") {
        make_uvgrid_image(img);
    } else if (type == "sky") {
        make_sunsky_image(img, js.value("sun_angle", pif / 4),
            js.value("turbidity", 3.0f), js.value("has_sun", false),
            js.value("sun_intensity", 1.0f), js.value("sun_temperature", 0.0f),
            js.value("ground_albedo", vec3f{0.7f, 0.7f, 0.7f}));
    } else if (type == "noise") {
        make_noise_image(img, js.value("c0", vec4f{0, 0, 0, 1}),
            js.value("c1", vec4f{1, 1, 1, 1}), js.value("scale", 1.0f),
            js.value("wrap", true));
    } else if (type == "fbm") {
        make_fbm_image(img, js.value("c0", vec4f{0, 0, 0, 1}),
            js.value("c1", vec4f{1, 1, 1, 1}), js.value("scale", 1.0f),
            js.value("lacunarity", 2.0f), js.value("gain", 0.5f),
            js.value("octaves", 6), js.value("wrap", true));
    } else if (type == "ridge") {
        make_ridge_image(img, js.value("c0", vec4f{0, 0, 0, 1}),
            js.value("c1", vec4f{1, 1, 1, 1}), js.value("scale", 1.0f),
            js.value("lacunarity", 2.0f), js.value("gain", 0.5f),
            js.value("offset", 1.0f), js.value("octaves", 6),
            js.value("wrap", true));
    } else if (type == "turbulence") {
        make_turbulence_image(img, js.value("c0", vec4f{0, 0, 0, 1}),
            js.value("c1", vec4f{1, 1, 1, 1}), js.value("scale", 1.0f),
            js.value("lacunarity", 2.0f), js.value("gain", 0.5f),
            js.value("octaves", 6), js.value("wrap", true));
    } else if (type == "montage") {
        auto sub_imgs = vector<image<vec4f>>(js.at("images").size());
        for (auto i = 0; i < sub_imgs.size(); i++) {
            apply_json_procedural(js.at("images").at(i), sub_imgs.at(i));
        }
        auto size = zero2i;
        for (auto& sub_img : sub_imgs) {
            size.x += sub_img.size().x;
            size.y = max(size.y, sub_img.size().y);
        }
        img.resize(size);
        auto pos = 0;
        for (auto& sub_img : sub_imgs) {
            set_image_region(img, sub_img, {pos, 0});
            pos += sub_img.size().x;
        }
    } else {
        throw std::invalid_argument("unknown image type" + type);
    }
    if (js.value("border", false)) {
        add_image_border(img, js.value("border_width", 2),
            js.value("border_color", vec4f{0, 0, 0, 1}));
    }
    if (js.value("bump_to_normal", false)) {
        auto buffer = img;
        bump_to_normal_map(img, buffer, js.value("bump_scale", 1.0f));
    }
}

void apply_json_procedural(const json& js, image<vec4b>& img) {
    auto imgf = image<vec4f>{};
    apply_json_procedural(js, imgf);
    auto srgb = js.value("srgb", true);
    if (srgb) {
        auto srgb = imgf;
        linear_to_srgb(srgb, imgf);
        imgf = srgb;
    }
    float_to_byte(img, imgf);
}

// load a JSON image
template <int N>
void load_json_image(const string& filename, image<vec<float, N>>& img) {
    if constexpr (N == 4) {
        auto js = json();
        load_json(filename, js);
        apply_json_procedural(js, img);
    } else {
        auto js = json();
        load_json(filename, js);
        auto img_rgba = image{img.size(), vec<float, 4>{}};
        apply_json_procedural(js, img_rgba);
        rgba_to_color(img, img_rgba);
    }
}
template <int N>
void load_json_image(const string& filename, image<vec<byte, N>>& img) {
    if constexpr (N == 4) {
        auto js = json();
        load_json(filename, js);
        apply_json_procedural(js, img);
    } else {
        auto js = json();
        load_json(filename, js);
        auto img_rgba = image{img.size(), vec<byte, 4>{}};
        apply_json_procedural(js, img_rgba);
        rgba_to_color(img, img_rgba);
    }
}
    
// Check if an image is a preset based on filename.
bool is_image_preset_filename(const string& filename) {
    return get_filename(filename).find("yocto::") == 0;
}
string get_image_preset_type(const string& filename) {
    return get_noextension(get_filename(filename).substr(7));
}

template<int N>
void load_image_preset(const string& filename, image<vec<float, N>>& img) {
    if constexpr(N == 4) {
        img.resize({1024, 1024});
        if(get_image_preset_type(filename) == "images2") img.resize({2048,1024});
        make_image_preset(img, get_image_preset_type(filename));
    } else {
        auto img4 = image<vec<float, 4>>({1024, 1024});
        if(get_image_preset_type(filename) == "images2") img4.resize({2048,1024});
        make_image_preset(img4, get_image_preset_type(filename));
        img.resize(img4.size());
        rgba_to_color(img, img4);
    }
}
template<int N>
void load_image_preset(const string& filename, image<vec<byte, N>>& img) {
    auto imgf = image<vec<float, N>>{};
    load_image_preset(filename,imgf);
    img.resize(imgf.size());
    linear_to_srgb8(img, imgf);
}
    
// check hdr extensions
bool is_hdr_filename(const string& filename) {
    auto ext = get_extension(filename);
    return ext == "hdr" || ext == "exr" || ext == "pfm";
}

// Loads an hdr image.
template <int N>
void load_image(const string& filename, image<vec<float, N>>& img) {
    if (is_image_preset_filename(filename)) {
        return load_image_preset(filename, img);
    }
    auto ext = get_extension(filename);
    if (ext == "exr" || ext == "EXR") {
        load_exr_image(filename, img);
    } else if (ext == "pfm" || ext == "PFM") {
        load_pfm_image(filename, img);
    } else if (ext == "hdr" || ext == "HDR") {
        load_stb_image(filename, img);
    } else if (ext == "png" || ext == "PNG") {
        auto img8 = image<vec<byte, N>>{};
        load_stb_image(filename, img8);
        img.resize(img8.size());
        srgb8_to_linear(img, img8);
    } else if (ext == "jpg" || ext == "JPG") {
        auto img8 = image<vec<byte, N>>{};
        load_stb_image(filename, img8);
        img.resize(img8.size());
        srgb8_to_linear(img, img8);
    } else if (ext == "tga" || ext == "TGA") {
        auto img8 = image<vec<byte, N>>{};
        load_stb_image(filename, img8);
        img.resize(img8.size());
        srgb8_to_linear(img, img8);
    } else if (ext == "bmp" || ext == "BMP") {
        auto img8 = image<vec<byte, N>>{};
        load_stb_image(filename, img8);
        img.resize(img8.size());
        srgb8_to_linear(img, img8);
    } else if (ext == "ppm" || ext == "PPM") {
        auto img8 = image<vec<byte, N>>{};
        load_pnm_image(filename, img8);
        // load_stb_image(filename, img8);
        img.resize(img8.size());
        srgb8_to_linear(img, img8);
    } else if (ext == "pgm" || ext == "PGM") {
        auto img8 = image<vec<byte, N>>{};
        load_pnm_image(filename, img8);
        // load_stb_image(filename, img8);
        img.resize(img8.size());
        srgb8_to_linear(img, img8);
    } else if (ext == "json" || ext == "JSON") {
        load_json_image(filename, img);
    } else {
        throw imageio_error("unsupported image format " + ext);
    }
}

// Saves an hdr image.
template <int N>
void save_image(const string& filename, const image<vec<float, N>>& img) {
    auto ext = get_extension(filename);
    if (ext == "png" || ext == "PNG") {
        auto img8 = image<vec<byte, N>>{img.size()};
        linear_to_srgb8(img8, img);
        save_png_image(filename, img8);
    } else if (ext == "jpg" || ext == "JPG") {
        auto img8 = image<vec<byte, N>>{img.size()};
        linear_to_srgb8(img8, img);
        save_jpg_image(filename, img8);
    } else if (ext == "tga" || ext == "TGA") {
        auto img8 = image<vec<byte, N>>{img.size()};
        linear_to_srgb8(img8, img);
        save_tga_image(filename, img8);
    } else if (ext == "bmp" || ext == "BMP") {
        auto img8 = image<vec<byte, N>>{img.size()};
        linear_to_srgb8(img8, img);
        save_bmp_image(filename, img8);
    } else if (ext == "ppm" || ext == "PPM") {
        auto img8 = image<vec<byte, N>>{img.size()};
        linear_to_srgb8(img8, img);
        save_pnm_image(filename, img8);
    } else if (ext == "pgm" || ext == "PGM") {
        auto img8 = image<vec<byte, N>>{img.size()};
        linear_to_srgb8(img8, img);
        save_pnm_image(filename, img8);
    } else if (ext == "hdr" || ext == "HDR") {
        auto img8 = image<vec<byte, N>>{img.size()};
        linear_to_srgb8(img8, img);
        save_hdr_image(filename, img);
    } else if (ext == "pfm" || ext == "PFM") {
        save_pfm_image(filename, img);
    } else if (ext == "exr" || ext == "EXR") {
        save_exr_image(filename, img);
    } else {
        throw imageio_error("unsupported image format " + ext);
    }
}

// Loads an hdr image.
template <int N>
void load_image_from_memory(
    const byte* data, int data_size, image<vec<float, N>>& img) {
    load_stb_image_from_memory(data, data_size, img);
}

// Loads an hdr image.
template <int N>
void load_image(const string& filename, image<vec<byte, N>>& img) {
    if (is_image_preset_filename(filename)) {
        return load_image_preset(filename, img);
    }
    auto ext = get_extension(filename);
    if (ext == "exr" || ext == "EXR") {
        auto imgf = image<vec<float, N>>{};
        load_exr_image(filename, imgf);
        img.resize(imgf.size());
        linear_to_srgb8(img, imgf);
    } else if (ext == "pfm" || ext == "PFM") {
        auto imgf = image<vec<float, N>>{};
        load_pfm_image(filename, imgf);
        img.resize(imgf.size());
        linear_to_srgb8(img, imgf);
    } else if (ext == "hdr" || ext == "HDR") {
        auto imgf = image<vec<float, N>>{};
        load_stb_image(filename, imgf);
        img.resize(imgf.size());
        linear_to_srgb8(img, imgf);
    } else if (ext == "png" || ext == "PNG") {
        load_stb_image(filename, img);
    } else if (ext == "jpg" || ext == "JPG") {
        load_stb_image(filename, img);
    } else if (ext == "tga" || ext == "TGA") {
        load_stb_image(filename, img);
    } else if (ext == "bmp" || ext == "BMP") {
        load_stb_image(filename, img);
    } else if (ext == "ppm" || ext == "PPM") {
        load_pnm_image(filename, img);
    } else if (ext == "pgm" || ext == "PGM") {
        load_pnm_image(filename, img);
    } else if (ext == "json" || ext == "JSON") {
        load_json_image(filename, img);
    } else {
        throw imageio_error("unsupported image format " + ext);
    }
}

// Saves an ldr image.
template <int N>
void save_image(const string& filename, const image<vec<byte, N>>& img) {
    auto ext = get_extension(filename);
    if (ext == "png" || ext == "PNG") {
        save_png_image(filename, img);
    } else if (ext == "jpg" || ext == "JPG") {
        save_jpg_image(filename, img);
    } else if (ext == "tga" || ext == "TGA") {
        save_tga_image(filename, img);
    } else if (ext == "bmp" || ext == "BMP") {
        save_bmp_image(filename, img);
    } else if (ext == "ppm" || ext == "PPM") {
        save_pnm_image(filename, img);
    } else if (ext == "pgm" || ext == "PGM") {
        save_pnm_image(filename, img);
    } else if (ext == "hdr" || ext == "HDR") {
        auto imgf = image<vec<float, N>>{img.size()};
        srgb8_to_linear(imgf, img);
        save_hdr_image(filename, imgf);
    } else if (ext == "pfm" || ext == "PFM") {
        auto imgf = image<vec<float, N>>{img.size()};
        srgb8_to_linear(imgf, img);
        save_pfm_image(filename, imgf);
    } else if (ext == "exr" || ext == "EXR") {
        auto imgf = image<vec<float, N>>{img.size()};
        srgb8_to_linear(imgf, img);
        save_exr_image(filename, imgf);
    } else {
        throw imageio_error("unsupported image format " + ext);
    }
}

// Loads an hdr image.
template <int N>
void load_image_from_memory(
    const byte* data, int data_size, image<vec<byte, N>>& img) {
    load_stb_image_from_memory(data, data_size, img);
}

// Specializations
template void load_image<1>(const string& filename, image<vec<float, 1>>& img);
template void load_image<2>(const string& filename, image<vec<float, 2>>& img);
template void load_image<3>(const string& filename, image<vec<float, 3>>& img);
template void load_image<4>(const string& filename, image<vec<float, 4>>& img);
template void save_image<1>(
    const string& filename, const image<vec<float, 1>>& img);
template void save_image<2>(
    const string& filename, const image<vec<float, 2>>& img);
template void save_image<3>(
    const string& filename, const image<vec<float, 3>>& img);
template void save_image<4>(
    const string& filename, const image<vec<float, 4>>& img);
template void load_image_from_memory<1>(
    const byte* data, int data_size, image<vec<float, 1>>& img);
template void load_image_from_memory<2>(
    const byte* data, int data_size, image<vec<float, 2>>& img);
template void load_image_from_memory<3>(
    const byte* data, int data_size, image<vec<float, 3>>& img);
template void load_image_from_memory<4>(
    const byte* data, int data_size, image<vec<float, 4>>& img);
template void load_image<1>(const string& filename, image<vec<byte, 1>>& img);
template void load_image<2>(const string& filename, image<vec<byte, 2>>& img);
template void load_image<3>(const string& filename, image<vec<byte, 3>>& img);
template void load_image<4>(const string& filename, image<vec<byte, 4>>& img);
template void save_image<1>(
    const string& filename, const image<vec<byte, 1>>& img);
template void save_image<2>(
    const string& filename, const image<vec<byte, 2>>& img);
template void save_image<3>(
    const string& filename, const image<vec<byte, 3>>& img);
template void save_image<4>(
    const string& filename, const image<vec<byte, 4>>& img);
template void load_image_from_memory<1>(
    const byte* data, int data_size, image<vec<byte, 1>>& img);
template void load_image_from_memory<2>(
    const byte* data, int data_size, image<vec<byte, 2>>& img);
template void load_image_from_memory<3>(
    const byte* data, int data_size, image<vec<byte, 3>>& img);
template void load_image_from_memory<4>(
    const byte* data, int data_size, image<vec<byte, 4>>& img);

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR VOLUME IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Volume load
float* load_yvol(
    const char* filename, int* w, int* h, int* d, int* nc, int req) {
    auto fs = fopen(filename, "rb");
    if (!fs) return nullptr;
    auto fs_guard = unique_ptr<FILE, void (*)(FILE*)>{
        fs, [](FILE* f) { fclose(f); }};

    // buffer
    char buffer[4096];
    auto toks = vector<string>();

    // read magic
    if (!fgets(buffer, sizeof(buffer), fs)) return nullptr;
    toks = split_string(buffer);
    if (toks[0] != "YVOL") return nullptr;

    // read w, h
    if (!fgets(buffer, sizeof(buffer), fs)) return nullptr;
    toks = split_string(buffer);
    *w   = atoi(toks[0].c_str());
    *h   = atoi(toks[1].c_str());
    *d   = atoi(toks[2].c_str());
    *nc  = atoi(toks[3].c_str());

    // read data
    auto nvoxels = (size_t)(*w) * (size_t)(*h) * (size_t)(*d);
    auto nvalues = nvoxels * (size_t)(*nc);
    auto voxels  = unique_ptr<float[]>(new float[nvalues]);
    if (fread(voxels.get(), sizeof(float), nvalues, fs) != nvalues)
        return nullptr;

    // proper number of channels
    if (!req || *nc == req) return voxels.release();

    // pack into channels
    if (req < 0 || req > 4) {
        return nullptr;
    }
    auto cvoxels = unique_ptr<float[]>(new float[req * nvoxels]);
    for (auto i = 0; i < nvoxels; i++) {
        auto vp = voxels.get() + i * (*nc);
        auto cp = cvoxels.get() + i * req;
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
        } else if (*nc == 2) {
            switch (req) {
                case 1: cp[0] = vp[0]; break;
                case 2:
                    cp[0] = vp[0];
                    cp[1] = vp[1];
                    break;
                case 3:
                    cp[0] = vp[0];
                    cp[1] = vp[1];
                    break;
                case 4:
                    cp[0] = vp[0];
                    cp[1] = vp[1];
                    break;
            }
        } else if (*nc == 3) {
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
        } else if (*nc == 4) {
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
                    cp[3] = vp[3];
                    break;
            }
        }
    }
    return cvoxels.release();
}

// save pfm
bool save_yvol(
    const char* filename, int w, int h, int d, int nc, const float* voxels) {
    auto fs = fopen(filename, "wb");
    if (!fs) return false;
    auto fs_guard = unique_ptr<FILE, void (*)(FILE*)>{
        fs, [](FILE* f) { fclose(f); }};

    if (fprintf(fs, "YVOL\n") < 0) return false;
    if (fprintf(fs, "%d %d %d %d\n", w, h, d, nc) < 0) return false;
    auto nvalues = (size_t)w * (size_t)h * (size_t)d * (size_t)nc;
    if (fwrite(voxels, sizeof(float), nvalues, fs) != nvalues) return false;

    return true;
}

// Loads volume data from binary format.
void load_volume(const string& filename, volume1f& vol) {
    auto width = 0, height = 0, depth = 0, ncomp = 0;
    auto voxels = load_yvol(
        filename.c_str(), &width, &height, &depth, &ncomp, 1);
    if (!voxels) {
        throw imageio_error("error loading volume " + filename);
    }
    vol = volume{{width, height, depth}, (const float*)voxels};
    delete[] voxels;
}

// Saves volume data in binary format.
void save_volume(const string& filename, const volume1f& vol) {
    if (!save_yvol(filename.c_str(), vol.size().x, vol.size().y, vol.size().z,
            1, vol.data())) {
        throw imageio_error("error saving volume " + filename);
    }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// BUILTIN IMAGES
// -----------------------------------------------------------------------------
namespace yocto {

// Loads/saves a 1-4 channel builtin image.
template <int N>
void load_builtin_image(const string& name, image<vec<byte, N>>& img) {
    static const char* logo_render = R"(
        P2
        144 28
        255
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 212 87 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 255 62 0 0 0 0 0 14 27 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 187 245 13 0 0 0 80 255 90 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 29 88 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 88 251 10 0 0 0 40 200 253 255 234 106 0 0 6 255 157 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 12 147 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 90 69 0 0 0 0 0 0 0 0 0 0 0 144 125 0 0 0 0 0 117 37 0 0 0 0 0 0 0 79 255 101 0 0 0 178 232 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 115 255 33 0 0 0 0 0 0 0 0 0 0 0 0 0 0 145 205 0 0 0 47 239 210 74 57 144 232 24 0 6 255 157 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 18 251 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 146 123 0 0 0 0 0 0 0 0 0 0 0 43 35 0 61 87 0 0 208 61 0 0 0 0 0 0 0 3 224 199 0 0 24 251 129 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 115 255 33 0 0 0 0 0 0 0 0 0 0 0 0 0 0 201 149 0 0 0 180 238 20 0 0 0 16 0 0 6 255 157 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 18 251 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 146 123 0 0 0 0 0 0 0 0 0 0 0 0 0 0 119 150 0 0 208 61 0 0 0 0 0 0 0 0 118 255 41 0 117 250 26 0 7 147 223 232 167 19 0 0 0 5 143 224 234 162 20 166 227 255 210 201 3 0 7 147 223 232 167 19 0 0 0 0 8 249 93 0 0 45 255 146 0 0 0 0 0 0 0 6 255 157 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 14 192 114 207 0 84 224 220 61 0 68 143 166 220 78 0 0 142 234 160 251 0 0 134 234 195 25 0 123 105 211 89 12 177 236 158 4 0 44 217 214 189 123 0 0 0 149 76 0 189 108 0 160 55 123 107 83 236 240 179 0 208 128 220 174 6 0 0 0 0 0 19 247 140 0 214 168 0 0 176 248 122 111 237 212 2 0 0 167 251 128 127 223 59 103 185 255 143 117 0 0 176 248 122 111 237 212 2 0 0 0 58 255 36 0 0 97 255 77 0 0 0 0 0 0 0 6 255 157 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 16 248 166 36 13 234 44 73 215 0 80 248 65 80 193 0 66 222 32 97 251 0 65 207 21 135 152 0 144 230 81 12 129 157 17 189 88 0 194 126 17 208 123 0 0 0 131 125 7 228 164 0 224 23 144 125 5 127 156 10 0 208 179 13 205 65 0 0 0 0 0 0 158 233 61 255 59 0 46 255 121 0 0 91 255 83 0 40 254 134 0 0 0 0 0 115 255 33 0 0 46 255 121 0 0 91 255 83 0 0 0 114 235 0 0 0 130 255 51 0 2 17 17 17 7 0 6 255 157 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 16 255 43 0 76 191 0 1 242 20 80 192 0 36 235 0 138 142 0 18 251 0 140 127 0 52 212 0 144 171 0 0 204 63 0 116 148 11 255 14 0 146 123 0 0 0 84 164 46 155 201 9 231 0 144 125 0 119 150 0 0 208 64 0 164 107 0 0 0 0 0 0 49 255 220 206 0 0 114 255 46 0 0 17 255 148 0 111 255 54 0 0 0 0 0 115 255 33 0 0 114 255 46 0 0 17 255 148 0 0 0 171 179 0 0 0 159 255 29 0 20 255 255 255 109 0 6 255 157 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 16 253 0 0 101 242 214 214 249 40 80 189 0 34 236 0 162 121 0 18 251 0 165 232 214 219 232 0 144 126 0 0 229 222 214 229 168 34 249 0 0 146 123 0 0 0 38 203 88 107 206 48 187 0 144 125 0 119 150 0 0 208 61 0 162 108 0 0 0 0 0 0 0 197 255 98 0 0 144 255 22 0 0 0 249 177 0 141 255 28 0 0 0 0 0 115 255 33 0 0 144 255 22 0 0 0 249 177 0 0 0 227 123 0 0 0 143 255 40 0 0 79 108 255 109 0 6 255 157 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 16 253 0 0 85 188 4 4 4 0 80 189 0 34 236 0 149 132 0 18 251 0 149 125 4 4 3 0 144 125 0 0 213 62 4 4 2 21 255 4 0 146 123 0 0 0 2 230 130 67 178 116 141 0 144 125 0 119 150 0 0 208 61 0 162 108 0 0 0 0 0 0 0 130 255 33 0 0 151 255 17 0 0 0 245 183 0 152 255 21 0 0 0 0 0 115 255 33 0 0 151 255 17 0 0 0 245 183 0 0 28 255 67 0 0 0 116 255 59 0 0 0 39 255 109 0 6 255 157 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 16 253 0 0 32 236 19 2 33 0 80 189 0 34 236 0 98 195 4 64 251 0 93 191 4 11 24 0 144 125 0 0 157 131 0 27 8 1 225 72 1 190 123 0 0 0 0 201 196 26 137 197 96 0 144 125 0 109 159 0 0 208 61 0 162 108 0 0 0 0 0 0 0 130 255 33 0 0 122 255 36 0 0 8 255 152 0 124 255 42 0 0 0 0 0 115 255 33 0 0 122 255 36 0 0 8 255 152 0 0 84 253 13 0 0 0 79 255 108 0 0 0 39 255 109 0 6 255 157 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 16 253 0 0 0 123 238 224 155 0 80 189 0 34 236 0 11 201 218 156 248 0 6 178 225 234 97 0 144 125 0 0 31 216 214 228 49 0 89 244 198 179 123 0 0 0 0 154 239 0 96 253 50 0 144 125 0 34 224 201 25 208 61 0 162 108 0 0 0 0 0 0 0 130 255 33 0 0 70 255 96 0 0 67 255 97 0 76 255 104 0 0 0 0 0 113 255 35 0 0 70 255 96 0 0 67 255 97 0 0 141 210 0 0 0 0 5 230 200 2 0 0 39 255 109 0 6 255 157 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 0 0 0 0 0 0 11 14 0 0 0 0 0 0 0 0 0 0 18 0 0 0 0 0 19 6 0 0 0 0 0 0 0 0 23 1 0 0 0 12 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 9 0 0 0 0 0 0 0 0 0 0 0 0 0 130 255 33 0 0 1 211 223 54 44 205 229 7 0 2 216 230 69 74 169 31 0 55 255 120 59 26 1 211 223 54 44 205 229 7 0 0 197 154 0 0 0 0 0 109 255 166 59 78 165 255 109 0 6 255 201 113 113 113 105 0 0 0 0 0 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 130 255 33 0 0 0 30 204 255 255 214 40 0 0 0 35 208 255 255 212 47 0 1 174 252 243 102 0 30 204 255 255 214 40 0 0 6 247 97 0 0 0 0 0 0 103 235 255 255 242 146 24 0 6 255 255 255 255 255 213 0 0 0 0 0 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 24 29 0 0 0 0 0 0 0 25 31 0 0 0 0 0 23 9 0 0 0 0 24 29 0 0 0 0 54 255 41 0 0 0 0 0 0 0 0 33 30 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 63 188 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
        )";
    if (name == "logo-render") {
        load_pnm_image_from_string(logo_render, img);
    } else {
        throw imageio_error("unknown builtin image " + name);
    }
}

// Loads/saves a 1-4 channel builtin image.
template <int N>
void load_builtin_image(const string& name, image<vec<float, N>>& img) {
    auto img8 = image<vec<byte, N>>();
    load_builtin_image(name, img8);
    img.resize(img8.size());
    srgb8_to_linear(img, img8);
}

// Specializations
template void load_builtin_image<1>(
    const string& filename, image<vec<float, 1>>& img);
template void load_builtin_image<2>(
    const string& filename, image<vec<float, 2>>& img);
template void load_builtin_image<3>(
    const string& filename, image<vec<float, 3>>& img);
template void load_builtin_image<4>(
    const string& filename, image<vec<float, 4>>& img);
template void load_builtin_image<1>(
    const string& filename, image<vec<byte, 1>>& img);
template void load_builtin_image<2>(
    const string& filename, image<vec<byte, 2>>& img);
template void load_builtin_image<3>(
    const string& filename, image<vec<byte, 3>>& img);
template void load_builtin_image<4>(
    const string& filename, image<vec<byte, 4>>& img);

}  // namespace yocto
