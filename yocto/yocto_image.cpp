//
// Implementation for Yocto/Image.
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

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_image.h"

#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wunused-variable"
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#endif

// #ifndef __clang_analyzer__

#define STB_IMAGE_IMPLEMENTATION
#include "ext/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "ext/stb_image_write.h"

#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "ext/stb_image_resize.h"

#define TINYEXR_IMPLEMENTATION
#include "ext/tinyexr.h"

// #endif

#ifndef _WIN32
#pragma GCC diagnostic pop
#endif

#include "ext/ArHosekSkyModel.cpp"

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGEIO
// -----------------------------------------------------------------------------
namespace yocto {

// Split a string
static inline vector<string> _image_split_string(const string& str) {
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
static inline float* load_pfm(
    const char* filename, int* w, int* h, int* nc, int req) {
    auto fs = fopen(filename, "rb");
    if (!fs) return nullptr;
    auto fs_guard = unique_ptr<FILE, void (*)(FILE*)>{
        fs, [](FILE* f) { fclose(f); }};

    // buffer
    char buffer[4096];
    auto toks = vector<string>();

    // read magic
    if (!fgets(buffer, sizeof(buffer), fs)) return nullptr;
    toks = _image_split_string(buffer);
    if (toks[0] == "Pf")
        *nc = 1;
    else if (toks[0] == "PF")
        *nc = 3;
    else
        return nullptr;

    // read w, h
    if (!fgets(buffer, sizeof(buffer), fs)) return nullptr;
    toks = _image_split_string(buffer);
    *w   = atoi(toks[0].c_str());
    *h   = atoi(toks[1].c_str());

    // read scale
    if (!fgets(buffer, sizeof(buffer), fs)) return nullptr;
    toks   = _image_split_string(buffer);
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
static inline bool save_pfm(
    const char* filename, int w, int h, int nc, const float* pixels) {
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

// load pfm image
template <int N>
static inline void load_pfm_image(
    const string& filename, image<vec<float, N>>& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = load_pfm(filename.c_str(), &width, &height, &ncomp, N);
    if (!pixels) {
        throw io_error("error loading image " + filename);
    }
    img = image{{width, height}, (const vec<float, N>*)pixels};
    delete[] pixels;
}
template <int N>
static inline void save_pfm_image(
    const string& filename, const image<vec<float, N>>& img) {
    if (!save_pfm(filename.c_str(), img.size().x, img.size().y, N,
            (float*)img.data())) {
        throw io_error("error saving image " + filename);
    }
}

// load exr image weith tiny exr
static inline const char* get_tinyexr_error(int error) {
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
        default: throw io_error("unknown tinyexr error");
    }
}

template <int N>
static inline void load_exr_image(
    const string& filename, image<vec<float, N>>& img) {
    // TODO
    if (N != 4) throw runtime_error("bad number of channels");
    auto width = 0, height = 0;
    auto pixels = (float*)nullptr;
    if (auto error = LoadEXR(
            &pixels, &width, &height, filename.c_str(), nullptr);
        error < 0) {
        throw io_error("error loading image " + filename + "("s +
                       get_tinyexr_error(error) + ")"s);
    }
    if (!pixels) {
        throw io_error("error loading image " + filename);
    }
    img = image{{width, height}, (const vec<float, N>*)pixels};
    free(pixels);
}
template <int N>
static inline void save_exr_image(
    const string& filename, const image<vec<float, N>>& img) {
    // TODO
    if (N != 4) throw runtime_error("bad number of channels");
    if (!SaveEXR((float*)img.data(), img.size().x, img.size().y, N,
            filename.c_str())) {
        throw io_error("error saving image " + filename);
    }
}

// load an image using stbi library
template <int N>
static inline void load_stb_image(
    const string& filename, image<vec<byte, N>>& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_load(filename.c_str(), &width, &height, &ncomp, N);
    if (!pixels) {
        throw io_error("error loading image " + filename);
    }
    img = image{{width, height}, (const vec<byte, N>*)pixels};
    free(pixels);
}
template <int N>
static inline void load_stb_image(
    const string& filename, image<vec<float, N>>& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_loadf(filename.c_str(), &width, &height, &ncomp, N);
    if (!pixels) {
        throw io_error("error loading image " + filename);
    }
    img = image{{width, height}, (const vec<float, N>*)pixels};
    free(pixels);
}

// save an image with stbi
template <int N>
static inline void save_png_image(
    const string& filename, const image<vec<byte, N>>& img) {
    if (!stbi_write_png(filename.c_str(), img.size().x, img.size().y, N,
            img.data(), img.size().x * 4)) {
        throw io_error("error saving image " + filename);
    }
}
template <int N>
static inline void save_jpg_image(
    const string& filename, const image<vec<byte, N>>& img) {
    if (!stbi_write_jpg(
            filename.c_str(), img.size().x, img.size().y, 4, img.data(), 75)) {
        throw io_error("error saving image " + filename);
    }
}
template <int N>
static inline void save_tga_image(
    const string& filename, const image<vec<byte, N>>& img) {
    if (!stbi_write_tga(
            filename.c_str(), img.size().x, img.size().y, 4, img.data())) {
        throw io_error("error saving image " + filename);
    }
}
template <int N>
static inline void save_bmp_image(
    const string& filename, const image<vec<byte, N>>& img) {
    if (!stbi_write_bmp(
            filename.c_str(), img.size().x, img.size().y, 4, img.data())) {
        throw io_error("error saving image " + filename);
    }
}
template <int N>
static inline void save_hdr_image(
    const string& filename, const image<vec<float, N>>& img) {
    if (!stbi_write_hdr(filename.c_str(), img.size().x, img.size().y, 4,
            (float*)img.data())) {
        throw io_error("error saving image " + filename);
    }
}

// load an image using stbi library
template <int N>
static inline void load_stb_image_from_memory(
    const byte* data, int data_size, image<vec<byte, N>>& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_load_from_memory(
        data, data_size, &width, &height, &ncomp, 4);
    if (!pixels) {
        throw io_error("error loading in-memory image");
    }
    img = image{{width, height}, (const vec<byte, N>*)pixels};
    free(pixels);
}
template <int N>
static inline void load_stb_image_from_memory(
    const byte* data, int data_size, image<vec<float, N>>& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_loadf_from_memory(
        data, data_size, &width, &height, &ncomp, 4);
    if (!pixels) {
        throw io_error("error loading in-memory image {}");
    }
    img = image{{width, height}, (const vec<float, N>*)pixels};
    free(pixels);
}

#if 0
static inline void apply_json_procedural(const json& js, image<vec4f>& img) {
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

#endif

template <int N>
inline void load_image_preset(
    const string& filename, image<vec<float, N>>& img) {
    auto [type, nfilename] = get_image_preset_type(filename);
    if constexpr (N == 4) {
        img.resize({1024, 1024});
        if (type == "images2") img.resize({2048, 1024});
        make_image_preset(img, type);
    } else {
        auto img4 = image<vec<float, 4>>({1024, 1024});
        if (type == "images2") img4.resize({2048, 1024});
        make_image_preset(img4, type);
        img.resize(img4.size());
        convert_color_channels(img, img4);
    }
}
template <int N>
inline void load_image_preset(
    const string& filename, image<vec<byte, N>>& img) {
    auto imgf = image<vec<float, N>>{};
    load_image_preset(filename, imgf);
    img.resize(imgf.size());
    linear_to_srgb(img, imgf);
}

// Loads an hdr image.
template <int N>
inline void load_channel_image(
    const string& filename, image<vec<float, N>>& img) {
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
    } else if (!is_hdr_filename(filename)) {
        auto img8 = image<vec<byte, N>>{};
        load_channel_image(filename, img8);
        srgb_to_linear(img, img8);
    } else {
        throw io_error("unsupported image format " + ext);
    }
}

// Saves an hdr image.
template <int N>
inline void save_channel_image(
    const string& filename, const image<vec<float, N>>& img) {
    auto ext = get_extension(filename);
    if (ext == "hdr" || ext == "HDR") {
        save_hdr_image(filename, img);
    } else if (ext == "pfm" || ext == "PFM") {
        save_pfm_image(filename, img);
    } else if (ext == "exr" || ext == "EXR") {
        save_exr_image(filename, img);
    } else if (!is_hdr_filename(filename)) {
        auto img8 = image<vec<byte, N>>{img.size()};
        linear_to_srgb(img8, img);
        save_channel_image(filename, img8);
    } else {
        throw io_error("unsupported image format " + ext);
    }
}

// Loads an hdr image.
template <int N>
inline void load_channel_image(
    const string& filename, image<vec<byte, N>>& img) {
    if (is_image_preset_filename(filename)) {
        return load_image_preset(filename, img);
    }
    auto ext = get_extension(filename);
    if (ext == "png" || ext == "PNG") {
        load_stb_image(filename, img);
    } else if (ext == "jpg" || ext == "JPG") {
        load_stb_image(filename, img);
    } else if (ext == "tga" || ext == "TGA") {
        load_stb_image(filename, img);
    } else if (ext == "bmp" || ext == "BMP") {
        load_stb_image(filename, img);
    } else if (is_hdr_filename(filename)) {
        auto imgf = image<vec<float, N>>{};
        load_channel_image(filename, imgf);
        linear_to_srgb(img, imgf);
    } else {
        throw io_error("unsupported image format " + ext);
    }
}

// Saves an ldr image.
template <int N>
inline void save_channel_image(
    const string& filename, const image<vec<byte, N>>& img) {
    auto ext = get_extension(filename);
    if (ext == "png" || ext == "PNG") {
        save_png_image(filename, img);
    } else if (ext == "jpg" || ext == "JPG") {
        save_jpg_image(filename, img);
    } else if (ext == "tga" || ext == "TGA") {
        save_tga_image(filename, img);
    } else if (ext == "bmp" || ext == "BMP") {
        save_bmp_image(filename, img);
    } else if (is_hdr_filename(filename)) {
        auto imgf = image<vec<float, N>>{img.size()};
        srgb_to_linear(imgf, img);
        save_channel_image(filename, imgf);
    } else {
        throw io_error("unsupported image format " + ext);
    }
}

// Loads/saves a 1-4 channels float image in linear color space.
void load_image(const string& filename, image<float>& img) {
    load_channel_image(filename, (image<vec1f>&)img);
}
void load_image(const string& filename, image<vec2f>& img) {
    load_channel_image(filename, img);
}
void load_image(const string& filename, image<vec3f>& img) {
    load_channel_image(filename, img);
}
void load_image(const string& filename, image<vec4f>& img) {
    load_channel_image(filename, img);
}
void save_image(const string& filename, const image<float>& img) {
    save_channel_image(filename, (const image<vec1f>&)img);
}
void save_image(const string& filename, const image<vec2f>& img) {
    save_channel_image(filename, img);
}
void save_image(const string& filename, const image<vec3f>& img) {
    save_channel_image(filename, img);
}
void save_image(const string& filename, const image<vec4f>& img) {
    save_channel_image(filename, img);
}

// Loads/saves a 1-4 byte image in sRGB color space.
void load_image(const string& filename, image<byte>& img) {
    load_channel_image(filename, (image<vec1b>&)img);
}
void load_image(const string& filename, image<vec2b>& img) {
    load_channel_image(filename, img);
}
void load_image(const string& filename, image<vec3b>& img) {
    load_channel_image(filename, img);
}
void load_image(const string& filename, image<vec4b>& img) {
    load_channel_image(filename, img);
}
void save_image(const string& filename, const image<byte>& img) {
    save_channel_image(filename, (const image<vec1b>&)img);
}
void save_image(const string& filename, const image<vec2b>& img) {
    save_channel_image(filename, img);
}
void save_image(const string& filename, const image<vec3b>& img) {
    save_channel_image(filename, img);
}
void save_image(const string& filename, const image<vec4b>& img) {
    save_channel_image(filename, img);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR VOLUME IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Volume load
static inline float* load_yvol(
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
    toks = _image_split_string(buffer);
    if (toks[0] != "YVOL") return nullptr;

    // read w, h
    if (!fgets(buffer, sizeof(buffer), fs)) return nullptr;
    toks = _image_split_string(buffer);
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
static inline bool save_yvol(
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
void load_volume(const string& filename, volume<float>& vol) {
    auto width = 0, height = 0, depth = 0, ncomp = 0;
    auto voxels = load_yvol(
        filename.c_str(), &width, &height, &depth, &ncomp, 1);
    if (!voxels) {
        throw io_error("error loading volume " + filename);
    }
    vol = volume{{width, height, depth}, (const float*)voxels};
    delete[] voxels;
}

// Saves volume data in binary format.
void save_volume(const string& filename, const volume<float>& vol) {
    if (!save_yvol(filename.c_str(), vol.size().x, vol.size().y, vol.size().z,
            1, vol.data())) {
        throw io_error("error saving volume " + filename);
    }
}

}  // namespace yocto
