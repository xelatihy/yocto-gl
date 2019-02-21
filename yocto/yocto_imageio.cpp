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
#include "yocto_json.h"
#include "yocto_utils.h"

#include <climits>
#include <cstdlib>

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
vector<float> load_pfm(const char* filename, int& w, int& h, int& nc, int req) {
    auto fs = input_file(filename, true);

    // buffer
    auto buffer = ""s;
    auto toks   = vector<string>();

    // read magic
    if (!read_line(fs, buffer)) return {};
    toks = split_string(buffer);
    if (toks[0] == "Pf")
        nc = 1;
    else if (toks[0] == "PF")
        nc = 3;
    else
        return {};

    // read w, h
    if (!read_line(fs, buffer)) return {};
    toks = split_string(buffer);
    w    = atoi(toks[0].c_str());
    h    = atoi(toks[1].c_str());

    // read scale
    if (!read_line(fs, buffer)) return {};
    toks   = split_string(buffer);
    auto s = atof(toks[0].c_str());

    // read the data (flip y)
    auto npixels = w * h;
    auto nvalues = w * h * nc;
    auto nrow    = w * nc;
    auto pixels  = vector<float>(nvalues);
    for (auto j = h - 1; j >= 0; j--) {
        read_values(fs, pixels.data() + j * nrow, nrow);
    }

    // endian conversion
    if (s > 0) {
        for (auto i = 0; i < nvalues; ++i) {
            auto dta = (uint8_t*)(pixels.data() + i);
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
    if (!req || nc == req) return pixels;

    // pack into channels
    if (req < 0 || req > 4) {
        return {};
    }
    auto cpixels = vector<float>(req * npixels);
    for (auto i = 0; i < npixels; i++) {
        auto vp = pixels.data() + i * nc;
        auto cp = cpixels.data() + i * req;
        if (nc == 1) {
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
    return cpixels;
}

// save pfm
bool save_pfm(const char* filename, int w, int h, int nc, const float* pixels) {
    auto fs = fopen(filename, "wb");
    if (!fs) return false;

    fprintf(fs, "%s\n", (nc == 1) ? "Pf" : "PF");
    fprintf(fs, "%d %d\n", w, h);
    fprintf(fs, "-1\n");
    if (nc == 1 || nc == 3) {
        fwrite(pixels, sizeof(float), w * h * nc, fs);
    } else {
        for (auto i = 0; i < w * h; i++) {
            auto vz = 0.0f;
            auto v  = pixels + i * nc;
            fwrite(v + 0, sizeof(float), 1, fs);
            fwrite(v + 1, sizeof(float), 1, fs);
            if (nc == 2)
                fwrite(&vz, sizeof(float), 1, fs);
            else
                fwrite(v + 2, sizeof(float), 1, fs);
        }
    }

    fclose(fs);

    return true;
}

// load pfm image
void load_pfm_image(const string& filename, image4f& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = load_pfm(filename.c_str(), width, height, ncomp, 4);
    if (pixels.empty()) {
        throw io_error("error loading image " + filename);
    }
    img = image{{width, height}, (const vec4f*)pixels.data()};
}
void save_pfm_image(const string& filename, const image4f& img) {
    if (!save_pfm(filename.c_str(), img.size().x, img.size().y, 4,
            (float*)img.data())) {
        throw io_error("error saving image " + filename);
    }
}

// load exr image weith tiny exr
void load_exr_image(const string& filename, image4f& img) {
    auto width = 0, height = 0;
    auto pixels = (float*)nullptr;
    if (LoadEXR(&pixels, &width, &height, filename.c_str(), nullptr) < 0) {
        throw io_error("error loading image " + filename);
    }
    if (!pixels) {
        throw io_error("error loading image " + filename);
    }
    img = image{{width, height}, (const vec4f*)pixels};
    free(pixels);
}
void save_exr_image(const string& filename, const image4f& img) {
    if (!SaveEXR((float*)img.data(), img.size().x, img.size().y, 4,
            filename.c_str())) {
        throw io_error("error saving image " + filename);
    }
}

// load an image using stbi library
void load_stb_image(const string& filename, image4b& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_load(filename.c_str(), &width, &height, &ncomp, 4);
    if (!pixels) {
        throw io_error("error loading image " + filename);
    }
    img = image{{width, height}, (const vec4b*)pixels};
    free(pixels);
}
void load_stb_image(const string& filename, image4f& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_loadf(filename.c_str(), &width, &height, &ncomp, 4);
    if (!pixels) {
        throw io_error("error loading image " + filename);
    }
    img = image{{width, height}, (const vec4f*)pixels};
    free(pixels);
}

// save an image with stbi
void save_png_image(const string& filename, const image4b& img) {
    if (!stbi_write_png(filename.c_str(), img.size().x, img.size().y, 4,
            img.data(), img.size().x * 4)) {
        throw io_error("error saving image " + filename);
    }
}
void save_jpg_image(const string& filename, const image4b& img) {
    if (!stbi_write_jpg(
            filename.c_str(), img.size().x, img.size().y, 4, img.data(), 75)) {
        throw io_error("error saving image " + filename);
    }
}
void save_tga_image(const string& filename, const image4b& img) {
    if (!stbi_write_tga(
            filename.c_str(), img.size().x, img.size().y, 4, img.data())) {
        throw io_error("error saving image " + filename);
    }
}
void save_bmp_image(const string& filename, const image4b& img) {
    if (!stbi_write_bmp(
            filename.c_str(), img.size().x, img.size().y, 4, img.data())) {
        throw io_error("error saving image " + filename);
    }
}
void save_hdr_image(const string& filename, const image4f& img) {
    if (!stbi_write_hdr(filename.c_str(), img.size().x, img.size().y, 4,
            (float*)img.data())) {
        throw io_error("error saving image " + filename);
    }
}

// load an image using stbi library
void load_stb_image_from_memory(const byte* data, int data_size, image4b& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_load_from_memory(
        data, data_size, &width, &height, &ncomp, 4);
    if (!pixels) {
        throw io_error("error loading in-memory image");
    }
    img = image{{width, height}, (const vec4b*)pixels};
    free(pixels);
}
void load_stbi_image_from_memory(
    const byte* data, int data_size, image4f& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_loadf_from_memory(
        data, data_size, &width, &height, &ncomp, 4);
    if (!pixels) {
        throw io_error("error loading in-memory image {}");
    }
    img = image{{width, height}, (const vec4f*)pixels};
    free(pixels);
}

void apply_json_procedural(const json& js, image4f& img) {
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
        make_bumpdimple_image(img, js.value("tile", 8));
    } else if (type == "uvramp") {
        make_uvramp_image(img);
    } else if (type == "gammaramp") {
        make_gammaramp_image(img);
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
        make_noise_image(img, js.value("scale", 1.0f), js.value("wrap", true));
    } else if (type == "fbm") {
        make_fbm_image(img, js.value("scale", 1.0f),
            js.value("lacunarity", 2.0f), js.value("gain", 0.5f),
            js.value("octaves", 6), js.value("wrap", true));
    } else if (type == "ridge") {
        make_ridge_image(img, js.value("scale", 1.0f),
            js.value("lacunarity", 2.0f), js.value("gain", 0.5f),
            js.value("offset", 1.0f), js.value("octaves", 6),
            js.value("wrap", true));
    } else if (type == "turbulence") {
        make_turbulence_image(img, js.value("scale", 1.0f),
            js.value("lacunarity", 2.0f), js.value("gain", 0.5f),
            js.value("octaves", 6), js.value("wrap", true));
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

void apply_json_procedural(const json& js, image4b& img) {
    auto imgf = image4f{};
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
void load_json_image(const string& filename, image4f& img) {
    auto js = json();
    load_json(filename, js);
    apply_json_procedural(js, img);
}
void load_json_image(const string& filename, image4b& img) {
    auto js = json();
    load_json(filename, js);
    apply_json_procedural(js, img);
}

// check hdr extensions
bool is_hdr_filename(const string& filename) {
    auto ext = get_extension(filename);
    return ext == "hdr" || ext == "exr" || ext == "pfm";
}

// Loads an hdr image.
void load_image(const string& filename, image4f& img) {
    auto ext = get_extension(filename);
    if (ext == "exr" || ext == "EXR") {
        load_exr_image(filename, img);
    } else if (ext == "pfm" || ext == "PFM") {
        load_pfm_image(filename, img);
    } else if (ext == "hdr" || ext == "HDR") {
        load_stb_image(filename, img);
    } else if (ext == "png" || ext == "PNG") {
        auto img8 = image4b{};
        load_stb_image(filename, img8);
        img.resize(img8.size());
        srgb_to_linear(img, img8);
    } else if (ext == "jpg" || ext == "JPG") {
        auto img8 = image4b{};
        load_stb_image(filename, img8);
        img.resize(img8.size());
        srgb_to_linear(img, img8);
    } else if (ext == "tga" || ext == "TGA") {
        auto img8 = image4b{};
        load_stb_image(filename, img8);
        img.resize(img8.size());
        srgb_to_linear(img, img8);
    } else if (ext == "bmp" || ext == "BMP") {
        auto img8 = image4b{};
        load_stb_image(filename, img8);
        img.resize(img8.size());
        srgb_to_linear(img, img8);
    } else if (ext == "json" || ext == "JSON") {
        load_json_image(filename, img);
    } else {
        throw io_error("unsupported image format " + ext);
    }
}

// Saves an hdr image.
void save_image(const string& filename, const image4f& img) {
    auto ext = get_extension(filename);
    if (ext == "png" || ext == "PNG") {
        auto img8 = image4b{img.size()};
        linear_to_srgb(img8, img);
        save_png_image(filename, img8);
    } else if (ext == "jpg" || ext == "JPG") {
        auto img8 = image4b{img.size()};
        linear_to_srgb(img8, img);
        save_jpg_image(filename, img8);
    } else if (ext == "tga" || ext == "TGA") {
        auto img8 = image4b{img.size()};
        linear_to_srgb(img8, img);
        save_tga_image(filename, img8);
    } else if (ext == "bmp" || ext == "BMP") {
        auto img8 = image4b{img.size()};
        linear_to_srgb(img8, img);
        save_bmp_image(filename, img8);
    } else if (ext == "hdr" || ext == "HDR") {
        auto img8 = image4b{img.size()};
        linear_to_srgb(img8, img);
        save_hdr_image(filename, img);
    } else if (ext == "pfm" || ext == "PFM") {
        save_pfm_image(filename, img);
    } else if (ext == "exr" || ext == "EXR") {
        save_exr_image(filename, img);
    } else {
        throw io_error("unsupported image format " + ext);
    }
}

// Loads an hdr image.
void load_image_from_memory(const byte* data, int data_size, image4f& img) {
    load_stbi_image_from_memory(data, data_size, img);
}

// Loads an hdr image.
void load_image(const string& filename, image4b& img) {
    auto ext = get_extension(filename);
    if (ext == "exr" || ext == "EXR") {
        auto imgf = image4f{};
        load_exr_image(filename, imgf);
        img.resize(imgf.size());
        linear_to_srgb(img, imgf);
    } else if (ext == "pfm" || ext == "PFM") {
        auto imgf = image4f{};
        load_pfm_image(filename, imgf);
        img.resize(imgf.size());
        linear_to_srgb(img, imgf);
    } else if (ext == "hdr" || ext == "HDR") {
        auto imgf = image4f{};
        load_stb_image(filename, imgf);
        img.resize(imgf.size());
        linear_to_srgb(img, imgf);
    } else if (ext == "png" || ext == "PNG") {
        load_stb_image(filename, img);
    } else if (ext == "jpg" || ext == "JPG") {
        load_stb_image(filename, img);
    } else if (ext == "tga" || ext == "TGA") {
        load_stb_image(filename, img);
    } else if (ext == "bmp" || ext == "BMP") {
        load_stb_image(filename, img);
    } else if (ext == "json" || ext == "JSON") {
        load_json_image(filename, img);
    } else {
        throw io_error("unsupported image format " + ext);
    }
}

// Saves an ldr image.
void save_image(const string& filename, const image4b& img) {
    auto ext = get_extension(filename);
    if (ext == "png" || ext == "PNG") {
        save_png_image(filename, img);
    } else if (ext == "jpg" || ext == "JPG") {
        save_jpg_image(filename, img);
    } else if (ext == "tga" || ext == "TGA") {
        save_tga_image(filename, img);
    } else if (ext == "bmp" || ext == "BMP") {
        save_bmp_image(filename, img);
    } else if (ext == "hdr" || ext == "HDR") {
        auto imgf = image4f{img.size()};
        srgb_to_linear(imgf, img);
        save_hdr_image(filename, imgf);
    } else if (ext == "pfm" || ext == "PFM") {
        auto imgf = image4f{img.size()};
        srgb_to_linear(imgf, img);
        save_pfm_image(filename, imgf);
    } else if (ext == "exr" || ext == "EXR") {
        auto imgf = image4f{img.size()};
        srgb_to_linear(imgf, img);
        save_exr_image(filename, imgf);
    } else {
        throw io_error("unsupported image format " + ext);
    }
}

// Loads an ldr image.
void load_image_from_memory(const byte* data, int data_size, image4b& img) {
    load_stb_image_from_memory(data, data_size, img);
}

// Convenience helper that saves an HDR images as wither a linear HDR file or
// a tonemapped LDR file depending on file name
void save_tonemapped_image(const string& filename, const image4f& hdr,
    float exposure, bool filmic, bool srgb) {
    if (is_hdr_filename(filename)) {
        save_image(filename, hdr);
    } else {
        auto ldr = image4b{hdr.size()};
        tonemap_image(ldr, hdr, exposure, filmic, srgb);
        save_image(filename, ldr);
    }
}

// Resize image.
void resize_image(image4f& res_img, const image4f& img) {
    stbir_resize_float_generic((float*)img.data(), img.size().x, img.size().y,
        sizeof(vec4f) * img.size().x, (float*)res_img.data(), res_img.size().x,
        res_img.size().y, sizeof(vec4f) * res_img.size().x, 4, 3, 0,
        STBIR_EDGE_CLAMP, STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR,
        nullptr);
}
void resize_image(image4f& res_img, const image4f& img, const vec2i& size_) {
    auto size = size_;
    if (size == zero2i) {
        throw std::invalid_argument("bad image size in resize_image");
    }
    if (size.y == 0) {
        size.y = (int)round(size.x * (float)img.size().y / (float)img.size().x);
    } else if (size.x == 0) {
        size.x = (int)round(size.y * (float)img.size().x / (float)img.size().y);
    }
    res_img = {size};
    resize_image(res_img, img);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR VOLUME IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Loads volume data from binary format.
void load_volume(const string& filename, volume1f& vol) {
    auto fs   = input_file(filename, true);
    auto size = zero3i;
    read_value(fs, size);
    vol.resize(size);
    read_values(fs, vol._voxels);
}

// Saves volume data in binary format.
void save_volume(const string& filename, const volume1f& vol) {
    auto fs = output_file(filename, true);
    write_value(fs, vol.size());
    write_values(fs, vol._voxels);
}

}  // namespace yocto
