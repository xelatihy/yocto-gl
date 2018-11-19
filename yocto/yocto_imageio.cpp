//
// Implementation for Yocto/ImageIO.
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

#include "yocto_imageio.h"
#include "yocto_utils.h"
#include "yocto_json.h"

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
// IMPLEMENTATION FOR IMAGEIO
// -----------------------------------------------------------------------------
namespace yocto {

// Split a string
vector<string> split_string(const string& str) {
    auto ret = vector<string>();
    if (empty(str)) return ret;
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
    auto fs = open(filename, "rb");
    if (!fs) return {};

    // buffer
    char buffer[256];
    auto toks = vector<string>();

    // read magic
    if (!fgets(buffer, 256, fs.fs)) return {};
    toks = split_string(buffer);
    if (toks[0] == "Pf")
        nc = 1;
    else if (toks[0] == "PF")
        nc = 3;
    else
        return {};

    // read w, h
    if (!fgets(buffer, 256, fs.fs)) return {};
    toks = split_string(buffer);
    w    = atoi(toks[0].c_str());
    h    = atoi(toks[1].c_str());

    // read scale
    if (!fgets(buffer, 256, fs.fs)) return {};
    toks   = split_string(buffer);
    auto s = atof(toks[0].c_str());

    // read the data (flip y)
    auto npixels = w * h;
    auto nvalues = w * h * nc;
    auto nrow    = w * nc;
    auto pixels  = vector<float>(nvalues);
    for (auto j = h - 1; j >= 0; j--) {
        if (fread(data(pixels) + j * nrow, sizeof(float), nrow, fs.fs) != nrow) {
            return {};
        }
    }

    // endian conversion
    if (s > 0) {
        for (auto i = 0; i < nvalues; ++i) {
            auto dta = (uint8_t*)(data(pixels) + i);
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
        auto vp = data(pixels) + i * nc;
        auto cp = data(cpixels) + i * req;
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
bool load_pfm_image(const string& filename, image4f& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = load_pfm(filename.c_str(), width, height, ncomp, 4);
    if (empty(pixels)) {
        log_io_error("error loading image {}", filename);
        return false;
    }
    img = make_image(width, height, (vec4f*)data(pixels));
    return true;
}
bool save_pfm_image(const string& filename, const image4f& img) {
    if (!save_pfm(filename.c_str(), img.width, img.height, 4, (float*)data(img))) {
        log_io_error("error saving image {}", filename);
        return false;
    }
    return true;
}

// load exr image weith tiny exr
bool load_exr_image(const string& filename, image4f& img) {
    auto width = 0, height = 0;
    auto pixels = (vec4f*)nullptr;
    if (LoadEXR((float**)&pixels, &width, &height, filename.c_str(), nullptr) <
        0) {
        log_io_error("error loading image {}", filename);
        return false;
    }
    if (!pixels) {
        log_io_error("error loading image {}", filename);
        return false;
    }
    img = make_image(width, height, pixels);
    free(pixels);
    return true;
}
bool save_exr_image(const string& filename, const image4f& img) {
    if (!SaveEXR((float*)data(img), img.width, img.height, 4, filename.c_str())) {
        log_io_error("error saving image {}", filename);
        return false;
    }
    return true;
}

// load an image using stbi library
bool load_stb_image(const string& filename, image4b& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = (vec4b*)stbi_load(filename.c_str(), &width, &height, &ncomp, 4);
    if (!pixels) {
        log_io_error("error loading image {}", filename);
        return false;
    }
    img = make_image(width, height, pixels);
    free(pixels);
    return true;
}
bool load_stb_image(const string& filename, image4f& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = (vec4f*)stbi_loadf(
        filename.c_str(), &width, &height, &ncomp, 4);
    if (!pixels) {
        log_io_error("error loading image {}", filename);
        return false;
    }
    img = make_image(width, height, pixels);
    free(pixels);
    return true;
}

// save an image with stbi
bool save_png_image(const string& filename, const image4b& img) {
    if (!stbi_write_png(filename.c_str(), img.width, img.height, 4, data(img),
            img.width * 4)) {
        log_io_error("error saving image {}", filename);
        return false;
    }
    return true;
}
bool save_jpg_image(const string& filename, const image4b& img) {
    if (!stbi_write_jpg(
            filename.c_str(), img.width, img.height, 4, data(img), 75)) {
        log_io_error("error saving image {}", filename);
        return false;
    }
    return true;
}
bool save_tga_image(const string& filename, const image4b& img) {
    if (!stbi_write_tga(filename.c_str(), img.width, img.height, 4, data(img))) {
        log_io_error("error saving image {}", filename);
        return false;
    }
    return true;
}
bool save_bmp_image(const string& filename, const image4b& img) {
    if (!stbi_write_bmp(filename.c_str(), img.width, img.height, 4, data(img))) {
        log_io_error("error saving image {}", filename);
        return false;
    }
    return true;
}
bool save_hdr_image(const string& filename, const image4f& img) {
    if (!stbi_write_hdr(
            filename.c_str(), img.width, img.height, 4, (float*)data(img))) {
        log_io_error("error saving image {}", filename);
        return false;
    }
    return true;
}

// load an image using stbi library
bool load_stb_image_from_memory(const byte* data, int data_size, image4b& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = (vec4b*)stbi_load_from_memory(
        data, data_size, &width, &height, &ncomp, 4);
    if (!pixels) {
        log_io_error("error loading in-memory image");
        return false;
    }
    img = make_image(width, height, pixels);
    free(pixels);
    return true;
}
bool load_stbi_image_from_memory(const byte* data, int data_size, image4f& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = (vec4f*)stbi_loadf_from_memory(
        data, data_size, &width, &height, &ncomp, 4);
    if (!pixels) {
        log_io_error("error loading in-memory image {}");
        return false;
    }
    img = make_image(width, height, pixels);
    free(pixels);
    return true;
}

bool apply_json_procedural(const json& js, image4f& img) {
    auto type = get_json_value(js, "type", ""s);
    auto width = get_json_value(js, "width", 512);
    auto height = get_json_value(js, "height", 512);
    if(type == "") {
        img = make_image(width, height, zero4f);
    } else if (type == "grid") {
        img = make_grid_image(width, height,
            get_json_value(js, "tile", 8),
            get_json_value(js, "c0", vec4f{0.2f, 0.2f, 0.2f, 1}),
            get_json_value(js, "c1", vec4f{0.8f, 0.8f, 0.8f, 1}));
    } else if (type == "checker") {
        img = make_checker_image(width, height,
            get_json_value(js, "tile", 8),
            get_json_value(js, "c0", vec4f{0.2f, 0.2f, 0.2f, 1}),
            get_json_value(js, "c1", vec4f{0.8f, 0.8f, 0.8f, 1}));
    } else if (type == "bump") {
        img = make_bumpdimple_image(
            width, height, get_json_value(js, "tile", 8));
    } else if (type == "uvramp") {
        img = make_uvramp_image(width, height);
    } else if (type == "gammaramp") {
        img = make_gammaramp_image(width, height);
    } else if (type == "blackbodyramp") {
        img = make_blackbodyramp_image(width, height);
    } else if (type == "uvgrid") {
        img = make_uvgrid_image(width, height);
    } else if (type == "sky") {
        if (width < height * 2) width = height * 2;
        img = make_sunsky_image(width, height,
            get_json_value(js, "sun_angle", pif / 4),
            get_json_value(js, "turbidity", 3.0f),
            get_json_value(js, "has_sun", false),
            get_json_value(js, "sun_angle_scale", 1.0f),
            get_json_value(js, "sun_emission_scale", 1.0f),
            get_json_value(js, "ground_albedo", vec3f{0.7f, 0.7f, 0.7f}));
    } else if (type == "noise") {
        img = make_noise_image(width, height,
            get_json_value(js, "scale", 1.0f), get_json_value(js, "wrap", true));
    } else if (type == "fbm") {
        img = make_fbm_image(width, height,
            get_json_value(js, "scale", 1.0f),
            get_json_value(js, "lacunarity", 2.0f),
            get_json_value(js, "gain", 0.5f), get_json_value(js, "octaves", 6),
            get_json_value(js, "wrap", true));
    } else if (type == "ridge") {
        img = make_ridge_image(width, height,
            get_json_value(js, "scale", 1.0f),
            get_json_value(js, "lacunarity", 2.0f),
            get_json_value(js, "gain", 0.5f), get_json_value(js, "offset", 1.0f),
            get_json_value(js, "octaves", 6), get_json_value(js, "wrap", true));
    } else if (type == "turbulence") {
        img = make_turbulence_image(width, height,
            get_json_value(js, "scale", 1.0f),
            get_json_value(js, "lacunarity", 2.0f),
            get_json_value(js, "gain", 0.5f), get_json_value(js, "octaves", 6),
            get_json_value(js, "wrap", true));
    } else {
        log_error("unknown image type {}", type);
        return false;
    }
    if (get_json_value(js, "bump_to_normal", false)) {
        img = bump_to_normal_map(
            img, get_json_value(js, "bump_scale", 1.0f));
    }
    return true;
}

bool apply_json_procedural(const json& js, image4b& img) {
    auto imgf = image4f{};
    if(!apply_json_procedural(js, imgf)) return false;
    auto srgb = get_json_value(js, "srgb", true);
    if(srgb) imgf = linear_to_srgb(imgf);
    img = float_to_byte(imgf);
    return true;
}

// load a JSON image
bool load_json_image(const string& filename, image4f& img) {
    auto js = json();
    if(!load_json(filename, js)) return false;
    if(!apply_json_procedural(js, img)) return false;
    return true;
}
bool load_json_image(const string& filename, image4b& img) {
    auto js = json();
    if(!load_json(filename, js)) return false;
    if(!apply_json_procedural(js, img)) return false;
    return true;
}

// check hdr extensions
bool is_hdr_filename(const string& filename) {
    auto ext = get_extension(filename);
    return ext == "hdr" || ext == "exr" || ext == "pfm";
}

// Loads an hdr image.
bool load_image_nolog(const string& filename, image4f& img) {
    auto ext = get_extension(filename);
    if (ext == "exr" || ext == "EXR") {
        return load_exr_image(filename, img);
    } else if (ext == "pfm" || ext == "PFM") {
        return load_pfm_image(filename, img);
    } else if (ext == "hdr" || ext == "HDR") {
        return load_stb_image(filename, img);
    } else if (ext == "png" || ext == "PNG") {
        auto img8 = image4b{};
        if (!load_stb_image(filename, img8)) return false;
        img = srgb_to_linear(byte_to_float(img8));
        return true;
    } else if (ext == "jpg" || ext == "JPG") {
        auto img8 = image4b{};
        if (!load_stb_image(filename, img8)) return false;
        img = srgb_to_linear(byte_to_float(img8));
        return true;
    } else if (ext == "tga" || ext == "TGA") {
        auto img8 = image4b{};
        if (!load_stb_image(filename, img8)) return false;
        img = srgb_to_linear(byte_to_float(img8));
        return true;
    } else if (ext == "bmp" || ext == "BMP") {
        auto img8 = image4b{};
        if (!load_stb_image(filename, img8)) return false;
        img = srgb_to_linear(byte_to_float(img8));
        return true;
    } else if (ext == "json" || ext == "JSON") {
        return load_json_image(filename, img);
    } else {
        log_io_error("unsupported image format {}", ext);
        return false;
    }
}
bool load_image(const string& filename, image4f& img) {
    auto scope = log_trace_scoped("loading image {}", filename);
    return load_image_nolog(filename, img);
}

// Saves an hdr image.
bool save_image_nolog(const string& filename, const image4f& img) {
    auto ext = get_extension(filename);
    if (ext == "png" || ext == "PNG") {
        return save_png_image(filename, float_to_byte(linear_to_srgb(img)));
    } else if (ext == "jpg" || ext == "JPG") {
        return save_jpg_image(filename, float_to_byte(linear_to_srgb(img)));
    } else if (ext == "tga" || ext == "TGA") {
        return save_tga_image(filename, float_to_byte(linear_to_srgb(img)));
    } else if (ext == "bmp" || ext == "BMP") {
        return save_bmp_image(filename, float_to_byte(linear_to_srgb(img)));
    } else if (ext == "hdr" || ext == "HDR") {
        return save_hdr_image(filename, img);
    } else if (ext == "pfm" || ext == "PFM") {
        return save_pfm_image(filename, img);
    } else if (ext == "exr" || ext == "EXR") {
        return save_exr_image(filename, img);
    } else {
        log_io_error("unsupported image format {}", ext);
        return false;
    }
}
bool save_image(const string& filename, const image4f& img) {
    auto scope = log_trace_scoped("saving image {}", filename);
    return save_image_nolog(filename, img);
}

// Loads an hdr image.
bool load_image_from_memory_nolog(const byte* data, int data_size, image4f& img) {
    return load_stbi_image_from_memory(data, data_size, img);
}
bool load_image_from_memory(const byte* data, int data_size, image4f& img) {
    auto scope = log_trace_scoped("loading image in memory");
    return load_image_from_memory_nolog(data, data_size, img);
}

// Loads an hdr image.
bool load_image_nolog(const string& filename, image4b& img) {
    auto ext = get_extension(filename);
    if (ext == "exr" || ext == "EXR") {
        auto imgf = image4f{};
        if (!load_exr_image(filename, imgf)) return false;
        img = float_to_byte(linear_to_srgb(imgf));
        return true;
    } else if (ext == "pfm" || ext == "PFM") {
        auto imgf = image4f{};
        if (!load_pfm_image(filename, imgf)) return false;
        img = float_to_byte(linear_to_srgb(imgf));
        return true;
    } else if (ext == "hdr" || ext == "HDR") {
        auto imgf = image4f{};
        if (!load_stb_image(filename, imgf)) return false;
        img = float_to_byte(linear_to_srgb(imgf));
        return true;
    } else if (ext == "png" || ext == "PNG") {
        return load_stb_image(filename, img);
    } else if (ext == "jpg" || ext == "JPG") {
        return load_stb_image(filename, img);
    } else if (ext == "tga" || ext == "TGA") {
        return load_stb_image(filename, img);
    } else if (ext == "bmp" || ext == "BMP") {
        return load_stb_image(filename, img);
    } else if (ext == "json" || ext == "JSON") {
        return load_json_image(filename, img);
    } else {
        log_io_error("unsupported image format {}", ext);
        return false;
    }
}
bool load_image(const string& filename, image4b& img) {
    auto scope = log_trace_scoped("loading image {}", filename);
    return load_image_nolog(filename, img);
}

// Saves an ldr image.
bool save_image_nolog(const string& filename, const image4b& img) {
    auto ext = get_extension(filename);
    if (ext == "png" || ext == "PNG") {
        return save_png_image(filename, img);
    } else if (ext == "jpg" || ext == "JPG") {
        return save_jpg_image(filename, img);
    } else if (ext == "tga" || ext == "TGA") {
        return save_tga_image(filename, img);
    } else if (ext == "bmp" || ext == "BMP") {
        return save_bmp_image(filename, img);
    } else if (ext == "hdr" || ext == "HDR") {
        return save_hdr_image(filename, srgb_to_linear(byte_to_float(img)));
    } else if (ext == "pfm" || ext == "PFM") {
        return save_pfm_image(filename, srgb_to_linear(byte_to_float(img)));
    } else if (ext == "exr" || ext == "EXR") {
        return save_exr_image(filename, srgb_to_linear(byte_to_float(img)));
    } else {
        log_io_error("unsupported image format {}", ext);
        return false;
    }
}
bool save_image(const string& filename, const image4b& img) {
    auto scope = log_trace_scoped("saving image {}", filename);
    return save_image_nolog(filename, img);
}

// Loads an ldr image.
bool load_image_from_memory_nolog(const byte* data, int data_size, image4b& img) {
    return load_stb_image_from_memory(data, data_size, img);
}
bool load_image_from_memory(const byte* data, int data_size, image4b& img) {
    auto scope = log_trace_scoped("loading image in memory");
    return load_image_from_memory_nolog(data, data_size, img);
}

// Convenience helper that saves an HDR images as wither a linear HDR file or
// a tonemapped LDR file depending on file name
bool save_tonemapped_image(const string& filename, const image4f& hdr,
    float exposure, bool filmic, bool srgb) {
    if (is_hdr_filename(filename)) {
        return save_image(filename, hdr);
    } else {
        return save_image(filename,
            float_to_byte(tonemap_image(hdr, exposure, filmic, srgb)));
    }
}

// Resize image.
image4f resize_image(const image4f& img, int width, int height) {
    if (width == 0 && height == 0) {
        log_error("bad image size in resize_image");
        return {};
    }
    if (height == 0) {
        height = (int)round(width * (float)img.height / (float)img.width);
    } else if (width == 0) {
        width = (int)round(height * (float)img.width / (float)img.height);
    }
    auto res_img = make_image(width, height, zero4f);
    stbir_resize_float_generic((float*)data(img), img.width, img.height,
        sizeof(vec4f) * img.width, (float*)data(res_img), res_img.width,
        res_img.height, sizeof(vec4f) * res_img.width, 4, 3, 0, STBIR_EDGE_CLAMP,
        STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
    return img;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR VOLUME IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Loads volume data from binary format.
bool load_volume_nolog(const string& filename, volume1f& vol) {
    auto fs = open(filename, "r");
    if (!fs) return false;
    if (!read_value(fs, vol.width)) return false;
    if (!read_value(fs, vol.height)) return false;
    if (!read_value(fs, vol.depth)) return false;
    vol.voxels.resize(vol.width * vol.height * vol.depth);
    if (!read_values(fs, size(vol), data(vol))) return false;
    return true;
}
bool load_volume(const string& filename, volume1f& vol) {
    auto scope = log_trace_scoped("loading volume {}", filename);
    return load_volume_nolog(filename, vol);
}

// Saves volume data in binary format.
bool save_volume_nolog(const string& filename, const volume1f& vol) {
    auto fs = open(filename, "w");
    if (!fs) return false;
    if (!write_value(fs, vol.width)) return false;
    if (!write_value(fs, vol.height)) return false;
    if (!write_value(fs, vol.depth)) return false;
    if (!write_values(fs, size(vol), data(vol))) return false;
    return true;
}
bool save_volume(const string& filename, const volume1f& vol) {
    auto scope = log_trace_scoped("saving volume {}", filename);
    return save_volume_nolog(filename, vol);
}

}  // namespace yocto
