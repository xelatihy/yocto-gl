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
        if (fread(pixels.data() + j * nrow, sizeof(float), nrow, fs.fs) != nrow) {
            return {};
        }
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
bool load_pfm_image(const string& filename, image<vec4f>& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = load_pfm(filename.c_str(), width, height, ncomp, 4);
    if (pixels.empty()) {
        log_io_error("error loading image {}", filename);
        return false;
    }
    img = {{width, height}, (vec4f*)pixels.data()};
    return true;
}
bool save_pfm_image(const string& filename, const image<vec4f>& img) {
    if (!save_pfm(filename.c_str(), img.size().x, img.size().y, 4,
            (float*)img.data())) {
        log_io_error("error saving image {}", filename);
        return false;
    }
    return true;
}

// load exr image weith tiny exr
bool load_exr_image(const string& filename, image<vec4f>& img) {
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
    img = {{width, height}, pixels};
    free(pixels);
    return true;
}
bool save_exr_image(const string& filename, const image<vec4f>& img) {
    if (!SaveEXR((float*)img.data(), img.size().x, img.size().y, 4,
            filename.c_str())) {
        log_io_error("error saving image {}", filename);
        return false;
    }
    return true;
}

// load an image using stbi library
bool load_stb_image(const string& filename, image<vec4b>& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = (vec4b*)stbi_load(filename.c_str(), &width, &height, &ncomp, 4);
    if (!pixels) {
        log_io_error("error loading image {}", filename);
        return false;
    }
    img = {{width, height}, pixels};
    free(pixels);
    return true;
}
bool load_stb_image(const string& filename, image<vec4f>& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = (vec4f*)stbi_loadf(
        filename.c_str(), &width, &height, &ncomp, 4);
    if (!pixels) {
        log_io_error("error loading image {}", filename);
        return false;
    }
    img = {{width, height}, pixels};
    free(pixels);
    return true;
}

// save an image with stbi
bool save_png_image(const string& filename, const image<vec4b>& img) {
    if (!stbi_write_png(filename.c_str(), img.size().x, img.size().y, 4,
            img.data(), img.size().x * 4)) {
        log_io_error("error saving image {}", filename);
        return false;
    }
    return true;
}
bool save_jpg_image(const string& filename, const image<vec4b>& img) {
    if (!stbi_write_jpg(
            filename.c_str(), img.size().x, img.size().y, 4, img.data(), 75)) {
        log_io_error("error saving image {}", filename);
        return false;
    }
    return true;
}
bool save_tga_image(const string& filename, const image<vec4b>& img) {
    if (!stbi_write_tga(
            filename.c_str(), img.size().x, img.size().y, 4, img.data())) {
        log_io_error("error saving image {}", filename);
        return false;
    }
    return true;
}
bool save_bmp_image(const string& filename, const image<vec4b>& img) {
    if (!stbi_write_bmp(
            filename.c_str(), img.size().x, img.size().y, 4, img.data())) {
        log_io_error("error saving image {}", filename);
        return false;
    }
    return true;
}
bool save_hdr_image(const string& filename, const image<vec4f>& img) {
    if (!stbi_write_hdr(filename.c_str(), img.size().x, img.size().y, 4,
            (float*)img.data())) {
        log_io_error("error saving image {}", filename);
        return false;
    }
    return true;
}

// load an image using stbi library
bool load_stb_image_from_memory(
    const byte* data, int data_size, image<vec4b>& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = (vec4b*)stbi_load_from_memory(
        data, data_size, &width, &height, &ncomp, 4);
    if (!pixels) {
        log_io_error("error loading in-memory image");
        return false;
    }
    img = {{width, height}, pixels};
    free(pixels);
    return true;
}
bool load_stbi_image_from_memory(
    const byte* data, int data_size, image<vec4f>& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = (vec4f*)stbi_loadf_from_memory(
        data, data_size, &width, &height, &ncomp, 4);
    if (!pixels) {
        log_io_error("error loading in-memory image {}");
        return false;
    }
    img = {{width, height}, pixels};
    free(pixels);
    return true;
}

// check hdr extensions
bool is_hdr_filename(const string& filename) {
    auto ext = get_extension(filename);
    return ext == "hdr" || ext == "exr" || ext == "pfm";
}

// Loads an hdr image.
bool load_image_nolog(const string& filename, image<vec4f>& img) {
    auto ext = get_extension(filename);
    if (ext == "exr" || ext == "EXR") {
        return load_exr_image(filename, img);
    } else if (ext == "pfm" || ext == "PFM") {
        return load_pfm_image(filename, img);
    } else if (ext == "hdr" || ext == "HDR") {
        return load_stb_image(filename, img);
    } else if (ext == "png" || ext == "PNG") {
        auto img8 = image<vec4b>{};
        if (!load_stb_image(filename, img8)) return false;
        img = srgb_to_linear(byte_to_float(img8));
        return true;
    } else if (ext == "jpg" || ext == "JPG") {
        auto img8 = image<vec4b>{};
        if (!load_stb_image(filename, img8)) return false;
        img = srgb_to_linear(byte_to_float(img8));
        return true;
    } else if (ext == "tga" || ext == "TGA") {
        auto img8 = image<vec4b>{};
        if (!load_stb_image(filename, img8)) return false;
        img = srgb_to_linear(byte_to_float(img8));
        return true;
    } else if (ext == "bmp" || ext == "BMP") {
        auto img8 = image<vec4b>{};
        if (!load_stb_image(filename, img8)) return false;
        img = srgb_to_linear(byte_to_float(img8));
        return true;
    } else {
        log_io_error("unsupported image format {}", ext);
        return false;
    }
}
bool load_image(const string& filename, image<vec4f>& img) {
    auto scope = log_trace_scoped("loading image {}", filename);
    return load_image_nolog(filename, img);
}

// Saves an hdr image.
bool save_image_nolog(const string& filename, const image<vec4f>& img) {
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
bool save_image(const string& filename, const image<vec4f>& img) {
    auto scope = log_trace_scoped("saving image {}", filename);
    return save_image_nolog(filename, img);
}

// Loads an hdr image.
bool load_image_from_memory_nolog(
    const byte* data, int data_size, image<vec4f>& img) {
    return load_stbi_image_from_memory(data, data_size, img);
}
bool load_image_from_memory(const byte* data, int data_size, image<vec4f>& img) {
    auto scope = log_trace_scoped("loading image in memory");
    return load_image_from_memory_nolog(data, data_size, img);
}

// Loads an hdr image.
bool load_image_nolog(const string& filename, image<vec4b>& img) {
    auto ext = get_extension(filename);
    if (ext == "exr" || ext == "EXR") {
        auto imgf = image<vec4f>{};
        if (!load_exr_image(filename, imgf)) return false;
        img = float_to_byte(linear_to_srgb(imgf));
        return true;
    } else if (ext == "pfm" || ext == "PFM") {
        auto imgf = image<vec4f>{};
        if (!load_pfm_image(filename, imgf)) return false;
        img = float_to_byte(linear_to_srgb(imgf));
        return true;
    } else if (ext == "hdr" || ext == "HDR") {
        auto imgf = image<vec4f>{};
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
    } else {
        log_io_error("unsupported image format {}", ext);
        return false;
    }
}
bool load_image(const string& filename, image<vec4b>& img) {
    auto scope = log_trace_scoped("loading image {}", filename);
    return load_image_nolog(filename, img);
}

// Saves an ldr image.
bool save_image_nolog(const string& filename, const image<vec4b>& img) {
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
bool save_image(const string& filename, const image<vec4b>& img) {
    auto scope = log_trace_scoped("saving image {}", filename);
    return save_image_nolog(filename, img);
}

// Loads an ldr image.
bool load_image_from_memory_nolog(
    const byte* data, int data_size, image<vec4b>& img) {
    return load_stb_image_from_memory(data, data_size, img);
}
bool load_image_from_memory(const byte* data, int data_size, image<vec4b>& img) {
    auto scope = log_trace_scoped("loading image in memory");
    return load_image_from_memory_nolog(data, data_size, img);
}

// Convenience helper that saves an HDR images as wither a linear HDR file or
// a tonemapped LDR file depending on file name
bool save_tonemapped_image(const string& filename, const image<vec4f>& hdr,
    float exposure, bool filmic, bool srgb) {
    if (is_hdr_filename(filename)) {
        return save_image(filename, hdr);
    } else {
        return save_image(filename,
            float_to_byte(tonemap_image(hdr, exposure, filmic, srgb)));
    }
}

// Resize image.
image<vec4f> resize_image(const image<vec4f>& img, const vec2i& size) {
    if (size == zero2i) {
        log_error("bad image size in resize_image");
    }
    auto res_img = image<vec4f>{get_image_size(size, get_image_aspect(img))};
    stbir_resize_float_generic((float*)img.data(), img.size().x, img.size().y,
        sizeof(vec4f) * img.size().x, (float*)res_img.data(), res_img.size().x,
        res_img.size().y, sizeof(vec4f) * res_img.size().x, 4, 3, 0,
        STBIR_EDGE_CLAMP, STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
    return img;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR VOLUME IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Loads volume data from binary format.
bool load_volume_nolog(const string& filename, volume<float>& vol) {
    auto fs = open(filename, "r");
    if (!fs) return false;
    auto size = zero3i;
    if (!read_value(fs, size)) return false;
    vol.resize(size);
    if (!read_values(fs, size.x * size.y * size.z, vol.data())) return false;
    return true;
}
bool load_volume(const string& filename, volume<float>& vol) {
    auto scope = log_trace_scoped("loading volume {}", filename);
    return load_volume_nolog(filename, vol);
}

// Saves volume data in binary format.
bool save_volume_nolog(const string& filename, const volume<float>& vol) {
    auto fs = open(filename, "w");
    if (!fs) return false;
    auto size = vol.size();
    if (!write_value(fs, size)) return false;
    if (!write_values(fs, size.x * size.y * size.z, vol.data())) return false;
    return true;
}
bool save_volume(const string& filename, const volume<float>& vol) {
    auto scope = log_trace_scoped("saving volume {}", filename);
    return save_volume_nolog(filename, vol);
}

}  // namespace yocto
