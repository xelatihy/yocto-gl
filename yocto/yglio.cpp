//
// Implementation for Yocto/GL Input and Output functions.
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

//
//
// LICENSE OF INCLUDED CODE FOR BASE64 (base64.h, base64.cpp)
//
// Copyright (C) 2004-2008 René Nyffenegger
//
// This source code is provided 'as-is', without any express or implied
// warranty. In no event will the author be held liable for any damages
// arising from the use of this software.
//
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
//
// 1. The origin of this source code must not be misrepresented; you must not
// claim that you wrote the original source code. If you use this source code
// in a product, an acknowledgment in the product documentation would be
// appreciated but is not required.
//
// 2. Altered source versions must be plainly marked as such, and must not be
// misrepresented as being the original source code.
//
// 3. This notice may not be removed or altered from any source distribution.
//
// René Nyffenegger rene.nyffenegger@adp-gmbh.ch
//
//

#include "yglio.h"

#include <cstdlib>
#include <deque>
#include <map>
#include <regex>

#include <array>
#include <climits>

#include "ext/json.hpp"

#if YGL_HAPPLY
#include "ext/happly.h"
#endif

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
// IMPLEMENTATION OF PATH UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

string normalize_path(const string& filename_) {
    auto filename = filename_;
    for (auto& c : filename)
        if (c == '\\') c = '/';
    if (filename.size() > 1 && filename[0] == '/' && filename[1] == '/') {
        log_error("absolute paths are not supported");
        return filename_;
    }
    if (filename.size() > 3 && filename[1] == ':' && filename[2] == '/' &&
        filename[3] == '/') {
        log_error("absolute paths are not supported");
        return filename_;
    }
    auto pos = (size_t)0;
    while ((pos = filename.find("//")) != filename.npos)
        filename = filename.substr(0, pos) + filename.substr(pos + 1);
    return filename;
}

// Get directory name (not including '/').
string get_dirname(const string& filename_) {
    auto filename = normalize_path(filename_);
    auto pos      = filename.rfind('/');
    if (pos == string::npos) return "";
    return filename.substr(0, pos);
}

// Get extension (not including '.').
string get_extension(const string& filename_) {
    auto filename = normalize_path(filename_);
    auto pos      = filename.rfind('.');
    if (pos == string::npos) return "";
    return filename.substr(pos + 1);
}

// Get filename without directory.
string get_filename(const string& filename_) {
    auto filename = normalize_path(filename_);
    auto pos      = filename.rfind('/');
    if (pos == string::npos) return "";
    return filename.substr(pos + 1);
}

// Replace extension.
string replace_extension(const string& filename_, const string& ext_) {
    auto filename = normalize_path(filename_);
    auto ext      = normalize_path(ext_);
    if (ext.at(0) == '.') ext = ext.substr(1);
    auto pos = filename.rfind('.');
    if (pos == string::npos) return filename;
    return filename.substr(0, pos) + "." + ext;
}

// Check if a file can be opened for reading.
bool exists_file(const string& filename) {
    auto f = fopen(filename.c_str(), "r");
    if (!f) return false;
    fclose(f);
    return true;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF FILE READING
// -----------------------------------------------------------------------------
namespace ygl {

// log io error
template <typename... Args>
void log_io_error(const string& fmt, const Args&... args) {
    log_error(fmt, args...);
}

// File stream wrapper
struct file_stream {
    string filename = "";
    string mode     = "";
    FILE*  fs       = nullptr;

    file_stream()                   = default;
    file_stream(const file_stream&) = delete;
    file_stream& operator=(const file_stream&) = delete;
    file_stream(file_stream&&)                 = default;
    file_stream& operator=(file_stream&&) = default;

    ~file_stream() {
        if (fs) {
            fclose(fs);
            fs = nullptr;
        }
    }

    operator bool() const { return fs; }
};

// Opens a file
file_stream open(const string& filename, const string& mode) {
    auto fs = fopen(filename.c_str(), mode.c_str());
    if (!fs) {
        log_io_error("cannot open {}", filename);
        return {};
    }
    return {filename, mode, fs};
}

// Close a file
bool close(file_stream& fs) {
    if (!fs) {
        log_io_error("cannot close {}", fs.filename);
        return false;
    }
    fclose(fs.fs);
    fs.fs = nullptr;
    return true;
}

// Gets the length of a file
size_t get_length(file_stream& fs) {
    if (!fs) return 0;
    fseek(fs.fs, 0, SEEK_END);
    auto fsize = ftell(fs.fs);
    fseek(fs.fs, 0, SEEK_SET);
    return fsize;
}

// Print to file
bool write_text(file_stream& fs, const string& str) {
    if (!fs) return false;
    if (fprintf(fs.fs, "%s", str.c_str()) < 0) {
        log_io_error("cannot write to {}", fs.filename);
        return false;
    }
    return true;
}

// Write to file
template <typename T>
bool write_value(file_stream& fs, const T& val) {
    if (!fs) return false;
    if (fwrite(&val, sizeof(T), 1, fs.fs) != 1) {
        log_io_error("cannot write to {}", fs.filename);
        return false;
    }
    return true;
}

// Write to file
template <typename T>
bool write_values(file_stream& fs, const vector<T>& vals) {
    if (!fs) return false;
    if (fwrite(vals.data(), sizeof(T), vals.size(), fs.fs) != vals.size()) {
        log_io_error("cannot write to {}", fs.filename);
        return false;
    }
    return true;
}

// Print shortcut
template <typename... Args>
bool print(file_stream& fs, const string& fmt, const Args&... args) {
    if (!fs) return false;
    return write_text(fs, format(fmt, args...));
}

// Read binary data to fill the whole buffer
bool read_line(file_stream& fs, string& val) {
    if (!fs) return false;
    // TODO: make lkne as large as possible
    val = "";
    char buf[4096];
    if (!fgets(buf, 4096, fs.fs)) return false;
    val = string(buf);
    return true;
}

// Read binary data to fill the whole buffer
template <typename T>
bool read_value(file_stream& fs, T& val) {
    if (!fs) return false;
    if (fread(&val, sizeof(T), 1, fs.fs) != 1) {
        log_io_error("cannot read from {}", fs.filename);
        return false;
    }
    return true;
}

// Read binary data to fill the whole buffer
template <typename T>
bool read_values(file_stream& fs, vector<T>& vals) {
    if (!fs) return false;
    if (fread(vals.data(), sizeof(T), vals.size(), fs.fs) != vals.size()) {
        log_io_error("cannot read from {}", fs.filename);
        return false;
    }
    return true;
}

// Load a text file
string load_text(const string& filename) {
    // https://stackoverflow.com/questions/2602013/read-whole-ascii-file-into-c-stdstring
    auto fs = open(filename, "rb");
    if (!fs) return {};
    auto buf = vector<char>(get_length(fs));
    if (!read_values(fs, buf)) return {};
    return string{buf.begin(), buf.end()};
}

// Save a text file
bool save_text(const string& filename, const string& str) {
    auto fs = open(filename, "wt");
    if (!fs) return false;
    if (!write_text(fs, str)) return false;
    return true;
}

// Load a binary file
vector<byte> load_binary(const string& filename) {
    // https://stackoverflow.com/questions/2602013/read-whole-ascii-file-into-c-stdstring
    auto fs = open(filename, "rb");
    if (!fs) return {};
    auto data = vector<byte>(get_length(fs));
    if (!read_values(fs, data)) return {};
    return data;
}

// Save a binary file
bool save_binary(const string& filename, const vector<byte>& data) {
    auto fs = open(filename.c_str(), "wb");
    if (!fs) return false;
    if (!write_values(fs, data)) return false;
    return true;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// JSON UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Json alias
using json = nlohmann::json;

// Load a JSON object
json load_json(const string& filename) {
    auto texture = load_text(filename);
    if (texture.empty()) return {};
    try {
        return json::parse(texture.begin(), texture.end());
    } catch (...) {
        log_io_error("could not parse json {}", filename);
        return {};
    }
}

// Save a JSON object
bool save_json(const string& filename, const json& js) {
    auto str = ""s;
    try {
        str = js.dump(4);
    } catch (...) {
        log_io_error("could not dump json {}", filename);
        return false;
    }
    return save_text(filename, str);
}

template <typename T, int N>
inline void to_json(json& js, const vec<T, N>& val) {
    js = (const std::array<T, N>&)val;
}
template <typename T, int N>
inline void from_json(const json& js, vec<T, N>& val) {
    (std::array<T, N>&)val = js.get<std::array<T, N>>();
}

template <typename T, int N>
inline void to_json(json& js, const frame<T, N>& val) {
    js = (const std::array<vec<T, N>, N + 1>&)val;
}
template <typename T, int N>
inline void from_json(const json& js, frame<T, N>& val) {
    (std::array<vec<T, N>, N + 1>&)val = js.get<std::array<vec<T, N>, N + 1>>();
}

template <typename T, int N, int M>
inline void to_json(json& js, const mat<T, N, M>& val) {
    js = (const std::array<vec<T, N>, M>&)val;
}
template <typename T, int N, int M>
inline void from_json(const json& js, mat<T, N, M>& val) {
    (std::array<vec<T, N>, M>&)val = js.get<std::array<vec<T, N>, M>>();
}

template <typename T>
inline void to_json(json& js, const bbox<T, 1>& val) {
    js = (const std::array<T, 2>&)val;
}
template <typename T>
inline void from_json(const json& js, bbox<T, 1>& val) {
    (std::array<T, 2>&)val = js.get<std::array<T, 2>>();
}
template <typename T, int N>
inline void to_json(json& js, const bbox<T, N>& val) {
    js = (const std::array<vec<T, N>, 2>&)val;
}
template <typename T, int N>
inline void from_json(const json& js, bbox<T, N>& val) {
    (std::array<vec<T, N>, 2>&)val = js.get<std::array<vec<T, N>, 2>>();
}

template <typename T>
inline void to_json(json& js, const image<T>& val) {
    js           = json::object();
    js["width"]  = val.width;
    js["height"] = val.height;
    js["pixels"] = val.pixels;
}
template <typename T>
inline void from_json(const json& js, image<T>& val) {
    auto width  = js.at("width").get<int>();
    auto height = js.at("height").get<int>();
    auto pixels = js.at("pixels").get<vector<T>>();
    val         = image<T>{width, height, pixels.data()};
}
template <typename T>
inline void to_json(json& js, const volume<T>& val) {
    js           = json::object();
    js["width"]  = val.width;
    js["height"] = val.height;
    js["depth"]  = val.depth;
    js["voxels"] = val.voxels;
}
template <typename T>
inline void from_json(const json& js, volume<T>& val) {
    auto width  = js.at("width").get<int>();
    auto height = js.at("height").get<int>();
    auto depth  = js.at("depth").get<int>();
    auto voxels = js.at("voxels").get<vector<T>>();
    val         = volume<T>{width, height, depth, voxels.data()};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGEIO
// -----------------------------------------------------------------------------
namespace ygl {

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

    // buffer
    char buf[256];
    auto toks = vector<string>();

    // read magic
    if (!fgets(buf, 256, fs)) return nullptr;
    toks = split_string(buf);
    if (toks[0] == "Pf")
        *nc = 1;
    else if (toks[0] == "PF")
        *nc = 3;
    else
        return nullptr;

    // read w, h
    if (!fgets(buf, 256, fs)) return nullptr;
    toks = split_string(buf);
    *w   = atoi(toks[0].c_str());
    *h   = atoi(toks[1].c_str());

    // read scale
    if (!fgets(buf, 256, fs)) return nullptr;
    toks   = split_string(buf);
    auto s = atof(toks[0].c_str());

    // read the data (flip y)
    auto npixels = (*w) * (*h);
    auto nvalues = (*w) * (*h) * (*nc);
    auto nrow    = (*w) * (*nc);
    auto pixels  = new float[nvalues];
    for (auto j = *h - 1; j >= 0; j--) {
        if (fread(pixels + j * nrow, sizeof(float), nrow, fs) != nrow) {
            delete[] pixels;
            return nullptr;
        }
    }

    // done reading
    fclose(fs);

    // endian conversion
    if (s > 0) {
        for (auto i = 0; i < nvalues; ++i) {
            auto dta = (uint8_t*)(pixels + i);
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
image<vec4f> load_pfm_image4f(const string& filename) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = (vec4f*)load_pfm(filename.c_str(), &width, &height, &ncomp, 4);
    if (!pixels) {
        log_io_error("error loading image {}", filename);
        return {};
    }
    auto img = image<vec4f>{width, height, pixels};
    delete[] pixels;
    return img;
}
bool save_pfm_image4f(const string& filename, const image<vec4f>& img) {
    if (!save_pfm(filename.c_str(), img.width, img.height, 4,
            (float*)img.pixels.data())) {
        log_io_error("error saving image {}", filename);
        return false;
    }
    return true;
}

// load exr image weith tiny exr
image<vec4f> load_exr_image4f(const string& filename) {
    auto width = 0, height = 0;
    auto pixels = (vec4f*)nullptr;
    if (LoadEXR((float**)&pixels, &width, &height, filename.c_str(), nullptr) <
        0) {
        log_io_error("error loading image {}", filename);
        return {};
    }
    if (!pixels) {
        log_io_error("error loading image {}", filename);
        return {};
    }
    auto img = image<vec4f>{width, height, pixels};
    free(pixels);
    return img;
}
bool save_exr_image4f(const string& filename, const image<vec4f>& img) {
    if (!SaveEXR((float*)img.pixels.data(), img.width, img.height, 4,
            filename.c_str())) {
        log_io_error("error saving image {}", filename);
        return false;
    }
    return true;
}

// load an image using stbi library
image<vec4b> load_stb_image4b(const string& filename) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = (vec4b*)stbi_load(filename.c_str(), &width, &height, &ncomp, 4);
    if (!pixels) {
        log_io_error("error loading image {}", filename);
        return {};
    }
    auto img = image<vec4b>{width, height, pixels};
    free(pixels);
    return img;
}
image<vec4f> load_stb_image4f(const string& filename) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = (vec4f*)stbi_loadf(
        filename.c_str(), &width, &height, &ncomp, 4);
    if (!pixels) {
        log_io_error("error loading image {}", filename);
        return {};
    }
    auto img = image<vec4f>{width, height, pixels};
    free(pixels);
    return img;
}

// save an image with stbi
bool save_png_image4b(const string& filename, const image<vec4b>& img) {
    if (!stbi_write_png(filename.c_str(), img.width, img.height, 4,
            img.pixels.data(), img.width * 4)) {
        log_io_error("error saving image {}", filename);
        return false;
    }
    return true;
}
bool save_jpg_image4b(const string& filename, const image<vec4b>& img) {
    if (!stbi_write_jpg(filename.c_str(), img.width, img.height, 4,
            img.pixels.data(), 75)) {
        log_io_error("error saving image {}", filename);
        return false;
    }
    return true;
}
bool save_tga_image4b(const string& filename, const image<vec4b>& img) {
    if (!stbi_write_tga(
            filename.c_str(), img.width, img.height, 4, img.pixels.data())) {
        log_io_error("error saving image {}", filename);
        return false;
    }
    return true;
}
bool save_bmp_image4b(const string& filename, const image<vec4b>& img) {
    if (!stbi_write_bmp(
            filename.c_str(), img.width, img.height, 4, img.pixels.data())) {
        log_io_error("error saving image {}", filename);
        return false;
    }
    return true;
}
bool save_hdr_image4f(const string& filename, const image<vec4f>& img) {
    if (!stbi_write_hdr(filename.c_str(), img.width, img.height, 4,
            (float*)img.pixels.data())) {
        log_io_error("error saving image {}", filename);
        return false;
    }
    return true;
}

// load an image using stbi library
image<vec4b> load_stb_image4b_from_memory(const byte* data, int data_size) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = (vec4b*)stbi_load_from_memory(
        data, data_size, &width, &height, &ncomp, 4);
    if (!pixels) {
        log_io_error("error loading in-memory image");
        return {};
    }
    auto img = image<vec4b>{width, height, pixels};
    free(pixels);
    return img;
}
image<vec4f> load_stbi_image4f_from_memory(const byte* data, int data_size) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = (vec4f*)stbi_loadf_from_memory(
        data, data_size, &width, &height, &ncomp, 4);
    if (!pixels) {
        log_io_error("error loading in-memory image {}");
        return {};
    }
    auto img = image<vec4f>{width, height, pixels};
    free(pixels);
    return img;
}

// check hdr extensions
bool is_hdr_filename(const string& filename) {
    auto ext = get_extension(filename);
    return ext == "hdr" || ext == "exr" || ext == "pfm";
}

// Loads an hdr image.
image<vec4f> load_image4f_nolog(const string& filename) {
    auto ext = get_extension(filename);
    if (ext == "exr" || ext == "EXR") {
        return load_exr_image4f(filename);
    } else if (ext == "pfm" || ext == "PFM") {
        return load_pfm_image4f(filename);
    } else if (ext == "hdr" || ext == "HDR") {
        return load_stb_image4f(filename);
    } else if (ext == "png" || ext == "PNG") {
        return srgb_to_linear(byte_to_float(load_stb_image4b(filename)));
    } else if (ext == "jpg" || ext == "JPG") {
        return srgb_to_linear(byte_to_float(load_stb_image4b(filename)));
    } else if (ext == "tga" || ext == "TGA") {
        return srgb_to_linear(byte_to_float(load_stb_image4b(filename)));
    } else if (ext == "bmp" || ext == "BMP") {
        return srgb_to_linear(byte_to_float(load_stb_image4b(filename)));
    } else {
        log_io_error("unsupported image format {}", ext);
        return {};
    }
}
image<vec4f> load_image4f(const string& filename) {
    auto scope = log_trace_scoped("loading image {}", filename);
    return load_image4f_nolog(filename);
}

// Saves an hdr image.
bool save_image4f_nolog(const string& filename, const image<vec4f>& img) {
    auto ext = get_extension(filename);
    if (ext == "png" || ext == "PNG") {
        return save_png_image4b(filename, float_to_byte(linear_to_srgb(img)));
    } else if (ext == "jpg" || ext == "JPG") {
        return save_jpg_image4b(filename, float_to_byte(linear_to_srgb(img)));
    } else if (ext == "tga" || ext == "TGA") {
        return save_tga_image4b(filename, float_to_byte(linear_to_srgb(img)));
    } else if (ext == "bmp" || ext == "BMP") {
        return save_bmp_image4b(filename, float_to_byte(linear_to_srgb(img)));
    } else if (ext == "hdr" || ext == "HDR") {
        return save_hdr_image4f(filename, img);
    } else if (ext == "pfm" || ext == "PFM") {
        return save_pfm_image4f(filename, img);
    } else if (ext == "exr" || ext == "EXR") {
        return save_exr_image4f(filename, img);
    } else {
        log_io_error("unsupported image format {}", ext);
        return false;
    }
}
bool save_image4f(const string& filename, const image<vec4f>& img) {
    auto scope = log_trace_scoped("saving image {}", filename);
    return save_image4f_nolog(filename, img);
}

// Loads an hdr image.
image<vec4f> load_image4f_from_memory_nolog(const byte* data, int data_size) {
    return load_stbi_image4f_from_memory(data, data_size);
}
image<vec4f> load_image4f_from_memory(const byte* data, int data_size) {
    auto scope = log_trace_scoped("loading image in memory");
    return load_image4f_from_memory_nolog(data, data_size);
}

// Loads an hdr image.
image<vec4b> load_image4b_nolog(const string& filename) {
    auto ext = get_extension(filename);
    if (ext == "exr" || ext == "EXR") {
        return float_to_byte(linear_to_srgb(load_exr_image4f(filename)));
    } else if (ext == "pfm" || ext == "PFM") {
        return float_to_byte(linear_to_srgb(load_pfm_image4f(filename)));
    } else if (ext == "hdr" || ext == "HDR") {
        return float_to_byte(linear_to_srgb(load_stb_image4f(filename)));
    } else if (ext == "png" || ext == "PNG") {
        return load_stb_image4b(filename);
    } else if (ext == "jpg" || ext == "JPG") {
        return load_stb_image4b(filename);
    } else if (ext == "tga" || ext == "TGA") {
        return load_stb_image4b(filename);
    } else if (ext == "bmp" || ext == "BMP") {
        return load_stb_image4b(filename);
    } else {
        log_io_error("unsupported image format {}", ext);
        return {};
    }
}
image<vec4b> load_image4b(const string& filename) {
    auto scope = log_trace_scoped("loading image {}", filename);
    return load_image4b_nolog(filename);
}

// Saves an ldr image.
bool save_image4b_nolog(const string& filename, const image<vec4b>& img) {
    auto ext = get_extension(filename);
    if (ext == "png" || ext == "PNG") {
        return save_png_image4b(filename, img);
    } else if (ext == "jpg" || ext == "JPG") {
        return save_jpg_image4b(filename, img);
    } else if (ext == "tga" || ext == "TGA") {
        return save_tga_image4b(filename, img);
    } else if (ext == "bmp" || ext == "BMP") {
        return save_bmp_image4b(filename, img);
    } else if (ext == "hdr" || ext == "HDR") {
        return save_hdr_image4f(filename, srgb_to_linear(byte_to_float(img)));
    } else if (ext == "pfm" || ext == "PFM") {
        return save_pfm_image4f(filename, srgb_to_linear(byte_to_float(img)));
    } else if (ext == "exr" || ext == "EXR") {
        return save_exr_image4f(filename, srgb_to_linear(byte_to_float(img)));
    } else {
        log_io_error("unsupported image format {}", ext);
        return false;
    }
}
bool save_image4b(const string& filename, const image<vec4b>& img) {
    auto scope = log_trace_scoped("saving image {}", filename);
    return save_image4b_nolog(filename, img);
}

// Loads an ldr image.
image<vec4b> load_image4b_from_memory_nolog(const byte* data, int data_size) {
    return load_stb_image4b_from_memory(data, data_size);
}
image<vec4b> load_image4b_from_memory(const byte* data, int data_size) {
    auto scope = log_trace_scoped("loading image in memory");
    return load_image4b_from_memory_nolog(data, data_size);
}

// Convenience helper that saves an HDR images as wither a linear HDR file or
// a tonemapped LDR file depending on file name
bool save_tonemapped_image(const string& filename, const image<vec4f>& hdr,
    float exposure, bool filmic, bool srgb) {
    if (is_hdr_filename(filename)) {
        return save_image4f(filename, hdr);
    } else {
        auto ldr = float_to_byte(tonemap_filmic(hdr, exposure, filmic, srgb));
        return save_image4b(filename, ldr);
    }
}

// Resize image.
image<vec4f> resize_image(const image<vec4f>& img, int width, int height) {
    if (!width && !height) {
        log_error("bad image size in resize_image");
        return img;
    }
    if (!width) width = (int)round(img.width * (height / (float)img.height));
    if (!height) height = (int)round(img.height * (width / (float)img.width));
    auto res_img = image<vec4f>{width, height};
    stbir_resize_float_generic((float*)img.pixels.data(), img.width, img.height,
        sizeof(vec4f) * img.width, (float*)res_img.pixels.data(), width, height,
        sizeof(vec4f) * width, 4, 3, 0, STBIR_EDGE_CLAMP, STBIR_FILTER_DEFAULT,
        STBIR_COLORSPACE_LINEAR, nullptr);
    return res_img;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR VOLUME IMAGE IO
// -----------------------------------------------------------------------------
namespace ygl {

// Loads volume data from binary format.
volume<float> load_volume1f_nolog(const string& filename) {
    auto fs = open(filename, "r");
    if (!fs) return {};
    auto size = zero3i;
    if (!read_value(fs, size)) return {};
    auto vol = volume<float>{size.x, size.y, size.z};
    if (!read_values(fs, vol.voxels)) return {};
    return vol;
}
volume<float> load_volume1f(const string& filename) {
    auto scope = log_trace_scoped("loading volume {}", filename);
    return load_volume1f_nolog(filename);
}

// Saves volume data in binary format.
bool save_volume1f_nolog(const string& filename, const volume<float>& vol) {
    auto fs = open(filename, "w");
    if (!fs) return false;
    auto size = vec3i{vol.width, vol.height, vol.depth};
    if (!write_value(fs, size)) return false;
    if (!write_values(fs, vol.voxels)) return false;
    return true;
}
bool save_volume1f(const string& filename, const volume<float>& vol) {
    auto scope = log_trace_scoped("saving volume {}", filename);
    return save_volume1f_nolog(filename, vol);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// GENERIC IMAGE LOADING
// -----------------------------------------------------------------------------
namespace ygl {

// Load a scene
yocto_scene* load_scene(
    const string& filename, bool load_textures, bool skip_missing) {
    auto ext = get_extension(filename);
    if (ext == "json" || ext == "JSON") {
        return load_json_scene(filename, load_textures, skip_missing);
    } else if (ext == "obj" || ext == "OBJ") {
        return load_obj_scene(filename, load_textures, skip_missing);
    } else if (ext == "gltf" || ext == "GLTF") {
        return load_gltf_scene(filename, load_textures, skip_missing);
    } else if (ext == "pbrt" || ext == "PBRT") {
        return load_pbrt_scene(filename, load_textures, skip_missing);
    } else if (ext == "ybin" || ext == "YBIN") {
        return load_ybin_scene(filename, load_textures, skip_missing);
    } else {
        log_io_error("unsupported scene format {}", ext);
        return nullptr;
    }
}

// Save a scene
bool save_scene(const string& filename, const yocto_scene* scene,
    bool save_textures, bool skip_missing) {
    auto ext = get_extension(filename);
    if (ext == "json" || ext == "JSON") {
        return save_json_scene(filename, scene, save_textures, skip_missing);
    } else if (ext == "obj" || ext == "OBJ") {
        return save_obj_scene(filename, scene, save_textures, skip_missing);
    } else if (ext == "gltf" || ext == "GLTF") {
        return save_gltf_scene(filename, scene, save_textures, skip_missing);
    } else if (ext == "pbrt" || ext == "PBRT") {
        return save_pbrt_scene(filename, scene, save_textures, skip_missing);
    } else if (ext == "ybin" || ext == "YBIN") {
        return save_ybin_scene(filename, scene, save_textures, skip_missing);
    } else {
        log_io_error("unsupported scene format {}", ext);
        return false;
    }
}

bool load_scene_textures(yocto_scene* scene, const string& dirname,
    bool skip_missing, bool assign_opacity) {
    // load images
    for (auto texture : scene->textures) {
        if (texture->filename == "" || !texture->hdr_image.pixels.empty() ||
            !texture->ldr_image.pixels.empty())
            continue;
        auto filename = normalize_path(dirname + "/" + texture->filename);
        if (is_hdr_filename(filename)) {
            texture->hdr_image = load_image4f_nolog(filename);
        } else {
            texture->ldr_image = load_image4b_nolog(filename);
        }
        if (texture->hdr_image.pixels.empty() &&
            texture->ldr_image.pixels.empty()) {
            if (!skip_missing) return false;
        }
    }

    // load volumes
    for (auto texture : scene->voltextures) {
        if (texture->filename == "" || !texture->volume_data.voxels.empty())
            continue;
        auto filename = normalize_path(dirname + "/" + texture->filename);
        texture->volume_data = load_volume1f_nolog(filename);
        if (texture->volume_data.voxels.empty()) {
            if (!skip_missing) return false;
        }
    }

    // assign opacity texture if needed
    if (assign_opacity) {
        auto has_opacity = unordered_map<yocto_texture*, bool>();
        for (auto& texture : scene->textures) {
            has_opacity[texture] = false;
            for (auto& p : texture->hdr_image.pixels) {
                if (p.w < 0.999f) {
                    has_opacity[texture] = true;
                    break;
                }
            }
            for (auto& p : texture->ldr_image.pixels) {
                if (p.w < 255) {
                    has_opacity[texture] = true;
                    break;
                }
            }
        }
        for (auto& mat : scene->materials) {
            if (mat->diffuse_texture && !mat->opacity_texture &&
                has_opacity.at(mat->diffuse_texture))
                mat->opacity_texture = mat->diffuse_texture;
        }
    }

    // done
    return true;
}

// helper to save textures
bool save_scene_textures(
    const yocto_scene* scene, const string& dirname, bool skip_missing) {
    // save images
    for (auto texture : scene->textures) {
        if (texture->hdr_image.pixels.empty() && texture->ldr_image.pixels.empty())
            continue;
        auto filename = normalize_path(dirname + "/" + texture->filename);
        if (is_hdr_filename(filename)) {
            if (!save_image4f_nolog(filename, texture->hdr_image)) {
                if (!skip_missing) return false;
            }
        } else {
            if (!save_image4b_nolog(filename, texture->ldr_image)) {
                if (!skip_missing) return false;
            }
        }
    }

    // save volumes
    for (auto texture : scene->voltextures) {
        if (texture->volume_data.voxels.empty()) continue;
        auto filename = normalize_path(dirname + "/" + texture->filename);
        if (!save_volume1f_nolog(filename, texture->volume_data)) {
            if (!skip_missing) return false;
        }
    }

    // done
    return true;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IO UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Encode in base64
string base64_encode(unsigned char const* bytes_to_encode, unsigned int in_len) {
    static const string base64_chars =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz"
        "0123456789+/";

    string        ret;
    int           i = 0;
    int           j = 0;
    unsigned char char_array_3[3];
    unsigned char char_array_4[4];

    while (in_len--) {
        char_array_3[i++] = *(bytes_to_encode++);
        if (i == 3) {
            char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
            char_array_4[1] = ((char_array_3[0] & 0x03) << 4) +
                              ((char_array_3[1] & 0xf0) >> 4);
            char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) +
                              ((char_array_3[2] & 0xc0) >> 6);
            char_array_4[3] = char_array_3[2] & 0x3f;

            for (i = 0; (i < 4); i++) ret += base64_chars[char_array_4[i]];
            i = 0;
        }
    }

    if (i) {
        for (j = i; j < 3; j++) char_array_3[j] = '\0';

        char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
        char_array_4[1] = ((char_array_3[0] & 0x03) << 4) +
                          ((char_array_3[1] & 0xf0) >> 4);
        char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) +
                          ((char_array_3[2] & 0xc0) >> 6);
        char_array_4[3] = char_array_3[2] & 0x3f;

        for (j = 0; (j < i + 1); j++) ret += base64_chars[char_array_4[j]];

        while ((i++ < 3)) ret += '=';
    }

    return ret;
}

// Decode from base64
string base64_decode(string const& encoded_string) {
    static const string base64_chars =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz"
        "0123456789+/";

    auto is_base64 = [](unsigned char c) -> bool {
        return (isalnum(c) || (c == '+') || (c == '/'));
    };

    int           in_len = (int)encoded_string.size();
    int           i      = 0;
    int           j      = 0;
    int           in_    = 0;
    unsigned char char_array_4[4], char_array_3[3];
    string        ret;

    while (in_len-- && (encoded_string[in_] != '=') &&
           is_base64(encoded_string[in_])) {
        char_array_4[i++] = encoded_string[in_];
        in_++;
        if (i == 4) {
            for (i = 0; i < 4; i++)
                char_array_4[i] = base64_chars.find(char_array_4[i]);

            char_array_3[0] = (char_array_4[0] << 2) +
                              ((char_array_4[1] & 0x30) >> 4);
            char_array_3[1] = ((char_array_4[1] & 0xf) << 4) +
                              ((char_array_4[2] & 0x3c) >> 2);
            char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

            for (i = 0; (i < 3); i++) ret += char_array_3[i];
            i = 0;
        }
    }

    if (i) {
        for (j = i; j < 4; j++) char_array_4[j] = 0;

        for (j = 0; j < 4; j++)
            char_array_4[j] = base64_chars.find(char_array_4[j]);

        char_array_3[0] = (char_array_4[0] << 2) +
                          ((char_array_4[1] & 0x30) >> 4);
        char_array_3[1] = ((char_array_4[1] & 0xf) << 4) +
                          ((char_array_4[2] & 0x3c) >> 2);
        char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

        for (j = 0; (j < i - 1); j++) ret += char_array_3[j];
    }

    return ret;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// BUILTIN JSON FORMAT
// -----------------------------------------------------------------------------
namespace ygl {

template <typename T>
bool operator==(const image<T>& a, const image<T>& b) {
    return a.width == b.width && a.height == b.height && a.pixels == b.pixels;
}
template <typename T>
bool operator==(const volume<T>& a, const volume<T>& b) {
    return a._.size == b.extents && a.data == b.data;
}

// Dumps a json value
template <typename T>
bool dump_json_value(json& js, const T& val) {
    try {
        js = val;
        return true;
    } catch (...) {
        return false;
    }
}

// Dumps a json value
template <typename T>
bool dump_json_value(json& js, const T& val, const char* name, const T& def) {
    if (val == def) return true;
    return dump_json_value(js[name], val);
}

// Dumps a json value
template <typename T>
bool parse_json_value(const json& js, T& val) {
    try {
        val = js.get<T>();
        return true;
    } catch (...) {
        return false;
    }
}

// Dumps a json value
template <typename T>
bool parse_json_value(const json& js, T& val, const char* name, const T& def) {
    if (!js.count(name)) return true;
    val = def;
    return parse_json_value(js.at(name), val);
}

// Dumps a json value
template <typename T>
bool dump_json_objref(json& js, T* val, const vector<T*>& refs) {
    return dump_json_value(js, val ? val->name : ""s);
}

// Dumps a json value
template <typename T>
bool dump_json_objref(json& js, T* val, const char* name, const vector<T*>& refs) {
    if (!val) return true;
    return dump_json_objref(js[name], val, refs);
}

// Dumps a json value
template <typename T>
bool dump_json_objref(json& js, const vector<T*>& val, const vector<T*>& refs) {
    js = json::array();
    for (auto v : val) {
        js.push_back({});
        if (!dump_json_objref(js.back(), v, refs)) return false;
    }
    return true;
}

// Dumps a json value
template <typename T>
bool dump_json_objref(
    json& js, const vector<T*>& val, const char* name, const vector<T*>& refs) {
    if (val.empty()) return true;
    return dump_json_objref(js[name], val, refs);
}

// Dumps a json value
template <typename T>
bool parse_json_objref(const json& js, T*& val, const vector<T*>& refs) {
    if (!js.is_string()) return false;
    auto name = ""s;
    if (!parse_json_value(js, name)) return false;
    val = nullptr;
    for (auto ref : refs) {
        if (ref->name == name) {
            val = ref;
            break;
        }
    }
    if (!val) return false;
    return true;
}

// Dumps a json value
template <typename T>
bool parse_json_objref(
    const json& js, T*& val, const char* name, const vector<T*>& refs) {
    if (!js.count(name)) return true;
    val = nullptr;
    return parse_json_objref(js.at(name), val, refs);
}

// Dumps a json value
template <typename T>
bool parse_json_objref(const json& js, vector<T*>& val, const vector<T*>& refs) {
    if (!js.is_array()) return false;
    for (auto& j : js) {
        val.push_back(nullptr);
        if (!parse_json_objref(j, val.back(), refs)) return false;
    }
    return true;
}

// Dumps a json value
template <typename T>
bool parse_json_objref(
    const json& js, vector<T*>& val, const char* name, const vector<T*>& refs) {
    if (!js.count(name)) return true;
    val = {};
    return parse_json_objref(js.at(name), val, refs);
}

// Starts a json object
bool dump_json_objbegin(json& js) {
    js = json::object();
    return true;
}

// Starts a json object
bool parse_json_objbegin(const json& js) { return js.is_object(); }

// Dumps a json value
template <typename T>
bool dump_json_objarray(
    json& js, const vector<T*>& val, const yocto_scene* scene) {
    js = json::array();
    for (auto& v : val) {
        js.push_back({});
        if (!dump_json_object(js.back(), v, scene)) return false;
    }
    return true;
}

// Dumps a json value
template <typename T>
bool dump_json_objarray(json& js, const vector<T*>& val, const char* name,
    const yocto_scene* scene) {
    if (val.empty()) return true;
    return dump_json_objarray(js[name], val, scene);
}

// Dumps a json value
template <typename T>
bool parse_json_objarray(
    const json& js, vector<T*>& val, const yocto_scene* scene) {
    if (!js.is_array()) return false;
    for (auto& j : js) {
        val.push_back(new T());
        if (!parse_json_object(j, val.back(), scene)) return false;
    }
    return true;
}

// Dumps a json value
template <typename T>
bool parse_json_objarray(const json& js, vector<T*>& val, const char* name,
    const yocto_scene* scene) {
    if (!js.count(name)) return true;
    val = {};
    return parse_json_objarray(js.at(name), val, scene);
}

// Parses and applied a JSON procedural
template <typename T>
bool parse_json_procedural(
    const json& js, T* val, const char* name, const yocto_scene* scene) {
    if (!js.count(name)) return true;
    return apply_json_procedural(js.at(name), val, scene);
}

// Procedural commands for cameras
bool apply_json_procedural(
    const json& js, yocto_camera* val, const yocto_scene* scene) {
    if (!parse_json_objbegin(js)) return false;
    if (js.count("from") || js.count("to")) {
        auto from           = js.value("from", zero3f);
        auto to             = js.value("to", zero3f);
        auto up             = js.value("up", vec3f{0, 1, 0});
        val->frame          = lookat_frame(from, to, up);
        val->focus_distance = length(from - to);
    }
    return true;
}

// Serialize struct
bool dump_json_object(
    json& js, const yocto_camera* val, const yocto_scene* scene) {
    static const auto def = yocto_camera();
    if (!dump_json_objbegin(js)) return false;
    if (!dump_json_value(js, val->name, "name", def.name)) return false;
    if (!dump_json_value(js, val->frame, "frame", def.frame)) return false;
    if (!dump_json_value(js, val->orthographic, "orthographic", def.orthographic))
        return false;
    if (!dump_json_value(js, val->film_size, "film_size", def.film_size))
        return false;
    if (!dump_json_value(js, val->focal_length, "focal_length", def.focal_length))
        return false;
    if (!dump_json_value(
            js, val->focus_distance, "focus_distance", def.focus_distance))
        return false;
    if (!dump_json_value(
            js, val->lens_aperture, "lens_aperture", def.lens_aperture))
        return false;
    return true;
}

// Serialize struct
bool parse_json_object(
    const json& js, yocto_camera* val, const yocto_scene* scene) {
    static const auto def = yocto_camera();
    if (!parse_json_objbegin(js)) return false;
    if (!parse_json_value(js, val->name, "name", def.name)) return false;
    if (!parse_json_value(js, val->frame, "frame", def.frame)) return false;
    if (!parse_json_value(js, val->film_size, "film_size", def.film_size))
        return false;
    if (!parse_json_value(js, val->focal_length, "focal_length", def.focal_length))
        return false;
    if (!parse_json_value(
            js, val->focus_distance, "focus_distance", def.focus_distance))
        return false;
    if (!parse_json_value(
            js, val->lens_aperture, "lens_aperture", def.lens_aperture))
        return false;
    if (!parse_json_procedural(js, val, "!!proc", scene)) return false;
    return true;
}

// Serialize struct
bool dump_json_object(
    json& js, const yocto_texture* val, const yocto_scene* scene) {
    static const auto def = yocto_texture();
    if (!dump_json_objbegin(js)) return false;
    if (!dump_json_value(js, val->name, "name", def.name)) return false;
    if (!dump_json_value(js, val->filename, "filename", def.filename))
        return false;
    if (!dump_json_value(
            js, val->clamp_to_edge, "clamp_to_edge", def.clamp_to_edge))
        return false;
    if (!dump_json_value(js, val->height_scale, "height_scale", def.height_scale))
        return false;
    if (!dump_json_value(js, val->no_interpolation, "no_interpolation",
            def.no_interpolation))
        return false;
    if (!dump_json_value(
            js, val->ldr_as_linear, "ldr_as_linear", def.ldr_as_linear))
        return false;
    if (val->filename == "") {
        if (!dump_json_value(js, val->hdr_image, "hdr_image", def.hdr_image))
            return false;
        if (!dump_json_value(js, val->ldr_image, "ldr_image", def.ldr_image))
            return false;
    }
    return true;
}

// Procedural commands for textures
bool apply_json_procedural(
    const json& js, yocto_texture* val, const yocto_scene* scene) {
    if (!parse_json_objbegin(js)) return false;
    auto type = js.value("type", ""s);
    if (type == "") return true;
    auto is_hdr = false;
    auto width  = js.value("width", 512);
    auto height = js.value("height", 512);
    if (type == "grid") {
        val->hdr_image = make_grid_image4f(width, height, js.value("tile", 8),
            js.value("c0", vec4f{0.2f, 0.2f, 0.2f, 1}),
            js.value("c1", vec4f{0.8f, 0.8f, 0.8f, 1}));
    } else if (type == "checker") {
        val->hdr_image = make_checker_image4f(width, height,
            js.value("tile", 8), js.value("c0", vec4f{0.2f, 0.2f, 0.2f, 1}),
            js.value("c1", vec4f{0.8f, 0.8f, 0.8f, 1}));
    } else if (type == "bump") {
        val->hdr_image = make_bumpdimple_image4f(
            width, height, js.value("tile", 8));
    } else if (type == "uvramp") {
        val->hdr_image = make_uvramp_image4f(width, height);
    } else if (type == "uvgrid") {
        val->hdr_image = make_uvgrid_image4f(width, height);
    } else if (type == "sky") {
        if (width < height * 2) width = height * 2;
        val->hdr_image = make_sunsky_image4f(width, height,
            js.value("sun_angle", pif / 4), js.value("turbidity", 3.0f),
            js.value("has_sun", false),
            js.value("ground_albedo", vec3f{0.7f, 0.7f, 0.7f}));
        is_hdr         = true;
    } else if (type == "noise") {
        val->hdr_image = make_noise_image4f(
            width, height, js.value("scale", 1.0f), js.value("wrap", true));
    } else if (type == "fbm") {
        val->hdr_image = make_fbm_image4f(width, height, js.value("scale", 1.0f),
            js.value("lacunarity", 2.0f), js.value("gain", 0.5f),
            js.value("octaves", 6), js.value("wrap", true));
    } else if (type == "ridge") {
        val->hdr_image = make_ridge_image4f(width, height,
            js.value("scale", 1.0f), js.value("lacunarity", 2.0f),
            js.value("gain", 0.5f), js.value("offset", 1.0f),
            js.value("octaves", 6), js.value("wrap", true));
    } else if (type == "turbulence") {
        val->hdr_image = make_turbulence_image4f(width, height,
            js.value("scale", 1.0f), js.value("lacunarity", 2.0f),
            js.value("gain", 0.5f), js.value("octaves", 6),
            js.value("wrap", true));
    } else {
        log_error("unknown texture type {}", type);
        return false;
    }
    if (js.value("bump_to_normal", false)) {
        val->hdr_image = bump_to_normal_map(
            val->hdr_image, js.value("bump_scale", 1.0f));
        val->ldr_as_linear = true;
    }
    if (!is_hdr) {
        if (!val->ldr_as_linear) {
            val->ldr_image = float_to_byte(linear_to_srgb(val->hdr_image));
        } else {
            val->ldr_image = float_to_byte(val->hdr_image);
        }
        val->hdr_image = {};
    }
    if (val->filename == "") {
        auto ext      = (is_hdr) ? string("hdr") : string("png");
        val->filename = "textures/" + val->name + "." + ext;
    }
    return true;
}

// Serialize struct
bool parse_json_object(
    const json& js, yocto_texture* val, const yocto_scene* scene) {
    static const auto def = yocto_texture();
    if (!parse_json_objbegin(js)) return false;
    if (!parse_json_value(js, val->name, "name", def.name)) return false;
    if (!parse_json_value(js, val->filename, "filename", def.filename))
        return false;
    if (!parse_json_value(
            js, val->clamp_to_edge, "clamp_to_edge", def.clamp_to_edge))
        return false;
    if (!parse_json_value(js, val->height_scale, "height_scale", def.height_scale))
        return false;
    if (!parse_json_value(js, val->no_interpolation, "no_interpolation",
            def.no_interpolation))
        return false;
    if (!parse_json_value(
            js, val->ldr_as_linear, "ldr_as_linear", def.ldr_as_linear))
        return false;
    if (!parse_json_value(js, val->hdr_image, "hdr_image", def.hdr_image))
        return false;
    if (!parse_json_value(js, val->ldr_image, "ldr_image", def.ldr_image))
        return false;
    if (!parse_json_procedural(js, val, "!!proc", scene)) return false;
    return true;
}

// Serialize struct
bool dump_json_object(
    json& js, const yocto_voltexture* val, const yocto_scene* scene) {
    static const auto def = yocto_voltexture();
    if (!dump_json_objbegin(js)) return false;
    if (!dump_json_value(js, val->name, "name", def.name)) return false;
    if (!dump_json_value(js, val->filename, "filename", def.filename))
        return false;
    if (!dump_json_value(
            js, val->clamp_to_edge, "clamp_to_edge", def.clamp_to_edge))
        return false;
    if (!dump_json_value(js, val->no_interpolation, "no_interpolation",
            def.no_interpolation))
        return false;
    if (val->filename == "") {
        if (!val->volume_data.voxels.empty()) js["vol"] = val->volume_data;
    }
    return true;
}

// Procedural commands for textures
bool apply_json_procedural(
    const json& js, yocto_voltexture* val, const yocto_scene* scene) {
    if (!parse_json_objbegin(js)) return false;
    auto type = js.value("type", ""s);
    if (type == "") return true;
    auto width  = js.value("width", 512);
    auto height = js.value("height", 512);
    auto depth  = js.value("depth", 512);
    if (type == "test_volume") {
        val->volume_data = make_test_volume1f(width, height, depth,
            js.value("scale", 10.0f), js.value("exponent", 6.0f));
    } else {
        log_error("unknown texture type {}", type);
        return false;
    }
    if (val->filename == "") {
        auto ext      = string("vol");
        val->filename = "textures/" + val->name + "." + ext;
    }
    return true;
}

// Serialize struct
bool parse_json_object(
    const json& js, yocto_voltexture* val, const yocto_scene* scene) {
    static const auto def = yocto_voltexture();
    if (!parse_json_objbegin(js)) return false;
    if (!parse_json_value(js, val->name, "name", def.name)) return false;
    if (!parse_json_value(js, val->filename, "filename", def.filename))
        return false;
    if (!parse_json_value(
            js, val->clamp_to_edge, "clamp_to_edge", def.clamp_to_edge))
        return false;
    if (!parse_json_value(js, val->no_interpolation, "no_interpolation",
            def.no_interpolation))
        return false;
    if (!parse_json_value(js, val->volume_data, "vol", def.volume_data))
        return false;
    if (!parse_json_procedural(js, val, "!!proc", scene)) return false;
    return true;
}

// Serialize struct
bool dump_json_object(
    json& js, const yocto_material* val, const yocto_scene* scene) {
    static const auto def = yocto_material();
    if (!dump_json_objbegin(js)) return false;
    if (!dump_json_value(js, val->name, "name", def.name)) return false;
    if (val->base_metallic != def.base_metallic)
        js["base_metallic"] = val->base_metallic;
    if (val->gltf_textures != def.gltf_textures)
        js["gltf_textures"] = val->gltf_textures;
    if (val->double_sided != def.double_sided)
        js["double_sided"] = val->double_sided;
    if (!dump_json_value(js, val->emission, "emission", def.emission))
        return false;
    if (!dump_json_value(js, val->diffuse, "diffuse", def.diffuse))
        return false;
    if (!dump_json_value(js, val->specular, "specular", def.specular))
        return false;
    if (!dump_json_value(js, val->transmission, "transmission", def.transmission))
        return false;
    if (!dump_json_value(js, val->roughness, "roughness", def.roughness))
        return false;
    if (!dump_json_value(js, val->opacity, "opacity", def.opacity))
        return false;
    if (!dump_json_value(js, val->fresnel, "fresnel", def.fresnel))
        return false;
    if (!dump_json_value(js, val->refract, "refract", def.refract))
        return false;
    if (!dump_json_objref(
            js, val->emission_texture, "emission_texture", scene->textures))
        return false;
    if (!dump_json_objref(
            js, val->diffuse_texture, "diffuse_texture", scene->textures))
        return false;
    if (!dump_json_objref(
            js, val->specular_texture, "specular_texture", scene->textures))
        return false;
    if (!dump_json_objref(js, val->transmission_texture, "transmission_texture",
            scene->textures))
        return false;
    if (!dump_json_objref(
            js, val->roughness_texture, "roughness_texture", scene->textures))
        return false;
    if (!dump_json_objref(
            js, val->opacity_texture, "opacity_texture", scene->textures))
        return false;
    if (!dump_json_objref(
            js, val->occlusion_texture, "occlusion_texture", scene->textures))
        return false;
    if (!dump_json_objref(js, val->bump_texture, "bump_texture", scene->textures))
        return false;
    if (!dump_json_objref(js, val->displacement_texture, "displacement_texture",
            scene->textures))
        return false;
    if (!dump_json_objref(
            js, val->normal_texture, "normal_texture", scene->textures))
        return false;
    if (!dump_json_objref(js, val->volume_density_texture,
            "volume_density_texture", scene->voltextures))
        return false;
    return true;
}

// Procedural commands for materials
bool apply_json_procedural(
    const json& js, yocto_material* val, const yocto_scene* scene) {
    if (!parse_json_objbegin(js)) return false;
    return true;
}

// Serialize struct
bool parse_json_object(
    const json& js, yocto_material* val, const yocto_scene* scene) {
    static const auto def = yocto_material();
    if (!parse_json_objbegin(js)) return false;
    if (!parse_json_value(js, val->name, "name", def.name)) return false;
    if (!parse_json_value(
            js, val->base_metallic, "base_metallic", def.base_metallic))
        return false;
    if (!parse_json_value(
            js, val->gltf_textures, "gltf_textures", def.gltf_textures))
        return false;
    if (!parse_json_value(js, val->double_sided, "double_sided", def.double_sided))
        return false;
    if (!parse_json_value(js, val->emission, "emission", def.emission))
        return false;
    if (!parse_json_value(js, val->diffuse, "diffuse", def.diffuse))
        return false;
    if (!parse_json_value(js, val->specular, "specular", def.specular))
        return false;
    if (!parse_json_value(js, val->transmission, "transmission", def.transmission))
        return false;
    if (!parse_json_value(js, val->roughness, "roughness", def.roughness))
        return false;
    if (!parse_json_value(js, val->opacity, "opacity", def.opacity))
        return false;
    if (!parse_json_value(
            js, val->volume_emission, "volume_emission", def.volume_emission))
        return false;
    if (!parse_json_value(
            js, val->volume_albedo, "volume_albedo", def.volume_albedo))
        return false;
    if (!parse_json_value(
            js, val->volume_density, "volume_density", def.volume_density))
        return false;
    if (!parse_json_value(
            js, val->volume_phaseg, "volume_phaseg", def.volume_phaseg))
        return false;
    if (!parse_json_value(js, val->fresnel, "fresnel", def.fresnel))
        return false;
    if (!parse_json_value(js, val->refract, "refract", def.refract))
        return false;
    if (!parse_json_objref(
            js, val->emission_texture, "emission_texture", scene->textures))
        return false;
    if (!parse_json_objref(
            js, val->diffuse_texture, "diffuse_texture", scene->textures))
        return false;
    if (!parse_json_objref(
            js, val->specular_texture, "specular_texture", scene->textures))
        return false;
    if (!parse_json_objref(js, val->transmission_texture,
            "transmission_texture", scene->textures))
        return false;
    if (!parse_json_objref(
            js, val->roughness_texture, "roughness_texture", scene->textures))
        return false;
    if (!parse_json_objref(
            js, val->opacity_texture, "opacity_texture", scene->textures))
        return false;
    if (!parse_json_objref(
            js, val->occlusion_texture, "occlusion_texture", scene->textures))
        return false;
    if (!parse_json_objref(js, val->bump_texture, "bump_texture", scene->textures))
        return false;
    if (!parse_json_objref(js, val->displacement_texture,
            "displacement_texture", scene->textures))
        return false;
    if (!parse_json_objref(
            js, val->normal_texture, "normal_texture", scene->textures))
        return false;
    if (!parse_json_objref(js, val->volume_density_texture,
            "volume_density_texture", scene->voltextures))
        return false;
    if (!parse_json_procedural(js, val, "!!proc", scene)) return false;
    return true;
}

// Serialize struct
bool dump_json_object(json& js, const yocto_shape* val, const yocto_scene* scene) {
    static const auto def = yocto_shape();
    if (!dump_json_objbegin(js)) return false;
    if (!dump_json_value(js, val->name, "name", def.name)) return false;
    if (!dump_json_value(js, val->filename, "filename", def.filename))
        return false;
    if (val->filename == "") {
        if (!dump_json_value(js, val->points, "points", def.points))
            return false;
        if (!dump_json_value(js, val->lines, "lines", def.lines)) return false;
        if (!dump_json_value(js, val->triangles, "triangles", def.triangles))
            return false;
        if (!dump_json_value(js, val->quads, "quads", def.quads)) return false;
        if (!dump_json_value(js, val->positions, "positions", def.positions))
            return false;
        if (!dump_json_value(js, val->normals, "normals", def.normals))
            return false;
        if (!dump_json_value(
                js, val->texturecoords, "texturecoords", def.texturecoords))
            return false;
        if (!dump_json_value(js, val->colors, "colors", def.colors))
            return false;
        if (!dump_json_value(js, val->radius, "radius", def.radius))
            return false;
        if (!dump_json_value(
                js, val->tangentspaces, "tangent_spaces", def.tangentspaces))
            return false;
    }
    return true;
}

// Procedural commands for materials
bool apply_json_procedural(
    const json& js, yocto_shape* val, const yocto_scene* scene) {
    if (!parse_json_objbegin(js)) return false;
    auto type = js.value("type", ""s);
    if (type == "") return true;
    auto shape        = make_shape_data();
    auto as_triangles = js.value("as_triangles", false);
    if (type == "quad") {
        shape = make_quad_shape(js.value("steps", vec2i{1, 1}),
            js.value("size", vec2f{2, 2}), js.value("uvsize", vec2f{1, 1}),
            as_triangles);
    } else if (type == "quady") {
        shape = make_quad_shape(js.value("steps", vec2i{1, 1}),
            js.value("size", vec2f{2, 2}), js.value("uvsize", vec2f{1, 1}),
            as_triangles);
    } else if (type == "quad_stack") {
        shape = make_quad_stack_shape(js.value("steps", vec3i{1, 1, 1}),
            js.value("size", vec3f{2, 2, 2}), js.value("uvsize", vec2f{1, 1}),
            as_triangles);
    } else if (type == "cube") {
        shape = make_cube_shape(js.value("steps", vec3i{1, 1, 1}),
            js.value("size", vec3f{2, 2, 2}),
            js.value("uvsize", vec3f{1, 1, 1}), as_triangles);
    } else if (type == "cube_rounded") {
        shape = make_cube_rounded_shape(js.value("steps", vec3i{32, 32, 32}),
            js.value("size", vec3f{2, 2, 2}), js.value("uvsize", vec3f{1, 1, 1}),
            js.value("radius", 0.3f), as_triangles);
    } else if (type == "sphere") {
        shape = make_sphere_shape(js.value("steps", vec2i{64, 32}),
            js.value("size", 2.0f), js.value("uvsize", vec2f{1, 1}),
            as_triangles);
    } else if (type == "sphere_cube") {
        shape = make_sphere_cube_shape(js.value("steps", 32),
            js.value("size", 2.0f), js.value("uvsize", 1.0f), as_triangles);
    } else if (type == "sphere_flipcap") {
        shape = make_sphere_flipcap_shape(js.value("steps", vec2i{64, 32}),
            js.value("size", 2.0f), js.value("uvsize", vec2f{1, 1}),
            js.value("zflip", vec2f{-0.75f, +0.75f}), as_triangles);
    } else if (type == "disk") {
        shape = make_disk_shape(js.value("steps", vec2i{32, 16}),
            js.value("size", 2.0f), js.value("uvsize", vec2f{1, 1}),
            as_triangles);
    } else if (type == "disk_quad") {
        shape = make_disk_quad_shape(js.value("steps", 32),
            js.value("size", 2.0f), js.value("uvsize", 1.0f), as_triangles);
    } else if (type == "disk_bulged") {
        shape = make_disk_bulged_shape(js.value("steps", 32),
            js.value("size", 2.0f), js.value("uvsize", 1.0f),
            js.value("height", 0.25f), as_triangles);
    } else if (type == "cylinder_side") {
        shape = make_cylinder_side_shape(js.value("steps", vec2i{64, 32}),
            js.value("size", vec2f{2.0f, 2.0f}),
            js.value("uvsize", vec2f{1, 1}), as_triangles);
    } else if (type == "cylinder") {
        shape = make_cylinder_shape(js.value("steps", vec3i{64, 32, 16}),
            js.value("size", vec2f{2.0f, 2.0f}),
            js.value("uvsize", vec3f{1, 1, 1}), as_triangles);
    } else if (type == "cylinder_rounded") {
        shape = make_cylinder_rounded_shape(js.value("steps", vec3i{64, 32, 16}),
            js.value("size", vec2f{2.0f, 2.0f}),
            js.value("uvsize", vec3f{1, 1, 1}), js.value("radius", 0.15f),
            as_triangles);
    } else if (type == "sphere_geodesic") {
        shape = make_geodesic_sphere_shape(
            js.value("tesselation", 4), js.value("size", 2.0f), as_triangles);
    } else if (type == "floor") {
        shape = make_floor_shape(js.value("steps", vec2i{1, 1}),
            js.value("size", vec2f{40, 40}), js.value("uvsize", vec2f{20, 20}),
            as_triangles);
    } else if (type == "matball") {
        shape = make_sphere_shape(js.value("steps", vec2i{64, 32}),
            js.value("size", 2.0f), js.value("uvsize", vec2f{1, 1}),
            as_triangles);
    } else if (type == "hairball") {
        auto base = make_sphere_cube_shape(
            32, js.value("size", 2.0f) * 0.8f, 1, as_triangles);
        shape = make_hair_shape(js.value("steps", vec2i{4, 65536}),
            base.triangles, base.quads, base.positions, base.normals,
            base.texturecoords, js.value("length", vec2f{0.2f, 0.2f}),
            js.value("radius", vec2f{0.001f, 0.001f}),
            js.value("noise", vec2f{0, 0}), js.value("clump", vec2f{0, 0}));
    } else if (type == "hairball_interior") {
        shape = make_sphere_cube_shape(
            32, js.value("size", 2.0f) * 0.8f, 1, as_triangles);
    } else if (type == "suzanne") {
        shape = make_suzanne_shape(js.value("size", 2.0f), as_triangles);
    } else {
        log_error("unknown shape type {}", type);
        return false;
    }
    if (js.value("flipyz", false)) {
        for (auto& p : shape.positions) p = {p.x, p.z, p.y};
        for (auto& n : shape.normals) n = {n.x, n.z, n.y};
    }
    val->points        = shape.points;
    val->lines         = shape.lines;
    val->triangles     = shape.triangles;
    val->quads         = shape.quads;
    val->positions     = shape.positions;
    val->normals       = shape.normals;
    val->texturecoords = shape.texturecoords;
    val->radius        = shape.radius;
    return true;
}

// Serialize struct
bool parse_json_object(
    const json& js, yocto_shape* val, const yocto_scene* scene) {
    static const auto def = yocto_shape();
    if (!parse_json_objbegin(js)) return false;
    if (!parse_json_value(js, val->name, "name", def.name)) return false;
    if (!parse_json_value(js, val->filename, "filename", def.filename))
        return false;
    if (!parse_json_value(js, val->points, "points", def.points)) return false;
    if (!parse_json_value(js, val->lines, "lines", def.lines)) return false;
    if (!parse_json_value(js, val->triangles, "triangles", def.triangles))
        return false;
    if (!parse_json_value(js, val->quads, "quads", def.quads)) return false;
    if (!parse_json_value(js, val->positions, "positions", def.positions))
        return false;
    if (!parse_json_value(js, val->normals, "normals", def.normals))
        return false;
    if (!parse_json_value(
            js, val->texturecoords, "texturecoords", def.texturecoords))
        return false;
    if (!parse_json_value(js, val->colors, "color", def.colors)) return false;
    if (!parse_json_value(js, val->radius, "radius", def.radius)) return false;
    if (!parse_json_value(
            js, val->tangentspaces, "tangent_spaces", def.tangentspaces))
        return false;
    if (!parse_json_procedural(js, val, "!!proc", scene)) return false;
    return true;
}

// Serialize struct
bool dump_json_object(
    json& js, const yocto_surface* val, const yocto_scene* scene) {
    static const auto def = yocto_surface();
    if (!dump_json_objbegin(js)) return false;
    if (!dump_json_value(js, val->name, "name", def.name)) return false;
    if (!dump_json_value(js, val->filename, "filename", def.filename))
        return false;
    if (!dump_json_value(
            js, val->subdivision_level, "level", def.subdivision_level))
        return false;
    if (val->catmull_clark != def.catmull_clark)
        js["catmull_clark"] = val->catmull_clark;
    if (val->compute_vertex_normals != def.compute_vertex_normals)
        js["compute_vertex_normals"] = val->compute_vertex_normals;
    if (val->filename == "") {
        if (!dump_json_value(js, val->positions_quads, "positions_quads",
                def.positions_quads))
            return false;
        if (!dump_json_value(js, val->texturecoords_quads,
                "texturecoords_quads", def.texturecoords_quads))
            return false;
        if (!dump_json_value(
                js, val->colors_quads, "colors_quads", def.colors_quads))
            return false;
        if (!dump_json_value(js, val->positions, "positions", def.positions))
            return false;
        if (!dump_json_value(
                js, val->texturecoords, "texturecoords", def.texturecoords))
            return false;
        if (!dump_json_value(js, val->colors, "colors", def.colors))
            return false;
    }
    return true;
}

// Procedural commands for subdivs
bool apply_json_procedural(
    const json& js, yocto_surface* val, const yocto_scene* scene) {
    if (!parse_json_objbegin(js)) return false;
    auto type = js.value("type", ""s);
    if (type == "") return true;
    auto shape = make_fvshape_data();
    if (type == "cube") {
        shape = make_cube_fvshape(js.value("steps", vec3i{1, 1, 1}),
            js.value("size", vec3f{2, 2, 2}), js.value("uvsize", vec3f{1, 1, 1}));
    } else if (type == "cube_open") {
        shape = make_cube_fvshape(js.value("steps", vec3i{1, 1, 1}),
            js.value("size", vec3f{2, 2, 2}), js.value("uvsize", vec3f{1, 1, 1}));
        shape.positions_quads.pop_back();
        shape.normals_quads.pop_back();
        shape.quads_texturecoords.pop_back();
    } else if (type == "suzanne") {
        auto qshp = make_suzanne_shape(js.value("size", 2.0f), false);
        shape.positions_quads = qshp.quads;
        shape.positions       = qshp.positions;
    } else {
        log_error("unknown shape type {}", type);
        return false;
    }
    val->positions_quads     = shape.positions_quads;
    val->positions           = shape.positions;
    val->texturecoords_quads = shape.quads_texturecoords;
    val->texturecoords       = shape.texturecoords;
    if (val->filename == "") val->filename = "meshes/" + val->name + ".obj";
    return true;
}

// Serialize struct
bool parse_json_object(
    const json& js, yocto_surface* val, const yocto_scene* scene) {
    static const auto def = yocto_surface();
    if (!parse_json_objbegin(js)) return false;
    if (!parse_json_value(js, val->name, "name", def.name)) return false;
    if (!parse_json_value(js, val->filename, "filename", def.filename))
        return false;
    if (!parse_json_value(js, val->subdivision_level, "subdivision_level",
            def.subdivision_level))
        return false;
    if (!parse_json_value(
            js, val->catmull_clark, "catmull_clark", def.catmull_clark))
        return false;
    if (!parse_json_value(js, val->compute_vertex_normals,
            "compute_vertex_normals", def.compute_vertex_normals))
        return false;
    if (!parse_json_value(
            js, val->positions_quads, "positions_quads", def.positions_quads))
        return false;
    if (!parse_json_value(js, val->texturecoords_quads, "texturecoords_quads",
            def.texturecoords_quads))
        return false;
    if (!parse_json_value(js, val->colors_quads, "colors_quads", def.colors_quads))
        return false;
    if (!parse_json_value(js, val->positions, "positions", def.positions))
        return false;
    if (!parse_json_value(
            js, val->texturecoords, "texturecoords", def.texturecoords))
        return false;
    if (!parse_json_value(js, val->colors, "colors", def.colors)) return false;
    if (!parse_json_procedural(js, val, "!!proc", scene)) return false;
    return true;
}

// Serialize struct
bool dump_json_object(
    json& js, const yocto_instance* val, const yocto_scene* scene) {
    static const auto def = yocto_instance();
    if (!dump_json_objbegin(js)) return false;
    if (!dump_json_value(js, val->name, "name", def.name)) return false;
    if (!dump_json_value(js, val->frame, "frame", def.frame)) return false;
    if (!dump_json_objref(js, val->shape, "shape", scene->shapes)) return false;
    if (!dump_json_objref(js, val->material, "material", scene->materials))
        return false;
    if (!dump_json_objref(js, val->surface, "surface", scene->surfaces))
        return false;
    return true;
}

// Procedural commands for instances
bool apply_json_procedural(
    const json& js, yocto_instance* val, const yocto_scene* scene) {
    if (!parse_json_objbegin(js)) return false;
    if (js.count("from")) {
        auto from  = js.value("from", zero3f);
        auto to    = js.value("to", zero3f);
        auto up    = js.value("up", vec3f{0, 1, 0});
        val->frame = lookat_frame(from, to, up, true);
    }
    if (js.count("translation") || js.count("rotation") || js.count("scale")) {
        auto translation = js.value("translation", zero3f);
        auto rotation    = js.value("rotation", zero4f);
        auto scaling     = js.value("scale", vec3f{1, 1, 1});
        val->frame = translation_frame(translation) * scaling_frame(scaling) *
                     rotation_frame(xyz(rotation), rotation.w);
    }
    return true;
}

// Serialize struct
bool parse_json_object(
    const json& js, yocto_instance* val, const yocto_scene* scene) {
    static const auto def = yocto_instance();
    if (!parse_json_objbegin(js)) return false;
    if (!parse_json_value(js, val->name, "name", def.name)) return false;
    if (!parse_json_value(js, val->frame, "frame", def.frame)) return false;
    if (!parse_json_objref(js, val->shape, "shape", scene->shapes))
        return false;
    if (!parse_json_objref(js, val->surface, "surface", scene->surfaces))
        return false;
    if (!parse_json_objref(js, val->material, "material", scene->materials))
        return false;
    if (!parse_json_procedural(js, val, "!!proc", scene)) return false;
    return true;
}

// Serialize struct
bool dump_json_object(
    json& js, const yocto_environment* val, const yocto_scene* scene) {
    static const auto def = yocto_environment();
    if (!dump_json_objbegin(js)) return false;
    if (!dump_json_value(js, val->name, "name", def.name)) return false;
    if (!dump_json_value(js, val->frame, "frame", def.frame)) return false;
    if (!dump_json_value(js, val->emission, "emission", def.emission))
        return false;
    if (!dump_json_objref(
            js, val->emission_texture, "emission_texture", scene->textures))
        return false;
    return true;
}

// Procedural commands for materials
bool apply_json_procedural(
    const json& js, yocto_environment* val, const yocto_scene* scene) {
    if (!parse_json_objbegin(js)) return false;
    if (js.count("rotation")) {
        auto rotation = js.value("rotation", zero4f);
        val->frame    = rotation_frame(xyz(rotation), rotation.w);
    }
    return true;
}

// Serialize struct
bool parse_json_object(
    const json& js, yocto_environment* val, const yocto_scene* scene) {
    static const auto def = yocto_environment();
    if (!parse_json_objbegin(js)) return false;
    if (!parse_json_value(js, val->name, "name", def.name)) return false;
    if (!parse_json_value(js, val->frame, "frame", def.frame)) return false;
    if (!parse_json_value(js, val->emission, "emission", def.emission))
        return false;
    if (!parse_json_objref(
            js, val->emission_texture, "emission_texture", scene->textures))
        return false;
    if (!parse_json_procedural(js, val, "!!proc", scene)) return false;
    return true;
}

// Serialize struct
bool dump_json_object(
    json& js, const yocto_scene_node* val, const yocto_scene* scene) {
    static const auto def = yocto_scene_node();
    if (!dump_json_objbegin(js)) return false;
    if (!dump_json_value(js, val->name, "name", def.name)) return false;
    if (!dump_json_value(js, val->local, "local", def.local)) return false;
    if (!dump_json_value(js, val->translation, "translation", def.translation))
        return false;
    if (!dump_json_value(js, val->rotation, "rotation", def.rotation))
        return false;
    if (!dump_json_value(js, val->scale, "scale", def.scale)) return false;
    if (!dump_json_value(js, val->weights, "weights", def.weights))
        return false;
    if (!dump_json_objref(js, val->parent, "parent", scene->nodes))
        return false;
    if (!dump_json_objref(js, val->camera, "camera", scene->cameras))
        return false;
    if (!dump_json_objref(js, val->instance, "instance", scene->instances))
        return false;
    if (!dump_json_objref(
            js, val->environment, "environment", scene->environments))
        return false;
    return true;
}

// Procedural commands for nodes
bool apply_json_procedural(
    const json& js, yocto_scene_node* val, const yocto_scene* scene) {
    if (!parse_json_objbegin(js)) return false;
    if (js.count("from")) {
        auto from  = js.value("from", zero3f);
        auto to    = js.value("to", zero3f);
        auto up    = js.value("up", vec3f{0, 1, 0});
        val->local = lookat_frame(from, to, up, true);
    }
    return true;
}

// Serialize struct
bool parse_json_object(
    const json& js, yocto_scene_node* val, const yocto_scene* scene) {
    static const auto def = yocto_scene_node();
    if (!parse_json_objbegin(js)) return false;
    if (!parse_json_value(js, val->name, "name", def.name)) return false;
    if (!parse_json_value(js, val->local, "local", def.local)) return false;
    if (!parse_json_value(js, val->translation, "translation", def.translation))
        return false;
    if (!parse_json_value(js, val->rotation, "rotation", def.rotation))
        return false;
    if (!parse_json_value(js, val->scale, "scale", def.scale)) return false;
    if (!parse_json_value(js, val->weights, "weights", def.weights))
        return false;
    if (!parse_json_objref(js, val->parent, "parent", scene->nodes))
        return false;
    if (!parse_json_objref(js, val->instance, "instance", scene->instances))
        return false;
    if (!parse_json_objref(js, val->camera, "camera", scene->cameras))
        return false;
    if (!parse_json_objref(
            js, val->environment, "environment", scene->environments))
        return false;
    if (!parse_json_procedural(js, val, "!!proc", scene)) return false;
    return true;
}

// Serialize enum
bool dump_json_value(json& js, const yocto_interpolation_type& val) {
    static auto names = map<yocto_interpolation_type, string>{
        {yocto_interpolation_type::linear, "linear"},
        {yocto_interpolation_type::step, "step"},
        {yocto_interpolation_type::bezier, "bezier"},
    };
    auto vals = names.at(val);
    return dump_json_value(js, vals);
}

// Serialize enum
bool parse_json_value(const json& js, yocto_interpolation_type& val) {
    static auto names = map<string, yocto_interpolation_type>{
        {"linear", yocto_interpolation_type::linear},
        {"step", yocto_interpolation_type::step},
        {"bezier", yocto_interpolation_type::bezier},
    };
    auto vals = ""s;
    if (!parse_json_value(js, vals)) return false;
    try {
        val = names.at(js.get<string>());
    } catch (...) {
        return false;
    }
    return true;
}

// Serialize struct
bool dump_json_object(
    json& js, const yocto_animation* val, const yocto_scene* scene) {
    static const auto def = yocto_animation();
    if (!dump_json_objbegin(js)) return false;
    if (!dump_json_value(js, val->name, "name", def.name)) return false;
    if (!dump_json_value(js, val->filename, "filename", def.filename))
        return false;
    if (!dump_json_value(
            js, val->animation_group, "animation_group", def.animation_group))
        return false;
    if (!dump_json_value(
            js, val->interpolation_type, "type", def.interpolation_type))
        return false;
    if (val->filename == "") {
        if (!dump_json_value(js, val->keyframes_times, "keyframes_times",
                def.keyframes_times))
            return false;
        if (!dump_json_value(js, val->translation_keyframes,
                "translation_keyframes", def.translation_keyframes))
            return false;
        if (!dump_json_value(js, val->rotation_keyframes, "rotation_keyframes",
                def.rotation_keyframes))
            return false;
        if (!dump_json_value(js, val->scale_keyframes, "scale_keyframes",
                def.scale_keyframes))
            return false;
    }
    if (!dump_json_objref(js, val->node_targets, "node_targets", scene->nodes))
        return false;
    return true;
}

// Procedural commands for animations
bool apply_json_procedural(
    const json& js, yocto_animation* val, const yocto_scene* scene) {
    if (!parse_json_objbegin(js)) return false;
    if (js.count("rotation_axisangle")) {
        for (auto& j : js.at("rotation_axisangle")) {
            val->rotation_keyframes.push_back(rotation_quat(j.get<vec4f>()));
        }
    }
    return true;
}

// Serialize struct
bool parse_json_object(
    const json& js, yocto_animation* val, const yocto_scene* scene) {
    static const auto def = yocto_animation();
    if (!parse_json_objbegin(js)) return false;
    if (!parse_json_value(js, val->name, "name", def.name)) return false;
    if (!parse_json_value(js, val->filename, "filename", def.filename))
        return false;
    if (!parse_json_value(
            js, val->animation_group, "animation_group", def.animation_group))
        return false;
    if (!parse_json_value(js, val->interpolation_type, "interpolation_type",
            def.interpolation_type))
        return false;
    if (!parse_json_value(
            js, val->keyframes_times, "keyframes_times", def.keyframes_times))
        return false;
    if (!parse_json_value(js, val->translation_keyframes,
            "translation_keyframes", def.translation_keyframes))
        return false;
    if (!parse_json_value(js, val->rotation_keyframes, "rotation_keyframes",
            def.rotation_keyframes))
        return false;
    if (!parse_json_value(
            js, val->scale_keyframes, "scale_keyframes", def.scale_keyframes))
        return false;
    if (!parse_json_objref(js, val->node_targets, "node_targets", scene->nodes))
        return false;
    if (!parse_json_procedural(js, val, "!!proc", scene)) return false;
    return true;
}

// Serialize struct
bool dump_json_object(json& js, const yocto_scene* val, const yocto_scene* scene) {
    static const auto def = yocto_scene();
    if (!dump_json_objbegin(js)) return false;
    if (!dump_json_value(js, val->name, "name", def.name)) return false;
    if (!dump_json_objarray(js, val->cameras, "cameras", scene)) return false;
    if (!dump_json_objarray(js, val->textures, "textures", scene)) return false;
    if (!dump_json_objarray(js, val->materials, "materials", scene))
        return false;
    if (!dump_json_objarray(js, val->shapes, "shapes", scene)) return false;
    if (!dump_json_objarray(js, val->surfaces, "surfaces", scene)) return false;
    if (!dump_json_objarray(js, val->instances, "instances", scene))
        return false;
    if (!dump_json_objarray(js, val->environments, "environments", scene))
        return false;
    if (!dump_json_objarray(js, val->nodes, "nodes", scene)) return false;
    if (!dump_json_objarray(js, val->animations, "animations", scene))
        return false;
    return true;
}

// Procedural commands for scenes
bool apply_json_procedural(
    const json& js, yocto_scene* val, const yocto_scene* scene) {
    if (!parse_json_objbegin(js)) return false;
    if (js.count("random_instances")) {
        auto& jjs  = js.at("random_instances");
        auto  num  = jjs.value("num", 100);
        auto  seed = jjs.value("seed", 13);
        auto  base = make_unique<yocto_instance>();
        parse_json_object(jjs.at("base"), base.get(), scene);
        auto ists = vector<unique_ptr<yocto_instance>>();
        for (auto& j : jjs.at("instances")) {
            ists.push_back(make_unique<yocto_instance>());
            parse_json_object(j, ists.back().get(), scene);
        }

        auto pos                 = vector<vec3f>();
        auto norm                = vector<vec3f>();
        auto texcoord            = vector<vec2f>();
        tie(pos, norm, texcoord) = sample_triangles_points(
            base->shape->triangles, base->shape->positions,
            base->shape->normals, base->shape->texturecoords, num, seed);

        auto nmap = unordered_map<yocto_instance*, int>();
        for (auto& instance : ists) nmap[instance.get()] = 0;
        auto rng = make_rng(seed, 17);
        for (auto i = 0; i < num; i++) {
            auto instance = ists.at(get_random_int(rng, (int)ists.size() - 1)).get();
            nmap[instance] += 1;
            val->instances.push_back(new yocto_instance());
            val->instances.back()->name = instance->name +
                                          std::to_string(nmap[instance]);
            val->instances.back()->frame = base->frame *
                                           translation_frame(pos[i]) *
                                           instance->frame;
            val->instances.back()->shape    = instance->shape;
            val->instances.back()->material = instance->material;
            val->instances.back()->surface  = instance->surface;
        }
    }
    return true;
}

// Json to scene
bool parse_json_object(
    const json& js, yocto_scene* val, const yocto_scene* scene) {
    static const auto def = yocto_scene();
    if (!parse_json_objbegin(js)) return false;
    if (!parse_json_value(js, val->name, "name", def.name)) return false;
    if (!parse_json_objarray(js, val->cameras, "cameras", scene)) return false;
    if (!parse_json_objarray(js, val->textures, "textures", scene))
        return false;
    if (!parse_json_objarray(js, val->voltextures, "voltextures", scene))
        return false;
    if (!parse_json_objarray(js, val->materials, "materials", scene))
        return false;
    if (!parse_json_objarray(js, val->shapes, "shapes", scene)) return false;
    if (!parse_json_objarray(js, val->surfaces, "surfaces", scene))
        return false;
    if (!parse_json_objarray(js, val->instances, "instances", scene))
        return false;
    if (!parse_json_objarray(js, val->environments, "environments", scene))
        return false;
    if (!parse_json_objarray(js, val->nodes, "nodes", scene)) return false;
    if (!parse_json_objarray(js, val->animations, "animations", scene))
        return false;
    if (!parse_json_procedural(js, val, "!!proc", scene)) return false;
    return true;
}

bool dump_json_object(json& js, const yocto_scene* val) {
    return dump_json_object(js, val, val);
}
bool parse_json_object(const json& js, yocto_scene* val) {
    return parse_json_object(js, val, val);
}

// Load a scene in the builtin JSON format.
yocto_scene* load_json_scene(
    const string& filename, bool load_textures, bool skip_missing) {
    auto scope = log_trace_scoped("loading scene {}", filename);
    // initialize
    auto scene = make_unique<yocto_scene>();

    // load jsonz
    auto js = load_json(filename);
    if (js.empty()) return nullptr;

    // deserialize json
    try {
        if (!parse_json_object(js, scene.get())) {
            log_io_error("could not deserialize json {}", filename);
            return nullptr;
        }
    } catch (...) {
        log_io_error("could not deserialize json {}", filename);
        return nullptr;
    }

    // load meshes
    auto dirname = get_dirname(filename);
    for (auto shape : scene->shapes) {
        if (shape->filename == "" || !shape->positions.empty()) continue;
        auto filename = normalize_path(dirname + "/" + shape->filename);
        if (!load_mesh(filename, shape->points, shape->lines, shape->triangles,
                shape->quads, shape->positions, shape->normals,
                shape->texturecoords, shape->colors, shape->radius, false)) {
            if (!skip_missing) return nullptr;
        }
    }

    // load suddivs
    for (auto surface : scene->surfaces) {
        if (surface->filename == "" || !surface->positions.empty()) continue;
        auto filename      = normalize_path(dirname + "/" + surface->filename);
        auto quads_normals = vector<vec4i>();
        auto norm          = vector<vec3f>();
        if (!load_fvmesh(filename, surface->positions_quads, surface->positions,
                quads_normals, norm, surface->texturecoords_quads,
                surface->texturecoords, surface->colors_quads, surface->colors)) {
            if (!skip_missing) return nullptr;
        }
    }

    // skip textures
    if (load_textures) {
        if (!load_scene_textures(scene.get(), dirname, skip_missing, false))
            return nullptr;
    }

    // fix scene
    if (scene->name == "") scene->name = get_filename(filename);
    add_missing_cameras(scene.get());
    add_missing_materials(scene.get());
    add_missing_names(scene.get());
    update_transforms(scene.get());

    // done
    return scene.release();
}

// Save a scene in the builtin JSON format.
bool save_json_scene(const string& filename, const yocto_scene* scene,
    bool save_textures, bool skip_missing) {
    // save json
    auto js = json();
    try {
        if (!dump_json_object(js, scene)) {
            log_io_error("could not serialize json {}", filename);
            return false;
        }
    } catch (...) {
        log_io_error("could not serialize json {}", filename);
        return false;
    }
    if (!save_json(filename, js)) return false;

    // save meshes
    auto dirname = get_dirname(filename);
    for (auto& shape : scene->shapes) {
        if (shape->filename == "") continue;
        auto filename = normalize_path(dirname + "/" + shape->filename);
        if (!save_mesh(filename, shape->points, shape->lines, shape->triangles,
                shape->quads, shape->positions, shape->normals,
                shape->texturecoords, shape->colors, shape->radius)) {
            if (!skip_missing) return false;
        }
    }

    // save subdivs
    for (auto& surface : scene->surfaces) {
        if (surface->filename == "") continue;
        auto filename = normalize_path(dirname + "/" + surface->filename);
        if (!save_fvmesh(filename, surface->positions_quads, surface->positions,
                {}, {}, surface->texturecoords_quads, surface->texturecoords,
                surface->colors_quads, surface->colors)) {
            if (!skip_missing) return false;
        }
    }

    // skip textures
    if (save_textures) {
        if (!save_scene_textures(scene, dirname, skip_missing)) return false;
    }

    // done
    return true;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// OBJ CONVERSION
// -----------------------------------------------------------------------------
namespace ygl {

inline bool operator==(obj_vertex a, obj_vertex b) {
    return a.position == b.position && a.texturecoord == b.texturecoord &&
           a.normal == b.normal;
}

struct obj_vertex_hash {
    size_t operator()(const obj_vertex& v) const {
        auto vh = std::hash<int>();
        auto h  = (size_t)0;
        for (auto i = 0; i < 3; i++)
            h ^= vh((&v.position)[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

// prepare obj line (remove comments and normalize whitespace)
void normalize_obj_line(char* s) {
    while (*s) {
        if (*s == '\t' || *s == '\r' || *s == '\n') {
            *s++ = ' ';
        } else if (*s == '#') {
            *s = 0;
            break;
        } else {
            s++;
        }
    }
}

// parse stream
inline int parse_int(char*& s) {
    if (!*s) return 0;
    while (*s == ' ') s++;
    if (!*s) return 0;
    auto val = 0;
    auto sn  = (*s == '-') ? -1 : 1;
    if (*s == '-' || *s == '+') s++;
    while (*s >= '0' && *s <= '9') val = val * 10 + (*s++ - '0');
    val *= sn;
    return val;
}
inline bool   parse_bool(char*& s) { return (bool)parse_int(s); }
inline double parse_double(char*& s) {
    if (!*s) return 0;
    while (*s == ' ') s++;
    auto val      = 0.0;
    auto mantissa = 0, fractional = 0, fractional_length = 0, exponent = 0;
    auto sn = (*s == '-') ? -1 : 1;
    if (*s == '-' || *s == '+') s++;
    while (*s >= '0' && *s <= '9') mantissa = mantissa * 10 + (*s++ - '0');
    if (*s == '.') {
        s++;
        while (*s >= '0' && *s <= '9') {
            fractional = fractional * 10 + (*s++ - '0');
            fractional_length++;
        }
    }
    mantissa *= sn;
    fractional *= sn;
    if (*s == 'e' || *s == 'E') {
        s++;
        auto en = (*s == '-') ? -1 : 1;
        if (*s == '-' || *s == '+') s++;
        while (*s >= '0' && *s <= '9') exponent = exponent * 10 + (*s++ - '0');
        exponent *= en;
    }
    val = (double)mantissa;
    if (fractional)
        val += fractional * std::pow(10.0, -(double)fractional_length);
    if (exponent) val *= std::pow(10.0, (double)exponent);
    return val;
}
inline float  parse_float(char*& s) { return (float)parse_double(s); }
inline string parse_string(char*& s) {
    if (!*s) return "";
    char buf[4096];
    auto valb = buf;
    while (*s == ' ') s++;
    while (*s && *s != ' ') *valb++ = *s++;
    *valb = 0;
    return buf;
}
inline vec2f parse_vec2f(char*& s) { return {parse_float(s), parse_float(s)}; }
inline vec3f parse_vec3f(char*& s) {
    return {parse_float(s), parse_float(s), parse_float(s)};
}
inline vec4f parse_vec4f(char*& s) {
    return {parse_float(s), parse_float(s), parse_float(s), parse_float(s)};
}
inline vec2i parse_vec2i(char*& s) { return {parse_int(s), parse_int(s)}; }
inline vec3i parse_vec3i(char*& s) {
    return {parse_int(s), parse_int(s), parse_int(s)};
}
inline vec4i parse_vec4i(char*& s) {
    return {parse_int(s), parse_int(s), parse_int(s), parse_int(s)};
}
inline frame3f parse_frame3f(char*& s) {
    if (!*s) return identity_frame3f;
    return {parse_vec3f(s), parse_vec3f(s), parse_vec3f(s), parse_vec3f(s)};
}

inline obj_vertex parse_obj_vertex(char*& s) {
    auto val     = obj_vertex{0, 0, 0};
    val.position = parse_int(s);
    if (*s == '/') {
        s++;
        if (*s == '/') {
            s++;
            val.normal = parse_int(s);
        } else {
            val.texturecoord = parse_int(s);
            if (*s == '/') {
                s++;
                val.normal = parse_int(s);
            }
        }
    }
    return val;
}

// Input for OBJ textures
inline obj_texture_info parse_obj_texture_info(char*& s) {
    // initialize
    auto info = obj_texture_info();

    // get tokens
    auto tokens = vector<string>();
    while (true) {
        auto v = parse_string(s);
        if (v == "") break;
        tokens.push_back(v);
    }
    if (tokens.empty()) return info;

    // texture name
    info.path = normalize_path(tokens.back());

    // texture options
    auto last = string();
    for (auto i = 0; i < tokens.size() - 1; i++) {
        if (tokens[i] == "-bm") info.scale = atof(tokens[i + 1].c_str());
        if (tokens[i] == "-clamp") info.clamp = true;
    }

    return info;
}

// Load obj materials
bool load_mtl(const string& filename, const obj_callbacks& cb, bool flip_tr) {
    // open file
    auto fs = open(filename, "rt");
    if (!fs) return false;

    // currently parsed material
    auto mat   = obj_material();
    auto first = true;

    // read the file line by line
    auto line = ""s;
    char buf[4096];
    while (read_line(fs, line)) {
        // line
        assert(line.size() < 4096);
        memcpy(buf, line.c_str(), line.size() + 1);
        normalize_obj_line(buf);
        auto ss = buf;

        // get command
        auto cmd = parse_string(ss);
        if (cmd == "") continue;

        // possible token values
        if (cmd == "newmtl") {
            if (!first && cb.material) cb.material(mat);
            first    = false;
            mat      = obj_material();
            mat.name = parse_string(ss);
        } else if (cmd == "illum") {
            mat.illum = parse_int(ss);
        } else if (cmd == "Ke") {
            mat.ke = parse_vec3f(ss);
        } else if (cmd == "Kd") {
            mat.kd = parse_vec3f(ss);
        } else if (cmd == "Ks") {
            mat.ks = parse_vec3f(ss);
        } else if (cmd == "Kt") {
            mat.kt = parse_vec3f(ss);
        } else if (cmd == "Tf") {
            mat.kt = {-1, -1, -1};
            mat.kt = parse_vec3f(ss);
            if (mat.kt.y < 0) mat.kt = {mat.kt.x, mat.kt.x, mat.kt.x};
            if (flip_tr) mat.kt = vec3f{1, 1, 1} - mat.kt;
        } else if (cmd == "Tr") {
            auto tr = vec3f{-1, -1, -1};
            tr      = parse_vec3f(ss);
            if (tr.y < 0) tr = {tr.x, tr.x, tr.x};
            mat.op = (tr.x + tr.y + tr.z) / 3;
            if (flip_tr) mat.op = 1 - mat.op;
        } else if (cmd == "Ns") {
            mat.ns = parse_float(ss);
            mat.rs = pow(2 / (mat.ns + 2), 1 / 4.0f);
            if (mat.rs < 0.01f) mat.rs = 0;
            if (mat.rs > 0.99f) mat.rs = 1;
        } else if (cmd == "d") {
            mat.op = parse_float(ss);
        } else if (cmd == "Pr" || cmd == "rs") {
            mat.rs = parse_float(ss);
        } else if (cmd == "map_Ke") {
            mat.ke_txt = parse_obj_texture_info(ss);
        } else if (cmd == "map_Kd") {
            mat.kd_txt = parse_obj_texture_info(ss);
        } else if (cmd == "map_Ks") {
            mat.ks_txt = parse_obj_texture_info(ss);
        } else if (cmd == "map_Tr") {
            mat.kt_txt = parse_obj_texture_info(ss);
        } else if (cmd == "map_d" || cmd == "map_Tr") {
            mat.op_txt = parse_obj_texture_info(ss);
        } else if (cmd == "map_Pr" || cmd == "map_rs") {
            mat.rs_txt = parse_obj_texture_info(ss);
        } else if (cmd == "map_occ" || cmd == "occ") {
            mat.occ_txt = parse_obj_texture_info(ss);
        } else if (cmd == "map_bump" || cmd == "bump") {
            mat.bump_txt = parse_obj_texture_info(ss);
        } else if (cmd == "map_disp" || cmd == "disp") {
            mat.disp_txt = parse_obj_texture_info(ss);
        } else if (cmd == "map_norm" || cmd == "norm") {
            mat.norm_txt = parse_obj_texture_info(ss);
        } else if (cmd == "Ve") {
            mat.ve = parse_vec3f(ss);
        } else if (cmd == "Va") {
            mat.va = parse_vec3f(ss);
        } else if (cmd == "Vd") {
            mat.vd = parse_vec3f(ss);
        } else if (cmd == "Vg") {
            mat.vg = parse_float(ss);
        } else if (cmd == "map_Vd") {
            mat.vd_txt = parse_obj_texture_info(ss);
        }
    }

    // issue current material
    if (!first && cb.material) cb.material(mat);

    // done
    return true;
}

// Load obj extensions
bool load_objx(const string& filename, const obj_callbacks& cb) {
    // open file
    auto fs = open(filename, "rt");
    if (!fs) return false;

    // read the file line by line
    char buf[4096];
    auto line = ""s;
    while (read_line(fs, line)) {
        // line
        memcpy(buf, line.c_str(), line.size() + 1);
        normalize_obj_line(buf);
        auto ss = buf;

        // get command
        auto cmd = parse_string(ss);
        if (cmd == "") continue;

        // possible token values
        if (cmd == "c") {
            auto camera     = obj_camera();
            camera.name     = parse_string(ss);
            camera.ortho    = parse_bool(ss);
            camera.film     = parse_vec2f(ss);
            camera.focal    = parse_float(ss);
            camera.focus    = parse_float(ss);
            camera.aperture = parse_float(ss);
            camera.frame    = parse_frame3f(ss);
            if (cb.camera) cb.camera(camera);
        } else if (cmd == "e") {
            auto environment        = obj_environment();
            environment.name        = parse_string(ss);
            environment.ke          = parse_vec3f(ss);
            environment.ke_txt.path = parse_string(ss);
            if (environment.ke_txt.path == "\"\"") environment.ke_txt.path = "";
            if (cb.environmnet) cb.environmnet(environment);
        } else {
            // unused
        }
    }

    // done
    return true;
}

// Load obj scene
bool load_obj(const string& filename, const obj_callbacks& cb,
    bool geometry_only, bool skip_missing, bool flip_texcoord, bool flip_tr) {
    // open file
    auto fs = open(filename, "rt");
    if (!fs) return false;

    // track vertex size
    auto vert_size = obj_vertex();
    auto verts     = vector<obj_vertex>();  // buffer to avoid reallocation

    // read the file line by line
    char buf[4096];
    auto line = ""s;
    while (read_line(fs, line)) {
        // line
        memcpy(buf, line.c_str(), line.size() + 1);
        normalize_obj_line(buf);
        auto ss = buf;

        // get command
        auto cmd = parse_string(ss);
        if (cmd == "") continue;

        // possible token values
        if (cmd == "v") {
            if (cb.vert) cb.vert(parse_vec3f(ss));
            vert_size.position += 1;
        } else if (cmd == "vn") {
            if (cb.norm) cb.norm(parse_vec3f(ss));
            vert_size.normal += 1;
        } else if (cmd == "vt") {
            auto v = parse_vec2f(ss);
            if (flip_texcoord) v.y = 1 - v.y;
            if (cb.texcoord) cb.texcoord(v);
            vert_size.texturecoord += 1;
        } else if (cmd == "f" || cmd == "l" || cmd == "p") {
            verts.clear();
            while (true) {
                auto vert = parse_obj_vertex(ss);
                if (!vert.position) break;
                if (vert.position < 0)
                    vert.position = vert_size.position + vert.position + 1;
                if (vert.texturecoord < 0)
                    vert.texturecoord = vert_size.texturecoord +
                                        vert.texturecoord + 1;
                if (vert.normal < 0)
                    vert.normal = vert_size.normal + vert.normal + 1;
                verts.push_back(vert);
            }
            if (cmd == "f" && cb.face) cb.face(verts);
            if (cmd == "l" && cb.line) cb.line(verts);
            if (cmd == "p" && cb.point) cb.point(verts);
        } else if (cmd == "o") {
            if (cb.object) cb.object(parse_string(ss));
        } else if (cmd == "usemtl") {
            if (cb.usemtl) cb.usemtl(parse_string(ss));
        } else if (cmd == "g") {
            if (cb.group) cb.group(parse_string(ss));
        } else if (cmd == "s") {
            if (cb.smoothing) cb.smoothing(parse_string(ss));
        } else if (cmd == "mtllib") {
            if (geometry_only) continue;
            auto mtlname = parse_string(ss);
            if (cb.mtllib) cb.mtllib(mtlname);
            auto mtlpath = get_dirname(filename) + "/" + mtlname;
            if (!load_mtl(mtlpath, cb, flip_tr)) {
                if (!skip_missing) return false;
            }
        } else {
            // unused
        }
    }

    // parse extensions if presents
    if (!geometry_only) {
        auto extname    = replace_extension(filename, "objx");
        auto ext_exists = exists_file(extname);
        if (ext_exists) {
            if (!load_objx(extname, cb)) return false;
        }
    }

    // done
    return true;
}

// Loads an OBJ
yocto_scene* load_obj_scene(const string& filename, bool load_textures,
    bool skip_missing, bool split_shapes) {
    auto scope = log_trace_scoped("loading scene {}", filename);
    auto scene = make_unique<yocto_scene>();

    // splitting policy
    auto split_material  = split_shapes;
    auto split_group     = split_shapes;
    auto split_smoothing = split_shapes;

    // current parsing values
    auto matname   = string();
    auto oname     = string();
    auto gname     = string();
    auto smoothing = true;
    auto instance  = (yocto_instance*)nullptr;

    // vertices
    auto opos      = std::deque<vec3f>();
    auto onorm     = std::deque<vec3f>();
    auto otexcoord = std::deque<vec2f>();

    // object maps
    auto tmap = unordered_map<string, yocto_texture*>();
    auto vmap = unordered_map<string, yocto_voltexture*>();
    auto mmap = unordered_map<string, yocto_material*>();

    // vertex maps
    auto name_map     = unordered_map<string, int>();
    auto vert_map     = unordered_map<obj_vertex, int, obj_vertex_hash>();
    auto pos_map      = unordered_map<int, int>();
    auto norm_map     = unordered_map<int, int>();
    auto texcoord_map = unordered_map<int, int>();

    // add object if needed
    auto is_instance_empty = [](yocto_instance* instance) {
        if (instance->surface) {
            return instance->surface->positions.empty();
        } else if (instance->shape) {
            return instance->shape->positions.empty();
        } else {
            return true;
        }
    };
    auto add_instance = [&](yocto_scene* scene, const string& objname,
                            const string& matname, const string& groupname,
                            bool smoothing) {
        if (scene->instances.empty() || objname != scene->instances.back()->name ||
            !is_instance_empty(scene->instances.back())) {
            auto instance = new yocto_instance();
            scene->instances.push_back(instance);
            instance->shape = new yocto_shape();
            scene->shapes.push_back(instance->shape);
        }
        name_map[objname] += 1;
        auto name = (name_map[objname] == 1) ?
                        objname :
                        objname + "_" + std::to_string(name_map[objname] - 1);
        if (objname == "") name = "object" + name;
        auto instance  = scene->instances.back();
        instance->name = name;
        if (instance->shape) instance->shape->name = instance->name;
        if (instance->surface) instance->surface->name = instance->name;
        if (matname != "") {
            auto it = mmap.find(matname);
            if (it == mmap.end()) {
                log_error("missing material {}", matname);
            }
            instance->material = it->second;
        }
        vert_map.clear();
        pos_map.clear();
        norm_map.clear();
        texcoord_map.clear();
        return instance;
    };
    // Parse texture options and name
    auto add_texture = [&scene, &tmap](
                           const obj_texture_info& info, bool force_linear) {
        if (info.path == "") return (yocto_texture*)nullptr;
        if (tmap.find(info.path) != tmap.end()) {
            return tmap.at(info.path);
        }

        // create texture
        auto texture           = new yocto_texture();
        texture->name          = info.path;
        texture->filename      = info.path;
        texture->clamp_to_edge = info.clamp;
        texture->height_scale  = info.scale;
        texture->ldr_as_linear = force_linear || is_hdr_filename(info.path);
        scene->textures.push_back(texture);
        tmap[info.path] = texture;

        return texture;
    };
    // Parse texture options and name
    auto add_voltexture = [&scene, &vmap](
                              const obj_texture_info& info, bool srgb) {
        if (info.path == "") return (yocto_voltexture*)nullptr;
        if (vmap.find(info.path) != vmap.end()) {
            return vmap.at(info.path);
        }

        // create texture
        auto texture      = new yocto_voltexture();
        texture->name     = info.path;
        texture->filename = info.path;
        scene->voltextures.push_back(texture);
        vmap[info.path] = texture;

        return texture;
    };
    // Add  vertices to the current shape
    auto add_verts = [&](const vector<obj_vertex>& verts) {
        for (auto& vert : verts) {
            auto it = vert_map.find(vert);
            if (it != vert_map.end()) continue;
            auto nverts = (int)instance->shape->positions.size();
            vert_map.insert(it, {vert, nverts});
            if (vert.position)
                instance->shape->positions.push_back(opos.at(vert.position - 1));
            if (vert.texturecoord)
                instance->shape->texturecoords.push_back(
                    otexcoord.at(vert.texturecoord - 1));
            if (vert.normal)
                instance->shape->normals.push_back(onorm.at(vert.normal - 1));
        }
    };

    // current objet
    instance = add_instance(scene.get(), "", "", "", true);

    // callbacks
    auto cb     = obj_callbacks();
    cb.vert     = [&](vec3f v) { opos.push_back(v); };
    cb.norm     = [&](vec3f v) { onorm.push_back(v); };
    cb.texcoord = [&](vec2f v) { otexcoord.push_back(v); };
    cb.face     = [&](const vector<obj_vertex>& verts) {
        add_verts(verts);
        if (verts.size() == 4) {
            instance->shape->quads.push_back(
                {vert_map.at(verts[0]), vert_map.at(verts[1]),
                    vert_map.at(verts[2]), vert_map.at(verts[3])});
        } else {
            for (auto i = 2; i < verts.size(); i++)
                instance->shape->triangles.push_back({vert_map.at(verts[0]),
                    vert_map.at(verts[i - 1]), vert_map.at(verts[i])});
        }
    };
    cb.line = [&](const vector<obj_vertex>& verts) {
        add_verts(verts);
        for (auto i = 1; i < verts.size(); i++)
            instance->shape->lines.push_back(
                {vert_map.at(verts[i - 1]), vert_map.at(verts[i])});
    };
    cb.point = [&](const vector<obj_vertex>& verts) {
        add_verts(verts);
        for (auto i = 0; i < verts.size(); i++)
            instance->shape->points.push_back(vert_map.at(verts[i]));
    };
    cb.object = [&](const string& name) {
        oname     = name;
        gname     = "";
        matname   = "";
        smoothing = true;
        instance  = add_instance(scene.get(), oname, matname, gname, smoothing);
    };
    cb.group = [&](const string& name) {
        gname = name;
        if (split_group) {
            instance = add_instance(
                scene.get(), oname, matname, gname, smoothing);
        }
    };
    cb.smoothing = [&](const string& name) {
        smoothing = (name == "on");
        if (split_smoothing) {
            instance = add_instance(
                scene.get(), oname, matname, gname, smoothing);
        }
    };
    cb.usemtl = [&](const string& name) {
        matname = name;
        if (split_material) {
            instance = add_instance(
                scene.get(), oname, matname, gname, smoothing);
        } else {
            if (matname != "") instance->material = mmap.at(matname);
        }
    };
    cb.material = [&](const obj_material& omat) {
        auto mat                    = new yocto_material();
        mat->name                   = omat.name;
        mat->emission               = omat.ke;
        mat->diffuse                = omat.kd;
        mat->specular               = omat.ks;
        mat->transmission           = omat.kt;
        mat->roughness              = omat.rs;
        mat->opacity                = omat.op;
        mat->emission_texture       = add_texture(omat.ke_txt, false);
        mat->diffuse_texture        = add_texture(omat.kd_txt, false);
        mat->specular_texture       = add_texture(omat.ks_txt, false);
        mat->transmission_texture   = add_texture(omat.kt_txt, false);
        mat->opacity_texture        = add_texture(omat.op_txt, true);
        mat->roughness_texture      = add_texture(omat.rs_txt, true);
        mat->occlusion_texture      = add_texture(omat.occ_txt, true);
        mat->bump_texture           = add_texture(omat.bump_txt, true);
        mat->displacement_texture   = add_texture(omat.disp_txt, true);
        mat->normal_texture         = add_texture(omat.norm_txt, true);
        mat->volume_emission        = omat.ve;
        mat->volume_albedo          = omat.va;
        mat->volume_density         = omat.vd;
        mat->volume_phaseg          = omat.vg;
        mat->volume_density_texture = add_voltexture(omat.vd_txt, false);
        scene->materials.push_back(mat);
        mmap[mat->name] = mat;
    };
    cb.camera = [&](const obj_camera& ocam) {
        auto camera            = new yocto_camera();
        camera->name           = ocam.name;
        camera->orthographic   = ocam.ortho;
        camera->film_size      = ocam.film;
        camera->focal_length   = ocam.focal;
        camera->focus_distance = ocam.focus;
        camera->lens_aperture  = ocam.aperture;
        camera->frame          = ocam.frame;
        scene->cameras.push_back(camera);
    };
    cb.environmnet = [&](const obj_environment& oenv) {
        auto environment              = new yocto_environment();
        environment->name             = oenv.name;
        environment->emission         = oenv.ke;
        environment->emission_texture = add_texture(oenv.ke_txt, true);
        scene->environments.push_back(environment);
    };

    // Parse obj
    if (!load_obj(filename, cb, false, skip_missing)) return nullptr;

    // cleanup empty
    // TODO: delete unused
    for (auto idx = 0; idx < scene->instances.size(); idx++) {
        if (!is_instance_empty(scene->instances[idx])) continue;
        auto instance = scene->instances[idx];
        if (instance->shape) {
            scene->shapes.erase(std::find(
                scene->shapes.begin(), scene->shapes.end(), instance->shape));
        }
        if (instance->surface) {
            scene->surfaces.erase(std::find(scene->surfaces.begin(),
                scene->surfaces.end(), instance->surface));
        }
        scene->instances.erase(scene->instances.begin() + idx);
        idx--;
    }

    // load textures
    auto dirname = get_dirname(filename);
    if (load_textures) {
        if (!load_scene_textures(scene.get(), dirname, skip_missing, false))
            return nullptr;
    }

    // fix scene
    scene->name = get_filename(filename);
    add_missing_cameras(scene.get());
    add_missing_materials(scene.get());
    add_missing_names(scene.get());
    update_transforms(scene.get());

    // done
    return scene.release();
}

bool save_mtl(
    const string& filename, const yocto_scene* scene, bool flip_tr = true) {
    // open
    auto fs = open(filename, "wt");
    if (!fs) return false;

    // for each material, dump all the values
    for (auto mat : scene->materials) {
        print(fs, "newmtl {}\n", mat->name);
        print(fs, "  illum 2\n");
        if (mat->emission != zero3f) print(fs, "  Ke {}\n", mat->emission);
        if (mat->diffuse != zero3f) print(fs, "  Kd {}\n", mat->diffuse);
        if (mat->specular != zero3f) print(fs, "  Ks {}\n", mat->specular);
        if (mat->transmission != zero3f)
            print(fs, "  Kt {}\n", mat->transmission);
        if (mat->roughness != 1.0f)
            print(fs, "  Ns {}\n",
                (int)clamp(
                    2 / pow(mat->roughness + 1e-10f, 4.0f) - 2, 0.0f, 1.0e12f));
        if (mat->opacity != 1.0f) print(fs, "  d {}\n", mat->opacity);
        if (mat->roughness != -1.0f) print(fs, "  Pr {}\n", mat->roughness);
        if (mat->emission_texture)
            print(fs, "  map_Ke {}\n", mat->emission_texture->filename);
        if (mat->diffuse_texture)
            print(fs, "  map_Kd {}\n", mat->diffuse_texture->filename);
        if (mat->specular_texture)
            print(fs, "  map_Ks {}\n", mat->specular_texture->filename);
        if (mat->transmission_texture)
            print(fs, "  map_Kt {}\n", mat->transmission_texture->filename);
        if (mat->opacity_texture && mat->opacity_texture != mat->diffuse_texture)
            print(fs, "  map_d  {}\n", mat->opacity_texture->filename);
        if (mat->roughness_texture)
            print(fs, "  map_Pr {}\n", mat->roughness_texture->filename);
        if (mat->occlusion_texture)
            print(fs, "  map_occ {}\n", mat->occlusion_texture->filename);
        if (mat->bump_texture)
            print(fs, "  map_bump {}\n", mat->bump_texture->filename);
        if (mat->displacement_texture)
            print(fs, "  map_disp {}\n", mat->displacement_texture->filename);
        if (mat->normal_texture)
            print(fs, "  map_norm {}\n", mat->normal_texture->filename);
        if (mat->volume_emission != zero3f)
            print(fs, "  Ve {}\n", mat->volume_emission);
        if (mat->volume_density != zero3f)
            print(fs, "  Vd {}\n", mat->volume_density);
        if (mat->volume_albedo != zero3f)
            print(fs, "  Va {}\n", mat->volume_albedo);
        if (mat->volume_phaseg != 0) print(fs, "  Vg {}\n", mat->volume_phaseg);
        if (mat->volume_density_texture)
            print(fs, "  map_Vd {}\n", mat->volume_density_texture->filename);
        print(fs, "\n");
    }

    // done
    return true;
}

bool save_objx(const string& filename, const yocto_scene* scene) {
    // scene
    auto fs = open(filename, "wt");
    if (!fs) return false;

    // cameras
    for (auto camera : scene->cameras) {
        print(fs, "c {} {} {} {} {} {} {}\n", camera->name,
            (int)camera->orthographic, camera->film_size, camera->focal_length,
            camera->focus_distance, camera->lens_aperture, camera->frame);
    }

    // environments
    for (auto environment : scene->environments) {
        print(fs, "e {} {} {} {}\n", environment->name.c_str(),
            environment->emission,
            ((environment->emission_texture) ?
                    environment->emission_texture->filename.c_str() :
                    "\"\""),
            environment->frame);
    }

    // done
    return true;
}

string to_string(const obj_vertex& v) {
    auto s = std::to_string(v.position);
    if (v.texturecoord) {
        s += "/" + std::to_string(v.texturecoord);
        if (v.normal) s += "/" + std::to_string(v.normal);
    } else {
        if (v.normal) s += "//" + std::to_string(v.normal);
    }
    return s;
}

bool save_obj(const string& filename, const yocto_scene* scene,
    bool flip_texcoord = true) {
    // scene
    auto fs = open(filename, "wt");
    if (!fs) return false;

    // material library
    if (!scene->materials.empty()) {
        auto mtlname = replace_extension(get_filename(filename), "mtl");
        print(fs, "mtllib {}\n", mtlname);
    }

    // shapes
    auto offset = obj_vertex{0, 0, 0};
    for (auto instance : scene->instances) {
        if (!instance->surface) {
            print(fs, "o {}\n", instance->name);
            if (instance->material)
                print(fs, "usemtl {}\n", instance->material->name);
            if (instance->frame == identity_frame3f) {
                for (auto& p : instance->shape->positions)
                    print(fs, "v {}\n", p);
                for (auto& n : instance->shape->normals)
                    print(fs, "vn {}\n", n);
                for (auto& t : instance->shape->texturecoords)
                    print(fs, "vt {}\n",
                        vec2f{t.x, (flip_texcoord) ? 1 - t.y : t.y});
            } else {
                for (auto& pp : instance->shape->positions) {
                    print(fs, "v {}\n", transform_point(instance->frame, pp));
                }
                for (auto& nn : instance->shape->normals) {
                    print(fs, "vn {}\n",
                        transform_direction(instance->frame, nn));
                }
                for (auto& t : instance->shape->texturecoords)
                    print(fs, "vt {}\n",
                        vec2f{t.x, (flip_texcoord) ? 1 - t.y : t.y});
            }
            auto mask = obj_vertex{1,
                instance->shape->texturecoords.empty() ? 0 : 1,
                instance->shape->normals.empty() ? 0 : 1};
            auto vert = [mask, offset](int i) {
                return obj_vertex{(i + offset.position + 1) * mask.position,
                    (i + offset.texturecoord + 1) * mask.texturecoord,
                    (i + offset.normal + 1) * mask.normal};
            };
            for (auto& p : instance->shape->points) {
                print(fs, "p {}\n", to_string(vert(p)));
            }
            for (auto& l : instance->shape->lines) {
                print(fs, "l {} {}\n", to_string(vert(l.x)),
                    to_string(vert(l.y)));
            }
            for (auto& t : instance->shape->triangles) {
                print(fs, "f {} {} {}\n", to_string(vert(t.x)),
                    to_string(vert(t.y)), to_string(vert(t.z)));
            }
            for (auto& q : instance->shape->quads) {
                if (q.z == q.w) {
                    print(fs, "f {} {} {}\n", to_string(vert(q.x)),
                        to_string(vert(q.y)), to_string(vert(q.z)));
                } else {
                    print(fs, "f {} {} {} {}\n", to_string(vert(q.x)),
                        to_string(vert(q.y)), to_string(vert(q.z)),
                        to_string(vert(q.w)));
                }
            }
            offset.position += instance->shape->positions.size();
            offset.texturecoord += instance->shape->texturecoords.size();
            offset.normal += instance->shape->normals.size();
        } else {
            print(fs, "o {}\n", instance->name);
            if (instance->material)
                print(fs, "usemtl {}\n", instance->material->name);
            if (instance->frame == identity_frame3f) {
                for (auto& p : instance->surface->positions)
                    print(fs, "v {}\n", p);
                for (auto& t : instance->surface->texturecoords)
                    print(fs, "vt {}\n",
                        vec2f{t.x, (flip_texcoord) ? 1 - t.y : t.y});
            } else {
                for (auto& pp : instance->surface->positions) {
                    auto p = transform_point(instance->frame, pp);
                    print(fs, "v {}\n", p);
                }
                for (auto& t : instance->surface->texturecoords)
                    print(fs, "vt {}\n",
                        vec2f{t.x, (flip_texcoord) ? 1 - t.y : t.y});
            }
            if (!instance->surface->texturecoords.empty()) {
                auto vert = [offset](int ip, int it) {
                    return obj_vertex{ip + offset.position + 1,
                        it + offset.texturecoord + 1, 0};
                };
                for (auto i = 0; i < instance->surface->positions_quads.size();
                     i++) {
                    auto qp = instance->surface->positions_quads[i];
                    auto qt = instance->surface->texturecoords_quads[i];
                    if (qp.z == qp.w) {
                        print(fs, "f {} {} {}\n", to_string(vert(qp.x, qt.x)),
                            to_string(vert(qp.y, qt.y)),
                            to_string(vert(qp.z, qt.z)));
                    } else {
                        print(fs, "f {} {} {} {}\n", to_string(vert(qp.x, qt.x)),
                            to_string(vert(qp.y, qt.y)),
                            to_string(vert(qp.z, qt.z)),
                            to_string(vert(qp.w, qt.w)));
                    }
                }
            } else {
                auto vert = [offset](int ip) {
                    return obj_vertex{ip + offset.position + 1, 0, 0};
                };
                for (auto& q : instance->surface->positions_quads) {
                    if (q.z == q.w) {
                        print(fs, "f {} {} {}\n", to_string(vert(q.x)),
                            to_string(vert(q.y)), to_string(vert(q.z)));
                    } else {
                        print(fs, "f {} {} {} {}\n", to_string(vert(q.x)),
                            to_string(vert(q.y)), to_string(vert(q.z)),
                            to_string(vert(q.w)));
                    }
                }
            }
            offset.position += instance->surface->positions.size();
            offset.texturecoord += instance->surface->texturecoords.size();
        }
    }

    return true;
}

bool save_obj_scene(const string& filename, const yocto_scene* scene,
    bool save_textures, bool skip_missing) {
    if (!save_obj(filename, scene, true)) return false;
    if (!scene->materials.empty()) {
        if (!save_mtl(replace_extension(filename, ".mtl"), scene, true))
            return false;
    }
    if (!scene->cameras.empty() || !scene->environments.empty()) {
        if (!save_objx(replace_extension(filename, ".objx"), scene))
            return false;
    }

    // skip textures if needed
    auto dirname = get_dirname(filename);
    if (save_textures) {
        if (!save_scene_textures(scene, dirname, skip_missing)) return false;
    }

    // done
    return true;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// GLTF CONVESION
// -----------------------------------------------------------------------------
namespace ygl {

static bool startswith(const string& str, const string& substr) {
    if (str.length() < substr.length()) return false;
    for (auto i = 0; i < substr.length(); i++)
        if (str[i] != substr[i]) return false;
    return true;
}

// convert gltf to scene
bool gltf_to_scene(yocto_scene* scene, const json& gltf, const string& dirname) {
    // convert textures
    if (gltf.count("images")) {
        for (auto iid = 0; iid < gltf.at("images").size(); iid++) {
            auto& gimg        = gltf.at("images").at(iid);
            auto  texture     = new yocto_texture();
            texture->name     = gimg.value("name", ""s);
            texture->filename = (startswith(gimg.value("uri", ""s), "data:")) ?
                                    string("[glTF-inline].png") :
                                    gimg.value("uri", ""s);
            scene->textures.push_back(texture);
        }
    }

    // load buffers
    auto bmap = vector<vector<byte>>();
    if (gltf.count("buffers")) {
        bmap.resize(gltf.at("buffers").size());
        for (auto bid = 0; bid < gltf.at("buffers").size(); bid++) {
            auto& gbuf = gltf.at("buffers").at(bid);
            auto& data = bmap.at(bid);
            auto  uri  = gbuf.value("uri", ""s);
            if (uri == "") continue;
            if (startswith(uri, "data:")) {
                // assume it is base64 and find ','
                auto pos = uri.find(',');
                if (pos == uri.npos) {
                    return false;
                }
                // decode
                auto data_char = base64_decode(uri.substr(pos + 1));
                data = vector<unsigned char>((unsigned char*)data_char.c_str(),
                    (unsigned char*)data_char.c_str() + data_char.length());
            } else {
                auto filename = normalize_path(dirname + "/" + uri);
                data          = load_binary(filename);
                if (data.empty()) return false;
            }
            if (gbuf.value("byteLength", -1) != data.size()) {
                return false;
            }
        }
    }

    // add a texture
    auto add_texture = [scene, &gltf](const json& ginfo, bool force_linear) {
        if (!gltf.count("images") || !gltf.count("textures"))
            return (yocto_texture*)nullptr;
        if (ginfo.is_null() || ginfo.empty()) return (yocto_texture*)nullptr;
        if (ginfo.value("index", -1) < 0) return (yocto_texture*)nullptr;
        auto& gtxt = gltf.at("textures").at(ginfo.value("index", -1));
        if (gtxt.empty() || gtxt.value("source", -1) < 0)
            return (yocto_texture*)nullptr;
        auto texture = scene->textures.at(gtxt.value("source", -1));
        if (!gltf.count("samplers") || gtxt.value("sampler", -1) < 0)
            return texture;
        auto& gsmp = gltf.at("samplers").at(gtxt.value("sampler", -1));
        texture->clamp_to_edge = gsmp.value("wrapS", ""s) == "ClampToEdge" ||
                                 gsmp.value("wrapT", ""s) == "ClampToEdge";
        texture->height_scale = gsmp.value("scale", 1.0f) *
                                gsmp.value("strength", 1.0f);
        texture->ldr_as_linear = force_linear ||
                                 is_hdr_filename(texture->filename);
        return texture;
    };

    // convert materials
    if (gltf.count("materials")) {
        for (auto mid = 0; mid < gltf.at("materials").size(); mid++) {
            auto& gmat    = gltf.at("materials").at(mid);
            auto  mat     = new yocto_material();
            mat->name     = gmat.value("name", ""s);
            mat->emission = gmat.value("emissiveFactor", zero3f);
            if (gmat.count("emissiveTexture"))
                mat->emission_texture = add_texture(
                    gmat.at("emissiveTexture"), false);
            if (gmat.count("extensions") &&
                gmat.at("extensions").count("KHR_materials_pbrSpecularGlossiness")) {
                mat->base_metallic = false;
                mat->gltf_textures = true;
                auto& gsg          = gmat.at("extensions")
                                .at("KHR_materials_pbrSpecularGlossiness");
                auto kb        = gsg.value("diffuseFactor", vec4f{1, 1, 1, 1});
                mat->diffuse   = {kb.x, kb.y, kb.z};
                mat->opacity   = kb.w;
                mat->specular  = gsg.value("specularFactor", vec3f{1, 1, 1});
                mat->roughness = 1 - gsg.value("glossinessFactor", 1.0f);
                if (gsg.count("diffuseTexture"))
                    mat->diffuse_texture = add_texture(
                        gsg.at("diffuseTexture"), false);
                if (gsg.count("specularGlossinessTexture"))
                    mat->specular_texture = add_texture(
                        gsg.at("specularGlossinessTexture"), false);
                mat->roughness_texture = mat->specular_texture;
            } else if (gmat.count("pbrMetallicRoughness")) {
                mat->base_metallic = true;
                mat->gltf_textures = true;
                auto& gmr          = gmat.at("pbrMetallicRoughness");
                auto  kb      = gmr.value("baseColorFactor", vec4f{1, 1, 1, 1});
                mat->diffuse  = {kb.x, kb.y, kb.z};
                mat->opacity  = kb.w;
                auto km       = gmr.value("metallicFactor", 1.0f);
                mat->specular = {km, km, km};
                mat->roughness = gmr.value("roughnessFactor", 1.0f);
                if (gmr.count("baseColorTexture"))
                    mat->diffuse_texture = add_texture(
                        gmr.at("baseColorTexture"), false);
                if (gmr.count("metallicRoughnessTexture"))
                    mat->specular_texture = add_texture(
                        gmr.at("metallicRoughnessTexture"), true);
                mat->roughness_texture = mat->specular_texture;
            }
            if (gmat.count("occlusionTexture"))
                mat->occlusion_texture = add_texture(
                    gmat.at("occlusionTexture"), true);
            if (gmat.count("normalTexture"))
                mat->normal_texture = add_texture(gmat.at("normalTexture"), true);
            mat->double_sided = gmat.value("doubleSided", true);
            scene->materials.push_back(mat);
        }
    }

    // get values from accessors
    auto accessor_values =
        [&gltf, &bmap](const json& gacc,
            bool normalize = false) -> vector<std::array<double, 4>> {
        auto gview  = gltf.at("bufferViews").at(gacc.value("bufferView", -1));
        auto data   = bmap.at(gview.value("buffer", -1)).data();
        auto offset = gacc.value("byteOffset", 0) + gview.value("byteOffset", 0);
        auto stride      = gview.value("byteStride", 0);
        auto compTypeNum = gacc.value("componentType", 5123);
        auto count       = gacc.value("count", -1);
        auto type        = gacc.value("type", ""s);
        auto ncomp       = 0;
        if (type == "SCALAR") ncomp = 1;
        if (type == "VEC2") ncomp = 2;
        if (type == "VEC3") ncomp = 3;
        if (type == "VEC4") ncomp = 4;
        auto compSize = 1;
        if (compTypeNum == 5122 || compTypeNum == 5123) {
            compSize = 2;
        }
        if (compTypeNum == 5124 || compTypeNum == 5125 || compTypeNum == 5126) {
            compSize = 4;
        }
        if (!stride) stride = compSize * ncomp;
        auto vals = vector<std::array<double, 4>>(count, {{0.0, 0.0, 0.0, 1.0}});
        for (auto i = 0; i < count; i++) {
            auto d = data + offset + i * stride;
            for (auto c = 0; c < ncomp; c++) {
                if (compTypeNum == 5120) {  // char
                    vals[i][c] = (double)(*(char*)d);
                    if (normalize) vals[i][c] /= SCHAR_MAX;
                } else if (compTypeNum == 5121) {  // byte
                    vals[i][c] = (double)(*(byte*)d);
                    if (normalize) vals[i][c] /= UCHAR_MAX;
                } else if (compTypeNum == 5122) {  // short
                    vals[i][c] = (double)(*(short*)d);
                    if (normalize) vals[i][c] /= SHRT_MAX;
                } else if (compTypeNum == 5123) {  // unsigned short
                    vals[i][c] = (double)(*(unsigned short*)d);
                    if (normalize) vals[i][c] /= USHRT_MAX;
                } else if (compTypeNum == 5124) {  // int
                    vals[i][c] = (double)(*(int*)d);
                    if (normalize) vals[i][c] /= INT_MAX;
                } else if (compTypeNum == 5125) {  // unsigned int
                    vals[i][c] = (double)(*(unsigned int*)d);
                    if (normalize) vals[i][c] /= UINT_MAX;
                } else if (compTypeNum == 5126) {  // float
                    vals[i][c] = (*(float*)d);
                }
                d += compSize;
            }
        }
        return vals;
    };

    // convert meshes
    auto meshes = vector<vector<tuple<yocto_shape*, yocto_material*>>>();
    if (gltf.count("meshes")) {
        for (auto mid = 0; mid < gltf.at("meshes").size(); mid++) {
            auto& gmesh = gltf.at("meshes").at(mid);
            meshes.push_back({});
            auto sid = 0;
            for (auto& gprim : gmesh.value("primitives", json::array())) {
                if (!gprim.count("attributes")) continue;
                auto shape  = new yocto_shape();
                shape->name = gmesh.value("name", ""s) +
                              ((sid) ? std::to_string(sid) : string());
                sid++;
                for (json::iterator gattr_it = gprim.at("attributes").begin();
                     gattr_it != gprim.at("attributes").end(); ++gattr_it) {
                    auto  semantic = gattr_it.key();
                    auto& gacc     = gltf.at("accessors")
                                     .at(gattr_it.value().get<int>());
                    auto vals = accessor_values(gacc);
                    if (semantic == "POSITION") {
                        shape->positions.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shape->positions.push_back({(float)vals[i][0],
                                (float)vals[i][1], (float)vals[i][2]});
                    } else if (semantic == "NORMAL") {
                        shape->normals.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shape->normals.push_back({(float)vals[i][0],
                                (float)vals[i][1], (float)vals[i][2]});
                    } else if (semantic == "TEXCOORD" ||
                               semantic == "TEXCOORD_0") {
                        shape->texturecoords.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shape->texturecoords.push_back(
                                {(float)vals[i][0], (float)vals[i][1]});
                    } else if (semantic == "COLOR" || semantic == "COLOR_0") {
                        shape->colors.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shape->colors.push_back(
                                {(float)vals[i][0], (float)vals[i][1],
                                    (float)vals[i][2], (float)vals[i][3]});
                    } else if (semantic == "TANGENT") {
                        shape->tangentspaces.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shape->tangentspaces.push_back(
                                {(float)vals[i][0], (float)vals[i][1],
                                    (float)vals[i][2], (float)vals[i][3]});
                        for (auto& t : shape->tangentspaces) t.w = -t.w;
                    } else if (semantic == "RADIUS") {
                        shape->radius.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shape->radius.push_back((float)vals[i][0]);
                    } else {
                        // ignore
                    }
                }
                // indices
                auto mode = gprim.value("mode", 4);
                if (!gprim.count("indices")) {
                    if (mode == 4) {
                        // triangles
                        shape->triangles.reserve(shape->positions.size() / 3);
                        for (auto i = 0; i < shape->positions.size() / 3; i++)
                            shape->triangles.push_back(
                                {i * 3 + 0, i * 3 + 1, i * 3 + 2});
                    } else if (mode == 6) {
                        // triangle fan
                        shape->triangles.reserve(shape->positions.size() - 2);
                        for (auto i = 2; i < shape->positions.size(); i++)
                            shape->triangles.push_back({0, i - 1, i});
                    } else if (mode == 5) {
                        // triangle strip
                        shape->triangles.reserve(shape->positions.size() - 2);
                        for (auto i = 2; i < shape->positions.size(); i++)
                            shape->triangles.push_back({i - 2, i - 1, i});
                    } else if (mode == 1) {
                        // lines
                        shape->lines.reserve(shape->positions.size() / 2);
                        for (auto i = 0; i < shape->positions.size() / 2; i++)
                            shape->lines.push_back({i * 2 + 0, i * 2 + 1});
                    } else if (mode == 2) {
                        // line loop
                        shape->lines.reserve(shape->positions.size());
                        for (auto i = 1; i < shape->positions.size(); i++)
                            shape->lines.push_back({i - 1, i});
                        shape->lines.back() = {
                            (int)shape->positions.size() - 1, 0};
                    } else if (mode == 3) {
                        // line strip
                        shape->lines.reserve(shape->positions.size() - 1);
                        for (auto i = 1; i < shape->positions.size(); i++)
                            shape->lines.push_back({i - 1, i});
                    } else if (mode == -1 || mode == 0) {
                        // points
                        printf("points not supported\n");
                    } else {
                        log_error("unknown primitive type");
                    }
                } else {
                    auto indices = accessor_values(
                        gltf.at("accessors").at(gprim.value("indices", -1)),
                        false);
                    if (mode == 4) {
                        // triangles
                        shape->triangles.reserve(indices.size());
                        for (auto i = 0; i < indices.size() / 3; i++)
                            shape->triangles.push_back(
                                {(int)indices[i * 3 + 0][0],
                                    (int)indices[i * 3 + 1][0],
                                    (int)indices[i * 3 + 2][0]});
                    } else if (mode == 6) {
                        // triangle fan
                        shape->triangles.reserve(indices.size() - 2);
                        for (auto i = 2; i < indices.size(); i++)
                            shape->triangles.push_back({(int)indices[0][0],
                                (int)indices[i - 1][0], (int)indices[i][0]});
                    } else if (mode == 5) {
                        // triangle strip
                        shape->triangles.reserve(indices.size() - 2);
                        for (auto i = 2; i < indices.size(); i++)
                            shape->triangles.push_back({(int)indices[i - 2][0],
                                (int)indices[i - 1][0], (int)indices[i][0]});
                    } else if (mode == 1) {
                        // lines
                        shape->lines.reserve(indices.size() / 2);
                        for (auto i = 0; i < indices.size() / 2; i++)
                            shape->lines.push_back({(int)indices[i * 2 + 0][0],
                                (int)indices[i * 2 + 1][0]});
                    } else if (mode == 2) {
                        // line loop
                        shape->lines.reserve(indices.size());
                        for (auto i = 1; i < indices.size(); i++)
                            shape->lines.push_back(
                                {(int)indices[i - 1][0], (int)indices[i][0]});
                        shape->lines.back() = {
                            (int)indices[indices.size() - 1][0],
                            (int)indices[0][0]};
                    } else if (mode == 3) {
                        // line strip
                        shape->lines.reserve(indices.size() - 1);
                        for (auto i = 1; i < indices.size(); i++)
                            shape->lines.push_back(
                                {(int)indices[i - 1][0], (int)indices[i][0]});
                    } else if (mode == -1 || mode == 0) {
                        // points
                        printf("points not supported\n");
                    } else {
                        log_error("unknown primitive type");
                    }
                }
                auto mat = (gprim.count("material")) ?
                               scene->materials.at(gprim.value("material", -1)) :
                               nullptr;
                meshes.back().push_back({shape, mat});
                scene->shapes.push_back(shape);
            }
        }
    }

    // convert cameras
    if (gltf.count("cameras")) {
        for (auto cid = 0; cid < gltf.at("cameras").size(); cid++) {
            auto& gcam           = gltf.at("cameras").at(cid);
            auto  camera         = new yocto_camera();
            camera->name         = gcam.value("name", ""s);
            camera->orthographic = gcam.value("type", ""s) == "orthographic";
            if (camera->orthographic) {
                printf("orthographic not supported well\n");
                auto ortho = gcam.value("orthographic", json::object());
                set_camera_fovy(camera, ortho.value("ymag", 0.0f),
                    ortho.value("xmag", 0.0f) / ortho.value("ymag", 0.0f));
                camera->focus_distance = maxf;
                camera->lens_aperture  = 0;
            } else {
                auto persp = gcam.value("perspective", json::object());
                set_camera_fovy(camera, persp.value("yfov", 1.0f),
                    persp.value("aspectRatio", 1.0f));
                camera->focus_distance = maxf;
                camera->lens_aperture  = 0;
            }
            scene->cameras.push_back(camera);
        }
    }

    // convert nodes
    if (gltf.count("nodes")) {
        for (auto nid = 0; nid < gltf.at("nodes").size(); nid++) {
            auto& gnde = gltf.at("nodes").at(nid);
            auto  node = new yocto_scene_node();
            node->name = gnde.value("name", ""s);
            if (gnde.count("camera"))
                node->camera = scene->cameras[gnde.value("camera", 0)];
            node->translation = gnde.value("translation", zero3f);
            node->rotation    = gnde.value("rotation", vec4f{0, 0, 0, 1});
            node->scale       = gnde.value("scale", vec3f{1, 1, 1});
            node->local = mat_to_frame(gnde.value("matrix", identity_mat4f));
            scene->nodes.push_back(node);
        }

        // set up parent pointers
        for (auto nid = 0; nid < gltf.at("nodes").size(); nid++) {
            auto& gnde = gltf.at("nodes").at(nid);
            if (!gnde.count("children")) continue;
            auto node = scene->nodes[nid];
            for (auto& cid : gnde.at("children"))
                scene->nodes[cid.get<int>()]->parent = node;
        }

        // set up instances
        for (auto nid = 0; nid < gltf.at("nodes").size(); nid++) {
            auto& gnde = gltf.at("nodes").at(nid);
            if (!gnde.count("mesh")) continue;
            auto  node = scene->nodes[nid];
            auto& shps = meshes.at(gnde.value("mesh", 0));
            if (shps.empty()) continue;
            if (shps.size() == 1) {
                node->instance           = new yocto_instance();
                node->instance->name     = node->name;
                node->instance->shape    = get<0>(shps[0]);
                node->instance->material = get<1>(shps[0]);
                scene->instances.push_back(node->instance);
            } else {
                for (auto shape : shps) {
                    auto child      = new yocto_scene_node();
                    child->name     = node->name + "_" + get<0>(shape)->name;
                    child->parent   = node;
                    child->instance = new yocto_instance();
                    child->instance->name     = child->name;
                    child->instance->shape    = get<0>(shape);
                    child->instance->material = get<1>(shape);
                    scene->instances.push_back(child->instance);
                }
            }
        }
    }

    // convert animations
    if (gltf.count("animations")) {
        for (auto& ganm : gltf.at("animations")) {
            auto aid         = 0;
            auto sampler_map = unordered_map<vec2i, int>();
            for (auto& gchannel : ganm.at("channels")) {
                auto path_ = gchannel.at("target").at("path").get<string>();
                auto path  = -1;
                if (path_ == "translation") path = 0;
                if (path_ == "rotation") path = 1;
                if (path_ == "scale") path = 2;
                if (path_ == "weights") path = 3;
                if (sampler_map.find({gchannel.at("sampler").get<int>(), path}) ==
                    sampler_map.end()) {
                    auto& gsampler = ganm.at("samplers")
                                         .at(gchannel.at("sampler").get<int>());
                    auto animation  = new yocto_animation();
                    animation->name = (ganm.count("name") ?
                                              ganm.value("name", ""s) :
                                              "anim") +
                                      std::to_string(aid++);
                    animation->animation_group = ganm.value("name", ""s);
                    auto input_view            = accessor_values(
                        gltf.at("accessors").at(gsampler.value("input", -1)));
                    animation->keyframes_times.resize(input_view.size());
                    for (auto i = 0; i < input_view.size(); i++)
                        animation->keyframes_times[i] = input_view[i][0];
                    auto type = gsampler.value("interpolation", "LINEAR");
                    if (type == "LINEAR")
                        animation->interpolation_type = yocto_interpolation_type::linear;
                    if (type == "STEP")
                        animation->interpolation_type = yocto_interpolation_type::step;
                    if (type == "CUBICSPLINE")
                        animation->interpolation_type = yocto_interpolation_type::bezier;
                    auto output_view = accessor_values(
                        gltf.at("accessors").at(gsampler.value("output", -1)));
                    switch (path) {
                        case 0: {  // translation
                            animation->translation_keyframes.reserve(
                                output_view.size());
                            for (auto i = 0; i < output_view.size(); i++)
                                animation->translation_keyframes.push_back(
                                    {(float)output_view[i][0],
                                        (float)output_view[i][1],
                                        (float)output_view[i][2]});
                        } break;
                        case 1: {  // rotation
                            animation->rotation_keyframes.reserve(
                                output_view.size());
                            for (auto i = 0; i < output_view.size(); i++)
                                animation->rotation_keyframes.push_back(
                                    {(float)output_view[i][0],
                                        (float)output_view[i][1],
                                        (float)output_view[i][2],
                                        (float)output_view[i][3]});
                        } break;
                        case 2: {  // scale
                            animation->scale_keyframes.reserve(
                                output_view.size());
                            for (auto i = 0; i < output_view.size(); i++)
                                animation->scale_keyframes.push_back(
                                    {(float)output_view[i][0],
                                        (float)output_view[i][1],
                                        (float)output_view[i][2]});
                        } break;
                        case 3: {  // weights
                            printf("weights not supported for now\n");
#if 0
                        // get a node that it refers to
                        auto ncomp = 0;
                        auto gnode = gltf->get(gchannel->target->node);
                        auto gmesh = gltf->get(gnode->mesh);
                        if (gmesh) {
                            for (auto gshp : gmesh->primitives) {
                                ncomp = max((int)gshp->targets.size(), ncomp);
                            }
                        }
                        if (ncomp) {
                            auto values = vector<float>();
                            values.reserve(output_view.size());
                            for (auto i = 0; i < output_view.size(); i++)
                                values.push_back(output_view.get(i));
                            animation->weights.resize(values.size() / ncomp);
                            for (auto i = 0; i < animation->weights.size(); i++) {
                                animation->weights[i].resize(ncomp);
                                for (auto j = 0; j < ncomp; j++)
                                    animation->weights[i][j] = values[i * ncomp + j];
                            }
                        }
#endif
                        } break;
                        default: { return false; }
                    }
                    sampler_map[{gchannel.at("sampler").get<int>(),
                        path}] = (int)scene->animations.size();
                    scene->animations.push_back(animation);
                }
                scene
                    ->animations[sampler_map.at(
                        {gchannel.at("sampler").get<int>(), path})]
                    ->node_targets.push_back(
                        scene->nodes[(int)gchannel.at("target").at("node").get<int>()]);
            }
        }
    }

    return true;
}

// Load a scene
yocto_scene* load_gltf_scene(
    const string& filename, bool load_textures, bool skip_missing) {
    // initialization
    auto scene = make_unique<yocto_scene>();

    // convert json
    auto js = load_json(filename);
    if (js.empty()) return nullptr;
    try {
        if (!gltf_to_scene(scene.get(), js, get_dirname(filename)))
            return nullptr;
    } catch (...) {
        return nullptr;
    }

    // load textures
    auto dirname = get_dirname(filename);
    if (load_textures) {
        if (!load_scene_textures(scene.get(), dirname, skip_missing, false))
            return nullptr;
    }

    // fix scene
    scene->name = get_filename(filename);
    add_missing_cameras(scene.get());
    add_missing_materials(scene.get());
    add_missing_names(scene.get());
    update_transforms(scene.get());

    // fix cameras
    auto bbox = compute_scene_bounds(scene.get());
    for (auto camera : scene->cameras) {
        auto center = (bbox.min + bbox.max) / 2;
        auto dist   = dot(-camera->frame.z, center - camera->frame.o);
        if (dist > 0) camera->focus_distance = dist;
    }

    // done
    return scene.release();
}

// convert gltf scene to json
bool scene_to_gltf(const yocto_scene* scene, json& js) {
    // init to emprt object
    js = json::object();

    // start creating json
    js["asset"]["version"] = "2.0";

    // prepare top level nodes
    if (!scene->cameras.empty()) js["cameras"] = json::array();
    if (!scene->textures.empty()) {
        js["textures"] = json::array();
        js["images"]   = json::array();
    }
    if (!scene->materials.empty()) js["materials"] = json::array();
    if (!scene->shapes.empty()) {
        js["meshes"]      = json::array();
        js["buffers"]     = json::array();
        js["bufferViews"] = json::array();
        js["accessors"]   = json::array();
    }
    if (!scene->instances.empty()) js["nodes"] = json::array();
    if (!scene->nodes.empty()) js["nodes"] = json::array();

    // convert cameras
    auto cmap = unordered_map<yocto_camera*, int>();
    for (auto camera : scene->cameras) {
        auto cjs    = json();
        cjs["name"] = camera->name;
        if (!camera->orthographic) {
            cjs["type"]                       = "perspective";
            cjs["perspective"]["aspectRatio"] = camera->film_size.x /
                                                camera->film_size.y;
            cjs["perspective"]["znear"] = 0.01f;
        } else {
            cjs["type"]                  = "orthographic";
            cjs["orthographic"]["xmag"]  = camera->film_size.x / 2;
            cjs["orthographic"]["ymag"]  = camera->film_size.y / 2;
            cjs["orthographic"]["znear"] = 0.01f;
        }
        cmap[camera] = (int)js["cameras"].size();
        js["cameras"].push_back(cjs);
    }

    // textures
    auto tmap = unordered_map<yocto_texture*, int>();
    for (auto& texture : scene->textures) {
        auto tjs = json(), ijs = json();
        tjs["source"] = (int)js["images"].size();
        ijs["uri"]    = texture->filename;
        js["images"].push_back(ijs);
        js["textures"].push_back(tjs);
        tmap[texture] = (int)js["textures"].size() - 1;
    }

    // material
    auto mmap = unordered_map<yocto_material*, int>();
    for (auto mat : scene->materials) {
        auto mjs           = json();
        mjs["name"]        = mat->name;
        mjs["doubleSided"] = mat->double_sided;
        if (mat->emission != zero3f) mjs["emissiveFactor"] = mat->emission;
        if (mat->emission_texture)
            mjs["emissiveTexture"]["index"] = tmap.at(mat->emission_texture);
        auto kd = vec4f{
            mat->diffuse.x, mat->diffuse.y, mat->diffuse.z, mat->opacity};
        if (mat->base_metallic) {
            auto mmjs               = json();
            mmjs["baseColorFactor"] = kd;
            mmjs["metallicFactor"]  = mat->specular.x;
            mmjs["roughnessFactor"] = mat->roughness;
            if (mat->diffuse_texture)
                mmjs["baseColorTexture"]["index"] = tmap.at(mat->diffuse_texture);
            if (mat->specular_texture)
                mmjs["metallicRoughnessTexture"]["index"] = tmap.at(
                    mat->specular_texture);
            mjs["pbrMetallicRoughness"] = mmjs;
        } else {
            auto mmjs                = json();
            mmjs["diffuseFactor"]    = kd;
            mmjs["specularFactor"]   = mat->specular;
            mmjs["glossinessFactor"] = 1 - mat->roughness;
            if (mat->diffuse_texture)
                mmjs["diffuseTexture"]["index"] = tmap.at(mat->diffuse_texture);
            if (mat->specular_texture)
                mmjs["specularGlossinessTexture"]["index"] = tmap.at(
                    mat->specular_texture);
            mjs["extensions"]["KHR_materials_pbrSpecularGlossiness"] = mmjs;
        }
        if (mat->normal_texture)
            mjs["normalTexture"]["index"] = tmap.at(mat->normal_texture);
        if (mat->occlusion_texture)
            mjs["occlusionTexture"]["index"] = tmap.at(mat->occlusion_texture);
        js["materials"].push_back(mjs);
        mmap[mat] = (int)js["materials"].size() - 1;
    }

    // determine shape materials
    auto shape_mats = unordered_map<yocto_shape*, int>();
    for (auto instance : scene->instances)
        if (instance->material)
            shape_mats[instance->shape] = mmap.at(instance->material);

    // shapes
    auto smap = unordered_map<yocto_shape*, int>();
    for (auto shape : scene->shapes) {
        auto mjs = json(), bjs = json(), pjs = json();
        auto bid          = js["buffers"].size();
        mjs["name"]       = shape->name;
        mjs["primitives"] = json::array();
        bjs["name"]       = shape->name;
        bjs["byteLength"] = 0;
        bjs["uri"]        = replace_extension(shape->filename, ".bin");
        auto mat_it       = shape_mats.find(shape);
        if (mat_it != shape_mats.end()) pjs["material"] = mat_it->second;
        auto add_accessor = [&js, &bjs, bid](
                                int count, string type, bool indices = false) {
            auto bytes = count * 4;
            if (type == "VEC2") bytes *= 2;
            if (type == "VEC3") bytes *= 3;
            if (type == "VEC4") bytes *= 4;
            auto ajs = json(), vjs = json();
            vjs["buffer"]        = bid;
            vjs["byteLength"]    = bytes;
            vjs["byteOffset"]    = bjs["byteLength"].get<int>();
            vjs["target"]        = (!indices) ? 34962 : 34963;
            bjs["byteLength"]    = bjs["byteLength"].get<int>() + bytes;
            ajs["bufferView"]    = (int)js["bufferViews"].size();
            ajs["byteOffset"]    = 0;
            ajs["componentType"] = (!indices) ? 5126 : 5125;
            ajs["count"]         = count;
            ajs["type"]          = type;
            js["accessors"].push_back(ajs);
            js["bufferViews"].push_back(vjs);
            return (int)js["accessors"].size() - 1;
        };
        auto nverts = (int)shape->positions.size();
        if (!shape->positions.empty())
            pjs["attributes"]["POSITION"] = add_accessor(nverts, "VEC3");
        if (!shape->normals.empty())
            pjs["attributes"]["NORMAL"] = add_accessor(nverts, "VEC3");
        if (!shape->texturecoords.empty())
            pjs["attributes"]["TEXCOORD_0"] = add_accessor(nverts, "VEC2");
        if (!shape->colors.empty())
            pjs["attributes"]["COLOR_0"] = add_accessor(nverts, "VEC4");
        if (!shape->radius.empty())
            pjs["attributes"]["RADIUS"] = add_accessor(nverts, "SCALAR");
        if (!shape->points.empty()) {
            pjs["indices"] = add_accessor(
                (int)shape->points.size(), "SCALAR", true);
            pjs["mode"] = 1;
        }
        if (!shape->lines.empty()) {
            pjs["indices"] = add_accessor(
                (int)shape->lines.size() * 2, "SCALAR", true);
            pjs["mode"] = 1;
        }
        if (!shape->triangles.empty()) {
            pjs["indices"] = add_accessor(
                (int)shape->triangles.size() * 3, "SCALAR", true);
            pjs["mode"] = 4;
        }
        if (!shape->quads.empty()) {
            auto triangles = convert_quads_to_triangles(shape->quads);
            pjs["indices"] = add_accessor(
                (int)triangles.size() * 3, "SCALAR", true);
            pjs["mode"] = 4;
        }
        mjs["primitives"].push_back(pjs);
        js["meshes"].push_back(mjs);
        js["buffers"].push_back(bjs);
        smap[shape] = (int)js["meshes"].size() - 1;
    }

    // nodes
    auto nmap = unordered_map<yocto_scene_node*, int>();
    for (auto& node : scene->nodes) {
        auto njs           = json();
        njs["name"]        = node->name;
        njs["matrix"]      = frame_to_mat(node->local);
        njs["translation"] = node->translation;
        njs["rotation"]    = node->rotation;
        njs["scale"]       = node->scale;
        if (node->camera) njs["camera"] = cmap.at(node->camera);
        if (node->instance) njs["mesh"] = smap.at(node->instance->shape);
        if (!node->children.empty()) {
            njs["children"] = json::array();
            for (auto& c : node->children)
                njs["children"].push_back(nmap.at(c));
        }
        js["nodes"].push_back(njs);
        nmap[node] = (int)js["nodes"].size() - 1;
    }

    // animations not supported yet
    if (!scene->animations.empty()) printf("animation not supported yet\n");

    // nodes from instances
    if (scene->nodes.empty()) {
        for (auto camera : scene->cameras) {
            auto njs      = json();
            njs["name"]   = camera->name;
            njs["camera"] = cmap.at(camera);
            njs["matrix"] = frame_to_mat(camera->frame);
            js["nodes"].push_back(njs);
        }
        for (auto instance : scene->instances) {
            auto njs      = json();
            njs["name"]   = instance->name;
            njs["mesh"]   = smap.at(instance->shape);
            njs["matrix"] = frame_to_mat(instance->frame);
            js["nodes"].push_back(njs);
        }
    }

    // done
    return true;
}

// save gltf mesh
bool save_gltf_mesh(const string& filename, const yocto_shape* shape) {
    auto fs = open(filename, "wb");
    if (!fs) return false;

    if (!write_values(fs, shape->positions)) return false;
    if (!write_values(fs, shape->normals)) return false;
    if (!write_values(fs, shape->texturecoords)) return false;
    if (!write_values(fs, shape->colors)) return false;
    if (!write_values(fs, shape->radius)) return false;
    if (!write_values(fs, shape->points)) return false;
    if (!write_values(fs, shape->lines)) return false;
    if (!write_values(fs, shape->triangles)) return false;
    if (!write_values(fs, convert_quads_to_triangles(shape->quads)))
        return false;

    return true;
}

// Save gltf json
bool save_gltf_scene(const string& filename, const yocto_scene* scene,
    bool save_textures, bool skip_missing) {
    // save json
    auto js = json();
    try {
        if (!scene_to_gltf(scene, js)) return false;
    } catch (...) {
        return false;
    }
    if (!save_json(filename, js)) return false;

    // meshes
    auto dirname = get_dirname(filename);
    for (auto& shape : scene->shapes) {
        if (shape->filename == "") continue;
        auto filename = normalize_path(dirname + "/" + shape->filename);
        filename      = replace_extension(filename, ".bin");
        if (!save_gltf_mesh(filename, shape)) {
            if (!skip_missing) return false;
        }
    }

    // save textures
    if (save_textures) {
        if (!save_scene_textures(scene, dirname, skip_missing)) return false;
    }

    // done
    return true;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF PBRT
// -----------------------------------------------------------------------------
namespace ygl {

// convert pbrt to json
bool pbrt_to_json(const string& filename, json& js) {
    auto split = [](const string& str) {
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
    };

    auto is_cmd = [](const vector<string>& tokens, int i) -> bool {
        auto& tok = tokens.at(i);
        return !(tok[0] == '[' || tok[0] == ']' || tok[0] == '\"' ||
                 tok[0] == '-' || tok[0] == '+' || tok[0] == '.' ||
                 std::isdigit(tok[0]));
    };
    auto is_number = [](const vector<string>& tokens, int i) -> bool {
        auto& tok = tokens.at(i);
        return tok[0] == '-' || tok[0] == '+' || tok[0] == '.' ||
               std::isdigit(tok[0]);
    };
    auto parse_string = [](const vector<string>& tokens, int& i) -> string {
        if (tokens[i][0] != '"') {
            log_error("string expected");
            return "";
        }
        auto tok = tokens[i++];
        tok      = tok.substr(1, tok.size() - 2);
        if (tok.find('|') != tok.npos) tok = tok.substr(tok.find('|') + 1);
        return tok;
    };
    auto parse_param = [&](const vector<string>& tokens, int& i,
                           json& js) -> void {
        auto list = false, first = true;
        while (i < tokens.size()) {
            if (is_cmd(tokens, i)) {
                break;
            } else if (tokens[i][0] == '[') {
                list = true;
                i++;
            } else if (tokens[i][0] == ']') {
                list = false;
                i++;
                break;
            } else if (tokens[i][0] == '"') {
                if (!first && !list) break;
                js.push_back(tokens[i].substr(1, tokens[i].size() - 2));
                i++;
                if (!list) break;
            } else {
                if (!first && !list) {
                    log_error("bad params");
                    break;
                }
                js.push_back(atof(tokens[i].c_str()));
                i++;
                if (!list) break;
            }
        }
    };
    auto parse_param_list = [&](const vector<string>& tokens, int& i,
                                json& js) -> void {
        while (i < tokens.size()) {
            if (is_cmd(tokens, i)) break;
            auto name = parse_string(tokens, i);
            js[name]  = json::array();
            parse_param(tokens, i, js.at(name));
            if (js.at(name).size() == 1) {
                js.at(name) = js.at(name).at(0);
            }
        }
    };
    auto parse_param_numbers = [&](const vector<string>& tokens, int& i,
                                   json& js) -> void {
        js["values"] = json::array();
        if (tokens[i][0] == '[') i++;
        while (is_number(tokens, i)) {
            js.at("values").push_back((float)atof(tokens[i++].c_str()));
        }
        if (tokens[i][0] == ']') i++;
    };

    auto fs = open(filename, "rt");
    if (!fs) return false;

    auto pbrt = ""s;
    auto line = ""s;
    while (read_line(fs, line)) {
        if (line.find('#') == line.npos)
            pbrt += line + "\n";
        else
            pbrt += line.substr(0, line.find('#')) + "\n";
    }

    auto re = std::regex("\"(\\w+)\\s+(\\w+)\"");
    pbrt    = std::regex_replace(pbrt, re, "\"$1|$2\"");
    pbrt    = std::regex_replace(pbrt, std::regex("\\["), " [ ");
    pbrt    = std::regex_replace(pbrt, std::regex("\\]"), " ] ");
    js      = json::array();

    auto tokens = split(pbrt);
    auto i      = 0;
    while (i < tokens.size()) {
        if (!is_cmd(tokens, i)) {
            runtime_error("command expected");
            break;
        }
        auto& tok   = tokens[i++];
        auto  jcmd  = json::object();
        jcmd["cmd"] = tok;
        if (tok == "Transform" || tok == "LookAt" || tok == "Scale" ||
            tok == "Rotate" || tok == "Translate" || tok == "ConcatTransform") {
            parse_param_numbers(tokens, i, jcmd);
        } else if (tok == "Integrator" || tok == "Sampler" ||
                   tok == "PixelFilter" || tok == "Film" || tok == "Camera" ||
                   tok == "Shape" || tok == "AreaLightSource" ||
                   tok == "LightSource") {
            jcmd["type"] = parse_string(tokens, i);
            parse_param_list(tokens, i, jcmd);
        } else if (tok == "Texture") {
            jcmd["name"]       = parse_string(tokens, i);
            jcmd["value_type"] = parse_string(tokens, i);
            jcmd["type"]       = parse_string(tokens, i);
            parse_param_list(tokens, i, jcmd);
        } else if (tok == "MakeNamedMaterial") {
            jcmd["name"] = parse_string(tokens, i);
            parse_param_list(tokens, i, jcmd);
        } else if (tok == "Material") {
            jcmd["type"] = parse_string(tokens, i);
            parse_param_list(tokens, i, jcmd);
        } else if (tok == "NamedMaterial" || tok == "ObjectBegin" ||
                   tok == "ObjectInstance") {
            jcmd["name"] = parse_string(tokens, i);
        } else if (tok == "WorldBegin" || tok == "AttributeBegin" ||
                   tok == "TransformBegin" || tok == "WorldEnd" ||
                   tok == "AttributeEnd" || tok == "TransformEnd" ||
                   tok == "ObjectEnd" || tok == "ReverseOrientation") {
        } else {
            log_error("unsupported command {}", tok);
        }
        js.push_back(jcmd);
    }
    // auto fstr = std::fstream(filename + ".json");
    // fstr << js;
    return true;
}

// load pbrt scenes
yocto_scene* load_pbrt_scene(
    const string& filename, bool load_textures, bool skip_missing) {
    auto scope = log_trace_scoped("loading scene {}", filename);
    // convert to json
    auto js = json();
    try {
        if (!pbrt_to_json(filename, js)) return nullptr;
    } catch (...) {
        return nullptr;
    }

    auto dirname_ = get_dirname(filename);

    struct stack_item {
        frame3f         frame     = identity_frame3f;
        yocto_material* mat       = nullptr;
        yocto_material* light_mat = nullptr;
        float           focus = 1, aspect = 1;
        bool            reverse = false;
    };

    // parse
    auto scene = make_unique<yocto_scene>();
    auto stack = vector<stack_item>();
    stack.push_back(stack_item());

    auto txt_map = map<string, yocto_texture*>();
    auto mat_map = map<string, yocto_material*>();
    auto mid     = 0;

    auto get_vec3f = [](const json& js) -> vec3f {
        if (js.is_number())
            return {js.get<float>(), js.get<float>(), js.get<float>()};
        if (js.is_array() && js.size() == 1)
            return {js.at(0).get<float>(), js.at(0).get<float>(),
                js.at(0).get<float>()};
        if (js.is_array() && js.size() == 3)
            return {js.at(0).get<float>(), js.at(1).get<float>(),
                js.at(2).get<float>()};
        printf("cannot handle vec3f\n");
        return zero3f;
    };

    auto get_vec4f = [](const json& js) -> vec4f {
        if (js.is_number())
            return {js.get<float>(), js.get<float>(), js.get<float>(),
                js.get<float>()};
        if (js.is_array() && js.size() == 4)
            return {js.at(0).get<float>(), js.at(1).get<float>(),
                js.at(2).get<float>(), js.at(3).get<float>()};
        printf("cannot handle vec4f\n");
        return zero4f;
    };

    auto get_mat4f = [](const json& js) -> frame3f {
        if (!js.is_array() || js.size() != 16) {
            printf("cannot handle vec4f\n");
            return identity_frame3f;
        }
        float m[16] = {0};
        for (auto i = 0; i < 16; i++) m[i] = js.at(i).get<float>();
        return {{m[0], m[1], m[2]}, {m[4], m[5], m[6]}, {m[8], m[9], m[10]},
            {m[12], m[13], m[14]}};
    };

    auto get_mat3f = [](const json& js) -> frame3f {
        if (!js.is_array() || js.size() != 9) {
            printf("cannot handle mat3f\n");
            return identity_frame3f;
        }
        auto m = identity_frame3f;
        for (auto i = 0; i < 9; i++) (&m.x.x)[i] = js.at(i).get<float>();
        return m;
    };

    auto get_vector_vec3i = [](const json& js) -> vector<vec3i> {
        if (!js.is_array() || js.size() % 3) {
            printf("cannot handle vector<vec3f>\n");
            return {};
        }
        auto vals = vector<vec3i>(js.size() / 3);
        for (auto i = 0; i < vals.size(); i++) {
            vals[i].x = (int)std::round(js.at(i * 3 + 0).get<float>());
            vals[i].y = (int)std::round(js.at(i * 3 + 1).get<float>());
            vals[i].z = (int)std::round(js.at(i * 3 + 2).get<float>());
        }
        return vals;
    };

    auto get_vector_vec3f = [](const json& js) -> vector<vec3f> {
        if (!js.is_array() || js.size() % 3) {
            printf("cannot handle vector<vec3f>\n");
            return {};
        }
        auto vals = vector<vec3f>(js.size() / 3);
        for (auto i = 0; i < vals.size(); i++) {
            vals[i].x = js.at(i * 3 + 0).get<float>();
            vals[i].y = js.at(i * 3 + 1).get<float>();
            vals[i].z = js.at(i * 3 + 2).get<float>();
        }
        return vals;
    };

    auto get_vector_vec2f = [](const json& js) -> vector<vec2f> {
        if (!js.is_array() || js.size() % 2) {
            printf("cannot handle vector<vec3f>\n");
            return {};
        }
        auto vals = vector<vec2f>(js.size() / 2);
        for (auto i = 0; i < vals.size(); i++) {
            vals[i].x = js.at(i * 2 + 0).get<float>();
            vals[i].y = js.at(i * 2 + 1).get<float>();
        }
        return vals;
    };

    auto get_scaled_texture =
        [&txt_map, &get_vec3f](const json& js) -> tuple<vec3f, yocto_texture*> {
        if (js.is_string()) return {{1, 1, 1}, txt_map.at(js.get<string>())};
        return {get_vec3f(js), nullptr};
    };

    auto use_hierarchy = false;

    map<string, vector<yocto_instance*>> objects;
    for (auto& jcmd : js) {
        auto cmd = jcmd.at("cmd").get<string>();
        if (cmd == "ObjectInstance") {
            use_hierarchy = true;
            break;
        }
    }

    auto lid = 0, sid = 0, cid = 0;
    auto cur_object = ""s;
    for (auto& jcmd : js) {
        auto cmd = jcmd.at("cmd").get<string>();
        if (cmd == "Integrator" || cmd == "Sampler" || cmd == "PixelFilter") {
        } else if (cmd == "Transform") {
            stack.back().frame = get_mat4f(jcmd.at("values"));
        } else if (cmd == " ConcatTransform") {
            stack.back().frame = stack.back().frame * get_mat4f(jcmd.at("value"
                                                                        "s"));
        } else if (cmd == "Scale") {
            auto v             = get_vec3f(jcmd.at("values"));
            stack.back().frame = stack.back().frame * scaling_frame(v);
        } else if (cmd == "Translate") {
            auto v             = get_vec3f(jcmd.at("values"));
            stack.back().frame = stack.back().frame * translation_frame(v);
        } else if (cmd == "Rotate") {
            auto v             = get_vec4f(jcmd.at("values"));
            stack.back().frame = stack.back().frame *
                                 rotation_frame(
                                     vec3f{v.y, v.z, v.w}, v.x * pif / 180);
        } else if (cmd == "LookAt") {
            auto m             = get_mat3f(jcmd.at("values"));
            stack.back().frame = stack.back().frame *
                                 inverse(lookat_frame(m.x, m.y, m.z, true));
            stack.back().focus = length(m.x - m.y);
        } else if (cmd == "ReverseOrientation") {
            stack.back().reverse = !stack.back().reverse;
        } else if (cmd == "Film") {
            stack.back().aspect = jcmd.at("xresolution").get<float>() /
                                  jcmd.at("yresolution").get<float>();
        } else if (cmd == "Camera") {
            auto camera            = new yocto_camera();
            camera->name           = "camera" + std::to_string(cid++);
            camera->frame          = inverse(stack.back().frame);
            camera->frame.z        = -camera->frame.z;
            camera->focus_distance = stack.back().focus;
            auto aspect            = stack.back().aspect;
            auto fovy              = 1.0f;
            auto type              = jcmd.at("type").get<string>();
            if (type == "perspective") {
                fovy = jcmd.at("fov").get<float>() * pif / 180;
            } else {
                printf("%s camera not supported\n", type.c_str());
            }
            set_camera_fovy(camera, fovy, aspect);
            scene->cameras.push_back(camera);
        } else if (cmd == "Texture") {
            auto found = false;
            auto name  = jcmd.at("name").get<string>();
            for (auto& texture : scene->textures) {
                if (texture->name == name) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                auto texture = new yocto_texture();
                scene->textures.push_back(texture);
                texture->name          = jcmd.at("name").get<string>();
                txt_map[texture->name] = texture;
                auto type              = jcmd.at("type").get<string>();
                if (type == "imagemap") {
                    texture->filename = jcmd.at("filename").get<string>();
                    if (get_extension(texture->filename) == "pfm")
                        texture->filename = replace_extension(texture->filename,
                            ".hd"
                            "r");
                } else {
                    printf("%s texture not supported\n", type.c_str());
                }
            }
        } else if (cmd == "MakeNamedMaterial" || cmd == "Material") {
            auto found = false;
            if (cmd == "MakeNamedMaterial") {
                auto name = jcmd.at("name").get<string>();
                for (auto mat : scene->materials) {
                    if (mat->name == name) {
                        found = true;
                        break;
                    }
                }
            }
            if (!found) {
                auto mat = new yocto_material();
                scene->materials.push_back(mat);
                if (cmd == "Material") {
                    mat->name        = "unnamed_mat" + std::to_string(mid++);
                    stack.back().mat = mat;
                } else {
                    mat->name          = jcmd.at("name").get<string>();
                    mat_map[mat->name] = mat;
                }
                auto type = "uber"s;
                if (jcmd.count("type")) type = jcmd.at("type").get<string>();
                if (type == "uber") {
                    if (jcmd.count("Kd"))
                        tie(mat->diffuse, mat->diffuse_texture) = get_scaled_texture(
                            jcmd.at("Kd"));
                    if (jcmd.count("Ks"))
                        tie(mat->specular, mat->specular_texture) = get_scaled_texture(
                            jcmd.at("Ks"));
                    if (jcmd.count("Kt"))
                        tie(mat->transmission, mat->transmission_texture) = get_scaled_texture(
                            jcmd.at("Kt"));
                    if (jcmd.count("opacity")) {
                        auto op         = vec3f{0, 0, 0};
                        auto op_txt     = (yocto_texture*)nullptr;
                        tie(op, op_txt) = get_scaled_texture(jcmd.at("opacity"));
                        mat->opacity    = (op.x + op.y + op.z) / 3;
                        mat->opacity_texture = op_txt;
                    }
                    mat->roughness = 0;
                } else if (type == "matte") {
                    mat->diffuse = {1, 1, 1};
                    if (jcmd.count("Kd"))
                        tie(mat->diffuse, mat->diffuse_texture) = get_scaled_texture(
                            jcmd.at("Kd"));
                    mat->roughness = 1;
                } else if (type == "mirror") {
                    mat->diffuse   = {0, 0, 0};
                    mat->specular  = {1, 1, 1};
                    mat->roughness = 0;
                } else if (type == "metal") {
                    auto eta       = get_vec3f(jcmd.at("eta"));
                    auto k         = get_vec3f(jcmd.at("k"));
                    mat->specular  = fresnel_metal(1, eta, k);
                    mat->roughness = 0;
                } else if (type == "substrate") {
                    if (jcmd.count("Kd"))
                        tie(mat->diffuse, mat->diffuse_texture) = get_scaled_texture(
                            jcmd.at("Kd"));
                    mat->specular = {0.04f, 0.04f, 0.04f};
                    if (jcmd.count("Ks"))
                        tie(mat->specular, mat->specular_texture) = get_scaled_texture(
                            jcmd.at("Ks"));
                    mat->roughness = 0;
                } else if (type == "glass") {
                    mat->specular     = {0.04f, 0.04f, 0.04f};
                    mat->transmission = {1, 1, 1};
                    if (jcmd.count("Ks"))
                        tie(mat->specular, mat->specular_texture) = get_scaled_texture(
                            jcmd.at("Ks"));
                    if (jcmd.count("Kt"))
                        tie(mat->transmission, mat->transmission_texture) = get_scaled_texture(
                            jcmd.at("Kt"));
                    mat->roughness = 0;
                } else if (type == "mix") {
                    printf("mix material not properly supported\n");
                    if (jcmd.count("namedmaterial1")) {
                        auto mat1 = jcmd.at("namedmaterial1").get<string>();
                        auto saved_name = mat->name;
                        *mat            = *mat_map.at(mat1);
                        mat->name       = saved_name;
                    } else {
                        printf("mix material missing front material\n");
                    }
                } else {
                    mat->diffuse = {1, 0, 0};
                    printf("%s material not supported\n", type.c_str());
                }
                if (jcmd.count("uroughness")) {
                    auto remap = js.count("remaproughness") &&
                                 js.at("remaproughness").get<bool>();
                    if (jcmd.count("uroughness"))
                        mat->roughness = jcmd.at("uroughness").get<float>();
                    // if (!remap) mat->rs = mat->rs * mat->rs;
                    if (remap) printf("remap roughness not supported\n");
                }
                if (jcmd.count("roughness")) {
                    auto remap = js.count("remaproughness") &&
                                 js.at("remaproughness").get<bool>();
                    if (jcmd.count("roughness"))
                        mat->roughness = jcmd.at("roughness").get<float>();
                    // if (!remap) mat->rs = mat->rs * mat->rs;
                    if (remap) printf("remap roughness not supported\n");
                }
                if (stack.back().light_mat) {
                    mat->emission         = stack.back().light_mat->emission;
                    mat->emission_texture = stack.back().light_mat->emission_texture;
                }
            }
        } else if (cmd == "NamedMaterial") {
            stack.back().mat = mat_map.at(jcmd.at("name").get<string>());
            if (stack.back().light_mat) {
                auto mat = new yocto_material(*stack.back().mat);
                mat->name += "_" + std::to_string(lid++);
                mat->emission         = stack.back().light_mat->emission;
                mat->emission_texture = stack.back().light_mat->emission_texture;
                scene->materials.push_back(mat);
                stack.back().mat = mat;
            }
        } else if (cmd == "Shape") {
            auto shape = new yocto_shape();
            auto type  = jcmd.at("type").get<string>();
            if (type == "plymesh") {
                auto filename   = jcmd.at("filename").get<string>();
                shape->name     = get_filename(filename);
                shape->filename = filename;
                if (!load_ply_mesh(dirname_ + "/" + filename, shape->points,
                        shape->lines, shape->triangles, shape->quads,
                        shape->positions, shape->normals, shape->texturecoords,
                        shape->colors, shape->radius, false))
                    return nullptr;
            } else if (type == "trianglemesh") {
                shape->name     = "mesh" + std::to_string(sid++);
                shape->filename = "models/" + shape->name + ".ply";
                if (jcmd.count("indices"))
                    shape->triangles = get_vector_vec3i(jcmd.at("indices"));
                if (jcmd.count("P"))
                    shape->positions = get_vector_vec3f(jcmd.at("P"));
                if (jcmd.count("N"))
                    shape->normals = get_vector_vec3f(jcmd.at("N"));
                if (jcmd.count("uv"))
                    shape->texturecoords = get_vector_vec2f(jcmd.at("uv"));
            } else if (type == "sphere") {
                shape->name     = "sphere" + std::to_string(sid++);
                shape->filename = "models/" + shape->name + ".ply";
                auto radius     = 1.0f;
                if (jcmd.count("radius"))
                    radius = jcmd.at("radius").get<float>();
                auto sshp = make_sphere_shape({64, 32}, 2 * radius, {1, 1}, true);
                shape->positions     = sshp.positions;
                shape->normals       = sshp.normals;
                shape->texturecoords = sshp.texturecoords;
                shape->triangles     = sshp.triangles;
            } else if (type == "disk") {
                shape->name     = "disk" + std::to_string(sid++);
                shape->filename = "models/" + shape->name + ".ply";
                auto radius     = 1.0f;
                if (jcmd.count("radius"))
                    radius = jcmd.at("radius").get<float>();
                auto sshp = make_disk_shape({32, 16}, 2 * radius, {1, 1}, true);
                shape->positions     = sshp.positions;
                shape->normals       = sshp.normals;
                shape->texturecoords = sshp.texturecoords;
                shape->triangles     = sshp.triangles;
            } else {
                printf("%s shape not supported\n", type.c_str());
            }
            auto frame = stack.back().frame;
            auto scl = vec3f{length(frame.x), length(frame.y), length(frame.z)};
            for (auto& p : shape->positions) p *= scl;
            frame = {normalize(frame.x), normalize(frame.y), normalize(frame.z),
                frame.o};
            if (stack.back().reverse) {
                for (auto& t : shape->triangles) swap(t.y, t.z);
            }
            scene->shapes.push_back(shape);
            auto instance      = new yocto_instance();
            instance->name     = shape->name;
            instance->frame    = frame;
            instance->shape    = shape;
            instance->material = stack.back().mat;
            if (cur_object != "") {
                objects[cur_object].push_back(instance);
            } else {
                scene->instances.push_back(instance);
            }
        } else if (cmd == "ObjectInstance") {
            static auto instances = map<string, int>();
            auto        name      = jcmd.at("name").get<string>();
            auto&       object    = objects.at(name);
            for (auto shape : object) {
                instances[shape->name] += 1;
                auto instance  = new yocto_instance();
                instance->name = shape->name + "_ist" +
                                 std::to_string(instances[shape->name]);
                instance->frame = stack.back().frame * shape->frame;
                instance->shape = shape->shape;
                scene->instances.push_back(instance);
            }
        } else if (cmd == "AreaLightSource") {
            auto type = jcmd.at("type").get<string>();
            if (type == "diffuse") {
                auto lmat              = new yocto_material();
                lmat->emission         = get_vec3f(jcmd.at("L"));
                stack.back().light_mat = lmat;
            } else {
                printf("%s area light not supported\n", type.c_str());
            }
        } else if (cmd == "LightSource") {
            auto type = jcmd.at("type").get<string>();
            if (type == "infinite") {
                auto environment  = new yocto_environment();
                environment->name = "environment" + std::to_string(lid++);
                // environment->frame =
                // frame3f{{1,0,0},{0,0,-1},{0,-1,0},{0,0,0}}
                // * stack.back().frame;
                environment->frame = stack.back().frame *
                                     frame3f{{0, 0, 1}, {0, 1, 0}, {1, 0, 0},
                                         {0, 0, 0}};
                environment->emission = {1, 1, 1};
                if (jcmd.count("scale"))
                    environment->emission *= get_vec3f(jcmd.at("scale"));
                if (jcmd.count("mapname")) {
                    auto texture      = new yocto_texture();
                    texture->filename = jcmd.at("mapname").get<string>();
                    texture->name     = environment->name;
                    scene->textures.push_back(texture);
                    environment->emission_texture = texture;
                }
                scene->environments.push_back(environment);
            } else if (type == "distant") {
                auto distant_dist = 100;
                auto shape        = new yocto_shape();
                shape->name       = "distant" + std::to_string(lid++);
                auto from = vec3f{0, 0, 0}, to = vec3f{0, 0, 0};
                if (jcmd.count("from")) from = get_vec3f(jcmd.at("from"));
                if (jcmd.count("to")) to = get_vec3f(jcmd.at("to"));
                auto dir  = normalize(from - to);
                auto size = distant_dist * sin(5 * pif / 180);
                auto sshp = make_quad_shape({1, 1}, {size, size}, {1, 1}, true);
                shape->positions     = sshp.positions;
                shape->normals       = sshp.normals;
                shape->texturecoords = sshp.texturecoords;
                shape->triangles     = sshp.triangles;
                scene->shapes.push_back(shape);
                auto mat      = new yocto_material();
                mat->name     = shape->name;
                mat->emission = {1, 1, 1};
                if (jcmd.count("L")) mat->emission *= get_vec3f(jcmd.at("L"));
                if (jcmd.count("scale"))
                    mat->emission *= get_vec3f(jcmd.at("scale"));
                mat->emission *= (distant_dist * distant_dist) / (size * size);
                scene->materials.push_back(mat);
                auto instance      = new yocto_instance();
                instance->name     = shape->name;
                instance->shape    = shape;
                instance->material = mat;
                instance->frame    = stack.back().frame *
                                  lookat_frame(dir * distant_dist, zero3f,
                                      {0, 1, 0}, true);
                scene->instances.push_back(instance);
                printf("%s light not properly supported\n", type.c_str());
            } else {
                printf("%s light not supported\n", type.c_str());
            }
        } else if (cmd == "WorldBegin") {
            stack.push_back(stack_item());
        } else if (cmd == "AttributeBegin") {
            stack.push_back(stack.back());
        } else if (cmd == "ObjectBegin") {
            auto name     = jcmd.at("name").get<string>();
            cur_object    = name;
            objects[name] = {};
        } else if (cmd == "ObjectEnd") {
            cur_object = "";
        } else if (cmd == "TransformBegin") {
            stack.push_back(stack.back());
        } else if (cmd == "WorldEnd" || cmd == "AttributeEnd" ||
                   cmd == "TransformEnd") {
            stack.pop_back();
        } else {
            printf("%s command not supported\n", cmd.c_str());
        }
    }
    if (use_hierarchy) {
        for (auto camera : scene->cameras) {
            auto node    = new yocto_scene_node();
            node->name   = camera->name;
            node->local  = camera->frame;
            node->camera = camera;
            scene->nodes.insert(scene->nodes.begin(), node);
        }
        for (auto environment : scene->environments) {
            auto node         = new yocto_scene_node();
            node->name        = environment->name;
            node->local       = environment->frame;
            node->environment = environment;
            scene->nodes.push_back(node);
        }
    }

    // load textures
    auto dirname = get_dirname(filename);
    if (load_textures) {
        if (!load_scene_textures(scene.get(), dirname, skip_missing, false))
            return nullptr;
    }

    // fix scene
    scene->name = get_filename(filename);
    add_missing_cameras(scene.get());
    add_missing_materials(scene.get());
    add_missing_names(scene.get());
    update_transforms(scene.get());

    return scene.release();
}

// Convert a scene to pbrt format
bool save_pbrt(const string& filename, const yocto_scene* scene) {
    auto fs = open(filename, "wt");
    if (!fs) return false;

#if 0
WorldBegin

#uniform blue - ish illumination from all directions
LightSource "infinite" "rgb L" [.4 .45 .5]

#approximate the sun
LightSource "distant"  "point from" [ -30 40  100 ]
   "blackbody L" [3000 1.5]

AttributeBegin
  Material "glass"
  Shape "sphere" "float radius" 1
AttributeEnd

AttributeBegin
  Texture "checks" "spectrum" "checkerboard"
          "float uscale" [8] "float vscale" [8]
          "rgb tex1" [.1 .1 .1] "rgb tex2" [.8 .8 .8]
  Material "matte" "texture Kd" "checks"
  Translate 0 0 -1
  Shape "trianglemesh"
      "integer indices" [0 1 2 0 2 3]
      "point P" [ -20 -20 0   20 -20 0   20 20 0   -20 20 0 ]
      "float st" [ 0 0   1 0    1 1   0 1 ]
AttributeEnd

WorldEnd
#endif

    // convert camera and settings
    auto camera = scene->cameras.front();
    auto from   = camera->frame.o;
    auto to     = camera->frame.o - camera->frame.z;
    auto up     = camera->frame.y;
    print(fs, "LookAt {} {} {}\n", from, to, up);
    print(fs, "Camera \"perspective\" \"float fov\" {}\n",
        evaluate_camera_fovy(camera) * 180 / pif);

    // save renderer
    print(fs, "Sampler \"random\" \"integer pixelsamples\" [64]\n");
    // fprintf(f, "Sampler \"sobol\" \"interger pixelsamples\" [64]\n");
    print(fs, "Integrator \"path\"\n");
    print(fs,
        "Film \"image\" \"string filename\" [\"{}\"] "
        "\"integer xresolution\" [{}] \"integer yresolution\" [{}]\n",
        replace_extension(filename, "exr"), evaluate_image_size(camera, 512).x,
        evaluate_image_size(camera, 512).y);

    // start world
    print(fs, "WorldBegin\n");

    // convert textures
    for (auto texture : scene->textures) {
        print(fs,
            "Texture \"{}\" \"spectrum\" \"imagemap\" "
            "\"string filename\" [\"{}\"]\n",
            texture->name, texture->filename);
    }

    // convert materials
    for (auto mat : scene->materials) {
        print(fs, "MakeNamedMaterial \"{}\" ", mat->name);
        print(fs, "\"string type\" \"{}\" ", "uber");
        if (mat->diffuse_texture)
            print(fs, "\"texture Kd\" [\"{}\"] ", mat->diffuse_texture->name);
        else
            print(fs, "\"rgb Kd\" [{}] ", mat->diffuse);
        if (mat->specular_texture)
            print(fs, "\"texture Ks\" [\"{}\"] ", mat->specular_texture->name);
        else
            print(fs, "\"rgb Ks\" [{}] ", mat->specular);
        print(fs, "\"float roughness\" [{}] ", mat->roughness);
        print(fs, "\n");
    }

    // convert instances
    for (auto instance : scene->instances) {
        print(fs, "AttributeBegin\n");
        print(fs, "TransformBegin\n");
        print(fs, "ConcatTransform [{}]\n", frame_to_mat(instance->frame));
        if (instance->material->emission != zero3f)
            print(fs, "AreaLightSource \"diffuse\" \"rgb L\" [ {} ]\n",
                instance->material->emission);
        print(fs, "NamedMaterial \"{}\"\n", instance->material->name);
        print(fs, "Shape \"plymesh\" \"string filename\" [\"{}\"]\n",
            instance->shape->filename.c_str());
        print(fs, "TransformEnd\n");
        print(fs, "AttributeEnd\n");
    }

    // end world
    print(fs, "WorldEnd\n");

    // done
    return true;
}

// Save a pbrt scene
bool save_pbrt_scene(const string& filename, const yocto_scene* scene,
    bool save_textures, bool skip_missing) {
    // save json
    if (!save_pbrt(filename, scene)) return false;

    // save meshes
    auto dirname = get_dirname(filename);
    for (auto& shape : scene->shapes) {
        if (shape->filename == "") continue;
        auto filename = normalize_path(dirname + "/" + shape->filename);
        if (!save_mesh(filename, shape->points, shape->lines, shape->triangles,
                shape->quads, shape->positions, shape->normals,
                shape->texturecoords, shape->colors, shape->radius)) {
            if (!skip_missing) return false;
        }
    }

    // skip textures
    if (save_textures) {
        if (!save_scene_textures(scene, dirname, skip_missing)) return false;
    }

    // done
    return true;
}

// Attempt to fix pbrt z-up.
void pbrt_flipyz_scene(const yocto_scene* scene) {
    // flip meshes
    for (auto shape : scene->shapes) {
        for (auto& p : shape->positions) swap(p.y, p.z);
        for (auto& n : shape->normals) swap(n.y, n.z);
    }
    for (auto instance : scene->instances) {
        instance->frame = instance->frame *
                          frame3f{{1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 0, 0}};
    }
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF BINARY SCENE FORMAT
// -----------------------------------------------------------------------------
namespace ygl {

// serialize_bin( ) can both save/load data to/from a binary file. The behaviour
// is set by the boolean 'save'. serialize_bin(var, file, true) : writes var as
// binary into file serialize_bin(var, file, false): read file as binary and set
// var

// Serialize type or struct with no allocated resource
template <typename T>
bool serialize_bin_value(T& val, file_stream& fs, bool save) {
    if (save) {
        return write_value(fs, val);
    } else {
        return read_value(fs, val);
    }
}

// Serialize vector
template <typename T>
bool serialize_bin_value(vector<T>& vec, file_stream& fs, bool save) {
    if (save) {
        auto count = (size_t)vec.size();
        if (!write_value(fs, count)) return false;
        if (!write_values(fs, vec)) return false;
        return true;
    } else {
        auto count = (size_t)0;
        if (!read_value(fs, count)) return false;
        vec = vector<T>(count);
        if (!read_values(fs, vec)) return false;
        return true;
    }
}

// Serialize string
bool serialize_bin_value(string& str, file_stream& fs, bool save) {
    if (save) {
        auto count = (size_t)str.size();
        if (!write_value(fs, count)) return false;
        auto vec = vector<char>(str.begin(), str.end());
        if (!write_values(fs, vec)) return false;
        return true;
    } else {
        auto count = (size_t)0;
        if (!read_value(fs, count)) return false;
        auto vec = vector<char>(count);
        if (!read_values(fs, vec)) return false;
        str = {vec.begin(), vec.end()};
        return true;
    }
}

// Serialize image
template <typename T>
bool serialize_bin_value(image<T>& img, file_stream& fs, bool save) {
    if (save) {
        if (!write_value(fs, img.width)) return false;
        if (!write_value(fs, img.height)) return false;
        if (!write_values(fs, img.pixels)) return false;
        return true;
    } else {
        auto width = 0, height = 0;
        if (!read_value(fs, width)) return false;
        if (!read_value(fs, height)) return false;
        img = image<T>{width, height};
        if (!read_values(fs, img.pixels)) return false;
        return true;
    }
}

// Serialize image
template <typename T>
bool serialize_bin_value(volume<T>& vol, file_stream& fs, bool save) {
    if (save) {
        if (!write_value(fs, vol.width)) return false;
        if (!write_value(fs, vol.height)) return false;
        if (!write_value(fs, vol.depth)) return false;
        if (!write_values(fs, vol.voxels)) return false;
        return true;
    } else {
        auto width = 0, height = 0, depth = 0;
        if (!read_value(fs, width)) return false;
        if (!read_value(fs, height)) return false;
        if (!read_value(fs, depth)) return false;
        vol = volume<T>{width, height, depth};
        if (!read_values(fs, vol.voxels)) return false;
        return true;
    }
}

// Serialize vector of pointers
template <typename T>
bool serialize_bin_object(vector<T*>& vec, file_stream& fs, bool save) {
    if (save) {
        auto count = (size_t)vec.size();
        if (!serialize_bin_value(count, fs, true)) return false;
        for (auto i = 0; i < vec.size(); ++i) {
            if (!serialize_bin_object(vec[i], fs, true)) return false;
        }
        return true;
    } else {
        auto count = (size_t)0;
        if (!serialize_bin_value(count, fs, false)) return false;
        vec = vector<T*>(count);
        for (auto i = 0; i < vec.size(); ++i) {
            vec[i] = new T();
            if (!serialize_bin_object(vec[i], fs, false)) return false;
        }
        return true;
    }
}

// Serialize vector of pointers
template <typename T>
bool serialize_bin_object(
    vector<T*>& vec, const yocto_scene* scene, file_stream& fs, bool save) {
    if (save) {
        auto count = (size_t)vec.size();
        if (!serialize_bin_value(count, fs, true)) return false;
        for (auto i = 0; i < vec.size(); ++i) {
            if (!serialize_bin_object(vec[i], scene, fs, true)) return false;
        }
        return true;
    } else {
        auto count = (size_t)0;
        if (!serialize_bin_value(count, fs, false)) return false;
        vec = vector<T*>(count);
        for (auto i = 0; i < vec.size(); ++i) {
            vec[i] = new T();
            if (!serialize_bin_object(vec[i], scene, fs, false)) return false;
        }
        return true;
    }
}

// Serialize a pointer. It is saved as an integer index (handle) of the array of
// pointers vec. On loading, the handle is converted back into a pointer.
template <typename T>
bool serialize_bin_handle(
    T*& val, const vector<T*>& vec, file_stream& fs, bool save) {
    if (save) {
        auto handle = -1;
        for (auto i = 0; i < vec.size(); ++i)
            if (vec[i] == val) {
                handle = i;
                break;
            }
        if (!serialize_bin_value(handle, fs, true)) return false;
        return true;
    } else {
        auto handle = -1;
        if (!serialize_bin_value(handle, fs, false)) return false;
        val = (handle == -1) ? nullptr : vec[handle];
        return true;
    }
}

// Serialize a pointer. It is saved as an integer index (handle) of the array of
// pointers vec. On loading, the handle is converted back into a pointer.
template <typename T>
bool serialize_bin_handle(
    vector<T*>& vals, const vector<T*>& vec_, file_stream& fs, bool save) {
    if (save) {
        auto count = (size_t)vals.size();
        if (!serialize_bin_value(count, fs, true)) return false;
        for (auto i = 0; i < vals.size(); ++i) {
            if (!serialize_bin_handle(vals[i], vec_, fs, true)) return false;
        }
    } else {
        auto count = (size_t)0;
        if (!serialize_bin_value(count, fs, false)) return false;
        vals = vector<T*>(count);
        for (auto i = 0; i < vals.size(); ++i) {
            if (!serialize_bin_handle(vals[i], vec_, fs, false)) return false;
        }
    }
}

// Serialize yocto types. This is mostly boiler plate code.
bool serialize_bin_object(yocto_camera* camera, file_stream& fs, bool save) {
    if (!serialize_bin_value(camera->name, fs, save)) return false;
    if (!serialize_bin_value(camera->frame, fs, save)) return false;
    if (!serialize_bin_value(camera->orthographic, fs, save)) return false;
    if (!serialize_bin_value(camera->film_size, fs, save)) return false;
    if (!serialize_bin_value(camera->focal_length, fs, save)) return false;
    if (!serialize_bin_value(camera->focus_distance, fs, save)) return false;
    if (!serialize_bin_value(camera->lens_aperture, fs, save)) return false;
    return true;
}

bool serialize_bin_object(bvh_tree* bvh, file_stream& fs, bool save) {
    if (!serialize_bin_value(bvh->positions, fs, save)) return false;
    if (!serialize_bin_value(bvh->radius, fs, save)) return false;
    if (!serialize_bin_value(bvh->points, fs, save)) return false;
    if (!serialize_bin_value(bvh->lines, fs, save)) return false;
    if (!serialize_bin_value(bvh->triangles, fs, save)) return false;
    if (!serialize_bin_value(bvh->quads, fs, save)) return false;
    if (!serialize_bin_value(bvh->nodes, fs, save)) return false;
    if (!serialize_bin_value(bvh->instances, fs, save)) return false;
    if (!serialize_bin_object(bvh->shape_bvhs, fs, save)) return false;
    if (!serialize_bin_value(bvh->nodes, fs, save)) return false;
    return true;
}

bool serialize_bin_object(
    yocto_shape* shape, const yocto_scene* scene, file_stream& fs, bool save) {
    if (!serialize_bin_value(shape->name, fs, save)) return false;
    if (!serialize_bin_value(shape->filename, fs, save)) return false;
    if (!serialize_bin_value(shape->points, fs, save)) return false;
    if (!serialize_bin_value(shape->lines, fs, save)) return false;
    if (!serialize_bin_value(shape->triangles, fs, save)) return false;
    if (!serialize_bin_value(shape->positions, fs, save)) return false;
    if (!serialize_bin_value(shape->normals, fs, save)) return false;
    if (!serialize_bin_value(shape->texturecoords, fs, save)) return false;
    if (!serialize_bin_value(shape->colors, fs, save)) return false;
    if (!serialize_bin_value(shape->radius, fs, save)) return false;
    if (!serialize_bin_value(shape->tangentspaces, fs, save)) return false;
    return true;
}

bool serialize_bin_object(yocto_surface* surface, file_stream& fs, bool save) {
    if (!serialize_bin_value(surface->name, fs, save)) return false;
    if (!serialize_bin_value(surface->filename, fs, save)) return false;
    if (!serialize_bin_value(surface->subdivision_level, fs, save))
        return false;
    if (!serialize_bin_value(surface->catmull_clark, fs, save)) return false;
    if (!serialize_bin_value(surface->compute_vertex_normals, fs, save))
        return false;
    if (!serialize_bin_value(surface->positions_quads, fs, save)) return false;
    if (!serialize_bin_value(surface->texturecoords_quads, fs, save))
        return false;
    if (!serialize_bin_value(surface->colors_quads, fs, save)) return false;
    if (!serialize_bin_value(surface->positions_creases, fs, save))
        return false;
    if (!serialize_bin_value(surface->texturecoords_quads, fs, save))
        return false;
    if (!serialize_bin_value(surface->positions, fs, save)) return false;
    if (!serialize_bin_value(surface->texturecoords, fs, save)) return false;
    if (!serialize_bin_value(surface->colors, fs, save)) return false;
    return true;
}

bool serialize_bin_object(yocto_texture* tex, file_stream& fs, bool save) {
    if (!serialize_bin_value(tex->name, fs, save)) return false;
    if (!serialize_bin_value(tex->filename, fs, save)) return false;
    if (!serialize_bin_value(tex->hdr_image, fs, save)) return false;
    if (!serialize_bin_value(tex->ldr_image, fs, save)) return false;
    if (!serialize_bin_value(tex->clamp_to_edge, fs, save)) return false;
    if (!serialize_bin_value(tex->height_scale, fs, save)) return false;
    if (!serialize_bin_value(tex->no_interpolation, fs, save)) return false;
    if (!serialize_bin_value(tex->ldr_as_linear, fs, save)) return false;
    if (!serialize_bin_value(tex->has_opacity, fs, save)) return false;
    return true;
}

bool serialize_bin_object(yocto_voltexture* tex, file_stream& fs, bool save) {
    if (!serialize_bin_value(tex->name, fs, save)) return false;
    if (!serialize_bin_value(tex->filename, fs, save)) return false;
    if (!serialize_bin_value(tex->volume_data, fs, save)) return false;
    if (!serialize_bin_value(tex->clamp_to_edge, fs, save)) return false;
    return true;
}

bool serialize_bin_object(yocto_environment* environment,
    const yocto_scene* scene, file_stream& fs, bool save) {
    if (!serialize_bin_value(environment->name, fs, save)) return false;
    if (!serialize_bin_value(environment->frame, fs, save)) return false;
    if (!serialize_bin_value(environment->emission, fs, save)) return false;
    if (!serialize_bin_handle(
            environment->emission_texture, scene->textures, fs, save))
        return false;
    return true;
}

bool serialize_bin_object(
    yocto_material* mat, const yocto_scene* scene, file_stream& fs, bool save) {
    if (!serialize_bin_value(mat->name, fs, save)) return false;
    if (!serialize_bin_value(mat->base_metallic, fs, save)) return false;
    if (!serialize_bin_value(mat->gltf_textures, fs, save)) return false;
    if (!serialize_bin_value(mat->double_sided, fs, save)) return false;
    if (!serialize_bin_value(mat->emission, fs, save)) return false;
    if (!serialize_bin_value(mat->diffuse, fs, save)) return false;
    if (!serialize_bin_value(mat->specular, fs, save)) return false;
    if (!serialize_bin_value(mat->transmission, fs, save)) return false;
    if (!serialize_bin_value(mat->roughness, fs, save)) return false;
    if (!serialize_bin_value(mat->opacity, fs, save)) return false;
    if (!serialize_bin_value(mat->fresnel, fs, save)) return false;
    if (!serialize_bin_value(mat->refract, fs, save)) return false;
    if (!serialize_bin_handle(mat->emission_texture, scene->textures, fs, save))
        return false;
    if (!serialize_bin_handle(mat->diffuse_texture, scene->textures, fs, save))
        return false;
    if (!serialize_bin_handle(mat->specular_texture, scene->textures, fs, save))
        return false;
    if (!serialize_bin_handle(
            mat->transmission_texture, scene->textures, fs, save))
        return false;
    if (!serialize_bin_handle(mat->roughness_texture, scene->textures, fs, save))
        return false;
    if (!serialize_bin_handle(mat->opacity_texture, scene->textures, fs, save))
        return false;
    if (!serialize_bin_handle(mat->occlusion_texture, scene->textures, fs, save))
        return false;
    if (!serialize_bin_handle(mat->bump_texture, scene->textures, fs, save))
        return false;
    if (!serialize_bin_handle(
            mat->displacement_texture, scene->textures, fs, save))
        return false;
    if (!serialize_bin_handle(mat->normal_texture, scene->textures, fs, save))
        return false;
    if (!serialize_bin_value(mat->volume_emission, fs, save)) return false;
    if (!serialize_bin_value(mat->volume_albedo, fs, save)) return false;
    if (!serialize_bin_value(mat->volume_density, fs, save)) return false;
    if (!serialize_bin_value(mat->volume_phaseg, fs, save)) return false;
    if (!serialize_bin_handle(
            mat->volume_density_texture, scene->voltextures, fs, save))
        return false;
    return true;
};

bool serialize_bin_object(yocto_instance* instance, const yocto_scene* scene,
    file_stream& fs, bool save) {
    if (!serialize_bin_value(instance->name, fs, save)) return false;
    if (!serialize_bin_value(instance->frame, fs, save)) return false;
    if (!serialize_bin_handle(instance->shape, scene->shapes, fs, save))
        return false;
    if (!serialize_bin_handle(instance->material, scene->materials, fs, save))
        return false;
    if (!serialize_bin_handle(instance->surface, scene->surfaces, fs, save))
        return false;
    return true;
};

bool serialize_scene(yocto_scene* scene, file_stream& fs, bool save) {
    if (!serialize_bin_value(scene->name, fs, save)) return false;
    if (!serialize_bin_object(scene->cameras, fs, save)) return false;
    if (!serialize_bin_object(scene->shapes, scene, fs, save)) return false;
    if (!serialize_bin_object(scene->surfaces, fs, save)) return false;
    if (!serialize_bin_object(scene->textures, fs, save)) return false;
    if (!serialize_bin_object(scene->voltextures, fs, save)) return false;
    if (!serialize_bin_object(scene->materials, scene, fs, save)) return false;
    if (!serialize_bin_object(scene->instances, scene, fs, save)) return false;
    if (!serialize_bin_object(scene->environments, scene, fs, save))
        return false;
    return true;
}

// Load/save a binary dump useful for very fast scene IO.
yocto_scene* load_ybin_scene(
    const string& filename, bool load_textures, bool skip_missing) {
    auto scope = log_trace_scoped("loading scene {}", filename);
    auto fs    = open(filename, "rb");
    if (!fs) return nullptr;
    auto scene = make_unique<yocto_scene>();
    if (!serialize_scene(scene.get(), fs, false)) return nullptr;
    return scene.release();
}

// Load/save a binary dump useful for very fast scene IO.
bool save_ybin_scene(const string& filename, const yocto_scene* scene,
    bool save_textures, bool skip_missing) {
    auto fs = open(filename, "wb");
    if (!fs) return false;
    if (!serialize_scene((yocto_scene*)scene, fs, true)) return false;
    return true;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE IO
// -----------------------------------------------------------------------------
namespace ygl {

// Reset mesh data
void reset_mesh_data(vector<int>& points, vector<vec2i>& lines,
    vector<vec3i>& triangles, vector<vec4i>& quads, vector<vec3f>& pos,
    vector<vec3f>& norm, vector<vec2f>& texcoord, vector<vec4f>& color,
    vector<float>& radius) {
    points    = {};
    lines     = {};
    triangles = {};
    quads     = {};
    pos       = {};
    norm      = {};
    texcoord  = {};
    color     = {};
    radius    = {};
}

// merge quads and triangles
void merge_triangles_and_quads(
    vector<vec3i>& triangles, vector<vec4i>& quads, bool force_triangles) {
    if (quads.empty()) return;
    if (force_triangles) {
        auto qtriangles = convert_quads_to_triangles(quads);
        triangles.insert(triangles.end(), qtriangles.begin(), qtriangles.end());
        quads = {};
    } else {
        auto tquads = convert_triangles_to_quads(triangles);
        quads.insert(quads.end(), tquads.begin(), tquads.end());
        triangles = {};
    }
}

// Load ply mesh
bool load_mesh(const string& filename, vector<int>& points,
    vector<vec2i>& lines, vector<vec3i>& triangles, vector<vec4i>& quads,
    vector<vec3f>& pos, vector<vec3f>& norm, vector<vec2f>& texcoord,
    vector<vec4f>& color, vector<float>& radius, bool force_triangles) {
    auto ext = get_extension(filename);
    if (ext == "ply" || ext == "PLY") {
        return load_ply_mesh(filename, points, lines, triangles, quads, pos,
            norm, texcoord, color, radius, force_triangles);
    } else if (ext == "obj" || ext == "OBJ") {
        return load_obj_mesh(filename, points, lines, triangles, quads, pos,
            norm, texcoord, force_triangles);
    } else {
        reset_mesh_data(points, lines, triangles, quads, pos, norm, texcoord,
            color, radius);
        return false;
    }
}

// Save ply mesh
bool save_mesh(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    const vector<vec4f>& colors, const vector<float>& radius, bool ascii) {
    auto ext = get_extension(filename);
    if (ext == "ply" || ext == "PLY") {
        return save_ply_mesh(filename, points, lines, triangles, quads,
            positions, normals, texturecoords, colors, radius, ascii);
    } else if (ext == "obj" || ext == "OBJ") {
        return save_obj_mesh(filename, points, lines, triangles, quads,
            positions, normals, texturecoords);
    } else {
        return false;
    }
}

// prepare obj line (remove comments and normalize whitespace)
void normalize_ply_line(char* s) {
    while (*s) {
        if (*s == '\t' || *s == '\r' || *s == '\n') {
            *s++ = ' ';
        } else {
            s++;
        }
    }
}

// Load ply mesh
ply_data* load_ply(const string& filename) {
    // open file
    auto fs = open(filename, "rb");
    if (!fs) return {};

    // parse header
    auto ply   = make_unique<ply_data>();
    auto ascii = false;
    char buf[4096];
    auto line = ""s;
    while (read_line(fs, line)) {
        memcpy(buf, line.c_str(), line.size() + 1);
        normalize_ply_line(buf);
        auto ss  = buf;
        auto cmd = parse_string(ss);
        if (cmd == "") continue;
        if (cmd == "ply") {
        } else if (cmd == "comment") {
        } else if (cmd == "format") {
            auto fmt = parse_string(ss);
            if (fmt != "ascii" && fmt != "binary_little_endian") return {};
            ascii = fmt == "ascii";
        } else if (cmd == "element") {
            auto elem  = ply_element();
            elem.name  = parse_string(ss);
            elem.count = parse_int(ss);
            ply->elements.push_back(elem);
        } else if (cmd == "property") {
            auto prop = ply_property();
            auto type = parse_string(ss);
            if (type == "list") {
                auto count_type = parse_string(ss);
                auto elem_type  = parse_string(ss);
                if (count_type != "uchar" && count_type != "uint8")
                    log_error("unsupported ply list type");
                if (elem_type != "int") log_error("unsupported ply list type");
                prop.type = ply_type::ply_int_list;
            } else if (type == "float") {
                prop.type = ply_type::ply_float;
            } else if (type == "uchar" || type == "uint8") {
                prop.type = ply_type::ply_uchar;
            } else if (type == "int") {
                prop.type = ply_type::ply_int;
            } else {
                return {};
            }
            prop.name = parse_string(ss);
            prop.scalars.resize(ply->elements.back().count);
            if (prop.type == ply_type::ply_int_list)
                prop.lists.resize(ply->elements.back().count);
            ply->elements.back().properties.push_back(prop);
        } else if (cmd == "end_header") {
            break;
        } else {
            return {};
        }
    }

    // parse content
    for (auto& elem : ply->elements) {
        for (auto vid = 0; vid < elem.count; vid++) {
            auto ss = (char*)nullptr;
            if (ascii) {
                if (!read_line(fs, line)) return {};
                memcpy(buf, line.c_str(), line.size() + 1);
                ss = buf;
            }
            for (auto pid = 0; pid < elem.properties.size(); pid++) {
                auto& prop = elem.properties[pid];
                if (prop.type == ply_type::ply_float) {
                    auto v = 0.0f;
                    if (ascii) {
                        v = parse_float(ss);
                    } else {
                        if (!read_value(fs, v)) return {};
                    }
                    prop.scalars[vid] = v;
                } else if (prop.type == ply_type::ply_int) {
                    auto v = 0;
                    if (ascii) {
                        v = parse_int(ss);
                    } else {
                        if (!read_value(fs, v)) return {};
                    }
                    prop.scalars[vid] = v;
                } else if (prop.type == ply_type::ply_uchar) {
                    auto vc = (unsigned char)0;
                    if (ascii) {
                        auto v = parse_int(ss);
                        vc     = (unsigned char)v;
                    } else {
                        if (!read_value(fs, vc)) return {};
                    }
                    prop.scalars[vid] = vc / 255.0f;
                } else if (prop.type == ply_type::ply_int_list) {
                    auto vc = (unsigned char)0;
                    if (ascii) {
                        auto v = parse_int(ss);
                        vc     = (unsigned char)v;
                    } else {
                        if (!read_value(fs, vc)) return {};
                    }
                    prop.scalars[vid] = vc;
                    for (auto i = 0; i < (int)prop.scalars[vid]; i++)
                        if (ascii) {
                            prop.lists[vid][i] = parse_int(ss);
                        } else {
                            if (!read_value(fs, prop.lists[vid][i])) return {};
                        }
                } else {
                    return {};
                }
            }
        }
    }

    return ply.release();
}

#if YGL_HAPPLY

// Load ply mesh
bool load_ply_mesh(const string& filename, vector<int>& points,
    vector<vec2i>& lines, vector<vec3i>& triangles, vector<vec4i>& quads,
    vector<vec3f>& pos, vector<vec3f>& norm, vector<vec2f>& texcoord,
    vector<vec4f>& color, vector<float>& radius, bool force_triangles) {
    // reset data
    reset_mesh_data(
        points, lines, triangles, quads, pos, norm, texcoord, color, radius);

    // open ply
    auto ply = unique_ptr<happly::PLYData>();
    try {
        ply = make_unique<happly::PLYData>(filename);
    } catch (...) {
        log_io_error("could not open PLY {}", filename);
        return false;
    }

    // get vertex data
    if (ply->hasElement("vertex")) {
        auto& vertex_element = ply->getElement("vertex");
        if (vertex_element.hasProperty("x") &&
            vertex_element.hasProperty("y") && vertex_element.hasProperty("z")) {
            auto x = vertex_element.getProperty<float>("x");
            auto y = vertex_element.getProperty<float>("y");
            auto z = vertex_element.getProperty<float>("z");
            pos.resize(x.size());
            for (auto i = 0; i < x.size(); i++) pos[i] = {x[i], y[i], z[i]};
        }
        if (vertex_element.hasProperty("nx") && vertex_element.hasProperty("ny") &&
            vertex_element.hasProperty("nz")) {
            auto x = vertex_element.getProperty<float>("nx");
            auto y = vertex_element.getProperty<float>("ny");
            auto z = vertex_element.getProperty<float>("nz");
            norm.resize(x.size());
            for (auto i = 0; i < x.size(); i++) norm[i] = {x[i], y[i], z[i]};
        }
        if (vertex_element.hasProperty("u") && vertex_element.hasProperty("v")) {
            auto x = vertex_element.getProperty<float>("u");
            auto y = vertex_element.getProperty<float>("v");
            texcoord.resize(x.size());
            for (auto i = 0; i < x.size(); i++) texcoord[i] = {x[i], y[i]};
        }
        if (vertex_element.hasProperty("red") &&
            vertex_element.hasProperty("green") &&
            vertex_element.hasProperty("blue")) {
            auto x = vertex_element.getProperty<float>("red");
            auto y = vertex_element.getProperty<float>("green");
            auto z = vertex_element.getProperty<float>("blue");
            for (auto i = 0; i < x.size(); i++)
                color[i] = {x[i], y[i], z[i], 1};
            if (vertex_element.hasProperty("alpha")) {
                auto w = vertex_element.getProperty<float>("alpha");
                for (auto i = 0; i < x.size(); i++) color[i].w = w[i];
            }
        }
        if (vertex_element.hasProperty("radius")) {
            auto x = vertex_element.getProperty<float>("radius");
            radius.resize(x.size());
            for (auto i = 0; i < x.size(); i++) radius[i] = x[i];
        }
    }

    // faces
    if (ply->hasElement("face")) {
        auto& face_element = ply->getElement("face");
        if (face_element.hasProperty("vertex_indices")) {
            auto indices = face_element.getListProperty<int>("vertex_indices");
            for (auto& f : indices) {
                if (f.size() == 4) {
                    quads.push_back({f[0], f[1], f[2], f[3]});
                } else {
                    for (auto i = 2; i < f.size(); i++)
                        triangles.push_back({f[0], f[i - 1], f[i]});
                }
            }
        }
    }

    // lines
    if (ply->hasElement("line")) {
        auto& face_element = ply->getElement("line");
        if (face_element.hasProperty("vertex_indices")) {
            auto indices = face_element.getListProperty<int>("vertex_indices");
            for (auto& l : indices) {
                for (auto i = 1; i < l.size(); i++)
                    lines.push_back({l[i - 1], l[i]});
            }
        }
    }

    // done
    return true;
}

// Save ply mesh
bool save_ply_mesh(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    const vector<vec4f>& colors, const vector<float>& radius, bool ascii) {
    using namespace happly;

    PLYData ply;

    auto get_channel = [](const auto& data, int channel_index) -> vector<float> {
        auto channel = vector<float>(data.size());
        for (auto i = 0; i < channel.size(); i++)
            channel[i] = (&data[i].x)[channel_index];
        return channel;
    };

    if (!positions.empty()) {
        if (!ply.hasElement("vertex"))
            ply.addElement("vertex", positions.size());
        ply.getElement("vertex").addProperty("x", get_channel(positions, 0));
        ply.getElement("vertex").addProperty("y", get_channel(positions, 1));
        ply.getElement("vertex").addProperty("z", get_channel(positions, 2));
    }

    if (!normals.empty()) {
        if (!ply.hasElement("vertex")) ply.addElement("vertex", normals.size());
        ply.getElement("vertex").addProperty("nx", get_channel(normals, 0));
        ply.getElement("vertex").addProperty("ny", get_channel(normals, 1));
        ply.getElement("vertex").addProperty("nz", get_channel(normals, 2));
    }

    if (!texturecoords.empty()) {
        if (!ply.hasElement("vertex"))
            ply.addElement("vertex", texturecoords.size());
        ply.getElement("vertex").addProperty("u", get_channel(texturecoords, 0));
        ply.getElement("vertex").addProperty("v", get_channel(texturecoords, 1));
    }

    if (!colors.empty()) {
        if (!ply.hasElement("vertex")) ply.addElement("vertex", colors.size());
        ply.getElement("vertex").addProperty("red", get_channel(colors, 0));
        ply.getElement("vertex").addProperty("green", get_channel(colors, 1));
        ply.getElement("vertex").addProperty("blue", get_channel(colors, 2));
        ply.getElement("vertex").addProperty("alpha", get_channel(colors, 3));
    }

    if (!radius.empty()) {
        if (!ply.hasElement("vertex")) ply.addElement("vertex", radius.size());
        ply.getElement("vertex").addProperty("radius", radius);
    }

    if (!triangles.empty() || !quads.empty()) {
        if (!ply.hasElement("face"))
            ply.addElement("face", triangles.size() + quads.size());
        auto face_property = vector<vector<int>>();
        for (auto t : triangles) face_property.push_back({t.x, t.y, t.z});
        for (auto q : quads) {
            if (q.z == q.w) {
                face_property.push_back({q.x, q.y, q.z});
            } else {
                face_property.push_back({q.x, q.y, q.z, q.w});
            }
        }
        ply.getElement("face").addListProperty("vertex_indices", face_property);
    }

    if (!lines.empty()) {
        if (!ply.hasElement("line")) ply.addElement("line", lines.size());
        auto line_property = vector<vector<int>>();
        for (auto l : lines) line_property.push_back({l.x, l.y});
        ply.getElement("line").addListProperty("vertex_indices", line_property);
    }

    // save
    try {
        ply.write(filename, ascii ? DataFormat::ASCII : DataFormat::Binary);
    } catch (...) {
        return false;
    }

    return true;
}

#else

// Load ply mesh
bool load_ply_mesh(const string& filename, vector<int>& points,
    vector<vec2i>& lines, vector<vec3i>& triangles, vector<vec4i>& quads,
    vector<vec3f>& pos, vector<vec3f>& norm, vector<vec2f>& texcoord,
    vector<vec4f>& color, vector<float>& radius, bool force_triangles) {
    // clear
    reset_mesh_data(
        points, lines, triangles, quads, pos, norm, texcoord, color, radius);

    // load ply
    auto ply = unique_ptr<ply_data>{load_ply(filename)};
    if (ply->elements.empty()) {
        log_io_error("empty ply file {}", filename);
        return false;
    }

    // copy vertex data
    for (auto& elem : ply->elements) {
        if (elem.name != "vertex") continue;
        auto count = elem.count;
        for (auto& prop : elem.properties) {
            auto vals        = prop.scalars.data();
            auto copy_floats = [vals, count](auto& vert, const auto& def,
                                   int stride, int offset) {
                if (vert.size() != count) vert.resize(count, def);
                auto dst = (float*)vert.data();
                for (auto i = 0; i < count; i++)
                    dst[i * stride + offset] = vals[i];
            };
            if (prop.name == "x") copy_floats(pos, zero3f, 3, 0);
            if (prop.name == "y") copy_floats(pos, zero3f, 3, 1);
            if (prop.name == "z") copy_floats(pos, zero3f, 3, 2);
            if (prop.name == "nx") copy_floats(norm, zero3f, 3, 0);
            if (prop.name == "ny") copy_floats(norm, zero3f, 3, 1);
            if (prop.name == "nz") copy_floats(norm, zero3f, 3, 2);
            if (prop.name == "u") copy_floats(texcoord, zero2f, 2, 0);
            if (prop.name == "v") copy_floats(texcoord, zero2f, 2, 1);
            if (prop.name == "red") copy_floats(color, vec4f{0, 0, 0, 1}, 4, 0);
            if (prop.name == "green")
                copy_floats(color, vec4f{0, 0, 0, 1}, 4, 1);
            if (prop.name == "blue")
                copy_floats(color, vec4f{0, 0, 0, 1}, 4, 2);
            if (prop.name == "alpha")
                copy_floats(color, vec4f{0, 0, 0, 1}, 4, 3);
            if (prop.name == "radius") copy_floats(radius, 0.0f, 1, 0);
        }
    }

    // copy triangle data
    for (auto& elem : ply->elements) {
        if (elem.name != "face") continue;
        auto count = elem.count;
        for (auto& prop : elem.properties) {
            if (prop.name == "vertex_indices") {
                for (auto fid = 0; fid < count; fid++) {
                    auto& list = prop.lists[fid];
                    auto  num  = (int)prop.scalars[fid];
                    if (num == 4) {
                        quads.push_back({list[0], list[1], list[2], list[3]});
                    } else {
                        for (auto i = 2; i < num; i++)
                            triangles.push_back({list[0], list[i - 1], list[i]});
                    }
                }
            }
        }
    }

    merge_triangles_and_quads(triangles, quads, force_triangles);

    // done
    return true;
}

// Save ply mesh
bool save_ply_mesh(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    const vector<vec4f>& colors, const vector<float>& radius, bool ascii) {
    auto fs = open(filename, "wb");
    if (!fs) return false;

    // header
    print(fs, "ply\n");
    if (ascii)
        print(fs, "format ascii 1.0\n");
    else
        print(fs, "format binary_little_endian 1.0\n");
    print(fs, "element vertex {}\n", (int)pos.size());
    if (!pos.empty())
        print(fs, "property float x\nproperty float y\nproperty float z\n");
    if (!norm.empty())
        print(fs,
            "property float nx\nproperty float ny\nproperty float "
            "nz\n");
    if (!texcoord.empty()) print(fs, "property float u\nproperty float v\n");
    if (!color.empty())
        print(fs,
            "property float red\nproperty float green\nproperty float "
            "blue\nproperty float alpha\n");
    if (!radius.empty()) print(fs, "property float radius\n");
    if (!lines.empty()) {
        print(fs, "element line {}\n", (int)lines.size());
        print(fs, "property list uchar int vertex_indices\n");
    }
    if (!triangles.empty() || !quads.empty()) {
        print(fs, "element face {}\n", (int)triangles.size() + (int)quads.size());
        print(fs, "property list uchar int vertex_indices\n");
    }
    print(fs, "end_header\n");

    // body
    if (ascii) {
        // write vertex data
        for (auto i = 0; i < pos.size(); i++) {
            if (!pos.empty()) print(fs, "{} ", pos[i]);
            if (!norm.empty()) print(fs, "{} ", norm[i]);
            if (!texcoord.empty()) print(fs, "{} ", texcoord[i]);
            if (!color.empty()) print(fs, "{} ", color[i]);
            if (!radius.empty()) print(fs, "{} ", radius[i]);
            print(fs, "\n");
        }

        // write face data
        for (auto& l : lines) print(fs, "2 {}\n", l);
        for (auto& t : triangles) print(fs, "3 {}\n", t);
        for (auto& q : quads) {
            if (q.z == q.w)
                print(fs, "3 {}\n", xyz(q));
            else
                print(fs, "4 {}\n", q);
        }
    } else {
        // write vertex data
        for (auto i = 0; i < pos.size(); i++) {
            if (!pos.empty()) write_value(fs, pos[i]);
            if (!norm.empty()) write_value(fs, norm[i]);
            if (!texcoord.empty()) write_value(fs, texcoord[i]);
            if (!color.empty()) write_value(fs, color[i]);
            if (!radius.empty()) write_value(fs, radius[i]);
        }

        // write face data
        for (auto& l : lines) {
            auto n = (byte)2;
            write_value(fs, n);
            write_value(fs, l);
        }
        for (auto& t : triangles) {
            auto n = (byte)3;
            write_value(fs, n);
            write_value(fs, t);
        }
        for (auto& q : quads) {
            if (q.z == q.w) {
                auto n = (byte)3;
                write_value(fs, n);
                write_value(fs, xyz(q));
            } else {
                auto n = (byte)4;
                write_value(fs, n);
                write_value(fs, q);
            }
        }
    }

    // done
    return true;
}

#endif

// Load ply mesh
bool load_obj_mesh(const string& filename, vector<int>& points,
    vector<vec2i>& lines, vector<vec3i>& triangles, vector<vec4i>& quads,
    vector<vec3f>& pos, vector<vec3f>& norm, vector<vec2f>& texcoord,
    bool flip_texcoord, bool force_triangles) {
    // clear
    auto color  = vector<vec4f>{};
    auto radius = vector<float>{};
    reset_mesh_data(
        points, lines, triangles, quads, pos, norm, texcoord, color, radius);

    // obj vertices
    auto opos      = std::deque<vec3f>();
    auto onorm     = std::deque<vec3f>();
    auto otexcoord = std::deque<vec2f>();

    // vertex maps
    auto vert_map = unordered_map<obj_vertex, int, obj_vertex_hash>();

    // Add  vertices to the current shape
    auto add_verts = [&](const vector<obj_vertex>& verts) {
        for (auto& vert : verts) {
            auto it = vert_map.find(vert);
            if (it != vert_map.end()) continue;
            auto nverts = (int)pos.size();
            vert_map.insert(it, {vert, nverts});
            if (vert.position) pos.push_back(opos.at(vert.position - 1));
            if (vert.texturecoord)
                texcoord.push_back(otexcoord.at(vert.texturecoord - 1));
            if (vert.normal) norm.push_back(onorm.at(vert.normal - 1));
        }
    };

    auto cb     = obj_callbacks();
    cb.vert     = [&](vec3f v) { opos.push_back(v); };
    cb.norm     = [&](vec3f v) { onorm.push_back(v); };
    cb.texcoord = [&](vec2f v) { otexcoord.push_back(v); };
    cb.face     = [&](const vector<obj_vertex>& verts) {
        add_verts(verts);
        if (verts.size() == 4) {
            quads.push_back({vert_map.at(verts[0]), vert_map.at(verts[1]),
                vert_map.at(verts[2]), vert_map.at(verts[3])});
        } else {
            for (auto i = 2; i < verts.size(); i++)
                triangles.push_back({vert_map.at(verts[0]),
                    vert_map.at(verts[i - 1]), vert_map.at(verts[i])});
        }
    };
    cb.line = [&](const vector<obj_vertex>& verts) {
        add_verts(verts);
        for (auto i = 1; i < verts.size(); i++)
            lines.push_back({vert_map.at(verts[i - 1]), vert_map.at(verts[i])});
    };
    cb.point = [&](const vector<obj_vertex>& verts) {
        add_verts(verts);
        for (auto i = 0; i < verts.size(); i++)
            points.push_back(vert_map.at(verts[i]));
    };

    // merging quads and triangles
    merge_triangles_and_quads(triangles, quads, force_triangles);

    // load obj
    return load_obj(filename, cb, true, true, flip_texcoord);
}

// Load ply mesh
bool save_obj_mesh(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    bool flip_texcoord) {
    auto fs = open(filename, "wt");
    if (!fs) return false;

    for (auto& p : positions) print(fs, "v {}\n", p);
    for (auto& n : normals) print(fs, "vn {}\n", n);
    for (auto& t : texturecoords)
        print(fs, "vt {}\n", vec2f{t.x, (flip_texcoord) ? 1 - t.y : t.y});
    auto mask = obj_vertex{
        1, texturecoords.empty() ? 0 : 1, normals.empty() ? 0 : 1};
    auto vert = [mask](int i) {
        return obj_vertex{(i + 1) * mask.position, (i + 1) * mask.texturecoord,
            (i + 1) * mask.normal};
    };
    for (auto& p : points) {
        print(fs, "p {}\n", to_string(vert(p)).c_str());
    }
    for (auto& l : lines) {
        print(fs, "l {} {}\n", to_string(vert(l.x)).c_str(),
            to_string(vert(l.y)).c_str());
    }
    for (auto& t : triangles) {
        print(fs, "f {} {} {}\n", to_string(vert(t.x)).c_str(),
            to_string(vert(t.y)).c_str(), to_string(vert(t.z)).c_str());
    }
    for (auto& q : quads) {
        if (q.z == q.w) {
            print(fs, "f {} {} {}\n", to_string(vert(q.x)).c_str(),
                to_string(vert(q.y)).c_str(), to_string(vert(q.z)).c_str());
        } else {
            print(fs, "f {} {} {} {}\n", to_string(vert(q.x)).c_str(),
                to_string(vert(q.y)).c_str(), to_string(vert(q.z)).c_str(),
                to_string(vert(q.w)).c_str());
        }
    }

    return true;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE IO
// -----------------------------------------------------------------------------
namespace ygl {

// Reset mesh data
void reset_fvmesh_data(vector<vec4i>& quads_positions, vector<vec3f>& pos,
    vector<vec4i>& quads_normals, vector<vec3f>& norm,
    vector<vec4i>& quads_texturecoords, vector<vec2f>& texcoord,
    vector<vec4i>& quads_colors, vector<vec4f>& color) {
    quads_positions     = {};
    pos                 = {};
    quads_normals       = {};
    norm                = {};
    quads_texturecoords = {};
    texcoord            = {};
    quads_colors        = {};
    color               = {};
}

// Load mesh
bool load_fvmesh(const string& filename, vector<vec4i>& quads_positions,
    vector<vec3f>& pos, vector<vec4i>& quads_normals, vector<vec3f>& norm,
    vector<vec4i>& quads_texturecoords, vector<vec2f>& texcoord,
    vector<vec4i>& quads_colors, vector<vec4f>& color) {
    auto ext = get_extension(filename);
    if (ext == "obj" || ext == "OBJ") {
        return load_obj_fvmesh(filename, quads_positions, pos, quads_normals,
            norm, quads_texturecoords, texcoord);
    } else {
        reset_fvmesh_data(quads_positions, pos, quads_normals, norm,
            quads_texturecoords, texcoord, quads_colors, color);
        log_io_error("unsupported mesh format {}", ext);
        return false;
    }
}

// Save mesh
bool save_fvmesh(const string& filename, const vector<vec4i>& quads_positions,
    const vector<vec3f>& positions, const vector<vec4i>& quads_normals,
    const vector<vec3f>& normals, const vector<vec4i>& quads_texturecoords,
    const vector<vec2f>& texturecoords, const vector<vec4i>& quads_colors,
    const vector<vec4f>& colors, bool ascii) {
    auto ext = get_extension(filename);
    if (ext == "obj" || ext == "OBJ") {
        return save_obj_fvmesh(filename, quads_positions, positions,
            quads_normals, normals, quads_texturecoords, texturecoords);
    } else {
        log_io_error("unsupported mesh format {}", ext);
        return false;
    }
}

// Load obj mesh
bool load_obj_fvmesh(const string& filename, vector<vec4i>& quads_positions,
    vector<vec3f>& pos, vector<vec4i>& quads_normals, vector<vec3f>& norm,
    vector<vec4i>& quads_texturecoords, vector<vec2f>& texcoord,
    bool flip_texcoord) {
    // clear
    vector<vec4i> quads_colors;
    vector<vec4f> color;
    reset_fvmesh_data(quads_positions, pos, quads_normals, norm,
        quads_texturecoords, texcoord, quads_colors, color);

    // obj vertex
    auto opos      = std::deque<vec3f>();
    auto onorm     = std::deque<vec3f>();
    auto otexcoord = std::deque<vec2f>();

    // vertex maps
    auto pos_map      = unordered_map<int, int>();
    auto texcoord_map = unordered_map<int, int>();
    auto norm_map     = unordered_map<int, int>();

    // add vertex
    auto add_verts = [&](const vector<obj_vertex>& verts) {
        for (auto& vert : verts) {
            if (!vert.position) continue;
            auto pos_it = pos_map.find(vert.position);
            if (pos_it != pos_map.end()) continue;
            auto nverts = (int)pos.size();
            pos_map.insert(pos_it, {vert.position, nverts});
            pos.push_back(opos.at(vert.position - 1));
        }
        for (auto& vert : verts) {
            if (!vert.texturecoord) continue;
            auto texcoord_it = texcoord_map.find(vert.texturecoord);
            if (texcoord_it != texcoord_map.end()) continue;
            auto nverts = (int)texcoord.size();
            texcoord_map.insert(texcoord_it, {vert.texturecoord, nverts});
            texcoord.push_back(otexcoord.at(vert.texturecoord - 1));
        }
        for (auto& vert : verts) {
            if (!vert.normal) continue;
            auto norm_it = norm_map.find(vert.normal);
            if (norm_it != norm_map.end()) continue;
            auto nverts = (int)norm.size();
            norm_map.insert(norm_it, {vert.normal, nverts});
            norm.push_back(onorm.at(vert.normal - 1));
        }
    };

    auto cb     = obj_callbacks();
    cb.vert     = [&](vec3f v) { opos.push_back(v); };
    cb.norm     = [&](vec3f v) { onorm.push_back(v); };
    cb.texcoord = [&](vec2f v) { otexcoord.push_back(v); };
    cb.face     = [&](const vector<obj_vertex>& verts) {
        add_verts(verts);
        if (verts.size() == 4) {
            if (verts[0].position) {
                quads_positions.push_back({pos_map.at(verts[0].position),
                    pos_map.at(verts[1].position), pos_map.at(verts[2].position),
                    pos_map.at(verts[3].position)});
            }
            if (verts[0].texturecoord) {
                quads_texturecoords.push_back(
                    {texcoord_map.at(verts[0].texturecoord),
                        texcoord_map.at(verts[1].texturecoord),
                        texcoord_map.at(verts[2].texturecoord),
                        texcoord_map.at(verts[3].texturecoord)});
            }
            if (verts[0].normal) {
                quads_normals.push_back({norm_map.at(verts[0].normal),
                    norm_map.at(verts[1].normal), norm_map.at(verts[2].normal),
                    norm_map.at(verts[3].normal)});
            }
        } else {
            if (verts[0].position) {
                for (auto i = 2; i < verts.size(); i++)
                    quads_positions.push_back({pos_map.at(verts[0].position),
                        pos_map.at(verts[1].position),
                        pos_map.at(verts[i].position),
                        pos_map.at(verts[i].position)});
            }
            if (verts[0].texturecoord) {
                for (auto i = 2; i < verts.size(); i++)
                    quads_texturecoords.push_back(
                        {texcoord_map.at(verts[0].texturecoord),
                            texcoord_map.at(verts[1].texturecoord),
                            texcoord_map.at(verts[i].texturecoord),
                            texcoord_map.at(verts[i].texturecoord)});
            }
            if (verts[0].normal) {
                for (auto i = 2; i < verts.size(); i++)
                    quads_normals.push_back({norm_map.at(verts[0].normal),
                        norm_map.at(verts[1].normal), norm_map.at(verts[i].normal),
                        norm_map.at(verts[i].normal)});
            }
        }
    };

    // load obj
    return load_obj(filename, cb, true, true, flip_texcoord);
}

// Load ply mesh
bool save_obj_fvmesh(const string& filename, const vector<vec4i>& quads_positions,
    const vector<vec3f>& positions, const vector<vec4i>& quads_normals,
    const vector<vec3f>& normals, const vector<vec4i>& quads_texturecoords,
    const vector<vec2f>& texturecoords, bool flip_texcoord) {
    auto fs = open(filename, "wt");
    if (!fs) return false;

    for (auto& p : positions) print(fs, "v {}\n", p);
    for (auto& n : normals) print(fs, "vn {}\n", n);
    for (auto& t : texturecoords)
        print(fs, "vt {}\n", vec2f{t.x, (flip_texcoord) ? 1 - t.y : t.y});
    auto mask = obj_vertex{
        1, texturecoords.empty() ? 0 : 1, normals.empty() ? 0 : 1};
    auto vert = [mask](int pif, int ti, int ni) {
        return obj_vertex{(pif + 1) * mask.position,
            (ti + 1) * mask.texturecoord, (ni + 1) * mask.normal};
    };
    for (auto i = 0; i < quads_positions.size(); i++) {
        auto qp = quads_positions.at(i);
        auto qt = !quads_texturecoords.empty() ? quads_texturecoords.at(i) :
                                                 vec4i{-1, -1, -1, -1};
        auto qn = !quads_normals.empty() ? quads_normals.at(i) :
                                           vec4i{-1, -1, -1, -1};
        if (qp.z != qp.w)
            print(fs, "f {} {} {} {}\n",
                to_string(vert(qp.x, qt.x, qn.x)).c_str(),
                to_string(vert(qp.y, qt.y, qn.y)).c_str(),
                to_string(vert(qp.z, qt.z, qn.z)).c_str(),
                to_string(vert(qp.w, qt.w, qn.w)).c_str());
        else
            print(fs, "f {} {} {}\n", to_string(vert(qp.x, qt.x, qn.x)).c_str(),
                to_string(vert(qp.y, qt.y, qn.y)).c_str(),
                to_string(vert(qp.z, qt.z, qn.z)).c_str());
    }

    return true;
}

}  // namespace ygl
