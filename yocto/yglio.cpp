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
    if (filename.size() > 1 && filename[0] == '/' && filename[1] == '/')
        throw runtime_error("no absolute paths");
    if (filename.size() > 3 && filename[1] == ':' && filename[2] == '/' &&
        filename[3] == '/')
        throw runtime_error("no absolute paths");
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

// Write to file
bool write_values(file_stream& fs, const string& vals) {
    if (!fs) return false;
    if (fwrite(vals.data(), 1, vals.size(), fs.fs) != vals.size()) {
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

// Read binary data to fill the whole buffer
bool read_values(file_stream& fs, string& vals) {
    if (!fs) return false;
    if (fread(vals.data(), 1, vals.size(), fs.fs) != vals.size()) {
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
    auto txt = load_text(filename);
    if (txt.empty()) return {};
    try {
        return json::parse(txt.begin(), txt.end());
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
image<vec4f> load_image4f(const string& filename) {
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

// Saves an hdr image.
bool save_image4f(const string& filename, const image<vec4f>& img) {
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

// Loads an hdr image.
image<vec4f> load_image4f_from_memory(const byte* data, int data_size) {
    return load_stbi_image4f_from_memory(data, data_size);
}

// Loads an hdr image.
image<vec4b> load_image4b(const string& filename) {
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

// Saves an ldr image.
bool save_image4b(const string& filename, const image<vec4b>& img) {
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

// Loads an ldr image.
image<vec4b> load_image4b_from_memory(const byte* data, int data_size) {
    return load_stb_image4b_from_memory(data, data_size);
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
    if (!width && !height) throw runtime_error("bad image size");
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
volume<float> load_volume1f(const string& filename) {
    auto fs = open(filename, "r");
    if (!fs) return {};
    auto size = zero3i;
    if (!read_value(fs, size)) return {};
    auto vol = volume<float>{size.x, size.y, size.z};
    if (!read_values(fs, vol.voxels)) return {};
    return vol;
}

// Saves volume data in binary format.
bool save_volume1f(const string& filename, const volume<float>& vol) {
    auto fs = open(filename, "w");
    if (!fs) return false;
    auto size = vec3i{vol.width, vol.height, vol.depth};
    if (!write_value(fs, size)) return false;
    if (!write_values(fs, vol.voxels)) return false;
    return true;
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
bool save_scene(const string& filename, const yocto_scene* scn,
    bool save_textures, bool skip_missing) {
    auto ext = get_extension(filename);
    if (ext == "json" || ext == "JSON") {
        return save_json_scene(filename, scn, save_textures, skip_missing);
    } else if (ext == "obj" || ext == "OBJ") {
        return save_obj_scene(filename, scn, save_textures, skip_missing);
    } else if (ext == "gltf" || ext == "GLTF") {
        return save_gltf_scene(filename, scn, save_textures, skip_missing);
    } else if (ext == "pbrt" || ext == "PBRT") {
        return save_pbrt_scene(filename, scn, save_textures, skip_missing);
    } else if (ext == "ybin" || ext == "YBIN") {
        return save_ybin_scene(filename, scn, save_textures, skip_missing);
    } else {
        log_io_error("unsupported scene format {}", ext);
        return false;
    }
}

bool load_scene_textures(yocto_scene* scn, const string& dirname,
    bool skip_missing, bool assign_opacity) {
    // load images
    for (auto txt : scn->textures) {
        if (txt->filename == "" || !txt->hdr_image.pixels.empty() ||
            !txt->ldr_image.pixels.empty())
            continue;
        auto filename = normalize_path(dirname + "/" + txt->filename);
        if (is_hdr_filename(filename)) {
            txt->hdr_image = load_image4f(filename);
        } else {
            txt->ldr_image = load_image4b(filename);
        }
        if (txt->hdr_image.pixels.empty() && txt->ldr_image.pixels.empty()) {
            if (!skip_missing) return false;
        }
    }

    // load volumes
    for (auto txt : scn->voltextures) {
        if (txt->filename == "" || !txt->volume_data.voxels.empty()) continue;
        auto filename    = normalize_path(dirname + "/" + txt->filename);
        txt->volume_data = load_volume1f(filename);
        if (txt->volume_data.voxels.empty()) {
            if (!skip_missing) return false;
        }
    }

    // assign opacity texture if needed
    if (assign_opacity) {
        auto has_opacity = unordered_map<yocto_texture*, bool>();
        for (auto& txt : scn->textures) {
            has_opacity[txt] = false;
            for (auto& p : txt->hdr_image.pixels) {
                if (p.w < 0.999f) {
                    has_opacity[txt] = true;
                    break;
                }
            }
            for (auto& p : txt->ldr_image.pixels) {
                if (p.w < 255) {
                    has_opacity[txt] = true;
                    break;
                }
            }
        }
        for (auto& mat : scn->materials) {
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
    const yocto_scene* scn, const string& dirname, bool skip_missing) {
    // save images
    for (auto txt : scn->textures) {
        if (txt->hdr_image.pixels.empty() && txt->ldr_image.pixels.empty())
            continue;
        auto filename = normalize_path(dirname + "/" + txt->filename);
        if (is_hdr_filename(filename)) {
            if (!save_image4f(filename, txt->hdr_image)) {
                if (!skip_missing) return false;
            }
        } else {
            if (!save_image4b(filename, txt->ldr_image)) {
                if (!skip_missing) return false;
            }
        }
    }

    // save volumes
    for (auto txt : scn->voltextures) {
        if (txt->volume_data.voxels.empty()) continue;
        auto filename = normalize_path(dirname + "/" + txt->filename);
        if (!save_volume1f(filename, txt->volume_data)) {
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
    } catch (...) { return false; }
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
    } catch (...) { return false; }
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
bool dump_json_objarray(json& js, const vector<T*>& val, const yocto_scene* scn) {
    js = json::array();
    for (auto& v : val) {
        js.push_back({});
        if (!dump_json_object(js.back(), v, scn)) return false;
    }
    return true;
}

// Dumps a json value
template <typename T>
bool dump_json_objarray(
    json& js, const vector<T*>& val, const char* name, const yocto_scene* scn) {
    if (val.empty()) return true;
    return dump_json_objarray(js[name], val, scn);
}

// Dumps a json value
template <typename T>
bool parse_json_objarray(const json& js, vector<T*>& val, const yocto_scene* scn) {
    if (!js.is_array()) return false;
    for (auto& j : js) {
        val.push_back(new T());
        if (!parse_json_object(j, val.back(), scn)) return false;
    }
    return true;
}

// Dumps a json value
template <typename T>
bool parse_json_objarray(
    const json& js, vector<T*>& val, const char* name, const yocto_scene* scn) {
    if (!js.count(name)) return true;
    val = {};
    return parse_json_objarray(js.at(name), val, scn);
}

// Parses and applied a JSON procedural
template <typename T>
bool parse_json_procedural(
    const json& js, T* val, const char* name, const yocto_scene* scn) {
    if (!js.count(name)) return true;
    return apply_json_procedural(js.at(name), val, scn);
}

// Procedural commands for cameras
bool apply_json_procedural(
    const json& js, yocto_camera* val, const yocto_scene* scn) {
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
bool dump_json_object(json& js, const yocto_camera* val, const yocto_scene* scn) {
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
bool parse_json_object(const json& js, yocto_camera* val, const yocto_scene* scn) {
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
    if (!parse_json_procedural(js, val, "!!proc", scn)) return false;
    return true;
}

// Serialize struct
bool dump_json_object(json& js, const yocto_texture* val, const yocto_scene* scn) {
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
    const json& js, yocto_texture* val, const yocto_scene* scn) {
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
        throw runtime_error("unknown texture type " + type);
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
    const json& js, yocto_texture* val, const yocto_scene* scn) {
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
    if (!parse_json_procedural(js, val, "!!proc", scn)) return false;
    return true;
}

// Serialize struct
bool dump_json_object(
    json& js, const yocto_voltexture* val, const yocto_scene* scn) {
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
    const json& js, yocto_voltexture* val, const yocto_scene* scn) {
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
        throw runtime_error("unknown texture type " + type);
    }
    if (val->filename == "") {
        auto ext      = string("vol");
        val->filename = "textures/" + val->name + "." + ext;
    }
    return true;
}

// Serialize struct
bool parse_json_object(
    const json& js, yocto_voltexture* val, const yocto_scene* scn) {
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
    if (!parse_json_procedural(js, val, "!!proc", scn)) return false;
    return true;
}

// Serialize struct
bool dump_json_object(
    json& js, const yocto_material* val, const yocto_scene* scn) {
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
            js, val->emission_texture, "emission_texture", scn->textures))
        return false;
    if (!dump_json_objref(
            js, val->diffuse_texture, "diffuse_texture", scn->textures))
        return false;
    if (!dump_json_objref(
            js, val->specular_texture, "specular_texture", scn->textures))
        return false;
    if (!dump_json_objref(js, val->transmission_texture, "transmission_texture",
            scn->textures))
        return false;
    if (!dump_json_objref(
            js, val->roughness_texture, "roughness_texture", scn->textures))
        return false;
    if (!dump_json_objref(
            js, val->opacity_texture, "opacity_texture", scn->textures))
        return false;
    if (!dump_json_objref(
            js, val->occlusion_texture, "occlusion_texture", scn->textures))
        return false;
    if (!dump_json_objref(js, val->bump_texture, "bump_texture", scn->textures))
        return false;
    if (!dump_json_objref(js, val->displacement_texture, "displacement_texture",
            scn->textures))
        return false;
    if (!dump_json_objref(
            js, val->normal_texture, "normal_texture", scn->textures))
        return false;
    if (!dump_json_objref(js, val->volume_density_texture,
            "volume_density_texture", scn->voltextures))
        return false;
    return true;
}

// Procedural commands for materials
bool apply_json_procedural(
    const json& js, yocto_material* val, const yocto_scene* scn) {
    if (!parse_json_objbegin(js)) return false;
    return true;
}

// Serialize struct
bool parse_json_object(
    const json& js, yocto_material* val, const yocto_scene* scn) {
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
            js, val->emission_texture, "emission_texture", scn->textures))
        return false;
    if (!parse_json_objref(
            js, val->diffuse_texture, "diffuse_texture", scn->textures))
        return false;
    if (!parse_json_objref(
            js, val->specular_texture, "specular_texture", scn->textures))
        return false;
    if (!parse_json_objref(js, val->transmission_texture,
            "transmission_texture", scn->textures))
        return false;
    if (!parse_json_objref(
            js, val->roughness_texture, "roughness_texture", scn->textures))
        return false;
    if (!parse_json_objref(
            js, val->opacity_texture, "opacity_texture", scn->textures))
        return false;
    if (!parse_json_objref(
            js, val->occlusion_texture, "occlusion_texture", scn->textures))
        return false;
    if (!parse_json_objref(js, val->bump_texture, "bump_texture", scn->textures))
        return false;
    if (!parse_json_objref(js, val->displacement_texture,
            "displacement_texture", scn->textures))
        return false;
    if (!parse_json_objref(
            js, val->normal_texture, "normal_texture", scn->textures))
        return false;
    if (!parse_json_objref(js, val->volume_density_texture,
            "volume_density_texture", scn->voltextures))
        return false;
    if (!parse_json_procedural(js, val, "!!proc", scn)) return false;
    return true;
}

// Serialize struct
bool dump_json_object(json& js, const yocto_shape* val, const yocto_scene* scn) {
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
                js, val->tangent_spaces, "tangent_spaces", def.tangent_spaces))
            return false;
    }
    return true;
}

// Procedural commands for materials
bool apply_json_procedural(
    const json& js, yocto_shape* val, const yocto_scene* scn) {
    if (!parse_json_objbegin(js)) return false;
    auto type = js.value("type", ""s);
    if (type == "") return true;
    auto shp = make_shape_data();
    if (type == "quad") {
        shp = make_quad(js.value("steps", vec2i{1, 1}),
            js.value("size", vec2f{2, 2}), js.value("uvsize", vec2f{1, 1}), true);
    } else if (type == "quady") {
        shp = make_quad(js.value("steps", vec2i{1, 1}),
            js.value("size", vec2f{2, 2}), js.value("uvsize", vec2f{1, 1}), true);
    } else if (type == "quad_stack") {
        shp = make_quad_stack(js.value("steps", vec3i{1, 1, 1}),
            js.value("size", vec3f{2, 2, 2}), js.value("uvsize", vec2f{1, 1}),
            true);
    } else if (type == "cube") {
        shp = make_cube(js.value("steps", vec3i{1, 1, 1}),
            js.value("size", vec3f{2, 2, 2}),
            js.value("uvsize", vec3f{1, 1, 1}), true);
    } else if (type == "cube_rounded") {
        shp = make_cube_rounded(js.value("steps", vec3i{32, 32, 32}),
            js.value("size", vec3f{2, 2, 2}),
            js.value("uvsize", vec3f{1, 1, 1}), js.value("radius", 0.3f), true);
    } else if (type == "sphere") {
        shp = make_sphere(js.value("steps", vec2i{64, 32}),
            js.value("size", 2.0f), js.value("uvsize", vec2f{1, 1}), true);
    } else if (type == "sphere_cube") {
        shp = make_sphere_cube(js.value("steps", 32), js.value("size", 2.0f),
            js.value("uvsize", 1.0f), true);
    } else if (type == "sphere_flipcap") {
        shp = make_sphere_flipcap(js.value("steps", vec2i{64, 32}),
            js.value("size", 2.0f), js.value("uvsize", vec2f{1, 1}),
            js.value("zflip", vec2f{-0.75f, +0.75f}), true);
    } else if (type == "disk") {
        shp = make_disk(js.value("steps", vec2i{32, 16}),
            js.value("size", 2.0f), js.value("uvsize", vec2f{1, 1}), true);
    } else if (type == "disk_quad") {
        shp = make_disk_quad(js.value("steps", 32), js.value("size", 2.0f),
            js.value("uvsize", 1.0f), true);
    } else if (type == "disk_bulged") {
        shp = make_disk_bulged(js.value("steps", 32), js.value("size", 2.0f),
            js.value("uvsize", 1.0f), js.value("height", 0.25f), true);
    } else if (type == "cylinder_side") {
        shp = make_cylinder_side(js.value("steps", vec2i{64, 32}),
            js.value("size", vec2f{2.0f, 2.0f}),
            js.value("uvsize", vec2f{1, 1}), true);
    } else if (type == "cylinder") {
        shp = make_cylinder(js.value("steps", vec3i{64, 32, 16}),
            js.value("size", vec2f{2.0f, 2.0f}),
            js.value("uvsize", vec3f{1, 1, 1}), true);
    } else if (type == "cylinder_rounded") {
        shp = make_cylinder_rounded(js.value("steps", vec3i{64, 32, 16}),
            js.value("size", vec2f{2.0f, 2.0f}),
            js.value("uvsize", vec3f{1, 1, 1}), js.value("radius", 0.15f), true);
    } else if (type == "sphere_geodesic") {
        shp = make_geodesic_sphere(
            js.value("tesselation", 4), js.value("size", 2.0f), true);
    } else if (type == "floor") {
        shp = make_floor(js.value("steps", vec2i{1, 1}),
            js.value("size", vec2f{40, 40}), js.value("uvsize", vec2f{20, 20}),
            true);
    } else if (type == "matball") {
        shp = make_sphere(js.value("steps", vec2i{64, 32}),
            js.value("size", 2.0f), js.value("uvsize", vec2f{1, 1}), true);
    } else if (type == "hairball") {
        auto base = make_sphere_cube(32, js.value("size", 2.0f) * 0.8f, 1, true);
        shp = make_hair(js.value("steps", vec2i{4, 65536}), base.triangles,
            base.positions, base.normals, base.texturecoords,
            js.value("length", vec2f{0.2f, 0.2f}),
            js.value("radius", vec2f{0.001f, 0.001f}),
            js.value("noise", vec2f{0, 0}), js.value("clump", vec2f{0, 0}));
    } else if (type == "hairball_interior") {
        shp = make_sphere_cube(32, js.value("size", 2.0f) * 0.8f, 1, true);
    } else if (type == "suzanne") {
        shp = make_suzanne(js.value("size", 2.0f), true);
    } else {
        throw runtime_error("unknown shape type " + type);
    }
    if (js.value("flipyz", false)) {
        for (auto& p : shp.positions) p = {p.x, p.z, p.y};
        for (auto& n : shp.normals) n = {n.x, n.z, n.y};
    }
    val->points        = shp.points;
    val->lines         = shp.lines;
    val->triangles     = shp.triangles;
    val->positions     = shp.positions;
    val->normals       = shp.normals;
    val->texturecoords = shp.texturecoords;
    val->radius        = shp.radius;
    return true;
}

// Serialize struct
bool parse_json_object(const json& js, yocto_shape* val, const yocto_scene* scn) {
    static const auto def = yocto_shape();
    if (!parse_json_objbegin(js)) return false;
    if (!parse_json_value(js, val->name, "name", def.name)) return false;
    if (!parse_json_value(js, val->filename, "filename", def.filename))
        return false;
    if (!parse_json_value(js, val->points, "points", def.points)) return false;
    if (!parse_json_value(js, val->lines, "lines", def.lines)) return false;
    if (!parse_json_value(js, val->triangles, "triangles", def.triangles))
        return false;
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
            js, val->tangent_spaces, "tangent_spaces", def.tangent_spaces))
        return false;
    if (!parse_json_procedural(js, val, "!!proc", scn)) return false;
    return true;
}

// Serialize struct
bool dump_json_object(json& js, const yocto_surface* val, const yocto_scene* scn) {
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
    if (val->compute_normals != def.compute_normals)
        js["compute_normals"] = val->compute_normals;
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
    const json& js, yocto_surface* val, const yocto_scene* scn) {
    if (!parse_json_objbegin(js)) return false;
    auto type = js.value("type", ""s);
    if (type == "") return true;
    auto shp = make_fvshape_data();
    if (type == "cube") {
        shp = make_fvcube(js.value("steps", vec3i{1, 1, 1}),
            js.value("size", vec3f{2, 2, 2}), js.value("uvsize", vec3f{1, 1, 1}));
    } else if (type == "cube_open") {
        shp = make_fvcube(js.value("steps", vec3i{1, 1, 1}),
            js.value("size", vec3f{2, 2, 2}), js.value("uvsize", vec3f{1, 1, 1}));
        shp.positions_quads.pop_back();
        shp.normals_quads.pop_back();
        shp.quads_texcoord.pop_back();
    } else if (type == "suzanne") {
        auto qshp           = make_suzanne(js.value("size", 2.0f), false);
        shp.positions_quads = qshp.quads;
        shp.positions       = qshp.positions;
    } else {
        throw runtime_error("unknown shape type " + type);
    }
    val->positions_quads     = shp.positions_quads;
    val->positions           = shp.positions;
    val->texturecoords_quads = shp.quads_texcoord;
    val->texturecoords       = shp.texturecoords;
    if (val->filename == "") val->filename = "meshes/" + val->name + ".obj";
    return true;
}

// Serialize struct
bool parse_json_object(
    const json& js, yocto_surface* val, const yocto_scene* scn) {
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
    if (!parse_json_value(
            js, val->compute_normals, "compute_normals", def.compute_normals))
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
    if (!parse_json_procedural(js, val, "!!proc", scn)) return false;
    return true;
}

// Serialize struct
bool dump_json_object(
    json& js, const yocto_instance* val, const yocto_scene* scn) {
    static const auto def = yocto_instance();
    if (!dump_json_objbegin(js)) return false;
    if (!dump_json_value(js, val->name, "name", def.name)) return false;
    if (!dump_json_value(js, val->frame, "frame", def.frame)) return false;
    if (!dump_json_objref(js, val->shape, "shape", scn->shapes)) return false;
    if (!dump_json_objref(js, val->material, "material", scn->materials))
        return false;
    if (!dump_json_objref(js, val->surface, "surface", scn->surfaces))
        return false;
    return true;
}

// Procedural commands for instances
bool apply_json_procedural(
    const json& js, yocto_instance* val, const yocto_scene* scn) {
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
    const json& js, yocto_instance* val, const yocto_scene* scn) {
    static const auto def = yocto_instance();
    if (!parse_json_objbegin(js)) return false;
    if (!parse_json_value(js, val->name, "name", def.name)) return false;
    if (!parse_json_value(js, val->frame, "frame", def.frame)) return false;
    if (!parse_json_objref(js, val->shape, "shape", scn->shapes)) return false;
    if (!parse_json_objref(js, val->surface, "surface", scn->surfaces))
        return false;
    if (!parse_json_objref(js, val->material, "material", scn->materials))
        return false;
    if (!parse_json_procedural(js, val, "!!proc", scn)) return false;
    return true;
}

// Serialize struct
bool dump_json_object(
    json& js, const yocto_environment* val, const yocto_scene* scn) {
    static const auto def = yocto_environment();
    if (!dump_json_objbegin(js)) return false;
    if (!dump_json_value(js, val->name, "name", def.name)) return false;
    if (!dump_json_value(js, val->frame, "frame", def.frame)) return false;
    if (!dump_json_value(js, val->emission, "emission", def.emission))
        return false;
    if (!dump_json_objref(
            js, val->emission_texture, "emission_txt", scn->textures))
        return false;
    return true;
}

// Procedural commands for materials
bool apply_json_procedural(
    const json& js, yocto_environment* val, const yocto_scene* scn) {
    if (!parse_json_objbegin(js)) return false;
    if (js.count("rotation")) {
        auto rotation = js.value("rotation", zero4f);
        val->frame    = rotation_frame(xyz(rotation), rotation.w);
    }
    return true;
}

// Serialize struct
bool parse_json_object(
    const json& js, yocto_environment* val, const yocto_scene* scn) {
    static const auto def = yocto_environment();
    if (!parse_json_objbegin(js)) return false;
    if (!parse_json_value(js, val->name, "name", def.name)) return false;
    if (!parse_json_value(js, val->frame, "frame", def.frame)) return false;
    if (!parse_json_value(js, val->emission, "emission", def.emission))
        return false;
    if (!parse_json_objref(
            js, val->emission_texture, "emission_txt", scn->textures))
        return false;
    if (!parse_json_procedural(js, val, "!!proc", scn)) return false;
    return true;
}

// Serialize struct
bool dump_json_object(
    json& js, const yocto_scene_node* val, const yocto_scene* scn) {
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
    if (!dump_json_objref(js, val->parent, "parent", scn->nodes)) return false;
    if (!dump_json_objref(js, val->camera, "camera", scn->cameras))
        return false;
    if (!dump_json_objref(js, val->instance, "instance", scn->instances))
        return false;
    if (!dump_json_objref(js, val->environment, "environment", scn->environments))
        return false;
    return true;
}

// Procedural commands for nodes
bool apply_json_procedural(
    const json& js, yocto_scene_node* val, const yocto_scene* scn) {
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
    const json& js, yocto_scene_node* val, const yocto_scene* scn) {
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
    if (!parse_json_objref(js, val->parent, "parent", scn->nodes)) return false;
    if (!parse_json_objref(js, val->instance, "instance", scn->instances))
        return false;
    if (!parse_json_objref(js, val->camera, "camera", scn->cameras))
        return false;
    if (!parse_json_objref(js, val->environment, "environment", scn->environments))
        return false;
    if (!parse_json_procedural(js, val, "!!proc", scn)) return false;
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
    } catch (...) { return false; }
    return true;
}

// Serialize struct
bool dump_json_object(
    json& js, const yocto_animation* val, const yocto_scene* scn) {
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
    if (!dump_json_objref(js, val->node_targets, "node_targets", scn->nodes))
        return false;
    return true;
}

// Procedural commands for animations
bool apply_json_procedural(
    const json& js, yocto_animation* val, const yocto_scene* scn) {
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
    const json& js, yocto_animation* val, const yocto_scene* scn) {
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
    if (!parse_json_objref(js, val->node_targets, "node_targets", scn->nodes))
        return false;
    if (!parse_json_procedural(js, val, "!!proc", scn)) return false;
    return true;
}

// Serialize struct
bool dump_json_object(json& js, const yocto_scene* val, const yocto_scene* scn) {
    static const auto def = yocto_scene();
    if (!dump_json_objbegin(js)) return false;
    if (!dump_json_value(js, val->name, "name", def.name)) return false;
    if (!dump_json_objarray(js, val->cameras, "cameras", scn)) return false;
    if (!dump_json_objarray(js, val->textures, "textures", scn)) return false;
    if (!dump_json_objarray(js, val->materials, "materials", scn)) return false;
    if (!dump_json_objarray(js, val->shapes, "shapes", scn)) return false;
    if (!dump_json_objarray(js, val->surfaces, "surfaces", scn)) return false;
    if (!dump_json_objarray(js, val->instances, "instances", scn)) return false;
    if (!dump_json_objarray(js, val->environments, "environments", scn))
        return false;
    if (!dump_json_objarray(js, val->nodes, "nodes", scn)) return false;
    if (!dump_json_objarray(js, val->animations, "animations", scn))
        return false;
    return true;
}

// Procedural commands for scenes
bool apply_json_procedural(
    const json& js, yocto_scene* val, const yocto_scene* scn) {
    if (!parse_json_objbegin(js)) return false;
    if (js.count("random_instances")) {
        auto& jjs  = js.at("random_instances");
        auto  num  = jjs.value("num", 100);
        auto  seed = jjs.value("seed", 13);
        auto  base = make_unique<yocto_instance>();
        parse_json_object(jjs.at("base"), base.get(), scn);
        auto ists = vector<unique_ptr<yocto_instance>>();
        for (auto& j : jjs.at("instances")) {
            ists.push_back(make_unique<yocto_instance>());
            parse_json_object(j, ists.back().get(), scn);
        }

        auto pos                 = vector<vec3f>();
        auto norm                = vector<vec3f>();
        auto texcoord            = vector<vec2f>();
        tie(pos, norm, texcoord) = sample_triangles_points(
            base->shape->triangles, base->shape->positions,
            base->shape->normals, base->shape->texturecoords, num, seed);

        auto nmap = unordered_map<yocto_instance*, int>();
        for (auto& ist : ists) nmap[ist.get()] = 0;
        auto rng = make_rng(seed, 17);
        for (auto i = 0; i < num; i++) {
            auto ist = ists.at(rand1i(rng, (int)ists.size() - 1)).get();
            nmap[ist] += 1;
            val->instances.push_back(new yocto_instance());
            val->instances.back()->name = ist->name + std::to_string(nmap[ist]);
            val->instances.back()->frame = base->frame *
                                           translation_frame(pos[i]) * ist->frame;
            val->instances.back()->shape    = ist->shape;
            val->instances.back()->material = ist->material;
            val->instances.back()->surface  = ist->surface;
        }
    }
    return true;
}

// Json to scene
bool parse_json_object(const json& js, yocto_scene* val, const yocto_scene* scn) {
    static const auto def = yocto_scene();
    if (!parse_json_objbegin(js)) return false;
    if (!parse_json_value(js, val->name, "name", def.name)) return false;
    if (!parse_json_objarray(js, val->cameras, "cameras", scn)) return false;
    if (!parse_json_objarray(js, val->textures, "textures", scn)) return false;
    if (!parse_json_objarray(js, val->voltextures, "voltextures", scn))
        return false;
    if (!parse_json_objarray(js, val->materials, "materials", scn))
        return false;
    if (!parse_json_objarray(js, val->shapes, "shapes", scn)) return false;
    if (!parse_json_objarray(js, val->surfaces, "surfaces", scn)) return false;
    if (!parse_json_objarray(js, val->instances, "instances", scn))
        return false;
    if (!parse_json_objarray(js, val->environments, "environments", scn))
        return false;
    if (!parse_json_objarray(js, val->nodes, "nodes", scn)) return false;
    if (!parse_json_objarray(js, val->animations, "animations", scn))
        return false;
    if (!parse_json_procedural(js, val, "!!proc", scn)) return false;
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
    // initialize
    auto scn = make_unique<yocto_scene>();

    // load jsonz
    auto js = load_json(filename);
    if (js.empty()) return nullptr;

    // deserialize json
    try {
        if (!parse_json_object(js, scn.get())) {
            log_io_error("could not deserialize json {}", filename);
            return nullptr;
        }
    } catch (...) {
        log_io_error("could not deserialize json {}", filename);
        return nullptr;
    }

    // load meshes
    auto dirname = get_dirname(filename);
    for (auto shp : scn->shapes) {
        if (shp->filename == "" || !shp->positions.empty()) continue;
        auto filename = normalize_path(dirname + "/" + shp->filename);
        if (!load_mesh(filename, shp->points, shp->lines, shp->triangles,
                shp->positions, shp->normals, shp->texturecoords, shp->colors,
                shp->radius)) {
            if (!skip_missing) return nullptr;
        }
    }

    // load suddivs
    for (auto sbd : scn->surfaces) {
        if (sbd->filename == "" || !sbd->positions.empty()) continue;
        auto filename   = normalize_path(dirname + "/" + sbd->filename);
        auto quads_norm = vector<vec4i>();
        auto norm       = vector<vec3f>();
        if (!load_fvmesh(filename, sbd->positions_quads, sbd->positions,
                quads_norm, norm, sbd->texturecoords_quads, sbd->texturecoords,
                sbd->colors_quads, sbd->colors)) {
            if (!skip_missing) return nullptr;
        }
    }

    // skip textures
    if (load_textures) {
        if (!load_scene_textures(scn.get(), dirname, skip_missing, false))
            return nullptr;
    }

    // fix scene
    if (scn->name == "") scn->name = get_filename(filename);
    add_missing_cameras(scn.get());
    add_missing_materials(scn.get());
    add_missing_names(scn.get());
    update_transforms(scn.get());

    // done
    return scn.release();
}

// Save a scene in the builtin JSON format.
bool save_json_scene(const string& filename, const yocto_scene* scn,
    bool save_textures, bool skip_missing) {
    // save json
    auto js = json();
    try {
        if (!dump_json_object(js, scn)) {
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
    for (auto& shp : scn->shapes) {
        if (shp->filename == "") continue;
        auto filename = normalize_path(dirname + "/" + shp->filename);
        if (!save_mesh(filename, shp->points, shp->lines, shp->triangles,
                shp->positions, shp->normals, shp->texturecoords, shp->colors,
                shp->radius)) {
            if (!skip_missing) return false;
        }
    }

    // save subdivs
    for (auto& sbd : scn->surfaces) {
        if (sbd->filename == "") continue;
        auto filename = normalize_path(dirname + "/" + sbd->filename);
        if (!save_fvmesh(filename, sbd->positions_quads, sbd->positions, {}, {},
                sbd->texturecoords_quads, sbd->texturecoords, sbd->colors_quads,
                sbd->colors)) {
            if (!skip_missing) return false;
        }
    }

    // skip textures
    if (save_textures) {
        if (!save_scene_textures(scn, dirname, skip_missing)) return false;
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
    return a.pos == b.pos && a.texcoord == b.texcoord && a.norm == b.norm;
}

struct obj_vertex_hash {
    size_t operator()(const obj_vertex& v) const {
        auto vh = std::hash<int>();
        auto h  = (size_t)0;
        for (auto i = 0; i < 3; i++)
            h ^= vh((&v.pos)[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
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
    auto val = obj_vertex{0, 0, 0};
    val.pos  = parse_int(s);
    if (*s == '/') {
        s++;
        if (*s == '/') {
            s++;
            val.norm = parse_int(s);
        } else {
            val.texcoord = parse_int(s);
            if (*s == '/') {
                s++;
                val.norm = parse_int(s);
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
            auto cam     = obj_camera();
            cam.name     = parse_string(ss);
            cam.ortho    = parse_bool(ss);
            cam.film     = parse_vec2f(ss);
            cam.focal    = parse_float(ss);
            cam.focus    = parse_float(ss);
            cam.aperture = parse_float(ss);
            cam.frame    = parse_frame3f(ss);
            if (cb.camera) cb.camera(cam);
        } else if (cmd == "e") {
            auto env        = obj_environment();
            env.name        = parse_string(ss);
            env.ke          = parse_vec3f(ss);
            env.ke_txt.path = parse_string(ss);
            if (env.ke_txt.path == "\"\"") env.ke_txt.path = "";
            if (cb.environmnet) cb.environmnet(env);
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
            vert_size.pos += 1;
        } else if (cmd == "vn") {
            if (cb.norm) cb.norm(parse_vec3f(ss));
            vert_size.norm += 1;
        } else if (cmd == "vt") {
            auto v = parse_vec2f(ss);
            if (flip_texcoord) v.y = 1 - v.y;
            if (cb.texcoord) cb.texcoord(v);
            vert_size.texcoord += 1;
        } else if (cmd == "f" || cmd == "l" || cmd == "p") {
            verts.clear();
            while (true) {
                auto vert = parse_obj_vertex(ss);
                if (!vert.pos) break;
                if (vert.pos < 0) vert.pos = vert_size.pos + vert.pos + 1;
                if (vert.texcoord < 0)
                    vert.texcoord = vert_size.texcoord + vert.texcoord + 1;
                if (vert.norm < 0) vert.norm = vert_size.norm + vert.norm + 1;
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
    auto scn = make_unique<yocto_scene>();

    // splitting policy
    auto split_material  = split_shapes;
    auto split_group     = split_shapes;
    auto split_smoothing = split_shapes;

    // current parsing values
    auto matname   = string();
    auto oname     = string();
    auto gname     = string();
    auto smoothing = true;
    auto ist       = (yocto_instance*)nullptr;

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
    auto is_instance_empty = [](yocto_instance* ist) {
        if (ist->surface) {
            return ist->surface->positions.empty();
        } else if (ist->shape) {
            return ist->shape->positions.empty();
        } else {
            return true;
        }
    };
    auto add_instance = [&](yocto_scene* scn, const string& objname,
                            const string& matname, const string& groupname,
                            bool smoothing) {
        if (scn->instances.empty() || objname != scn->instances.back()->name ||
            !is_instance_empty(scn->instances.back())) {
            auto ist = new yocto_instance();
            scn->instances.push_back(ist);
            ist->shape = new yocto_shape();
            scn->shapes.push_back(ist->shape);
        }
        name_map[objname] += 1;
        auto name = (name_map[objname] == 1) ?
                        objname :
                        objname + "_" + std::to_string(name_map[objname] - 1);
        if (objname == "") name = "object" + name;
        auto ist  = scn->instances.back();
        ist->name = name;
        if (ist->shape) ist->shape->name = ist->name;
        if (ist->surface) ist->surface->name = ist->name;
        if (matname != "") {
            auto it = mmap.find(matname);
            if (it == mmap.end())
                throw runtime_error("missing material " + matname);
            ist->material = it->second;
        }
        vert_map.clear();
        pos_map.clear();
        norm_map.clear();
        texcoord_map.clear();
        return ist;
    };
    // Parse texture options and name
    auto add_texture = [&scn, &tmap](
                           const obj_texture_info& info, bool force_linear) {
        if (info.path == "") return (yocto_texture*)nullptr;
        if (tmap.find(info.path) != tmap.end()) { return tmap.at(info.path); }

        // create texture
        auto txt           = new yocto_texture();
        txt->name          = info.path;
        txt->filename      = info.path;
        txt->clamp_to_edge = info.clamp;
        txt->height_scale  = info.scale;
        txt->ldr_as_linear = force_linear || is_hdr_filename(info.path);
        scn->textures.push_back(txt);
        tmap[info.path] = txt;

        return txt;
    };
    // Parse texture options and name
    auto add_voltexture = [&scn, &vmap](const obj_texture_info& info, bool srgb) {
        if (info.path == "") return (yocto_voltexture*)nullptr;
        if (vmap.find(info.path) != vmap.end()) { return vmap.at(info.path); }

        // create texture
        auto txt      = new yocto_voltexture();
        txt->name     = info.path;
        txt->filename = info.path;
        scn->voltextures.push_back(txt);
        vmap[info.path] = txt;

        return txt;
    };
    // Add  vertices to the current shape
    auto add_verts = [&](const vector<obj_vertex>& verts) {
        for (auto& vert : verts) {
            auto it = vert_map.find(vert);
            if (it != vert_map.end()) continue;
            auto nverts = (int)ist->shape->positions.size();
            vert_map.insert(it, {vert, nverts});
            if (vert.pos)
                ist->shape->positions.push_back(opos.at(vert.pos - 1));
            if (vert.texcoord)
                ist->shape->texturecoords.push_back(
                    otexcoord.at(vert.texcoord - 1));
            if (vert.norm)
                ist->shape->normals.push_back(onorm.at(vert.norm - 1));
        }
    };

    // current objet
    ist = add_instance(scn.get(), "", "", "", true);

    // callbacks
    auto cb     = obj_callbacks();
    cb.vert     = [&](vec3f v) { opos.push_back(v); };
    cb.norm     = [&](vec3f v) { onorm.push_back(v); };
    cb.texcoord = [&](vec2f v) { otexcoord.push_back(v); };
    cb.face     = [&](const vector<obj_vertex>& verts) {
        add_verts(verts);
        for (auto i = 2; i < verts.size(); i++)
            ist->shape->triangles.push_back({vert_map.at(verts[0]),
                vert_map.at(verts[i - 1]), vert_map.at(verts[i])});
    };
    cb.line = [&](const vector<obj_vertex>& verts) {
        add_verts(verts);
        for (auto i = 1; i < verts.size(); i++)
            ist->shape->lines.push_back(
                {vert_map.at(verts[i - 1]), vert_map.at(verts[i])});
    };
    cb.point = [&](const vector<obj_vertex>& verts) {
        add_verts(verts);
        for (auto i = 0; i < verts.size(); i++)
            ist->shape->points.push_back(vert_map.at(verts[i]));
    };
    cb.object = [&](const string& name) {
        oname     = name;
        gname     = "";
        matname   = "";
        smoothing = true;
        ist       = add_instance(scn.get(), oname, matname, gname, smoothing);
    };
    cb.group = [&](const string& name) {
        gname = name;
        if (split_group) {
            ist = add_instance(scn.get(), oname, matname, gname, smoothing);
        }
    };
    cb.smoothing = [&](const string& name) {
        smoothing = (name == "on");
        if (split_smoothing) {
            ist = add_instance(scn.get(), oname, matname, gname, smoothing);
        }
    };
    cb.usemtl = [&](const string& name) {
        matname = name;
        if (split_material) {
            ist = add_instance(scn.get(), oname, matname, gname, smoothing);
        } else {
            if (matname != "") ist->material = mmap.at(matname);
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
        scn->materials.push_back(mat);
        mmap[mat->name] = mat;
    };
    cb.camera = [&](const obj_camera& ocam) {
        auto cam            = new yocto_camera();
        cam->name           = ocam.name;
        cam->orthographic   = ocam.ortho;
        cam->film_size      = ocam.film;
        cam->focal_length   = ocam.focal;
        cam->focus_distance = ocam.focus;
        cam->lens_aperture  = ocam.aperture;
        cam->frame          = ocam.frame;
        scn->cameras.push_back(cam);
    };
    cb.environmnet = [&](const obj_environment& oenv) {
        auto env              = new yocto_environment();
        env->name             = oenv.name;
        env->emission         = oenv.ke;
        env->emission_texture = add_texture(oenv.ke_txt, true);
        scn->environments.push_back(env);
    };

    // Parse obj
    if (!load_obj(filename, cb, false, skip_missing)) return nullptr;

    // cleanup empty
    // TODO: delete unused
    for (auto idx = 0; idx < scn->instances.size(); idx++) {
        if (!is_instance_empty(scn->instances[idx])) continue;
        auto ist = scn->instances[idx];
        if (ist->shape) {
            scn->shapes.erase(
                std::find(scn->shapes.begin(), scn->shapes.end(), ist->shape));
        }
        if (ist->surface) {
            scn->surfaces.erase(std::find(
                scn->surfaces.begin(), scn->surfaces.end(), ist->surface));
        }
        scn->instances.erase(scn->instances.begin() + idx);
        idx--;
    }

    // load textures
    auto dirname = get_dirname(filename);
    if (load_textures) {
        if (!load_scene_textures(scn.get(), dirname, skip_missing, false))
            return nullptr;
    }

    // fix scene
    scn->name = get_filename(filename);
    add_missing_cameras(scn.get());
    add_missing_materials(scn.get());
    add_missing_names(scn.get());
    update_transforms(scn.get());

    // done
    return scn.release();
}

bool save_mtl(
    const string& filename, const yocto_scene* scn, bool flip_tr = true) {
    // open
    auto fs = open(filename, "wt");
    if (!fs) return false;

    // for each material, dump all the values
    for (auto mat : scn->materials) {
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

bool save_objx(const string& filename, const yocto_scene* scn) {
    // scene
    auto fs = open(filename, "wt");
    if (!fs) return false;

    // cameras
    for (auto cam : scn->cameras) {
        print(fs, "c {} {} {} {} {} {} {}\n", cam->name, (int)cam->orthographic,
            cam->film_size, cam->focal_length, cam->focus_distance,
            cam->lens_aperture, cam->frame);
    }

    // environments
    for (auto env : scn->environments) {
        print(fs, "e {} {} {} {}\n", env->name.c_str(), env->emission,
            ((env->emission_texture) ? env->emission_texture->filename.c_str() :
                                       "\"\""),
            env->frame);
    }

    // done
    return true;
}

string to_string(const obj_vertex& v) {
    auto s = std::to_string(v.pos);
    if (v.texcoord) {
        s += "/" + std::to_string(v.texcoord);
        if (v.norm) s += "/" + std::to_string(v.norm);
    } else {
        if (v.norm) s += "//" + std::to_string(v.norm);
    }
    return s;
}

bool save_obj(
    const string& filename, const yocto_scene* scn, bool flip_texcoord = true) {
    // scene
    auto fs = open(filename, "wt");
    if (!fs) return false;

    // material library
    if (!scn->materials.empty()) {
        auto mtlname = replace_extension(get_filename(filename), "mtl");
        print(fs, "mtllib {}\n", mtlname);
    }

    // shapes
    auto offset = obj_vertex{0, 0, 0};
    for (auto ist : scn->instances) {
        if (!ist->surface) {
            print(fs, "o {}\n", ist->name);
            if (ist->material) print(fs, "usemtl {}\n", ist->material->name);
            if (ist->frame == identity_frame3f) {
                for (auto& p : ist->shape->positions) print(fs, "v {}\n", p);
                for (auto& n : ist->shape->normals) print(fs, "vn {}\n", n);
                for (auto& t : ist->shape->texturecoords)
                    print(fs, "vt {}\n",
                        vec2f{t.x, (flip_texcoord) ? 1 - t.y : t.y});
            } else {
                for (auto& pp : ist->shape->positions) {
                    print(fs, "v {}\n", transform_point(ist->frame, pp));
                }
                for (auto& nn : ist->shape->normals) {
                    print(fs, "vn {}\n", transform_direction(ist->frame, nn));
                }
                for (auto& t : ist->shape->texturecoords)
                    print(fs, "vt {}\n",
                        vec2f{t.x, (flip_texcoord) ? 1 - t.y : t.y});
            }
            auto mask = obj_vertex{1, ist->shape->texturecoords.empty() ? 0 : 1,
                ist->shape->normals.empty() ? 0 : 1};
            auto vert = [mask, offset](int i) {
                return obj_vertex{(i + offset.pos + 1) * mask.pos,
                    (i + offset.texcoord + 1) * mask.texcoord,
                    (i + offset.norm + 1) * mask.norm};
            };
            for (auto& t : ist->shape->triangles) {
                print(fs, "f {} {} {}\n", to_string(vert(t.x)),
                    to_string(vert(t.y)), to_string(vert(t.z)));
            }
            for (auto& l : ist->shape->lines) {
                print(fs, "l {} {}\n", to_string(vert(l.x)),
                    to_string(vert(l.y)));
            }
            offset.pos += ist->shape->positions.size();
            offset.texcoord += ist->shape->texturecoords.size();
            offset.norm += ist->shape->normals.size();
        } else {
            print(fs, "o {}\n", ist->name);
            if (ist->material) print(fs, "usemtl {}\n", ist->material->name);
            if (ist->frame == identity_frame3f) {
                for (auto& p : ist->surface->positions) print(fs, "v {}\n", p);
                for (auto& t : ist->surface->texturecoords)
                    print(fs, "vt {}\n",
                        vec2f{t.x, (flip_texcoord) ? 1 - t.y : t.y});
            } else {
                for (auto& pp : ist->surface->positions) {
                    auto p = transform_point(ist->frame, pp);
                    print(fs, "v {}\n", p);
                }
                for (auto& t : ist->surface->texturecoords)
                    print(fs, "vt {}\n",
                        vec2f{t.x, (flip_texcoord) ? 1 - t.y : t.y});
            }
            if (!ist->surface->texturecoords.empty()) {
                auto vert = [offset](int ip, int it) {
                    return obj_vertex{
                        ip + offset.pos + 1, it + offset.texcoord + 1, 0};
                };
                for (auto i = 0; i < ist->surface->positions_quads.size(); i++) {
                    auto qp = ist->surface->positions_quads[i];
                    auto qt = ist->surface->texturecoords_quads[i];
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
                    return obj_vertex{ip + offset.pos + 1, 0, 0};
                };
                for (auto& q : ist->surface->positions_quads) {
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
            offset.pos += ist->surface->positions.size();
            offset.texcoord += ist->surface->texturecoords.size();
        }
    }

    return true;
}

bool save_obj_scene(const string& filename, const yocto_scene* scn,
    bool save_textures, bool skip_missing) {
    if (!save_obj(filename, scn, true)) return false;
    if (!scn->materials.empty()) {
        if (!save_mtl(replace_extension(filename, ".mtl"), scn, true))
            return false;
    }
    if (!scn->cameras.empty() || !scn->environments.empty()) {
        if (!save_objx(replace_extension(filename, ".objx"), scn)) return false;
    }

    // skip textures if needed
    auto dirname = get_dirname(filename);
    if (save_textures) {
        if (!save_scene_textures(scn, dirname, skip_missing)) return false;
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
bool gltf_to_scene(yocto_scene* scn, const json& gltf, const string& dirname) {
    // convert textures
    if (gltf.count("images")) {
        for (auto iid = 0; iid < gltf.at("images").size(); iid++) {
            auto& gimg    = gltf.at("images").at(iid);
            auto  txt     = new yocto_texture();
            txt->name     = gimg.value("name", ""s);
            txt->filename = (startswith(gimg.value("uri", ""s), "data:")) ?
                                string("[glTF-inline].png") :
                                gimg.value("uri", ""s);
            scn->textures.push_back(txt);
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
                if (pos == uri.npos) { return false; }
                // decode
                auto data_char = base64_decode(uri.substr(pos + 1));
                data = vector<unsigned char>((unsigned char*)data_char.c_str(),
                    (unsigned char*)data_char.c_str() + data_char.length());
            } else {
                auto filename = normalize_path(dirname + "/" + uri);
                data          = load_binary(filename);
                if (data.empty()) return false;
            }
            if (gbuf.value("byteLength", -1) != data.size()) { return false; }
        }
    }

    // add a texture
    auto add_texture = [scn, &gltf](const json& ginfo, bool force_linear) {
        if (!gltf.count("images") || !gltf.count("textures"))
            return (yocto_texture*)nullptr;
        if (ginfo.is_null() || ginfo.empty()) return (yocto_texture*)nullptr;
        if (ginfo.value("index", -1) < 0) return (yocto_texture*)nullptr;
        auto& gtxt = gltf.at("textures").at(ginfo.value("index", -1));
        if (gtxt.empty() || gtxt.value("source", -1) < 0)
            return (yocto_texture*)nullptr;
        auto txt = scn->textures.at(gtxt.value("source", -1));
        if (!gltf.count("samplers") || gtxt.value("sampler", -1) < 0)
            return txt;
        auto& gsmp         = gltf.at("samplers").at(gtxt.value("sampler", -1));
        txt->clamp_to_edge = gsmp.value("wrapS", ""s) == "ClampToEdge" ||
                             gsmp.value("wrapT", ""s) == "ClampToEdge";
        txt->height_scale = gsmp.value("scale", 1.0f) *
                            gsmp.value("strength", 1.0f);
        txt->ldr_as_linear = force_linear || is_hdr_filename(txt->filename);
        return txt;
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
            scn->materials.push_back(mat);
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
        if (compTypeNum == 5122 || compTypeNum == 5123) { compSize = 2; }
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
                auto shp  = new yocto_shape();
                shp->name = gmesh.value("name", ""s) +
                            ((sid) ? std::to_string(sid) : string());
                sid++;
                for (json::iterator gattr_it = gprim.at("attributes").begin();
                     gattr_it != gprim.at("attributes").end(); ++gattr_it) {
                    auto  semantic = gattr_it.key();
                    auto& gacc     = gltf.at("accessors")
                                     .at(gattr_it.value().get<int>());
                    auto vals = accessor_values(gacc);
                    if (semantic == "POSITION") {
                        shp->positions.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shp->positions.push_back({(float)vals[i][0],
                                (float)vals[i][1], (float)vals[i][2]});
                    } else if (semantic == "NORMAL") {
                        shp->normals.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shp->normals.push_back({(float)vals[i][0],
                                (float)vals[i][1], (float)vals[i][2]});
                    } else if (semantic == "TEXCOORD" ||
                               semantic == "TEXCOORD_0") {
                        shp->texturecoords.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shp->texturecoords.push_back(
                                {(float)vals[i][0], (float)vals[i][1]});
                    } else if (semantic == "COLOR" || semantic == "COLOR_0") {
                        shp->colors.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shp->colors.push_back(
                                {(float)vals[i][0], (float)vals[i][1],
                                    (float)vals[i][2], (float)vals[i][3]});
                    } else if (semantic == "TANGENT") {
                        shp->tangent_spaces.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shp->tangent_spaces.push_back(
                                {(float)vals[i][0], (float)vals[i][1],
                                    (float)vals[i][2], (float)vals[i][3]});
                        for (auto& t : shp->tangent_spaces) t.w = -t.w;
                    } else if (semantic == "RADIUS") {
                        shp->radius.reserve(vals.size());
                        for (auto i = 0; i < vals.size(); i++)
                            shp->radius.push_back((float)vals[i][0]);
                    } else {
                        // ignore
                    }
                }
                // indices
                auto mode = gprim.value("mode", 4);
                if (!gprim.count("indices")) {
                    if (mode == 4) {
                        // triangles
                        shp->triangles.reserve(shp->positions.size() / 3);
                        for (auto i = 0; i < shp->positions.size() / 3; i++)
                            shp->triangles.push_back(
                                {i * 3 + 0, i * 3 + 1, i * 3 + 2});
                    } else if (mode == 6) {
                        // triangle fan
                        shp->triangles.reserve(shp->positions.size() - 2);
                        for (auto i = 2; i < shp->positions.size(); i++)
                            shp->triangles.push_back({0, i - 1, i});
                    } else if (mode == 5) {
                        // triangle strip
                        shp->triangles.reserve(shp->positions.size() - 2);
                        for (auto i = 2; i < shp->positions.size(); i++)
                            shp->triangles.push_back({i - 2, i - 1, i});
                    } else if (mode == 1) {
                        // lines
                        shp->lines.reserve(shp->positions.size() / 2);
                        for (auto i = 0; i < shp->positions.size() / 2; i++)
                            shp->lines.push_back({i * 2 + 0, i * 2 + 1});
                    } else if (mode == 2) {
                        // line loop
                        shp->lines.reserve(shp->positions.size());
                        for (auto i = 1; i < shp->positions.size(); i++)
                            shp->lines.push_back({i - 1, i});
                        shp->lines.back() = {(int)shp->positions.size() - 1, 0};
                    } else if (mode == 3) {
                        // line strip
                        shp->lines.reserve(shp->positions.size() - 1);
                        for (auto i = 1; i < shp->positions.size(); i++)
                            shp->lines.push_back({i - 1, i});
                    } else if (mode == -1 || mode == 0) {
                        // points
                        printf("points not supported\n");
                    } else {
                        throw runtime_error("unknown primitive type");
                    }
                } else {
                    auto indices = accessor_values(
                        gltf.at("accessors").at(gprim.value("indices", -1)),
                        false);
                    if (mode == 4) {
                        // triangles
                        shp->triangles.reserve(indices.size());
                        for (auto i = 0; i < indices.size() / 3; i++)
                            shp->triangles.push_back({(int)indices[i * 3 + 0][0],
                                (int)indices[i * 3 + 1][0],
                                (int)indices[i * 3 + 2][0]});
                    } else if (mode == 6) {
                        // triangle fan
                        shp->triangles.reserve(indices.size() - 2);
                        for (auto i = 2; i < indices.size(); i++)
                            shp->triangles.push_back({(int)indices[0][0],
                                (int)indices[i - 1][0], (int)indices[i][0]});
                    } else if (mode == 5) {
                        // triangle strip
                        shp->triangles.reserve(indices.size() - 2);
                        for (auto i = 2; i < indices.size(); i++)
                            shp->triangles.push_back({(int)indices[i - 2][0],
                                (int)indices[i - 1][0], (int)indices[i][0]});
                    } else if (mode == 1) {
                        // lines
                        shp->lines.reserve(indices.size() / 2);
                        for (auto i = 0; i < indices.size() / 2; i++)
                            shp->lines.push_back({(int)indices[i * 2 + 0][0],
                                (int)indices[i * 2 + 1][0]});
                    } else if (mode == 2) {
                        // line loop
                        shp->lines.reserve(indices.size());
                        for (auto i = 1; i < indices.size(); i++)
                            shp->lines.push_back(
                                {(int)indices[i - 1][0], (int)indices[i][0]});
                        shp->lines.back() = {(int)indices[indices.size() - 1][0],
                            (int)indices[0][0]};
                    } else if (mode == 3) {
                        // line strip
                        shp->lines.reserve(indices.size() - 1);
                        for (auto i = 1; i < indices.size(); i++)
                            shp->lines.push_back(
                                {(int)indices[i - 1][0], (int)indices[i][0]});
                    } else if (mode == -1 || mode == 0) {
                        // points
                        printf("points not supported\n");
                    } else {
                        throw runtime_error("unknown primitive type");
                    }
                }
                auto mat = (gprim.count("material")) ?
                               scn->materials.at(gprim.value("material", -1)) :
                               nullptr;
                meshes.back().push_back({shp, mat});
                scn->shapes.push_back(shp);
            }
        }
    }

    // convert cameras
    if (gltf.count("cameras")) {
        for (auto cid = 0; cid < gltf.at("cameras").size(); cid++) {
            auto& gcam        = gltf.at("cameras").at(cid);
            auto  cam         = new yocto_camera();
            cam->name         = gcam.value("name", ""s);
            cam->orthographic = gcam.value("type", ""s) == "orthographic";
            if (cam->orthographic) {
                printf("orthographic not supported well\n");
                auto ortho = gcam.value("orthographic", json::object());
                set_camera_fovy(cam, ortho.value("ymag", 0.0f),
                    ortho.value("xmag", 0.0f) / ortho.value("ymag", 0.0f));
                cam->focus_distance = maxf;
                cam->lens_aperture  = 0;
            } else {
                auto persp = gcam.value("perspective", json::object());
                set_camera_fovy(cam, persp.value("yfov", 1.0f),
                    persp.value("aspectRatio", 1.0f));
                cam->focus_distance = maxf;
                cam->lens_aperture  = 0;
            }
            scn->cameras.push_back(cam);
        }
    }

    // convert nodes
    if (gltf.count("nodes")) {
        for (auto nid = 0; nid < gltf.at("nodes").size(); nid++) {
            auto& gnde = gltf.at("nodes").at(nid);
            auto  nde  = new yocto_scene_node();
            nde->name  = gnde.value("name", ""s);
            if (gnde.count("camera"))
                nde->camera = scn->cameras[gnde.value("camera", 0)];
            nde->translation = gnde.value("translation", zero3f);
            nde->rotation    = gnde.value("rotation", vec4f{0, 0, 0, 1});
            nde->scale       = gnde.value("scale", vec3f{1, 1, 1});
            nde->local = mat_to_frame(gnde.value("matrix", identity_mat4f));
            scn->nodes.push_back(nde);
        }

        // set up parent pointers
        for (auto nid = 0; nid < gltf.at("nodes").size(); nid++) {
            auto& gnde = gltf.at("nodes").at(nid);
            if (!gnde.count("children")) continue;
            auto nde = scn->nodes[nid];
            for (auto& cid : gnde.at("children"))
                scn->nodes[cid.get<int>()]->parent = nde;
        }

        // set up instances
        for (auto nid = 0; nid < gltf.at("nodes").size(); nid++) {
            auto& gnde = gltf.at("nodes").at(nid);
            if (!gnde.count("mesh")) continue;
            auto  nde  = scn->nodes[nid];
            auto& shps = meshes.at(gnde.value("mesh", 0));
            if (shps.empty()) continue;
            if (shps.size() == 1) {
                nde->instance           = new yocto_instance();
                nde->instance->name     = nde->name;
                nde->instance->shape    = get<0>(shps[0]);
                nde->instance->material = get<1>(shps[0]);
                scn->instances.push_back(nde->instance);
            } else {
                for (auto shp : shps) {
                    auto child            = new yocto_scene_node();
                    child->name           = nde->name + "_" + get<0>(shp)->name;
                    child->parent         = nde;
                    child->instance       = new yocto_instance();
                    child->instance->name = child->name;
                    child->instance->shape    = get<0>(shp);
                    child->instance->material = get<1>(shp);
                    scn->instances.push_back(child->instance);
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
                    auto anm  = new yocto_animation();
                    anm->name = (ganm.count("name") ? ganm.value("name", ""s) :
                                                      "anim") +
                                std::to_string(aid++);
                    anm->animation_group = ganm.value("name", ""s);
                    auto input_view      = accessor_values(
                        gltf.at("accessors").at(gsampler.value("input", -1)));
                    anm->keyframes_times.resize(input_view.size());
                    for (auto i = 0; i < input_view.size(); i++)
                        anm->keyframes_times[i] = input_view[i][0];
                    auto type = gsampler.value("interpolation", "LINEAR");
                    if (type == "LINEAR")
                        anm->interpolation_type = yocto_interpolation_type::linear;
                    if (type == "STEP")
                        anm->interpolation_type = yocto_interpolation_type::step;
                    if (type == "CUBICSPLINE")
                        anm->interpolation_type = yocto_interpolation_type::bezier;
                    auto output_view = accessor_values(
                        gltf.at("accessors").at(gsampler.value("output", -1)));
                    switch (path) {
                        case 0: {  // translation
                            anm->translation_keyframes.reserve(
                                output_view.size());
                            for (auto i = 0; i < output_view.size(); i++)
                                anm->translation_keyframes.push_back(
                                    {(float)output_view[i][0],
                                        (float)output_view[i][1],
                                        (float)output_view[i][2]});
                        } break;
                        case 1: {  // rotation
                            anm->rotation_keyframes.reserve(output_view.size());
                            for (auto i = 0; i < output_view.size(); i++)
                                anm->rotation_keyframes.push_back(
                                    {(float)output_view[i][0],
                                        (float)output_view[i][1],
                                        (float)output_view[i][2],
                                        (float)output_view[i][3]});
                        } break;
                        case 2: {  // scale
                            anm->scale_keyframes.reserve(output_view.size());
                            for (auto i = 0; i < output_view.size(); i++)
                                anm->scale_keyframes.push_back(
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
                            anm->weights.resize(values.size() / ncomp);
                            for (auto i = 0; i < anm->weights.size(); i++) {
                                anm->weights[i].resize(ncomp);
                                for (auto j = 0; j < ncomp; j++)
                                    anm->weights[i][j] = values[i * ncomp + j];
                            }
                        }
#endif
                        } break;
                        default: { return false; }
                    }
                    sampler_map[{gchannel.at("sampler").get<int>(),
                        path}] = (int)scn->animations.size();
                    scn->animations.push_back(anm);
                }
                scn->animations[sampler_map.at(
                                    {gchannel.at("sampler").get<int>(), path})]
                    ->node_targets.push_back(
                        scn->nodes[(int)gchannel.at("target").at("node").get<int>()]);
            }
        }
    }

    return true;
}

// Load a scene
yocto_scene* load_gltf_scene(
    const string& filename, bool load_textures, bool skip_missing) {
    // initialization
    auto scn = make_unique<yocto_scene>();

    // convert json
    auto js = load_json(filename);
    if (js.empty()) return nullptr;
    try {
        if (!gltf_to_scene(scn.get(), js, get_dirname(filename)))
            return nullptr;
    } catch (...) { return nullptr; }

    // load textures
    auto dirname = get_dirname(filename);
    if (load_textures) {
        if (!load_scene_textures(scn.get(), dirname, skip_missing, false))
            return nullptr;
    }

    // fix scene
    scn->name = get_filename(filename);
    add_missing_cameras(scn.get());
    add_missing_materials(scn.get());
    add_missing_names(scn.get());
    update_transforms(scn.get());

    // fix cameras
    auto bbox = compute_bbox(scn.get());
    for (auto cam : scn->cameras) {
        auto center = (bbox.min + bbox.max) / 2;
        auto dist   = dot(-cam->frame.z, center - cam->frame.o);
        if (dist > 0) cam->focus_distance = dist;
    }

    // done
    return scn.release();
}

// convert gltf scene to json
bool scene_to_gltf(const yocto_scene* scn, json& js) {
    // init to emprt object
    js = json::object();

    // start creating json
    js["asset"]["version"] = "2.0";

    // prepare top level nodes
    if (!scn->cameras.empty()) js["cameras"] = json::array();
    if (!scn->textures.empty()) {
        js["textures"] = json::array();
        js["images"]   = json::array();
    }
    if (!scn->materials.empty()) js["materials"] = json::array();
    if (!scn->shapes.empty()) {
        js["meshes"]      = json::array();
        js["buffers"]     = json::array();
        js["bufferViews"] = json::array();
        js["accessors"]   = json::array();
    }
    if (!scn->instances.empty()) js["nodes"] = json::array();
    if (!scn->nodes.empty()) js["nodes"] = json::array();

    // convert cameras
    auto cmap = unordered_map<yocto_camera*, int>();
    for (auto cam : scn->cameras) {
        auto cjs    = json();
        cjs["name"] = cam->name;
        if (!cam->orthographic) {
            cjs["type"]                       = "perspective";
            cjs["perspective"]["aspectRatio"] = cam->film_size.x /
                                                cam->film_size.y;
            cjs["perspective"]["znear"] = 0.01f;
        } else {
            cjs["type"]                  = "orthographic";
            cjs["orthographic"]["xmag"]  = cam->film_size.x / 2;
            cjs["orthographic"]["ymag"]  = cam->film_size.y / 2;
            cjs["orthographic"]["znear"] = 0.01f;
        }
        cmap[cam] = (int)js["cameras"].size();
        js["cameras"].push_back(cjs);
    }

    // textures
    auto tmap = unordered_map<yocto_texture*, int>();
    for (auto& txt : scn->textures) {
        auto tjs = json(), ijs = json();
        tjs["source"] = (int)js["images"].size();
        ijs["uri"]    = txt->filename;
        js["images"].push_back(ijs);
        js["textures"].push_back(tjs);
        tmap[txt] = (int)js["textures"].size() - 1;
    }

    // material
    auto mmap = unordered_map<yocto_material*, int>();
    for (auto mat : scn->materials) {
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
    for (auto ist : scn->instances)
        if (ist->material) shape_mats[ist->shape] = mmap.at(ist->material);

    // shapes
    auto smap = unordered_map<yocto_shape*, int>();
    for (auto shp : scn->shapes) {
        auto mjs = json(), bjs = json(), pjs = json();
        auto bid          = js["buffers"].size();
        mjs["name"]       = shp->name;
        mjs["primitives"] = json::array();
        bjs["name"]       = shp->name;
        bjs["byteLength"] = 0;
        bjs["uri"]        = replace_extension(shp->filename, ".bin");
        auto mat_it       = shape_mats.find(shp);
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
        auto nverts = (int)shp->positions.size();
        if (!shp->positions.empty())
            pjs["attributes"]["POSITION"] = add_accessor(nverts, "VEC3");
        if (!shp->normals.empty())
            pjs["attributes"]["NORMAL"] = add_accessor(nverts, "VEC3");
        if (!shp->texturecoords.empty())
            pjs["attributes"]["TEXCOORD_0"] = add_accessor(nverts, "VEC2");
        if (!shp->colors.empty())
            pjs["attributes"]["COLOR_0"] = add_accessor(nverts, "VEC4");
        if (!shp->radius.empty())
            pjs["attributes"]["RADIUS"] = add_accessor(nverts, "SCALAR");
        if (!shp->lines.empty()) {
            pjs["indices"] = add_accessor(
                (int)shp->lines.size() * 2, "SCALAR", true);
            pjs["mode"] = 1;
        }
        if (!shp->triangles.empty()) {
            pjs["indices"] = add_accessor(
                (int)shp->triangles.size() * 3, "SCALAR", true);
            pjs["mode"] = 4;
        }
        mjs["primitives"].push_back(pjs);
        js["meshes"].push_back(mjs);
        js["buffers"].push_back(bjs);
        smap[shp] = (int)js["meshes"].size() - 1;
    }

    // nodes
    auto nmap = unordered_map<yocto_scene_node*, int>();
    for (auto& nde : scn->nodes) {
        auto njs           = json();
        njs["name"]        = nde->name;
        njs["matrix"]      = frame_to_mat(nde->local);
        njs["translation"] = nde->translation;
        njs["rotation"]    = nde->rotation;
        njs["scale"]       = nde->scale;
        if (nde->camera) njs["camera"] = cmap.at(nde->camera);
        if (nde->instance) njs["mesh"] = smap.at(nde->instance->shape);
        if (!nde->children.empty()) {
            njs["children"] = json::array();
            for (auto& c : nde->children) njs["children"].push_back(nmap.at(c));
        }
        js["nodes"].push_back(njs);
        nmap[nde] = (int)js["nodes"].size() - 1;
    }

    // animations not supported yet
    if (!scn->animations.empty()) printf("animation not supported yet\n");

    // nodes from instances
    if (scn->nodes.empty()) {
        for (auto cam : scn->cameras) {
            auto njs      = json();
            njs["name"]   = cam->name;
            njs["camera"] = cmap.at(cam);
            njs["matrix"] = frame_to_mat(cam->frame);
            js["nodes"].push_back(njs);
        }
        for (auto ist : scn->instances) {
            auto njs      = json();
            njs["name"]   = ist->name;
            njs["mesh"]   = smap.at(ist->shape);
            njs["matrix"] = frame_to_mat(ist->frame);
            js["nodes"].push_back(njs);
        }
    }

    // done
    return true;
}

// save gltf mesh
bool save_gltf_mesh(const string& filename, const yocto_shape* shp) {
    auto fs = open(filename, "wb");
    if (!fs) return false;

    if (!write_values(fs, shp->positions)) return false;
    if (!write_values(fs, shp->normals)) return false;
    if (!write_values(fs, shp->texturecoords)) return false;
    if (!write_values(fs, shp->colors)) return false;
    if (!write_values(fs, shp->radius)) return false;
    if (!write_values(fs, shp->lines)) return false;
    if (!write_values(fs, shp->triangles)) return false;

    return true;
}

// Save gltf json
bool save_gltf_scene(const string& filename, const yocto_scene* scn,
    bool save_textures, bool skip_missing) {
    // save json
    auto js = json();
    try {
        if (!scene_to_gltf(scn, js)) return false;
    } catch (...) { return false; }
    if (!save_json(filename, js)) return false;

    // meshes
    auto dirname = get_dirname(filename);
    for (auto& shp : scn->shapes) {
        if (shp->filename == "") continue;
        auto filename = normalize_path(dirname + "/" + shp->filename);
        filename      = replace_extension(filename, ".bin");
        if (!save_gltf_mesh(filename, shp)) {
            if (!skip_missing) return false;
        }
    }

    // save textures
    if (save_textures) {
        if (!save_scene_textures(scn, dirname, skip_missing)) return false;
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
        if (tokens[i][0] != '"') throw runtime_error("string expected");
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
                if (!first && !list) throw runtime_error("bad params");
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
            if (js.at(name).size() == 1) { js.at(name) = js.at(name).at(0); }
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
        if (!is_cmd(tokens, i)) throw runtime_error("command expected");
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
            throw runtime_error("unsupported command " + tok);
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
    // convert to json
    auto js = json();
    try {
        if (!pbrt_to_json(filename, js)) return nullptr;
    } catch (...) { return nullptr; }

    auto dirname_ = get_dirname(filename);

    struct stack_item {
        frame3f         frame     = identity_frame3f;
        yocto_material* mat       = nullptr;
        yocto_material* light_mat = nullptr;
        float           focus = 1, aspect = 1;
        bool            reverse = false;
    };

    // parse
    auto scn   = make_unique<yocto_scene>();
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
            auto cam            = new yocto_camera();
            cam->name           = "cam" + std::to_string(cid++);
            cam->frame          = inverse(stack.back().frame);
            cam->frame.z        = -cam->frame.z;
            cam->focus_distance = stack.back().focus;
            auto aspect         = stack.back().aspect;
            auto fovy           = 1.0f;
            auto type           = jcmd.at("type").get<string>();
            if (type == "perspective") {
                fovy = jcmd.at("fov").get<float>() * pif / 180;
            } else {
                printf("%s camera not supported\n", type.c_str());
            }
            set_camera_fovy(cam, fovy, aspect);
            scn->cameras.push_back(cam);
        } else if (cmd == "Texture") {
            auto found = false;
            auto name  = jcmd.at("name").get<string>();
            for (auto& txt : scn->textures) {
                if (txt->name == name) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                auto txt = new yocto_texture();
                scn->textures.push_back(txt);
                txt->name          = jcmd.at("name").get<string>();
                txt_map[txt->name] = txt;
                auto type          = jcmd.at("type").get<string>();
                if (type == "imagemap") {
                    txt->filename = jcmd.at("filename").get<string>();
                    if (get_extension(txt->filename) == "pfm")
                        txt->filename = replace_extension(txt->filename,
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
                for (auto mat : scn->materials) {
                    if (mat->name == name) {
                        found = true;
                        break;
                    }
                }
            }
            if (!found) {
                auto mat = new yocto_material();
                scn->materials.push_back(mat);
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
                scn->materials.push_back(mat);
                stack.back().mat = mat;
            }
        } else if (cmd == "Shape") {
            auto shp  = new yocto_shape();
            auto type = jcmd.at("type").get<string>();
            if (type == "plymesh") {
                auto filename = jcmd.at("filename").get<string>();
                shp->name     = get_filename(filename);
                shp->filename = filename;
                if (!load_ply_mesh(dirname_ + "/" + filename, shp->points,
                        shp->lines, shp->triangles, shp->positions, shp->normals,
                        shp->texturecoords, shp->colors, shp->radius))
                    return nullptr;
            } else if (type == "trianglemesh") {
                shp->name     = "mesh" + std::to_string(sid++);
                shp->filename = "models/" + shp->name + ".ply";
                if (jcmd.count("indices"))
                    shp->triangles = get_vector_vec3i(jcmd.at("indices"));
                if (jcmd.count("P"))
                    shp->positions = get_vector_vec3f(jcmd.at("P"));
                if (jcmd.count("N"))
                    shp->normals = get_vector_vec3f(jcmd.at("N"));
                if (jcmd.count("uv"))
                    shp->texturecoords = get_vector_vec2f(jcmd.at("uv"));
            } else if (type == "sphere") {
                shp->name     = "sphere" + std::to_string(sid++);
                shp->filename = "models/" + shp->name + ".ply";
                auto radius   = 1.0f;
                if (jcmd.count("radius"))
                    radius = jcmd.at("radius").get<float>();
                auto sshp = make_sphere({64, 32}, 2 * radius, {1, 1}, true);
                shp->positions     = sshp.positions;
                shp->normals       = sshp.normals;
                shp->texturecoords = sshp.texturecoords;
                shp->triangles     = sshp.triangles;
            } else if (type == "disk") {
                shp->name     = "disk" + std::to_string(sid++);
                shp->filename = "models/" + shp->name + ".ply";
                auto radius   = 1.0f;
                if (jcmd.count("radius"))
                    radius = jcmd.at("radius").get<float>();
                auto sshp      = make_disk({32, 16}, 2 * radius, {1, 1}, true);
                shp->positions = sshp.positions;
                shp->normals   = sshp.normals;
                shp->texturecoords = sshp.texturecoords;
                shp->triangles     = sshp.triangles;
            } else {
                printf("%s shape not supported\n", type.c_str());
            }
            auto frame = stack.back().frame;
            auto scl = vec3f{length(frame.x), length(frame.y), length(frame.z)};
            for (auto& p : shp->positions) p *= scl;
            frame = {normalize(frame.x), normalize(frame.y), normalize(frame.z),
                frame.o};
            if (stack.back().reverse) {
                for (auto& t : shp->triangles) swap(t.y, t.z);
            }
            scn->shapes.push_back(shp);
            auto ist      = new yocto_instance();
            ist->name     = shp->name;
            ist->frame    = frame;
            ist->shape    = shp;
            ist->material = stack.back().mat;
            if (cur_object != "") {
                objects[cur_object].push_back(ist);
            } else {
                scn->instances.push_back(ist);
            }
        } else if (cmd == "ObjectInstance") {
            static auto instances = map<string, int>();
            auto        name      = jcmd.at("name").get<string>();
            auto&       object    = objects.at(name);
            for (auto shp : object) {
                instances[shp->name] += 1;
                auto ist  = new yocto_instance();
                ist->name = shp->name + "_ist" +
                            std::to_string(instances[shp->name]);
                ist->frame = stack.back().frame * shp->frame;
                ist->shape = shp->shape;
                scn->instances.push_back(ist);
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
                auto env  = new yocto_environment();
                env->name = "env" + std::to_string(lid++);
                // env->frame = frame3f{{1,0,0},{0,0,-1},{0,-1,0},{0,0,0}} *
                // stack.back().frame;
                env->frame = stack.back().frame * frame3f{{0, 0, 1}, {0, 1, 0},
                                                      {1, 0, 0}, {0, 0, 0}};
                env->emission = {1, 1, 1};
                if (jcmd.count("scale"))
                    env->emission *= get_vec3f(jcmd.at("scale"));
                if (jcmd.count("mapname")) {
                    auto txt      = new yocto_texture();
                    txt->filename = jcmd.at("mapname").get<string>();
                    txt->name     = env->name;
                    scn->textures.push_back(txt);
                    env->emission_texture = txt;
                }
                scn->environments.push_back(env);
            } else if (type == "distant") {
                auto distant_dist = 100;
                auto shp          = new yocto_shape();
                shp->name         = "distant" + std::to_string(lid++);
                auto from = vec3f{0, 0, 0}, to = vec3f{0, 0, 0};
                if (jcmd.count("from")) from = get_vec3f(jcmd.at("from"));
                if (jcmd.count("to")) to = get_vec3f(jcmd.at("to"));
                auto dir       = normalize(from - to);
                auto size      = distant_dist * sin(5 * pif / 180);
                auto sshp      = make_quad({1, 1}, {size, size}, {1, 1}, true);
                shp->positions = sshp.positions;
                shp->normals   = sshp.normals;
                shp->texturecoords = sshp.texturecoords;
                shp->triangles     = sshp.triangles;
                scn->shapes.push_back(shp);
                auto mat      = new yocto_material();
                mat->name     = shp->name;
                mat->emission = {1, 1, 1};
                if (jcmd.count("L")) mat->emission *= get_vec3f(jcmd.at("L"));
                if (jcmd.count("scale"))
                    mat->emission *= get_vec3f(jcmd.at("scale"));
                mat->emission *= (distant_dist * distant_dist) / (size * size);
                scn->materials.push_back(mat);
                auto ist      = new yocto_instance();
                ist->name     = shp->name;
                ist->shape    = shp;
                ist->material = mat;
                ist->frame = stack.back().frame * lookat_frame(dir * distant_dist,
                                                      zero3f, {0, 1, 0}, true);
                scn->instances.push_back(ist);
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
        for (auto cam : scn->cameras) {
            auto nde    = new yocto_scene_node();
            nde->name   = cam->name;
            nde->local  = cam->frame;
            nde->camera = cam;
            scn->nodes.insert(scn->nodes.begin(), nde);
        }
        for (auto env : scn->environments) {
            auto nde         = new yocto_scene_node();
            nde->name        = env->name;
            nde->local       = env->frame;
            nde->environment = env;
            scn->nodes.push_back(nde);
        }
    }

    // load textures
    auto dirname = get_dirname(filename);
    if (load_textures) {
        if (!load_scene_textures(scn.get(), dirname, skip_missing, false))
            return nullptr;
    }

    // fix scene
    scn->name = get_filename(filename);
    add_missing_cameras(scn.get());
    add_missing_materials(scn.get());
    add_missing_names(scn.get());
    update_transforms(scn.get());

    return scn.release();
}

// Convert a scene to pbrt format
bool save_pbrt(const string& filename, const yocto_scene* scn) {
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
    auto cam  = scn->cameras.front();
    auto from = cam->frame.o;
    auto to   = cam->frame.o - cam->frame.z;
    auto up   = cam->frame.y;
    print(fs, "LookAt {} {} {}\n", from, to, up);
    print(fs, "Camera \"perspective\" \"float fov\" {}\n",
        eval_camera_fovy(cam) * 180 / pif);

    // save renderer
    print(fs, "Sampler \"random\" \"integer pixelsamples\" [64]\n");
    // fprintf(f, "Sampler \"sobol\" \"interger pixelsamples\" [64]\n");
    print(fs, "Integrator \"path\"\n");
    print(fs,
        "Film \"image\" \"string filename\" [\"{}\"] "
        "\"integer xresolution\" [{}] \"integer yresolution\" [{}]\n",
        replace_extension(filename, "exr"), eval_image_size(cam, 512).x,
        eval_image_size(cam, 512).y);

    // start world
    print(fs, "WorldBegin\n");

    // convert textures
    for (auto txt : scn->textures) {
        print(fs,
            "Texture \"{}\" \"spectrum\" \"imagemap\" "
            "\"string filename\" [\"{}\"]\n",
            txt->name, txt->filename);
    }

    // convert materials
    for (auto mat : scn->materials) {
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
    for (auto ist : scn->instances) {
        print(fs, "AttributeBegin\n");
        print(fs, "TransformBegin\n");
        print(fs, "ConcatTransform [{}]\n", frame_to_mat(ist->frame));
        if (ist->material->emission != zero3f)
            print(fs, "AreaLightSource \"diffuse\" \"rgb L\" [ {} ]\n",
                ist->material->emission);
        print(fs, "NamedMaterial \"{}\"\n", ist->material->name);
        print(fs, "Shape \"plymesh\" \"string filename\" [\"{}\"]\n",
            ist->shape->filename.c_str());
        print(fs, "TransformEnd\n");
        print(fs, "AttributeEnd\n");
    }

    // end world
    print(fs, "WorldEnd\n");

    // done
    return true;
}

// Save a pbrt scene
bool save_pbrt_scene(const string& filename, const yocto_scene* scn,
    bool save_textures, bool skip_missing) {
    // save json
    if (!save_pbrt(filename, scn)) return false;

    // save meshes
    auto dirname = get_dirname(filename);
    for (auto& shp : scn->shapes) {
        if (shp->filename == "") continue;
        auto filename = normalize_path(dirname + "/" + shp->filename);
        if (!save_mesh(filename, shp->points, shp->lines, shp->triangles,
                shp->positions, shp->normals, shp->texturecoords, shp->colors,
                shp->radius)) {
            if (!skip_missing) return false;
        }
    }

    // skip textures
    if (save_textures) {
        if (!save_scene_textures(scn, dirname, skip_missing)) return false;
    }

    // done
    return true;
}

// Attempt to fix pbrt z-up.
void pbrt_flipyz_scene(const yocto_scene* scn) {
    // flip meshes
    for (auto shp : scn->shapes) {
        for (auto& p : shp->positions) swap(p.y, p.z);
        for (auto& n : shp->normals) swap(n.y, n.z);
    }
    for (auto ist : scn->instances) {
        ist->frame = ist->frame *
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
bool serialize_bin_value(string& vec, file_stream& fs, bool save) {
    if (save) {
        auto count = (size_t)vec.size();
        if (!write_value(fs, count)) return false;
        if (!write_values(fs, vec)) return false;
        return true;
    } else {
        auto count = (size_t)0;
        if (!read_value(fs, count)) return false;
        vec = string(count, ' ');
        if (!read_values(fs, vec)) return false;
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
    vector<T*>& vec, const yocto_scene* scn, file_stream& fs, bool save) {
    if (save) {
        auto count = (size_t)vec.size();
        if (!serialize_bin_value(count, fs, true)) return false;
        for (auto i = 0; i < vec.size(); ++i) {
            if (!serialize_bin_object(vec[i], scn, fs, true)) return false;
        }
        return true;
    } else {
        auto count = (size_t)0;
        if (!serialize_bin_value(count, fs, false)) return false;
        vec = vector<T*>(count);
        for (auto i = 0; i < vec.size(); ++i) {
            vec[i] = new T();
            if (!serialize_bin_object(vec[i], scn, fs, false)) return false;
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
bool serialize_bin_object(yocto_camera* cam, file_stream& fs, bool save) {
    if (!serialize_bin_value(cam->name, fs, save)) return false;
    if (!serialize_bin_value(cam->frame, fs, save)) return false;
    if (!serialize_bin_value(cam->orthographic, fs, save)) return false;
    if (!serialize_bin_value(cam->film_size, fs, save)) return false;
    if (!serialize_bin_value(cam->focal_length, fs, save)) return false;
    if (!serialize_bin_value(cam->focus_distance, fs, save)) return false;
    if (!serialize_bin_value(cam->lens_aperture, fs, save)) return false;
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
    yocto_shape* shp, const yocto_scene* scn, file_stream& fs, bool save) {
    if (!serialize_bin_value(shp->name, fs, save)) return false;
    if (!serialize_bin_value(shp->filename, fs, save)) return false;
    if (!serialize_bin_value(shp->points, fs, save)) return false;
    if (!serialize_bin_value(shp->lines, fs, save)) return false;
    if (!serialize_bin_value(shp->triangles, fs, save)) return false;
    if (!serialize_bin_value(shp->positions, fs, save)) return false;
    if (!serialize_bin_value(shp->normals, fs, save)) return false;
    if (!serialize_bin_value(shp->texturecoords, fs, save)) return false;
    if (!serialize_bin_value(shp->colors, fs, save)) return false;
    if (!serialize_bin_value(shp->radius, fs, save)) return false;
    if (!serialize_bin_value(shp->tangent_spaces, fs, save)) return false;
    return true;
}

bool serialize_bin_object(yocto_surface* sbd, file_stream& fs, bool save) {
    if (!serialize_bin_value(sbd->name, fs, save)) return false;
    if (!serialize_bin_value(sbd->filename, fs, save)) return false;
    if (!serialize_bin_value(sbd->subdivision_level, fs, save)) return false;
    if (!serialize_bin_value(sbd->catmull_clark, fs, save)) return false;
    if (!serialize_bin_value(sbd->compute_normals, fs, save)) return false;
    if (!serialize_bin_value(sbd->positions_quads, fs, save)) return false;
    if (!serialize_bin_value(sbd->texturecoords_quads, fs, save)) return false;
    if (!serialize_bin_value(sbd->colors_quads, fs, save)) return false;
    if (!serialize_bin_value(sbd->positions_creases, fs, save)) return false;
    if (!serialize_bin_value(sbd->texturecoords_quads, fs, save)) return false;
    if (!serialize_bin_value(sbd->positions, fs, save)) return false;
    if (!serialize_bin_value(sbd->texturecoords, fs, save)) return false;
    if (!serialize_bin_value(sbd->colors, fs, save)) return false;
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

bool serialize_bin_object(yocto_environment* env, const yocto_scene* scn,
    file_stream& fs, bool save) {
    if (!serialize_bin_value(env->name, fs, save)) return false;
    if (!serialize_bin_value(env->frame, fs, save)) return false;
    if (!serialize_bin_value(env->emission, fs, save)) return false;
    if (!serialize_bin_handle(env->emission_texture, scn->textures, fs, save))
        return false;
    return true;
}

bool serialize_bin_object(
    yocto_material* mat, const yocto_scene* scn, file_stream& fs, bool save) {
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
    if (!serialize_bin_handle(mat->emission_texture, scn->textures, fs, save))
        return false;
    if (!serialize_bin_handle(mat->diffuse_texture, scn->textures, fs, save))
        return false;
    if (!serialize_bin_handle(mat->specular_texture, scn->textures, fs, save))
        return false;
    if (!serialize_bin_handle(mat->transmission_texture, scn->textures, fs, save))
        return false;
    if (!serialize_bin_handle(mat->roughness_texture, scn->textures, fs, save))
        return false;
    if (!serialize_bin_handle(mat->opacity_texture, scn->textures, fs, save))
        return false;
    if (!serialize_bin_handle(mat->occlusion_texture, scn->textures, fs, save))
        return false;
    if (!serialize_bin_handle(mat->bump_texture, scn->textures, fs, save))
        return false;
    if (!serialize_bin_handle(mat->displacement_texture, scn->textures, fs, save))
        return false;
    if (!serialize_bin_handle(mat->normal_texture, scn->textures, fs, save))
        return false;
    if (!serialize_bin_value(mat->volume_emission, fs, save)) return false;
    if (!serialize_bin_value(mat->volume_albedo, fs, save)) return false;
    if (!serialize_bin_value(mat->volume_density, fs, save)) return false;
    if (!serialize_bin_value(mat->volume_phaseg, fs, save)) return false;
    if (!serialize_bin_handle(
            mat->volume_density_texture, scn->voltextures, fs, save))
        return false;
    return true;
};

bool serialize_bin_object(
    yocto_instance* ist, const yocto_scene* scn, file_stream& fs, bool save) {
    if (!serialize_bin_value(ist->name, fs, save)) return false;
    if (!serialize_bin_value(ist->frame, fs, save)) return false;
    if (!serialize_bin_handle(ist->shape, scn->shapes, fs, save)) return false;
    if (!serialize_bin_handle(ist->material, scn->materials, fs, save))
        return false;
    if (!serialize_bin_handle(ist->surface, scn->surfaces, fs, save))
        return false;
    return true;
};

bool serialize_scene(yocto_scene* scn, file_stream& fs, bool save) {
    if (!serialize_bin_value(scn->name, fs, save)) return false;
    if (!serialize_bin_object(scn->cameras, fs, save)) return false;
    if (!serialize_bin_object(scn->shapes, scn, fs, save)) return false;
    if (!serialize_bin_object(scn->surfaces, fs, save)) return false;
    if (!serialize_bin_object(scn->textures, fs, save)) return false;
    if (!serialize_bin_object(scn->voltextures, fs, save)) return false;
    if (!serialize_bin_object(scn->materials, scn, fs, save)) return false;
    if (!serialize_bin_object(scn->instances, scn, fs, save)) return false;
    if (!serialize_bin_object(scn->environments, scn, fs, save)) return false;
    return true;
}

// Load/save a binary dump useful for very fast scene IO.
yocto_scene* load_ybin_scene(
    const string& filename, bool load_textures, bool skip_missing) {
    auto fs = open(filename, "rb");
    if (!fs) return nullptr;
    auto scn = make_unique<yocto_scene>();
    if (!serialize_scene(scn.get(), fs, false)) return nullptr;
    return scn.release();
}

// Load/save a binary dump useful for very fast scene IO.
bool save_ybin_scene(const string& filename, const yocto_scene* scn,
    bool save_textures, bool skip_missing) {
    auto fs = open(filename, "wb");
    if (!fs) return false;
    if (!serialize_scene((yocto_scene*)scn, fs, true)) return false;
    return true;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE IO
// -----------------------------------------------------------------------------
namespace ygl {

// Reset mesh data
void reset_mesh_data(vector<int>& points, vector<vec2i>& lines,
    vector<vec3i>& triangles, vector<vec3f>& pos, vector<vec3f>& norm,
    vector<vec2f>& texcoord, vector<vec4f>& color, vector<float>& radius) {
    points    = {};
    lines     = {};
    triangles = {};
    pos       = {};
    norm      = {};
    texcoord  = {};
    color     = {};
    radius    = {};
}

// Load ply mesh
bool load_mesh(const string& filename, vector<int>& points, vector<vec2i>& lines,
    vector<vec3i>& triangles, vector<vec3f>& pos, vector<vec3f>& norm,
    vector<vec2f>& texcoord, vector<vec4f>& color, vector<float>& radius) {
    auto ext = get_extension(filename);
    if (ext == "ply" || ext == "PLY") {
        return load_ply_mesh(filename, points, lines, triangles, pos, norm,
            texcoord, color, radius);
    } else if (ext == "obj" || ext == "OBJ") {
        return load_obj_mesh(
            filename, points, lines, triangles, pos, norm, texcoord);
    } else {
        reset_mesh_data(
            points, lines, triangles, pos, norm, texcoord, color, radius);
        return false;
    }
}

// Save ply mesh
bool save_mesh(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec3f>& pos, const vector<vec3f>& norm,
    const vector<vec2f>& texcoord, const vector<vec4f>& color,
    const vector<float>& radius, bool ascii) {
    auto ext = get_extension(filename);
    if (ext == "ply" || ext == "PLY") {
        return save_ply_mesh(filename, points, lines, triangles, pos, norm,
            texcoord, color, radius, ascii);
    } else if (ext == "obj" || ext == "OBJ") {
        return save_obj_mesh(
            filename, points, lines, triangles, pos, norm, texcoord);
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
                    throw runtime_error("unsupported ply list type");
                if (elem_type != "int")
                    throw runtime_error("unsupported ply list type");
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

// Load ply mesh
bool load_ply_mesh(const string& filename, vector<int>& points,
    vector<vec2i>& lines, vector<vec3i>& triangles, vector<vec3f>& pos,
    vector<vec3f>& norm, vector<vec2f>& texcoord, vector<vec4f>& color,
    vector<float>& radius) {
    // clear
    reset_mesh_data(points, lines, triangles, pos, norm, texcoord, color, radius);

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
                    for (auto i = 2; i < (int)prop.scalars[fid]; i++)
                        triangles.push_back({list[0], list[i - 1], list[i]});
                }
            }
        }
    }

    // done
    return true;
}

// Save ply mesh
bool save_ply_mesh(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec3f>& pos, const vector<vec3f>& norm,
    const vector<vec2f>& texcoord, const vector<vec4f>& color,
    const vector<float>& radius, bool ascii) {
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
    if (!triangles.empty()) {
        print(fs, "element face {}\n", (int)triangles.size());
        print(fs, "property list uchar int vertex_indices\n");
    }
    if (!lines.empty()) {
        print(fs, "element line {}\n", (int)lines.size());
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
        for (auto i = 0; i < triangles.size(); i++)
            print(fs, "3 {}\n", triangles[i]);
        for (auto i = 0; i < lines.size(); i++) print(fs, "2 {}\n", lines[i]);
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
        for (auto i = 0; i < triangles.size(); i++) {
            auto n = (byte)3;
            write_value(fs, n);
            write_value(fs, triangles[i]);
        }
        for (auto i = 0; i < lines.size(); i++) {
            auto n = (byte)2;
            write_value(fs, n);
            write_value(fs, lines[i]);
        }
    }

    // done
    return true;
}

// Load ply mesh
bool load_obj_mesh(const string& filename, vector<int>& points,
    vector<vec2i>& lines, vector<vec3i>& triangles, vector<vec3f>& pos,
    vector<vec3f>& norm, vector<vec2f>& texcoord, bool flip_texcoord) {
    // clear
    auto color  = vector<vec4f>{};
    auto radius = vector<float>{};
    reset_mesh_data(points, lines, triangles, pos, norm, texcoord, color, radius);

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
            if (vert.pos) pos.push_back(opos.at(vert.pos - 1));
            if (vert.texcoord)
                texcoord.push_back(otexcoord.at(vert.texcoord - 1));
            if (vert.norm) norm.push_back(onorm.at(vert.norm - 1));
        }
    };

    auto cb     = obj_callbacks();
    cb.vert     = [&](vec3f v) { opos.push_back(v); };
    cb.norm     = [&](vec3f v) { onorm.push_back(v); };
    cb.texcoord = [&](vec2f v) { otexcoord.push_back(v); };
    cb.face     = [&](const vector<obj_vertex>& verts) {
        add_verts(verts);
        for (auto i = 2; i < verts.size(); i++)
            triangles.push_back({vert_map.at(verts[0]),
                vert_map.at(verts[i - 1]), vert_map.at(verts[i])});
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

    // load obj
    return load_obj(filename, cb, true, true, flip_texcoord);
}

// Load ply mesh
bool save_obj_mesh(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec3f>& pos, const vector<vec3f>& norm,
    const vector<vec2f>& texcoord, bool flip_texcoord) {
    auto fs = open(filename, "wt");
    if (!fs) return false;

    for (auto& p : pos) print(fs, "v {}\n", p);
    for (auto& n : norm) print(fs, "vn {}\n", n);
    for (auto& t : texcoord)
        print(fs, "vt {}\n", vec2f{t.x, (flip_texcoord) ? 1 - t.y : t.y});
    auto mask = obj_vertex{1, texcoord.empty() ? 0 : 1, norm.empty() ? 0 : 1};
    auto vert = [mask](int i) {
        return obj_vertex{
            (i + 1) * mask.pos, (i + 1) * mask.texcoord, (i + 1) * mask.norm};
    };
    for (auto& t : triangles)
        print(fs, "f {} {} {}\n", to_string(vert(t.x)).c_str(),
            to_string(vert(t.y)).c_str(), to_string(vert(t.z)).c_str());
    for (auto& l : lines)
        print(fs, "l {} {}\n", to_string(vert(l.x)).c_str(),
            to_string(vert(l.y)).c_str());
    for (auto& p : points) print(fs, "p {}\n", to_string(vert(p)).c_str());

    return true;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE IO
// -----------------------------------------------------------------------------
namespace ygl {

// Reset mesh data
void reset_fvmesh_data(vector<vec4i>& quads_pos, vector<vec3f>& pos,
    vector<vec4i>& quads_norm, vector<vec3f>& norm, vector<vec4i>& quads_texcoord,
    vector<vec2f>& texcoord, vector<vec4i>& quads_color, vector<vec4f>& color) {
    quads_pos      = {};
    pos            = {};
    quads_norm     = {};
    norm           = {};
    quads_texcoord = {};
    texcoord       = {};
    quads_color    = {};
    color          = {};
}

// Load mesh
bool load_fvmesh(const string& filename, vector<vec4i>& quads_pos,
    vector<vec3f>& pos, vector<vec4i>& quads_norm, vector<vec3f>& norm,
    vector<vec4i>& quads_texcoord, vector<vec2f>& texcoord,
    vector<vec4i>& quads_color, vector<vec4f>& color) {
    auto ext = get_extension(filename);
    if (ext == "obj" || ext == "OBJ") {
        return load_obj_fvmesh(filename, quads_pos, pos, quads_norm, norm,
            quads_texcoord, texcoord);
    } else {
        reset_fvmesh_data(quads_pos, pos, quads_norm, norm, quads_texcoord,
            texcoord, quads_color, color);
        log_io_error("unsupported mesh format {}", ext);
        return false;
    }
}

// Save mesh
bool save_fvmesh(const string& filename, const vector<vec4i>& quads_pos,
    const vector<vec3f>& pos, const vector<vec4i>& quads_norm,
    const vector<vec3f>& norm, const vector<vec4i>& quads_texcoord,
    const vector<vec2f>& texcoord, const vector<vec4i>& quads_color,
    const vector<vec4f>& color, bool ascii) {
    auto ext = get_extension(filename);
    if (ext == "obj" || ext == "OBJ") {
        return save_obj_fvmesh(filename, quads_pos, pos, quads_norm, norm,
            quads_texcoord, texcoord);
    } else {
        log_io_error("unsupported mesh format {}", ext);
        return false;
    }
}

// Load obj mesh
bool load_obj_fvmesh(const string& filename, vector<vec4i>& quads_pos,
    vector<vec3f>& pos, vector<vec4i>& quads_norm, vector<vec3f>& norm,
    vector<vec4i>& quads_texcoord, vector<vec2f>& texcoord, bool flip_texcoord) {
    // clear
    vector<vec4i> quads_color;
    vector<vec4f> color;
    reset_fvmesh_data(quads_pos, pos, quads_norm, norm, quads_texcoord,
        texcoord, quads_color, color);

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
            if (!vert.pos) continue;
            auto pos_it = pos_map.find(vert.pos);
            if (pos_it != pos_map.end()) continue;
            auto nverts = (int)pos.size();
            pos_map.insert(pos_it, {vert.pos, nverts});
            pos.push_back(opos.at(vert.pos - 1));
        }
        for (auto& vert : verts) {
            if (!vert.texcoord) continue;
            auto texcoord_it = texcoord_map.find(vert.texcoord);
            if (texcoord_it != texcoord_map.end()) continue;
            auto nverts = (int)texcoord.size();
            texcoord_map.insert(texcoord_it, {vert.texcoord, nverts});
            texcoord.push_back(otexcoord.at(vert.texcoord - 1));
        }
        for (auto& vert : verts) {
            if (!vert.norm) continue;
            auto norm_it = norm_map.find(vert.norm);
            if (norm_it != norm_map.end()) continue;
            auto nverts = (int)norm.size();
            norm_map.insert(norm_it, {vert.norm, nverts});
            norm.push_back(onorm.at(vert.norm - 1));
        }
    };

    auto cb     = obj_callbacks();
    cb.vert     = [&](vec3f v) { opos.push_back(v); };
    cb.norm     = [&](vec3f v) { onorm.push_back(v); };
    cb.texcoord = [&](vec2f v) { otexcoord.push_back(v); };
    cb.face     = [&](const vector<obj_vertex>& verts) {
        add_verts(verts);
        if (verts.size() == 4) {
            if (verts[0].pos) {
                quads_pos.push_back(
                    {pos_map.at(verts[0].pos), pos_map.at(verts[1].pos),
                        pos_map.at(verts[2].pos), pos_map.at(verts[3].pos)});
            }
            if (verts[0].texcoord) {
                quads_texcoord.push_back({texcoord_map.at(verts[0].texcoord),
                    texcoord_map.at(verts[1].texcoord),
                    texcoord_map.at(verts[2].texcoord),
                    texcoord_map.at(verts[3].texcoord)});
            }
            if (verts[0].norm) {
                quads_norm.push_back({norm_map.at(verts[0].norm),
                    norm_map.at(verts[1].norm), norm_map.at(verts[2].norm),
                    norm_map.at(verts[3].norm)});
            }
        } else {
            if (verts[0].pos) {
                for (auto i = 2; i < verts.size(); i++)
                    quads_pos.push_back({pos_map.at(verts[0].pos),
                        pos_map.at(verts[1].pos), pos_map.at(verts[i].pos),
                        pos_map.at(verts[i].pos)});
            }
            if (verts[0].texcoord) {
                for (auto i = 2; i < verts.size(); i++)
                    quads_texcoord.push_back({texcoord_map.at(verts[0].texcoord),
                        texcoord_map.at(verts[1].texcoord),
                        texcoord_map.at(verts[i].texcoord),
                        texcoord_map.at(verts[i].texcoord)});
            }
            if (verts[0].norm) {
                for (auto i = 2; i < verts.size(); i++)
                    quads_norm.push_back({norm_map.at(verts[0].norm),
                        norm_map.at(verts[1].norm), norm_map.at(verts[i].norm),
                        norm_map.at(verts[i].norm)});
            }
        }
    };

    // load obj
    return load_obj(filename, cb, true, true, flip_texcoord);
}

// Load ply mesh
bool save_obj_fvmesh(const string& filename, const vector<vec4i>& quads_pos,
    const vector<vec3f>& pos, const vector<vec4i>& quads_norm,
    const vector<vec3f>& norm, const vector<vec4i>& quads_texcoord,
    const vector<vec2f>& texcoord, bool flip_texcoord) {
    auto fs = open(filename, "wt");
    if (!fs) return false;

    for (auto& p : pos) print(fs, "v {}\n", p);
    for (auto& n : norm) print(fs, "vn {}\n", n);
    for (auto& t : texcoord)
        print(fs, "vt {}\n", vec2f{t.x, (flip_texcoord) ? 1 - t.y : t.y});
    auto mask = obj_vertex{1, texcoord.empty() ? 0 : 1, norm.empty() ? 0 : 1};
    auto vert = [mask](int pif, int ti, int ni) {
        return obj_vertex{(pif + 1) * mask.pos, (ti + 1) * mask.texcoord,
            (ni + 1) * mask.norm};
    };
    for (auto i = 0; i < quads_pos.size(); i++) {
        auto qp = quads_pos.at(i);
        auto qt = !quads_texcoord.empty() ? quads_texcoord.at(i) :
                                            vec4i{-1, -1, -1, -1};
        auto qn = !quads_norm.empty() ? quads_norm.at(i) : vec4i{-1, -1, -1, -1};
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
