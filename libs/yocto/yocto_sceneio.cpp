//
// Implementation for Yocto/Scene Input and Output functions.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
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

#include "yocto_sceneio.h"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <unordered_map>

#define JSON_USE_IMPLICIT_CONVERSIONS 0
#include "ext/json.hpp"
#include "ext/stb_image.h"
#include "ext/stb_image_resize.h"
#include "ext/stb_image_write.h"
#include "ext/tinyexr.h"
#include "yocto_color.h"
#include "yocto_geometry.h"
#include "yocto_image.h"
#include "yocto_modelio.h"
#include "yocto_parallel.h"
#include "yocto_shading.h"
#include "yocto_shape.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::unique_ptr;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// PATH UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Make a path from a utf8 string
static std::filesystem::path make_path(const string& filename) {
  return std::filesystem::u8path(filename);
}

// Normalize path
string normalize_path(const string& filename) {
  return make_path(filename).generic_u8string();
}

// Get directory name (not including /)
string path_dirname(const string& filename) {
  return make_path(filename).parent_path().generic_u8string();
}

// Get extension (including .)
string path_extension(const string& filename) {
  return make_path(filename).extension().u8string();
}

// Get filename without directory.
string path_filename(const string& filename) {
  return make_path(filename).filename().u8string();
}

// Get filename without directory and extension.
string path_basename(const string& filename) {
  return make_path(filename).stem().u8string();
}

// Joins paths
string path_join(const string& patha, const string& pathb) {
  return (make_path(patha) / make_path(pathb)).generic_u8string();
}
string path_join(
    const string& patha, const string& pathb, const string& pathc) {
  return (make_path(patha) / make_path(pathb) / make_path(pathc))
      .generic_u8string();
}

// Replaces extensions
string replace_extension(const string& filename, const string& ext) {
  return make_path(filename).replace_extension(ext).u8string();
}

// Check if a file can be opened for reading.
bool path_exists(const string& filename) { return exists(make_path(filename)); }

// Check if a file is a directory
bool path_isdir(const string& filename) {
  return is_directory(make_path(filename));
}

// Check if a file is a file
bool path_isfile(const string& filename) {
  return is_regular_file(make_path(filename));
}

// Get the current directory
string path_current() { return std::filesystem::current_path().u8string(); }

// Create a directory and all missing parent directories if needed
bool make_directory(const string& dirname, string& error) {
  if (path_exists(dirname)) return true;
  try {
    create_directories(make_path(dirname));
    return true;
  } catch (...) {
    error = dirname + ": cannot create directory";
    return false;
  }
}

// List the contents of a directory
bool list_directory(
    const string& dirname, vector<string>& entries, string& error) {
  entries.clear();
  try {
    for (auto entry : std::filesystem::directory_iterator(make_path(dirname))) {
      entries.push_back(entry.path().generic_u8string());
    }
    return true;
  } catch (...) {
    error = dirname + ": cannot list directory";
    return false;
  }
}

// Create a directory and all missing parent directories if needed
io_status make_directory(const string& dirname) {
  auto error = string{};
  if (!make_directory(dirname, error)) return io_status{error};
  return io_status{};
}

// List the contents of a directory
pair<io_status, vector<string>> list_directory(const string& dirname) {
  auto error   = string{};
  auto entries = vector<string>{};
  if (!list_directory(dirname, entries, error)) return {io_status{error}, {}};
  return {io_status{}, std::move(entries)};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// FILE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Opens a file with a utf8 file name
static FILE* fopen_utf8(const char* filename, const char* mode) {
#ifdef _WIN32
  auto path8    = std::filesystem::u8path(filename);
  auto str_mode = string{mode};
  auto wmode    = std::wstring(str_mode.begin(), str_mode.end());
  return _wfopen(path8.c_str(), wmode.c_str());
#else
  return fopen(filename, mode);
#endif
}

// Opens a file with utf8 filename
FILE* fopen_utf8(const string& filename, const string& mode) {
#ifdef _WIN32
  auto path8 = std::filesystem::u8path(filename);
  auto wmode = std::wstring(mode.begin(), mode.end());
  return _wfopen(path8.c_str(), wmode.c_str());
#else
  return fopen(filename.c_str(), mode.c_str());
#endif
}

// Load a text file
pair<io_status, string> load_text(const string& filename) {
  auto error = string{};
  auto str   = string{};
  if (!load_text(filename, str, error)) return {io_status{error}, {}};
  return {io_status{}, std::move(str)};
}
io_status load_text(const string& filename, string& text) {
  auto error = string{};
  if (!load_text(filename, text, error)) return io_status{error};
  return io_status{};
}

// Save a text file
io_status save_text(const string& filename, const string& text) {
  auto error = string{};
  if (!save_text(filename, text, error)) return io_status{error};
  return io_status{};
}

// Load a binary file
pair<io_status, vector<byte>> load_binary(const string& filename) {
  auto error = string{};
  auto data  = vector<byte>{};
  if (!load_binary(filename, data, error)) return {io_status{error}, {}};
  return {io_status{}, std::move(data)};
}
io_status load_binary(const string& filename, vector<byte>& data) {
  auto error = string{};
  if (!load_binary(filename, data, error)) return io_status{error};
  return io_status{};
}

// Save a binary file
io_status save_binary(const string& filename, const vector<byte>& data) {
  auto error = string{};
  if (!save_binary(filename, data, error)) return io_status{error};
  return io_status{};
}

// Load a text file
bool load_text(const string& filename, string& str, string& error) {
  // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
  auto fs = fopen_utf8(filename.c_str(), "rb");
  if (!fs) {
    error = filename + ": file not found";
    return false;
  }
  fseek(fs, 0, SEEK_END);
  auto length = ftell(fs);
  fseek(fs, 0, SEEK_SET);
  str.resize(length);
  if (fread(str.data(), 1, length, fs) != length) {
    fclose(fs);
    error = filename + ": read error";
    return false;
  }
  fclose(fs);
  return true;
}

// Save a text file
bool save_text(const string& filename, const string& str, string& error) {
  auto fs = fopen_utf8(filename.c_str(), "wt");
  if (!fs) {
    error = filename + ": file not found";
    return false;
  }
  if (fprintf(fs, "%s", str.c_str()) < 0) {
    fclose(fs);
    error = filename + ": write error";
    return false;
  }
  fclose(fs);
  return true;
}

// Load a binary file
bool load_binary(const string& filename, vector<byte>& data, string& error) {
  // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
  auto fs = fopen_utf8(filename.c_str(), "rb");
  if (!fs) {
    error = filename + ": file not found";
    return false;
  }
  fseek(fs, 0, SEEK_END);
  auto length = ftell(fs);
  fseek(fs, 0, SEEK_SET);
  data.resize(length);
  if (fread(data.data(), 1, length, fs) != length) {
    fclose(fs);
    error = filename + ": read error";
    return false;
  }
  fclose(fs);
  return true;
}

// Save a binary file
bool save_binary(
    const string& filename, const vector<byte>& data, string& error) {
  auto fs = fopen_utf8(filename.c_str(), "wb");
  if (!fs) {
    error = filename + ": file not found";
    return false;
  }
  if (fwrite(data.data(), 1, data.size(), fs) != data.size()) {
    fclose(fs);
    error = filename + ": write error";
    return false;
  }
  fclose(fs);
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON SUPPORT
// -----------------------------------------------------------------------------
namespace yocto {

// nlohmann json
using nlohmann::ordered_json;

// Load/save json
static bool load_json(
    const string& filename, ordered_json& json, string& error) {
  auto text = string{};
  if (!load_text(filename, text, error)) return false;
  try {
    json = ordered_json::parse(text);
    return true;
  } catch (...) {
    error = filename + ": json parse error";
    return false;
  }
}
static bool save_json(
    const string& filename, const ordered_json& json, string& error) {
  return save_text(filename, json.dump(2), error);
}

// conversions
inline void to_json(ordered_json& json, const vec2f& value) {
  nlohmann::to_json(json, (const array<float, 2>&)value);
}
inline void to_json(ordered_json& json, const vec3f& value) {
  nlohmann::to_json(json, (const array<float, 3>&)value);
}
inline void to_json(ordered_json& json, const vec4f& value) {
  nlohmann::to_json(json, (const array<float, 4>&)value);
}
inline void to_json(ordered_json& json, const mat3f& value) {
  nlohmann::to_json(json, (const array<float, 9>&)value);
}
inline void to_json(ordered_json& json, const mat4f& value) {
  nlohmann::to_json(json, (const array<float, 16>&)value);
}
inline void to_json(ordered_json& json, const frame3f& value) {
  nlohmann::to_json(json, (const array<float, 12>&)value);
}
inline void from_json(const ordered_json& json, const vec2f& value) {
  nlohmann::from_json(json, (array<float, 2>&)value);
}
inline void from_json(const ordered_json& json, const vec3f& value) {
  nlohmann::from_json(json, (array<float, 3>&)value);
}
inline void from_json(const ordered_json& json, const vec4f& value) {
  nlohmann::from_json(json, (array<float, 4>&)value);
}
inline void from_json(const ordered_json& json, const mat3f& value) {
  nlohmann::from_json(json, (array<float, 12>&)value);
}
inline void from_json(const ordered_json& json, const mat4f& value) {
  nlohmann::from_json(json, (array<float, 16>&)value);
}
inline void from_json(const ordered_json& json, const frame3f& value) {
  nlohmann::from_json(json, (array<float, 12>&)value);
}

// setup json value type
using json_value = nlohmann::ordered_json;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Check if an image is HDR based on filename.
bool is_hdr_filename(const string& filename) {
  auto ext = path_extension(filename);
  return ext == ".hdr" || ext == ".exr" || ext == ".pfm";
}

bool is_ldr_filename(const string& filename) {
  auto ext = path_extension(filename);
  return ext == ".png" || ext == ".jpg" || ext == ".jpeg" || ext == ".bmp" ||
         ext == ".tga";
}

// Loads/saves an image. Chooses hdr or ldr based on file name.
bool load_image(const string& filename, image_data& image, string& error) {
  auto read_error = [&]() {
    error = filename + ": read error";
    return false;
  };

  // conversion helpers
  auto from_linear = [](const float* pixels, int width, int height) {
    return vector<vec4f>{
        (vec4f*)pixels, (vec4f*)pixels + (size_t)width * (size_t)height};
  };
  auto from_srgb = [](const byte* pixels, int width, int height) {
    auto pixelsf = vector<vec4f>((size_t)width * (size_t)height);
    for (auto idx = (size_t)0; idx < pixelsf.size(); idx++) {
      pixelsf[idx] = byte_to_float(((vec4b*)pixels)[idx]);
    }
    return pixelsf;
  };

  auto ext = path_extension(filename);
  if (ext == ".exr" || ext == ".EXR") {
    auto buffer = vector<byte>{};
    if (!load_binary(filename, buffer, error)) return false;
    auto pixels = (float*)nullptr;
    if (LoadEXRFromMemory(&pixels, &image.width, &image.height, buffer.data(),
            buffer.size(), nullptr) != 0)
      return read_error();
    image.linear = true;
    image.pixels = from_linear(pixels, image.width, image.height);
    free(pixels);
    return true;
  } else if (ext == ".hdr" || ext == ".HDR") {
    auto buffer = vector<byte>{};
    if (!load_binary(filename, buffer, error)) return false;
    auto ncomp  = 0;
    auto pixels = stbi_loadf_from_memory(buffer.data(), (int)buffer.size(),
        &image.width, &image.height, &ncomp, 4);
    if (!pixels) return read_error();
    image.linear = true;
    image.pixels = from_linear(pixels, image.width, image.height);
    free(pixels);
    return true;
  } else if (ext == ".png" || ext == ".PNG") {
    auto buffer = vector<byte>{};
    if (!load_binary(filename, buffer, error)) return false;
    auto ncomp  = 0;
    auto pixels = stbi_load_from_memory(buffer.data(), (int)buffer.size(),
        &image.width, &image.height, &ncomp, 4);
    if (!pixels) return read_error();
    image.linear = false;
    image.pixels = from_srgb(pixels, image.width, image.height);
    free(pixels);
    return true;
  } else if (ext == ".jpg" || ext == ".JPG" || ext == ".jpeg" ||
             ext == ".JPEG") {
    auto buffer = vector<byte>{};
    if (!load_binary(filename, buffer, error)) return false;
    auto ncomp  = 0;
    auto pixels = stbi_load_from_memory(buffer.data(), (int)buffer.size(),
        &image.width, &image.height, &ncomp, 4);
    if (!pixels) return read_error();
    image.linear = false;
    image.pixels = from_srgb(pixels, image.width, image.height);
    free(pixels);
    return true;
  } else if (ext == ".tga" || ext == ".TGA") {
    auto buffer = vector<byte>{};
    if (!load_binary(filename, buffer, error)) return false;
    auto ncomp  = 0;
    auto pixels = stbi_load_from_memory(buffer.data(), (int)buffer.size(),
        &image.width, &image.height, &ncomp, 4);
    if (!pixels) return read_error();
    image.linear = false;
    image.pixels = from_srgb(pixels, image.width, image.height);
    free(pixels);
    return true;
  } else if (ext == ".bmp" || ext == ".BMP") {
    auto buffer = vector<byte>{};
    if (!load_binary(filename, buffer, error)) return false;
    auto ncomp  = 0;
    auto pixels = stbi_load_from_memory(buffer.data(), (int)buffer.size(),
        &image.width, &image.height, &ncomp, 4);
    if (!pixels) return read_error();
    image.linear = false;
    image.pixels = from_srgb(pixels, image.width, image.height);
    free(pixels);
    return true;
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    // create preset
    if (!make_image_preset(filename, image, error)) return false;
    return true;
  } else {
    error = filename + ": unknown format";
    return false;
  }
}

// Saves an hdr image.
bool save_image(
    const string& filename, const image_data& image, string& error) {
  auto write_error = [&]() {
    error = filename + ": write error";
    return false;
  };

  // conversion helpers
  auto to_linear = [](const image_data& image) {
    if (image.linear) return image.pixels;
    auto pixelsf = vector<vec4f>(image.pixels.size());
    srgb_to_rgb(pixelsf, image.pixels);
    return pixelsf;
  };
  auto to_srgb = [](const image_data& image) {
    auto pixelsb = vector<vec4b>(image.pixels.size());
    if (image.linear) {
      rgb_to_srgb(pixelsb, image.pixels);
    } else {
      float_to_byte(pixelsb, image.pixels);
    }
    return pixelsb;
  };

  // write data
  auto stbi_write_data = [](void* context, void* data, int size) {
    auto& buffer = *(vector<byte>*)context;
    buffer.insert(buffer.end(), (byte*)data, (byte*)data + size);
  };

  auto ext = path_extension(filename);
  if (ext == ".hdr" || ext == ".HDR") {
    auto buffer = vector<byte>{};
    if (!stbi_write_hdr_to_func(stbi_write_data, &buffer, (int)image.width,
            (int)image.height, 4, (const float*)to_linear(image).data()))
      return write_error();
    if (!save_binary(filename, buffer, error)) return false;
    return true;
  } else if (ext == ".exr" || ext == ".EXR") {
    auto data = (byte*)nullptr;
    auto size = (size_t)0;
    if (SaveEXRToMemory((const float*)to_linear(image).data(), (int)image.width,
            (int)image.height, 4, 1, &data, &size, nullptr) < 0)
      return write_error();
    auto buffer = vector<byte>{data, data + size};
    free(data);
    if (!save_binary(filename, buffer, error)) return false;
    return true;
  } else if (ext == ".png" || ext == ".PNG") {
    auto buffer = vector<byte>{};
    if (!stbi_write_png_to_func(stbi_write_data, &buffer, (int)image.width,
            (int)image.height, 4, (const byte*)to_srgb(image).data(),
            (int)image.width * 4))
      return write_error();
    if (!save_binary(filename, buffer, error)) return false;
    return true;
  } else if (ext == ".jpg" || ext == ".JPG" || ext == ".jpeg" ||
             ext == ".JPEG") {
    auto buffer = vector<byte>{};
    if (!stbi_write_jpg_to_func(stbi_write_data, &buffer, (int)image.width,
            (int)image.height, 4, (const byte*)to_srgb(image).data(), 75))
      return write_error();
    if (!save_binary(filename, buffer, error)) return false;
    return true;
  } else if (ext == ".tga" || ext == ".TGA") {
    auto buffer = vector<byte>{};
    if (!stbi_write_tga_to_func(stbi_write_data, &buffer, (int)image.width,
            (int)image.height, 4, (const byte*)to_srgb(image).data()))
      return write_error();
    if (!save_binary(filename, buffer, error)) return false;
    return true;
  } else if (ext == ".bmp" || ext == ".BMP") {
    auto buffer = vector<byte>{};
    if (!stbi_write_bmp_to_func(stbi_write_data, &buffer, (int)image.width,
            (int)image.height, 4, (const byte*)to_srgb(image).data()))
      return write_error();
    if (!save_binary(filename, buffer, error)) return false;
    return true;
  } else {
    error = filename + ": unknown format";
    return false;
  }
}

image_data make_image_preset(const string& type_) {
  auto type  = path_basename(type_);
  auto width = 1024, height = 1024;
  if (type.find("sky") != type.npos) width = 2048;
  if (type.find("images2") != type.npos) width = 2048;
  if (type == "grid") {
    return make_grid(width, height);
  } else if (type == "checker") {
    return make_checker(width, height);
  } else if (type == "bumps") {
    return make_bumps(width, height);
  } else if (type == "uvramp") {
    return make_uvramp(width, height);
  } else if (type == "gammaramp") {
    return make_gammaramp(width, height);
  } else if (type == "blackbodyramp") {
    return make_blackbodyramp(width, height);
  } else if (type == "uvgrid") {
    return make_uvgrid(width, height);
  } else if (type == "colormapramp") {
    return make_colormapramp(width, height);
  } else if (type == "sky") {
    return make_sunsky(width, height, pif / 4, 3.0f, false, 1.0f, 1.0f,
        vec3f{0.7f, 0.7f, 0.7f});
  } else if (type == "sunsky") {
    return make_sunsky(width, height, pif / 4, 3.0f, true, 1.0f, 1.0f,
        vec3f{0.7f, 0.7f, 0.7f});
  } else if (type == "noise") {
    return make_noisemap(width, height, 1);
  } else if (type == "fbm") {
    return make_fbmmap(width, height, 1);
  } else if (type == "ridge") {
    return make_ridgemap(width, height, 1);
  } else if (type == "turbulence") {
    return make_turbulencemap(width, height, 1);
  } else if (type == "bump-normal") {
    return make_bumps(width, height);
    // TODO(fabio): fix color space
    // img   = srgb_to_rgb(bump_to_normal(img, 0.05f));
  } else if (type == "images1") {
    auto sub_types  = vector<string>{"grid", "uvgrid", "checker", "gammaramp",
        "bumps", "bump-normal", "noise", "fbm", "blackbodyramp"};
    auto sub_images = vector<image_data>();
    for (auto& sub_type : sub_types)
      sub_images.push_back(make_image_preset(sub_type));
    auto montage_size = zero2i;
    for (auto& sub_image : sub_images) {
      montage_size.x += sub_image.width;
      montage_size.y = max(montage_size.y, sub_image.height);
    }
    auto image = make_image(
        montage_size.x, montage_size.y, sub_images[0].linear);
    auto pos = 0;
    for (auto& sub_image : sub_images) {
      set_region(image, sub_image, pos, 0);
      pos += sub_image.width;
    }
    return image;
  } else if (type == "images2") {
    auto sub_types  = vector<string>{"sky", "sunsky"};
    auto sub_images = vector<image_data>();
    for (auto& sub_type : sub_types)
      sub_images.push_back(make_image_preset(sub_type));
    auto montage_size = zero2i;
    for (auto& sub_image : sub_images) {
      montage_size.x += sub_image.width;
      montage_size.y = max(montage_size.y, sub_image.height);
    }
    auto image = make_image(
        montage_size.x, montage_size.y, sub_images[0].linear);
    auto pos = 0;
    for (auto& sub_image : sub_images) {
      set_region(image, sub_image, pos, 0);
      pos += sub_image.width;
    }
    return image;
  } else if (type == "test-floor") {
    return add_border(make_grid(width, height), 0.0025f);
  } else if (type == "test-grid") {
    return make_grid(width, height);
  } else if (type == "test-checker") {
    return make_checker(width, height);
  } else if (type == "test-bumps") {
    return make_bumps(width, height);
  } else if (type == "test-uvramp") {
    return make_uvramp(width, height);
  } else if (type == "test-gammaramp") {
    return make_gammaramp(width, height);
  } else if (type == "test-blackbodyramp") {
    return make_blackbodyramp(width, height);
  } else if (type == "test-colormapramp") {
    return make_colormapramp(width, height);
    // TODO(fabio): fix color space
    // img   = srgb_to_rgb(img);
  } else if (type == "test-uvgrid") {
    return make_uvgrid(width, height);
  } else if (type == "test-sky") {
    return make_sunsky(width, height, pif / 4, 3.0f, false, 1.0f, 1.0f,
        vec3f{0.7f, 0.7f, 0.7f});
  } else if (type == "test-sunsky") {
    return make_sunsky(width, height, pif / 4, 3.0f, true, 1.0f, 1.0f,
        vec3f{0.7f, 0.7f, 0.7f});
  } else if (type == "test-noise") {
    return make_noisemap(width, height);
  } else if (type == "test-fbm") {
    return make_noisemap(width, height);
  } else if (type == "test-bumps-normal") {
    return bump_to_normal(make_bumps(width, height), 0.05f);
  } else if (type == "test-bumps-displacement") {
    return make_bumps(width, height);
    // TODO(fabio): fix color space
    // img   = srgb_to_rgb(img);
  } else if (type == "test-fbm-displacement") {
    return make_fbmmap(width, height);
    // TODO(fabio): fix color space
    // img   = srgb_to_rgb(img);
  } else if (type == "test-checker-opacity") {
    return make_checker(width, height, 1, {1, 1, 1, 1}, {0, 0, 0, 0});
  } else if (type == "test-grid-opacity") {
    return make_grid(width, height, 1, {1, 1, 1, 1}, {0, 0, 0, 0});
  } else {
    return {};
  }
}

// Loads/saves an image. Chooses hdr or ldr based on file name.
image_data load_image(const string& filename, string& error) {
  auto image = image_data{};
  if (!load_image(filename, image, error)) return image_data{};
  return image;
}
pair<io_status, image_data> load_image(const string& filename) {
  auto error = string{};
  auto image = image_data{};
  if (!load_image(filename, image, error)) return {io_status{error}, {}};
  return {io_status{}, std::move(image)};
}
io_status load_image(const string& filename, image_data& image) {
  auto error = string{};
  if (!load_image(filename, image, error)) return io_status{error};
  return io_status{};
}
io_status save_image(const string& filename, const image_data& image) {
  auto error = string{};
  if (!save_image(filename, image, error)) return io_status{error};
  return io_status{};
}

bool make_image_preset(
    const string& filename, image_data& image, string& error) {
  image = make_image_preset(path_basename(filename));
  if (image.pixels.empty()) {
    error = "unknown preset";
    return false;
  }
  return true;
}

}  // namespace yocto

#if 0

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR VOLUME IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Volume load
static bool load_yvol(const string& filename, int& width, int& height,
    int& depth, int& components, vector<float>& voxels, string& error) {
  // error helpers
  auto open_error = [filename, &error]() {
    error = filename + ": file not found";
    return false;
  };
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error";
    return false;
  };
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };

  // Split a string
  auto split_string = [](const string& str) -> vector<string> {
    auto ret = vector<string>();
    if (str.empty()) return ret;
    auto lpos = (size_t)0;
    while (lpos != string::npos) {
      auto pos = str.find_first_of(" \t\n\r", lpos);
      if (pos != string::npos) {
        if (pos > lpos) ret.push_back(str.substr(lpos, pos - lpos));
        lpos = pos + 1;
      } else {
        if (lpos < str.size()) ret.push_back(str.substr(lpos));
        lpos = pos;
      }
    }
    return ret;
  };

  auto fs       = fopen_utf8(filename.c_str(), "rb");
  auto fs_guard = unique_ptr<FILE, int (*)(FILE*)>(fs, &fclose);
  if (!fs) return open_error();

  // buffer
  auto buffer = array<char, 4096>{};
  auto toks   = vector<string>();

  // read magic
  if (!fgets(buffer.data(), (int)buffer.size(), fs)) return parse_error();
  toks = split_string(buffer.data());
  if (toks[0] != "YVOL") return parse_error();

  // read width, height
  if (!fgets(buffer.data(), (int)buffer.size(), fs)) return parse_error();
  toks       = split_string(buffer.data());
  width      = atoi(toks[0].c_str());
  height     = atoi(toks[1].c_str());
  depth      = atoi(toks[2].c_str());
  components = atoi(toks[3].c_str());

  // read data
  auto nvoxels = (size_t)width * (size_t)height * (size_t)depth;
  auto nvalues = nvoxels * (size_t)components;
  voxels       = vector<float>(nvalues);
  if (!read_values(fs, voxels.data(), nvalues)) return read_error();

  // done
  return true;
}

// save pfm
static bool save_yvol(const string& filename, int width, int height, int depth,
    int components, const vector<float>& voxels, string& error) {
  // error helpers
  auto open_error = [filename, &error]() {
    error = filename + ": file not found";
    return false;
  };
  auto write_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };

  auto fs       = fopen_utf8(filename.c_str(), "wb");
  auto fs_guard = unique_ptr<FILE, int (*)(FILE*)>(fs, &fclose);
  if (!fs) return open_error();

  if (!write_text(fs, "YVOL\n")) return write_error();
  if (!write_text(fs, std::to_string(width) + " " + std::to_string(height) +
                          " " + std::to_string(depth) + " " +
                          std::to_string(components) + "\n"))
    return write_error();
  auto nvalues = (size_t)width * (size_t)height * (size_t)depth *
                 (size_t)components;
  if (!write_values(fs, voxels.data(), nvalues)) return write_error();
  return true;
}

// Loads volume data from binary format.
bool load_volume(const string& filename, volume<float>& vol, string& error) {
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };
  auto width = 0, height = 0, depth = 0, ncomp = 0;
  auto voxels = vector<float>{};
  if (!load_yvol(filename, width, height, depth, ncomp, voxels, error))
    return false;
  if (ncomp != 1) voxels = convert_components(voxels, ncomp, 1);
  vol = volume{{width, height, depth}, (const float*)voxels.data()};
  return true;
}

// Saves volume data in binary format.
bool save_volume(
    const string& filename, const volume<float>& vol, string& error) {
  return save_yvol(filename, vol.width(), vol.height(), vol.depth(), 1,
      {vol.data(), vol.data() + vol.count()}, error);
}

}  // namespace yocto

#endif

// -----------------------------------------------------------------------------
// SHAPE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load mesh
bool load_shape(const string& filename, shape_data& shape, string& error,
    bool flip_texcoord) {
  auto shape_error = [&]() {
    error = filename + ": empty shape";
    return false;
  };

  shape = {};

  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = ply_model{};
    if (!load_ply(filename, ply, error)) return false;
    get_positions(ply, shape.positions);
    get_normals(ply, shape.normals);
    get_texcoords(ply, shape.texcoords, flip_texcoord);
    get_colors(ply, shape.colors);
    get_radius(ply, shape.radius);
    get_faces(ply, shape.triangles, shape.quads);
    get_lines(ply, shape.lines);
    get_points(ply, shape.points);
    if (shape.points.empty() && shape.lines.empty() &&
        shape.triangles.empty() && shape.quads.empty())
      return shape_error();
    return true;
  } else if (ext == ".obj" || ext == ".OBJ") {
    auto obj = obj_shape{};
    if (!load_obj(filename, obj, error, false)) return false;
    auto materials = vector<int>{};
    get_positions(obj, shape.positions);
    get_normals(obj, shape.normals);
    get_texcoords(obj, shape.texcoords, flip_texcoord);
    get_faces(obj, shape.triangles, shape.quads, materials);
    get_lines(obj, shape.lines, materials);
    get_points(obj, shape.points, materials);
    if (shape.points.empty() && shape.lines.empty() &&
        shape.triangles.empty() && shape.quads.empty())
      return shape_error();
    return true;
  } else if (ext == ".stl" || ext == ".STL") {
    auto stl = stl_model{};
    if (!load_stl(filename, stl, error, true)) return false;
    if (stl.shapes.size() != 1) return shape_error();
    auto fnormals = vector<vec3f>{};
    if (!get_triangles(stl, 0, shape.triangles, shape.positions, fnormals))
      return shape_error();
    return true;
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    if (!make_shape_preset(filename, shape, error)) return false;
    return true;
  } else {
    error = filename + ": unknown format";
    return false;
  }
}

// Save ply mesh
bool save_shape(const string& filename, const shape_data& shape, string& error,
    bool flip_texcoord, bool ascii) {
  auto shape_error = [&]() {
    error = filename + ": empty shape";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = ply_model{};
    add_positions(ply, shape.positions);
    add_normals(ply, shape.normals);
    add_texcoords(ply, shape.texcoords, flip_texcoord);
    add_colors(ply, shape.colors);
    add_radius(ply, shape.radius);
    add_faces(ply, shape.triangles, shape.quads);
    add_lines(ply, shape.lines);
    add_points(ply, shape.points);
    if (!save_ply(filename, ply, error)) return false;
    return true;
  } else if (ext == ".obj" || ext == ".OBJ") {
    auto obj = obj_shape{};
    add_positions(obj, shape.positions);
    add_normals(obj, shape.normals);
    add_texcoords(obj, shape.texcoords, flip_texcoord);
    add_triangles(obj, shape.triangles, 0, !shape.normals.empty(),
        !shape.texcoords.empty());
    add_quads(
        obj, shape.quads, 0, !shape.normals.empty(), !shape.texcoords.empty());
    add_lines(
        obj, shape.lines, 0, !shape.normals.empty(), !shape.texcoords.empty());
    add_points(
        obj, shape.points, 0, !shape.normals.empty(), !shape.texcoords.empty());
    if (!save_obj(filename, obj, error)) return false;
    return true;
  } else if (ext == ".stl" || ext == ".STL") {
    auto stl = stl_model{};
    if (!shape.lines.empty()) return shape_error();
    if (!shape.points.empty()) return shape_error();
    if (!shape.triangles.empty()) {
      add_triangles(stl, shape.triangles, shape.positions, {});
    } else if (!shape.quads.empty()) {
      add_triangles(stl, quads_to_triangles(shape.quads), shape.positions, {});
    } else {
      return shape_error();
    }
    if (!save_stl(filename, stl, error)) return false;
    return true;
  } else if (ext == ".cpp" || ext == ".CPP") {
    auto to_cpp = [](const string& name, const string& vname,
                      const auto& values) -> string {
      using T = typename std::remove_const_t<
          std::remove_reference_t<decltype(values)>>::value_type;
      if (values.empty()) return ""s;
      auto str = "auto " + name + "_" + vname + " = ";
      if constexpr (std::is_same_v<int, T>) str += "vector<int>{\n";
      if constexpr (std::is_same_v<float, T>) str += "vector<float>{\n";
      if constexpr (std::is_same_v<vec2i, T>) str += "vector<vec2i>{\n";
      if constexpr (std::is_same_v<vec2f, T>) str += "vector<vec2f>{\n";
      if constexpr (std::is_same_v<vec3i, T>) str += "vector<vec3i>{\n";
      if constexpr (std::is_same_v<vec3f, T>) str += "vector<vec3f>{\n";
      if constexpr (std::is_same_v<vec4i, T>) str += "vector<vec4i>{\n";
      if constexpr (std::is_same_v<vec4f, T>) str += "vector<vec4f>{\n";
      for (auto& value : values) {
        if constexpr (std::is_same_v<int, T> || std::is_same_v<float, T>) {
          str += std::to_string(value) + ",\n";
        } else if constexpr (std::is_same_v<vec2i, T> ||
                             std::is_same_v<vec2f, T>) {
          str += "{" + std::to_string(value.x) + "," + std::to_string(value.y) +
                 "},\n";
        } else if constexpr (std::is_same_v<vec3i, T> ||
                             std::is_same_v<vec3f, T>) {
          str += "{" + std::to_string(value.x) + "," + std::to_string(value.y) +
                 "," + std::to_string(value.z) + "},\n";
        } else if constexpr (std::is_same_v<vec4i, T> ||
                             std::is_same_v<vec4f, T>) {
          str += "{" + std::to_string(value.x) + "," + std::to_string(value.y) +
                 "," + std::to_string(value.z) + "," + std::to_string(value.w) +
                 "},\n";
        } else {
          throw std::invalid_argument{"cannot print this"};
        }
      }
      str += "};\n\n";
      return str;
    };

    auto name = string{"shape"};
    auto str  = ""s;
    str += to_cpp(name, "positions", shape.positions);
    str += to_cpp(name, "normals", shape.normals);
    str += to_cpp(name, "texcoords", shape.texcoords);
    str += to_cpp(name, "colors", shape.colors);
    str += to_cpp(name, "radius", shape.radius);
    str += to_cpp(name, "points", shape.points);
    str += to_cpp(name, "lines", shape.lines);
    str += to_cpp(name, "triangles", shape.triangles);
    str += to_cpp(name, "quads", shape.quads);
    if (!save_text(filename, str, error)) return false;
    return true;
  } else {
    error = filename + ": unknown format";
    return false;
  }
}

// Load face-varying mesh
bool load_fvshape(const string& filename, fvshape_data& shape, string& error,
    bool flip_texcoord) {
  auto shape_error = [&]() {
    error = filename + ": empty shape";
    return false;
  };

  shape = {};

  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = ply_model{};
    if (!load_ply(filename, ply, error)) return false;
    get_positions(ply, shape.positions);
    get_normals(ply, shape.normals);
    get_texcoords(ply, shape.texcoords, flip_texcoord);
    get_quads(ply, shape.quadspos);
    if (!shape.normals.empty()) shape.quadsnorm = shape.quadspos;
    if (!shape.texcoords.empty()) shape.quadstexcoord = shape.quadspos;
    if (shape.quadspos.empty()) return shape_error();
    return true;
  } else if (ext == ".obj" || ext == ".OBJ") {
    auto obj = obj_shape{};
    if (!load_obj(filename, obj, error, true)) return false;
    auto materials = vector<int>{};
    get_positions(obj, shape.positions);
    get_normals(obj, shape.normals);
    get_texcoords(obj, shape.texcoords, flip_texcoord);
    get_fvquads(
        obj, shape.quadspos, shape.quadsnorm, shape.quadstexcoord, materials);
    if (shape.quadspos.empty()) return shape_error();
    return true;
  } else if (ext == ".stl" || ext == ".STL") {
    auto stl = stl_model{};
    if (!load_stl(filename, stl, error, true)) return false;
    if (stl.shapes.empty()) return shape_error();
    if (stl.shapes.size() > 1) return shape_error();
    auto fnormals  = vector<vec3f>{};
    auto triangles = vector<vec3i>{};
    if (!get_triangles(stl, 0, triangles, shape.positions, fnormals))
      return shape_error();
    shape.quadspos = triangles_to_quads(triangles);
    return true;
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    if (!make_fvshape_preset(filename, shape, error)) return false;
    return true;
  } else {
    error = filename + ": unknown format";
    return false;
  }
}

// Save ply mesh
bool save_fvshape(const string& filename, const fvshape_data& shape,
    string& error, bool flip_texcoord, bool ascii) {
  auto shape_error = [&]() {
    error = filename + ": empty shape";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply             = ply_model{};
    auto split_quads     = vector<vec4i>{};
    auto split_positions = vector<vec3f>{};
    auto split_normals   = vector<vec3f>{};
    auto split_texcoords = vector<vec2f>{};
    split_facevarying(split_quads, split_positions, split_normals,
        split_texcoords, shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
        shape.positions, shape.normals, shape.texcoords);
    add_positions(ply, split_positions);
    add_normals(ply, split_normals);
    add_texcoords(ply, split_texcoords, flip_texcoord);
    add_faces(ply, {}, split_quads);
    if (!save_ply(filename, ply, error)) return false;
    return true;
  } else if (ext == ".obj" || ext == ".OBJ") {
    auto obj = obj_shape{};
    add_positions(obj, shape.positions);
    add_normals(obj, shape.positions);
    add_texcoords(obj, shape.texcoords, flip_texcoord);
    add_fvquads(obj, shape.quadspos, shape.quadsnorm, shape.quadstexcoord, 0);
    if (!save_obj(filename, obj, error)) return false;
    return true;
  } else if (ext == ".stl" || ext == ".STL") {
    auto stl = stl_model{};
    if (!shape.quadspos.empty()) {
      auto split_quads     = vector<vec4i>{};
      auto split_positions = vector<vec3f>{};
      auto split_normals   = vector<vec3f>{};
      auto split_texcoords = vector<vec2f>{};
      split_facevarying(split_quads, split_positions, split_normals,
          split_texcoords, shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
          shape.positions, shape.normals, shape.texcoords);
      add_triangles(stl, quads_to_triangles(split_quads), split_positions, {});
    } else {
      return shape_error();
    }
    if (!save_stl(filename, stl, error)) return false;
    return true;
  } else if (ext == ".cpp" || ext == ".CPP") {
    auto to_cpp = [](const string& name, const string& vname,
                      const auto& values) -> string {
      using T = typename std::remove_const_t<
          std::remove_reference_t<decltype(values)>>::value_type;
      if (values.empty()) return ""s;
      auto str = "auto " + name + "_" + vname + " = ";
      if constexpr (std::is_same_v<int, T>) str += "vector<int>{\n";
      if constexpr (std::is_same_v<float, T>) str += "vector<float>{\n";
      if constexpr (std::is_same_v<vec2i, T>) str += "vector<vec2i>{\n";
      if constexpr (std::is_same_v<vec2f, T>) str += "vector<vec2f>{\n";
      if constexpr (std::is_same_v<vec3i, T>) str += "vector<vec3i>{\n";
      if constexpr (std::is_same_v<vec3f, T>) str += "vector<vec3f>{\n";
      if constexpr (std::is_same_v<vec4i, T>) str += "vector<vec4i>{\n";
      if constexpr (std::is_same_v<vec4f, T>) str += "vector<vec4f>{\n";
      for (auto& value : values) {
        if constexpr (std::is_same_v<int, T> || std::is_same_v<float, T>) {
          str += std::to_string(value) + ",\n";
        } else if constexpr (std::is_same_v<vec2i, T> ||
                             std::is_same_v<vec2f, T>) {
          str += "{" + std::to_string(value.x) + "," + std::to_string(value.y) +
                 "},\n";
        } else if constexpr (std::is_same_v<vec3i, T> ||
                             std::is_same_v<vec3f, T>) {
          str += "{" + std::to_string(value.x) + "," + std::to_string(value.y) +
                 "," + std::to_string(value.z) + "},\n";
        } else if constexpr (std::is_same_v<vec4i, T> ||
                             std::is_same_v<vec4f, T>) {
          str += "{" + std::to_string(value.x) + "," + std::to_string(value.y) +
                 "," + std::to_string(value.z) + "," + std::to_string(value.w) +
                 "},\n";
        } else {
          throw std::invalid_argument{"cannot print this"};
        }
      }
      str += "};\n\n";
      return str;
    };
    auto name = string{"shape"};
    auto str  = ""s;
    str += to_cpp(name, "positions", shape.positions);
    str += to_cpp(name, "normals", shape.normals);
    str += to_cpp(name, "texcoords", shape.texcoords);
    str += to_cpp(name, "quadspos", shape.quadspos);
    str += to_cpp(name, "quadsnorm", shape.quadsnorm);
    str += to_cpp(name, "quadstexcoord", shape.quadstexcoord);
    if (!save_text(filename, str, error)) return false;
    return true;
  } else {
    error = filename + ": unknown format";
    return false;
  }
}

// Shape presets used for testing.
shape_data make_shape_preset(const string& type) {
  if (type == "default-quad") {
    return make_rect();
  } else if (type == "default-quady") {
    return make_recty();
  } else if (type == "default-cube") {
    return make_box();
  } else if (type == "default-cube-rounded") {
    return make_rounded_box();
  } else if (type == "default-sphere") {
    return make_sphere();
  } else if (type == "default-matcube") {
    return make_rounded_box();
  } else if (type == "default-matsphere") {
    return make_uvspherey();
  } else if (type == "default-disk") {
    return make_disk();
  } else if (type == "default-disk-bulged") {
    return make_bulged_disk();
  } else if (type == "default-quad-bulged") {
    return make_bulged_rect();
  } else if (type == "default-uvsphere") {
    return make_uvsphere();
  } else if (type == "default-uvsphere-flipcap") {
    return make_capped_uvsphere();
  } else if (type == "default-uvspherey") {
    return make_uvspherey();
  } else if (type == "default-uvspherey-flipcap") {
    return make_capped_uvspherey();
  } else if (type == "default-uvdisk") {
    return make_uvdisk();
  } else if (type == "default-uvcylinder") {
    return make_uvcylinder();
  } else if (type == "default-uvcylinder-rounded") {
    return make_rounded_uvcylinder({32, 32, 32});
  } else if (type == "default-geosphere") {
    return make_geosphere();
  } else if (type == "default-floor") {
    return make_floor();
  } else if (type == "default-floor-bent") {
    return make_bent_floor();
  } else if (type == "default-matball") {
    return make_sphere();
  } else if (type == "default-hairball") {
    auto base = make_sphere(pow2(5), 0.8f);
    return make_hair(base, {4, 65536}, {0.2f, 0.2f}, {0.002f, 0.001f});
  } else if (type == "default-hairball-interior") {
    return make_sphere(pow2(5), 0.8f);
  } else if (type == "default-suzanne") {
    return make_monkey();
  } else if (type == "default-cube-facevarying") {
    return fvshape_to_shape(make_fvbox());
  } else if (type == "default-sphere-facevarying") {
    return fvshape_to_shape(make_fvsphere());
  } else if (type == "default-quady-displaced") {
    return make_recty({256, 256});
  } else if (type == "default-sphere-displaced") {
    return make_sphere(128);
  } else if (type == "test-cube") {
    auto shape = make_rounded_box(
        {32, 32, 32}, {0.075f, 0.075f, 0.075f}, {1, 1, 1}, 0.3f * 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-uvsphere") {
    auto shape = make_uvsphere({32, 32}, 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-uvsphere-flipcap") {
    auto shape = make_capped_uvsphere({32, 32}, 0.075f, {1, 1}, 0.3f * 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-uvspherey") {
    auto shape = make_uvspherey({32, 32}, 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-uvspherey-flipcap") {
    auto shape = make_capped_uvspherey({32, 32}, 0.075f, {1, 1}, 0.3f * 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-sphere") {
    auto shape = make_sphere(32, 0.075f, 1);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-matcube") {
    auto shape = make_rounded_box(
        {32, 32, 32}, {0.075f, 0.075f, 0.075f}, {1, 1, 1}, 0.3f * 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-matsphere") {
    auto shape = make_uvspherey({32, 32}, 0.075f, {2, 1});
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-sphere-displaced") {
    auto shape = make_sphere(128, 0.075f, 1);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-smallsphere") {
    auto shape = make_sphere(32, 0.015f, 1);
    for (auto& p : shape.positions) p += {0, 0.015f, 0};
    return shape;
  } else if (type == "test-disk") {
    auto shape = make_disk(32, 0.075f, 1);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-uvcylinder") {
    auto shape = make_rounded_uvcylinder(
        {32, 32, 32}, {0.075f, 0.075f}, {1, 1, 1}, 0.3f * 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-floor") {
    return make_floor({1, 1}, {2, 2}, {20, 20});
  } else if (type == "test-smallfloor") {
    return make_floor({1, 1}, {0.5f, 0.5f}, {1, 1});
  } else if (type == "test-quad") {
    return make_rect({1, 1}, {0.075f, 0.075f}, {1, 1});
  } else if (type == "test-quady") {
    return make_recty({1, 1}, {0.075f, 0.075f}, {1, 1});
  } else if (type == "test-quad-displaced") {
    return make_rect({256, 256}, {0.075f, 0.075f}, {1, 1});
  } else if (type == "test-quady-displaced") {
    return make_recty({256, 256}, {0.075f, 0.075f}, {1, 1});
  } else if (type == "test-matball") {
    auto shape = make_sphere(32, 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-geosphere") {
    auto shape = make_geosphere(0.075f, 3);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-geosphere-flat") {
    auto shape = make_geosphere(0.075f, 3);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    shape.normals = {};
    return shape;
  } else if (type == "test-geosphere-subdivided") {
    auto shape = make_geosphere(0.075f, 6);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-hairball1") {
    auto base = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : base.positions) p += {0, 0.075f, 0};
    return make_hair(base, {4, 65536}, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f}, {0.03f, 100});
  } else if (type == "test-hairball2") {
    auto base = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : base.positions) p += {0, 0.075f, 0};
    return make_hair(base, {4, 65536}, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f});
  } else if (type == "test-hairball3") {
    auto base = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : base.positions) p += {0, 0.075f, 0};
    return make_hair(base, {4, 65536}, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f}, {0, 0}, {0.5, 128});
  } else if (type == "test-hairball-interior") {
    auto shape = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-suzanne-subdiv") {
    auto shape = make_monkey(0.075f * 0.8f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-cube-subdiv") {
    auto fvshape    = make_fvcube(0.075f);
    auto shape      = shape_data{};
    shape.quads     = fvshape.quadspos;
    shape.positions = fvshape.positions;
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-arealight1") {
    return make_rect({1, 1}, {0.2f, 0.2f});
  } else if (type == "test-arealight2") {
    return make_rect({1, 1}, {0.2f, 0.2f});
  } else if (type == "test-largearealight1") {
    return make_rect({1, 1}, {0.4f, 0.4f});
  } else if (type == "test-largearealight2") {
    return make_rect({1, 1}, {0.4f, 0.4f});
  } else if (type == "test-pointlight1") {
    return make_point(0);
  } else if (type == "test-pointlight2") {
    return make_point(0);
  } else if (type == "test-point") {
    return make_points(1);
  } else if (type == "test-points") {
    return make_points(4096);
  } else if (type == "test-points-random") {
    auto shape = make_random_points(4096, {0.075f, 0.075f, 0.075f});
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-points-grid") {
    auto shape = make_points({256, 256}, {0.075f, 0.075f});
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    for (auto& r : shape.radius) r *= 0.075f;
    return shape;
  } else if (type == "test-lines-grid") {
    auto shape = make_lines({256, 256}, {0.075f, 0.075f});
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    for (auto& r : shape.radius) r *= 0.075f;
    return shape;
  } else if (type == "test-thickpoints-grid") {
    auto shape = make_points({16, 16}, {0.075f, 0.075f});
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    for (auto& r : shape.radius) r *= 0.075f * 10;
    return shape;
  } else if (type == "test-thicklines-grid") {
    auto shape = make_lines({16, 16}, {0.075f, 0.075f});
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    for (auto& r : shape.radius) r *= 0.075f * 10;
    return shape;
  } else if (type == "test-particles") {
    return make_points(4096);
  } else if (type == "test-cloth") {
    return make_rect({64, 64}, {0.2f, 0.2f});
  } else if (type == "test-clothy") {
    return make_recty({64, 64}, {0.2f, 0.2f});
  } else {
    return {};
  }
}

// Shape presets used for testing.
fvshape_data make_fvshape_preset(const string& type) {
  if (type == "default-quad") {
    return shape_to_fvshape(make_rect());
  } else if (type == "default-quady") {
    return shape_to_fvshape(make_recty());
  } else if (type == "default-cube") {
    return shape_to_fvshape(make_box());
  } else if (type == "default-cube-rounded") {
    return shape_to_fvshape(make_rounded_box());
  } else if (type == "default-sphere") {
    return shape_to_fvshape(make_sphere());
  } else if (type == "default-matcube") {
    return shape_to_fvshape(make_rounded_box());
  } else if (type == "default-matsphere") {
    return shape_to_fvshape(make_uvspherey());
  } else if (type == "default-disk") {
    return shape_to_fvshape(make_disk());
  } else if (type == "default-disk-bulged") {
    return shape_to_fvshape(make_bulged_disk());
  } else if (type == "default-quad-bulged") {
    return shape_to_fvshape(make_bulged_rect());
  } else if (type == "default-uvsphere") {
    return shape_to_fvshape(make_uvsphere());
  } else if (type == "default-uvsphere-flipcap") {
    return shape_to_fvshape(make_capped_uvsphere());
  } else if (type == "default-uvspherey") {
    return shape_to_fvshape(make_uvspherey());
  } else if (type == "default-uvspherey-flipcap") {
    return shape_to_fvshape(make_capped_uvspherey());
  } else if (type == "default-uvdisk") {
    return shape_to_fvshape(make_uvdisk());
  } else if (type == "default-uvcylinder") {
    return shape_to_fvshape(make_uvcylinder());
  } else if (type == "default-uvcylinder-rounded") {
    return shape_to_fvshape(make_rounded_uvcylinder({32, 32, 32}));
  } else if (type == "default-geosphere") {
    return shape_to_fvshape(make_geosphere());
  } else if (type == "default-floor") {
    return shape_to_fvshape(make_floor());
  } else if (type == "default-floor-bent") {
    return shape_to_fvshape(make_bent_floor());
  } else if (type == "default-matball") {
    return shape_to_fvshape(make_sphere());
  } else if (type == "default-hairball-interior") {
    return shape_to_fvshape(make_sphere(pow2(5), 0.8f));
  } else if (type == "default-suzanne") {
    return shape_to_fvshape(make_monkey());
  } else if (type == "default-cube-facevarying") {
    return make_fvbox();
  } else if (type == "default-sphere-facevarying") {
    return make_fvsphere();
  } else if (type == "default-quady-displaced") {
    return shape_to_fvshape(make_recty({256, 256}));
  } else if (type == "default-sphere-displaced") {
    return shape_to_fvshape(make_sphere(128));
  } else if (type == "test-cube") {
    auto shape = make_rounded_box(
        {32, 32, 32}, {0.075f, 0.075f, 0.075f}, {1, 1, 1}, 0.3f * 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-matsphere") {
    auto shape = make_uvspherey({32, 32}, 0.075f, {2, 1});
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-uvsphere") {
    auto shape = make_uvsphere({32, 32}, 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-uvsphere-flipcap") {
    auto shape = make_capped_uvsphere({32, 32}, 0.075f, {1, 1}, 0.3f * 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-uvspherey") {
    auto shape = make_uvspherey({32, 32}, 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-uvspherey-flipcap") {
    auto shape = make_capped_uvspherey({32, 32}, 0.075f, {1, 1}, 0.3f * 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-sphere") {
    auto shape = make_sphere(32, 0.075f, 1);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-sphere-displaced") {
    auto shape = make_sphere(128, 0.075f, 1);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-matcube") {
    auto shape = make_rounded_box(
        {32, 32, 32}, {0.075f, 0.075f, 0.075f}, {1, 1, 1}, 0.3f * 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-disk") {
    auto shape = make_disk(32, 0.075f, 1);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-uvcylinder") {
    auto shape = make_rounded_uvcylinder(
        {32, 32, 32}, {0.075f, 0.075f}, {1, 1, 1}, 0.3f * 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-floor") {
    return shape_to_fvshape(make_floor({1, 1}, {2, 2}, {20, 20}));
  } else if (type == "test-smallfloor") {
    return shape_to_fvshape(make_floor({1, 1}, {0.5f, 0.5f}, {1, 1}));
  } else if (type == "test-quad") {
    return shape_to_fvshape(make_rect({1, 1}, {0.075f, 0.075f}, {1, 1}));
  } else if (type == "test-quady") {
    return shape_to_fvshape(make_recty({1, 1}, {0.075f, 0.075f}, {1, 1}));
  } else if (type == "test-quad-displaced") {
    return shape_to_fvshape(make_rect({256, 256}, {0.075f, 0.075f}, {1, 1}));
  } else if (type == "test-quady-displaced") {
    return shape_to_fvshape(make_recty({256, 256}, {0.075f, 0.075f}, {1, 1}));
  } else if (type == "test-matball") {
    auto shape = make_sphere(32, 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-suzanne-subdiv") {
    auto shape = make_monkey(0.075f * 0.8f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-cube-subdiv") {
    auto fvshape = make_fvcube(0.075f);
    for (auto& p : fvshape.positions) p += {0, 0.075f, 0};
    return fvshape;
  } else if (type == "test-arealight1") {
    return shape_to_fvshape(make_rect({1, 1}, {0.2f, 0.2f}));
  } else if (type == "test-arealight2") {
    return shape_to_fvshape(make_rect({1, 1}, {0.2f, 0.2f}));
  } else if (type == "test-largearealight1") {
    return shape_to_fvshape(make_rect({1, 1}, {0.4f, 0.4f}));
  } else if (type == "test-largearealight2") {
    return shape_to_fvshape(make_rect({1, 1}, {0.4f, 0.4f}));
  } else if (type == "test-cloth") {
    return shape_to_fvshape(make_rect({64, 64}, {0.2f, 0.2f}));
  } else if (type == "test-clothy") {
    return shape_to_fvshape(make_recty({64, 64}, {0.2f, 0.2f}));
  } else {
    return {};
  }
}

// Load mesh
shape_data load_shape(
    const string& filename, string& error, bool flip_texcoord) {
  auto shape = shape_data{};
  if (!load_shape(filename, shape, error, flip_texcoord)) return shape_data{};
  return shape;
}
pair<io_status, shape_data> load_shape(
    const string& filename, bool flip_texcoord) {
  auto error = string{};
  auto shape = shape_data{};
  if (!load_shape(filename, shape, error, flip_texcoord))
    return {io_status{error}, {}};
  return {io_status{}, std::move(shape)};
}
io_status load_shape(
    const string& filename, shape_data& shape, bool flip_texcoord) {
  auto error = string{};
  if (!load_shape(filename, shape, error, flip_texcoord))
    return io_status{error};
  return io_status{};
}
io_status save_shape(const string& filename, const shape_data& shape,
    bool flip_texcoord, bool ascii) {
  auto error = string{};
  if (!save_shape(filename, shape, error, flip_texcoord, ascii))
    return io_status{error};
  return io_status{};
}

// Load mesh
pair<io_status, fvshape_data> load_fvshape(
    const string& filename, bool flip_texcoord) {
  auto error = string{};
  auto shape = fvshape_data{};
  if (!load_fvshape(filename, shape, error, flip_texcoord))
    return {io_status{error}, {}};
  return {io_status{}, std::move(shape)};
}
io_status load_fvshape(
    const string& filename, fvshape_data& fvshape, bool flip_texcoord) {
  auto error = string{};
  if (!load_fvshape(filename, fvshape, error, flip_texcoord))
    return io_status{error};
  return io_status{};
}
io_status save_fvshape(const string& filename, const fvshape_data& fvshape,
    bool flip_texcoord, bool ascii) {
  auto error = string{};
  if (!save_fvshape(filename, fvshape, error, flip_texcoord, ascii))
    return io_status{error};
  return io_status{};
}

// Shape presets used ofr testing.
bool make_shape_preset(
    const string& filename, shape_data& shape, string& error) {
  shape = make_shape_preset(path_basename(filename));
  if (shape.positions.empty()) {
    error = "unknown preset";
    return false;
  }
  return true;
}

// Shape presets used for testing.
bool make_fvshape_preset(
    const string& filename, fvshape_data& fvshape, string& error) {
  fvshape = make_fvshape_preset(path_basename(filename));
  if (fvshape.positions.empty()) {
    error = "unknown preset";
    return false;
  }
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// TEXTURE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Loads/saves an image. Chooses hdr or ldr based on file name.
bool load_texture(
    const string& filename, texture_data& texture, string& error) {
  auto read_error = [&]() {
    error = filename + ": rad error";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".exr" || ext == ".EXR") {
    auto pixels = (float*)nullptr;
    if (LoadEXR(&pixels, &texture.width, &texture.height, filename.c_str(),
            nullptr) != 0)
      return read_error();
    texture.linear  = true;
    texture.pixelsf = vector<vec4f>{
        (vec4f*)pixels, (vec4f*)pixels + texture.width * texture.height};
    free(pixels);
    return true;
  } else if (ext == ".hdr" || ext == ".HDR") {
    auto buffer = vector<byte>{};
    if (!load_binary(filename, buffer, error)) return false;
    auto ncomp  = 0;
    auto pixels = stbi_loadf_from_memory(buffer.data(), (int)buffer.size(),
        &texture.width, &texture.height, &ncomp, 4);
    if (!pixels) return read_error();
    texture.linear  = true;
    texture.pixelsf = vector<vec4f>{
        (vec4f*)pixels, (vec4f*)pixels + texture.width * texture.height};
    free(pixels);
    return true;
  } else if (ext == ".png" || ext == ".PNG") {
    auto buffer = vector<byte>{};
    if (!load_binary(filename, buffer, error)) return false;
    auto ncomp  = 0;
    auto pixels = stbi_load_from_memory(buffer.data(), (int)buffer.size(),
        &texture.width, &texture.height, &ncomp, 4);
    if (!pixels) return read_error();
    texture.linear  = false;
    texture.pixelsb = vector<vec4b>{
        (vec4b*)pixels, (vec4b*)pixels + texture.width * texture.height};
    free(pixels);
    return true;
  } else if (ext == ".jpg" || ext == ".JPG" || ext == ".jpeg" ||
             ext == ".JPEG") {
    auto buffer = vector<byte>{};
    if (!load_binary(filename, buffer, error)) return false;
    auto ncomp  = 0;
    auto pixels = stbi_load_from_memory(buffer.data(), (int)buffer.size(),
        &texture.width, &texture.height, &ncomp, 4);
    if (!pixels) return read_error();
    texture.linear  = false;
    texture.pixelsb = vector<vec4b>{
        (vec4b*)pixels, (vec4b*)pixels + texture.width * texture.height};
    free(pixels);
    return true;
  } else if (ext == ".tga" || ext == ".TGA") {
    auto buffer = vector<byte>{};
    if (!load_binary(filename, buffer, error)) return false;
    auto ncomp  = 0;
    auto pixels = stbi_load_from_memory(buffer.data(), (int)buffer.size(),
        &texture.width, &texture.height, &ncomp, 4);
    if (!pixels) return read_error();
    texture.linear  = false;
    texture.pixelsb = vector<vec4b>{
        (vec4b*)pixels, (vec4b*)pixels + texture.width * texture.height};
    free(pixels);
    return true;
  } else if (ext == ".bmp" || ext == ".BMP") {
    auto buffer = vector<byte>{};
    if (!load_binary(filename, buffer, error)) return false;
    auto ncomp  = 0;
    auto pixels = stbi_load_from_memory(buffer.data(), (int)buffer.size(),
        &texture.width, &texture.height, &ncomp, 4);
    if (!pixels) return read_error();
    texture.linear  = false;
    texture.pixelsb = vector<vec4b>{
        (vec4b*)pixels, (vec4b*)pixels + texture.width * texture.height};
    free(pixels);
    return true;
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    if (!make_texture_preset(filename, texture, error)) return false;
    return true;
  } else {
    error = filename + ": unknown format";
    return false;
  }
}

// Saves an hdr image.
bool save_texture(
    const string& filename, const texture_data& texture, string& error) {
  auto write_error = [&]() {
    error = filename + ": write error";
    return false;
  };

  // check for correct handling
  if (!texture.pixelsf.empty() && is_ldr_filename(filename))
    throw std::invalid_argument(
        filename + ": cannot save hdr texture to ldr file");
  if (!texture.pixelsb.empty() && is_hdr_filename(filename))
    throw std::invalid_argument(
        filename + ": cannot save ldr texture to hdr file");

  // write data
  auto stbi_write_data = [](void* context, void* data, int size) {
    auto& buffer = *(vector<byte>*)context;
    buffer.insert(buffer.end(), (byte*)data, (byte*)data + size);
  };

  auto ext = path_extension(filename);
  if (ext == ".hdr" || ext == ".HDR") {
    auto buffer = vector<byte>{};
    if (!stbi_write_hdr_to_func(stbi_write_data, &buffer, (int)texture.width,
            (int)texture.height, 4, (const float*)texture.pixelsf.data()))
      return write_error();
    if (!save_binary(filename, buffer, error)) return false;
    return true;
  } else if (ext == ".exr" || ext == ".EXR") {
    auto data = (byte*)nullptr;
    auto size = (size_t)0;
    if (SaveEXRToMemory((const float*)texture.pixelsf.data(),
            (int)texture.width, (int)texture.height, 4, 1, &data, &size,
            nullptr) < 0)
      return write_error();
    auto buffer = vector<byte>{data, data + size};
    free(data);
    if (!save_binary(filename, buffer, error)) return false;
    return true;
  } else if (ext == ".png" || ext == ".PNG") {
    auto buffer = vector<byte>{};
    if (!stbi_write_png_to_func(stbi_write_data, &buffer, (int)texture.width,
            (int)texture.height, 4, (const byte*)texture.pixelsb.data(),
            (int)texture.width * 4))
      return write_error();
    if (!save_binary(filename, buffer, error)) return false;
    return true;
  } else if (ext == ".jpg" || ext == ".JPG" || ext == ".jpeg" ||
             ext == ".JPEG") {
    auto buffer = vector<byte>{};
    if (!stbi_write_jpg_to_func(stbi_write_data, &buffer, (int)texture.width,
            (int)texture.height, 4, (const byte*)texture.pixelsb.data(), 75))
      return write_error();
    if (!save_binary(filename, buffer, error)) return false;
    return true;
  } else if (ext == ".tga" || ext == ".TGA") {
    auto buffer = vector<byte>{};
    if (!stbi_write_tga_to_func(stbi_write_data, &buffer, (int)texture.width,
            (int)texture.height, 4, (const byte*)texture.pixelsb.data()))
      return write_error();
    if (!save_binary(filename, buffer, error)) return false;
    return true;
  } else if (ext == ".bmp" || ext == ".BMP") {
    auto buffer = vector<byte>{};
    if (!stbi_write_bmp_to_func(stbi_write_data, &buffer, (int)texture.width,
            (int)texture.height, 4, (const byte*)texture.pixelsb.data()))
      return write_error();
    if (!save_binary(filename, buffer, error)) return false;
    return true;
  } else {
    error = filename + ": unknown format";
    return false;
  }
}

texture_data make_texture_preset(const string& type) {
  return image_to_texture(make_image_preset(type));
}

// Loads/saves an image. Chooses hdr or ldr based on file name.
pair<io_status, texture_data> load_texture(const string& filename) {
  auto error   = string{};
  auto texture = texture_data{};
  if (!load_texture(filename, texture, error)) return {io_status{error}, {}};
  return {io_status{}, std::move(texture)};
}
io_status load_texture(const string& filename, texture_data& texture) {
  auto error = string{};
  if (!load_texture(filename, texture, error)) return io_status{error};
  return io_status{};
}
io_status save_texture(const string& filename, const texture_data& texture) {
  auto error = string{};
  if (!save_texture(filename, texture, error)) return io_status{error};
  return io_status{};
}

bool make_texture_preset(
    const string& filename, texture_data& texture, string& error) {
  texture = make_texture_preset(path_basename(filename));
  if (texture.width == 0 || texture.height == 0) {
    error = "unknown preset";
    return false;
  }
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// make element name
[[maybe_unused]] static string get_element_name(
    const string& name, int idx, size_t size) {
  // there are much better ways to do this, but fine for now
  auto num_str  = std::to_string(idx + 1);
  auto size_str = std::to_string(size + 1);
  while (num_str.size() < size_str.size()) num_str = "0" + num_str;
  return name + num_str;
}

// get name
[[maybe_unused]] static string get_camera_name(
    const scene_data& scene, int idx) {
  if (idx < 0) return "";
  if (scene.camera_names.empty())
    return get_element_name("camera", idx, scene.cameras.size());
  return scene.camera_names[idx];
}
[[maybe_unused]] static string get_environment_name(
    const scene_data& scene, int idx) {
  if (idx < 0) return "";
  if (scene.environment_names.empty())
    return get_element_name("environment", idx, scene.environments.size());
  return scene.environment_names[idx];
}
[[maybe_unused]] static string get_shape_name(
    const scene_data& scene, int idx) {
  if (idx < 0) return "";
  if (scene.shape_names.empty())
    return get_element_name("shape", idx, scene.shapes.size());
  return scene.shape_names[idx];
}
[[maybe_unused]] static string get_texture_name(
    const scene_data& scene, int idx) {
  if (idx < 0) return "";
  if (scene.texture_names.empty())
    return get_element_name("texture", idx, scene.textures.size());
  return scene.texture_names[idx];
}
[[maybe_unused]] static string get_instance_name(
    const scene_data& scene, int idx) {
  if (idx < 0) return "";
  if (scene.instance_names.empty())
    return get_element_name("instance", idx, scene.instances.size());
  return scene.instance_names[idx];
}
[[maybe_unused]] static string get_material_name(
    const scene_data& scene, int idx) {
  if (idx < 0) return "";
  if (scene.material_names.empty())
    return get_element_name("material", idx, scene.materials.size());
  return scene.material_names[idx];
}
[[maybe_unused]] static string get_subdiv_name(
    const scene_data& scene, int idx) {
  if (idx < 0) return "";
  if (scene.subdiv_names.empty())
    return get_element_name("subdiv", idx, scene.subdivs.size());
  return scene.subdiv_names[idx];
}

[[maybe_unused]] static string get_camera_name(
    const scene_data& scene, const camera_data& camera) {
  return get_camera_name(scene, (int)(&camera - scene.cameras.data()));
}
[[maybe_unused]] static string get_environment_name(
    const scene_data& scene, const environment_data& environment) {
  return get_environment_name(
      scene, (int)(&environment - scene.environments.data()));
}
[[maybe_unused]] static string get_shape_name(
    const scene_data& scene, const shape_data& shape) {
  return get_shape_name(scene, (int)(&shape - scene.shapes.data()));
}
[[maybe_unused]] static string get_texture_name(
    const scene_data& scene, const texture_data& texture) {
  return get_texture_name(scene, (int)(&texture - scene.textures.data()));
}
[[maybe_unused]] static string get_instance_name(
    const scene_data& scene, const instance_data& instance) {
  return get_instance_name(scene, (int)(&instance - scene.instances.data()));
}
[[maybe_unused]] static string get_material_name(
    const scene_data& scene, const material_data& material) {
  return get_material_name(scene, (int)(&material - scene.materials.data()));
}
[[maybe_unused]] static string get_subdiv_name(
    const scene_data& scene, const subdiv_data& subdiv) {
  return get_subdiv_name(scene, (int)(&subdiv - scene.subdivs.data()));
}

template <typename T>
static vector<string> make_names(const vector<T>& elements,
    const vector<string>& names, const string& prefix) {
  if (names.size() == elements.size()) return names;
  auto nnames = vector<string>(elements.size());
  for (auto idx : range(elements.size())) {
    // there are much better ways to do this, but fine for now
    auto num_str  = std::to_string(idx + 1);
    auto size_str = std::to_string(elements.size());
    while (num_str.size() < size_str.size()) num_str = "0" + num_str;
    nnames[idx] = prefix + num_str;
  }
  return nnames;
}

// Add missing cameras.
void add_missing_camera(scene_data& scene) {
  if (!scene.cameras.empty()) return;
  scene.camera_names.emplace_back("camera");
  auto& camera        = scene.cameras.emplace_back();
  camera.orthographic = false;
  camera.film         = 0.036f;
  camera.aspect       = (float)16 / (float)9;
  camera.aperture     = 0;
  camera.lens         = 0.050f;
  auto bbox           = compute_bounds(scene);
  auto center         = (bbox.max + bbox.min) / 2;
  auto bbox_radius    = length(bbox.max - bbox.min) / 2;
  auto camera_dir     = vec3f{0, 0, 1};
  auto camera_dist = bbox_radius * camera.lens / (camera.film / camera.aspect);
  camera_dist *= 2.0f;  // correction for tracer camera implementation
  auto from    = camera_dir * camera_dist + center;
  auto to      = center;
  auto up      = vec3f{0, 1, 0};
  camera.frame = lookat_frame(from, to, up);
  camera.focus = length(from - to);
}

// Add missing radius.
static void add_missing_radius(scene_data& scene, float radius = 0.001f) {
  for (auto& shape : scene.shapes) {
    if (shape.points.empty() && shape.lines.empty()) continue;
    if (!shape.radius.empty()) continue;
    shape.radius.assign(shape.positions.size(), radius);
  }
}

// Add missing cameras.
void add_missing_material(scene_data& scene) {
  auto default_material = invalidid;
  for (auto& instance : scene.instances) {
    if (instance.material >= 0) continue;
    if (default_material == invalidid) {
      auto& material   = scene.materials.emplace_back();
      material.color   = {0.8f, 0.8f, 0.8f};
      default_material = (int)scene.materials.size() - 1;
    }
    instance.material = default_material;
  }
}

// Reduce memory usage
static void trim_memory(scene_data& scene) {
  for (auto& shape : scene.shapes) {
    shape.points.shrink_to_fit();
    shape.lines.shrink_to_fit();
    shape.triangles.shrink_to_fit();
    shape.quads.shrink_to_fit();
    shape.positions.shrink_to_fit();
    shape.normals.shrink_to_fit();
    shape.texcoords.shrink_to_fit();
    shape.colors.shrink_to_fit();
    shape.radius.shrink_to_fit();
    shape.tangents.shrink_to_fit();
  }
  for (auto& subdiv : scene.subdivs) {
    subdiv.positions.shrink_to_fit();
    subdiv.normals.shrink_to_fit();
    subdiv.texcoords.shrink_to_fit();
    subdiv.quadspos.shrink_to_fit();
    subdiv.quadsnorm.shrink_to_fit();
    subdiv.quadstexcoord.shrink_to_fit();
  }
  for (auto& texture : scene.textures) {
    texture.pixelsf.shrink_to_fit();
    texture.pixelsb.shrink_to_fit();
  }
  scene.cameras.shrink_to_fit();
  scene.shapes.shrink_to_fit();
  scene.subdivs.shrink_to_fit();
  scene.instances.shrink_to_fit();
  scene.materials.shrink_to_fit();
  scene.textures.shrink_to_fit();
  scene.environments.shrink_to_fit();
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// TEST SCENES
// -----------------------------------------------------------------------------
namespace yocto {

enum struct test_cameras_type { standard, wide };
enum struct test_environments_type { none, sky, sunsky };
enum struct test_arealights_type { none, standard, large };
enum struct test_floor_type { none, standard };
enum struct test_instance_name_type { material, shape };
enum struct test_shapes_type {
  // clang-format off
  features1, features2, rows, bunny_sphere,
  shapes1, shapes2, shapes3
  // clang-format off
};
enum struct test_materials_type {
  // clang-format off
  features1, features2, uvgrid, hair, plastic_metal,
  materials1, materials2, materials3, materials4, materials5,
  // clang-format on
};

struct test_params {
  test_cameras_type       cameras       = test_cameras_type::standard;
  test_environments_type  environments  = test_environments_type::sky;
  test_arealights_type    arealights    = test_arealights_type::standard;
  test_floor_type         floor         = test_floor_type::standard;
  test_shapes_type        shapes        = test_shapes_type::features1;
  test_materials_type     materials     = test_materials_type::features1;
  test_instance_name_type instance_name = test_instance_name_type::material;
};

// Scene test
scene_data make_test(const test_params& params) {
  return {};
  // // cameras
  // switch (params.cameras) {
  //   case test_cameras_type::standard: {
  //     add_camera(scene, "default", {-0.75, 0.4, 0.9}, {-0.075, 0.05, -0.05},
  //         {0, 1, 0}, 0.05, 2.4, 0);
  //   } break;
  //   // TODO(fabio): fix wide camera
  //   case test_cameras_type::wide: {
  //     add_camera(scene, "default", {-0.75, 0.4, 0.9}, {-0.075, 0.05, -0.05},
  //         {0, 1, 0}, 0.05, 2.4, 0);
  //   } break;
  // }
  // // TODO(fabio): port other cameras
  // switch (params.environments) {
  //   case test_environments_type::none: break;
  //   case test_environments_type::sky: {
  //     add_environment(scene, "sky", identity3x4f, {0.5, 0.5, 0.5},
  //         add_texture(scene, "sky",
  //             make_sunsky(
  //                 {2048, 1024}, pif / 4, 3.0, false, 1.0, 1.0, {0.7, 0.7,
  //                 0.7}),
  //             true));
  //   } break;
  //   case test_environments_type::sunsky: {
  //     add_environment(scene, "sunsky", identity3x4f, {0.5, 0.5, 0.5},
  //         add_texture(scene, "sky",
  //             make_sunsky(
  //                 {2048, 1024}, pif / 4, 3.0, true, 1.0, 1.0, {0.7, 0.7,
  //                 0.7}),
  //             true));
  //   } break;
  // }
  // switch (params.arealights) {
  //   case test_arealights_type::none: break;
  //   case test_arealights_type::standard: {
  //     add_instance(scene, "arealight1",
  //         lookat_frame({-0.4, 0.8, 0.8}, {0, 0.1, 0}, {0, 1, 0}, true),
  //         add_shape(scene, "arealight1", make_rect({1, 1}, {0.2, 0.2})),
  //         add_emission_material(
  //             scene, "arealight1", {20, 20, 20}, invalidid));
  //     add_instance(scene, "arealight2",
  //         lookat_frame({+0.4, 0.8, 0.8}, {0, 0.1, 0}, {0, 1, 0}, true),
  //         add_shape(scene, "arealight2", make_rect({1, 1}, {0.2, 0.2})),
  //         add_emission_material(
  //             scene, "arealight2", {20, 20, 20}, invalidid));
  //   } break;
  //   case test_arealights_type::large: {
  //     add_instance(scene, "largearealight1",
  //         lookat_frame({-0.8, 1.6, 1.6}, {0, 0.1, 0}, {0, 1, 0}, true),
  //         add_shape(scene, "largearealight1", make_rect({1, 1}, {0.4, 0.4})),
  //         add_emission_material(
  //             scene, "largearealight1", {10, 10, 10}, invalidid));
  //     add_instance(scene, "largearealight2",
  //         lookat_frame({+0.8, 1.6, 1.6}, {0, 0.1, 0}, {0, 1, 0}, true),
  //         add_shape(scene, "largearealight2", make_rect({1, 1}, {0.4, 0.4})),
  //         add_emission_material(
  //             scene, "largearealight2", {10, 10, 10}, invalidid));
  //   } break;
  // }
  // switch (params.floor) {
  //   case test_floor_type::none: break;
  //   case test_floor_type::standard: {
  //     add_instance(scene, "floor", identity3x4f,
  //         add_shape(scene, "floor", make_floor({1, 1}, {2, 2}, {20, 20})),
  //         add_matte_material(scene, "floor", {1, 1, 1},
  //             add_texture(scene, "floor", make_grid({1024, 1024}))));
  //   } break;
  // }
  // auto shapes = vector<int>{}, shapesi =
  // vector<int>{}; auto subdivs   =
  // vector<int>{}; auto materials =
  // vector<int>{}; switch (params.shapes) {
  //   case test_shapes_type::features1: {
  //     auto bunny  = add_shape(scene, "sphere", make_sphere(32, 0.075, 1));
  //     auto sphere = add_shape(scene, "sphere", make_sphere(32, 0.075, 1));
  //     shapes      = {bunny, sphere, bunny, sphere, bunny};
  //   } break;
  //   case test_shapes_type::features2: {
  //     shapes  = {add_shape(scene, "sphere", make_sphere(32, 0.075, 1)),
  //         add_shape(scene, "suzanne", make_monkey(0.075f * 0.8f)),
  //         add_shape(scene, "hair",
  //             make_hair(make_sphere(32, 0.075f * 0.8f, 1), {4, 65536},
  //                 {0.1f * 0.15f, 0.1f * 0.15f},
  //                 {0.001f * 0.15f, 0.0005f * 0.15f}, {0.03, 100})),
  //         add_shape(scene, "displaced", make_sphere(128, 0.075f, 1)),
  //         add_shape(scene, "cube",
  //             make_rounded_box({32, 32, 32}, {0.075, 0.075, 0.075}, {1, 1,
  //             1},
  //                 0.3 * 0.075f))};
  //     shapesi = {invalidid, invalidid,
  //         add_shape(scene, "hairi", make_sphere(32, 0.075f * 0.8f, 1)),
  //         invalidid, invalidid};
  //     subdivs = {add_subdiv(scene, "suzanne", make_monkey(0.075f * 0.8f),
  //                    shapes[1], 2),
  //         add_subdiv(scene, "displaced", make_sphere(128, 0.075f, 1),
  //         shapes[3],
  //             0, 0.025,
  //             add_texture(scene, "bumps-displacement", make_bumps({1024,
  //             1024}),
  //                 false, true))};
  //   } break;
  //   case test_shapes_type::rows: {
  //     auto bunny  = add_shape(scene, "bunny", make_sphere(32, 0.075, 1));
  //     auto sphere = add_shape(scene, "sphere", make_sphere(32, 0.075, 1));
  //     shapes      = {bunny, bunny, bunny, bunny, bunny, sphere, sphere,
  //     sphere,
  //         sphere, sphere};
  //   } break;
  //   case test_shapes_type::bunny_sphere: {
  //     auto bunny  = add_shape(scene, "bunny", make_sphere(32, 0.075, 1));
  //     auto sphere = add_shape(scene, "sphere", make_sphere(32, 0.075, 1));
  //     shapes      = {bunny, sphere, bunny, sphere, bunny};
  //   } break;
  //   case test_shapes_type::shapes1: {
  //     shapes = {
  //         add_shape(scene, "sphere", make_sphere(32, 0.075, 1)),
  //         add_shape(scene, "uvsphere-flipcap",
  //             make_capped_uvsphere({32, 32}, 0.075, {1, 1}, 0.3 * 0.075)),
  //         add_shape(scene, "disk", make_disk(32, 0.075f, 1)),
  //         add_shape(scene, "uvcylinder",
  //             make_rounded_uvcylinder(
  //                 {32, 32, 32}, {0.075, 0.075}, {1, 1, 1}, 0.3 * 0.075)),
  //         add_shape(scene, "cube",
  //             make_rounded_box({32, 32, 32}, {0.075, 0.075, 0.075}, {1, 1,
  //             1},
  //                 0.3 * 0.075f)),
  //     };
  //   } break;
  //   case test_shapes_type::shapes2: {
  //     shapes = {
  //         add_shape(scene, "cube-subdiv", make_fvcube(0.075)),
  //         add_shape(scene, "suzanne-subdiv", make_monkey(0.075)),
  //         add_shape(scene, "displaced", make_sphere(128, 0.075f, 1)),
  //         add_shape(scene, "bunny", make_sphere(32, 0.075, 1)),
  //         add_shape(scene, "teapot", make_sphere(32, 0.075, 1)),
  //     };
  //     subdivs = {
  //         add_subdiv(scene, "cube-subdiv", make_fvcube(0.075), shapes[0], 4),
  //         add_subdiv(scene, "suzanne-subdiv", make_monkey(0.075), shapes[1],
  //         2), add_subdiv(scene, "displaced", make_sphere(128, 0.075f, 1),
  //         shapes[2],
  //             0, 0.025,
  //             add_texture(scene, "bumps-displacement", make_bumps({1024,
  //             1024}),
  //                 false, true))};
  //   } break;
  //   case test_shapes_type::shapes3: {
  //     shapes = {
  //         invalidid,
  //         add_shape(scene, "hair1",
  //             make_hair(make_sphere(32, 0.075f * 0.8f, 1), {4, 65536},
  //                 {0.1f * 0.15f, 0.1f * 0.15f},
  //                 {0.001f * 0.15f, 0.0005f * 0.15f}, {0.03, 100})),
  //         add_shape(scene, "hair2",
  //             make_hair(make_sphere(32, 0.075f * 0.8f, 1), {4, 65536},
  //                 {0.1f * 0.15f, 0.1f * 0.15f},
  //                 {0.001f * 0.15f, 0.0005f * 0.15f})),
  //         add_shape(scene, "hair3",
  //             make_hair(make_sphere(32, 0.075f * 0.8f, 1), {4, 65536},
  //                 {0.1f * 0.15f, 0.1f * 0.15f},
  //                 {0.001f * 0.15f, 0.0005f * 0.15f}, {0, 0}, {0.5, 128})),
  //         invalidid,
  //     };
  //   } break;
  // }
  // switch (params.materials) {
  //   case test_materials_type::features1: {
  //     materials = {
  //         add_plastic_material(scene, "coated", {1, 1, 1}, 0.2,
  //             add_texture(scene, "uvgrid", make_uvgrid({1024, 1024}))),
  //         add_glass_material(scene, "glass", {1, 0.5, 0.5}, 0),
  //         add_glass_material(
  //             scene, "jade", {0.5, 0.5, 0.5}, 0, {0.3, 0.6, 0.3}),
  //         add_plastic_material(scene, "bumped", {0.5, 0.7, 0.5}, 0.2,
  //             invalidid, invalidid,
  //             add_texture(scene, "bumps-normal",
  //                 bump_to_normal(make_bumps({1024, 1024}), 0.05), false,
  //                 true)),
  //         add_metal_material(scene, "metal", {0.66, 0.45, 0.34}, 0.2),
  //     };
  //   } break;
  //   case test_materials_type::features2: {
  //     auto uvgrid  = add_plastic_material(scene, "uvgrid", {1, 1, 1}, 0.2,
  //         add_texture(scene, "uvgrid", make_uvgrid({1024, 1024})));
  //     auto plastic = add_plastic_material(
  //         scene, "plastic", {0.5, 0.7, 0.5}, 0.2);
  //     auto hair = add_matte_material(scene, "hair", {0.7, 0.7, 0.7});
  //     materials = {uvgrid, plastic, hair, plastic, uvgrid};
  //   } break;
  //   case test_materials_type::uvgrid: {
  //     auto uvgrid = add_plastic_material(scene, "uvgrid", {1, 1, 1}, 0.2,
  //         add_texture(scene, "uvgrid", make_uvgrid({1024, 1024})));
  //     materials   = {uvgrid, uvgrid, uvgrid, uvgrid, uvgrid};
  //   } break;
  //   case test_materials_type::hair: {
  //     auto hair = add_matte_material(scene, "hair", {0.7, 0.7, 0.7});
  //     materials = {hair, hair, hair, hair, hair};
  //   } break;
  //   case test_materials_type::plastic_metal: {
  //     materials = {
  //         add_plastic_material(scene, "plastic1", {0.5, 0.5, 0.7}, 0.01),
  //         add_plastic_material(scene, "plastic2", {0.5, 0.7, 0.5}, 0.2),
  //         add_matte_material(scene, "matte", {0.7, 0.7, 0.7}),
  //         add_metal_material(scene, "metal1", {0.7, 0.7, 0.7}, 0),
  //         add_metal_material(scene, "metal2", {0.66, 0.45, 0.34}, 0.2),
  //     };
  //   } break;
  //   case test_materials_type::materials1: {
  //     materials = {
  //         add_plastic_material(scene, "plastic1", {0.5, 0.5, 0.7}, 0.01),
  //         add_plastic_material(scene, "plastic2", {0.5, 0.7, 0.5}, 0.2),
  //         add_matte_material(scene, "matte", {0.7, 0.7, 0.7}),
  //         add_plastic_material(scene, "metal1", {0.7, 0.7, 0.7}, 0),
  //         add_plastic_material(scene, "metal2", {0.66, 0.45, 0.34}, 0.2),
  //     };
  //   } break;
  //   case test_materials_type::materials2: {
  //     materials = {
  //         add_glass_material(scene, "glass1", {1, 1, 1}, 0),
  //         add_glass_material(scene, "glass2", {1, 0.7, 0.7}, 0.1),
  //         add_transparent_material(scene, "transparent", {0.7, 0.5, 0.5},
  //         0.2), add_thinglass_material(scene, "tglass1", {1, 1, 1}, 0),
  //         add_thinglass_material(scene, "tglass2", {1, 0.7, 0.7}, 0.1),
  //     };
  //   } break;
  //   case test_materials_type::materials3: {
  //     auto bumps_normal = add_texture(scene, "bumps-normal",
  //         bump_to_normal(make_bumps({1024, 1024}), 0.05), false, true);
  //     materials         = {
  //         add_plastic_material(scene, "plastic1", {0.5, 0.5, 0.7}, 0.01,
  //             invalidid, invalidid, bumps_normal),
  //         add_plastic_material(scene, "plastic2", {0.5, 0.7, 0.5}, 0.2),
  //         add_metal_material(scene, "metal1", {0.7, 0.7, 0.7}, 0,
  //             invalidid, invalidid, bumps_normal),
  //         add_metal_material(scene, "metal2", {0.66, 0.45, 0.34}, 0.2),
  //         add_metal_material(scene, "metal3", {0.66, 0.45, 0.34}, 0.2),
  //     };
  //   } break;
  //   case test_materials_type::materials4: {
  //     materials = {
  //         add_volume_material(
  //             scene, "cloud", {0.65, 0.65, 0.65}, {0.9, 0.9, 0.9}, 1),
  //         add_glass_material(scene, "glass", {1, 0.5, 0.5}, 0),
  //         add_glass_material(
  //             scene, "jade", {0.5, 0.5, 0.5}, 0, {0.3, 0.6, 0.3}),
  //         add_glass_material(
  //             scene, "jade2", {0.5, 0.5, 0.5}, 0, {0.3, 0.6, 0.3}),
  //         add_volume_material(scene, "smoke", {0.5, 0.5, 0.5}, {0.2, 0.2,
  //         0.2}),
  //     };
  //   } break;
  //   case test_materials_type::materials5: {
  //     materials = {
  //         add_glass_material(scene, "skin1a", {0.76, 0.48, 0.23}, 0.25,
  //             {0.436, 0.227, 0.131}, invalidid, invalidid,
  //             invalidid, 1.5, -0.8, 0.001),
  //         add_glass_material(scene, "skin2a", {0.82, 0.55, 0.4}, 0.25,
  //             {0.623, 0.433, 0.343}, invalidid, invalidid,
  //             invalidid, 1.5, -0.8, 0.001),
  //         add_glass_material(scene, "skins", {0.76, 0.48, 0.23}, 0,
  //             {0.436, 0.227, 0.131}, invalidid, invalidid,
  //             invalidid, 1.5, -0.8, 0.001),
  //         add_glass_material(scene, "skin1b", {0.76, 0.48, 0.23}, 0.25,
  //             {0.436, 0.227, 0.131}, invalidid, invalidid,
  //             invalidid, 1.5, -0.8, 0.001),
  //         add_glass_material(scene, "skin2b", {0.82, 0.55, 0.4}, 0.25,
  //             {0.623, 0.433, 0.343}, invalidid, invalidid,
  //             invalidid, 1.5, -0.8, 0.001),
  //     };
  //   } break;
  // }
  // for (auto idx = 0; idx < shapes.size(); idx++) {
  //   if (!shapes[idx]) continue;
  //   if (shapes.size() > 5) {
  //     add_instance(scene,
  //         scene.shape_names[idx] + "-" + scene.shape_names[idx % 5],
  //         {{1, 0, 0}, {0, 1, 0}, {0, 0, 1},
  //             {0.2f * (idx % 5 - 2), 0.075, -0.4f * (idx / 5)}},
  //         shapes[idx], materials[idx % 5]);
  //   } else {
  //     auto name = params.instance_name == test_instance_name_type::material
  //                     ? scene.material_names[idx]
  //                     : scene.shape_names[idx];
  //     add_instance(scene, name,
  //         {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.2f * (idx % 5 - 2), 0.075,
  //         0}}, shapes[idx], materials[idx]);
  //   }
  //   if (!shapesi.empty() && shapesi[idx]) {
  //     // TODO(fabio): fix name
  //     add_instance(scene, "",
  //         {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.2f * (idx - 2), 0.075, 0}},
  //         shapesi[idx], materials[idx]);
  //   }
  // }
}

// Scene presets used for testing.
scene_data make_scene_preset(const string& type) {
  if (type == "cornellbox") {
    return make_cornellbox();
  } else if (type == "features1") {
    return make_test({test_cameras_type::standard, test_environments_type::sky,
        test_arealights_type::standard, test_floor_type::standard,
        test_shapes_type::features1, test_materials_type::features1,
        test_instance_name_type::material});
  } else if (type == "features2") {
    return make_test({test_cameras_type::standard, test_environments_type::sky,
        test_arealights_type::standard, test_floor_type::standard,
        test_shapes_type::features2, test_materials_type::features2,
        test_instance_name_type::shape});
  } else if (type == "materials1") {
    return make_test({test_cameras_type::wide, test_environments_type::sky,
        test_arealights_type::large, test_floor_type::standard,
        test_shapes_type::rows, test_materials_type::materials1,
        test_instance_name_type::material});
  } else if (type == "materials2") {
    return make_test({test_cameras_type::wide, test_environments_type::sky,
        test_arealights_type::large, test_floor_type::standard,
        test_shapes_type::rows, test_materials_type::materials2,
        test_instance_name_type::material});
  } else if (type == "materials3") {
    return make_test({test_cameras_type::wide, test_environments_type::sky,
        test_arealights_type::large, test_floor_type::standard,
        test_shapes_type::rows, test_materials_type::materials3,
        test_instance_name_type::material});
  } else if (type == "materials4") {
    return make_test({test_cameras_type::wide, test_environments_type::sky,
        test_arealights_type::large, test_floor_type::standard,
        test_shapes_type::rows, test_materials_type::materials4,
        test_instance_name_type::material});
  } else if (type == "materials5") {
    return make_test({test_cameras_type::wide, test_environments_type::sky,
        test_arealights_type::large, test_floor_type::standard,
        test_shapes_type::rows, test_materials_type::materials5,
        test_instance_name_type::material});
  } else if (type == "shapes1") {
    return make_test({test_cameras_type::standard, test_environments_type::sky,
        test_arealights_type::large, test_floor_type::standard,
        test_shapes_type::shapes1, test_materials_type::uvgrid,
        test_instance_name_type::shape});
  } else if (type == "shapes2") {
    return make_test({test_cameras_type::standard, test_environments_type::sky,
        test_arealights_type::large, test_floor_type::standard,
        test_shapes_type::shapes2, test_materials_type::uvgrid,
        test_instance_name_type::shape});
  } else if (type == "shapes3") {
    return make_test({test_cameras_type::standard, test_environments_type::sky,
        test_arealights_type::large, test_floor_type::standard,
        test_shapes_type::shapes3, test_materials_type::hair,
        test_instance_name_type::shape});
  } else if (type == "environments1") {
    return make_test({test_cameras_type::standard, test_environments_type::sky,
        test_arealights_type::none, test_floor_type::standard,
        test_shapes_type::bunny_sphere, test_materials_type::plastic_metal,
        test_instance_name_type::material});
  } else if (type == "environments2") {
    return make_test({test_cameras_type::standard,
        test_environments_type::sunsky, test_arealights_type::none,
        test_floor_type::standard, test_shapes_type::bunny_sphere,
        test_materials_type::plastic_metal, test_instance_name_type::material});
  } else if (type == "arealights1") {
    return make_test({test_cameras_type::standard, test_environments_type::none,
        test_arealights_type::standard, test_floor_type::standard,
        test_shapes_type::bunny_sphere, test_materials_type::plastic_metal,
        test_instance_name_type::material});
  } else {
    return {};
  }
}

// Scene presets used for testing.
bool make_scene_preset(
    const string& filename, scene_data& scene, string& error) {
  auto type = path_basename(filename);
  if (type == "cornellbox") {
    scene = make_cornellbox();
    return true;
  } else if (type == "features1") {
    scene = make_test({test_cameras_type::standard, test_environments_type::sky,
        test_arealights_type::standard, test_floor_type::standard,
        test_shapes_type::features1, test_materials_type::features1,
        test_instance_name_type::material});
    return true;
  } else if (type == "features2") {
    scene = make_test({test_cameras_type::standard, test_environments_type::sky,
        test_arealights_type::standard, test_floor_type::standard,
        test_shapes_type::features2, test_materials_type::features2,
        test_instance_name_type::shape});
    return true;
  } else if (type == "materials1") {
    scene = make_test({test_cameras_type::wide, test_environments_type::sky,
        test_arealights_type::large, test_floor_type::standard,
        test_shapes_type::rows, test_materials_type::materials1,
        test_instance_name_type::material});
    return true;
  } else if (type == "materials2") {
    scene = make_test({test_cameras_type::wide, test_environments_type::sky,
        test_arealights_type::large, test_floor_type::standard,
        test_shapes_type::rows, test_materials_type::materials2,
        test_instance_name_type::material});
    return true;
  } else if (type == "materials3") {
    scene = make_test({test_cameras_type::wide, test_environments_type::sky,
        test_arealights_type::large, test_floor_type::standard,
        test_shapes_type::rows, test_materials_type::materials3,
        test_instance_name_type::material});
    return true;
  } else if (type == "materials4") {
    scene = make_test({test_cameras_type::wide, test_environments_type::sky,
        test_arealights_type::large, test_floor_type::standard,
        test_shapes_type::rows, test_materials_type::materials4,
        test_instance_name_type::material});
    return true;
  } else if (type == "materials5") {
    scene = make_test({test_cameras_type::wide, test_environments_type::sky,
        test_arealights_type::large, test_floor_type::standard,
        test_shapes_type::rows, test_materials_type::materials5,
        test_instance_name_type::material});
    return true;
  } else if (type == "shapes1") {
    scene = make_test({test_cameras_type::standard, test_environments_type::sky,
        test_arealights_type::large, test_floor_type::standard,
        test_shapes_type::shapes1, test_materials_type::uvgrid,
        test_instance_name_type::shape});
    return true;
  } else if (type == "shapes2") {
    scene = make_test({test_cameras_type::standard, test_environments_type::sky,
        test_arealights_type::large, test_floor_type::standard,
        test_shapes_type::shapes2, test_materials_type::uvgrid,
        test_instance_name_type::shape});
    return true;
  } else if (type == "shapes3") {
    scene = make_test({test_cameras_type::standard, test_environments_type::sky,
        test_arealights_type::large, test_floor_type::standard,
        test_shapes_type::shapes3, test_materials_type::hair,
        test_instance_name_type::shape});
    return true;
  } else if (type == "environments1") {
    scene = make_test({test_cameras_type::standard, test_environments_type::sky,
        test_arealights_type::none, test_floor_type::standard,
        test_shapes_type::bunny_sphere, test_materials_type::plastic_metal,
        test_instance_name_type::material});
    return true;
  } else if (type == "environments2") {
    scene = make_test({test_cameras_type::standard,
        test_environments_type::sunsky, test_arealights_type::none,
        test_floor_type::standard, test_shapes_type::bunny_sphere,
        test_materials_type::plastic_metal, test_instance_name_type::material});
    return true;
  } else if (type == "arealights1") {
    scene = make_test({test_cameras_type::standard,
        test_environments_type::none, test_arealights_type::standard,
        test_floor_type::standard, test_shapes_type::bunny_sphere,
        test_materials_type::plastic_metal, test_instance_name_type::material});
    return true;
  } else {
    error = "unknown preset";
    return false;
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GENERIC SCENE LOADING
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a scene in the builtin JSON format.
static bool load_json_scene(
    const string& filename, scene_data& scene, string& error, bool noparallel);
static bool save_json_scene(const string& filename, const scene_data& scene,
    string& error, bool noparallel);

// Load/save a scene from/to OBJ.
static bool load_obj_scene(
    const string& filename, scene_data& scene, string& error, bool noparallel);
static bool save_obj_scene(const string& filename, const scene_data& scene,
    string& error, bool noparallel);

// Load/save a scene from/to PLY. Loads/saves only one mesh with no other data.
static bool load_ply_scene(
    const string& filename, scene_data& scene, string& error, bool noparallel);
static bool save_ply_scene(const string& filename, const scene_data& scene,
    string& error, bool noparallel);

// Load/save a scene from/to STL. Loads/saves only one mesh with no other data.
static bool load_stl_scene(
    const string& filename, scene_data& scene, string& error, bool noparallel);
static bool save_stl_scene(const string& filename, const scene_data& scene,
    string& error, bool noparallel);

// Load/save a scene from/to glTF.
static bool load_gltf_scene(
    const string& filename, scene_data& scene, string& error, bool noparallel);
static bool save_gltf_scene(const string& filename, const scene_data& scene,
    string& error, bool noparallel);

// Load/save a scene from/to pbrt-> This is not robust at all and only
// works on scene that have been previously adapted since the two renderers
// are too different to match.
static bool load_pbrt_scene(
    const string& filename, scene_data& scene, string& error, bool noparallel);
static bool save_pbrt_scene(const string& filename, const scene_data& scene,
    string& error, bool noparallel);

// Load a scene
bool load_scene(
    const string& filename, scene_data& scene, string& error, bool noparallel) {
  auto ext = path_extension(filename);
  if (ext == ".json" || ext == ".JSON") {
    return load_json_scene(filename, scene, error, noparallel);
  } else if (ext == ".obj" || ext == ".OBJ") {
    return load_obj_scene(filename, scene, error, noparallel);
  } else if (ext == ".gltf" || ext == ".GLTF") {
    return load_gltf_scene(filename, scene, error, noparallel);
  } else if (ext == ".pbrt" || ext == ".PBRT") {
    return load_pbrt_scene(filename, scene, error, noparallel);
  } else if (ext == ".ply" || ext == ".PLY") {
    return load_ply_scene(filename, scene, error, noparallel);
  } else if (ext == ".stl" || ext == ".STL") {
    return load_stl_scene(filename, scene, error, noparallel);
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    return make_scene_preset(filename, scene, error);
  } else {
    error = filename + ": unknown format";
    return false;
  }
}

// Save a scene
bool save_scene(const string& filename, const scene_data& scene, string& error,
    bool noparallel) {
  auto ext = path_extension(filename);
  if (ext == ".json" || ext == ".JSON") {
    return save_json_scene(filename, scene, error, noparallel);
  } else if (ext == ".obj" || ext == ".OBJ") {
    return save_obj_scene(filename, scene, error, noparallel);
  } else if (ext == ".pbrt" || ext == ".PBRT") {
    return save_pbrt_scene(filename, scene, error, noparallel);
  } else if (ext == ".gltf" || ext == ".GLTF") {
    return save_gltf_scene(filename, scene, error, noparallel);
  } else if (ext == ".ply" || ext == ".PLY") {
    return save_ply_scene(filename, scene, error, noparallel);
  } else if (ext == ".stl" || ext == ".STL") {
    return save_stl_scene(filename, scene, error, noparallel);
  } else {
    error = filename + ": unknown format";
    return false;
  }
}

// Load/save a scene
pair<io_status, scene_data> load_scene(
    const string& filename, bool noparallel) {
  auto error = string{};
  auto scene = scene_data{};
  if (!load_scene(filename, scene, error, noparallel))
    return {io_status{error}, {}};
  return {io_status{}, std::move(scene)};
}
io_status load_scene(
    const string& filename, scene_data& scene, bool noparallel) {
  auto error = string{};
  if (!load_scene(filename, scene, error, noparallel)) return io_status{error};
  return io_status{};
}
io_status save_scene(
    const string& filename, const scene_data& scene, bool noparallel) {
  auto error = string{};
  if (!save_scene(filename, scene, error, noparallel)) return io_status{error};
  return io_status{};
}

// Make missing scene directories
bool make_scene_directories(
    const string& filename, const scene_data& scene, string& error) {
  // make a directory if needed
  if (!make_directory(path_dirname(filename), error)) return false;
  if (!scene.shapes.empty())
    if (!make_directory(path_join(path_dirname(filename), "shapes"), error))
      return false;
  if (!scene.textures.empty())
    if (!make_directory(path_join(path_dirname(filename), "textures"), error))
      return false;
  if (!scene.subdivs.empty())
    if (!make_directory(path_join(path_dirname(filename), "subdivs"), error))
      return false;
  return true;
}

// Add environment
bool add_environment(scene_data& scene, const string& filename, string& error) {
  auto texture = texture_data{};
  if (!load_texture(filename, texture, error)) return false;
  scene.textures.push_back(std::move(texture));
  scene.environments.push_back(
      {identity3x4f, {1, 1, 1}, (int)scene.textures.size() - 1});
  return true;
}

// Make missing scene directories
io_status make_scene_directories(
    const string& filename, const scene_data& scene) {
  auto error = string{};
  if (!make_scene_directories(filename, scene, error)) return io_status{error};
  return io_status{};
}

// Add environment
io_status add_environment(scene_data& scene, const string& filename) {
  auto error = string{};
  if (!add_environment(scene, filename, error)) return io_status{error};
  return io_status{};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// INDIVIDUAL ELEMENTS
// -----------------------------------------------------------------------------
namespace yocto {

// load instances
static bool load_instance(
    const string& filename, vector<frame3f>& frames, string& error) {
  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = ply_model{};
    if (!load_ply(filename, ply, error)) return false;
    if (!get_values(ply, "instance",
            {"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz", "ox", "oy",
                "oz"},
            frames)) {
      error = filename + ": parse error";
      return false;
    }
    return true;
  } else {
    error = filename + ": unknown format";
    return false;
  }
}

// save instances
[[maybe_unused]] static bool save_instance(const string& filename,
    const vector<frame3f>& frames, string& error, bool ascii = false) {
  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = ply_model{};
    add_values(ply, "instance",
        {"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz", "ox", "oy",
            "oz"},
        frames);
    if (!save_ply(filename, ply, error)) return false;
    return true;
  } else {
    error = filename + ": unknown format";
    return false;
  }
}

// load subdiv
bool load_subdiv(const string& filename, subdiv_data& subdiv, string& error) {
  auto lsubdiv = fvshape_data{};
  if (!load_fvshape(filename, lsubdiv, error, true)) return false;
  subdiv.quadspos      = lsubdiv.quadspos;
  subdiv.quadsnorm     = lsubdiv.quadsnorm;
  subdiv.quadstexcoord = lsubdiv.quadstexcoord;
  subdiv.positions     = lsubdiv.positions;
  subdiv.normals       = lsubdiv.normals;
  subdiv.texcoords     = lsubdiv.texcoords;
  return true;
}

// save subdiv
bool save_subdiv(
    const string& filename, const subdiv_data& subdiv, string& error) {
  auto ssubdiv          = fvshape_data{};
  ssubdiv.quadspos      = subdiv.quadspos;
  ssubdiv.quadsnorm     = subdiv.quadsnorm;
  ssubdiv.quadstexcoord = subdiv.quadstexcoord;
  ssubdiv.positions     = subdiv.positions;
  ssubdiv.normals       = subdiv.normals;
  ssubdiv.texcoords     = subdiv.texcoords;
  if (!save_fvshape(filename, ssubdiv, error, true)) return false;
  return true;
}

// load/save subdiv
pair<io_status, subdiv_data> load_subdiv(const string& filename) {
  auto error  = string{};
  auto subdiv = subdiv_data{};
  if (!load_subdiv(filename, subdiv, error)) return {io_status{error}, {}};
  return {io_status{}, std::move(subdiv)};
}
io_status load_subdiv(const string& filename, subdiv_data& subdiv) {
  auto error = string{};
  if (!load_subdiv(filename, subdiv, error)) return io_status{error};
  return io_status{};
}
io_status save_subdiv(const string& filename, const subdiv_data& subdiv) {
  auto error = string{};
  if (!save_subdiv(filename, subdiv, error)) return io_status{error};
  return io_status{};
}

// save binary shape
static bool save_binshape(
    const string& filename, const shape_data& shape, string& error) {
  auto write_values = [](vector<byte>& buffer, const auto& values) {
    if (values.empty()) return;
    buffer.insert(buffer.end(), (byte*)values.data(),
        (byte*)values.data() + values.size() * sizeof(values.front()));
  };

  auto buffer = vector<byte>{};

  write_values(buffer, shape.positions);
  write_values(buffer, shape.normals);
  write_values(buffer, shape.texcoords);
  write_values(buffer, shape.colors);
  write_values(buffer, shape.radius);
  write_values(buffer, shape.points);
  write_values(buffer, shape.lines);
  write_values(buffer, shape.triangles);
  write_values(buffer, quads_to_triangles(shape.quads));

  if (!save_binary(filename, buffer, error)) return false;
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON IO
// -----------------------------------------------------------------------------
namespace yocto {

// Material type
enum struct material_type40 {
  // clang-format off
  matte, glossy, metallic, transparent, refractive, subsurface, volume, gltfpbr
  // clang-format on
};

// Enum labels
static const auto material_type40_names = std::vector<std::string>{"matte",
    "glossy", "metallic", "transparent", "refractive", "subsurface", "volume",
    "gltfpbr"};

NLOHMANN_JSON_SERIALIZE_ENUM(
    material_type40, {
                         {material_type40::matte, "matte"},
                         {material_type40::glossy, "glossy"},
                         {material_type40::metallic, "metallic"},
                         {material_type40::transparent, "transparent"},
                         {material_type40::refractive, "refractive"},
                         {material_type40::subsurface, "subsurface"},
                         {material_type40::volume, "volume"},
                         {material_type40::gltfpbr, "gltfpbr"},
                     })
NLOHMANN_JSON_SERIALIZE_ENUM(
    material_type, {
                       {material_type::matte, "matte"},
                       {material_type::glossy, "glossy"},
                       {material_type::reflective, "reflective"},
                       {material_type::transparent, "transparent"},
                       {material_type::refractive, "refractive"},
                       {material_type::subsurface, "subsurface"},
                       {material_type::volumetric, "volumetric"},
                       {material_type::gltfpbr, "gltfpbr"},
                   })

// Load a scene in the builtin JSON format.
static bool load_json_scene_version40(const string& filename,
    const json_value& json, scene_data& scene, string& error, bool noparallel) {
  auto parse_error = [filename, &error](const string& patha,
                         const string& pathb = "", const string& pathc = "") {
    auto path = patha;
    if (!pathb.empty()) path += "/" + pathb;
    if (!pathc.empty()) path += "/" + pathc;
    error = filename + ": parse error at " + path;
    return false;
  };
  auto key_error = [filename, &error](const string& patha,
                       const string& pathb = "", const string& pathc = "") {
    auto path = patha;
    if (!pathb.empty()) path += "/" + pathb;
    if (!pathc.empty()) path += "/" + pathc;
    error = filename + "; unknow key at " + path;
    return false;
  };

  // parse json value
  auto get_opt = [](const json_value& json, const string& key, auto& value) {
    value = json.value(key, value);
  };
  auto get_of3 = [](const json_value& json, const string& key, auto& value) {
    auto valuea = json.value(key, (array<float, 12>&)value);
    value       = *(frame3f*)&valuea;
  };
  auto get_om3 = [](const json_value& json, const string& key, auto& value) {
    auto valuea = json.value(key, (array<float, 9>&)value);
    value       = *(mat3f*)&valuea;
  };

  // parse json reference
  auto shape_map = unordered_map<string, int>{};
  auto get_shp   = [&scene, &shape_map](
                     const json_value& json, const string& key, int& value) {
    auto name = json.value(key, string{});
    if (name.empty()) return;
    auto it = shape_map.find(name);
    if (it != shape_map.end()) {
      value = it->second;
    } else {
      scene.shape_names.emplace_back(name);
      scene.shapes.emplace_back();
      auto shape_id   = (int)scene.shapes.size() - 1;
      shape_map[name] = shape_id;
      value           = shape_id;
    }
  };

  // parse json reference
  auto material_map = unordered_map<string, int>{};
  auto get_mat      = [&material_map](
                     const json_value& json, const string& key, int& value) {
    auto name = json.value(key, string{});
    if (name.empty()) return;
    auto it = material_map.find(name);
    if (it != material_map.end()) {
      value = it->second;
    } else {
      throw std::out_of_range{"missing key"};
    }
  };

  // parse json reference
  auto texture_map = unordered_map<string, int>{};
  auto get_tex     = [&scene, &texture_map](
                     const json_value& json, const string& key, int& value) {
    auto name = json.value(key, string{});
    if (name.empty()) return;
    auto it = texture_map.find(name);
    if (it != texture_map.end()) {
      value = it->second;
    } else {
      scene.texture_names.emplace_back(name);
      scene.textures.emplace_back();
      auto texture_id   = (int)scene.textures.size() - 1;
      texture_map[name] = texture_id;
      value             = texture_id;
    }
  };

  // load json instance
  struct ply_instance {
    vector<frame3f> frames = {};
  };
  using ply_instance_handle = int;
  auto ply_instances        = vector<ply_instance>{};
  auto ply_instances_names  = vector<string>{};
  auto ply_instance_map     = unordered_map<string, ply_instance_handle>{
      {"", invalidid}};
  auto instance_ply = unordered_map<int, ply_instance_handle>{};
  auto get_ist      = [&scene, &ply_instances, &ply_instances_names,
                     &ply_instance_map, &instance_ply](const json_value& json,
                     const string& key, const instance_data& instance) {
    auto name = json.value(key, string{});
    if (name.empty()) return;
    auto instance_id = (int)(&instance - scene.instances.data());
    auto it          = ply_instance_map.find(name);
    if (it != ply_instance_map.end()) {
      instance_ply[instance_id] = it->second;
    } else {
      ply_instances_names.emplace_back(name);
      ply_instances.emplace_back(ply_instance());
      auto ply_instance_id      = (int)ply_instances.size() - 1;
      ply_instance_map[name]    = ply_instance_id;
      instance_ply[instance_id] = ply_instance_id;
    }
  };
  auto get_ply_instance_name = [&ply_instances, &ply_instances_names](
                                   const scene_data&   scene,
                                   const ply_instance& instance) -> string {
    return ply_instances_names[&instance - ply_instances.data()];
  };

  // parsing values
  try {
    if (json.contains("asset")) {
      auto& element = json.at("asset");
      get_opt(element, "copyright", scene.copyright);
    }
    if (json.contains("cameras")) {
      for (auto& [key, element] : json.at("cameras").items()) {
        auto& camera = scene.cameras.emplace_back();
        scene.camera_names.emplace_back(key);
        get_of3(element, "frame", camera.frame);
        get_opt(element, "orthographic", camera.orthographic);
        get_opt(element, "ortho", camera.orthographic);
        get_opt(element, "lens", camera.lens);
        get_opt(element, "aspect", camera.aspect);
        get_opt(element, "film", camera.film);
        get_opt(element, "focus", camera.focus);
        get_opt(element, "aperture", camera.aperture);
        if (element.contains("lookat")) {
          get_om3(element, "lookat", (mat3f&)camera.frame);
          camera.focus = length(camera.frame.x - camera.frame.y);
          camera.frame = lookat_frame(
              camera.frame.x, camera.frame.y, camera.frame.z);
        }
      }
    }
    if (json.contains("environments")) {
      for (auto& [key, element] : json.at("environments").items()) {
        auto& environment = scene.environments.emplace_back();
        scene.environment_names.emplace_back(key);
        get_of3(element, "frame", environment.frame);
        get_opt(element, "emission", environment.emission);
        get_tex(element, "emission_tex", environment.emission_tex);
        if (element.contains("lookat")) {
          get_om3(element, "lookat", (mat3f&)environment.frame);
          environment.frame = lookat_frame(environment.frame.x,
              environment.frame.y, environment.frame.z, false);
        }
      }
    }
    if (json.contains("materials")) {
      for (auto& [key, element] : json.at("materials").items()) {
        auto& material = scene.materials.emplace_back();
        scene.material_names.emplace_back(key);
        material_map[key] = (int)scene.materials.size() - 1;
        auto type40       = material_type40::matte;
        get_opt(element, "type", type40);
        material.type = (material_type)type40;
        get_opt(element, "emission", material.emission);
        get_opt(element, "color", material.color);
        get_opt(element, "metallic", material.metallic);
        get_opt(element, "roughness", material.roughness);
        get_opt(element, "ior", material.ior);
        get_opt(element, "trdepth", material.trdepth);
        get_opt(element, "scattering", material.scattering);
        get_opt(element, "scanisotropy", material.scanisotropy);
        get_opt(element, "opacity", material.opacity);
        get_tex(element, "emission_tex", material.emission_tex);
        get_tex(element, "color_tex", material.color_tex);
        get_tex(element, "roughness_tex", material.roughness_tex);
        get_tex(element, "scattering_tex", material.scattering_tex);
        get_tex(element, "normal_tex", material.normal_tex);
      }
    }
    if (json.contains("instances")) {
      for (auto& [key, element] : json.at("instances").items()) {
        auto& instance = scene.instances.emplace_back();
        scene.instance_names.emplace_back(key);
        get_of3(element, "frame", instance.frame);
        get_shp(element, "shape", instance.shape);
        get_mat(element, "material", instance.material);
        if (element.contains("lookat")) {
          get_om3(element, "lookat", (mat3f&)instance.frame);
          instance.frame = lookat_frame(
              instance.frame.x, instance.frame.y, instance.frame.z, false);
        }
      }
    }
    if (json.contains("objects")) {
      for (auto& [key, element] : json.at("objects").items()) {
        auto& instance = scene.instances.emplace_back();
        scene.instance_names.emplace_back(key);
        get_of3(element, "frame", instance.frame);
        get_shp(element, "shape", instance.shape);
        get_mat(element, "material", instance.material);
        if (element.contains("lookat")) {
          get_om3(element, "lookat", (mat3f&)instance.frame);
          instance.frame = lookat_frame(
              instance.frame.x, instance.frame.y, instance.frame.z, false);
        }
        if (element.contains("instance")) {
          get_ist(element, "instance", instance);
        }
      }
    }
    if (json.contains("subdivs")) {
      for (auto& [key, element] : json.at("subdivs").items()) {
        auto& subdiv = scene.subdivs.emplace_back();
        scene.subdiv_names.emplace_back(key);
        get_shp(element, "shape", subdiv.shape);
        get_opt(element, "subdivisions", subdiv.subdivisions);
        get_opt(element, "catmullclark", subdiv.catmullclark);
        get_opt(element, "smooth", subdiv.smooth);
        get_opt(element, "displacement", subdiv.displacement);
        get_tex(element, "displacement_tex", subdiv.displacement_tex);
      }
    }
  } catch (...) {
    error = filename + ": parse error";
    return false;
  }

  // dirname
  auto dirname         = path_dirname(filename);
  auto dependent_error = [&filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // get filename from name
  auto find_path = [dirname](const string& name, const string& group,
                       const vector<string>& extensions) {
    for (auto& extension : extensions) {
      auto path = path_join(dirname, group, name + extension);
      if (path_exists(path)) return path_join(group, name + extension);
    }
    return path_join(group, name + extensions.front());
  };

  // load resources
  if (noparallel) {
    auto error = string{};
    // load shapes
    for (auto& shape : scene.shapes) {
      auto path = find_path(
          get_shape_name(scene, shape), "shapes", {".ply", ".obj"});
      if (!load_shape(path_join(dirname, path), shape, error, true))
        return dependent_error();
    }
    // load subdivs
    for (auto& subdiv : scene.subdivs) {
      auto path = find_path(
          get_subdiv_name(scene, subdiv), "subdivs", {".ply", ".obj"});
      if (!load_subdiv(path_join(dirname, path), subdiv, error))
        return dependent_error();
    }
    // load textures
    for (auto& texture : scene.textures) {
      auto path = find_path(get_texture_name(scene, texture), "textures",
          {".hdr", ".exr", ".png", ".jpg"});
      if (!load_texture(path_join(dirname, path), texture, error))
        return dependent_error();
    }
    // load instances
    for (auto& ply_instance : ply_instances) {
      auto path = find_path(
          get_ply_instance_name(scene, ply_instance), "instances", {".ply"});
      if (!load_instance(path_join(dirname, path), ply_instance.frames, error))
        return dependent_error();
    }
  } else {
    // load shapes
    if (!parallel_foreach(scene.shapes, error, [&](auto& shape, string& error) {
          auto path = find_path(
              get_shape_name(scene, shape), "shapes", {".ply", ".obj"});
          return load_shape(path_join(dirname, path), shape, error, true);
        }))
      return dependent_error();
    // load subdivs
    if (!parallel_foreach(
            scene.subdivs, error, [&](auto& subdiv, string& error) {
              auto path = find_path(
                  get_subdiv_name(scene, subdiv), "subdivs", {".ply", ".obj"});
              return load_subdiv(path_join(dirname, path), subdiv, error);
            }))
      return dependent_error();
    // load textures
    if (!parallel_foreach(
            scene.textures, error, [&](auto& texture, string& error) {
              auto path = find_path(get_texture_name(scene, texture),
                  "textures", {".hdr", ".exr", ".png", ".jpg"});
              return load_texture(path_join(dirname, path), texture, error);
            }))
      return dependent_error();
    // load instances
    if (!parallel_foreach(
            ply_instances, error, [&](auto& ply_instance, string& error) {
              auto path = find_path(get_ply_instance_name(scene, ply_instance),
                  "instances", {".ply"});
              return load_instance(
                  path_join(dirname, path), ply_instance.frames, error);
            }))
      return dependent_error();
  }

  // apply instances
  if (!ply_instances.empty()) {
    auto instances      = scene.instances;
    auto instance_names = scene.instance_names;
    scene.instances.clear();
    scene.instance_names.clear();
    for (auto& instance : instances) {
      auto it = instance_ply.find((int)(&instance - instances.data()));
      if (it == instance_ply.end()) {
        auto& ninstance = scene.instances.emplace_back();
        scene.instance_names.emplace_back(
            instance_names[&instance - instances.data()]);
        ninstance.frame    = instance.frame;
        ninstance.shape    = instance.shape;
        ninstance.material = instance.material;
      } else {
        auto& ply_instance = ply_instances[it->second];
        auto  instance_id  = 0;
        for (auto& frame : ply_instance.frames) {
          auto& ninstance = scene.instances.emplace_back();
          scene.instance_names.emplace_back(
              instance_names[&instance - instances.data()] + "_" +
              std::to_string(instance_id++));
          ninstance.frame    = frame * instance.frame;
          ninstance.shape    = instance.shape;
          ninstance.material = instance.material;
        }
      }
    }
  }

  // fix scene
  add_missing_camera(scene);
  add_missing_radius(scene);
  trim_memory(scene);

  // done
  return true;
}

// Load a scene in the builtin JSON format.
static bool load_json_scene_version41(const string& filename, json_value& json,
    scene_data& scene, string& error, bool noparallel) {
  // check version
  if (!json.contains("asset") || !json.at("asset").contains("version"))
    return load_json_scene_version40(filename, json, scene, error, noparallel);

  // parse json value
  auto get_opt = [](const json_value& json, const string& key, auto& value) {
    value = json.value(key, value);
  };
  auto get_ref = [](const json_value& json, const string& key, int& value,
                     const unordered_map<string, int>& map) {
    auto values = json.value(key, string{});
    value       = values.empty() ? -1 : map.at(values);
  };

  // references
  auto shape_map    = unordered_map<string, int>{};
  auto texture_map  = unordered_map<string, int>{};
  auto material_map = unordered_map<string, int>{};

  // filenames
  auto shape_filenames   = vector<string>{};
  auto texture_filenames = vector<string>{};
  auto subdiv_filenames  = vector<string>{};

  // parsing values
  try {
    if (json.contains("asset")) {
      auto& element = json.at("asset");
      get_opt(element, "copyright", scene.copyright);
    }
    if (json.contains("cameras")) {
      auto& group = json.at("cameras");
      scene.cameras.reserve(group.size());
      scene.camera_names.reserve(group.size());
      for (auto& [key, element] : group.items()) {
        auto& camera = scene.cameras.emplace_back();
        scene.camera_names.push_back(key);
        get_opt(element, "frame", camera.frame);
        get_opt(element, "orthographic", camera.orthographic);
        get_opt(element, "ortho", camera.orthographic);
        get_opt(element, "lens", camera.lens);
        get_opt(element, "aspect", camera.aspect);
        get_opt(element, "film", camera.film);
        get_opt(element, "focus", camera.focus);
        get_opt(element, "aperture", camera.aperture);
        if (element.contains("lookat")) {
          get_opt(element, "lookat", (mat3f&)camera.frame);
          camera.focus = length(camera.frame.x - camera.frame.y);
          camera.frame = lookat_frame(
              camera.frame.x, camera.frame.y, camera.frame.z);
        }
      }
    }
    if (json.contains("textures")) {
      auto& group = json.at("textures");
      scene.textures.reserve(group.size());
      scene.texture_names.reserve(group.size());
      texture_filenames.reserve(group.size());
      for (auto& [key, element] : group.items()) {
        [[maybe_unused]] auto& texture = scene.textures.emplace_back();
        scene.texture_names.push_back(key);
        auto& datafile   = texture_filenames.emplace_back();
        texture_map[key] = (int)scene.textures.size() - 1;
        if (element.is_string()) {
          auto filename       = element.get<string>();
          element             = json_value::object();
          element["datafile"] = filename;
        }
        get_opt(element, "datafile", datafile);
      }
    }
    if (json.contains("materials")) {
      auto& group = json.at("materials");
      scene.materials.reserve(group.size());
      scene.material_names.reserve(group.size());
      for (auto& [key, element] : json.at("materials").items()) {
        auto& material = scene.materials.emplace_back();
        scene.material_names.push_back(key);
        material_map[key] = (int)scene.materials.size() - 1;
        get_opt(element, "type", material.type);
        get_opt(element, "emission", material.emission);
        get_opt(element, "color", material.color);
        get_opt(element, "metallic", material.metallic);
        get_opt(element, "roughness", material.roughness);
        get_opt(element, "ior", material.ior);
        get_opt(element, "trdepth", material.trdepth);
        get_opt(element, "scattering", material.scattering);
        get_opt(element, "scanisotropy", material.scanisotropy);
        get_opt(element, "opacity", material.opacity);
        get_ref(element, "emission_tex", material.emission_tex, texture_map);
        get_ref(element, "color_tex", material.color_tex, texture_map);
        get_ref(element, "roughness_tex", material.roughness_tex, texture_map);
        get_ref(
            element, "scattering_tex", material.scattering_tex, texture_map);
        get_ref(element, "normal_tex", material.normal_tex, texture_map);
      }
    }
    if (json.contains("shapes")) {
      auto& group = json.at("shapes");
      scene.shapes.reserve(group.size());
      scene.shape_names.reserve(group.size());
      shape_filenames.reserve(group.size());
      for (auto& [key, element] : group.items()) {
        [[maybe_unused]] auto& shape = scene.shapes.emplace_back();
        scene.shape_names.push_back(key);
        auto& datafile = shape_filenames.emplace_back();
        shape_map[key] = (int)scene.shapes.size() - 1;
        if (element.is_string()) {
          auto filename       = element.get<string>();
          element             = json_value::object();
          element["datafile"] = filename;
        }
        get_opt(element, "datafile", datafile);
      }
    }
    if (json.contains("subdivs")) {
      auto& group = json.at("subdivs");
      scene.subdivs.reserve(group.size());
      scene.subdiv_names.reserve(group.size());
      subdiv_filenames.reserve(group.size());
      for (auto& [key, element] : group.items()) {
        auto& subdiv = scene.subdivs.emplace_back();
        scene.subdiv_names.emplace_back(key);
        auto& datafile = subdiv_filenames.emplace_back();
        get_opt(element, "datafile", datafile);
        get_ref(element, "shape", subdiv.shape, shape_map);
        get_opt(element, "subdivisions", subdiv.subdivisions);
        get_opt(element, "catmullclark", subdiv.catmullclark);
        get_opt(element, "smooth", subdiv.smooth);
        get_opt(element, "displacement", subdiv.displacement);
        get_ref(
            element, "displacement_tex", subdiv.displacement_tex, texture_map);
      }
    }
    if (json.contains("instances")) {
      auto& group = json.at("instances");
      scene.instances.reserve(group.size());
      scene.instance_names.reserve(group.size());
      for (auto& [key, element] : group.items()) {
        auto& instance = scene.instances.emplace_back();
        scene.instance_names.emplace_back(key);
        get_opt(element, "frame", instance.frame);
        get_ref(element, "shape", instance.shape, shape_map);
        get_ref(element, "material", instance.material, material_map);
        if (element.contains("lookat")) {
          get_opt(element, "lookat", (mat3f&)instance.frame);
          instance.frame = lookat_frame(
              instance.frame.x, instance.frame.y, instance.frame.z, false);
        }
      }
    }
    if (json.contains("environments")) {
      auto& group = json.at("environments");
      scene.instances.reserve(group.size());
      scene.instance_names.reserve(group.size());
      for (auto& [key, element] : group.items()) {
        auto& environment = scene.environments.emplace_back();
        scene.environment_names.push_back(key);
        get_opt(element, "frame", environment.frame);
        get_opt(element, "emission", environment.emission);
        get_ref(element, "emission_tex", environment.emission_tex, texture_map);
        if (element.contains("lookat")) {
          get_opt(element, "lookat", (mat3f&)environment.frame);
          environment.frame = lookat_frame(environment.frame.x,
              environment.frame.y, environment.frame.z, false);
        }
      }
    }
  } catch (...) {
    error = filename + ": parse error";
    return false;
  }

  // prepare data
  auto dirname         = path_dirname(filename);
  auto dependent_error = [&filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // fix paths
  for (auto& datafile : shape_filenames)
    datafile = path_join(dirname, "shapes", datafile);
  for (auto& datafile : texture_filenames)
    datafile = path_join(dirname, "textures", datafile);
  for (auto& datafile : subdiv_filenames)
    datafile = path_join(dirname, "subdivs", datafile);

  // load resources
  if (noparallel) {
    auto error = string{};
    // load shapes
    for (auto idx : range(scene.shapes.size())) {
      if (!load_shape(shape_filenames[idx], scene.shapes[idx], error, true))
        return dependent_error();
    }
    // load subdivs
    for (auto idx : range(scene.subdivs.size())) {
      if (!load_subdiv(subdiv_filenames[idx], scene.subdivs[idx], error))
        return dependent_error();
    }
    // load textures
    for (auto idx : range(scene.textures.size())) {
      if (!load_texture(texture_filenames[idx], scene.textures[idx], error))
        return dependent_error();
    }
  } else {
    // load shapes
    if (!parallel_for(
            scene.shapes.size(), error, [&](size_t idx, string& error) {
              return load_shape(
                  shape_filenames[idx], scene.shapes[idx], error, true);
            }))
      return dependent_error();
    // load subdivs
    if (!parallel_for(
            scene.subdivs.size(), error, [&](size_t idx, string& error) {
              return load_subdiv(
                  subdiv_filenames[idx], scene.subdivs[idx], error);
            }))
      return dependent_error();
    // load textures
    if (!parallel_for(
            scene.textures.size(), error, [&](size_t idx, string& error) {
              return load_texture(
                  texture_filenames[idx], scene.textures[idx], error);
            }))
      return dependent_error();
  }

  // fix scene
  add_missing_camera(scene);
  add_missing_radius(scene);
  trim_memory(scene);

  // done
  return false;
}

// Load a scene in the builtin JSON format.
static bool load_json_scene(
    const string& filename, scene_data& scene, string& error, bool noparallel) {
  // open file
  auto json = json_value{};
  if (!load_json(filename, json, error)) return false;

  // check version
  if (!json.contains("asset") || !json.at("asset").contains("version"))
    return load_json_scene_version40(filename, json, scene, error, noparallel);
  if (json.contains("asset") && json.at("asset").contains("version") &&
      json.at("asset").at("version") == "4.1")
    return load_json_scene_version41(filename, json, scene, error, noparallel);

  // parse json value
  auto get_opt = [](const json_value& json, const string& key, auto& value) {
    value = json.value(key, value);
  };

  // filenames
  auto shape_filenames   = vector<string>{};
  auto texture_filenames = vector<string>{};
  auto subdiv_filenames  = vector<string>{};

  // errors
  auto parse_error = [&filename, &error]() {
    error = filename + ": parse error";
    return false;
  };

  // parsing values
  try {
    if (json.contains("asset")) {
      auto& element = json.at("asset");
      get_opt(element, "copyright", scene.copyright);
      auto version = string{};
      get_opt(element, "version", version);
      if (version != "4.2" && version != "5.0") return parse_error();
    }
    if (json.contains("cameras")) {
      auto& group = json.at("cameras");
      scene.cameras.reserve(group.size());
      scene.camera_names.reserve(group.size());
      for (auto& element : group) {
        auto& camera = scene.cameras.emplace_back();
        auto& name   = scene.camera_names.emplace_back();
        get_opt(element, "name", name);
        get_opt(element, "frame", camera.frame);
        get_opt(element, "orthographic", camera.orthographic);
        get_opt(element, "lens", camera.lens);
        get_opt(element, "aspect", camera.aspect);
        get_opt(element, "film", camera.film);
        get_opt(element, "focus", camera.focus);
        get_opt(element, "aperture", camera.aperture);
      }
    }
    if (json.contains("textures")) {
      auto& group = json.at("textures");
      scene.textures.reserve(group.size());
      scene.texture_names.reserve(group.size());
      texture_filenames.reserve(group.size());
      for (auto& element : group) {
        [[maybe_unused]] auto& texture = scene.textures.emplace_back();
        auto&                  name    = scene.texture_names.emplace_back();
        auto&                  uri     = texture_filenames.emplace_back();
        get_opt(element, "name", name);
        get_opt(element, "uri", uri);
      }
    }
    if (json.contains("materials")) {
      auto& group = json.at("materials");
      scene.materials.reserve(group.size());
      scene.material_names.reserve(group.size());
      for (auto& element : json.at("materials")) {
        auto& material = scene.materials.emplace_back();
        auto& name     = scene.material_names.emplace_back();
        get_opt(element, "name", name);
        get_opt(element, "type", material.type);
        get_opt(element, "emission", material.emission);
        get_opt(element, "color", material.color);
        get_opt(element, "metallic", material.metallic);
        get_opt(element, "roughness", material.roughness);
        get_opt(element, "ior", material.ior);
        get_opt(element, "trdepth", material.trdepth);
        get_opt(element, "scattering", material.scattering);
        get_opt(element, "scanisotropy", material.scanisotropy);
        get_opt(element, "opacity", material.opacity);
        get_opt(element, "emission_tex", material.emission_tex);
        get_opt(element, "color_tex", material.color_tex);
        get_opt(element, "roughness_tex", material.roughness_tex);
        get_opt(element, "scattering_tex", material.scattering_tex);
        get_opt(element, "normal_tex", material.normal_tex);
      }
    }
    if (json.contains("shapes")) {
      auto& group = json.at("shapes");
      scene.shapes.reserve(group.size());
      scene.shape_names.reserve(group.size());
      shape_filenames.reserve(group.size());
      for (auto& element : group) {
        [[maybe_unused]] auto& shape = scene.shapes.emplace_back();
        auto&                  name  = scene.shape_names.emplace_back();
        auto&                  uri   = shape_filenames.emplace_back();
        get_opt(element, "name", name);
        get_opt(element, "uri", uri);
      }
    }
    if (json.contains("subdivs")) {
      auto& group = json.at("subdivs");
      scene.subdivs.reserve(group.size());
      scene.subdiv_names.reserve(group.size());
      subdiv_filenames.reserve(group.size());
      for (auto& element : group) {
        auto& subdiv = scene.subdivs.emplace_back();
        auto& name   = scene.subdiv_names.emplace_back();
        auto& uri    = subdiv_filenames.emplace_back();
        get_opt(element, "name", name);
        get_opt(element, "uri", uri);
        get_opt(element, "shape", subdiv.shape);
        get_opt(element, "subdivisions", subdiv.subdivisions);
        get_opt(element, "catmullclark", subdiv.catmullclark);
        get_opt(element, "smooth", subdiv.smooth);
        get_opt(element, "displacement", subdiv.displacement);
        get_opt(element, "displacement_tex", subdiv.displacement_tex);
      }
    }
    if (json.contains("instances")) {
      auto& group = json.at("instances");
      scene.instances.reserve(group.size());
      scene.instance_names.reserve(group.size());
      for (auto& element : group) {
        auto& instance = scene.instances.emplace_back();
        auto& name     = scene.instance_names.emplace_back();
        get_opt(element, "name", name);
        get_opt(element, "frame", instance.frame);
        get_opt(element, "shape", instance.shape);
        get_opt(element, "material", instance.material);
      }
    }
    if (json.contains("environments")) {
      auto& group = json.at("environments");
      scene.instances.reserve(group.size());
      scene.instance_names.reserve(group.size());
      for (auto& element : group) {
        auto& environment = scene.environments.emplace_back();
        auto& name        = scene.environment_names.emplace_back();
        get_opt(element, "name", name);
        get_opt(element, "frame", environment.frame);
        get_opt(element, "emission", environment.emission);
        get_opt(element, "emission_tex", environment.emission_tex);
      }
    }
  } catch (...) {
    return parse_error();
  }

  // prepare data
  auto dirname         = path_dirname(filename);
  auto dependent_error = [&filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // load resources
  if (noparallel) {
    // load shapes
    for (auto idx : range(scene.shapes.size())) {
      if (!load_shape(path_join(dirname, shape_filenames[idx]),
              scene.shapes[idx], error, true))
        return dependent_error();
    }
    // load subdivs
    for (auto idx : range(scene.subdivs.size())) {
      if (!load_subdiv(path_join(dirname, subdiv_filenames[idx]),
              scene.subdivs[idx], error))
        return dependent_error();
    }
    // load textures
    for (auto idx : range(scene.textures.size())) {
      if (!load_texture(path_join(dirname, texture_filenames[idx]),
              scene.textures[idx], error))
        return dependent_error();
    }
  } else {
    // load shapes
    if (!parallel_for(
            scene.shapes.size(), error, [&](size_t idx, string& error) {
              return load_shape(path_join(dirname, shape_filenames[idx]),
                  scene.shapes[idx], error, true);
            }))
      return dependent_error();
    // load subdivs
    if (!parallel_for(
            scene.subdivs.size(), error, [&](size_t idx, string& error) {
              return load_subdiv(path_join(dirname, subdiv_filenames[idx]),
                  scene.subdivs[idx], error);
            }))
      return dependent_error();
    // load textures
    if (!parallel_for(
            scene.textures.size(), error, [&](size_t idx, string& error) {
              return load_texture(path_join(dirname, texture_filenames[idx]),
                  scene.textures[idx], error);
            }))
      return dependent_error();
  }

  // fix scene
  add_missing_camera(scene);
  add_missing_radius(scene);
  trim_memory(scene);

  // done
  return true;
}

// Save a scene in the builtin JSON format.
static bool save_json_scene(const string& filename, const scene_data& scene,
    string& error, bool noparallel) {
  // helpers to handel old code paths
  auto add_object = [](json_value& json, const string& name) -> json_value& {
    auto& item = json[name];
    item       = json_value::object();
    return item;
  };
  auto add_array = [](json_value& json, const string& name) -> json_value& {
    auto& item = json[name];
    item       = json_value::array();
    return item;
  };
  auto append_object = [](json_value& json) -> json_value& {
    auto& item = json.emplace_back();
    item       = json_value::object();
    return item;
  };
  auto set_val = [](json_value& json, const string& name, const auto& value,
                     const auto& def) {
    if (value == def) return;
    json[name] = value;
  };
  auto set_ref = [](json_value& json, const string& name, int value) {
    if (value < 0) return;
    json[name] = value;
  };
  auto reserve_values = [](json_value& json, size_t size) {
    json.get_ptr<json_value::array_t*>()->reserve(size);
  };

  // names
  auto get_name = [](const vector<string>& names, size_t idx) -> string {
    return (idx < names.size()) ? names[idx] : "";
  };
  auto get_filename = [](const vector<string>& names, size_t idx,
                          const string& basename,
                          const string& extension) -> string {
    if (idx < names.size()) {
      return basename + "s/" + names[idx] + extension;
    } else {
      return basename + "s/" + basename + std::to_string(idx) + extension;
    }
  };

  // filenames
  auto shape_filenames   = vector<string>(scene.shapes.size());
  auto texture_filenames = vector<string>(scene.textures.size());
  auto subdiv_filenames  = vector<string>(scene.subdivs.size());
  for (auto idx : range(shape_filenames.size())) {
    shape_filenames[idx] = get_filename(
        scene.shape_names, idx, "shape", ".ply");
  }
  for (auto idx : range(texture_filenames.size())) {
    texture_filenames[idx] = get_filename(scene.texture_names, idx, "texture",
        (scene.textures[idx].pixelsf.empty() ? ".png" : ".hdr"));
  }
  for (auto idx : range(subdiv_filenames.size())) {
    subdiv_filenames[idx] = get_filename(
        scene.subdiv_names, idx, "subdiv", ".obj");
  }

  // save json file
  auto json = json_value::object();

  // asset
  {
    auto& element = add_object(json, "asset");
    set_val(element, "copyright", scene.copyright, "");
    set_val(element, "generator",
        "Yocto/GL - https://github.com/xelatihy/yocto-gl"s, ""s);
    set_val(element, "version", "4.2"s, ""s);
  }

  if (!scene.cameras.empty()) {
    auto  default_ = camera_data{};
    auto& group    = add_array(json, "cameras");
    reserve_values(group, scene.cameras.size());
    for (auto&& [idx, camera] : enumerate(scene.cameras)) {
      auto& element = append_object(group);
      set_val(element, "name", get_name(scene.camera_names, idx), "");
      set_val(element, "frame", camera.frame, default_.frame);
      set_val(
          element, "orthographic", camera.orthographic, default_.orthographic);
      set_val(element, "lens", camera.lens, default_.lens);
      set_val(element, "aspect", camera.aspect, default_.aspect);
      set_val(element, "film", camera.film, default_.film);
      set_val(element, "focus", camera.focus, default_.focus);
      set_val(element, "aperture", camera.aperture, default_.aperture);
    }
  }

  if (!scene.textures.empty()) {
    auto& group = add_array(json, "textures");
    reserve_values(group, scene.textures.size());
    for (auto&& [idx, texture] : enumerate(scene.textures)) {
      auto& element = append_object(group);
      set_val(element, "name", get_name(scene.texture_names, idx), "");
      set_val(element, "uri", texture_filenames[idx], ""s);
    }
  }

  if (!scene.materials.empty()) {
    auto  default_ = material_data{};
    auto& group    = add_array(json, "materials");
    reserve_values(group, scene.materials.size());
    for (auto&& [idx, material] : enumerate(scene.materials)) {
      auto& element = append_object(group);
      set_val(element, "name", get_name(scene.material_names, idx), "");
      set_val(element, "type", material.type, default_.type);
      set_val(element, "emission", material.emission, default_.emission);
      set_val(element, "color", material.color, default_.color);
      set_val(element, "metallic", material.metallic, default_.metallic);
      set_val(element, "roughness", material.roughness, default_.roughness);
      set_val(element, "ior", material.ior, default_.ior);
      set_val(element, "trdepth", material.trdepth, default_.trdepth);
      set_val(element, "scattering", material.scattering, default_.scattering);
      set_val(element, "scanisotropy", material.scanisotropy,
          default_.scanisotropy);
      set_val(element, "opacity", material.opacity, default_.opacity);
      set_val(element, "emission_tex", material.emission_tex,
          default_.emission_tex);
      set_val(element, "color_tex", material.color_tex, default_.color_tex);
      set_val(element, "roughness_tex", material.roughness_tex,
          default_.roughness_tex);
      set_val(element, "scattering_tex", material.scattering_tex,
          default_.scattering_tex);
      set_val(element, "normal_tex", material.normal_tex, default_.normal_tex);
    }
  }

  if (!scene.shapes.empty()) {
    auto& group = add_array(json, "shapes");
    reserve_values(group, scene.shapes.size());
    for (auto&& [idx, shape] : enumerate(scene.shapes)) {
      auto& element = append_object(group);
      set_val(element, "name", get_name(scene.shape_names, idx), "");
      set_val(element, "uri", shape_filenames[idx], "");
    }
  }

  if (!scene.subdivs.empty()) {
    auto  default_ = subdiv_data{};
    auto& group    = add_array(json, "subdivs");
    reserve_values(group, scene.subdivs.size());
    for (auto&& [idx, subdiv] : enumerate(scene.subdivs)) {
      auto& element = append_object(group);
      set_val(element, "name", get_name(scene.subdiv_names, idx), "");
      set_ref(element, "shape", subdiv.shape);
      set_val(element, "uri", subdiv_filenames[idx], "");
      set_val(
          element, "subdivisions", subdiv.subdivisions, default_.subdivisions);
      set_val(
          element, "catmullclark", subdiv.catmullclark, default_.subdivisions);
      set_val(element, "smooth", subdiv.smooth, default_.subdivisions);
      set_val(
          element, "displacement", subdiv.displacement, default_.subdivisions);
      set_val(element, "displacement_tex",
          get_texture_name(scene, subdiv.displacement_tex), "");
    }
  }

  if (!scene.instances.empty()) {
    auto  default_ = instance_data{};
    auto& group    = add_array(json, "instances");
    reserve_values(group, scene.instances.size());
    for (auto&& [idx, instance] : enumerate(scene.instances)) {
      auto& element = append_object(group);
      set_val(element, "name", get_name(scene.instance_names, idx), "");
      set_val(element, "frame", instance.frame, default_.frame);
      set_val(element, "shape", instance.shape, default_.shape);
      set_val(element, "material", instance.material, default_.material);
    }
  }

  if (!scene.environments.empty()) {
    auto  default_ = environment_data{};
    auto& group    = add_array(json, "environments");
    reserve_values(group, scene.environments.size());
    for (auto&& [idx, environment] : enumerate(scene.environments)) {
      auto& element = append_object(group);
      set_val(element, "name", get_name(scene.environment_names, idx), "");
      set_val(element, "frame", environment.frame, default_.frame);
      set_val(element, "emission", environment.emission, default_.emission);
      set_val(element, "emission_tex", environment.emission_tex,
          default_.emission_tex);
    }
  }

  // save json
  if (!save_json(filename, json, error)) return false;

  // prepare data
  auto dirname         = path_dirname(filename);
  auto dependent_error = [&filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // dirname
  if (noparallel) {
    // save shapes
    for (auto idx : range(scene.shapes.size())) {
      if (!save_shape(path_join(dirname, shape_filenames[idx]),
              scene.shapes[idx], error, true))
        return dependent_error();
    }
    // save subdiv
    for (auto idx : range(scene.subdivs.size())) {
      if (!save_subdiv(path_join(dirname, subdiv_filenames[idx]),
              scene.subdivs[idx], error))
        return dependent_error();
    }
    // save textures
    for (auto idx : range(scene.textures.size())) {
      if (!save_texture(path_join(dirname, texture_filenames[idx]),
              scene.textures[idx], error))
        return dependent_error();
    }
  } else {
    // save shapes
    if (!parallel_for(scene.shapes.size(), error, [&](auto idx, string& error) {
          return save_shape(path_join(dirname, shape_filenames[idx]),
              scene.shapes[idx], error, true);
        }))
      return dependent_error();
    // save subdivs
    if (!parallel_for(
            scene.subdivs.size(), error, [&](auto idx, string& error) {
              return save_subdiv(path_join(dirname, subdiv_filenames[idx]),
                  scene.subdivs[idx], error);
            }))
      return dependent_error();
    // save textures
    if (!parallel_for(
            scene.textures.size(), error, [&](auto idx, string& error) {
              return save_texture(path_join(dirname, texture_filenames[idx]),
                  scene.textures[idx], error);
            }))
      return dependent_error();
  }

  // done
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// OBJ CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

// Loads an OBJ
static bool load_obj_scene(
    const string& filename, scene_data& scene, string& error, bool noparallel) {
  // load obj
  auto obj = obj_model{};
  if (!load_obj(filename, obj, error, false, true)) return false;

  // convert cameras
  scene.cameras.reserve(obj.cameras.size());
  for (auto& ocamera : obj.cameras) {
    auto& camera        = scene.cameras.emplace_back();
    camera.frame        = ocamera.frame;
    camera.orthographic = ocamera.ortho;
    camera.film         = ocamera.film;
    camera.aspect       = ocamera.aspect;
    camera.focus        = ocamera.focus;
    camera.lens         = ocamera.lens;
    camera.aperture     = ocamera.aperture;
  }

  // convert between roughness and exponent
  auto exponent_to_roughness = [](float exponent) {
    if (exponent >= 1000) return 0.0f;
    auto roughness = exponent;
    roughness      = pow(2 / (roughness + 2), 1 / 4.0f);
    if (roughness < 0.01f) roughness = 0;
    if (roughness > 0.99f) roughness = 1;
    return roughness;
  };

  // handler for textures
  auto texture_paths = vector<string>{};
  for (auto& otexture : obj.textures) {
    scene.textures.emplace_back();
    texture_paths.emplace_back(otexture.path);
  }

  // handler for materials
  scene.materials.reserve(obj.materials.size());
  for (auto& omaterial : obj.materials) {
    auto& material        = scene.materials.emplace_back();
    material.type         = material_type::gltfpbr;
    material.emission     = omaterial.emission;
    material.emission_tex = omaterial.emission_tex;
    if (max(omaterial.transmission) > 0.1) {
      material.type      = material_type::transparent;
      material.color     = omaterial.transmission;
      material.color_tex = omaterial.transmission_tex;
    } else if (max(omaterial.specular) > 0.2) {
      material.type      = material_type::reflective;
      material.color     = omaterial.specular;
      material.color_tex = omaterial.specular_tex;
    } else if (max(omaterial.specular) > 0) {
      material.type      = material_type::glossy;
      material.color     = omaterial.diffuse;
      material.color_tex = omaterial.diffuse_tex;
    } else {
      material.type      = material_type::matte;
      material.color     = omaterial.diffuse;
      material.color_tex = omaterial.diffuse_tex;
    }
    material.roughness  = exponent_to_roughness(omaterial.exponent);
    material.ior        = omaterial.ior;
    material.metallic   = 0;
    material.opacity    = omaterial.opacity;
    material.normal_tex = omaterial.normal_tex;
  }

  // convert shapes
  scene.shapes.reserve(obj.shapes.size());
  scene.instances.reserve(obj.shapes.size());
  for (auto& oshape : obj.shapes) {
    if (oshape.elements.empty()) continue;
    auto& shape       = scene.shapes.emplace_back();
    auto& instance    = scene.instances.emplace_back();
    instance.shape    = (int)scene.shapes.size() - 1;
    instance.material = oshape.elements.front().material;
    get_positions(oshape, shape.positions);
    get_normals(oshape, shape.normals);
    get_texcoords(oshape, shape.texcoords, true);
    get_faces(oshape, instance.material, shape.triangles, shape.quads);
    get_lines(oshape, instance.material, shape.lines);
    get_points(oshape, instance.material, shape.points);
  }

  // convert environments
  scene.environments.reserve(obj.environments.size());
  for (auto& oenvironment : obj.environments) {
    auto& environment        = scene.environments.emplace_back();
    environment.frame        = oenvironment.frame;
    environment.emission     = oenvironment.emission;
    environment.emission_tex = oenvironment.emission_tex;
  }

  // names
  scene.camera_names   = make_names(scene.cameras, {}, "camera");
  scene.texture_names  = make_names(scene.textures, {}, "texture");
  scene.material_names = make_names(scene.materials, {}, "material");
  scene.shape_names    = make_names(scene.shapes, {}, "shape");
  scene.subdiv_names   = make_names(scene.subdivs, {}, "subdiv");
  scene.instance_names = make_names(scene.instances, {}, "instance");

  // dirname
  auto dirname         = path_dirname(filename);
  auto dependent_error = [&filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  if (noparallel) {
    // load textures
    for (auto& texture : scene.textures) {
      auto& path = texture_paths[&texture - &scene.textures.front()];
      if (!load_texture(path_join(dirname, path), texture, error))
        return dependent_error();
    }
  } else {
    // load textures
    if (!parallel_foreach(
            scene.textures, error, [&](auto& texture, string& error) {
              auto& path = texture_paths[&texture - &scene.textures.front()];
              return load_texture(path_join(dirname, path), texture, error);
            }))
      return dependent_error();
  }

  // fix scene
  add_missing_camera(scene);
  add_missing_radius(scene);

  // done
  return true;
}

static bool save_obj_scene(const string& filename, const scene_data& scene,
    string& error, bool noparallel) {
  // build obj
  auto obj = obj_model{};

  // convert cameras
  for (auto& camera : scene.cameras) {
    auto& ocamera    = obj.cameras.emplace_back();
    ocamera.name     = get_camera_name(scene, camera);
    ocamera.frame    = camera.frame;
    ocamera.ortho    = camera.orthographic;
    ocamera.film     = camera.film;
    ocamera.aspect   = camera.aspect;
    ocamera.focus    = camera.focus;
    ocamera.lens     = camera.lens;
    ocamera.aperture = camera.aperture;
  }

  // helper
  auto roughness_to_exponent = [](float roughness) -> float {
    if (roughness < 0.01f) return 10000;
    if (roughness > 0.99f) return 10;
    return 2 / pow(roughness, 4.0f) - 2;
  };

  // convert textures
  for (auto& texture : scene.textures) {
    auto& otexture = obj.textures.emplace_back();
    otexture.path  = "textures/" + get_texture_name(scene, texture) +
                    (!texture.pixelsf.empty() ? ".hdr"s : ".png"s);
  }

  // convert materials
  for (auto& material : scene.materials) {
    auto& omaterial        = obj.materials.emplace_back();
    omaterial.name         = get_material_name(scene, material);
    omaterial.illum        = 2;
    omaterial.emission     = material.emission;
    omaterial.diffuse      = material.color;
    omaterial.specular     = {0, 0, 0};
    omaterial.exponent     = roughness_to_exponent(material.roughness);
    omaterial.opacity      = material.opacity;
    omaterial.emission_tex = material.emission_tex;
    omaterial.diffuse_tex  = material.color_tex;
    omaterial.normal_tex   = material.normal_tex;
  }

  // convert objects
  for (auto& instance : scene.instances) {
    auto& shape     = scene.shapes[instance.shape];
    auto  positions = shape.positions, normals = shape.normals;
    for (auto& p : positions) p = transform_point(instance.frame, p);
    for (auto& n : normals) n = transform_normal(instance.frame, n);
    auto& oshape = obj.shapes.emplace_back();
    oshape.name  = get_shape_name(scene, shape);
    add_positions(oshape, positions);
    add_normals(oshape, normals);
    add_texcoords(oshape, shape.texcoords, true);
    add_triangles(oshape, shape.triangles, instance.material,
        !shape.normals.empty(), !shape.texcoords.empty());
    add_quads(oshape, shape.quads, instance.material, !shape.normals.empty(),
        !shape.texcoords.empty());
    add_lines(oshape, shape.lines, instance.material, !shape.normals.empty(),
        !shape.texcoords.empty());
    add_points(oshape, shape.points, instance.material, !shape.normals.empty(),
        !shape.texcoords.empty());
  }

  // convert environments
  for (auto& environment : scene.environments) {
    auto& oenvironment        = obj.environments.emplace_back();
    oenvironment.name         = get_environment_name(scene, environment);
    oenvironment.frame        = environment.frame;
    oenvironment.emission     = environment.emission;
    oenvironment.emission_tex = environment.emission_tex;
  }

  // save obj
  if (!save_obj(filename, obj, error)) return false;

  // dirname
  auto dirname         = path_dirname(filename);
  auto dependent_error = [&filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  if (noparallel) {
    // save textures
    for (auto& texture : scene.textures) {
      auto path = "textures/" + get_texture_name(scene, texture) +
                  (!texture.pixelsf.empty() ? ".hdr"s : ".png"s);
      if (!save_texture(path_join(dirname, path), texture, error))
        return dependent_error();
    }
  } else {
    // save textures
    if (!parallel_foreach(
            scene.textures, error, [&](auto& texture, string& error) {
              auto path = "textures/" + get_texture_name(scene, texture) +
                          (!texture.pixelsf.empty() ? ".hdr"s : ".png"s);
              return save_texture(path_join(dirname, path), texture, error);
            }))
      return dependent_error();
  }

  // done
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PLY CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static bool load_ply_scene(
    const string& filename, scene_data& scene, string& error, bool noparallel) {
  // load ply mesh and make instance
  auto shape = shape_data{};
  if (!load_shape(filename, shape, error, true)) return false;
  scene.shapes.push_back(shape);
  scene.instances.push_back({identity3x4f, (int)scene.shapes.size() - 1, -1});

  // fix scene
  add_missing_material(scene);
  add_missing_camera(scene);
  add_missing_radius(scene);

  // done
  return true;
}

static bool save_ply_scene(const string& filename, const scene_data& scene,
    string& error, bool noparallel) {
  // save shape
  if (scene.shapes.empty()) throw std::invalid_argument{"empty shape"};
  return save_shape(filename, scene.shapes.front(), error, false);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// STL CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static bool load_stl_scene(
    const string& filename, scene_data& scene, string& error, bool noparallel) {
  // load ply mesh and make instance
  auto shape = shape_data{};
  if (!load_shape(filename, shape, error, true)) return false;
  scene.instances.push_back({identity3x4f, (int)scene.shapes.size() - 1, -1});

  // fix scene
  add_missing_material(scene);
  add_missing_camera(scene);
  add_missing_radius(scene);

  // done
  return true;
}

static bool save_stl_scene(const string& filename, const scene_data& scene,
    string& error, bool noparallel) {
  // save shape
  if (scene.shapes.empty()) throw std::invalid_argument{"empty shape"};
  return save_shape(filename, scene.shapes.front(), error, false);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GLTF CONVESION
// -----------------------------------------------------------------------------
namespace yocto {

// Load a scene
static bool load_gltf_scene(
    const string& filename, scene_data& scene, string& error, bool noparallel) {
  // load gltf
  auto gltf = json_value{};
  if (!load_json(filename, gltf, error)) return false;

  // errors
  auto parse_error = [&filename, &error]() {
    error = filename + ": parse error";
    return false;
  };

  // parse buffers
  auto buffers_paths = vector<string>{};
  auto buffers       = vector<vector<byte>>();
  try {
    if (gltf.contains("buffers")) {
      for (auto& gbuffer : gltf.at("buffers")) {
        if (!gbuffer.contains("uri")) return parse_error();
        buffers_paths.push_back(gbuffer.value("uri", ""));
        buffers.emplace_back();
      }
    }
  } catch (...) {
    return parse_error();
  }

  // dirname
  auto dirname         = path_dirname(filename);
  auto dependent_error = [&filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  if (noparallel) {
    // load buffers
    for (auto& buffer : buffers) {
      auto& path = buffers_paths[&buffer - &buffers.front()];
      if (!load_binary(path_join(dirname, path), buffer, error))
        return dependent_error();
    }
  } else {
    // load buffers
    if (!parallel_foreach(buffers, error, [&](auto& buffer, string& error) {
          auto& path = buffers_paths[&buffer - &buffers.front()];
          return load_binary(path_join(dirname, path), buffer, error);
        }))
      return dependent_error();
  }

  // convert asset
  if (gltf.contains("asset")) {
    try {
      scene.copyright = gltf.value("copyright", ""s);
    } catch (...) {
      return parse_error();
    }
  }

  // convert cameras
  auto cameras = vector<camera_data>{};
  if (gltf.contains("cameras")) {
    try {
      for (auto& gcamera : gltf.at("cameras")) {
        auto& camera = cameras.emplace_back();
        auto  type   = gcamera.value("type", "perspective");
        if (type == "orthographic") {
          auto& gortho  = gcamera.at("orthographic");
          auto  xmag    = gortho.value("xmag", 1.0f);
          auto  ymag    = gortho.value("ymag", 1.0f);
          camera.aspect = xmag / ymag;
          camera.lens   = ymag;  // this is probably bogus
          camera.film   = 0.036f;
        } else if (type == "perspective") {
          auto& gpersp  = gcamera.at("perspective");
          camera.aspect = gpersp.value("aspectRatio", 0.0f);
          auto yfov     = gpersp.value("yfov", radians(45.0f));
          if (camera.aspect == 0) camera.aspect = 16.0f / 9.0f;
          camera.film = 0.036f;
          if (camera.aspect >= 1) {
            camera.lens = (camera.film / camera.aspect) / (2 * tan(yfov / 2));
          } else {
            camera.lens = camera.film / (2 * tan(yfov / 2));
          }
          camera.focus = 1;
        } else {
          return parse_error();
        }
      }
    } catch (...) {
      return parse_error();
    }
  }

  // convert color textures
  auto get_texture = [&gltf](
                         const json_value& json, const string& name) -> int {
    if (!json.contains(name)) return invalidid;
    auto& ginfo    = json.at(name);
    auto& gtexture = gltf.at("textures").at(ginfo.value("index", -1));
    return gtexture.value("source", -1);
  };

  // https://stackoverflow.com/questions/3418231/replace-part-of-a-string-with-another-string
  auto replace = [](const std::string& str_, const std::string& from,
                     const std::string& to) -> string {
    auto str = str_;
    if (from.empty()) return str;
    size_t start_pos = 0;
    while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
      str.replace(start_pos, from.length(), to);
      start_pos += to.length();
    }
    return str;
  };

  // convert textures
  auto texture_paths = vector<string>{};
  if (gltf.contains("images")) {
    try {
      for (auto& gimage : gltf.at("images")) {
        scene.textures.emplace_back();
        texture_paths.push_back(replace(gimage.value("uri", ""), "%20", " "));
      }
    } catch (...) {
      return parse_error();
    }
  }

  // convert materials
  if (gltf.contains("materials")) {
    try {
      for (auto& gmaterial : gltf.at("materials")) {
        auto& material    = scene.materials.emplace_back();
        material.type     = material_type::gltfpbr;
        material.emission = gmaterial.value("emissiveFactor", vec3f{0, 0, 0});
        material.emission_tex = get_texture(gmaterial, "emissiveTexture");
        material.normal_tex   = get_texture(gmaterial, "normalTexture");
        if (gmaterial.contains("pbrMetallicRoughness")) {
          auto& gpbr         = gmaterial.at("pbrMetallicRoughness");
          auto  base         = gpbr.value("baseColorFactor", vec4f{1, 1, 1, 1});
          material.color     = xyz(base);
          material.opacity   = base.w;
          material.metallic  = gpbr.value("metallicFactor", 1.0f);
          material.roughness = gpbr.value("roughnessFactor", 1.0f);
          material.color_tex = get_texture(gpbr, "baseColorTexture");
          material.roughness_tex = get_texture(
              gpbr, "metallicRoughnessTexture");
        }
        if (gmaterial.contains("extensions") &&
            gmaterial.at("extensions").contains("KHR_materials_transmission")) {
          auto& gtransmission =
              gmaterial.at("extensions").at("KHR_materials_transmission");
          auto transmission = gtransmission.value("transmissionFactor", 1.0f);
          if (transmission > 0) {
            material.type      = material_type::transparent;
            material.color     = {transmission, transmission, transmission};
            material.color_tex = get_texture(gmaterial, "transmissionTexture");
            // material.roughness = 0; // leave it set from before
          }
        }
      }
    } catch (...) {
      return parse_error();
    }
  }

  // convert meshes
  auto mesh_primitives = vector<vector<instance_data>>{};
  if (gltf.contains("meshes")) {
    try {
      auto type_components = unordered_map<string, int>{
          {"SCALAR", 1}, {"VEC2", 2}, {"VEC3", 3}, {"VEC4", 4}};
      for (auto& gmesh : gltf.at("meshes")) {
        auto& primitives = mesh_primitives.emplace_back();
        if (!gmesh.contains("primitives")) continue;
        for (auto& gprimitive : gmesh.at("primitives")) {
          if (!gprimitive.contains("attributes")) continue;
          auto& shape       = scene.shapes.emplace_back();
          auto& instance    = primitives.emplace_back();
          instance.shape    = (int)scene.shapes.size() - 1;
          instance.material = gprimitive.value("material", -1);
          for (auto& [gname, gattribute] :
              gprimitive.at("attributes").items()) {
            auto& gaccessor = gltf.at("accessors").at(gattribute.get<int>());
            if (gaccessor.contains("sparse")) return parse_error();
            auto& gview =
                gltf.at("bufferViews").at(gaccessor.value("bufferView", -1));
            auto& buffer      = buffers.at(gview.value("buffer", 0));
            auto  components  = type_components.at(gaccessor.value("type", ""));
            auto  dcomponents = components;
            auto  count       = gaccessor.value("count", (size_t)0);
            auto  data        = (float*)nullptr;
            if (gname == "POSITION") {
              if (components != 3) return parse_error();
              shape.positions.resize(count);
              data = (float*)shape.positions.data();
            } else if (gname == "NORMAL") {
              if (components != 3) return parse_error();
              shape.normals.resize(count);
              data = (float*)shape.normals.data();
            } else if (gname == "TEXCOORD" || gname == "TEXCOORD_0") {
              if (components != 2) return parse_error();
              shape.texcoords.resize(count);
              data = (float*)shape.texcoords.data();
            } else if (gname == "COLOR" || gname == "COLOR_0") {
              if (components != 3 && components != 4) return parse_error();
              shape.colors.resize(count);
              data = (float*)shape.colors.data();
              if (components == 3) {
                dcomponents = 4;
                for (auto& c : shape.colors) c.w = 1;
              }
            } else if (gname == "TANGENT") {
              if (components != 4) return parse_error();
              shape.tangents.resize(count);
              data = (float*)shape.tangents.data();
            } else if (gname == "RADIUS") {
              if (components != 1) return parse_error();
              shape.radius.resize(count);
              data = (float*)shape.radius.data();
            } else {
              // ignore
              continue;
            }
            // convert values
            auto current = buffer.data() +
                           gaccessor.value("byteOffset", (size_t)0) +
                           gview.value("byteOffset", (size_t)0);
            auto stride = gaccessor.value("byteStride", (size_t)0);
            auto ctype  = gaccessor.value("componentType", 0);
            if (ctype == 5121) {
              if (stride == 0) stride = components * 1;
              for (auto idx = 0; idx < count; idx++, current += stride) {
                for (auto comp = 0; comp < components; comp++) {
                  data[idx * dcomponents + comp] =
                      *(byte*)(current + comp * 1) / 255.0f;
                }
              }
            } else if (ctype == 5123) {
              if (stride == 0) stride = components * 2;
              for (auto idx = 0; idx < count; idx++, current += stride) {
                for (auto comp = 0; comp < components; comp++) {
                  data[idx * dcomponents + comp] =
                      *(ushort*)(current + comp * 2) / 65535.0f;
                }
              }
            } else if (ctype == 5126) {
              if (stride == 0) stride = components * 4;
              for (auto idx = 0; idx < count; idx++, current += stride) {
                for (auto comp = 0; comp < components; comp++) {
                  data[idx * dcomponents + comp] = *(
                      float*)(current + comp * 4);
                }
              }
            } else {
              return parse_error();
            }
            // fixes
            if (gname == "TANGENT") {
              for (auto& t : shape.tangents) t.w = -t.w;
            }
          }
          // mode
          auto mode = gprimitive.value("mode", 4);
          // indices
          if (!gprimitive.contains("indices")) {
            if (mode == 4) {  // triangles
              shape.triangles.resize(shape.positions.size() / 3);
              for (auto i = 0; i < shape.positions.size() / 3; i++)
                shape.triangles[i] = {i * 3 + 0, i * 3 + 1, i * 3 + 2};
            } else if (mode == 6) {  // fans
              shape.triangles.resize(shape.positions.size() - 2);
              for (auto i = 2; i < shape.positions.size(); i++)
                shape.triangles[i - 2] = {0, i - 1, i};
            } else if (mode == 5) {  // strips
              shape.triangles.resize(shape.positions.size() - 2);
              for (auto i = 2; i < shape.positions.size(); i++)
                shape.triangles[i - 2] = {i - 2, i - 1, i};
            } else if (mode == 1) {  // lines
              shape.lines.resize(shape.positions.size() / 2);
              for (auto i = 0; i < shape.positions.size() / 2; i++)
                shape.lines[i] = {i * 2 + 0, i * 2 + 1};
            } else if (mode == 2) {  // lines loops
              shape.lines.resize(shape.positions.size());
              for (auto i = 1; i < shape.positions.size(); i++)
                shape.lines[i - 1] = {i - 1, i};
              shape.lines.back() = {(int)shape.positions.size() - 1, 0};
            } else if (mode == 3) {  // lines strips
              shape.lines.resize(shape.positions.size() - 1);
              for (auto i = 1; i < shape.positions.size(); i++)
                shape.lines[i - 1] = {i - 1, i};
            } else if (mode == 0) {  // points strips
              return parse_error();
            } else {
              return parse_error();
            }
          } else {
            auto& gaccessor =
                gltf.at("accessors").at(gprimitive.value("indices", -1));
            auto& gview =
                gltf.at("bufferViews").at(gaccessor.value("bufferView", -1));
            auto& buffer = buffers.at(gview.value("buffer", 0));
            if (gaccessor.value("type", "") != "SCALAR") return parse_error();
            auto count   = gaccessor.value("count", (size_t)0);
            auto indices = vector<int>(count);
            // convert values
            auto current = buffer.data() +
                           gaccessor.value("byteOffset", (size_t)0) +
                           gview.value("byteOffset", (size_t)0);
            auto stride = gaccessor.value("byteStride", (size_t)0);
            auto ctype  = gaccessor.value("componentType", 0);
            if (ctype == 5121) {
              if (stride == 0) stride = 1;
              for (auto idx = 0; idx < count; idx++, current += stride) {
                indices[idx] = (int)*(byte*)current;
              }
            } else if (ctype == 5123) {
              if (stride == 0) stride = 2;
              for (auto idx = 0; idx < count; idx++, current += stride) {
                indices[idx] = (int)*(ushort*)current;
              }
            } else if (ctype == 5125) {
              if (stride == 0) stride = 4;
              for (auto idx = 0; idx < count; idx++, current += stride) {
                indices[idx] = (int)*(uint*)current;
              }
            } else {
              return parse_error();
            }
            if (mode == 4) {  // triangles
              shape.triangles.resize(indices.size() / 3);
              for (auto i = 0; i < (int)indices.size() / 3; i++) {
                shape.triangles[i] = {
                    indices[i * 3 + 0], indices[i * 3 + 1], indices[i * 3 + 2]};
              }
            } else if (mode == 6) {  // fans
              shape.triangles.resize(indices.size() - 2);
              for (auto i = 2; i < (int)indices.size(); i++) {
                shape.triangles[i - 2] = {
                    indices[0], indices[i - 1], indices[i + 0]};
              }
            } else if (mode == 5) {  // strips
              shape.triangles.resize(indices.size() - 2);
              for (auto i = 2; i < (int)indices.size(); i++) {
                shape.triangles[i - 2] = {
                    indices[i - 2], indices[i - 1], indices[i + 0]};
              }
            } else if (mode == 1) {  // lines
              shape.lines.resize(indices.size() / 2);
              for (auto i = 0; i < (int)indices.size() / 2; i++) {
                shape.lines[i] = {indices[i * 2 + 0], indices[i * 2 + 1]};
              }
            } else if (mode == 2) {  // lines loops
              shape.lines.resize(indices.size());
              for (auto i = 0; i < (int)indices.size(); i++) {
                shape.lines[i] = {
                    indices[i + 0], indices[i + 1] % (int)indices.size()};
              }
            } else if (mode == 3) {  // lines strips
              shape.lines.resize(indices.size() - 1);
              for (auto i = 0; i < (int)indices.size() - 1; i++) {
                shape.lines[i] = {indices[i + 0], indices[i + 1]};
              }
            } else if (mode == 0) {  // points strips
                                     // points
              return parse_error();
            } else {
              return parse_error();
            }
          }
        }
      }
    } catch (...) {
      return parse_error();
    }
  }

  // convert nodes
  if (gltf.contains("nodes")) {
    try {
      auto parents = vector<int>(gltf.at("nodes").size(), -1);
      auto lxforms = vector<frame3f>(gltf.at("nodes").size(), identity3x4f);
      auto node_id = 0;
      for (auto& gnode : gltf.at("nodes")) {
        auto& xform = lxforms.at(node_id);
        if (gnode.contains("matrix")) {
          xform = mat_to_frame(gnode.value("matrix", identity4x4f));
        }
        if (gnode.contains("scale")) {
          xform = scaling_frame(gnode.value("scale", vec3f{1, 1, 1})) * xform;
        }
        if (gnode.contains("rotation")) {
          xform = rotation_frame(gnode.value("rotation", vec4f{0, 0, 0, 1})) *
                  xform;
        }
        if (gnode.contains("translation")) {
          xform = translation_frame(
                      gnode.value("translation", vec3f{0, 0, 0})) *
                  xform;
        }
        if (gnode.contains("children")) {
          for (auto& gchild : gnode.at("children")) {
            parents.at(gchild.get<int>()) = node_id;
          }
        }
        node_id++;
      }
      auto xforms = vector<frame3f>(gltf.at("nodes").size(), identity3x4f);
      node_id     = 0;
      for (auto& gnode : gltf.at("nodes")) {
        if (!gnode.contains("camera") && !gnode.contains("mesh")) {
          node_id++;
          continue;
        }
        auto& xform = xforms.at(node_id);
        xform       = lxforms.at(node_id);
        auto parent = parents.at(node_id);
        while (parent >= 0) {
          xform  = lxforms.at(parent) * xform;
          parent = parents.at(parent);
        }
        if (gnode.contains("camera")) {
          auto& camera = scene.cameras.emplace_back();
          camera       = cameras.at(gnode.value("camera", -1));
          camera.frame = xform;
        }
        if (gnode.contains("mesh")) {
          for (auto& primitive : mesh_primitives.at(gnode.value("mesh", -1))) {
            auto& instance = scene.instances.emplace_back();
            instance       = primitive;
            instance.frame = xform;
          }
        }
        node_id++;
      }
    } catch (...) {
      return parse_error();
    }
  }

  if (noparallel) {
    // load texture
    for (auto& texture : scene.textures) {
      auto& path = texture_paths[&texture - &scene.textures.front()];
      if (!load_texture(path_join(dirname, path), texture, error))
        return dependent_error();
    }
  } else {
    // load textures
    if (!parallel_foreach(
            scene.textures, error, [&](auto& texture, string& error) {
              auto& path = texture_paths[&texture - &scene.textures.front()];
              return load_texture(path_join(dirname, path), texture, error);
            }))
      return dependent_error();
  }

  // fix scene
  add_missing_material(scene);
  add_missing_camera(scene);
  add_missing_radius(scene);

  // done
  return true;
}

// Load a scene
static bool save_gltf_scene(const string& filename, const scene_data& scene,
    string& error, bool noparallel) {
  // convert scene to json
  auto gltf = json_value::object();

  // asset
  {
    auto& gasset        = gltf["asset"];
    gasset              = json_value::object();
    gasset["version"]   = "2.0";
    gasset["generator"] = "Yocto/GL - https://github.com/xelatihy/yocto-gl";
    if (!scene.copyright.empty()) gasset["copyright"] = scene.copyright;
  }

  // cameras
  if (!scene.cameras.empty()) {
    auto& gcameras = gltf["cameras"];
    gcameras       = json_value::array();
    for (auto& camera : scene.cameras) {
      auto& gcamera               = gcameras.emplace_back();
      gcamera                     = json_value::object();
      gcamera["name"]             = get_camera_name(scene, camera);
      gcamera["type"]             = "perspective";
      auto& gperspective          = gcamera["perspective"];
      gperspective                = json_value::object();
      gperspective["aspectRatio"] = camera.aspect;
      gperspective["yfov"]        = 0.660593;  // TODO(fabio): yfov
      gperspective["znear"]       = 0.001;     // TODO(fabio): configurable?
    }
  }

  // textures
  if (!scene.textures.empty()) {
    auto& gtextures  = gltf["textures"];
    gtextures        = json_value::array();
    auto& gsamplers  = gltf["samplers"];
    gsamplers        = json_value::array();
    auto& gimages    = gltf["images"];
    gimages          = json_value::array();
    auto& gsampler   = gsamplers.emplace_back();
    gsampler         = json_value::object();
    gsampler["name"] = "sampler";
    for (auto& texture : scene.textures) {
      auto  name          = get_texture_name(scene, texture);
      auto& gimage        = gimages.emplace_back();
      gimage              = json_value::object();
      gimage["name"]      = name;
      gimage["uri"]       = "textures/" + name + ".png";
      auto& gtexture      = gtextures.emplace_back();
      gtexture            = json_value::object();
      gtexture["name"]    = name;
      gtexture["sampler"] = 0;
      gtexture["source"]  = (int)gimages.size() - 1;
    }
  }

  // materials
  if (!scene.materials.empty()) {
    auto& gmaterials = gltf["materials"];
    gmaterials       = json_value::array();
    for (auto& material : scene.materials) {
      auto& gmaterial             = gmaterials.emplace_back();
      gmaterial                   = json_value::object();
      gmaterial["name"]           = get_material_name(scene, material);
      gmaterial["emissiveFactor"] = material.emission;
      auto& gpbr                  = gmaterial["pbrMetallicRoughness"];
      gpbr                        = json_value::object();
      gpbr["baseColorFactor"]     = vec4f{material.color.x, material.color.y,
          material.color.z, material.opacity};
      gpbr["metallicFactor"]      = material.metallic;
      gpbr["roughnessFactor"]     = material.roughness;
      if (material.emission_tex != invalidid) {
        gmaterial["emissiveTexture"]          = json_value::object();
        gmaterial["emissiveTexture"]["index"] = material.emission_tex;
      }
      if (material.normal_tex != invalidid) {
        gmaterial["normalTexture"]          = json_value::object();
        gmaterial["normalTexture"]["index"] = material.normal_tex;
      }
      if (material.color_tex != invalidid) {
        gpbr["baseColorTexture"]          = json_value::object();
        gpbr["baseColorTexture"]["index"] = material.color_tex;
      }
      if (material.roughness_tex != invalidid) {
        gpbr["metallicRoughnessTexture"]          = json_value::object();
        gpbr["metallicRoughnessTexture"]["index"] = material.roughness_tex;
      }
    }
  }

  // add an accessor
  auto set_view = [](json_value& gview, json_value& gbuffer, const auto& data,
                      size_t buffer_id) {
    gview                 = json_value::object();
    gview["buffer"]       = buffer_id;
    gview["byteLength"]   = data.size() * sizeof(data.front());
    gview["byteOffset"]   = gbuffer["byteLength"];
    gbuffer["byteLength"] = gbuffer.value("byteLength", (size_t)0) +
                            data.size() * sizeof(data.front());
  };
  auto set_vaccessor = [](json_value& gaccessor, const auto& data,
                           size_t view_id, bool minmax = false) {
    static auto types = unordered_map<size_t, string>{
        {1, "SCALAR"}, {2, "VEC2"}, {3, "VEC3"}, {4, "VEC4"}};
    gaccessor                  = json_value::object();
    gaccessor["bufferView"]    = view_id;
    gaccessor["componentType"] = 5126;
    gaccessor["count"]         = data.size();
    gaccessor["type"]          = types.at(sizeof(data.front()) / sizeof(float));
    if constexpr (sizeof(data.front()) == sizeof(vec3f)) {
      if (minmax) {
        auto bbox = invalidb3f;
        for (auto& value : data) bbox = merge(bbox, value);
        gaccessor["min"] = bbox.min;
        gaccessor["max"] = bbox.max;
      }
    }
  };
  auto set_iaccessor = [](json_value& gaccessor, const auto& data,
                           size_t view_id, bool minmax = false) {
    gaccessor                  = json_value::object();
    gaccessor["bufferView"]    = view_id;
    gaccessor["componentType"] = 5125;
    gaccessor["count"] = data.size() * sizeof(data.front()) / sizeof(int);
    gaccessor["type"]  = "SCALAR";
  };

  // meshes
  auto shape_primitives = vector<json_value>();
  shape_primitives.reserve(scene.shapes.size());
  if (!scene.shapes.empty()) {
    auto& gaccessors = gltf["accessors"];
    gaccessors       = json_value::array();
    auto& gviews     = gltf["bufferViews"];
    gviews           = json_value::array();
    auto& gbuffers   = gltf["buffers"];
    gbuffers         = json_value::array();
    for (auto& shape : scene.shapes) {
      auto& gbuffer         = gbuffers.emplace_back();
      gbuffer["uri"]        = "shapes/" + get_shape_name(scene, shape) + ".bin";
      gbuffer["byteLength"] = (size_t)0;
      auto& gprimitive      = shape_primitives.emplace_back();
      gprimitive            = json_value::object();
      auto& gattributes     = gprimitive["attributes"];
      gattributes           = json_value::object();
      if (!shape.positions.empty()) {
        set_view(gviews.emplace_back(), gbuffer, shape.positions,
            gbuffers.size() - 1);
        set_vaccessor(gaccessors.emplace_back(), shape.positions,
            gviews.size() - 1, true);
        gattributes["POSITION"] = (int)gaccessors.size() - 1;
      }
      if (!shape.normals.empty()) {
        set_view(
            gviews.emplace_back(), gbuffer, shape.normals, gbuffers.size() - 1);
        set_vaccessor(
            gaccessors.emplace_back(), shape.normals, gviews.size() - 1);
        gattributes["NORMAL"] = (int)gaccessors.size() - 1;
      }
      if (!shape.texcoords.empty()) {
        set_view(gviews.emplace_back(), gbuffer, shape.texcoords,
            gbuffers.size() - 1);
        set_vaccessor(
            gaccessors.emplace_back(), shape.texcoords, gviews.size() - 1);
        gattributes["TEXCOORD_0"] = (int)gaccessors.size() - 1;
      }
      if (!shape.colors.empty()) {
        set_view(
            gviews.emplace_back(), gbuffer, shape.colors, gbuffers.size() - 1);
        set_vaccessor(
            gaccessors.emplace_back(), shape.colors, gviews.size() - 1);
        gattributes["COLOR_0"] = (int)gaccessors.size() - 1;
      }
      if (!shape.radius.empty()) {
        set_view(
            gviews.emplace_back(), gbuffer, shape.radius, gbuffers.size() - 1);
        set_vaccessor(
            gaccessors.emplace_back(), shape.radius, gviews.size() - 1);
        gattributes["RADIUS"] = (int)gaccessors.size() - 1;
      }
      if (!shape.points.empty()) {
        set_view(
            gviews.emplace_back(), gbuffer, shape.points, gbuffers.size() - 1);
        set_iaccessor(
            gaccessors.emplace_back(), shape.points, gviews.size() - 1);
        gprimitive["indices"] = (int)gaccessors.size() - 1;
        gprimitive["mode"]    = 0;
      } else if (!shape.lines.empty()) {
        set_view(
            gviews.emplace_back(), gbuffer, shape.lines, gbuffers.size() - 1);
        set_iaccessor(
            gaccessors.emplace_back(), shape.lines, gviews.size() - 1);
        gprimitive["indices"] = (int)gaccessors.size() - 1;
        gprimitive["mode"]    = 1;
      } else if (!shape.triangles.empty()) {
        set_view(gviews.emplace_back(), gbuffer, shape.triangles,
            gbuffers.size() - 1);
        set_iaccessor(
            gaccessors.emplace_back(), shape.triangles, gviews.size() - 1);
        gprimitive["indices"] = (int)gaccessors.size() - 1;
        gprimitive["mode"]    = 4;
      } else if (!shape.quads.empty()) {
        auto triangles = quads_to_triangles(shape.quads);
        set_view(
            gviews.emplace_back(), gbuffer, triangles, gbuffers.size() - 1);
        set_iaccessor(gaccessors.emplace_back(), triangles, gviews.size() - 1);
        gprimitive["indices"] = (int)gaccessors.size() - 1;
        gprimitive["mode"]    = 4;
      }
    }
  }

  // meshes
  using mesh_key = pair<int, int>;
  struct mesh_key_hash {
    size_t operator()(const mesh_key& v) const {
      const std::hash<int> hasher = std::hash<int>();
      auto                 h      = (size_t)0;
      h ^= hasher(v.first) + 0x9e3779b9 + (h << 6) + (h >> 2);
      h ^= hasher(v.second) + 0x9e3779b9 + (h << 6) + (h >> 2);
      return h;
    }
  };
  auto mesh_map = unordered_map<mesh_key, size_t, mesh_key_hash>{};
  if (!scene.instances.empty()) {
    auto& gmeshes = gltf["meshes"];
    gmeshes       = json_value::array();
    for (auto& instance : scene.instances) {
      auto key = mesh_key{instance.shape, instance.material};
      if (mesh_map.find(key) != mesh_map.end()) continue;
      auto& gmesh   = gmeshes.emplace_back();
      gmesh         = json_value::object();
      gmesh["name"] = get_shape_name(scene, instance.shape) + "_" +
                      get_material_name(scene, instance.material);
      gmesh["primitives"] = json_value::array();
      gmesh["primitives"].push_back(shape_primitives.at(instance.shape));
      gmesh["primitives"].back()["material"] = instance.material;
      mesh_map[key]                          = gmeshes.size() - 1;
    }
  } else if (!scene.shapes.empty()) {
    auto& gmeshes = gltf["meshes"];
    gmeshes       = json_value::array();
    auto shape_id = 0;
    for (auto& primitives : shape_primitives) {
      auto& gmesh         = gmeshes.emplace_back();
      gmesh               = json_value::object();
      gmesh["name"]       = get_shape_name(scene, shape_id++);
      gmesh["primitives"] = json_value::array();
      gmesh["primitives"].push_back(primitives);
    }
  }

  // nodes
  if (!scene.cameras.empty() || !scene.instances.empty()) {
    auto& gnodes   = gltf["nodes"];
    gnodes         = json_value::array();
    auto camera_id = 0;
    for (auto& camera : scene.cameras) {
      auto& gnode     = gnodes.emplace_back();
      gnode           = json_value::object();
      gnode["name"]   = get_camera_name(scene, camera);
      gnode["matrix"] = frame_to_mat(camera.frame);
      gnode["camera"] = camera_id++;
    }
    for (auto& instance : scene.instances) {
      auto& gnode     = gnodes.emplace_back();
      gnode           = json_value::object();
      gnode["name"]   = get_instance_name(scene, instance);
      gnode["matrix"] = frame_to_mat(instance.frame);
      gnode["mesh"]   = mesh_map.at({instance.shape, instance.material});
    }
    // root children
    auto& groot     = gnodes.emplace_back();
    groot           = json_value::object();
    groot["name"]   = "root";
    auto& gchildren = groot["children"];
    gchildren       = json_value::array();
    for (auto idx = (size_t)0; idx < gnodes.size() - 1; idx++)
      gchildren.push_back(idx);
    // scene
    auto& gscenes     = gltf["scenes"];
    gscenes           = json_value::array();
    auto& gscene      = gscenes.emplace_back();
    gscene            = json_value::object();
    auto& gscenenodes = gscene["nodes"];
    gscenenodes       = json_value::array();
    gscenenodes.push_back(gnodes.size() - 1);
    gltf["scene"] = 0;
  }

  // save json
  if (!save_json(filename, gltf, error)) return false;

  // get filename from name
  auto make_filename = [filename](const string& name, const string& group,
                           const string& extension) {
    return path_join(path_dirname(filename), group, name + extension);
  };

  // dirname
  auto dirname         = path_dirname(filename);
  auto dependent_error = [&filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  if (noparallel) {
    // save shapes
    for (auto& shape : scene.shapes) {
      auto path = "shapes/" + get_shape_name(scene, shape) + ".bin";
      if (!save_binshape(path_join(dirname, path), shape, error))
        return dependent_error();
    }
    // save textures
    for (auto& texture : scene.textures) {
      auto path = "textures/" + get_texture_name(scene, texture) +
                  (!texture.pixelsf.empty() ? ".hdr" : ".png");
      if (!save_texture(path_join(dirname, path), texture, error))
        return dependent_error();
    }
  } else {
    // save shapes
    if (!parallel_foreach(scene.shapes, error, [&](auto& shape, string& error) {
          auto path = "shapes/" + get_shape_name(scene, shape) + ".bin";
          return save_binshape(path_join(dirname, path), shape, error);
        }))
      return dependent_error();
    // save textures
    if (!parallel_foreach(
            scene.textures, error, [&](auto& texture, string& error) {
              auto path = "textures/" + get_texture_name(scene, texture) +
                          (!texture.pixelsf.empty() ? ".hdr"s : ".png"s);
              return save_texture(path_join(dirname, path), texture, error);
            }))
      return dependent_error();
  }

  // done
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF PBRT
// -----------------------------------------------------------------------------
namespace yocto {

// load pbrt scenes
static bool load_pbrt_scene(
    const string& filename, scene_data& scene, string& error, bool noparallel) {
  // load pbrt
  auto pbrt = pbrt_model{};
  if (!load_pbrt(filename, pbrt, error)) return false;

  // convert cameras
  for (auto& pcamera : pbrt.cameras) {
    auto& camera  = scene.cameras.emplace_back();
    camera.frame  = pcamera.frame;
    camera.aspect = pcamera.aspect;
    camera.film   = 0.036f;
    camera.lens   = pcamera.lens;
    camera.focus  = pcamera.focus;
  }

  // convert material
  auto texture_paths = vector<string>{};
  for (auto& ptexture : pbrt.textures) {
    scene.textures.emplace_back();
    texture_paths.push_back(ptexture.filename);
  }

  // material type map
  auto material_type_map = unordered_map<pbrt_mtype, material_type>{
      {pbrt_mtype::matte, material_type::matte},
      {pbrt_mtype::plastic, material_type::glossy},
      {pbrt_mtype::metal, material_type::reflective},
      {pbrt_mtype::glass, material_type::refractive},
      {pbrt_mtype::thinglass, material_type::transparent},
      {pbrt_mtype::subsurface, material_type::matte},
  };

  // convert material
  for (auto& pmaterial : pbrt.materials) {
    auto& material = scene.materials.emplace_back();
    material.type  = material_type_map.at(pmaterial.type);
    if (pmaterial.emission != vec3f{0, 0, 0}) {
      material.type = material_type::matte;
    }
    material.emission  = pmaterial.emission;
    material.color     = pmaterial.color;
    material.ior       = pmaterial.ior;
    material.roughness = pmaterial.roughness;
    material.opacity   = pmaterial.opacity;
    material.color_tex = pmaterial.color_tex;
  }

  // convert shapes
  auto shapes_paths = vector<string>{};
  for (auto& pshape : pbrt.shapes) {
    auto& shape = scene.shapes.emplace_back();
    shapes_paths.emplace_back(pshape.filename_);
    shape.positions = pshape.positions;
    shape.normals   = pshape.normals;
    shape.texcoords = pshape.texcoords;
    shape.triangles = pshape.triangles;
    for (auto& uv : shape.texcoords) uv.y = 1 - uv.y;
    if (!pshape.instanced) {
      auto& instance    = scene.instances.emplace_back();
      instance.frame    = pshape.frame;
      instance.shape    = (int)scene.shapes.size() - 1;
      instance.material = pshape.material;
    } else {
      for (auto& frame : pshape.instances) {
        auto& instance    = scene.instances.emplace_back();
        instance.frame    = frame * pshape.frame;
        instance.shape    = (int)scene.shapes.size() - 1;
        instance.material = pshape.material;
      }
    }
  }

  // convert environments
  for (auto& penvironment : pbrt.environments) {
    auto& environment        = scene.environments.emplace_back();
    environment.frame        = penvironment.frame;
    environment.emission     = penvironment.emission;
    environment.emission_tex = penvironment.emission_tex;
  }

  // lights
  for (auto& plight : pbrt.lights) {
    auto& shape = scene.shapes.emplace_back();
    shapes_paths.emplace_back();
    shape.triangles   = plight.area_triangles;
    shape.positions   = plight.area_positions;
    shape.normals     = plight.area_normals;
    auto& material    = scene.materials.emplace_back();
    material.emission = plight.area_emission;
    auto& instance    = scene.instances.emplace_back();
    instance.shape    = (int)scene.shapes.size() - 1;
    instance.material = (int)scene.materials.size() - 1;
    instance.frame    = plight.area_frame;
  }

  // dirname
  auto dirname         = path_dirname(filename);
  auto dependent_error = [&filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  if (noparallel) {
    // load shape
    for (auto& shape : scene.shapes) {
      auto& path = shapes_paths[&shape - &scene.shapes.front()];
      if (path.empty()) continue;
      if (!load_shape(path_join(dirname, path), shape, error, true))
        return dependent_error();
    }
    // load texture
    for (auto& texture : scene.textures) {
      auto& path = texture_paths[&texture - &scene.textures.front()];
      if (!load_texture(path_join(dirname, path), texture, error))
        return dependent_error();
    }
  } else {
    // load shapes
    if (!parallel_foreach(scene.shapes, error, [&](auto& shape, string& error) {
          auto& path = shapes_paths[&shape - &scene.shapes.front()];
          if (path.empty()) return true;
          return load_shape(path_join(dirname, path), shape, error, true);
        }))
      return dependent_error();
    // load textures
    if (!parallel_foreach(
            scene.textures, error, [&](auto& texture, string& error) {
              auto& path = texture_paths[&texture - &scene.textures.front()];
              return load_texture(path_join(dirname, path), texture, error);
            }))
      return dependent_error();
  }

  // fix scene
  add_missing_camera(scene);
  add_missing_radius(scene);

  // done
  return true;
}

// Save a pbrt scene
static bool save_pbrt_scene(const string& filename, const scene_data& scene,
    string& error, bool noparallel) {
  // save pbrt
  auto pbrt = pbrt_model{};

  // convert camera
  auto& camera       = scene.cameras.front();
  auto& pcamera      = pbrt.cameras.emplace_back();
  pcamera.frame      = camera.frame;
  pcamera.lens       = camera.lens;
  pcamera.aspect     = camera.aspect;
  pcamera.resolution = {1280, (int)(1280 / pcamera.aspect)};

  // convert textures
  for (auto& texture : scene.textures) {
    auto& ptexture    = pbrt.textures.emplace_back();
    ptexture.filename = "textures/" + get_texture_name(scene, texture) +
                        (!texture.pixelsf.empty() ? ".hdr" : ".png");
  }

  // material type map
  auto material_type_map = unordered_map<material_type, pbrt_mtype>{
      {material_type::matte, pbrt_mtype::matte},
      {material_type::glossy, pbrt_mtype::plastic},
      {material_type::reflective, pbrt_mtype::metal},
      {material_type::refractive, pbrt_mtype::glass},
      {material_type::transparent, pbrt_mtype::thinglass},
      {material_type::subsurface, pbrt_mtype::matte},
      {material_type::volumetric, pbrt_mtype::matte},
  };

  // convert materials
  for (auto& material : scene.materials) {
    auto& pmaterial     = pbrt.materials.emplace_back();
    pmaterial.name      = get_material_name(scene, material);
    pmaterial.type      = material_type_map.at(material.type);
    pmaterial.emission  = material.emission;
    pmaterial.color     = material.color;
    pmaterial.roughness = material.roughness;
    pmaterial.ior       = material.ior;
    pmaterial.opacity   = material.opacity;
    pmaterial.color_tex = material.color_tex;
  }

  // convert instances
  for (auto& instance : scene.instances) {
    auto& pshape     = pbrt.shapes.emplace_back();
    pshape.filename_ = get_shape_name(scene, instance.shape) + ".ply";
    pshape.frame     = instance.frame;
    pshape.frend     = instance.frame;
    pshape.material  = instance.material;
  }

  // convert environments
  for (auto& environment : scene.environments) {
    auto& penvironment        = pbrt.environments.emplace_back();
    penvironment.emission     = environment.emission;
    penvironment.emission_tex = environment.emission_tex;
  }

  // save pbrt
  if (!save_pbrt(filename, pbrt, error)) return false;

  // dirname
  auto dirname         = path_dirname(filename);
  auto dependent_error = [&filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  if (noparallel) {
    // save textures
    for (auto& shape : scene.shapes) {
      auto path = "shapes/" + get_shape_name(scene, shape) + ".ply";
      if (!save_shape(path_join(dirname, path), shape, error, true))
        return dependent_error();
    }
    // save shapes
    for (auto& texture : scene.textures) {
      auto path = "textures/" + get_texture_name(scene, texture) +
                  (!texture.pixelsf.empty() ? ".hdr" : ".png");
      if (!save_texture(path_join(dirname, path), texture, error))
        return dependent_error();
    }
  } else {
    // save shapes
    if (!parallel_foreach(scene.shapes, error, [&](auto& shape, string& error) {
          auto path = "shapes/" + get_shape_name(scene, shape) + ".ply";
          return save_shape(path_join(dirname, path), shape, error, true);
        }))
      return dependent_error();
    // save textures
    if (!parallel_foreach(
            scene.textures, error, [&](auto& texture, string& error) {
              auto path = "textures/" + get_texture_name(scene, texture) +
                          (!texture.pixelsf.empty() ? ".hdr"s : ".png"s);
              return save_texture(path_join(dirname, path), texture, error);
            }))
      return dependent_error();
  }

  // done
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON CLI
// -----------------------------------------------------------------------------
namespace yocto {

// Using directive
using ordered_json = nlohmann::ordered_json;

// Parse command line arguments to Json without schema
static bool cli_to_json_value(ordered_json& json, const string& arg) {
  if (arg.empty()) throw std::invalid_argument("should not have gotten here");
  json = ordered_json::parse(arg, nullptr, false);
  if (json.is_discarded()) json = arg;
  return true;
}
static pair<bool, int> cli_to_json_option(
    ordered_json& json, const vector<string>& args, int pos) {
  if (pos >= args.size() || args[pos].find("--") == 0) {
    json = true;
    return {true, pos};
  } else {
    while (pos < (int)args.size() && args[pos].find("--") != 0) {
      if (json.is_array()) {
        if (!cli_to_json_value(json.emplace_back(), args[pos++]))
          return {false, pos};
      } else if (json.is_null()) {
        if (!cli_to_json_value(json, args[pos++])) return {false, pos};
      } else {
        auto item = json;
        json      = ordered_json::array();
        json.push_back(item);
        if (!cli_to_json_value(json.emplace_back(), args[pos++]))
          return {false, pos};
      }
    }
    return {true, pos};
  }
}
static bool cli_to_json_command(
    ordered_json& json, const vector<string>& args, int pos) {
  if (pos >= args.size()) return true;
  if (args[pos].find("--") == 0) {
    while (pos < (int)args.size() && args[pos].find("--") == 0) {
      auto result = cli_to_json_option(
          json[args[pos].substr(2)], args, pos + 1);
      if (!result.first) return false;
      pos = result.second;
    }
    return true;
  } else {
    return cli_to_json_command(json[args[pos]], args, pos + 1);
  }
}
bool cli_to_json(ordered_json& json, const vector<string>& args) {
  return cli_to_json_command(json, args, 1);
}
bool cli_to_json(ordered_json& json, int argc, const char** argv) {
  return cli_to_json(json, vector<string>{argv, argv + argc});
}

// Validate Cli Json against a schema
bool validate_cli(const ordered_json& json, const ordered_json& schema);

// Get Cli usage from Json
string cli_usage(const ordered_json& json, const ordered_json& schema);

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON SUPPORT FOR YOCTO TYPES
// -----------------------------------------------------------------------------
namespace yocto {

// Validate Json against a schema
bool validate_json(const ordered_json& json, const ordered_json& schema);

}  // namespace yocto

#if 0

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR VOLUME IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Volume load
static bool load_yvol(const string& filename, int& width, int& height,
    int& depth, int& components, vector<float>& voxels, string& error) {
  // error helpers
  auto open_error = [filename, &error]() {
    error = filename + ": file not found";
    return false;
  };
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error";
    return false;
  };
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };

  // Split a string
  auto split_string = [](const string& str) -> vector<string> {
    auto ret = vector<string>();
    if (str.empty()) return ret;
    auto lpos = (size_t)0;
    while (lpos != string::npos) {
      auto pos = str.find_first_of(" \t\n\r", lpos);
      if (pos != string::npos) {
        if (pos > lpos) ret.push_back(str.substr(lpos, pos - lpos));
        lpos = pos + 1;
      } else {
        if (lpos < str.size()) ret.push_back(str.substr(lpos));
        lpos = pos;
      }
    }
    return ret;
  };

  auto fs       = fopen_utf8(filename.c_str(), "rb");
  auto fs_guard = unique_ptr<FILE, int (*)(FILE*)>(fs, &fclose);
  if (!fs) return open_error();

  // buffer
  auto buffer = array<char, 4096>{};
  auto toks   = vector<string>();

  // read magic
  if (!fgets(buffer.data(), (int)buffer.size(), fs)) return parse_error();
  toks = split_string(buffer.data());
  if (toks[0] != "YVOL") return parse_error();

  // read width, height
  if (!fgets(buffer.data(), (int)buffer.size(), fs)) return parse_error();
  toks       = split_string(buffer.data());
  width      = atoi(toks[0].c_str());
  height     = atoi(toks[1].c_str());
  depth      = atoi(toks[2].c_str());
  components = atoi(toks[3].c_str());

  // read data
  auto nvoxels = (size_t)width * (size_t)height * (size_t)depth;
  auto nvalues = nvoxels * (size_t)components;
  voxels       = vector<float>(nvalues);
  if (!read_values(fs, voxels.data(), nvalues)) return read_error();

  // done
  return true;
}

// save pfm
static bool save_yvol(const string& filename, int width, int height, int depth,
    int components, const vector<float>& voxels, string& error) {
  // error helpers
  auto open_error = [filename, &error]() {
    error = filename + ": file not found";
    return false;
  };
  auto write_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };

  auto fs       = fopen_utf8(filename.c_str(), "wb");
  auto fs_guard = unique_ptr<FILE, int (*)(FILE*)>(fs, &fclose);
  if (!fs) return open_error();

  if (!write_text(fs, "YVOL\n")) return write_error();
  if (!write_text(fs, std::to_string(width) + " " + std::to_string(height) +
                          " " + std::to_string(depth) + " " +
                          std::to_string(components) + "\n"))
    return write_error();
  auto nvalues = (size_t)width * (size_t)height * (size_t)depth *
                 (size_t)components;
  if (!write_values(fs, voxels.data(), nvalues)) return write_error();
  return true;
}

// Loads volume data from binary format.
bool load_volume(const string& filename, volume<float>& vol, string& error) {
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };
  auto width = 0, height = 0, depth = 0, ncomp = 0;
  auto voxels = vector<float>{};
  if (!load_yvol(filename, width, height, depth, ncomp, voxels, error))
    return false;
  if (ncomp != 1) voxels = convert_components(voxels, ncomp, 1);
  vol = volume{{width, height, depth}, (const float*)voxels.data()};
  return true;
}

// Saves volume data in binary format.
bool save_volume(
    const string& filename, const volume<float>& vol, string& error) {
  return save_yvol(filename, vol.width(), vol.height(), vol.depth(), 1,
      {vol.data(), vol.data() + vol.count()}, error);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR DEPRECATED CODE
// -----------------------------------------------------------------------------
namespace yocto {

image<vec4f> filter_bilateral(const image<vec4f>& img, float spatial_sigma,
    float range_sigma, const vector<image<vec4f>>& features,
    const vector<float>& features_sigma) {
  auto filtered     = image{img.imsize(), zero4f};
  auto filter_width = (int)ceil(2.57f * spatial_sigma);
  auto sw           = 1 / (2.0f * spatial_sigma * spatial_sigma);
  auto rw           = 1 / (2.0f * range_sigma * range_sigma);
  auto fw           = vector<float>();
  for (auto feature_sigma : features_sigma)
    fw.push_back(1 / (2.0f * feature_sigma * feature_sigma));
  for (auto j = 0; j < img.height(); j++) {
    for (auto i = 0; i < img.width(); i++) {
      auto av = zero4f;
      auto aw = 0.0f;
      for (auto fj = -filter_width; fj <= filter_width; fj++) {
        for (auto fi = -filter_width; fi <= filter_width; fi++) {
          auto ii = i + fi, jj = j + fj;
          if (ii < 0 || jj < 0) continue;
          if (ii >= img.width() || jj >= img.height()) continue;
          auto uv  = vec2f{float(i - ii), float(j - jj)};
          auto rgb = img[{i, j}] - img[{i, j}];
          auto w   = (float)exp(-dot(uv, uv) * sw) *
                   (float)exp(-dot(rgb, rgb) * rw);
          for (auto fi = 0; fi < features.size(); fi++) {
            auto feat = features[fi][{i, j}] - features[fi][{i, j}];
            w *= exp(-dot(feat, feat) * fw[fi]);
          }
          av += w * img[{ii, jj}];
          aw += w;
        }
      }
      filtered[{i, j}] = av / aw;
    }
  }
  return filtered;
}

image<vec4f> filter_bilateral(
    const image<vec4f>& img, float spatial_sigma, float range_sigma) {
  auto filtered = image{img.imsize(), zero4f};
  auto fwidth   = (int)ceil(2.57f * spatial_sigma);
  auto sw       = 1 / (2.0f * spatial_sigma * spatial_sigma);
  auto rw       = 1 / (2.0f * range_sigma * range_sigma);
  for (auto j = 0; j < img.height(); j++) {
    for (auto i = 0; i < img.width(); i++) {
      auto av = zero4f;
      auto aw = 0.0f;
      for (auto fj = -fwidth; fj <= fwidth; fj++) {
        for (auto fi = -fwidth; fi <= fwidth; fi++) {
          auto ii = i + fi, jj = j + fj;
          if (ii < 0 || jj < 0) continue;
          if (ii >= img.width() || jj >= img.height()) continue;
          auto uv  = vec2f{float(i - ii), float(j - jj)};
          auto rgb = img[{i, j}] - img[{ii, jj}];
          auto w   = exp(-dot(uv, uv) * sw) * exp(-dot(rgb, rgb) * rw);
          av += w * img[{ii, jj}];
          aw += w;
        }
      }
      filtered[{i, j}] = av / aw;
    }
  }
  return filtered;
}

}  // namespace yocto

#endif
