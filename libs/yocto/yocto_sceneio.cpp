//
// Implementation for Yocto/Scene Input and Output functions.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2022 Fabio Pellacini
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

#include <cgltf/cgltf.h>
#include <cgltf/cgltf_write.h>
#include <stb_image/stb_image.h>
#include <stb_image/stb_image_resize.h>
#include <stb_image/stb_image_write.h>
#include <tinyexr/tinyexr.h>

#include <algorithm>
#include <cassert>
#include <cctype>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <future>
#include <memory>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <thread>
#include <unordered_map>

#include "yocto_color.h"
#include "yocto_geometry.h"
#include "yocto_image.h"
#include "yocto_modelio.h"
#include "yocto_pbrtio.h"
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
// PARALLEL HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename T, typename Func>
inline bool parallel_for(T num, Func&& func) {
  auto              futures  = vector<std::future<void>>{};
  auto              nthreads = std::thread::hardware_concurrency();
  std::atomic<T>    next_idx(0);
  std::atomic<bool> has_error(false);
  for (auto thread_id = 0; thread_id < (int)nthreads; thread_id++) {
    futures.emplace_back(
        std::async(std::launch::async, [&func, &next_idx, &has_error, num]() {
          while (true) {
            if (has_error) break;
            auto idx = next_idx.fetch_add(1);
            if (idx >= num) break;
            try {
              func(idx);
            } catch (std::exception& error) {
              has_error = true;
              throw;
            }
          }
        }));
  }
  for (auto& f : futures) f.get();
  return !(bool)has_error;
}

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename Sequence1, typename Sequence2, typename Func>
inline void parallel_zip(
    Sequence1&& sequence1, Sequence2&& sequence2, Func&& func) {
  if (std::size(sequence1) != std::size(sequence2))
    throw std::out_of_range{"invalid sequence lengths"};
  auto                num      = std::size(sequence1);
  auto                futures  = vector<std::future<void>>{};
  auto                nthreads = std::thread::hardware_concurrency();
  std::atomic<size_t> next_idx(0);
  std::atomic<bool>   has_error(false);
  for (auto thread_id = 0; thread_id < (int)nthreads; thread_id++) {
    futures.emplace_back(std::async(std::launch::async,
        [&func, &next_idx, &has_error, num, &sequence1, &sequence2]() {
          try {
            while (true) {
              auto idx = next_idx.fetch_add(1);
              if (idx >= num) break;
              if (has_error) break;
              func(std::forward<Sequence1>(sequence1)[idx],
                  std::forward<Sequence2>(sequence2)[idx]);
            }
          } catch (...) {
            has_error = true;
            throw;
          }
        }));
  }
  for (auto& f : futures) f.get();
}

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes a reference to a `T`.
template <typename T, typename Func>
inline bool parallel_foreach(vector<T>& values, Func&& func) {
  return parallel_for(values.size(),
      [&func, &values](size_t idx) { return func(values[idx]); });
}
template <typename T, typename Func>
inline bool parallel_foreach(const vector<T>& values, Func&& func) {
  return parallel_for(values.size(),
      [&func, &values](size_t idx) { return func(values[idx]); });
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PATH HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

// Make a path from a utf8 string
static std::filesystem::path to_path(const string& filename) {
  auto filename8 = std::u8string((char8_t*)filename.data(), filename.size());
  return std::filesystem::path(filename8);
}

// Make a utf8 string from a path
static string to_string(const std::filesystem::path& path) {
  auto string8 = path.u8string();
  return string((char*)string8.data(), string8.size());
}

// Normalize a path
string path_normalized(const string& path) { return to_string(to_path(path)); }

// Get directory name (not including /)
string path_dirname(const string& path) {
  return to_string(to_path(path).parent_path());
}

// Get filename without directory and extension.
string path_basename(const string& path) {
  return to_string(to_path(path).stem());
}

// Get extension
string path_extension(const string& path) {
  return to_string(to_path(path).extension());
}

// Check if a file can be opened for reading.
bool path_exists(const string& path) { return exists(to_path(path)); }

// Replace the extension of a file
string replace_extension(const string& path, const string& extension) {
  auto ext = to_path(extension).extension();
  return to_string(to_path(path).replace_extension(ext));
}

// Create a directory and all missing parent directories if needed
void make_directory(const string& path) {
  if (path_exists(path)) return;
  try {
    create_directories(to_path(path));
  } catch (...) {
    throw io_error{"cannot create directory " + path};
  }
}

// Joins paths
static string path_join(const string& patha, const string& pathb) {
  return to_string(to_path(patha) / to_path(pathb));
}
static string path_join(
    const string& patha, const string& pathb, const string& pathc) {
  return to_string(to_path(patha) / to_path(pathb) / to_path(pathc));
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// FILE WATCHER
// -----------------------------------------------------------------------------
namespace yocto {

// Initialize file watcher
watch_context make_watch_context(const vector<string>& filenames, int delay) {
  return {{0}, {}, filenames, vector<int64_t>(filenames.size(), 0),
      (int64_t)delay, {false}};
}

// Start file watcher
void watch_start(watch_context& context) {
  // stop
  if (context.worker.valid()) context.worker.get();
  context.stop = false;

  // initialize file times
  for (auto index : range(context.filenames.size())) {
    auto time = std::filesystem::last_write_time(context.filenames[index])
                    .time_since_epoch()
                    .count();
    context.filetimes[index] = (int64_t)time;
  }

  // start watcher
  context.worker = std::async(std::launch::async, [&]() {
    // until done
    while (!context.stop) {
      // sleep
      std::this_thread::sleep_for(std::chrono::milliseconds(context.delay));

      // check times
      auto changed = false;
      for (auto index : range(context.filenames.size())) {
        auto time = std::filesystem::last_write_time(context.filenames[index])
                        .time_since_epoch()
                        .count();
        if ((int64_t)time != context.filetimes[index]) {
          changed                  = true;
          context.filetimes[index] = (int64_t)time;
        }
      }

      // update version
      if (changed) context.version++;
    }
  });
}

// Stop file watcher
void watch_stop(watch_context& context) {
  context.stop = true;
  if (context.worker.valid()) context.worker.get();
}

// Got file versions
int get_version(const watch_context& context) { return context.version; }

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

// Load a binary file
string load_text(const string& filename) {
  auto str = string{};
  load_text(filename, str);
  return str;
}

// Load a text file
void load_text(const string& filename, string& str) {
  // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
  auto fs = fopen_utf8(filename.c_str(), "rb");
  if (fs == nullptr) throw io_error("cannot open " + filename);
  fseek(fs, 0, SEEK_END);
  auto length = ftell(fs);
  fseek(fs, 0, SEEK_SET);
  str = string(length, '\0');
  if (fread(str.data(), 1, length, fs) != length) {
    fclose(fs);
    throw io_error("cannot read " + filename);
  }
  fclose(fs);
}

// Save a text file
void save_text(const string& filename, const string& str) {
  auto fs = fopen_utf8(filename.c_str(), "wt");
  if (fs == nullptr) throw io_error("cannot create " + filename);
  if (fprintf(fs, "%s", str.c_str()) < 0) {
    fclose(fs);
    throw io_error("cannot write " + filename);
  }
  fclose(fs);
}

// Load a binary file
vector<byte> load_binary(const string& filename) {
  auto data = vector<byte>{};
  load_binary(filename, data);
  return data;
}

// Load a binary file
void load_binary(const string& filename, vector<byte>& data) {
  // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
  auto fs = fopen_utf8(filename.c_str(), "rb");
  if (fs == nullptr) throw io_error("cannot open " + filename);
  fseek(fs, 0, SEEK_END);
  auto length = ftell(fs);
  fseek(fs, 0, SEEK_SET);
  data = vector<byte>(length);
  if (fread(data.data(), 1, length, fs) != length) {
    fclose(fs);
    throw io_error("cannot read " + filename);
  }
  fclose(fs);
}

// Save a binary file
void save_binary(const string& filename, const vector<byte>& data) {
  auto fs = fopen_utf8(filename.c_str(), "wb");
  if (fs == nullptr) throw io_error("cannot create " + filename);
  if (fwrite(data.data(), 1, data.size(), fs) != data.size()) {
    fclose(fs);
    throw io_error("cannot write " + filename);
  }
  fclose(fs);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON SUPPORT
// -----------------------------------------------------------------------------
namespace yocto {

// Json values
using json_value = nlohmann::ordered_json;

// Load/save json
[[maybe_unused]] static void load_json(
    const string& filename, json_value& json);
[[maybe_unused]] static json_value load_json(const string& filename) {
  auto json = json_value{};
  load_json(filename, json);
  return json;
}
[[maybe_unused]] static void load_json(
    const string& filename, json_value& json) {
  auto text = load_text(filename);
  try {
    json = json_value::parse(text);
  } catch (...) {
    throw io_error("cannot parse " + filename);
  }
}
[[maybe_unused]] static void save_json(
    const string& filename, const json_value& json) {
  return save_text(filename, json.dump(2));
}

// Json conversions
template <typename T, size_t N>
inline void to_json(json_value& json, const vec<T, N>& value) {
  nlohmann::to_json(json, (const array<T, N>&)value);
}
template <typename T, size_t N>
inline void to_json(json_value& json, const frame<T, N>& value) {
  nlohmann::to_json(json, (const array<array<T, N>, (N + 1)>&)value);
}
template <typename T, size_t N, size_t M>
inline void to_json(json_value& json, const mat<T, N, M>& value) {
  nlohmann::to_json(json, (const array<array<T, N>, M>&)value);
}
template <typename T, size_t N>
inline void from_json(const json_value& json, vec<T, N>& value) {
  nlohmann::from_json(json, (array<T, N>&)value);
}
template <typename T, size_t N>
inline void from_json(const json_value& json, frame<T, N>& value) {
  if (json.is_array() && !json.empty() && !json.front().is_array()) {
    nlohmann::from_json(json, (array<T, N*(N + 1)>&)value);
  } else {
    nlohmann::from_json(json, (array<array<T, N>, (N + 1)>&)value);
  }
}
template <typename T, size_t N, size_t M>
inline void from_json(const json_value& json, mat<T, N, M>& value) {
  if (json.is_array() && !json.empty() && !json.front().is_array()) {
    nlohmann::from_json(json, (array<T, N * M>&)value);
  } else {
    nlohmann::from_json(json, (array<array<T, N>, M>&)value);
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// MATH TYPE SUPPORT
// -----------------------------------------------------------------------------
namespace yocto {

static vec3f   to_math(const array<float, 3>& value) { return (vec3f&)value; }
static frame3f to_math(const array<float, 12>& value) {
  return (frame3f&)value;
}

static array<float, 3> to_array(const vec3f& value) {
  return (array<float, 3>&)value;
}
static array<float, 12> to_array(const frame3f& value) {
  return (array<float, 12>&)value;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Check if an image is HDR based on filename.
bool is_srgb_filename(const string& filename) {
  auto ext = path_extension(filename);
  if (ext == ".ypreset" || ext == ".YPRESET") return is_srgb_preset(filename);
  return ext == ".png" || ext == ".jpg" || ext == ".jpeg" || ext == ".bmp" ||
         ext == ".tga";
}

// Loads a float image.
array2d<vec4f> load_image(const string& filename, bool srgb) {
  auto image = array2d<vec4f>{};
  load_image(filename, image, srgb);
  return image;
}

// Loads a byte image.
array2d<vec4b> load_imageb(const string& filename, bool srgb) {
  auto image = array2d<vec4b>{};
  load_image(filename, image, srgb);
  return image;
}

// Loads a float image.
void load_image(const string& filename, array2d<vec4f>& image, bool srgb) {
  auto ext = path_extension(filename);
  if (ext == ".exr" || ext == ".EXR") {
    auto buffer = load_binary(filename);
    auto width = 0, height = 0;
    auto pixels = (float*)nullptr;
    if (LoadEXRFromMemory(&pixels, &width, &height, buffer.data(),
            buffer.size(), nullptr) != 0)
      throw io_error{"cannot read " + filename};
    image = {(vec4f*)pixels, {(size_t)width, (size_t)height}};
    if (srgb) image = rgb_to_srgb(image);
    free(pixels);
  } else if (ext == ".hdr" || ext == ".HDR") {
    auto buffer = load_binary(filename);
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_loadf_from_memory(
        buffer.data(), (int)buffer.size(), &width, &height, &ncomp, 4);
    if (pixels == nullptr) throw io_error{"cannot read " + filename};
    image = {(vec4f*)pixels, {(size_t)width, (size_t)height}};
    if (srgb) image = rgb_to_srgb(image);
    free(pixels);
  } else if (ext == ".png" || ext == ".PNG" || ext == ".jpg" || ext == ".JPG" ||
             ext == ".jpeg" || ext == ".JPEG" || ext == ".tga" ||
             ext == ".TGA" || ext == ".bmp" || ext == ".BMP") {
    auto imageb = load_imageb(filename);
    image       = srgb ? byte_to_float(imageb) : srgbb_to_rgb(imageb);
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    image = make_image_preset(filename);
    if (is_srgb_preset(filename) != srgb) {
      image = srgb ? rgb_to_srgb(image) : srgb_to_rgb(image);
    }
  } else {
    throw io_error{"unsupported format " + filename};
  }
}

// Loads a byte image.
void load_image(const string& filename, array2d<vec4b>& image, bool srgb) {
  auto ext = path_extension(filename);
  if (ext == ".exr" || ext == ".EXR" || ext == ".hdr" || ext == ".HDR") {
    auto imagef = load_image(filename);
    image       = srgb ? rgb_to_srgbb(imagef) : float_to_byte(imagef);
  } else if (ext == ".png" || ext == ".PNG" || ext == ".jpg" || ext == ".JPG" ||
             ext == ".jpeg" || ext == ".JPEG" || ext == ".tga" ||
             ext == ".TGA" || ext == ".bmp" || ext == ".BMP") {
    auto buffer = load_binary(filename);
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_load_from_memory(
        buffer.data(), (int)buffer.size(), &width, &height, &ncomp, 4);
    if (pixels == nullptr) throw io_error{"cannot read " + filename};
    image = {(vec4b*)pixels, {(size_t)width, (size_t)height}};
    if (!srgb) image = float_to_byte(srgbb_to_rgb(image));
    free(pixels);
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    auto imagef = load_image(filename);
    if (is_srgb_preset(filename) != srgb) {
      image = srgb ? rgb_to_srgbb(imagef) : float_to_byte(srgb_to_rgb(imagef));
    } else {
      image = float_to_byte(imagef);
    }
  } else {
    throw io_error{"unsupported format " + filename};
  }
}

// Saves a float image.
void save_image(
    const string& filename, const array2d<vec4f>& image, bool srgb) {
  // write data
  auto stbi_write_data = [](void* context, void* data, int size) {
    auto& buffer = *(vector<byte>*)context;
    buffer.insert(buffer.end(), (byte*)data, (byte*)data + size);
  };

  // grab data for low level apis
  auto [width, height] = (vec2i)image.extents();
  auto num_channels    = 4;

  auto ext = path_extension(filename);
  if (ext == ".hdr" || ext == ".HDR") {
    auto& image_ = srgb ? srgb_to_rgb(image) : image;
    auto  buffer = vector<byte>{};
    if (!(bool)stbi_write_hdr_to_func(stbi_write_data, &buffer, width, height,
            num_channels, (const float*)image_.data()))
      throw io_error{"cannot write " + filename};
    return save_binary(filename, buffer);
  } else if (ext == ".exr" || ext == ".EXR") {
    auto& image_ = srgb ? srgb_to_rgb(image) : image;
    auto  data   = (byte*)nullptr;
    auto  count  = (size_t)0;
    if (SaveEXRToMemory((const float*)image_.data(), width, height,
            num_channels, 1, &data, &count, nullptr) < 0)
      throw io_error{"cannot write " + filename};
    auto buffer = vector<byte>{data, data + count};
    free(data);
    return save_binary(filename, buffer);
  } else if (ext == ".png" || ext == ".PNG" || ext == ".jpg" || ext == ".JPG" ||
             ext == ".jpeg" || ext == ".JPEG" || ext == ".tga" ||
             ext == ".TGA" || ext == ".bmp" || ext == ".BMP") {
    save_image(filename, srgb ? float_to_byte(image) : rgb_to_srgbb(image));
  } else {
    throw io_error{"unsupported format " + filename};
  }
}

// Saves a byte image.
void save_image(
    const string& filename, const array2d<vec4b>& image, bool srgb) {
  // write data
  auto stbi_write_data = [](void* context, void* data, int size) {
    auto& buffer = *(vector<byte>*)context;
    buffer.insert(buffer.end(), (byte*)data, (byte*)data + size);
  };

  // grab data for low level apis
  auto [width, height] = (vec2i)image.extents();
  auto num_channels    = 4;

  auto ext = path_extension(filename);
  if (ext == ".hdr" || ext == ".HDR" || ext == ".exr" || ext == ".EXR") {
    return save_image(
        filename, srgb ? srgbb_to_rgb(image) : byte_to_float(image));
  } else if (ext == ".png" || ext == ".PNG") {
    auto& image_ = srgb ? image : rgb_to_srgbb(byte_to_float(image));
    auto  buffer = vector<byte>{};
    if (!(bool)stbi_write_png_to_func(stbi_write_data, &buffer, width, height,
            num_channels, (const byte*)image_.data(), width * 4))
      throw io_error{"cannot write " + filename};
    return save_binary(filename, buffer);
  } else if (ext == ".jpg" || ext == ".JPG" || ext == ".jpeg" ||
             ext == ".JPEG") {
    auto& image_ = srgb ? image : rgb_to_srgbb(byte_to_float(image));
    auto  buffer = vector<byte>{};
    if (!(bool)stbi_write_jpg_to_func(stbi_write_data, &buffer, width, height,
            num_channels, (const byte*)image_.data(), 75))
      throw io_error{"cannot write " + filename};
    return save_binary(filename, buffer);
  } else if (ext == ".tga" || ext == ".TGA") {
    auto& image_ = srgb ? image : rgb_to_srgbb(byte_to_float(image));
    auto  buffer = vector<byte>{};
    if (!(bool)stbi_write_tga_to_func(stbi_write_data, &buffer, width, height,
            num_channels, (const byte*)image_.data()))
      throw io_error{"cannot write " + filename};
    return save_binary(filename, buffer);
  } else if (ext == ".bmp" || ext == ".BMP") {
    auto& image_ = srgb ? image : rgb_to_srgbb(byte_to_float(image));
    auto  buffer = vector<byte>{};
    if (!(bool)stbi_write_bmp_to_func(stbi_write_data, &buffer, width, height,
            num_channels, (const byte*)image_.data()))
      throw io_error{"cannot write " + filename};
    return save_binary(filename, buffer);
  } else {
    throw io_error{"unsupported format " + filename};
  }
}

bool is_srgb_preset(const string& type_) {
  auto type = path_basename(type_);
  return type.find("sky") == string::npos;
}
array2d<vec4f> make_image_preset(const string& type_) {
  auto type    = path_basename(type_);
  auto extents = vec2s{1024, 1024};
  if (type.find("sky") != string::npos) extents = {2048, 1024};
  if (type.find("images2") != string::npos) extents = {2048, 1024};
  if (type == "grid") {
    return make_grid(extents);
  } else if (type == "checker") {
    return make_checker(extents);
  } else if (type == "bumps") {
    return make_bumps(extents);
  } else if (type == "uvramp") {
    return make_uvramp(extents);
  } else if (type == "gammaramp") {
    return make_gammaramp(extents);
  } else if (type == "uvgrid") {
    return make_uvgrid(extents);
  } else if (type == "colormapramp") {
    return make_colormapramp(extents);
  } else if (type == "sky") {
    return make_sunsky(
        extents, pif / 4, 3.0f, false, 1.0f, 1.0f, vec3f{0.7f, 0.7f, 0.7f});
  } else if (type == "sunsky") {
    return make_sunsky(
        extents, pif / 4, 3.0f, true, 1.0f, 1.0f, vec3f{0.7f, 0.7f, 0.7f});
  } else if (type == "gnoisemap") {
    return make_gnoisemap(extents, 1.0f);
  } else if (type == "vnoisemap") {
    return make_vnoisemap(extents, 1.0f);
  } else if (type == "fnoisemap") {
    return make_fnoisemap(extents, 1.0f);
  } else if (type == "rnoisemap") {
    return make_rnoisemap(extents, 1.0f);
  } else if (type == "tnoisemap") {
    return make_tnoisemap(extents, 1.0f);
  } else if (type == "bump-normal") {
    return make_bumps(extents);
    // TODO(fabio): fix color space
    // img   = srgb_to_rgb(bump_to_normal(img, 0.05f));
  } else if (type == "images1") {
    auto sub_types  = vector<string>{"grid", "uvgrid", "checker", "gammaramp",
         "bumps", "bump-normal", "noise", "fbm", "blackbodyramp"};
    auto sub_images = vector<array2d<vec4f>>();
    for (auto& sub_type : sub_types)
      sub_images.push_back(make_image_preset(sub_type));
    auto montage_size = vec2s{0, 0};
    for (auto& sub_image : sub_images) {
      montage_size = {montage_size.x + sub_image.extents().x,
          max(montage_size.y, sub_image.extents().y)};
    }
    auto image = array2d<vec4f>(montage_size);
    auto pos   = (size_t)0;
    for (auto& sub_image : sub_images) {
      set_region(image, sub_image, {pos, 0});
      pos += sub_image.extents().x;
    }
    return image;
  } else if (type == "images2") {
    auto sub_types  = vector<string>{"sky", "sunsky"};
    auto sub_images = vector<array2d<vec4f>>();
    for (auto& sub_type : sub_types)
      sub_images.push_back(make_image_preset(sub_type));
    auto montage_size = vec2s{0, 0};
    for (auto& sub_image : sub_images) {
      montage_size = {montage_size.x + sub_image.extents().x,
          max(montage_size.y, sub_image.extents().y)};
    }
    auto image = array2d<vec4f>(montage_size);
    auto pos   = (size_t)0;
    for (auto& sub_image : sub_images) {
      set_region(image, sub_image, {pos, 0});
      pos += sub_image.extents().x;
    }
    return image;
  } else if (type == "test-floor") {
    return add_border(make_grid(extents), 0.0025f);
  } else if (type == "test-grid") {
    return make_grid(extents);
  } else if (type == "test-checker") {
    return make_checker(extents);
  } else if (type == "test-bumps") {
    return make_bumps(extents);
  } else if (type == "test-uvramp") {
    return make_uvramp(extents);
  } else if (type == "test-gammaramp") {
    return make_gammaramp(extents);
  } else if (type == "test-colormapramp") {
    return make_colormapramp(extents);
    // TODO(fabio): fix color space
    // img   = srgb_to_rgb(img);
  } else if (type == "test-uvgrid") {
    return make_uvgrid(extents);
  } else if (type == "test-sky") {
    return make_sunsky(
        extents, pif / 4, 3.0f, false, 1.0f, 1.0f, vec3f{0.7f, 0.7f, 0.7f});
  } else if (type == "test-sunsky") {
    return make_sunsky(
        extents, pif / 4, 3.0f, true, 1.0f, 1.0f, vec3f{0.7f, 0.7f, 0.7f});
  } else if (type == "test-bumps-normal") {
    return bump_to_normal(make_bumps(extents), 0.05f);
  } else if (type == "test-bumps-displacement") {
    return make_bumps(extents);
    // TODO(fabio): fix color space
    // img   = srgb_to_rgb(img);
  } else if (type == "test-checker-opacity") {
    return make_checker(extents, 1.0f, vec4f{1, 1, 1, 1}, vec4f{0, 0, 0, 0});
  } else if (type == "test-grid-opacity") {
    return make_grid(extents, 1.0f, vec4f{1, 1, 1, 1}, vec4f{0, 0, 0, 0});
  } else {
    throw io_error{"unknown preset " + type};
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load mesh
void load_shape(const string& filename, shape_data& shape, bool flip_texcoord) {
  shape = {};

  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = load_ply(filename);
    // TODO: remove when all as arrays
    get_positions(ply, (vector<array<float, 3>>&)shape.positions);
    get_normals(ply, (vector<array<float, 3>>&)shape.normals);
    get_texcoords(
        ply, (vector<array<float, 2>>&)shape.texcoords, flip_texcoord);
    get_colors(ply, (vector<array<float, 4>>&)shape.colors);
    get_radius(ply, shape.radius);
    get_faces(ply, (vector<array<int, 3>>&)shape.triangles,
        (vector<array<int, 4>>&)shape.quads);
    get_lines(ply, (vector<array<int, 2>>&)shape.lines);
    get_points(ply, shape.points);
    if (shape.points.empty() && shape.lines.empty() &&
        shape.triangles.empty() && shape.quads.empty())
      throw io_error{"empty shape " + filename};
  } else if (ext == ".obj" || ext == ".OBJ") {
    auto obj       = load_sobj(filename, false);
    auto materials = vector<int>{};
    // TODO: remove when all as arrays
    get_positions(obj, (vector<array<float, 3>>&)shape.positions);
    get_normals(obj, (vector<array<float, 3>>&)shape.normals);
    get_texcoords(
        obj, (vector<array<float, 2>>&)shape.texcoords, flip_texcoord);
    get_faces(obj, (vector<array<int, 3>>&)shape.triangles,
        (vector<array<int, 4>>&)shape.quads, materials);
    get_lines(obj, (vector<array<int, 2>>&)shape.lines, materials);
    get_points(obj, shape.points, materials);
    if (shape.points.empty() && shape.lines.empty() &&
        shape.triangles.empty() && shape.quads.empty())
      throw io_error{"empty shape " + filename};
  } else if (ext == ".stl" || ext == ".STL") {
    auto stl = load_stl(filename, true);
    if (stl.shapes.size() != 1) throw io_error{"empty shape " + filename};
    auto fnormals = vector<vec3f>{};
    if (!get_triangles(stl, 0, (vector<array<int, 3>>&)shape.triangles,
            (vector<array<float, 3>>&)shape.positions,
            (vector<array<float, 3>>&)fnormals))
      throw io_error{"empty shape " + filename};
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    shape = make_shape_preset(filename);
  } else {
    throw io_error("unsupported format " + filename);
  }
}

// Save ply mesh
void save_shape(const string& filename, const shape_data& shape,
    bool flip_texcoord, bool ascii) {
  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = ply_model{};
    // TODO: remove when all as arrays
    add_positions(ply, (const vector<array<float, 3>>&)shape.positions);
    add_normals(ply, (const vector<array<float, 3>>&)shape.normals);
    add_texcoords(
        ply, (const vector<array<float, 2>>&)shape.texcoords, flip_texcoord);
    add_colors(ply, (const vector<array<float, 4>>&)shape.colors);
    add_radius(ply, shape.radius);
    add_faces(ply, (const vector<array<int, 3>>&)shape.triangles,
        (const vector<array<int, 4>>&)shape.quads);
    add_lines(ply, (const vector<array<int, 2>>&)shape.lines);
    add_points(ply, shape.points);
    save_ply(filename, ply);
  } else if (ext == ".obj" || ext == ".OBJ") {
    auto obj = obj_shape{};
    // TODO: remove when all as arrays
    add_positions(obj, (const vector<array<float, 3>>&)shape.positions);
    add_normals(obj, (const vector<array<float, 3>>&)shape.normals);
    add_texcoords(
        obj, (const vector<array<float, 2>>&)shape.texcoords, flip_texcoord);
    add_triangles(obj, (const vector<array<int, 3>>&)shape.triangles, 0,
        !shape.normals.empty(), !shape.texcoords.empty());
    add_quads(obj, (const vector<array<int, 4>>&)shape.quads, 0,
        !shape.normals.empty(), !shape.texcoords.empty());
    add_lines(obj, (const vector<array<int, 2>>&)shape.lines, 0,
        !shape.normals.empty(), !shape.texcoords.empty());
    add_points(
        obj, shape.points, 0, !shape.normals.empty(), !shape.texcoords.empty());
    save_obj(filename, obj);
  } else if (ext == ".stl" || ext == ".STL") {
    auto stl = stl_model{};
    if (is_lines(shape)) throw io_error{"empty shape " + filename};
    if (is_points(shape)) throw io_error{"empty shape " + filename};
    if (is_triangles(shape)) {
      add_triangles(stl, (const vector<array<int, 3>>&)shape.triangles,
          (const vector<array<float, 3>>&)shape.positions, {});
    } else if (is_quads(shape)) {
      auto triangles = quads_to_triangles(shape.quads);
      add_triangles(stl, (const vector<array<int, 3>>&)triangles,
          (const vector<array<float, 3>>&)shape.positions, {});
    } else {
      throw io_error{"empty shape " + filename};
    }
    save_stl(filename, stl);
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
    save_text(filename, str);
  } else {
    throw io_error("unsupported format " + filename);
  }
}

// Load face-varying mesh
void load_fvshape(
    const string& filename, fvshape_data& shape, bool flip_texcoord) {
  shape = {};

  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = load_ply(filename);
    // TODO: remove when all as arrays
    get_positions(ply, (vector<array<float, 3>>&)shape.positions);
    get_normals(ply, (vector<array<float, 3>>&)shape.normals);
    get_texcoords(
        ply, (vector<array<float, 2>>&)shape.texcoords, flip_texcoord);
    get_quads(ply, (vector<array<int, 4>>&)shape.quadspos);
    if (!shape.normals.empty()) shape.quadsnorm = shape.quadspos;
    if (!shape.texcoords.empty()) shape.quadstexcoord = shape.quadspos;
    if (shape.quadspos.empty()) throw io_error{"empty shape " + filename};
  } else if (ext == ".obj" || ext == ".OBJ") {
    auto obj = load_sobj(filename, true);
    // TODO: remove when all as arrays
    auto materials = vector<int>{};
    get_positions(obj, (vector<array<float, 3>>&)shape.positions);
    get_normals(obj, (vector<array<float, 3>>&)shape.normals);
    get_texcoords(
        obj, (vector<array<float, 2>>&)shape.texcoords, flip_texcoord);
    get_fvquads(obj, (vector<array<int, 4>>&)shape.quadspos,
        (vector<array<int, 4>>&)shape.quadsnorm,
        (vector<array<int, 4>>&)shape.quadstexcoord, materials);
    if (shape.quadspos.empty()) throw io_error{"empty shape " + filename};
  } else if (ext == ".stl" || ext == ".STL") {
    auto stl = load_stl(filename, true);
    if (stl.shapes.empty()) throw io_error{"empty shape " + filename};
    if (stl.shapes.size() > 1) throw io_error{"empty shape " + filename};
    auto fnormals  = vector<vec3f>{};
    auto triangles = vector<vec3i>{};
    if (!get_triangles(stl, 0, (vector<array<int, 3>>&)triangles,
            (vector<array<float, 3>>&)shape.positions,
            (vector<array<float, 3>>&)fnormals))
      throw io_error{"empty shape " + filename};
    shape.quadspos = triangles_to_quads(triangles);
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    shape = make_fvshape_preset(filename);
  } else {
    throw io_error("unsupported format " + filename);
  }
}

// Save ply mesh
void save_fvshape(const string& filename, const fvshape_data& shape,
    bool flip_texcoord, bool ascii) {
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
    // TODO: remove when all as arrays
    add_positions(ply, (const vector<array<float, 3>>&)split_positions);
    add_normals(ply, (const vector<array<float, 3>>&)split_normals);
    add_texcoords(
        ply, (const vector<array<float, 2>>&)split_texcoords, flip_texcoord);
    add_quads(ply, (const vector<array<int, 4>>&)split_quads);
    save_ply(filename, ply);
  } else if (ext == ".obj" || ext == ".OBJ") {
    auto obj = obj_shape{};
    // TODO: remove when all as arrays
    add_positions(obj, (const vector<array<float, 3>>&)shape.positions);
    add_normals(obj, (const vector<array<float, 3>>&)shape.normals);
    add_texcoords(
        obj, (const vector<array<float, 2>>&)shape.texcoords, flip_texcoord);
    add_fvquads(obj, (const vector<array<int, 4>>&)shape.quadspos,
        (const vector<array<int, 4>>&)shape.quadsnorm,
        (const vector<array<int, 4>>&)shape.quadstexcoord, 0);
    save_obj(filename, obj);
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
      auto triangles = quads_to_triangles(split_quads);
      add_triangles(stl, (const vector<array<int, 3>>&)triangles,
          (const vector<array<float, 3>>&)split_positions, {});
    } else {
      throw io_error{"empty shape " + filename};
    }
    save_stl(filename, stl);
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
    save_text(filename, str);
  } else {
    throw io_error("unsupported format " + filename);
  }
}

// Shape presets used for testing.
shape_data make_shape_preset(const string& type_) {
  auto type       = path_basename(type_);
  auto test_xform = translation_frame(vec3f{0, 0.75f, 0}) *
                    scaling_frame(0.75f);
  auto test_radius_scale = 0.75f;
  if (type == "quad") {
    return make_quad();
  } else if (type == "quady") {
    return make_quady();
  } else if (type == "cube") {
    return make_cube();
  } else if (type == "rounded_cube") {
    return make_rounded_box();
  } else if (type == "sphere") {
    return make_sphere();
  } else if (type == "matcube") {
    return make_rounded_box();
  } else if (type == "matsphere") {
    return make_uvspherey();
  } else if (type == "disk") {
    return make_disk();
  } else if (type == "bulged_disk") {
    return make_bulged_disk();
  } else if (type == "bulged_quad") {
    return make_bulged_rect();
  } else if (type == "rect") {
    return make_rect();
  } else if (type == "recty") {
    return make_recty();
  } else if (type == "box") {
    return make_box();
  } else if (type == "rounded_box") {
    return make_rounded_box();
  } else if (type == "tsphere") {
    return make_sphere();
  } else if (type == "uvsphere") {
    return make_uvsphere();
  } else if (type == "capped_uvsphere") {
    return make_capped_uvsphere();
  } else if (type == "uvspherey") {
    return make_uvspherey();
  } else if (type == "capped_uvspherey") {
    return make_capped_uvspherey();
  } else if (type == "uvdisk") {
    return make_uvdisk();
  } else if (type == "uvcylinder") {
    return make_uvcylinder();
  } else if (type == "rounded_uvcylinder") {
    return make_rounded_uvcylinder({32, 32, 32});
  } else if (type == "uvcapsule") {
    return make_uvcapsule();
  } else if (type == "uvcone") {
    return make_uvcone();
  } else if (type == "geosphere") {
    return make_geosphere();
  } else if (type == "floor") {
    return make_floor();
  } else if (type == "bent_floor") {
    return make_bent_floor();
  } else if (type == "matball") {
    return make_sphere();
  } else if (type == "lineball") {
    auto base = transform_shape(make_sphere(), scaling_frame(0.8f));
    return make_random_lines(base, pow2(16), 4, {0.2f, 0.2f}, {0.002f, 0.001f});
  } else if (type == "hairball") {
    auto base = transform_shape(make_sphere(), scaling_frame(0.8f));
    return make_random_hairs(base, pow2(16), 4, {0.2f, 0.2f}, {0.002f, 0.001f});
  } else if (type == "hairball_interior") {
    return transform_shape(make_sphere(), scaling_frame(0.8f));
  } else if (type == "monkey") {
    return make_monkey();
  } else if (type == "wtcube") {
    return make_wtcube();
  } else if (type == "opcube") {
    return make_opcube();
  } else if (type == "wtbox") {
    return make_wtbox();
  } else if (type == "fvcube") {
    return fvshape_to_shape(make_fvbox());
  } else if (type == "fvsphere") {
    return fvshape_to_shape(make_fvsphere());
  } else if (type == "monkey_subdiv") {
    return add_normals(subdivide_shape(make_monkey(0), 2, true));
  } else if (type == "wtcube_subdiv") {
    return add_normals(subdivide_shape(make_wtcube(), 4, true));
  } else if (type == "opcube_subdiv") {
    return add_normals(subdivide_shape(make_opcube(), 4, true));
  } else if (type == "wtbox_subdiv") {
    return add_normals(subdivide_shape(make_wtbox(), 4, true));
  } else if (type == "opbox_subdiv") {
    return add_normals(subdivide_shape(make_opbox(), 4, true));
  } else if (type == "fvcube_subdiv") {
    return fvshape_to_shape(
        add_normals(subdivide_fvshape(make_fvcube(), 4, true)));
  } else if (type == "displaced_quady") {
    return make_recty({pow2(8), pow2(8)});
  } else if (type == "displaced_sphere") {
    return make_sphere(pow2(7));
  } else if (type == "floor") {
    return make_floor();
  } else if (type == "test_cube") {
    return transform_shape(make_rounded_box(), test_xform);
  } else if (type == "test_uvsphere") {
    return transform_shape(make_uvsphere(), test_xform);
  } else if (type == "test_capped_uvsphere") {
    return transform_shape(make_capped_uvsphere(), test_xform);
  } else if (type == "test_uvspherey") {
    return transform_shape(make_uvspherey(), test_xform);
  } else if (type == "test_capped_uvspherey") {
    return transform_shape(make_capped_uvspherey(), test_xform);
  } else if (type == "test_sphere") {
    return transform_shape(make_sphere(), test_xform);
  } else if (type == "test_matcube") {
    return transform_shape(make_rounded_box(), test_xform);
  } else if (type == "test_matsphere") {
    return transform_shape(make_uvspherey(), test_xform);
  } else if (type == "test_displaced_sphere") {
    return transform_shape(make_sphere(pow2(7)), test_xform);
  } else if (type == "test_smallsphere") {
    return transform_shape(make_sphere(),
        translation_frame(vec3f{0, 0.015f, 0}) * scaling_frame(0.015f));
  } else if (type == "test_disk") {
    return transform_shape(make_disk(), test_xform);
  } else if (type == "test_uvcylinder") {
    return transform_shape(make_rounded_uvcylinder(), test_xform);
  } else if (type == "test_floor") {
    return make_floor();
  } else if (type == "test_quad") {
    return transform_shape(make_quad(), test_xform);
  } else if (type == "test_quady") {
    return transform_shape(make_quady(), test_xform);
  } else if (type == "test_displaced_quad") {
    return transform_shape(make_rect({pow2(8), pow2(8)}), test_xform);
  } else if (type == "test_displaced_quady") {
    return transform_shape(make_recty({pow2(8), pow2(8)}), test_xform);
  } else if (type == "test_matball") {
    return transform_shape(make_sphere(), test_xform);
  } else if (type == "test_geosphere") {
    return transform_shape(make_geosphere(), test_xform);
  } else if (type == "test_faceted_geosphere") {
    return transform_shape(remove_normals(make_geosphere()), test_xform);
  } else if (type == "test_subdivided_geosphere") {
    return transform_shape(make_geosphere(6), test_xform);
  } else if (type == "test_lineball") {
    auto base = transform_shape(
        make_sphere(), test_xform * scaling_frame(0.8f));
    return make_random_lines(base, pow2(16), 4, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f});
  } else if (type == "test_hairball") {
    auto base = transform_shape(
        make_sphere(), test_xform * scaling_frame(0.8f));
    return make_random_hairs(base, pow2(16), 4, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f});
  } else if (type == "test_hairball-interior") {
    return transform_shape(make_sphere(), test_xform * scaling_frame(0.8f));
  } else if (type == "test_suzanne_subdiv") {
    return transform_shape(make_monkey(0), test_xform * scaling_frame(0.8f));
  } else if (type == "test_wtcube") {
    return transform_shape(make_wtcube(), test_xform);
  } else if (type == "test_arealight1") {
    return transform_shape(make_quad(), scaling_frame(0.2f));
  } else if (type == "test_arealight2") {
    return transform_shape(make_quad(), scaling_frame(0.2f));
  } else if (type == "test_largearealight1") {
    return transform_shape(make_quad(), scaling_frame(0.4f));
  } else if (type == "test_largearealight2") {
    return transform_shape(make_quad(), scaling_frame(0.4f));
  } else if (type == "test_pointlight1") {
    return make_point(0);
  } else if (type == "test_pointlight2") {
    return make_point(0);
  } else if (type == "test_point") {
    return make_point();
  } else if (type == "test_points") {
    return make_random_points(4096);
  } else if (type == "test_points_random") {
    return transform_shape(make_random_points(4096), test_xform);
  } else if (type == "test_points_grid") {
    return transform_shape(make_point_grid(), test_xform, test_radius_scale);
  } else if (type == "test_lines_grid") {
    return transform_shape(make_lines(), test_xform, test_radius_scale);
  } else if (type == "test_thickpoints_grid") {
    return transform_shape(
        make_point_grid(), test_xform, test_radius_scale * 10);
  } else if (type == "test_thicklines_grid") {
    return transform_shape(make_lines(), test_xform, test_radius_scale);
  } else if (type == "test_particles") {
    return make_points(4096);
  } else if (type == "test_cloth") {
    return transform_shape(make_rect({pow2(6), pow2(6)}), scaling_frame(0.2f));
  } else if (type == "test_clothy") {
    return transform_shape(make_recty({pow2(6), pow2(6)}), scaling_frame(0.2f));
  } else {
    throw io_error{"unknown preset " + type};
  }
}

// Shape presets used for testing.
fvshape_data make_fvshape_preset(const string& type) {
  auto test_xform = translation_frame(vec3f{0, 0.75f, 0}) *
                    scaling_frame(0.75f);
  if (type == "fvcube") {
    return make_fvbox();
  } else if (type == "fvsphere") {
    return make_fvsphere();
  } else if (type == "test_facevarying_cube") {
    return transform_fvshape(make_fvbox(), test_xform);
  } else if (type == "test_facevarying_sphere") {
    return transform_fvshape(make_fvsphere(), test_xform);
  } else {
    return shape_to_fvshape(make_shape_preset(type));
  }
}

// Load mesh
shape_data load_shape(const string& filename, bool flip_texcoord) {
  auto shape = shape_data{};
  load_shape(filename, shape, flip_texcoord);
  return shape;
}

// Load mesh
fvshape_data load_fvshape(const string& filename, bool flip_texcoord) {
  auto shape = fvshape_data{};
  load_fvshape(filename, shape, flip_texcoord);
  return shape;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// TEXTURE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Loads/saves an image. Chooses hdr or ldr based on file name.
void load_texture(const string& filename, texture_data& texture) {
  auto ext = path_extension(filename);
  if (ext == ".exr" || ext == ".EXR" || ext == ".hdr" || ext == ".HDR") {
    texture.pixelsf = load_image(filename);
  } else if (ext == ".png" || ext == ".PNG" || ext == ".jpg" || ext == ".JPG" ||
             ext == ".jpeg" || ext == ".JPEG" || ext == ".tga" ||
             ext == ".TGA" || ext == ".bmp" || ext == ".BMP") {
    texture.pixelsb = load_imageb(filename);
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    texture = make_texture_preset(filename);
  } else {
    throw io_error("unsupported format " + filename);
  }
}

// Saves an hdr image.
void save_texture(const string& filename, const texture_data& texture) {
  if (!texture.pixelsf.empty()) {
    save_image(filename, texture.pixelsf);
  } else {
    save_image(filename, texture.pixelsb);
  }
}

texture_data make_texture_preset(const string& type) {
  return image_to_texture(make_image_preset(type), !is_srgb_preset(type));
}

// Loads/saves an image. Chooses hdr or ldr based on file name.
texture_data load_texture(const string& filename) {
  auto texture = texture_data{};
  load_texture(filename, texture);
  return texture;
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
  if (scene.shape_names.empty() || scene.shape_names[idx].empty())
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
void add_missing_material(scene_data& scene, bool glossy = false,
    const array2d<vec4f>& texture = {}) {
  auto default_material = invalidid;
  for (auto& instance : scene.instances) {
    if (instance.material != invalidid) continue;
    if (default_material == invalidid) {
      auto& material   = scene.materials.emplace_back();
      default_material = (int)scene.materials.size() - 1;
      material.type    = glossy ? material_type::glossy : material_type::matte;
      material.color   = texture.empty() ? vec3f{0.8f, 0.8f, 0.8f}
                                         : vec3f{1, 1, 1};
      material.roughness = 0.1f;
      if (!texture.empty()) {
        scene.textures.push_back(image_to_texture(texture, false));
        material.color_tex = (int)scene.textures.size() - 1;
      }
    }
    instance.material = default_material;
  }
}

// Add missing cameras.
void add_missing_lights(scene_data& scene) {
  if (has_lights(scene)) return;
  add_sky(scene);
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

// Make a test scene with a callback to add your own objects.
template <typename Func>
scene_data make_test_scene(Func&& func) {
  auto scene = scene_data{};

  add_camera(scene, "camera", {0, 10, 35}, {0, 0, 0}, {0, 1, 0}, 0.1, 3, 0);

  add_environment(scene, "sky", identity3x4f, {0.5, 0.5, 0.5}, make_sunsky());

  add_instance(scene, "arealight1",
      lookat_frame(vec3f{-10, 20, 20}, {0, 0, 0}, {0, 1, 0}, true),
      add_shape(scene, "arealight1", scale_shape(make_rect(), 5)),
      add_emission_material(scene, "arealight1", {10, 10, 10}));
  add_instance(scene, "arealight2",
      lookat_frame(vec3f{+10, 20, 20}, {0, 0, 0}, {0, 1, 0}, true),
      add_shape(scene, "arealight2", scale_shape(make_rect(), 5)),
      add_emission_material(scene, "arealight2", {10, 10, 10}));

  add_instance(scene, "floor", translation_frame(vec3f{0, -1, 0}), make_floor(),
      make_matte_material(), make_grid());

  auto num = 5;
  for (auto idx : range(num)) {
    func(scene, "object" + std::to_string(idx + 1),
        translation_frame(vec3f{2.5f * (idx - num / 2), 0.0f, 0.0f}), idx);
  }

  return scene;
}

// Scene test 1
scene_data make_features1_scene() {
  return make_test_scene(
      [](scene_data& scene, const string& name, const frame3f& frame, int idx) {
        if (idx == 0) {
          add_instance(scene, name, frame, make_sphere(),
              make_glossy_material(), make_uvgrid());
        } else if (idx == 1) {
          add_instance(scene, name, frame, make_sphere(),
              make_refractive_material({1.0, 0.5, 0.5}));
        } else if (idx == 2) {
          add_instance(scene, name, frame, make_sphere(),
              make_scattering_material({0.5, 0.5, 0.5}, {0.3, 0.6, 0.3}),
              make_uvgrid());
        } else if (idx == 3) {
          add_instance(scene, name, frame, make_sphere(),
              make_glossy_material({0.5, 0.7, 0.5}), {}, {},
              bump_to_normal(make_bumps(), 0.05));
        } else if (idx == 4) {
          add_instance(scene, name, frame, make_sphere(),
              make_reflective_material({0.66, 0.45, 0.34}, 0.2f));
        } else {
          throw std::out_of_range("unknown instance");
        }
      });
}

// Scene test 2
scene_data make_features2_scene() {
  return make_test_scene([](scene_data& scene, const string& name,
                             const frame3f& frame, int idx) {
    if (idx == 0) {
      add_instance(scene, name, frame, make_sphere(), make_glossy_material(),
          make_uvgrid());
    } else if (idx == 1) {
      add_instance(scene, name, frame, make_monkey(),
          make_glossy_material({1.0, 0.5, 0.5}));
      add_subdiv(scene, name, make_monkey(0), (int)scene.shapes.size() - 1, 2);
    } else if (idx == 2) {
      add_instance(scene, name, frame, scale_shape(make_sphere(), 0.65f),
          make_matte_material());
      add_instance(scene, name, frame, make_random_hairs(scene.shapes.back()),
          make_matte_material());
    } else if (idx == 3) {
      add_instance(scene, name, frame, make_sphere(128),
          make_glossy_material({0.5, 0.7, 0.5}));
      add_subdiv(scene, name, make_sphere(128), (int)scene.shapes.size() - 1, 0,
          0.25f, make_bumps());
    } else if (idx == 4) {
      add_instance(scene, name, frame, make_rounded_box(),
          make_glossy_material(), make_uvgrid());
    } else {
      throw std::out_of_range("unknown instance");
    }
  });
}

// Scene test materials 1
scene_data make_materials1_scene() {
  return make_test_scene([shape = invalidid](scene_data& scene,
                             const string& name, const frame3f& frame,
                             int idx) mutable {
    if (shape == invalidid) {
      shape = add_shape(scene, "sphere", make_sphere());
    }

    if (idx == 0) {
      add_instance(scene, name, frame, shape,
          add_glossy_material(scene, "glossys", {0.5, 0.5, 0.7}, 0));
    } else if (idx == 1) {
      add_instance(scene, name, frame, shape,
          add_glossy_material(scene, "glossyr", {0.5, 0.7, 0.5}, 0.2));
    } else if (idx == 2) {
      add_instance(scene, name, frame, shape,
          add_glossy_material(scene, "matte", {0.7, 0.7, 0.7}));
    } else if (idx == 3) {
      add_instance(scene, name, frame, shape,
          add_reflective_material(scene, "reflectives", {0.7, 0.7, 0.7}, 0));
    } else if (idx == 4) {
      add_instance(scene, name, frame, shape,
          add_reflective_material(
              scene, "reflectiver", {0.66, 0.45, 0.34}, 0.2));
    } else {
      throw std::out_of_range("unknown instance");
    }
  });
}

// Scene test materials 2
scene_data make_materials2_scene() {
  return make_test_scene([shape = invalidid](scene_data& scene,
                             const string& name, const frame3f& frame,
                             int idx) mutable {
    if (shape == invalidid) {
      shape = add_shape(scene, "sphere", make_sphere());
    }

    if (idx == 0) {
      add_instance(scene, name, frame, shape,
          add_refractive_material(scene, "refractives", {1.0, 1.0, 1.0}, 0));
    } else if (idx == 1) {
      add_instance(scene, name, frame, shape,
          add_refractive_material(scene, "refractiver", {1.0, 0.7, 0.7}, 0.1));
    } else if (idx == 2) {
      add_instance(scene, name, frame, shape,
          add_opacity_material(scene, "opacity", {0.7, 0.5, 0.5}, 0.2));
    } else if (idx == 3) {
      add_instance(scene, name, frame, shape,
          add_transparent_material(scene, "transparents", {1.0, 1.0, 1.0}, 0));
    } else if (idx == 4) {
      add_instance(scene, name, frame, shape,
          add_transparent_material(
              scene, "transparentr", {1.0, 0.7, 0.7}, 0.1));
    } else {
      throw std::out_of_range("unknown instance");
    }
  });
}

// Scene test shapes 1
scene_data make_shapes1_scene() {
  return make_test_scene(
      [material = invalidid](scene_data& scene, const string& name,
          const frame3f& frame, int idx) mutable {
        if (material == invalidid) {
          auto texture = add_texture(scene, "uvgrid", make_uvgrid());
          material     = add_material(
              scene, "material", make_glossy_material({1, 1, 1}, 0.2, texture));
        }

        if (idx == 0) {
          add_instance(scene, name, frame,
              add_shape(scene, "sphere", make_sphere()), material);
        } else if (idx == 1) {
          add_instance(scene, name, frame,
              add_shape(scene, "capped_uvspherey",
                  flipyz_shape(make_capped_uvsphere())),
              material);
        } else if (idx == 2) {
          add_instance(scene, name, frame,
              add_shape(scene, "disk", make_disk()), material);
        } else if (idx == 3) {
          add_instance(scene, name, frame,
              add_shape(scene, "uvcylindery", flipyz_shape(make_uvcylinder())),
              material);
        } else if (idx == 4) {
          add_instance(scene, name, frame,
              add_shape(scene, "rounded_box", make_rounded_box()), material);
        } else {
          throw std::out_of_range("unknown instance");
        }
      });
}

// Scene test shapes 2
scene_data make_shapes2_scene() {
  return make_test_scene(
      [material = invalidid](scene_data& scene, const string& name,
          const frame3f& frame, int idx) mutable {
        if (material == invalidid) {
          material = add_material(
              scene, "material", make_glossy_material({0.5, 0.5, 1}, 0.2));
        }

        if (idx == 0) {
          add_instance(scene, name, frame,
              add_shape(scene, "sdcube", make_sdcube()), material);
        } else if (idx == 1) {
          add_instance(scene, name, frame,
              add_shape(scene, "sdmonkey", make_monkey()), material);
        } else if (idx == 2) {
          add_instance(scene, name, frame,
              add_shape(scene, "displaced",
                  displace_shape(make_sphere(128), make_bumps(), 0.1f)),
              material);
        } else if (idx == 3) {
          add_instance(scene, name, frame,
              add_shape(scene, "bunny", make_sphere()), material);
        } else if (idx == 4) {
          add_instance(scene, name, frame,
              add_shape(scene, "teapot", make_sphere()), material);
        } else {
          throw std::out_of_range("unknown instance");
        }
      });
}

// Scene test materials 4
scene_data make_materials4_scene() {
  return make_test_scene(
      [](scene_data& scene, const string& name, const frame3f& frame, int idx) {
        if (idx == 0) {
          add_instance(scene, name, frame, make_sphere(),
              make_volumetric_material({0.5, 0.5, 0.5}, {0.9, 0.9, 0.9}));
        } else if (idx == 1) {
          add_instance(scene, name, frame, make_sphere(),
              make_refractive_material({1.0, 0.5, 0.5}, 0));
        } else if (idx == 2) {
          add_instance(scene, name, frame, make_sphere(),
              make_refractive_material({1.0, 1.0, 1.0}, 0));
        } else if (idx == 3) {
          add_instance(scene, name, frame, make_sphere(),
              make_scattering_material({0.5, 0.5, 0.5}, {0.3, 0.6, 0.3}, 0));
        } else if (idx == 4) {
          add_instance(scene, name, frame, make_sphere(),
              make_volumetric_material({0.65, 0.65, 0.65}, {0.2, 0.2, 0.2}));
        } else {
          throw std::out_of_range("unknown instance");
        }
      });
}

// Scene test
scene_data make_test(const test_params& params) {
  // scene
  auto scene = scene_data{};
  // cameras
  switch (params.cameras) {
    case test_cameras_type::standard: {
      add_camera(scene, "default", {-0.75, 0.4, 0.9}, {-0.075, 0.05, -0.05},
          {0, 1, 0}, 0.05, 2.4, 0);
    } break;
    // TODO(fabio): fix wide camera
    case test_cameras_type::wide: {
      add_camera(scene, "default", {-0.75, 0.4, 0.9}, {-0.075, 0.05, -0.05},
          {0, 1, 0}, 0.05, 2.4, 0);
    } break;
  }
  // TODO(fabio): port other cameras
  switch (params.environments) {
    case test_environments_type::none: break;
    case test_environments_type::sky: {
      add_environment(scene, "sky",
          {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}}, {0.5, 0.5, 0.5},
          add_texture(scene, "sky",
              make_sunsky({2048, 1024}, pif / 4, 3.0f, false, 1.0f, 1.0f,
                  {0.7f, 0.7f, 0.7f})));
    } break;
    case test_environments_type::sunsky: {
      add_environment(scene, "sunsky",
          {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}}, {0.5, 0.5, 0.5},
          add_texture(scene, "sky",
              make_sunsky({2048, 1024}, pif / 4, 3.0f, true, 1.0f, 1.0f,
                  {0.7f, 0.7f, 0.7f})));
    } break;
  }
  switch (params.arealights) {
    case test_arealights_type::none: break;
    case test_arealights_type::standard: {
      add_instance(scene, "arealight1",
          lookat_frame(vec3f{-0.4f, 0.8f, 0.8f}, {0.0f, 0.1f, 0.0f},
              {0.0f, 1.0f, 0.0f}, true),
          add_shape(scene, "arealight1", make_rect({1, 1}, {0.2f, 0.2f})),
          add_emission_material(scene, "arealight1", {20, 20, 20}));
      add_instance(scene, "arealight2",
          lookat_frame(vec3f{+0.4f, 0.8f, 0.8f}, {0.0f, 0.1f, 0.0f},
              {0.0f, 1.0f, 0.0f}, true),
          add_shape(scene, "arealight2", make_rect({1, 1}, {0.2, 0.2})),
          add_emission_material(scene, "arealight2", {20, 20, 20}));
    } break;
    case test_arealights_type::large: {
      add_instance(scene, "largearealight1",
          lookat_frame(vec3f{-0.8f, 1.6f, 1.6f}, {0.0f, 0.1f, 0.0f},
              {0.0f, 1.0f, 0.0f}, true),
          add_shape(scene, "largearealight1", make_rect({1, 1}, {0.4f, 0.4f})),
          add_emission_material(scene, "largearealight1", {10, 10, 10}));
      add_instance(scene, "largearealight2",
          lookat_frame(vec3f{+0.8f, 1.6f, 1.6f}, {0.0f, 0.1f, 0.0f},
              {0.0f, 1.0f, 0.0f}, true),
          add_shape(scene, "largearealight2", make_rect({1, 1}, {0.4f, 0.4f})),
          add_emission_material(scene, "largearealight2", {10, 10, 10}));
    } break;
  }
  switch (params.floor) {
    case test_floor_type::none: break;
    case test_floor_type::standard: {
      add_instance(scene, "floor", {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}},
          add_shape(scene, "floor", make_floor({1, 1}, {2, 2}, {20, 20})),
          add_matte_material(scene, "floor", {1, 1, 1},
              add_texture(scene, "floor", make_grid({1024, 1024}))));
    } break;
  }
  auto shapes = vector<int>{}, shapesi = vector<int>{};
  auto subdivs   = vector<int>{};
  auto materials = vector<int>{};
  switch (params.shapes) {
    case test_shapes_type::features1: {
      auto bunny  = add_shape(scene, "sphere", make_sphere(32, 0.075, 1));
      auto sphere = add_shape(scene, "sphere", make_sphere(32, 0.075, 1));
      shapes      = {bunny, sphere, bunny, sphere, bunny};
    } break;
    case test_shapes_type::features2: {
      shapes  = {add_shape(scene, "sphere", make_sphere(32, 0.075, 1)),
           add_shape(scene, "suzanne", make_monkey(0, 0.075f * 0.8f)),
           add_shape(scene, "hair",
               make_random_hairs(make_sphere(32, 0.075f * 0.8f, 1), 65536, 4,
                   {0.1f * 0.15f, 0.1f * 0.15f},
                   {0.001f * 0.15f, 0.0005f * 0.15f}, 0.03)),
           add_shape(scene, "displaced", make_sphere(128, 0.075f, 1)),
           add_shape(scene, "cube",
               make_rounded_box({32, 32, 32}, {0.075, 0.075, 0.075},
                   0.3 * 0.075f, {1, 1, 1}))};
      shapesi = {invalidid, invalidid,
          add_shape(scene, "hairi", make_sphere(32, 0.075f * 0.8f, 1)),
          invalidid, invalidid};
      subdivs = {add_subdiv(scene, "suzanne", make_monkey(0, 0.075f * 0.8f),
                     shapes[1], 2),
          add_subdiv(scene, "displaced", make_sphere(128, 0.075f, 1), shapes[3],
              0, 0.025,
              add_texture(
                  scene, "bumps-displacement", make_bumps({1024, 1024})))};
    } break;
    case test_shapes_type::rows: {
      auto bunny  = add_shape(scene, "bunny", make_sphere(32, 0.075, 1));
      auto sphere = add_shape(scene, "sphere", make_sphere(32, 0.075, 1));
      shapes      = {bunny, bunny, bunny, bunny, bunny, sphere, sphere, sphere,
               sphere, sphere};
    } break;
    case test_shapes_type::bunny_sphere: {
      auto bunny  = add_shape(scene, "bunny", make_sphere(32, 0.075, 1));
      auto sphere = add_shape(scene, "sphere", make_sphere(32, 0.075, 1));
      shapes      = {bunny, sphere, bunny, sphere, bunny};
    } break;
    case test_shapes_type::shapes1: {
      shapes = {
          add_shape(scene, "sphere", make_sphere(32, 0.075, 1)),
          add_shape(scene, "uvsphere-flipcap",
              make_capped_uvsphere({32, 32}, 0.075, 0.3 * 0.075, {1, 1})),
          add_shape(scene, "disk", make_disk(32, 0.075f, 1)),
          add_shape(scene, "uvcylinder",
              make_rounded_uvcylinder(
                  {32, 32, 32}, {0.075, 0.075}, 0.3 * 0.075, {1, 1, 1})),
          add_shape(scene, "cube",
              make_rounded_box({32, 32, 32}, {0.075, 0.075, 0.075},
                  0.3 * 0.075f, {1, 1, 1})),
      };
    } break;
    case test_shapes_type::shapes2: {
      shapes = {
          add_shape(scene, "cube-subdiv", make_wtcube(1, 0.075)),
          add_shape(scene, "suzanne-subdiv", make_monkey(0, 0.075)),
          add_shape(scene, "displaced", make_sphere(128, 0.075f, 1)),
          add_shape(scene, "bunny", make_sphere(32, 0.075, 1)),
          add_shape(scene, "teapot", make_sphere(32, 0.075, 1)),
      };
      subdivs = {
          add_subdiv(scene, "cube-subdiv", make_fvcube(1, 0.075), shapes[0], 4),
          add_subdiv(
              scene, "suzanne-subdiv", make_monkey(0, 0.075), shapes[1], 2),
          add_subdiv(scene, "displaced", make_sphere(128, 0.075f), shapes[2], 0,
              0.025,
              add_texture(
                  scene, "bumps-displacement", make_bumps({1024, 1024})))};
    } break;
    case test_shapes_type::shapes3: {
      shapes = {
          invalidid,
          add_shape(scene, "hair1",
              make_random_hairs(make_sphere(32, 0.075f * 0.8f), 65536, 4,
                  {0.1f * 0.15f, 0.1f * 0.15f},
                  {0.001f * 0.15f, 0.0005f * 0.15f}, 0.03)),
          add_shape(scene, "hair2",
              make_random_hairs(make_sphere(32, 0.075f * 0.8f), 65536, 4,
                  {0.1f * 0.15f, 0.1f * 0.15f},
                  {0.001f * 0.15f, 0.0005f * 0.15f})),
          add_shape(scene, "hair3",
              make_random_hairs(make_sphere(32, 0.075f * 0.8f), 65536, 4,
                  {0.1f * 0.15f, 0.1f * 0.15f},
                  {0.001f * 0.15f, 0.0005f * 0.15f})),
          invalidid,
      };
    } break;
  }
  switch (params.materials) {
    case test_materials_type::features1: {
      materials = {
          add_glossy_material(scene, "coated", {1, 1, 1}, 0.2,
              add_texture(scene, "uvgrid", make_uvgrid({1024, 1024}))),
          add_refractive_material(scene, "glass", {1, 0.5, 0.5}, 0),
          add_scattering_material(
              scene, "jade", {0.5, 0.5, 0.5}, {0.3, 0.6, 0.3}, 0),
          add_glossy_material(scene, "bumped", {0.5, 0.7, 0.5}, 0.2, invalidid,
              invalidid,
              add_texture(scene, "bumps-normal",
                  bump_to_normal(make_bumps({1024, 1024}), 0.05))),
          add_reflective_material(scene, "metal", {0.66, 0.45, 0.34}, 0.2),
      };
    } break;
    case test_materials_type::features2: {
      auto uvgrid  = add_glossy_material(scene, "uvgrid", {1, 1, 1}, 0.2,
           add_texture(scene, "uvgrid", make_uvgrid({1024, 1024})));
      auto plastic = add_glossy_material(
          scene, "plastic", {0.5, 0.7, 0.5}, 0.2);
      auto hair = add_matte_material(scene, "hair", {0.7, 0.7, 0.7});
      materials = {uvgrid, plastic, hair, plastic, uvgrid};
    } break;
    case test_materials_type::uvgrid: {
      auto uvgrid = add_glossy_material(scene, "uvgrid", {1, 1, 1}, 0.2,
          add_texture(scene, "uvgrid", make_uvgrid({1024, 1024})));
      materials   = {uvgrid, uvgrid, uvgrid, uvgrid, uvgrid};
    } break;
    case test_materials_type::hair: {
      auto hair = add_matte_material(scene, "hair", {0.7, 0.7, 0.7});
      materials = {hair, hair, hair, hair, hair};
    } break;
    case test_materials_type::plastic_metal: {
      materials = {
          add_glossy_material(scene, "plastic1", {0.5, 0.5, 0.7}, 0.01),
          add_glossy_material(scene, "plastic2", {0.5, 0.7, 0.5}, 0.2),
          add_matte_material(scene, "matte", {0.7, 0.7, 0.7}),
          add_reflective_material(scene, "metal1", {0.7, 0.7, 0.7}, 0),
          add_reflective_material(scene, "metal2", {0.66, 0.45, 0.34}, 0.2),
      };
    } break;
    case test_materials_type::materials1: {
      materials = {
          add_glossy_material(scene, "plastic1", {0.5, 0.5, 0.7}, 0.01),
          add_glossy_material(scene, "plastic2", {0.5, 0.7, 0.5}, 0.2),
          add_matte_material(scene, "matte", {0.7, 0.7, 0.7}),
          add_glossy_material(scene, "metal1", {0.7, 0.7, 0.7}, 0),
          add_glossy_material(scene, "metal2", {0.66, 0.45, 0.34}, 0.2),
      };
    } break;
    case test_materials_type::materials2: {
      materials = {
          add_glossy_material(scene, "glass1", {1, 1, 1}, 0),
          add_glossy_material(scene, "glass2", {1, 0.7, 0.7}, 0.1),
          add_transparent_material(scene, "transparent", {0.7, 0.5, 0.5}, 0.2),
          add_transparent_material(scene, "tglass1", {1, 1, 1}, 0),
          add_transparent_material(scene, "tglass2", {1, 0.7, 0.7}, 0.1),
      };
    } break;
    case test_materials_type::materials3: {
      auto bumps_normal = add_texture(scene, "bumps-normal",
          bump_to_normal(make_bumps({1024, 1024}), 0.05));
      materials         = {
          add_glossy_material(scene, "plastic1", {0.5, 0.5, 0.7}, 0.01,
                      invalidid, invalidid, bumps_normal),
          add_glossy_material(scene, "plastic2", {0.5, 0.7, 0.5}, 0.2),
          add_reflective_material(scene, "metal1", {0.7, 0.7, 0.7}, 0,
                      invalidid, invalidid, bumps_normal),
          add_reflective_material(scene, "metal2", {0.66, 0.45, 0.34}, 0.2),
          add_reflective_material(scene, "metal3", {0.66, 0.45, 0.34}, 0.2),
      };
    } break;
    case test_materials_type::materials4: {
      materials = {
          add_volumetric_material(
              scene, "cloud", {0.65, 0.65, 0.65}, {0.9, 0.9, 0.9}, 1),
          add_refractive_material(scene, "glass", {1, 0.5, 0.5}, 0),
          add_scattering_material(
              scene, "jade", {0.5, 0.5, 0.5}, {0.3, 0.6, 0.3}, 0),
          add_scattering_material(
              scene, "jade2", {0.5, 0.5, 0.5}, {0.3, 0.6, 0.3}, 0),
          add_volumetric_material(
              scene, "smoke", {0.5, 0.5, 0.5}, {0.2, 0.2, 0.2}),
      };
    } break;
    case test_materials_type::materials5: {
      materials = {
          add_scattering_material(scene, "skin1a", {0.76, 0.48, 0.23},
              {0.436, 0.227, 0.131}, 0.25, invalidid, invalidid, invalidid,
              invalidid, 1.5, -0.8, 0.001),
          add_scattering_material(scene, "skin2a", {0.82, 0.55, 0.4},
              {0.623, 0.433, 0.343}, 0.25, invalidid, invalidid, invalidid,
              invalidid, 1.5, -0.8, 0.001),
          add_scattering_material(scene, "skins", {0.76, 0.48, 0.23},
              {0.436, 0.227, 0.131}, 0, invalidid, invalidid, invalidid,
              invalidid, 1.5, -0.8, 0.001),
          add_scattering_material(scene, "skin1b", {0.76, 0.48, 0.23},
              {0.436, 0.227, 0.131}, 0.25, invalidid, invalidid, invalidid,
              invalidid, 1.5, -0.8, 0.001),
          add_scattering_material(scene, "skin2b", {0.82, 0.55, 0.4},
              {0.623, 0.433, 0.343}, 0.25, invalidid, invalidid, invalidid,
              invalidid, 1.5, -0.8, 0.001),
      };
    } break;
  }
  for (auto idx : range(shapes.size())) {
    if (!shapes[idx]) continue;
    if (shapes.size() > 5) {
      add_instance(scene,
          scene.shape_names[idx] + "-" + scene.shape_names[idx % 5],
          {{1, 0, 0}, {0, 1, 0}, {0, 0, 1},
              {0.2f * (idx % 5 - 2), 0.075, -0.4f * (idx / 5)}},
          shapes[idx], materials[idx % 5]);
    } else {
      auto name = params.instance_name == test_instance_name_type::material
                      ? scene.material_names[idx]
                      : scene.shape_names[idx];
      add_instance(scene, name,
          {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.2f * (idx % 5 - 2), 0.075, 0}},
          shapes[idx], materials[idx]);
    }
    if (!shapesi.empty() && shapesi[idx]) {
      // TODO(fabio): fix name
      add_instance(scene, "",
          {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.2f * (idx - 2), 0.075, 0}},
          shapesi[idx], materials[idx]);
    }
  }
  return scene;
}

// Scene presets used for testing.
scene_data make_scene_preset(const string& type_) {
  auto type = path_basename(type_);
  if (type == "cornellbox") {
    return make_cornellbox();
  } else if (type == "features1") {
    return make_features1_scene();
  } else if (type == "features2") {
    return make_features2_scene();
  } else if (type == "materials1") {
    return make_materials1_scene();
  } else if (type == "materials2") {
    return make_materials2_scene();
  } else if (type == "materials3") {
    return make_materials4_scene();
  } else if (type == "materials4") {
    return make_materials4_scene();
  } else if (type == "shapes1") {
    return make_shapes1_scene();
  } else if (type == "shapes2") {
    return make_shapes2_scene();
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
  } else if (type.starts_with("shape_")) {
    auto scene = scene_data{};
    scene.shapes.push_back(make_shape_preset(type.substr(6)));
    scene.instances.push_back({identity3x4f, 0, invalidid});
    add_missing_material(scene);
    add_missing_camera(scene);
    add_missing_radius(scene);
    add_missing_lights(scene);
    return scene;
  } else if (type.starts_with("tshape_")) {
    auto scene = scene_data{};
    scene.shapes.push_back(make_shape_preset(type.substr(7)));
    scene.instances.push_back({identity3x4f, 0, invalidid});
    add_missing_material(scene, true, make_uvgrid());
    add_missing_camera(scene);
    add_missing_radius(scene);
    add_missing_lights(scene);
    return scene;
  } else if (type.starts_with("uvshape_")) {
    auto scene = scene_data{};
    scene.shapes.push_back(make_shape_preset(type.substr(8)));
    scene.instances.push_back({identity3x4f, 0, invalidid});
    add_missing_material(scene, true, make_uvramp());
    add_missing_camera(scene);
    add_missing_radius(scene);
    add_missing_lights(scene);
    return scene;
  } else if (type.starts_with("orshape_")) {
    auto scene = scene_data{};
    scene.shapes.push_back(make_shape_preset(type.substr(8)));
    scene.instances.push_back({identity3x4f, 0, invalidid});
    add_missing_material(scene, true, make_orgrid());
    add_missing_camera(scene);
    add_missing_radius(scene);
    add_missing_lights(scene);
    return scene;
  } else {
    throw io_error{"unknown preset " + type};
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GENERIC SCENE LOADING
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a scene in the builtin JSON format.
static void load_json_scene(
    const string& filename, scene_data& scene, bool noparallel);
static void save_json_scene(
    const string& filename, const scene_data& scene, bool noparallel);

// Load/save a scene from/to OBJ.
static void load_obj_scene(
    const string& filename, scene_data& scene, bool noparallel);
static void save_obj_scene(
    const string& filename, const scene_data& scene, bool noparallel);

// Load/save a scene from/to PLY. Loads/saves only one mesh with no other
// data.
static void load_ply_scene(
    const string& filename, scene_data& scene, bool noparallel);
static void save_ply_scene(
    const string& filename, const scene_data& scene, bool noparallel);

// Load/save a scene from/to STL. Loads/saves only one mesh with no other
// data.
static void load_stl_scene(
    const string& filename, scene_data& scene, bool noparallel);
static void save_stl_scene(
    const string& filename, const scene_data& scene, bool noparallel);

// Load/save a scene from/to glTF.
static void load_gltf_scene(
    const string& filename, scene_data& scene, bool noparallel);
static void save_gltf_scene(
    const string& filename, const scene_data& scene, bool noparallel);

// Load/save a scene from/to pbrt. This is not robust at all and only
// works on scene that have been previously adapted since the two renderers
// are too different to match.
static void load_pbrt_scene(
    const string& filename, scene_data& scene, bool noparallel);
static void save_pbrt_scene(
    const string& filename, const scene_data& scene, bool noparallel);

// Load/save a scene from/to mitsuba. This is not robust at all and only
// works on scene that have been previously adapted since the two renderers
// are too different to match. For now, only saving is allowed.
static void load_mitsuba_scene(
    const string& filename, scene_data& scene, bool noparallel);
static void save_mitsuba_scene(
    const string& filename, const scene_data& scene, bool noparallel);

// Load a scene
void load_scene(const string& filename, scene_data& scene, bool noparallel) {
  auto ext = path_extension(filename);
  if (ext == ".json" || ext == ".JSON") {
    return load_json_scene(filename, scene, noparallel);
  } else if (ext == ".obj" || ext == ".OBJ") {
    return load_obj_scene(filename, scene, noparallel);
  } else if (ext == ".gltf" || ext == ".GLTF") {
    return load_gltf_scene(filename, scene, noparallel);
  } else if (ext == ".pbrt" || ext == ".PBRT") {
    return load_pbrt_scene(filename, scene, noparallel);
  } else if (ext == ".xml" || ext == ".XML") {
    return load_mitsuba_scene(filename, scene, noparallel);
  } else if (ext == ".ply" || ext == ".PLY") {
    return load_ply_scene(filename, scene, noparallel);
  } else if (ext == ".stl" || ext == ".STL") {
    return load_stl_scene(filename, scene, noparallel);
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    scene = make_scene_preset(filename);
  } else {
    throw io_error("unsupported format " + filename);
  }
}

// Save a scene
void save_scene(
    const string& filename, const scene_data& scene, bool noparallel) {
  auto ext = path_extension(filename);
  if (ext == ".json" || ext == ".JSON") {
    return save_json_scene(filename, scene, noparallel);
  } else if (ext == ".obj" || ext == ".OBJ") {
    return save_obj_scene(filename, scene, noparallel);
  } else if (ext == ".gltf" || ext == ".GLTF") {
    return save_gltf_scene(filename, scene, noparallel);
  } else if (ext == ".pbrt" || ext == ".PBRT") {
    return save_pbrt_scene(filename, scene, noparallel);
  } else if (ext == ".xml" || ext == ".XML") {
    return save_mitsuba_scene(filename, scene, noparallel);
  } else if (ext == ".ply" || ext == ".PLY") {
    return save_ply_scene(filename, scene, noparallel);
  } else if (ext == ".stl" || ext == ".STL") {
    return save_stl_scene(filename, scene, noparallel);
  } else {
    throw io_error("unsupported format " + filename);
  }
}

// Load/save a scene
scene_data load_scene(const string& filename, bool noparallel) {
  auto scene = scene_data{};
  load_scene(filename, scene, noparallel);
  return scene;
}

// Make missing scene directories
void make_scene_directories(const string& filename, const scene_data& scene) {
  // make a directory if needed
  make_directory(path_dirname(filename));
  if (!scene.shapes.empty())
    make_directory(path_join(path_dirname(filename), "shapes"));
  if (!scene.textures.empty())
    make_directory(path_join(path_dirname(filename), "textures"));
  if (!scene.subdivs.empty())
    make_directory(path_join(path_dirname(filename), "subdivs"));
}

// Add environment
void add_environment(
    scene_data& scene, const string& name, const string& filename) {
  auto texture = load_texture(filename);
  scene.textures.push_back(std::move(texture));
  scene.environments.push_back({{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}},
      {1, 1, 1}, (int)scene.textures.size() - 1});
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// INDIVIDUAL ELEMENTS
// -----------------------------------------------------------------------------
namespace yocto {

// load instances
static void load_instance(const string& filename, vector<frame3f>& frames) {
  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = load_ply(filename);
    // TODO: remove when all as arrays
    if (!get_values(ply, "instance",
            {"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz", "ox", "oy",
                "oz"},
            (vector<array<float, 12>>&)frames)) {
      throw io_error{"cannot parse " + filename};
    }
  } else {
    throw io_error("unsupported format " + filename);
  }
}

// save instances
[[maybe_unused]] static void save_instance(
    const string& filename, const vector<frame3f>& frames, bool ascii = false) {
  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = ply_model{};
    // TODO: remove when all as arrays
    add_values(ply, "instance",
        {"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz", "ox", "oy",
            "oz"},
        (const vector<array<float, 12>>&)frames);
    save_ply(filename, ply);
  } else {
    throw io_error("unsupported format " + filename);
  }
}

// load subdiv
void load_subdiv(const string& filename, subdiv_data& subdiv) {
  auto lsubdiv         = load_fvshape(filename, true);
  subdiv.quadspos      = lsubdiv.quadspos;
  subdiv.quadsnorm     = lsubdiv.quadsnorm;
  subdiv.quadstexcoord = lsubdiv.quadstexcoord;
  subdiv.positions     = lsubdiv.positions;
  subdiv.normals       = lsubdiv.normals;
  subdiv.texcoords     = lsubdiv.texcoords;
}

// save subdiv
void save_subdiv(const string& filename, const subdiv_data& subdiv) {
  auto ssubdiv          = fvshape_data{};
  ssubdiv.quadspos      = subdiv.quadspos;
  ssubdiv.quadsnorm     = subdiv.quadsnorm;
  ssubdiv.quadstexcoord = subdiv.quadstexcoord;
  ssubdiv.positions     = subdiv.positions;
  ssubdiv.normals       = subdiv.normals;
  ssubdiv.texcoords     = subdiv.texcoords;
  save_fvshape(filename, ssubdiv, true);
}

// load/save subdiv
subdiv_data load_subdiv(const string& filename) {
  auto subdiv = subdiv_data{};
  load_subdiv(filename, subdiv);
  return subdiv;
}

// save binary shape
static void save_binshape(const string& filename, const shape_data& shape) {
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

  save_binary(filename, buffer);
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
static void load_json_scene_version40(const string& filename,
    const json_value& json, scene_data& scene, bool noparallel) {
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
      ply_instances.push_back({});
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
          auto lookat = mat3f{};
          get_om3(element, "lookat", lookat);
          auto from = lookat[0], to = lookat[1], up = lookat[2];
          camera.focus = length(from - to);
          camera.frame = lookat_frame(from, to, up);
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
          auto lookat = mat3f{};
          get_om3(element, "lookat", lookat);
          auto from = lookat[0], to = lookat[1], up = lookat[2];
          environment.frame = lookat_frame(from, to, up, false);
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
          auto lookat = mat3f{};
          get_om3(element, "lookat", lookat);
          auto from = lookat[0], to = lookat[1], up = lookat[2];
          instance.frame = lookat_frame(from, to, up, false);
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
          auto lookat = mat3f{};
          get_om3(element, "lookat", lookat);
          auto from = lookat[0], to = lookat[1], up = lookat[2];
          instance.frame = lookat_frame(from, to, up, false);
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
    throw io_error("cannot parse " + filename);
  }

  // dirname
  auto dirname = path_dirname(filename);

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
  try {
    if (noparallel) {
      // load shapes
      for (auto& shape : scene.shapes) {
        auto path = find_path(
            get_shape_name(scene, shape), "shapes", {".ply", ".obj"});
        load_shape(path_join(dirname, path), shape, true);
      }
      // load subdivs
      for (auto& subdiv : scene.subdivs) {
        auto path = find_path(
            get_subdiv_name(scene, subdiv), "subdivs", {".ply", ".obj"});
        load_subdiv(path_join(dirname, path), subdiv);
      }
      // load textures
      for (auto& texture : scene.textures) {
        auto path = find_path(get_texture_name(scene, texture), "textures",
            {".hdr", ".exr", ".png", ".jpg"});
        load_texture(path_join(dirname, path), texture);
      }
      // load instances
      for (auto& ply_instance : ply_instances) {
        auto path = find_path(
            get_ply_instance_name(scene, ply_instance), "instances", {".ply"});
        load_instance(path_join(dirname, path), ply_instance.frames);
      }
    } else {
      // load shapes
      parallel_foreach(scene.shapes, [&](auto& shape) {
        auto path = find_path(
            get_shape_name(scene, shape), "shapes", {".ply", ".obj"});
        return load_shape(path_join(dirname, path), shape, true);
      });
      // load subdivs
      parallel_foreach(scene.subdivs, [&](auto& subdiv) {
        auto path = find_path(
            get_subdiv_name(scene, subdiv), "subdivs", {".ply", ".obj"});
        return load_subdiv(path_join(dirname, path), subdiv);
      });
      // load textures
      parallel_foreach(scene.textures, [&](auto& texture) {
        auto path = find_path(get_texture_name(scene, texture), "textures",
            {".hdr", ".exr", ".png", ".jpg"});
        return load_texture(path_join(dirname, path), texture);
      });
      // load instances
      parallel_foreach(ply_instances, [&](auto& ply_instance) {
        auto path = find_path(
            get_ply_instance_name(scene, ply_instance), "instances", {".ply"});
        return load_instance(path_join(dirname, path), ply_instance.frames);
      });
    }
  } catch (std::exception& except) {
    throw io_error(
        "cannot load " + filename + " since " + string(except.what()));
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
}

// Load a scene in the builtin JSON format.
static void load_json_scene_version41(const string& filename, json_value& json,
    scene_data& scene, bool noparallel) {
  // check version
  if (!json.contains("asset") || !json.at("asset").contains("version"))
    return load_json_scene_version40(filename, json, scene, noparallel);

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
          auto lookat = mat3f{};
          get_opt(element, "lookat", lookat);
          auto from = lookat[0], to = lookat[1], up = lookat[2];
          camera.focus = length(from - to);
          camera.frame = lookat_frame(from, to, up);
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
          auto lookat = mat3f{};
          get_opt(element, "lookat", lookat);
          auto from = lookat[0], to = lookat[1], up = lookat[2];
          instance.frame = lookat_frame(from, to, up, false);
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
          auto lookat = mat3f{};
          get_opt(element, "lookat", lookat);
          auto from = lookat[0], to = lookat[1], up = lookat[2];
          environment.frame = lookat_frame(from, to, up, false);
        }
      }
    }
  } catch (...) {
    throw io_error("cannot parse " + filename);
  }

  // prepare data
  auto dirname = path_dirname(filename);

  // fix paths
  for (auto& datafile : shape_filenames)
    datafile = path_join(dirname, "shapes", datafile);
  for (auto& datafile : texture_filenames)
    datafile = path_join(dirname, "textures", datafile);
  for (auto& datafile : subdiv_filenames)
    datafile = path_join(dirname, "subdivs", datafile);

  // load resources
  if (noparallel) {
    try {
      // load shapes
      for (auto&& [filename, shape] : zip(shape_filenames, scene.shapes)) {
        load_shape(filename, shape, true);
      }
      // load subdivs
      for (auto&& [filename, subdiv] : zip(subdiv_filenames, scene.subdivs)) {
        load_subdiv(filename, subdiv);
      }
      // load textures
      for (auto&& [filename, texture] :
          zip(texture_filenames, scene.textures)) {
        load_texture(filename, texture);
      }
    } catch (std::exception& except) {
      throw io_error(
          "cannot load " + filename + " since " + string(except.what()));
    }
  } else {
    try {
      // load shapes
      parallel_zip(
          shape_filenames, scene.shapes, [&](auto&& filename, auto&& shape) {
            return load_shape(filename, shape, true);
          });
      // load subdivs
      parallel_zip(
          subdiv_filenames, scene.subdivs, [&](auto&& filename, auto&& subdiv) {
            return load_subdiv(filename, subdiv);
          });
      // load textures
      parallel_zip(texture_filenames, scene.textures,
          [&](auto&& filename, auto&& texture) {
            return load_texture(filename, texture);
          });
    } catch (std::exception& except) {
      throw io_error(
          "cannot load " + filename + " since " + string(except.what()));
    }
  }

  // fix scene
  add_missing_camera(scene);
  add_missing_radius(scene);
  trim_memory(scene);
}

// Load a scene in the builtin JSON format.
static void load_json_scene(
    const string& filename, scene_data& scene, bool noparallel) {
  // open file
  auto json = load_json(filename);

  // check version
  if (!json.contains("asset") || !json.at("asset").contains("version"))
    load_json_scene_version40(filename, json, scene, noparallel);
  if (json.contains("asset") && json.at("asset").contains("version") &&
      json.at("asset").at("version") == "4.1")
    load_json_scene_version41(filename, json, scene, noparallel);

  // parse json value
  auto get_opt = [](const json_value& json, const string& key, auto& value) {
    value = json.value(key, value);
  };

  // filenames
  auto shape_filenames   = vector<string>{};
  auto texture_filenames = vector<string>{};
  auto subdiv_filenames  = vector<string>{};

  // parsing values
  try {
    if (json.contains("asset")) {
      auto& element = json.at("asset");
      get_opt(element, "copyright", scene.copyright);
      auto version = string{};
      get_opt(element, "version", version);
      if (version != "4.2" && version != "5.0")
        throw io_error("unsupported format version " + filename);
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
        if (element.contains("lookat")) {
          auto lookat = mat3f{};
          get_opt(element, "lookat", lookat);
          auto from = lookat[0], to = lookat[1], up = lookat[2];
          camera.focus = length(from - to);
          camera.frame = lookat_frame(from, to, up);
        }
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
        get_opt(element, "nearest", texture.nearest);
        get_opt(element, "clamp", texture.clamp);
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
        if (element.contains("lookat")) {
          auto lookat = mat3f{};
          get_opt(element, "lookat", lookat);
          auto from = lookat[0], to = lookat[1], up = lookat[2];
          instance.frame = lookat_frame(from, to, up, true);
        }
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
        if (element.contains("lookat")) {
          auto lookat = mat3f{};
          get_opt(element, "lookat", lookat);
          auto from = lookat[0], to = lookat[1], up = lookat[2];
          environment.frame = lookat_frame(from, to, up, true);
        }
      }
    }
  } catch (...) {
    throw io_error("cannot parse " + filename);
  }

  // prepare data
  auto dirname = path_dirname(filename);

  // load resources
  if (noparallel) {
    try {
      // load shapes
      for (auto&& [filename, shape] : zip(shape_filenames, scene.shapes)) {
        load_shape(path_join(dirname, filename), shape, true);
      }
      // load subdivs
      for (auto&& [filename, subdiv] : zip(subdiv_filenames, scene.subdivs)) {
        load_subdiv(path_join(dirname, filename), subdiv);
      }
      // load textures
      for (auto&& [filename, texture] :
          zip(texture_filenames, scene.textures)) {
        load_texture(path_join(dirname, filename), texture);
      }
    } catch (std::exception& except) {
      throw io_error(
          "cannot load " + filename + " since " + string(except.what()));
    }
  } else {
    try {
      // load shapes
      parallel_zip(
          shape_filenames, scene.shapes, [&](auto&& filename, auto&& shape) {
            return load_shape(path_join(dirname, filename), shape, true);
          });
      // load subdivs
      parallel_zip(
          subdiv_filenames, scene.subdivs, [&](auto&& filename, auto&& subdiv) {
            return load_subdiv(path_join(dirname, filename), subdiv);
          });
      // load textures
      parallel_zip(texture_filenames, scene.textures,
          [&](auto&& filename, auto&& texture) {
            return load_texture(path_join(dirname, filename), texture);
          });
    } catch (std::exception& except) {
      throw io_error(
          "cannot load " + filename + " since " + string(except.what()));
    }
  }

  // fix scene
  add_missing_camera(scene);
  add_missing_radius(scene);
  trim_memory(scene);
}

// Save a scene in the builtin JSON format.
static void save_json_scene(
    const string& filename, const scene_data& scene, bool noparallel) {
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
    auto  default_ = texture_data{};
    auto& group    = add_array(json, "textures");
    reserve_values(group, scene.textures.size());
    for (auto&& [idx, texture] : enumerate(scene.textures)) {
      auto& element = append_object(group);
      set_val(element, "name", get_name(scene.texture_names, idx), "");
      set_val(element, "uri", texture_filenames[idx], ""s);
      set_val(element, "nearest", texture.nearest, default_.nearest);
      set_val(element, "clamp", texture.clamp, default_.clamp);
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
          element, "catmullclark", subdiv.catmullclark, default_.catmullclark);
      set_val(element, "smooth", subdiv.smooth, default_.smooth);
      set_val(
          element, "displacement", subdiv.displacement, default_.displacement);
      set_val(element, "displacement_tex", subdiv.displacement_tex,
          default_.displacement_tex);
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
  save_json(filename, json);

  // prepare data
  auto dirname = path_dirname(filename);

  // dirname
  if (noparallel) {
    try {
      // save shapes
      for (auto&& [filename, shape] : zip(shape_filenames, scene.shapes)) {
        save_shape(path_join(dirname, filename), shape, true);
      }
      // save subdiv
      for (auto&& [filename, subdiv] : zip(subdiv_filenames, scene.subdivs)) {
        save_subdiv(path_join(dirname, filename), subdiv);
      }
      // save textures
      for (auto&& [filename, texture] :
          zip(texture_filenames, scene.textures)) {
        save_texture(path_join(dirname, filename), texture);
      }
    } catch (std::exception& except) {
      throw io_error(
          "cannot save " + filename + " since " + string(except.what()));
    }
  } else {
    try {
      // save shapes
      parallel_zip(
          shape_filenames, scene.shapes, [&](auto&& filename, auto&& shape) {
            return save_shape(path_join(dirname, filename), shape, true);
          });
      // save subdivs
      parallel_zip(
          subdiv_filenames, scene.subdivs, [&](auto&& filename, auto&& subdiv) {
            return save_subdiv(path_join(dirname, filename), subdiv);
          });
      // save textures
      parallel_zip(texture_filenames, scene.textures,
          [&](auto&& filename, auto&& texture) {
            return save_texture(path_join(dirname, filename), texture);
          });
    } catch (std::exception& except) {
      throw io_error(
          "cannot save " + filename + " since " + string(except.what()));
    }
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// OBJ CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

// Loads an OBJ
static void load_obj_scene(
    const string& filename, scene_data& scene, bool noparallel) {
  // load obj
  auto obj = load_obj(filename, false, true);

  // convert cameras
  scene.cameras.reserve(obj.cameras.size());
  for (auto& ocamera : obj.cameras) {
    auto& camera        = scene.cameras.emplace_back();
    camera.frame        = to_math(ocamera.frame);
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
    material.emission     = to_math(omaterial.emission);
    material.emission_tex = omaterial.emission_tex;
    if (max(to_math(omaterial.transmission)) > 0.1) {
      material.type      = material_type::transparent;
      material.color     = to_math(omaterial.transmission);
      material.color_tex = omaterial.transmission_tex;
    } else if (max(to_math(omaterial.specular)) > 0.2) {
      material.type      = material_type::reflective;
      material.color     = to_math(omaterial.specular);
      material.color_tex = omaterial.specular_tex;
    } else if (max(to_math(omaterial.specular)) > 0) {
      material.type      = material_type::glossy;
      material.color     = to_math(omaterial.diffuse);
      material.color_tex = omaterial.diffuse_tex;
    } else {
      material.type      = material_type::matte;
      material.color     = to_math(omaterial.diffuse);
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
    get_positions(oshape, (vector<array<float, 3>>&)shape.positions);
    get_normals(oshape, (vector<array<float, 3>>&)shape.normals);
    get_texcoords(oshape, (vector<array<float, 2>>&)shape.texcoords, true);
    get_faces(oshape, instance.material,
        (vector<array<int, 3>>&)shape.triangles,
        (vector<array<int, 4>>&)shape.quads);
    get_lines(oshape, instance.material, (vector<array<int, 2>>&)shape.lines);
    get_points(oshape, instance.material, shape.points);
  }

  // convert environments
  scene.environments.reserve(obj.environments.size());
  for (auto& oenvironment : obj.environments) {
    auto& environment        = scene.environments.emplace_back();
    environment.frame        = to_math(oenvironment.frame);
    environment.emission     = to_math(oenvironment.emission);
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
  auto dirname = path_dirname(filename);

  try {
    if (noparallel) {
      // load textures
      for (auto&& [path, texture] : zip(texture_paths, scene.textures)) {
        load_texture(path_join(dirname, path), texture);
      }
    } else {
      // load textures
      parallel_zip(
          texture_paths, scene.textures, [&](auto&& path, auto&& texture) {
            return load_texture(path_join(dirname, path), texture);
          });
    }
  } catch (std::exception& except) {
    throw io_error(
        "cannot load " + filename + " since " + string(except.what()));
  }

  // fix scene
  add_missing_camera(scene);
  add_missing_radius(scene);
}

static void save_obj_scene(
    const string& filename, const scene_data& scene, bool noparallel) {
  // build obj
  auto obj = obj_model{};

  // convert cameras
  for (auto& camera : scene.cameras) {
    auto& ocamera    = obj.cameras.emplace_back();
    ocamera.name     = get_camera_name(scene, camera);
    ocamera.frame    = to_array(camera.frame);
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
    omaterial.emission     = to_array(material.emission);
    omaterial.diffuse      = to_array(material.color);
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
    add_positions(oshape, (const vector<array<float, 3>>&)positions);
    add_normals(oshape, (const vector<array<float, 3>>&)normals);
    add_texcoords(
        oshape, (const vector<array<float, 2>>&)shape.texcoords, true);
    add_triangles(oshape, (const vector<array<int, 3>>&)shape.triangles,
        instance.material, !shape.normals.empty(), !shape.texcoords.empty());
    add_quads(oshape, (const vector<array<int, 4>>&)shape.quads,
        instance.material, !shape.normals.empty(), !shape.texcoords.empty());
    add_lines(oshape, (const vector<array<int, 2>>&)shape.lines,
        instance.material, !shape.normals.empty(), !shape.texcoords.empty());
    add_points(oshape, shape.points, instance.material, !shape.normals.empty(),
        !shape.texcoords.empty());
  }

  // convert environments
  for (auto& environment : scene.environments) {
    auto& oenvironment        = obj.environments.emplace_back();
    oenvironment.name         = get_environment_name(scene, environment);
    oenvironment.frame        = to_array(environment.frame);
    oenvironment.emission     = to_array(environment.emission);
    oenvironment.emission_tex = environment.emission_tex;
  }

  // save obj
  save_obj(filename, obj);

  // dirname
  auto dirname = path_dirname(filename);

  try {
    if (noparallel) {
      // save textures
      for (auto& texture : scene.textures) {
        auto path = "textures/" + get_texture_name(scene, texture) +
                    (!texture.pixelsf.empty() ? ".hdr"s : ".png"s);
        save_texture(path_join(dirname, path), texture);
      }
    } else {
      // save textures
      parallel_foreach(scene.textures, [&](auto& texture) {
        auto path = "textures/" + get_texture_name(scene, texture) +
                    (!texture.pixelsf.empty() ? ".hdr"s : ".png"s);
        return save_texture(path_join(dirname, path), texture);
      });
    }
  } catch (std::exception& except) {
    throw io_error(
        "cannot save " + filename + " since " + string(except.what()));
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PLY CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static void load_ply_scene(
    const string& filename, scene_data& scene, bool noparallel) {
  // load ply mesh and make instance
  auto shape = load_shape(filename, true);
  scene.shapes.push_back(shape);
  scene.instances.push_back({{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}},
      (int)scene.shapes.size() - 1, -1});

  // fix scene
  add_missing_material(scene);
  add_missing_camera(scene);
  add_missing_radius(scene);
  add_missing_lights(scene);
}

static void save_ply_scene(
    const string& filename, const scene_data& scene, bool noparallel) {
  // save shape
  if (scene.shapes.empty()) throw std::invalid_argument{"empty shape"};
  return save_shape(filename, scene.shapes.front(), false);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// STL CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static void load_stl_scene(
    const string& filename, scene_data& scene, bool noparallel) {
  // load ply mesh and make instance
  auto shape = load_shape(filename, true);
  scene.instances.push_back({{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}},
      (int)scene.shapes.size() - 1, -1});

  // fix scene
  add_missing_material(scene);
  add_missing_camera(scene);
  add_missing_radius(scene);
  add_missing_lights(scene);
}

static void save_stl_scene(
    const string& filename, const scene_data& scene, bool noparallel) {
  // save shape
  if (scene.shapes.empty()) throw std::invalid_argument{"empty shape"};
  return save_shape(filename, scene.shapes.front(), false);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GLTF CONVESION
// -----------------------------------------------------------------------------
namespace yocto {

// Load a scene
static void load_gltf_scene(
    const string& filename, scene_data& scene, bool noparallel) {
  // load gltf data
  auto data = load_binary(filename);

  // parse glTF
  auto options = cgltf_options{};
  memset(&options, 0, sizeof(options));
  auto cgltf_ptr    = (cgltf_data*)nullptr;
  auto cgltf_result = cgltf_parse(
      &options, data.data(), data.size(), &cgltf_ptr);
  if (cgltf_result != cgltf_result_success) {
    throw io_error{"cannot parse " + filename};
  }

  // load buffers
  auto dirname_      = path_dirname(filename);
  dirname_           = dirname_.empty() ? string{"./"} : (dirname_ + "/");
  auto buffer_result = cgltf_load_buffers(
      &options, cgltf_ptr, dirname_.c_str());
  if (buffer_result != cgltf_result_success) {
    cgltf_free(cgltf_ptr);
    throw io_error{"cannot load " + filename + " since cannot load buffers"};
  }

  // setup parsing
  auto& cgltf       = *cgltf_ptr;
  auto  cgltf_guard = std::unique_ptr<cgltf_data, void (*)(cgltf_data*)>(
      cgltf_ptr, cgltf_free);

  // convert cameras
  auto cameras = vector<camera_data>{};
  for (auto idx : range(cgltf.cameras_count)) {
    auto& gcamera = cgltf.cameras[idx];
    auto& camera  = cameras.emplace_back();
    if (gcamera.type == cgltf_camera_type_orthographic) {
      auto& gortho  = gcamera.data.orthographic;
      auto  xmag    = gortho.xmag;
      auto  ymag    = gortho.ymag;
      camera.aspect = xmag / ymag;
      camera.lens   = ymag;  // this is probably bogus
      camera.film   = 0.036f;
    } else if (gcamera.type == cgltf_camera_type_perspective) {
      auto& gpersp  = gcamera.data.perspective;
      camera.aspect = (bool)gpersp.has_aspect_ratio ? gpersp.aspect_ratio
                                                    : 0.0f;
      auto yfov     = gpersp.yfov;
      if (camera.aspect == 0) camera.aspect = 16.0f / 9.0f;
      camera.film = 0.036f;
      if (camera.aspect >= 1) {
        camera.lens = (camera.film / camera.aspect) / (2 * tan(yfov / 2));
      } else {
        camera.lens = camera.film / (2 * tan(yfov / 2));
      }
      camera.focus = 1;
    } else {
      throw io_error{
          "cannot load " + filename + " for unsupported camera type"};
    }
  }

  // convert color textures
  auto get_texture = [&cgltf](const cgltf_texture_view& gview) -> int {
    if (gview.texture == nullptr) return -1;
    auto& gtexture = *gview.texture;
    if (gtexture.image == nullptr) return -1;
    return (int)(gtexture.image - cgltf.images);
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
  for (auto idx : range(cgltf.images_count)) {
    auto& gimage = cgltf.images[idx];
    scene.textures.emplace_back();
    texture_paths.push_back(replace(gimage.uri, "%20", " "));
  }

  // convert materials
  for (auto idx : range(cgltf.materials_count)) {
    auto& gmaterial   = cgltf.materials[idx];
    auto& material    = scene.materials.emplace_back();
    material.type     = material_type::gltfpbr;
    material.emission = {gmaterial.emissive_factor[0],
        gmaterial.emissive_factor[1], gmaterial.emissive_factor[2]};
    if ((bool)gmaterial.has_emissive_strength)
      material.emission *= gmaterial.emissive_strength.emissive_strength;
    material.emission_tex = get_texture(gmaterial.emissive_texture);
    material.normal_tex   = get_texture(gmaterial.normal_texture);
    if ((bool)gmaterial.has_pbr_metallic_roughness) {
      auto& gpbr        = gmaterial.pbr_metallic_roughness;
      material.type     = material_type::gltfpbr;
      material.color    = {gpbr.base_color_factor[0], gpbr.base_color_factor[1],
             gpbr.base_color_factor[2]};
      material.opacity  = gpbr.base_color_factor[3];
      material.metallic = gpbr.metallic_factor;
      material.roughness     = gpbr.roughness_factor;
      material.color_tex     = get_texture(gpbr.base_color_texture);
      material.roughness_tex = get_texture(gpbr.metallic_roughness_texture);
    }
    if ((bool)gmaterial.has_transmission) {
      auto& gtransmission = gmaterial.transmission;
      auto  transmission  = gtransmission.transmission_factor;
      if (transmission > 0) {
        material.type      = material_type::transparent;
        material.color     = {transmission, transmission, transmission};
        material.color_tex = get_texture(gtransmission.transmission_texture);
        // material.roughness = 0; // leave it set from before
      }
    }
  }

  // convert meshes
  auto mesh_primitives = vector<vector<instance_data>>{};
  for (auto idx : range(cgltf.meshes_count)) {
    auto& gmesh      = cgltf.meshes[idx];
    auto& primitives = mesh_primitives.emplace_back();
    if (gmesh.primitives == nullptr) continue;
    for (auto idx : range(gmesh.primitives_count)) {
      auto& gprimitive = gmesh.primitives[idx];
      if (gprimitive.attributes == nullptr) continue;
      auto& shape       = scene.shapes.emplace_back();
      auto& instance    = primitives.emplace_back();
      instance.shape    = (int)scene.shapes.size() - 1;
      instance.material = gprimitive.material != nullptr
                              ? (int)(gprimitive.material - cgltf.materials)
                              : -1;
      for (auto idx : range(gprimitive.attributes_count)) {
        auto& gattribute = gprimitive.attributes[idx];
        auto& gaccessor  = *gattribute.data;
        if ((bool)gaccessor.is_sparse)
          throw io_error{
              "cannot load " + filename + " for unsupported sparse accessor"};
        auto gname       = string{gattribute.name};
        auto count       = gaccessor.count;
        auto components  = cgltf_num_components(gaccessor.type);
        auto dcomponents = components;
        auto data        = (float*)nullptr;
        if (gname == "POSITION") {
          if (components != 3)
            throw io_error{"cannot load " + filename +
                           " for unsupported position components"};

          shape.positions = vector<vec3f>(count);
          data            = (float*)shape.positions.data();
        } else if (gname == "NORMAL") {
          if (components != 3)
            throw io_error{"cannot load " + filename +
                           " for unsupported normal components"};
          shape.normals = vector<vec3f>(count);
          data          = (float*)shape.normals.data();
        } else if (gname == "TEXCOORD" || gname == "TEXCOORD_0") {
          if (components != 2)
            throw io_error{"cannot load " + filename +
                           " for unsupported texture components"};
          shape.texcoords = vector<vec2f>(count);
          data            = (float*)shape.texcoords.data();
        } else if (gname == "COLOR" || gname == "COLOR_0") {
          if (components != 3 && components != 4)
            throw io_error{"cannot load " + filename +
                           " for unsupported color components"};
          shape.colors = vector<vec4f>(count);
          data         = (float*)shape.colors.data();
          if (components == 3) {
            dcomponents = 4;
            for (auto& c : shape.colors) c.w = 1;
          }
        } else if (gname == "TANGENT") {
          if (components != 4)
            throw io_error{"cannot load " + filename +
                           " for unsupported tangent components"};
          shape.tangents = vector<vec4f>(count);
          data           = (float*)shape.tangents.data();
        } else if (gname == "RADIUS") {
          if (components != 1)
            throw io_error{"cannot load " + filename +
                           " for unsupported radius components"};
          shape.radius = vector<float>(count);
          data         = (float*)shape.radius.data();
        } else {
          // ignore
          continue;
        }
        // convert values
        for (auto idx : range(count)) {
          if (!(bool)cgltf_accessor_read_float(
                  &gaccessor, idx, &data[idx * dcomponents], components))
            throw io_error{"cannot load " + filename +
                           " for unsupported accessor conversion"};
        }
        // fixes
        if (gname == "TANGENT") {
          for (auto& t : shape.tangents) t.w = -t.w;
        }
      }
      // indices
      if (gprimitive.indices == nullptr) {
        if (gprimitive.type == cgltf_primitive_type_triangles) {
          shape.triangles = vector<vec3i>(shape.positions.size() / 3);
          for (auto i = 0; i < (int)shape.positions.size() / 3; i++)
            shape.triangles[i] = {i * 3 + 0, i * 3 + 1, i * 3 + 2};
        } else if (gprimitive.type == cgltf_primitive_type_triangle_fan) {
          shape.triangles = vector<vec3i>(shape.positions.size() - 2);
          for (auto i = 2; i < (int)shape.positions.size(); i++)
            shape.triangles[i - 2] = {0, i - 1, i};
        } else if (gprimitive.type == cgltf_primitive_type_triangle_strip) {
          shape.triangles = vector<vec3i>(shape.positions.size() - 2);
          for (auto i = 2; i < (int)shape.positions.size(); i++)
            shape.triangles[i - 2] = {i - 2, i - 1, i};
        } else if (gprimitive.type == cgltf_primitive_type_lines) {
          shape.lines = vector<vec2i>(shape.positions.size() / 2);
          for (auto i = 0; i < (int)shape.positions.size() / 2; i++)
            shape.lines[i] = {i * 2 + 0, i * 2 + 1};
        } else if (gprimitive.type == cgltf_primitive_type_line_loop) {
          shape.lines = vector<vec2i>(shape.positions.size());
          for (auto i = 1; i < (int)shape.positions.size(); i++)
            shape.lines[i - 1] = {i - 1, i};
          shape.lines.back() = {(int)shape.positions.size() - 1, 0};
        } else if (gprimitive.type == cgltf_primitive_type_line_strip) {
          shape.lines = vector<vec2i>(shape.positions.size() - 1);
          for (auto i = 1; i < (int)shape.positions.size(); i++)
            shape.lines[i - 1] = {i - 1, i};
        } else if (gprimitive.type == cgltf_primitive_type_points) {
          throw io_error{
              "cannot load " + filename + " for unsupported point primitive"};

        } else {
          throw io_error{
              "cannot load " + filename + " for unsupported primitive type"};
        }
      } else {
        auto& gaccessor = *gprimitive.indices;
        if (gaccessor.type != cgltf_type_scalar)
          throw io_error{"cannot load " + filename +
                         " for unsupported non-scalar indices"};
        auto indices = vector<int>(gaccessor.count);
        for (auto idx : range(gaccessor.count)) {
          if (!(bool)cgltf_accessor_read_uint(
                  &gaccessor, idx, (cgltf_uint*)&indices[idx], 1))
            throw io_error{
                "cannot load " + filename + " for unsupported accessor type"};
        }
        if (gprimitive.type == cgltf_primitive_type_triangles) {
          shape.triangles = vector<vec3i>(indices.size() / 3);
          for (auto i = 0; i < (int)indices.size() / 3; i++) {
            shape.triangles[i] = {
                indices[i * 3 + 0], indices[i * 3 + 1], indices[i * 3 + 2]};
          }
        } else if (gprimitive.type == cgltf_primitive_type_triangle_fan) {
          shape.triangles = vector<vec3i>(indices.size() - 2);
          for (auto i = 2; i < (int)indices.size(); i++) {
            shape.triangles[i - 2] = {
                indices[0], indices[i - 1], indices[i + 0]};
          }
        } else if (gprimitive.type == cgltf_primitive_type_triangle_strip) {
          shape.triangles = vector<vec3i>(indices.size() - 2);
          for (auto i = 2; i < (int)indices.size(); i++) {
            shape.triangles[i - 2] = {
                indices[i - 2], indices[i - 1], indices[i + 0]};
          }
        } else if (gprimitive.type == cgltf_primitive_type_lines) {
          shape.lines = vector<vec2i>(indices.size() / 2);
          for (auto i = 0; i < (int)indices.size() / 2; i++) {
            shape.lines[i] = {indices[i * 2 + 0], indices[i * 2 + 1]};
          }
        } else if (gprimitive.type == cgltf_primitive_type_line_loop) {
          shape.lines = vector<vec2i>(indices.size());
          for (auto i : range(indices.size())) {
            shape.lines[i] = {
                indices[i + 0], indices[i + 1] % (int)indices.size()};
          }
        } else if (gprimitive.type == cgltf_primitive_type_line_strip) {
          shape.lines = vector<vec2i>(indices.size() - 1);
          for (auto i = 0; i < (int)indices.size() - 1; i++) {
            shape.lines[i] = {indices[i + 0], indices[i + 1]};
          }
        } else if (gprimitive.type == cgltf_primitive_type_points) {
          throw io_error{
              "cannot load " + filename + " for unsupported points indices"};

        } else {
          throw io_error{
              "cannot load " + filename + " for unsupported primitive type"};
        }
      }
    }
  }

  // convert nodes
  for (auto idx : range(cgltf.nodes_count)) {
    auto& gnode = cgltf.nodes[idx];
    if (gnode.camera != nullptr) {
      auto& camera = scene.cameras.emplace_back();
      camera       = cameras.at(gnode.camera - cgltf.cameras);
      auto xform   = mat4f{
            {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
      cgltf_node_transform_world(&gnode, yocto::data(xform));
      camera.frame = to_frame(xform);
    }
    if (gnode.mesh != nullptr) {
      for (auto& primitive : mesh_primitives.at(gnode.mesh - cgltf.meshes)) {
        auto& instance = scene.instances.emplace_back();
        instance       = primitive;
        auto xform     = mat4f{
                {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
        cgltf_node_transform_world(&gnode, yocto::data(xform));
        instance.frame = to_frame(xform);
      }
    }
  }

  // dirname
  auto dirname = path_dirname(filename);

  try {
    if (noparallel) {
      // load texture
      for (auto& texture : scene.textures) {
        auto& path = texture_paths[&texture - &scene.textures.front()];
        load_texture(path_join(dirname, path), texture);
      }
    } else {
      // load textures
      parallel_foreach(scene.textures, [&](auto& texture) {
        auto& path = texture_paths[&texture - &scene.textures.front()];
        return load_texture(path_join(dirname, path), texture);
      });
    }
  } catch (std::exception& except) {
    throw io_error(
        "cannot load " + filename + " since " + string(except.what()));
  }

  // fix scene
  add_missing_material(scene);
  add_missing_camera(scene);
  add_missing_radius(scene);
  add_missing_lights(scene);
}

// Load a scene
static void save_gltf_scene(
    const string& filename, const scene_data& scene, bool noparallel) {
  // prepare data
  auto cgltf_ptr = (cgltf_data*)malloc(sizeof(cgltf_data));
  memset(cgltf_ptr, 0, sizeof(cgltf_data));
  auto& cgltf = *cgltf_ptr;

  // helpers
  auto copy_string = [](const string& str) -> char* {
    if (str.empty()) return nullptr;
    auto copy = (char*)malloc(str.size() + 1);
    if (copy == nullptr) return nullptr;
    memcpy(copy, str.c_str(), str.size() + 1);
    return copy;
  };
  auto alloc_array = [](cgltf_size& size, auto*& ptr, size_t count) {
    using type =
        typename std::remove_pointer_t<std::remove_reference_t<decltype(ptr)>>;
    size = count;
    ptr  = (type*)malloc(count * sizeof(type));
    memset(ptr, 0, count * sizeof(type));
  };

  // asset
  cgltf.asset.version   = copy_string("2.0");
  cgltf.asset.generator = copy_string(
      "Yocto/GL - https://github.com/xelatihy/yocto-gl");
  cgltf.asset.copyright = copy_string(scene.copyright);

  // cameras
  if (!scene.cameras.empty()) {
    alloc_array(cgltf.cameras_count, cgltf.cameras, scene.cameras.size());
    for (auto idx : range(scene.cameras.size())) {
      auto& camera       = scene.cameras[idx];
      auto& gcamera      = cgltf.cameras[idx];
      gcamera.name       = copy_string(get_camera_name(scene, camera));
      gcamera.type       = cgltf_camera_type_perspective;
      auto& gperspective = gcamera.data.perspective;
      gperspective.has_aspect_ratio = 1;
      gperspective.aspect_ratio     = camera.aspect;
      gperspective.yfov             = 0.660593;  // TODO(fabio): yfov
      gperspective.znear            = 0.001;     // TODO(fabio): configurable?
    }
  }

  // textures
  if (!scene.textures.empty()) {
    alloc_array(cgltf.textures_count, cgltf.textures, scene.textures.size());
    alloc_array(cgltf.images_count, cgltf.images, scene.textures.size());
    alloc_array(cgltf.samplers_count, cgltf.samplers, 1);
    auto& gsampler  = cgltf.samplers[0];
    gsampler.name   = copy_string("sampler");
    gsampler.wrap_s = 10497;
    gsampler.wrap_t = 10497;
    for (auto idx : range(scene.textures.size())) {
      auto& texture    = scene.textures[idx];
      auto& gtexture   = cgltf.textures[idx];
      auto& gimage     = cgltf.images[idx];
      auto  name       = get_texture_name(scene, texture);
      gimage.name      = copy_string(name);
      gimage.uri       = copy_string("textures/" + name + ".png");
      gtexture.name    = copy_string(name);
      gtexture.sampler = &gsampler;
      gtexture.image   = &gimage;
    }
  }

  // materials
  if (!scene.materials.empty()) {
    alloc_array(cgltf.materials_count, cgltf.materials, scene.materials.size());
    for (auto idx : range(scene.materials.size())) {
      auto& material  = scene.materials[idx];
      auto& gmaterial = cgltf.materials[idx];
      gmaterial.name  = copy_string(get_material_name(scene, material));
      auto emission_scale =
          max(material.emission) > 1.0f ? max(material.emission) : 1.0f;
      gmaterial.emissive_factor[0] = material.emission.x / emission_scale;
      gmaterial.emissive_factor[1] = material.emission.y / emission_scale;
      gmaterial.emissive_factor[2] = material.emission.z / emission_scale;
      if (emission_scale > 1.0f) {
        gmaterial.has_emissive_strength               = 1;
        gmaterial.emissive_strength.emissive_strength = emission_scale;
      }
      gmaterial.has_pbr_metallic_roughness = 1;
      auto& gpbr                           = gmaterial.pbr_metallic_roughness;
      gpbr.base_color_factor[0]            = material.color.x;
      gpbr.base_color_factor[1]            = material.color.y;
      gpbr.base_color_factor[2]            = material.color.z;
      gpbr.base_color_factor[3]            = material.opacity;
      gpbr.metallic_factor                 = material.metallic;
      gpbr.roughness_factor                = material.roughness;
      if (material.emission_tex != invalidid) {
        gmaterial.emissive_texture.texture = cgltf.textures +
                                             material.emission_tex;
        gmaterial.emissive_texture.scale = 1.0f;
      }
      if (material.normal_tex != invalidid) {
        gmaterial.normal_texture.texture = cgltf.textures + material.normal_tex;
        gmaterial.normal_texture.scale   = 1.0f;
      }
      if (material.color_tex != invalidid) {
        gpbr.base_color_texture.texture = cgltf.textures + material.color_tex;
        gpbr.base_color_texture.scale   = 1.0f;
      }
      if (material.roughness_tex != invalidid) {
        gpbr.metallic_roughness_texture.texture = cgltf.textures +
                                                  material.roughness_tex;
        gpbr.metallic_roughness_texture.scale = 1.0f;
      }
    }
  }

  // buffers
  auto shape_accessor_start = vector<int>{};
  if (!scene.shapes.empty()) {
    alloc_array(cgltf.buffers_count, cgltf.buffers, scene.shapes.size());
    alloc_array(
        cgltf.accessors_count, cgltf.accessors, scene.shapes.size() * 6);
    alloc_array(
        cgltf.buffer_views_count, cgltf.buffer_views, scene.shapes.size() * 6);
    shape_accessor_start     = vector<int>(scene.shapes.size(), 0);
    cgltf.accessors_count    = 0;
    cgltf.buffer_views_count = 0;
    auto add_vertex = [](cgltf_data& cgltf, cgltf_buffer& gbuffer, size_t count,
                          cgltf_type type, const float* data) {
      if (count == 0) return;
      auto  components         = cgltf_num_components(type);
      auto& gbufferview        = cgltf.buffer_views[cgltf.buffer_views_count++];
      gbufferview.buffer       = &gbuffer;
      gbufferview.offset       = gbuffer.size;
      gbufferview.size         = sizeof(float) * components * count;
      gbufferview.type         = cgltf_buffer_view_type_vertices;
      auto& gaccessor          = cgltf.accessors[cgltf.accessors_count++];
      gaccessor.buffer_view    = &gbufferview;
      gaccessor.count          = count;
      gaccessor.type           = type;
      gaccessor.component_type = cgltf_component_type_r_32f;
      gaccessor.has_min        = 1;
      gaccessor.has_max        = 1;
      for (auto component : range(components)) {
        gaccessor.min[component] = flt_max;
        gaccessor.max[component] = flt_min;
        for (auto idx : range(count)) {
          gaccessor.min[component] = min(
              gaccessor.min[component], data[idx * components + component]);
          gaccessor.max[component] = max(
              gaccessor.max[component], data[idx * components + component]);
        }
      }
      gbuffer.size += gbufferview.size;
    };
    auto add_element = [](cgltf_data& cgltf, cgltf_buffer& gbuffer,
                           size_t count, cgltf_type type) {
      if (count == 0) return;
      auto  components         = cgltf_num_components(type);
      auto& gbufferview        = cgltf.buffer_views[cgltf.buffer_views_count++];
      gbufferview.buffer       = &gbuffer;
      gbufferview.offset       = gbuffer.size;
      gbufferview.size         = sizeof(int) * components * count;
      gbufferview.type         = cgltf_buffer_view_type_indices;
      auto& gaccessor          = cgltf.accessors[cgltf.accessors_count++];
      gaccessor.buffer_view    = &gbufferview;
      gaccessor.count          = count * components;
      gaccessor.type           = cgltf_type_scalar;
      gaccessor.component_type = cgltf_component_type_r_32u;
      gbuffer.size += gbufferview.size;
    };
    for (auto idx : range(scene.shapes.size())) {
      auto& shape               = scene.shapes[idx];
      auto& gbuffer             = cgltf.buffers[idx];
      shape_accessor_start[idx] = (int)cgltf.accessors_count;
      gbuffer.uri               = copy_string(
          "shapes/" + get_shape_name(scene, shape) + ".bin");
      add_vertex(cgltf, gbuffer, shape.positions.size(), cgltf_type_vec3,
          (const float*)shape.positions.data());
      add_vertex(cgltf, gbuffer, shape.normals.size(), cgltf_type_vec3,
          (const float*)shape.normals.data());
      add_vertex(cgltf, gbuffer, shape.texcoords.size(), cgltf_type_vec2,
          (const float*)shape.texcoords.data());
      add_vertex(cgltf, gbuffer, shape.colors.size(), cgltf_type_vec4,
          (const float*)shape.colors.data());
      add_vertex(cgltf, gbuffer, shape.radius.size(), cgltf_type_scalar,
          (const float*)shape.radius.data());
      add_element(cgltf, gbuffer, shape.points.size(), cgltf_type_scalar);
      add_element(cgltf, gbuffer, shape.lines.size(), cgltf_type_vec2);
      add_element(cgltf, gbuffer, shape.triangles.size(), cgltf_type_vec3);
      add_element(cgltf, gbuffer, quads_to_triangles(shape.quads).size(),
          cgltf_type_vec3);
    }
  }

  // meshes
  using mesh_key = pair<int, int>;
  struct mesh_hash {
    auto operator()(const pair<int, int>& key) const {
      auto packed = ((uint64_t)key.first << 32) | (uint64_t)key.second;
      return std::hash<uint64_t>()(packed);
    }
  };
  auto mesh_map = unordered_map<mesh_key, size_t, mesh_hash>{};
  if (!scene.instances.empty()) {
    alloc_array(cgltf.meshes_count, cgltf.meshes, scene.instances.size());
    cgltf.meshes_count = 0;
    for (auto idx : range(scene.instances.size())) {
      auto& instance = scene.instances[idx];
      if (mesh_map.find({instance.shape, instance.material}) != mesh_map.end())
        continue;
      mesh_map[{instance.shape, instance.material}] = (int)cgltf.meshes_count;
      auto& shape = scene.shapes[instance.shape];
      auto& gmesh = cgltf.meshes[cgltf.meshes_count++];
      alloc_array(gmesh.primitives_count, gmesh.primitives, 1);
      auto& gprimitive    = gmesh.primitives[0];
      gprimitive.material = cgltf.materials + instance.material;
      alloc_array(gprimitive.attributes_count, gprimitive.attributes, 5);
      gprimitive.attributes_count = 0;
      auto cur_accessor           = shape_accessor_start[instance.shape];
      if (!shape.positions.empty()) {
        auto& gattribute = gprimitive.attributes[gprimitive.attributes_count++];
        gattribute.type  = cgltf_attribute_type_position;
        gattribute.name  = copy_string("POSITION");
        gattribute.data  = cgltf.accessors + cur_accessor++;
      }
      if (!shape.normals.empty()) {
        auto& gattribute = gprimitive.attributes[gprimitive.attributes_count++];
        gattribute.type  = cgltf_attribute_type_normal;
        gattribute.name  = copy_string("NORMAL");
        gattribute.data  = cgltf.accessors + cur_accessor++;
      }
      if (!shape.texcoords.empty()) {
        auto& gattribute = gprimitive.attributes[gprimitive.attributes_count++];
        gattribute.type  = cgltf_attribute_type_texcoord;
        gattribute.name  = copy_string("TEXCOORD_0");
        gattribute.data  = cgltf.accessors + cur_accessor++;
      }
      if (!shape.colors.empty()) {
        auto& gattribute = gprimitive.attributes[gprimitive.attributes_count++];
        gattribute.type  = cgltf_attribute_type_color;
        gattribute.name  = copy_string("COLOR_0");
        gattribute.data  = cgltf.accessors + cur_accessor++;
      }
      if (!shape.radius.empty()) {
        auto& gattribute = gprimitive.attributes[gprimitive.attributes_count++];
        gattribute.type  = cgltf_attribute_type_invalid;
        gattribute.name  = copy_string("RADIUS");
        gattribute.data  = cgltf.accessors + cur_accessor++;
      }
      if (is_points(shape)) {
        gprimitive.type    = cgltf_primitive_type_points;
        gprimitive.indices = cgltf.accessors + cur_accessor++;
      } else if (is_lines(shape)) {
        gprimitive.type    = cgltf_primitive_type_lines;
        gprimitive.indices = cgltf.accessors + cur_accessor++;
      } else if (is_triangles(shape)) {
        gprimitive.type    = cgltf_primitive_type_triangles;
        gprimitive.indices = cgltf.accessors + cur_accessor++;
      } else if (is_quads(shape)) {
        gprimitive.type    = cgltf_primitive_type_triangles;
        gprimitive.indices = cgltf.accessors + cur_accessor++;
      }
    }
  }

  // nodes
  if (!scene.cameras.empty() || !scene.instances.empty()) {
    alloc_array(cgltf.nodes_count, cgltf.nodes,
        scene.cameras.size() + scene.instances.size() + 1);
    for (auto idx : range(scene.cameras.size())) {
      auto& camera = scene.cameras[idx];
      auto& gnode  = cgltf.nodes[idx];
      gnode.name   = copy_string(get_camera_name(scene, camera));
      auto xform   = to_mat(camera.frame);
      memcpy(gnode.matrix, &xform, sizeof(mat4f));
      gnode.has_matrix = 1;
      gnode.camera     = cgltf.cameras + idx;
    }
    for (auto idx : range(scene.instances.size())) {
      auto& instance = scene.instances[idx];
      auto& gnode    = cgltf.nodes[idx + scene.cameras.size()];
      gnode.name     = copy_string(get_instance_name(scene, instance));
      auto xform     = to_mat(instance.frame);
      memcpy(gnode.matrix, &xform, sizeof(mat4f));
      gnode.has_matrix = 1;
      gnode.mesh       = cgltf.meshes +
                   mesh_map.at({instance.shape, instance.material});
    }
    // root children
    auto& groot = cgltf.nodes[cgltf.nodes_count - 1];
    groot.name  = copy_string("root");
    alloc_array(groot.children_count, groot.children, cgltf.nodes_count - 1);
    for (auto idx : range(cgltf.nodes_count - 1))
      groot.children[idx] = cgltf.nodes + idx;
    // scene
    alloc_array(cgltf.scenes_count, cgltf.scenes, 1);
    auto& gscene = cgltf.scenes[0];
    alloc_array(gscene.nodes_count, gscene.nodes, 1);
    gscene.nodes[0] = cgltf.nodes + cgltf.nodes_count - 1;
    cgltf.scene     = cgltf.scenes;
  }

  // save gltf
  auto options = cgltf_options{};
  memset(&options, 0, sizeof(options));
  options.memory.free = [](void*, void* ptr) { free(ptr); };
  cgltf.memory.free   = [](void*, void* ptr) { free(ptr); };
  auto result         = cgltf_write_file(&options, filename.c_str(), &cgltf);
  if (result != cgltf_result_success) {
    cgltf_free(&cgltf);
    throw io_error{"cannot save " + filename};
  }

  // cleanup
  cgltf_free(&cgltf);

  // get filename from name
  auto make_filename = [filename](const string& name, const string& group,
                           const string& extension) {
    return path_join(path_dirname(filename), group, name + extension);
  };

  // dirname
  auto dirname = path_dirname(filename);

  try {
    if (noparallel) {
      // save shapes
      for (auto& shape : scene.shapes) {
        auto path = "shapes/" + get_shape_name(scene, shape) + ".bin";
        save_binshape(path_join(dirname, path), shape);
      }
      // save textures
      for (auto& texture : scene.textures) {
        auto path = "textures/" + get_texture_name(scene, texture) + ".png";
        save_texture(path_join(dirname, path), texture);
      }
    } else {
      // save shapes
      parallel_foreach(scene.shapes, [&](auto& shape) {
        auto path = "shapes/" + get_shape_name(scene, shape) + ".bin";
        return save_binshape(path_join(dirname, path), shape);
      });
      // save textures
      parallel_foreach(scene.textures, [&](auto& texture) {
        auto path = "textures/" + get_texture_name(scene, texture) + ".png";
        return save_texture(path_join(dirname, path), texture);
      });
    }
  } catch (std::exception& except) {
    throw io_error(
        "cannot save " + filename + " since " + string(except.what()));
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF PBRT
// -----------------------------------------------------------------------------
namespace yocto {

// load pbrt scenes
static void load_pbrt_scene(
    const string& filename, scene_data& scene, bool noparallel) {
  // load pbrt
  auto pbrt = load_pbrt(filename);

  // convert cameras
  for (auto& pcamera : pbrt.cameras) {
    auto& camera  = scene.cameras.emplace_back();
    camera.frame  = to_math(pcamera.frame);
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
    if (to_math(pmaterial.emission) != vec3f{0, 0, 0}) {
      material.type = material_type::matte;
    }
    material.emission  = to_math(pmaterial.emission);
    material.color     = to_math(pmaterial.color);
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
    shape.positions = (const vector<vec3f>&)pshape.positions;
    shape.normals   = (const vector<vec3f>&)pshape.normals;
    shape.texcoords = (const vector<vec2f>&)pshape.texcoords;
    shape.triangles = (const vector<vec3i>&)pshape.triangles;
    for (auto& uv : shape.texcoords) uv = flip_v(uv);
    if (!pshape.instanced) {
      auto& instance    = scene.instances.emplace_back();
      instance.frame    = to_math(pshape.frame);
      instance.shape    = (int)scene.shapes.size() - 1;
      instance.material = pshape.material;
    } else {
      for (auto& frame : pshape.instances) {
        auto& instance    = scene.instances.emplace_back();
        instance.frame    = to_math(frame) * to_math(pshape.frame);
        instance.shape    = (int)scene.shapes.size() - 1;
        instance.material = pshape.material;
      }
    }
  }

  // convert environments
  for (auto& penvironment : pbrt.environments) {
    auto& environment        = scene.environments.emplace_back();
    environment.frame        = to_math(penvironment.frame);
    environment.emission     = to_math(penvironment.emission);
    environment.emission_tex = penvironment.emission_tex;
  }

  // lights
  for (auto& plight : pbrt.lights) {
    auto& shape = scene.shapes.emplace_back();
    shapes_paths.emplace_back();
    shape.triangles   = (const vector<vec3i>&)plight.area_triangles;
    shape.positions   = (const vector<vec3f>&)plight.area_positions;
    shape.normals     = (const vector<vec3f>&)plight.area_normals;
    auto& material    = scene.materials.emplace_back();
    material.emission = to_math(plight.area_emission);
    auto& instance    = scene.instances.emplace_back();
    instance.shape    = (int)scene.shapes.size() - 1;
    instance.material = (int)scene.materials.size() - 1;
    instance.frame    = to_math(plight.area_frame);
  }

  // dirname
  auto dirname = path_dirname(filename);

  try {
    if (noparallel) {
      // load shape
      for (auto& shape : scene.shapes) {
        auto& path = shapes_paths[&shape - &scene.shapes.front()];
        if (path.empty()) continue;
        load_shape(path_join(dirname, path), shape, true);
      }
      // load texture
      for (auto& texture : scene.textures) {
        auto& path = texture_paths[&texture - &scene.textures.front()];
        load_texture(path_join(dirname, path), texture);
      }
    } else {
      // load shapes
      parallel_foreach(scene.shapes, [&](auto& shape) {
        auto& path = shapes_paths[&shape - &scene.shapes.front()];
        if (path.empty()) return;
        load_shape(path_join(dirname, path), shape, true);
      });
      // load textures
      parallel_foreach(scene.textures, [&](auto& texture) {
        auto& path = texture_paths[&texture - &scene.textures.front()];
        return load_texture(path_join(dirname, path), texture);
      });
    }
  } catch (std::exception& except) {
    throw io_error(
        "cannot load " + filename + " since " + string(except.what()));
  }

  // fix scene
  add_missing_camera(scene);
  add_missing_radius(scene);
}

// Save a pbrt scene
static void save_pbrt_scene(
    const string& filename, const scene_data& scene, bool noparallel) {
  // save pbrt
  auto pbrt = pbrt_model{};

  // convert camera
  auto& camera       = scene.cameras.front();
  auto& pcamera      = pbrt.cameras.emplace_back();
  pcamera.frame      = to_array(camera.frame);
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
    pmaterial.emission  = to_array(material.emission);
    pmaterial.color     = to_array(material.color);
    pmaterial.roughness = material.roughness;
    pmaterial.ior       = material.ior;
    pmaterial.opacity   = material.opacity;
    pmaterial.color_tex = material.color_tex;
  }

  // convert instances
  for (auto& instance : scene.instances) {
    auto& pshape     = pbrt.shapes.emplace_back();
    pshape.filename_ = get_shape_name(scene, instance.shape) + ".ply";
    pshape.frame     = to_array(instance.frame);
    pshape.frend     = to_array(instance.frame);
    pshape.material  = instance.material;
  }

  // convert environments
  for (auto& environment : scene.environments) {
    auto& penvironment        = pbrt.environments.emplace_back();
    penvironment.emission     = to_array(environment.emission);
    penvironment.emission_tex = environment.emission_tex;
  }

  // save pbrt
  save_pbrt(filename, pbrt);

  // dirname
  auto dirname = path_dirname(filename);

  try {
    if (noparallel) {
      // save textures
      for (auto& shape : scene.shapes) {
        auto path = "shapes/" + get_shape_name(scene, shape) + ".ply";
        save_shape(path_join(dirname, path), shape, true);
      }
      // save shapes
      for (auto& texture : scene.textures) {
        auto path = "textures/" + get_texture_name(scene, texture) +
                    (!texture.pixelsf.empty() ? ".hdr" : ".png");
        save_texture(path_join(dirname, path), texture);
      }
    } else {
      // save shapes
      parallel_foreach(scene.shapes, [&](auto& shape) {
        auto path = "shapes/" + get_shape_name(scene, shape) + ".ply";
        return save_shape(path_join(dirname, path), shape, true);
      });
      // save textures
      parallel_foreach(scene.textures, [&](auto& texture) {
        auto path = "textures/" + get_texture_name(scene, texture) +
                    (!texture.pixelsf.empty() ? ".hdr"s : ".png"s);
        return save_texture(path_join(dirname, path), texture);
      });
    }
  } catch (std::exception& except) {
    throw io_error(
        "cannot save " + filename + " since " + string(except.what()));
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF MITSUBA
// -----------------------------------------------------------------------------
namespace yocto {

// load mitsuba scenes
static void load_mitsuba_scene(
    const string& filename, scene_data& scene, bool noparallel) {
  // not implemented
  throw io_error(
      "cannot load " + filename + " since format is not supported for reading");
}

// To xml helpers
static void xml_attribute(string& xml, const string& name, bool value) {
  xml += " " + name + "=\"" + (value ? string("true") : string("false")) + "\"";
}
static void xml_attribute(string& xml, const string& name, int value) {
  xml += " " + name + "=\"" + std::to_string(value) + "\"";
}
static void xml_attribute(string& xml, const string& name, float value) {
  xml += " " + name + "=\"" + std::to_string(value) + "\"";
}
static void xml_attribute(string& xml, const string& name, const char* value) {
  xml += " " + name + "=\"" + string(value) + "\"";
}
static void xml_attribute(
    string& xml, const string& name, const string& value) {
  xml += " " + name + "=\"" + value + "\"";
}
static void xml_attribute(string& xml, const string& name, const vec3f& value) {
  xml += " " + name + "=\"";
  for (auto v : value) {
    if (xml.back() != '\"') xml += ' ';
    xml += std::to_string(v);
  }
}
static void xml_attribute(
    string& xml, const string& name, const frame3f& value_) {
  xml += " " + name + "=\"";
  auto value = to_mat(value_);
  for (auto ij : range(vec2i(4, 4))) {
    xml += (xml.back() == '"' ? "" : " ") + std::to_string(value[ij.y][ij.x]);
  }
  xml += "\"";
}
template <typename T, typename... Ts>
static void xml_attributes(
    string& xml, const string& name, const T& value, const Ts&... attributes) {
  xml_attribute(xml, name, value);
  if constexpr (sizeof...(attributes) != 0) xml_attributes(xml, attributes...);
}
template <typename... Ts>
static void xml_element(string& xml, const string& indent, const string& name,
    const Ts&... attributes) {
  xml += indent + "<" + name;
  xml_attributes(xml, attributes...);
  xml += "/>\n";
}
template <typename... Ts>
static void xml_begin(
    string& xml, string& indent, const string& name, const Ts&... attributes) {
  xml += indent + "<" + name;
  xml_attributes(xml, attributes...);
  xml += ">\n";
  indent += "  ";
}
static void xml_end(string& xml, string& indent, const string& name) {
  indent = indent.substr(0, indent.size() - 2);
  xml += indent + "</" + name + ">\n";
}
template <typename T>
static void xml_default(
    string& xml, const string& indent, const string& name, const T& value) {
  xml_element(xml, indent, "default", "name", name, "value", value);
}
template <typename T>
static void xml_property(string& xml, const string& indent, const string& type,
    const string& name, const T& value, const string& ref) {
  if (ref.empty()) {
    if (name.empty()) {
      xml_element(xml, indent, type, "value", value);
    } else {
      xml_element(xml, indent, type, "name", name, "value", value);
    }
  } else {
    xml_element(xml, indent, type, "name", name, "value", ref);
  }
}
static void xml_property(string& xml, const string& indent, const string& name,
    int value, const string& ref = "") {
  xml_property(xml, indent, "integer", name, value, ref);
}
static void xml_property(string& xml, const string& indent, const string& name,
    float value, const string& ref = "") {
  xml_property(xml, indent, "float", name, value, ref);
}
static void xml_property(string& xml, const string& indent, const string& name,
    bool value, const string& ref = "") {
  xml_property(xml, indent, "boolean", name, value, ref);
}
static void xml_property(string& xml, const string& indent, const string& name,
    const char* value, const string& ref = "") {
  xml_property(xml, indent, "string", name, value, ref);
}
static void xml_property(string& xml, const string& indent, const string& name,
    const string& value, const string& ref = "") {
  xml_property(xml, indent, "string", name, value, ref);
}
static void xml_property(string& xml, const string& indent, const string& name,
    const frame3f& value, const string& ref = "") {
  xml_property(xml, indent, "matrix", name, value, ref);
}
static void xml_property(string& xml, const string& indent, const string& name,
    const vec3f& value, const string& ref = "") {
  xml_property(xml, indent, "rgb", name, value, ref);
}

// Save a mitsuba scene
static void save_mitsuba_scene(
    const string& filename, const scene_data& scene, bool noparallel) {
  // write xml directly
  auto xml = string{};

  // begin
  auto indent = string{};
  xml_begin(xml, indent, "scene", "version", "3.0.0");

  // defaults
  xml_default(xml, indent, "integrator", "path");
  xml_default(xml, indent, "spp", 64);
  xml_default(xml, indent, "resx", 1440);
  xml_default(xml, indent, "resy", 720);
  xml_default(xml, indent, "pixel_format", "rgb");
  xml_default(xml, indent, "max_depth", 8);
  xml_default(xml, indent, "rr_depth", 64);

  // integrator
  xml_begin(xml, indent, "integrator", "type", "$integrator");
  xml_property(xml, indent, "max_depth", 0, "$max_depth");
  xml_property(xml, indent, "rr_depth", 0, "$rr_depth");
  xml_property(xml, indent, "hide_emitters", false);
  xml_end(xml, indent, "integrator");

  // film
  xml_begin(xml, indent, "film", "type", "hdrfilm", "id", "film");
  xml_property(xml, indent, "width", 0, "$resx");
  xml_property(xml, indent, "height", 0, "$resy");
  xml_element(xml, indent, "rfilter", "type", "box");
  xml_property(xml, indent, "pixel_format", "", "$pixel_format");
  xml_end(xml, indent, "film");

  // sampler
  xml_begin(xml, indent, "sampler", "type", "independent", "id", "sampler");
  xml_property(xml, indent, "sample_count", 0, "$spp");
  xml_end(xml, indent, "sampler");

  // camera
  auto& camera = scene.cameras.at(0);
  xml_begin(xml, indent, "sensor", "type", "perspective");
  xml_property(xml, indent, "fov_axis", "smaller");
  xml_property(xml, indent, "fov", 20.0f);
  xml_begin(xml, indent, "transform", "name", "to_world");
  xml_element(xml, indent, "lookat", "origin", camera.frame.o, "target",
      camera.frame.o - camera.frame.z, "up", vec3f{0, 1, 0});
  xml_end(xml, indent, "transform");
  xml_element(xml, indent, "ref", "id", "sampler");
  xml_element(xml, indent, "ref", "id", "film");
  xml_end(xml, indent, "sensor");

  // textures
  auto tid = 0;
  for (auto& texture : scene.textures) {
    if (texture.pixelsf.empty()) {
      xml_begin(xml, indent, "texture", "type", "bitmap", "id",
          "texture" + std::to_string(tid));
      xml_property(xml, indent, "filename",
          "textures/" + get_texture_name(scene, texture) +
              (texture.pixelsf.empty() ? ".png" : ".hdr"));
      xml_end(xml, indent, "texture");
    }
    tid += 1;
  }

  // environments
  for (auto& environment : scene.environments) {
    if (environment.emission_tex != invalidid) {
      auto& texture = scene.textures.at(environment.emission_tex);
      xml_begin(xml, indent, "emitter", "type", "envmap");
      xml_property(xml, indent, "scale", mean(environment.emission));
      xml_property(xml, indent, "filename",
          "textures/" + get_texture_name(scene, texture) + ".hdr");
      xml_end(xml, indent, "emitter");
    } else {
      xml_begin(xml, indent, "emitter", "type", "constant");
      xml_property(xml, indent, "rgb", "radiance", environment.emission, "");
      xml_end(xml, indent, "emitter");
    }
  }

  // property or texture
  auto xml_property_or_texture = [](string& xml, const string& indent,
                                     const string& name, auto value,
                                     int texture) {
    if (texture == invalidid) {
      xml_property(xml, indent, "rgb", name, value, "");
    } else {
      xml_element(xml, indent, "ref", "id", "texture" + std::to_string(texture),
          "name", name);
    }
  };

  // materials
  auto mid = 0;
  for (auto& material : scene.materials) {
    // xml_begin(xml, indent, "bsdf", "type", "twosided", "id",
    //     "material" + std::to_string(mid));
    switch (material.type) {
      case material_type::matte: {
        xml_begin(xml, indent, "bsdf", "type", "diffuse", "id",
            "material" + std::to_string(mid));
        xml_property_or_texture(
            xml, indent, "reflectance", material.color, material.color_tex);
        xml_end(xml, indent, "bsdf");
      } break;
      case material_type::reflective: {
        xml_begin(xml, indent, "bsdf", "type",
            material.roughness < 0.03f ? "conductor" : "roughconductor", "id",
            "material" + std::to_string(mid));
        xml_property_or_texture(
            xml, indent, "eta", reflectivity_to_eta(material.color), invalidid);
        xml_property_or_texture(xml, indent, "k", vec3f{0, 0, 0}, invalidid);
        if (material.roughness >= 0.03f) {
          xml_property(
              xml, indent, "alpha", material.roughness * material.roughness);
        }
        xml_end(xml, indent, "bsdf");
      } break;
      case material_type::glossy: {
        xml_begin(xml, indent, "bsdf", "type",
            material.roughness < 0.03f ? "plastic" : "roughplastic", "id",
            "material" + std::to_string(mid));
        xml_property_or_texture(xml, indent, "diffuse_reflectance",
            material.color, material.color_tex);
        if (material.roughness >= 0.03f) {
          xml_property(
              xml, indent, "alpha", material.roughness * material.roughness);
        }
        xml_end(xml, indent, "bsdf");
      } break;
      case material_type::transparent: {
        xml_begin(xml, indent, "bsdf", "type",
            material.roughness < 0.03f ? "conductor" : "roughconductor", "id",
            "material" + std::to_string(mid));
        xml_property_or_texture(
            xml, indent, "eta", reflectivity_to_eta(material.color), invalidid);
        xml_property_or_texture(xml, indent, "k", vec3f{0, 0, 0}, invalidid);
        if (material.roughness >= 0.03f) {
          xml_property_or_texture(
              xml, indent, "alpha", material.roughness, invalidid);
        }
        xml_end(xml, indent, "bsdf");
      } break;
      case material_type::volumetric:
      case material_type::subsurface:
      case material_type::refractive: {
        xml_begin(xml, indent, "bsdf", "type",
            material.roughness < 0.03f ? "dielectric" : "roughdielectric", "id",
            "material" + std::to_string(mid));
        xml_property(xml, indent, "int_ior", 1.5f);
        if (material.roughness >= 0.03f) {
          xml_property_or_texture(
              xml, indent, "alpha", material.roughness, invalidid);
        }
        xml_end(xml, indent, "bsdf");
        if (material.color != vec3f{1, 1, 1}) {
          xml_begin(xml, indent, "medium", "type", "homogeneous", "id",
              "medium" + std::to_string(mid));
          xml_property(xml, indent, "albedo", material.scattering);
          auto density = -log(clamp(material.color, 0.0001f, 1.0f)) /
                         material.trdepth;
          xml_property(xml, indent, "sigma_t", mean(density));
          xml_end(xml, indent, "medium");
        }
      } break;
      case material_type::gltfpbr: {
        xml_begin(xml, indent, "bsdf", "type", "diffuse", "id",
            "material" + std::to_string(mid));
        xml_property_or_texture(
            xml, indent, "reflectance", material.color, material.color_tex);
        xml_end(xml, indent, "bsdf");
      } break;
    }
    // todo
    // xml_end(xml, indent, "bsdf");
    mid += 1;
  }

  // flattened instances
  for (auto& instance : scene.instances) {
    auto& shape    = scene.shapes[instance.shape];
    auto& material = scene.materials[instance.material];
    xml_begin(xml, indent, "shape", "type", "ply");
    xml_property(xml, indent, "filename",
        "shapes/" + get_shape_name(scene, shape) + ".ply");
    if (instance.frame != frame3f{}) {
      xml_begin(xml, indent, "transform", "name", "to_world");
      xml_property(xml, indent, "", instance.frame);
      xml_end(xml, indent, "transform");
    }
    if (material.emission != vec3f{0, 0, 0}) {
      xml_property(xml, indent, "flip_normals", true);
      xml_begin(xml, indent, "emitter", "type", "area");
      xml_property(xml, indent, "rgb", "radiance", material.emission, "");
      xml_end(xml, indent, "emitter");
    }
    xml_element(xml, indent, "ref", "id",
        "material" + std::to_string(instance.material));
    if (material.type == material_type::refractive &&
        material.color != vec3f{1, 1, 1}) {
      xml_element(xml, indent, "ref", "name", "interior", "id",
          "medium" + std::to_string(instance.material));
    }
    xml_end(xml, indent, "shape");
  }

  // close
  xml_end(xml, indent, "scene");

  // save xml
  save_text(filename, xml);

  // dirname
  auto dirname = path_dirname(filename);

  auto triangulate = [](const shape_data& shape) -> shape_data {
    if (shape.quads.empty()) return shape;
    auto tshape      = shape;
    tshape.triangles = quads_to_triangles(shape.quads);
    tshape.quads.clear();
    return tshape;
  };

  try {
    if (noparallel) {
      // save shapes
      for (auto& shape : scene.shapes) {
        auto path = "shapes/" + get_shape_name(scene, shape) + ".ply";
        save_shape(path_join(dirname, path), triangulate(shape), true);
      }
      // save textures
      for (auto& texture : scene.textures) {
        auto path = "textures/" + get_texture_name(scene, texture) +
                    (!texture.pixelsf.empty() ? ".hdr" : ".png");
        save_texture(path_join(dirname, path), texture);
      }
    } else {
      // save shapes
      parallel_foreach(scene.shapes, [&](auto& shape) {
        auto path = "shapes/" + get_shape_name(scene, shape) + ".ply";
        return save_shape(path_join(dirname, path), triangulate(shape), true);
      });
      // save textures
      parallel_foreach(scene.textures, [&](auto& texture) {
        auto path = "textures/" + get_texture_name(scene, texture) +
                    (!texture.pixelsf.empty() ? ".hdr"s : ".png"s);
        return save_texture(path_join(dirname, path), texture);
      });
    }
  } catch (std::exception& except) {
    throw io_error(
        "cannot save " + filename + " since " + string(except.what()));
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PARAMETER IO
// -----------------------------------------------------------------------------
namespace yocto {

// Conversion to/from json
template <typename T>
static void to_json(
    json_value& json, const T& value, const vector<pair<T, string>>& labels) {
  auto it = std::find_if(labels.begin(), labels.end(),
      [value](const pair<T, string>& kv) -> bool { return kv.first == value; });
  if (it == labels.end()) throw std::invalid_argument{"bad enum value"};
  json = it->second;
}
template <typename T>
static void from_json(
    const json_value& json, T& value, const vector<pair<T, string>>& labels) {
  auto label = json.get<string>();
  auto it    = std::find_if(labels.begin(), labels.end(),
         [&label](
          const pair<T, string>& kv) -> bool { return kv.second == label; });
  // TODO: better error type
  if (it == labels.end()) throw std::invalid_argument{"bad enum value"};
  value = it->first;
}

// Conversion to/from json
static void to_json(json_value& json, const trace_sampler_type& value) {
  to_json(json, value, trace_sampler_labels);
}
static void from_json(const json_value& json, trace_sampler_type& value) {
  from_json(json, value, trace_sampler_labels);
}
static void to_json(json_value& json, const trace_falsecolor_type& value) {
  to_json(json, value, trace_falsecolor_labels);
}
static void from_json(const json_value& json, trace_falsecolor_type& value) {
  from_json(json, value, trace_falsecolor_labels);
}

// Conversion to/from json
static void to_json(json_value& json, const trace_params& value) {
  json["camera"]         = value.camera;
  json["resolution"]     = value.resolution;
  json["sampler"]        = value.sampler;
  json["falsecolor"]     = value.falsecolor;
  json["samples"]        = value.samples;
  json["bounces"]        = value.bounces;
  json["clamp"]          = value.clamp;
  json["nocaustics"]     = value.nocaustics;
  json["envhidden"]      = value.envhidden;
  json["tentfilter"]     = value.tentfilter;
  json["seed"]           = value.seed;
  json["embreebvh"]      = value.embreebvh;
  json["highqualitybvh"] = value.highqualitybvh;
  json["noparallel"]     = value.noparallel;
  json["pratio"]         = value.pratio;
  json["denoise"]        = value.denoise;
  json["batch"]          = value.batch;
}
static void from_json(const json_value& json, trace_params& value) {
  value.camera         = json.value("camera", value.camera);
  value.resolution     = json.value("resolution", value.resolution);
  value.sampler        = json.value("sampler", value.sampler);
  value.falsecolor     = json.value("falsecolor", value.falsecolor);
  value.samples        = json.value("samples", value.samples);
  value.bounces        = json.value("bounces", value.bounces);
  value.clamp          = json.value("clamp", value.clamp);
  value.nocaustics     = json.value("nocaustics", value.nocaustics);
  value.envhidden      = json.value("envhidden", value.envhidden);
  value.tentfilter     = json.value("tentfilter", value.tentfilter);
  value.seed           = json.value("seed", value.seed);
  value.embreebvh      = json.value("embreebvh", value.embreebvh);
  value.highqualitybvh = json.value("highqualitybvh", value.highqualitybvh);
  value.noparallel     = json.value("noparallel", value.noparallel);
  value.pratio         = json.value("pratio", value.pratio);
  value.denoise        = json.value("denoise", value.denoise);
  value.batch          = json.value("batch", value.batch);
}

// Conversion to/from json
static void to_json(json_value& json, const colorgrade_params& value) {
  json["exposure"]         = value.exposure;
  json["tint"]             = value.tint;
  json["lincontrast"]      = value.lincontrast;
  json["logcontrast"]      = value.logcontrast;
  json["linsaturation"]    = value.linsaturation;
  json["filmic"]           = value.filmic;
  json["srgb"]             = value.srgb;
  json["contrast"]         = value.contrast;
  json["saturation"]       = value.saturation;
  json["shadows"]          = value.shadows;
  json["midtones"]         = value.midtones;
  json["highlights"]       = value.highlights;
  json["shadows_color"]    = value.shadows_color;
  json["midtones_color"]   = value.midtones_color;
  json["highlights_color"] = value.highlights_color;
}
static void from_json(const json_value& json, colorgrade_params& value) {
  value.exposure         = json.value("exposure", value.exposure);
  value.tint             = json.value("tint", value.tint);
  value.lincontrast      = json.value("lincontrast", value.lincontrast);
  value.logcontrast      = json.value("logcontrast", value.logcontrast);
  value.linsaturation    = json.value("linsaturation", value.linsaturation);
  value.filmic           = json.value("filmic", value.filmic);
  value.srgb             = json.value("srgb", value.srgb);
  value.contrast         = json.value("contrast", value.contrast);
  value.saturation       = json.value("saturation", value.saturation);
  value.shadows          = json.value("shadows", value.shadows);
  value.midtones         = json.value("midtones", value.midtones);
  value.highlights       = json.value("highlights", value.highlights);
  value.shadows_color    = json.value("shadows_color", value.shadows_color);
  value.midtones_color   = json.value("midtones_color", value.midtones_color);
  value.highlights_color = json.value(
      "highlights_color", value.highlights_color);
}

// Load/Save/Update trace params
template <typename Params>
static Params load_params(const string& filename) {
  auto params = Params{};
  load_params(filename, params);
  return params;
}
template <typename Params>
static void load_params(const string& filename, Params& params) {
  // json
  auto json = load_json(filename);

  // conversion
  try {
    params = {};
    from_json(json, params);
  } catch (...) {
    throw io_error{"error parsing params"};
  }
}
template <typename Params>
static void update_params(const string& filename, Params& params) {
  // json
  auto json = load_json(filename);

  // conversion
  try {
    from_json(json, params);
  } catch (...) {
    throw io_error{"error parsing params"};
  }
}
template <typename Params>
static void save_params(const string& filename, const Params& params) {
  auto json = json_value{};
  to_json(json, params);
  return save_json(filename, json);
}

// Load/Save/Update trace params
trace_params load_trace_params(const string& filename) {
  return load_params<trace_params>(filename);
}
void load_trace_params(const string& filename, trace_params& params) {
  return load_params(filename, params);
}
void update_trace_params(const string& filename, trace_params& params) {
  return update_params(filename, params);
}
void save_trace_params(const string& filename, const trace_params& params) {
  return save_params(filename, params);
}

// Load/Save/Update color grade params
colorgrade_params load_colorgrade_params(const string& filename) {
  return load_params<colorgrade_params>(filename);
}
void load_colorgrade_params(const string& filename, colorgrade_params& params) {
  return load_params(filename, params);
}
void update_colorgrade_params(
    const string& filename, colorgrade_params& params) {
  return update_params(filename, params);
}
void save_colorgrade_params(
    const string& filename, const colorgrade_params& params) {
  return save_params(filename, params);
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
  if (pos >= (int)args.size() || args[pos].find("--") == 0) {
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
  if (pos >= (int)args.size()) return true;
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
// HELPERS FOR JSON MANIPULATION
// -----------------------------------------------------------------------------
namespace yocto {

// Validate a Json value againt a schema. Returns the first error found.
void validate_json(const json_value& json, const json_value& schema);
bool validate_json(
    const json_value& json, const json_value& schema, string& error);

// Converts command line arguments to Json. Never throws since a conversion
// is always possible in our conventions. Validation is done using a schema.
json_value make_json_cli(const vector<string>& args) {
  // init json
  auto json = json_value{};
  if (args.size() < 2) return json;

  // split into commans and options; use spans when available for speed
  auto commands = vector<string>{};
  auto options  = vector<pair<string, vector<string>>>{};
  for (auto& arg : args) {
    if (arg.find("--") == 0) {
      // start option
      options.push_back({arg.substr(2), {}});
    } else if (!options.empty()) {
      // add value
      options.back().second.push_back(arg);
    } else {
      // add command
      commands.push_back(arg);
    }
  }

  // build commands
  auto current = &json;
  for (auto& command : commands) {
    auto& json      = *current;
    json["command"] = command;
    json[command]   = json_value::object();
    current         = &json[command];
  }

  // build options
  for (auto& [name, values] : options) {
    auto& json = *current;
    if (values.empty()) {
      json[name] = true;
    } else {
      json[name] = json_value::array();
      for (auto& value : values) {
        if (value == "true") {
          json[name].push_back(true);
        } else if (value == "false") {
          json[name].push_back(false);
        } else if (value == "null") {
          json[name].push_back(nullptr);
        } else if ((bool)std::isdigit((int)value[0]) || value[0] == '-' ||
                   value[0] == '+') {
          try {
            if (value.find('.') != string::npos) {
              json[name].push_back(std::stod(value));
            } else if (value.find('-') == 0) {
              json[name].push_back(std::stoll(value));
            } else {
              json[name].push_back(std::stoull(value));
            }
          } catch (...) {
            json[name].push_back(value);
          }
        } else {
          json[name].push_back(value);
        }
      }
      if (values.size() == 1) {
        json[name] = json[name].front();
      }
    }
  }

  // done
  return json;
}
json_value make_json_cli(int argc, const char** argv) {
  return make_json_cli(vector<string>{argv, argv + argc});
}

// Validates a JSON against a schema including CLI constraints.
json_value validate_json_cli(const vector<string>& args);
json_value validate_json_cli(int argc, const char** argv);

// Helpers for creating schemas

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
    error = "cannot open " + filename;
    return false;
  };
  auto parse_error = [filename, &error]() {
    error = "cannot parse " + filename;
    return false;
  };
  auto read_error = [filename, &error]() {
    error = "cannot read " + filename;
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
  if (!read_values(fs, voxels.data(), nvalues)) throw io_error("cannot read " + filename);

  // done
  return true;
}

// save pfm
static bool save_yvol(const string& filename, int width, int height, int depth,
    int components, const vector<float>& voxels, string& error) {
  // error helpers
  auto open_error = [filename, &error]() {
    error = "cannot create " + filename;
    return false;
  };
  auto write_error = [filename, &error]() {
    error = "cannot read " + filename;
    return false;
  };

  auto fs       = fopen_utf8(filename.c_str(), "wb");
  auto fs_guard = unique_ptr<FILE, int (*)(FILE*)>(fs, &fclose);
  if (!fs) return open_error();

  if (!write_text(fs, "YVOL\n")) throw io_error("cannot write " + filename);
  if (!write_text(fs, std::to_string(width) + " " + std::to_string(height) +
                          " " + std::to_string(depth) + " " +
                          std::to_string(components) + "\n"))
    throw io_error("cannot write " + filename);
  auto nvalues = (size_t)width * (size_t)height * (size_t)depth *
                 (size_t)components;
  if (!write_values(fs, voxels.data(), nvalues)) throw io_error("cannot write " + filename);
  return true;
}

// Loads volume data from binary format.
bool load_volume(const string& filename, volume<float>& vol, string& error) {
  auto read_error = [filename, &error]() {
    error = "cannot read " + filename;
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
  auto filtered     = image{img.imsize(), vec4f{0,0,0,0}};
  auto filter_width = (int)ceil(2.57f * spatial_sigma);
  auto sw           = 1 / (2.0f * spatial_sigma * spatial_sigma);
  auto rw           = 1 / (2.0f * range_sigma * range_sigma);
  auto fw           = vector<float>();
  for (auto feature_sigma : features_sigma)
    fw.push_back(1 / (2.0f * feature_sigma * feature_sigma));
  for (auto j : range(img.height())) {
    for (auto i : range(img.width())) {
      auto av = vec4f{0,0,0,0};
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
          for (auto fi : range(features.size())) {
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
  auto filtered = image{img.imsize(), vec4f{0,0,0,0}};
  auto fwidth   = (int)ceil(2.57f * spatial_sigma);
  auto sw       = 1 / (2.0f * spatial_sigma * spatial_sigma);
  auto rw       = 1 / (2.0f * range_sigma * range_sigma);
  for (auto j : range(img.height())) {
    for (auto i : range(img.width())) {
      auto av = vec4f{0,0,0,0};
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
