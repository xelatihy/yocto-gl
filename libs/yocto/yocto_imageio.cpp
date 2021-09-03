//
// Implementation for Yocto/Image Input and Output functions.
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

#include "yocto_imageio.h"

#include <cassert>
#include <cctype>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <memory>

#include "ext/stb_image.h"
#include "ext/stb_image_resize.h"
#include "ext/stb_image_write.h"
#include "ext/tinyexr.h"
#include "yocto_color.h"
#include "yocto_commonio.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::unique_ptr;
using namespace std::string_literals;

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
image_data load_image(const string& filename) {
  auto image = image_data{};
  load_image(filename, image);
  return image;
}
void load_image(const string& filename, image_data& image) {
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
    auto buffer = load_binary(filename);
    auto pixels = (float*)nullptr;
    if (LoadEXRFromMemory(&pixels, &image.width, &image.height, buffer.data(),
            buffer.size(), nullptr) != 0)
      throw io_error::read_error(filename);
    image.linear = true;
    image.pixels = from_linear(pixels, image.width, image.height);
    free(pixels);
  } else if (ext == ".hdr" || ext == ".HDR") {
    auto buffer = load_binary(filename);
    auto ncomp  = 0;
    auto pixels = stbi_loadf_from_memory(buffer.data(), (int)buffer.size(),
        &image.width, &image.height, &ncomp, 4);
    if (!pixels) throw io_error::read_error(filename);
    image.linear = true;
    image.pixels = from_linear(pixels, image.width, image.height);
    free(pixels);
  } else if (ext == ".png" || ext == ".PNG") {
    auto buffer = load_binary(filename);
    auto ncomp  = 0;
    auto pixels = stbi_load_from_memory(buffer.data(), (int)buffer.size(),
        &image.width, &image.height, &ncomp, 4);
    if (!pixels) throw io_error::read_error(filename);
    image.linear = false;
    image.pixels = from_srgb(pixels, image.width, image.height);
    free(pixels);
  } else if (ext == ".jpg" || ext == ".JPG" || ext == ".jpeg" ||
             ext == ".JPEG") {
    auto buffer = load_binary(filename);
    auto ncomp  = 0;
    auto pixels = stbi_load_from_memory(buffer.data(), (int)buffer.size(),
        &image.width, &image.height, &ncomp, 4);
    if (!pixels) throw io_error::read_error(filename);
    image.linear = false;
    image.pixels = from_srgb(pixels, image.width, image.height);
    free(pixels);
  } else if (ext == ".tga" || ext == ".TGA") {
    auto buffer = load_binary(filename);
    auto ncomp  = 0;
    auto pixels = stbi_load_from_memory(buffer.data(), (int)buffer.size(),
        &image.width, &image.height, &ncomp, 4);
    if (!pixels) throw io_error::read_error(filename);
    image.linear = false;
    image.pixels = from_srgb(pixels, image.width, image.height);
    free(pixels);
  } else if (ext == ".bmp" || ext == ".BMP") {
    auto buffer = load_binary(filename);
    auto ncomp  = 0;
    auto pixels = stbi_load_from_memory(buffer.data(), (int)buffer.size(),
        &image.width, &image.height, &ncomp, 4);
    if (!pixels) throw io_error::read_error(filename);
    image.linear = false;
    image.pixels = from_srgb(pixels, image.width, image.height);
    free(pixels);
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    // create preset
    image = make_image_preset(path_basename(filename));
  } else {
    throw io_error::format_error(filename);
  }
}

// Saves an hdr image.
void save_image(const string& filename, const image_data& image) {
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
      throw io_error::write_error(filename);
    save_binary(filename, buffer);
  } else if (ext == ".exr" || ext == ".EXR") {
    auto data = (byte*)nullptr;
    auto size = (size_t)0;
    if (SaveEXRToMemory((const float*)to_linear(image).data(), (int)image.width,
            (int)image.height, 4, 1, &data, &size, nullptr) < 0)
      throw io_error::write_error(filename);
    auto buffer = vector<byte>{data, data + size};
    free(data);
    save_binary(filename, buffer);
  } else if (ext == ".png" || ext == ".PNG") {
    auto buffer = vector<byte>{};
    if (!stbi_write_png_to_func(stbi_write_data, &buffer, (int)image.width,
            (int)image.height, 4, (const byte*)to_srgb(image).data(),
            (int)image.width * 4))
      throw io_error::write_error(filename);
    save_binary(filename, buffer);
  } else if (ext == ".jpg" || ext == ".JPG" || ext == ".jpeg" ||
             ext == ".JPEG") {
    auto buffer = vector<byte>{};
    if (!stbi_write_jpg_to_func(stbi_write_data, &buffer, (int)image.width,
            (int)image.height, 4, (const byte*)to_srgb(image).data(), 75))
      throw io_error::write_error(filename);
    save_binary(filename, buffer);
  } else if (ext == ".tga" || ext == ".TGA") {
    auto buffer = vector<byte>{};
    if (!stbi_write_tga_to_func(stbi_write_data, &buffer, (int)image.width,
            (int)image.height, 4, (const byte*)to_srgb(image).data()))
      throw io_error::write_error(filename);
    save_binary(filename, buffer);
  } else if (ext == ".bmp" || ext == ".BMP") {
    auto buffer = vector<byte>{};
    if (!stbi_write_bmp_to_func(stbi_write_data, &buffer, (int)image.width,
            (int)image.height, 4, (const byte*)to_srgb(image).data()))
      throw io_error::write_error(filename);
    save_binary(filename, buffer);
  } else {
    throw io_error::format_error(filename);
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
    throw io_error::preset_error(type_);
  }
}

// Loads/saves an image. Chooses hdr or ldr based on file name.
bool load_image(const string& filename, image_data& image, string& error) {
  try {
    load_image(filename, image);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// Saves an hdr image.
bool save_image(
    const string& filename, const image_data& image, string& error) {
  try {
    save_image(filename, image);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

bool make_image_preset(image_data& image, const string& type, string& error) {
  try {
    image = make_image_preset(type);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
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
