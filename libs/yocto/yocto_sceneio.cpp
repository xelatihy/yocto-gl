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
#include <memory>
#include <stdexcept>
#include <unordered_map>

#include "ext/stb_image.h"
#include "ext/stb_image_resize.h"
#include "ext/stb_image_write.h"
#include "ext/tinyexr.h"
#include "yocto_color.h"
#include "yocto_commonio.h"
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
// IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

template <typename T>
static T _load_pfm_swap_endian(T value) {
  // https://stackoverflow.com/questions/105252/how-do-i-convert-between-big-endian-and-little-endian-values-in-c
  static_assert(sizeof(char) == 1, "sizeof(char) == 1");
  union {
    T             value;
    unsigned char bytes[sizeof(T)];
  } source, dest;
  source.value = value;
  for (auto k = (size_t)0; k < sizeof(T); k++)
    dest.bytes[k] = source.bytes[sizeof(T) - k - 1];
  return dest.value;
}

// Split a string
static vector<string> _load_pfm_split_string(const string& str) {
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
}

// Pfm load
static float* load_pfm(
    const string& filename, int* width, int* height, int* components, int req) {
  auto fs       = fopen_utf8(filename, "rb");
  auto fs_guard = unique_ptr<FILE, int (*)(FILE*)>(fs, &fclose);
  if (!fs) return nullptr;

  // buffer
  auto buffer = array<char, 4096>{};
  auto toks   = vector<string>();

  // read magic
  if (!fgets(buffer.data(), (int)buffer.size(), fs)) return nullptr;
  toks = _load_pfm_split_string(buffer.data());
  if (toks[0] == "Pf") {
    *components = 1;
  } else if (toks[0] == "PF") {
    *components = 3;
  } else {
    return nullptr;
  }

  // read width, height
  if (!fgets(buffer.data(), (int)buffer.size(), fs)) return nullptr;
  toks    = _load_pfm_split_string(buffer.data());
  *width  = atoi(toks[0].c_str());
  *height = atoi(toks[1].c_str());

  // read scale
  if (!fgets(buffer.data(), (int)buffer.size(), fs)) return nullptr;
  toks   = _load_pfm_split_string(buffer.data());
  auto s = (float)atof(toks[0].c_str());

  // read the data (flip y)
  auto npixels = (size_t)*width * (size_t)*height;
  auto nvalues = npixels * (size_t)*components;
  auto nrow    = (size_t)*width * (size_t)*components;
  auto pixels  = unique_ptr<float[]>{new float[nvalues]};
  for (auto j = *height - 1; j >= 0; j--) {
    if (fread(pixels.get() + j * nrow, 4, nrow, fs) != nrow) return nullptr;
  }

  // endian conversion
  if (s > 0) {
    for (auto i = (size_t)0; i < nvalues; ++i) {
      pixels[i] = _load_pfm_swap_endian(pixels[i]);
    }
  }

  // scale
  auto scl = (s > 0) ? s : -s;
  if (scl != 1) {
    for (auto i = (size_t)0; i < nvalues; i++) pixels[i] *= scl;
  }

  // check convertions
  if (req == 0 || *components == req) return pixels.release();

  // convert channels
  auto cpixels = unique_ptr<float[]>{new float[npixels * req]};
  if (req == 1) {
    if (*components == 3) {
      for (auto i = (size_t)0; i < npixels; i++) {
        cpixels[i] =
            (pixels[i * 3 + 0] + pixels[i * 3 + 1] + pixels[i * 3 + 2]) / 3;
      }
    }
  } else if (req == 2) {
    if (*components == 1) {
      for (auto i = (size_t)0; i < npixels; i++) {
        cpixels[i * 2 + 0] = pixels[i];
        cpixels[i * 2 + 1] = 1;
      }
    }
    if (*components == 3) {
      for (auto i = (size_t)0; i < npixels; i++) {
        cpixels[i * 2 + 0] =
            (pixels[i * 3 + 0] + pixels[i * 3 + 1] + pixels[i * 3 + 2]) / 3;
        cpixels[i * 2 + 1] = 1;
      }
    }
  } else if (req == 3) {
    if (*components == 1) {
      for (auto i = (size_t)0; i < npixels; i++) {
        cpixels[i * 3 + 0] = pixels[i];
        cpixels[i * 3 + 1] = pixels[i];
        cpixels[i * 3 + 2] = pixels[i];
      }
    }
  } else if (req == 4) {
    if (*components == 1) {
      for (auto i = (size_t)0; i < npixels; i++) {
        cpixels[i * 4 + 0] = pixels[i];
        cpixels[i * 4 + 1] = pixels[i];
        cpixels[i * 4 + 2] = pixels[i];
        cpixels[i * 4 + 3] = 1;
      }
    }
    if (*components == 3) {
      for (auto i = (size_t)0; i < npixels; i++) {
        cpixels[i * 4 + 0] = pixels[i * 3 + 0];
        cpixels[i * 4 + 1] = pixels[i * 3 + 1];
        cpixels[i * 4 + 2] = pixels[i * 3 + 2];
        cpixels[i * 4 + 3] = 1;
      }
    }
  } else {
    return nullptr;
  }

  // done
  return cpixels.release();
}

// save pfm
[[maybe_unused]] static bool save_pfm(const char* filename, int width,
    int height, int components, const float* pixels) {
  auto fs       = fopen_utf8(filename, "wb");
  auto fs_guard = unique_ptr<FILE, int (*)(FILE*)>(fs, &fclose);
  if (!fs) return false;

  if (fprintf(fs, "%s\n", (components == 1) ? "Pf" : "PF") < 0) return false;
  if (fprintf(fs, "%d %d\n", width, height) < 0) return false;
  if (fprintf(fs, "-1\n") < 0) return false;
  if (components == 1 || components == 3) {
    for (auto j = height - 1; j >= 0; j--) {
      if (fwrite(pixels + j * width * components, 4, width * components, fs) !=
          width * components)
        return false;
    }
  } else {
    for (auto j = height - 1; j >= 0; j--) {
      for (auto i = 0; i < width; i++) {
        auto vz = 0.0f;
        auto v  = pixels + (j * width + i) * components;
        if (fwrite(&v[0], 4, 1, fs) != 1) return false;
        if (fwrite(&v[1], 4, 1, fs) != 1) return false;
        if (components == 2) {
          if (fwrite(&vz, 4, 1, fs) != 1) return false;
        } else {
          if (fwrite(&v[2], 4, 1, fs) != 1) return false;
        }
      }
    }
  }

  return true;
}

// save pfm
static bool save_pfm_to_func(
    void (*func)(void* context, const void* data, int size), void* context,
    int width, int height, int components, const float* pixels) {
  auto buffer = array<char, 512>{};
  func(context, (components == 1) ? "Pf\n" : "PF\n", 3);
  func(context, buffer.data(),
      snprintf(buffer.data(), buffer.size(), "%d %d\n", width, height));
  func(context, "-1\n", 3);
  if (components == 1 || components == 3) {
    for (auto j = height - 1; j >= 0; j--) {
      func(context, pixels + j * width * components, 4 * width * components);
    }
  } else {
    for (auto j = height - 1; j >= 0; j--) {
      for (auto i = 0; i < width; i++) {
        auto vz = 0.0f;
        auto v  = pixels + (j * width + i) * components;
        func(context, &v[0], 4);
        func(context, &v[1], 4);
        if (components == 2) {
          func(context, &vz, 4);
        } else {
          func(context, &v[2], 4);
        }
      }
    }
  }

  return true;
}

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
  } else if (ext == ".pfm" || ext == ".PFM") {
    auto ncomp  = 0;
    auto pixels = load_pfm(
        filename.c_str(), &image.width, &image.height, &ncomp, 4);
    if (!pixels) throw io_error::read_error(filename);
    image.linear = true;
    image.pixels = from_linear(pixels, image.width, image.height);
    delete[] pixels;
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
  auto pfm_write_data = [](void* context, const void* data, int size) {
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
  } else if (ext == ".pfm" || ext == ".PFM") {
    auto buffer = vector<byte>{};
    if (!save_pfm_to_func(pfm_write_data, &buffer, image.width, image.height, 4,
            (const float*)to_linear(image).data()))
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

// -----------------------------------------------------------------------------
// TEXTURE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Loads/saves an image. Chooses hdr or ldr based on file name.
texture_data load_texture(const string& filename) {
  auto texture = texture_data{};
  load_texture(filename, texture);
  return texture;
}
void load_texture(const string& filename, texture_data& texture) {
  auto ext = path_extension(filename);
  if (ext == ".exr" || ext == ".EXR") {
    auto pixels = (float*)nullptr;
    if (LoadEXR(&pixels, &texture.width, &texture.height, filename.c_str(),
            nullptr) != 0)
      throw io_error::read_error(filename);
    texture.linear  = true;
    texture.pixelsf = vector<vec4f>{
        (vec4f*)pixels, (vec4f*)pixels + texture.width * texture.height};
    free(pixels);
  } else if (ext == ".pfm" || ext == ".PFM") {
    auto ncomp  = 0;
    auto pixels = load_pfm(
        filename.c_str(), &texture.width, &texture.height, &ncomp, 4);
    if (!pixels) throw io_error::read_error(filename);
    texture.linear  = true;
    texture.pixelsf = vector<vec4f>{
        (vec4f*)pixels, (vec4f*)pixels + texture.width * texture.height};
    delete[] pixels;
  } else if (ext == ".hdr" || ext == ".HDR") {
    auto buffer = load_binary(filename);
    auto ncomp  = 0;
    auto pixels = stbi_loadf_from_memory(buffer.data(), (int)buffer.size(),
        &texture.width, &texture.height, &ncomp, 4);
    if (!pixels) throw io_error::read_error(filename);
    texture.linear  = true;
    texture.pixelsf = vector<vec4f>{
        (vec4f*)pixels, (vec4f*)pixels + texture.width * texture.height};
    free(pixels);
  } else if (ext == ".png" || ext == ".PNG") {
    auto buffer = load_binary(filename);
    auto ncomp  = 0;
    auto pixels = stbi_load_from_memory(buffer.data(), (int)buffer.size(),
        &texture.width, &texture.height, &ncomp, 4);
    if (!pixels) throw io_error::read_error(filename);
    texture.linear  = false;
    texture.pixelsb = vector<vec4b>{
        (vec4b*)pixels, (vec4b*)pixels + texture.width * texture.height};
    free(pixels);
  } else if (ext == ".jpg" || ext == ".JPG" || ext == ".jpeg" ||
             ext == ".JPEG") {
    auto buffer = load_binary(filename);
    auto ncomp  = 0;
    auto pixels = stbi_load_from_memory(buffer.data(), (int)buffer.size(),
        &texture.width, &texture.height, &ncomp, 4);
    if (!pixels) throw io_error::read_error(filename);
    texture.linear  = false;
    texture.pixelsb = vector<vec4b>{
        (vec4b*)pixels, (vec4b*)pixels + texture.width * texture.height};
    free(pixels);
  } else if (ext == ".tga" || ext == ".TGA") {
    auto buffer = load_binary(filename);
    auto ncomp  = 0;
    auto pixels = stbi_load_from_memory(buffer.data(), (int)buffer.size(),
        &texture.width, &texture.height, &ncomp, 4);
    if (!pixels) throw io_error::read_error(filename);
    texture.linear  = false;
    texture.pixelsb = vector<vec4b>{
        (vec4b*)pixels, (vec4b*)pixels + texture.width * texture.height};
    free(pixels);
  } else if (ext == ".bmp" || ext == ".BMP") {
    auto buffer = load_binary(filename);
    auto ncomp  = 0;
    auto pixels = stbi_load_from_memory(buffer.data(), (int)buffer.size(),
        &texture.width, &texture.height, &ncomp, 4);
    if (!pixels) throw io_error::read_error(filename);
    texture.linear  = false;
    texture.pixelsb = vector<vec4b>{
        (vec4b*)pixels, (vec4b*)pixels + texture.width * texture.height};
    free(pixels);
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    // create preset
    texture = make_texture_preset(path_basename(filename));
  } else {
    throw io_error::format_error(filename);
  }
}

// Saves an hdr image.
void save_texture(const string& filename, const texture_data& texture) {
  // check for correct handling
  if (!texture.pixelsf.empty() && is_ldr_filename(filename))
    throw io_error{filename, "cannot save hdr texture to ldr file"};
  if (!texture.pixelsb.empty() && is_hdr_filename(filename))
    throw io_error{filename, "cannot save ldr texture to hdr file"};

  // write data
  auto stbi_write_data = [](void* context, void* data, int size) {
    auto& buffer = *(vector<byte>*)context;
    buffer.insert(buffer.end(), (byte*)data, (byte*)data + size);
  };
  auto pfm_write_data = [](void* context, const void* data, int size) {
    auto& buffer = *(vector<byte>*)context;
    buffer.insert(buffer.end(), (byte*)data, (byte*)data + size);
  };

  auto ext = path_extension(filename);
  if (ext == ".hdr" || ext == ".HDR") {
    auto buffer = vector<byte>{};
    if (!stbi_write_hdr_to_func(stbi_write_data, &buffer, (int)texture.width,
            (int)texture.height, 4, (const float*)texture.pixelsf.data()))
      throw io_error::write_error(filename);
    save_binary(filename, buffer);
  } else if (ext == ".pfm" || ext == ".PFM") {
    auto buffer = vector<byte>{};
    if (!save_pfm_to_func(pfm_write_data, &buffer, texture.width,
            texture.height, 4, (const float*)texture.pixelsf.data()))
      throw io_error::write_error(filename);
    save_binary(filename, buffer);
  } else if (ext == ".exr" || ext == ".EXR") {
    auto data = (byte*)nullptr;
    auto size = (size_t)0;
    if (SaveEXRToMemory((const float*)texture.pixelsf.data(),
            (int)texture.width, (int)texture.height, 4, 1, &data, &size,
            nullptr) < 0)
      throw io_error::write_error(filename);
    auto buffer = vector<byte>{data, data + size};
    free(data);
    save_binary(filename, buffer);
  } else if (ext == ".png" || ext == ".PNG") {
    auto buffer = vector<byte>{};
    if (!stbi_write_png_to_func(stbi_write_data, &buffer, (int)texture.width,
            (int)texture.height, 4, (const byte*)texture.pixelsb.data(),
            (int)texture.width * 4))
      throw io_error::write_error(filename);
    save_binary(filename, buffer);
  } else if (ext == ".jpg" || ext == ".JPG" || ext == ".jpeg" ||
             ext == ".JPEG") {
    auto buffer = vector<byte>{};
    if (!stbi_write_jpg_to_func(stbi_write_data, &buffer, (int)texture.width,
            (int)texture.height, 4, (const byte*)texture.pixelsb.data(), 75))
      throw io_error::write_error(filename);
    save_binary(filename, buffer);
  } else if (ext == ".tga" || ext == ".TGA") {
    auto buffer = vector<byte>{};
    if (!stbi_write_tga_to_func(stbi_write_data, &buffer, (int)texture.width,
            (int)texture.height, 4, (const byte*)texture.pixelsb.data()))
      throw io_error::write_error(filename);
    save_binary(filename, buffer);
  } else if (ext == ".bmp" || ext == ".BMP") {
    auto buffer = vector<byte>{};
    if (!stbi_write_bmp_to_func(stbi_write_data, &buffer, (int)texture.width,
            (int)texture.height, 4, (const byte*)texture.pixelsb.data()))
      throw io_error::write_error(filename);
    save_binary(filename, buffer);
  } else {
    throw io_error::format_error(filename);
  }
}

texture_data make_texture_preset(const string& type) {
  return image_to_texture(make_image_preset(type));
}

// Loads/saves an image. Chooses hdr or ldr based on file name.
bool load_texture(
    const string& filename, texture_data& texture, string& error) {
  try {
    load_texture(filename, texture);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// Saves an hdr image.
bool save_texture(
    const string& filename, const texture_data& texture, string& error) {
  try {
    save_texture(filename, texture);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

bool make_texture_preset(
    texture_data& texture, const string& type, string& error) {
  try {
    texture = make_texture_preset(type);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load ply mesh
shape_data load_shape(const string& filename, bool flip_texcoord) {
  auto shape = shape_data{};
  load_shape(filename, shape, flip_texcoord);
  return shape;
}
void load_shape(const string& filename, shape_data& shape, bool flip_texcoord) {
  shape = {};

  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = ply_model{};
    load_ply(filename, ply);
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
      throw io_error::shape_error(filename);
  } else if (ext == ".obj" || ext == ".OBJ") {
    auto obj = obj_shape{};
    load_obj(filename, obj, false);
    auto materials = vector<int>{};
    get_positions(obj, shape.positions);
    get_normals(obj, shape.normals);
    get_texcoords(obj, shape.texcoords, flip_texcoord);
    get_faces(obj, shape.triangles, shape.quads, materials);
    get_lines(obj, shape.lines, materials);
    get_points(obj, shape.points, materials);
    if (shape.points.empty() && shape.lines.empty() &&
        shape.triangles.empty() && shape.quads.empty())
      throw io_error::shape_error(filename);
  } else if (ext == ".stl" || ext == ".STL") {
    auto stl = stl_model{};
    load_stl(filename, stl, true);
    if (stl.shapes.size() != 1) throw io_error::shape_error(filename);
    auto fnormals = vector<vec3f>{};
    if (!get_triangles(stl, 0, shape.triangles, shape.positions, fnormals))
      throw io_error::shape_error(filename);
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    shape = make_shape_preset(path_basename(filename));
  } else {
    throw io_error::format_error(filename);
  }
}

// Save ply mesh
void save_shape(const string& filename, const shape_data& shape,
    bool flip_texcoord, bool ascii) {
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
    save_ply(filename, ply);
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
    save_obj(filename, obj);
  } else if (ext == ".stl" || ext == ".STL") {
    auto stl = stl_model{};
    if (!shape.lines.empty()) throw io_error{filename, "lines not supported"};
    if (!shape.points.empty()) throw io_error{filename, "points not supported"};
    if (!shape.triangles.empty()) {
      add_triangles(stl, shape.triangles, shape.positions, {});
    } else if (!shape.quads.empty()) {
      add_triangles(stl, quads_to_triangles(shape.quads), shape.positions, {});
    } else {
      throw io_error::shape_error(filename);
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
    throw io_error::format_error(filename);
  }
}

// Load face-varying mesh
fvshape_data load_fvshape(const string& filename, bool flip_texcoord) {
  auto shape = fvshape_data{};
  load_fvshape(filename, shape, flip_texcoord);
  return shape;
}
void load_fvshape(
    const string& filename, fvshape_data& shape, bool flip_texcoord) {
  shape = {};

  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = ply_model{};
    load_ply(filename, ply);
    get_positions(ply, shape.positions);
    get_normals(ply, shape.normals);
    get_texcoords(ply, shape.texcoords, flip_texcoord);
    get_quads(ply, shape.quadspos);
    if (!shape.normals.empty()) shape.quadsnorm = shape.quadspos;
    if (!shape.texcoords.empty()) shape.quadstexcoord = shape.quadspos;
    if (shape.quadspos.empty()) throw io_error::shape_error(filename);
  } else if (ext == ".obj" || ext == ".OBJ") {
    auto obj = obj_shape{};
    load_obj(filename, obj, true);
    auto materials = vector<int>{};
    get_positions(obj, shape.positions);
    get_normals(obj, shape.normals);
    get_texcoords(obj, shape.texcoords, flip_texcoord);
    get_fvquads(
        obj, shape.quadspos, shape.quadsnorm, shape.quadstexcoord, materials);
    if (shape.quadspos.empty()) throw io_error::shape_error(filename);
  } else if (ext == ".stl" || ext == ".STL") {
    auto stl = stl_model{};
    load_stl(filename, stl, true);
    if (stl.shapes.empty()) throw io_error::shape_error(filename);
    if (stl.shapes.size() > 1) throw io_error::shape_error(filename);
    auto fnormals  = vector<vec3f>{};
    auto triangles = vector<vec3i>{};
    if (!get_triangles(stl, 0, triangles, shape.positions, fnormals))
      throw io_error::shape_error(filename);
    shape.quadspos = triangles_to_quads(triangles);
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    shape = make_fvshape_preset(path_basename(filename));
  } else {
    throw io_error::format_error(filename);
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
    add_positions(ply, split_positions);
    add_normals(ply, split_normals);
    add_texcoords(ply, split_texcoords, flip_texcoord);
    add_faces(ply, {}, split_quads);
    save_ply(filename, ply);
  } else if (ext == ".obj" || ext == ".OBJ") {
    auto obj = obj_shape{};
    add_positions(obj, shape.positions);
    add_normals(obj, shape.positions);
    add_texcoords(obj, shape.texcoords, flip_texcoord);
    add_fvquads(obj, shape.quadspos, shape.quadsnorm, shape.quadstexcoord, 0);
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
      add_triangles(stl, quads_to_triangles(split_quads), split_positions, {});
    } else {
      throw io_error::shape_error(filename);
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
    throw io_error::format_error(filename);
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
  } else if (type == "test-particles") {
    return make_points(4096);
  } else if (type == "test-cloth") {
    return make_rect({64, 64}, {0.2f, 0.2f});
  } else if (type == "test-clothy") {
    return make_recty({64, 64}, {0.2f, 0.2f});
  } else {
    throw io_error::preset_error(type);
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
    throw io_error::preset_error(type);
  }
}

// Load ply mesh
bool load_shape(const string& filename, shape_data& shape, string& error,
    bool flip_texcoord) {
  try {
    load_shape(filename, shape, flip_texcoord);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// Save ply mesh
bool save_shape(const string& filename, const shape_data& shape, string& error,
    bool flip_texcoord, bool ascii) {
  try {
    save_shape(filename, shape, flip_texcoord, ascii);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// Load ply mesh
bool load_fvshape(const string& filename, fvshape_data& fvshape, string& error,
    bool flip_texcoord) {
  try {
    load_fvshape(filename, fvshape, flip_texcoord);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// Save ply mesh
bool save_fvshape(const string& filename, const fvshape_data& fvshape,
    string& error, bool flip_texcoord, bool ascii) {
  try {
    save_fvshape(filename, fvshape, flip_texcoord, ascii);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// Shape presets used ofr testing.
bool make_shape_preset(shape_data& shape, const string& type, string& error) {
  try {
    shape = make_shape_preset(type);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// Shape presets used for testing.
bool make_fvshape_preset(
    fvshape_data& fvshape, const string& type, string& error) {
  try {
    fvshape = make_fvshape_preset(type);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
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
    throw io_error::preset_error(type);
  }
}

// Scene presets used for testing.
bool make_scene_preset(scene_data& scene, const string& type, string& error) {
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
static void load_json_scene(
    const string& filename, scene_data& scene, bool noparallel);
static void save_json_scene(
    const string& filename, const scene_data& scene, bool noparallel);

// Load/save a scene from/to OBJ.
static void load_obj_scene(
    const string& filename, scene_data& scene, bool noparallel);
static void save_obj_scene(
    const string& filename, const scene_data& scene, bool noparallel);

// Load/save a scene from/to PLY. Loads/saves only one mesh with no other data.
static void load_ply_scene(
    const string& filename, scene_data& scene, bool noparallel);
static void save_ply_scene(
    const string& filename, const scene_data& scene, bool noparallel);

// Load/save a scene from/to STL. Loads/saves only one mesh with no other data.
static void load_stl_scene(
    const string& filename, scene_data& scene, bool noparallel);
static void save_stl_scene(
    const string& filename, const scene_data& scene, bool noparallel);

// Load/save a scene from/to glTF.
static void load_gltf_scene(
    const string& filename, scene_data& scene, bool noparallel);
static void save_gltf_scene(
    const string& filename, const scene_data& scene, bool noparallel);

// Load/save a scene from/to pbrt-> This is not robust at all and only
// works on scene that have been previously adapted since the two renderers
// are too different to match.
static void load_pbrt_scene(
    const string& filename, scene_data& scene, bool noparallel);
static void save_pbrt_scene(
    const string& filename, const scene_data& scene, bool noparallel);

// Load a scene
scene_data load_scene(const string& filename, bool noparallel) {
  auto scene = scene_data{};
  load_scene(filename, scene, noparallel);
  return scene;
}

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
  } else if (ext == ".ply" || ext == ".PLY") {
    return load_ply_scene(filename, scene, noparallel);
  } else if (ext == ".stl" || ext == ".STL") {
    return load_stl_scene(filename, scene, noparallel);
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    scene = make_scene_preset(path_basename(filename));
  } else {
    throw io_error::format_error(filename);
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
  } else if (ext == ".pbrt" || ext == ".PBRT") {
    return save_pbrt_scene(filename, scene, noparallel);
  } else if (ext == ".gltf" || ext == ".GLTF") {
    return save_gltf_scene(filename, scene, noparallel);
  } else if (ext == ".ply" || ext == ".PLY") {
    return save_ply_scene(filename, scene, noparallel);
  } else if (ext == ".stl" || ext == ".STL") {
    return save_stl_scene(filename, scene, noparallel);
  } else {
    throw io_error::format_error(filename);
  }
}

// Load a scene
bool load_scene(
    const string& filename, scene_data& scene, string& error, bool noparallel) {
  try {
    load_scene(filename, scene, noparallel);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// Save a scene
bool save_scene(const string& filename, const scene_data& scene, string& error,
    bool noparallel) {
  try {
    save_scene(filename, scene, noparallel);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
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
void add_environment(scene_data& scene, const string& filename) {
  scene.textures.push_back(load_texture(filename));
  scene.environments.push_back(
      {identity3x4f, {1, 1, 1}, (int)scene.textures.size() - 1});
}

// Make missing scene directories
bool make_scene_directories(
    const string& filename, const scene_data& scene, string& error) {
  try {
    make_scene_directories(filename, scene);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// Add environment
bool add_environment(scene_data& scene, const string& filename, string& error) {
  try {
    add_environment(scene, filename);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
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
    get_values(ply, "instance",
        {"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz", "ox", "oy",
            "oz"},
        frames);
  } else {
    io_error::format_error(filename);
  }
}

// save instances
[[maybe_unused]] static void save_instance(
    const string& filename, const vector<frame3f>& frames, bool ascii = false) {
  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = ply_model{};
    add_values(ply, "instance",
        {"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz", "ox", "oy",
            "oz"},
        frames);
    save_ply(filename, ply);
  } else {
    throw io_error::format_error(filename);
  }
}

// load subdiv
subdiv_data load_subdiv(const string& filename) {
  auto subdiv = subdiv_data{};
  load_subdiv(filename, subdiv);
  return subdiv;
}
void load_subdiv(const string& filename, subdiv_data& subdiv) {
  auto lsubdiv = fvshape_data{};
  load_fvshape(filename, lsubdiv, true);
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

// load subdiv
bool load_subdiv(const string& filename, subdiv_data& subdiv, string& error) {
  try {
    load_subdiv(filename, subdiv);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// save subdiv
bool save_subdiv(
    const string& filename, const subdiv_data& subdiv, string& error) {
  try {
    save_subdiv(filename, subdiv);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
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

// Json specializations
template <>
struct json_enum_trait<material_type> {
  static const vector<string>& labels() { return material_type_names; }
};

// Load a scene in the builtin JSON format.
static void load_json_scene_version40(const string& filename,
    const json_value& json, scene_data& scene, bool noparallel) {
  auto parse_error = [filename](const string& patha, const string& pathb = "",
                         const string& pathc = "") {
    auto path = patha;
    if (!pathb.empty()) path += "/" + pathb;
    if (!pathc.empty()) path += "/" + pathc;
    throw io_error(filename, "parse error at " + path);
  };
  auto key_error = [filename](const string& patha, const string& pathb = "",
                       const string& pathc = "") {
    auto path = patha;
    if (!pathb.empty()) path += "/" + pathb;
    if (!pathc.empty()) path += "/" + pathc;
    throw io_error(filename, "unknow key at " + path);
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
  } catch (const json_error& error) {
    throw io_error::parse_error(
        filename, json_value::get_path(&json, error.where()));
  } catch (...) {
    throw io_error::parse_error(filename);
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
        load_shape(path_join(dirname, path), shape, true);
      });
      // load subdivs
      parallel_foreach(scene.subdivs, [&](auto& subdiv) {
        auto path = find_path(
            get_subdiv_name(scene, subdiv), "subdivs", {".ply", ".obj"});
        load_subdiv(path_join(dirname, path), subdiv);
      });
      // load textures
      parallel_foreach(scene.textures, [&](auto& texture) {
        auto path = find_path(get_texture_name(scene, texture), "textures",
            {".hdr", ".exr", ".png", ".jpg"});
        load_texture(path_join(dirname, path), texture);
      });
      // load instances
      parallel_foreach(ply_instances, [&](auto& ply_instance) {
        auto path = find_path(
            get_ply_instance_name(scene, ply_instance), "instances", {".ply"});
        load_instance(path_join(dirname, path), ply_instance.frames);
      });
    }
  } catch (const io_error& exception) {
    throw io_error::dependent_error(filename, exception);
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
static void load_json_scene(
    const string& filename, scene_data& scene, bool noparallel) {
  // ends with
  auto ends_with = [](const string& str, const string& end) {
    if (str.size() < end.size()) return false;
    return str.substr(str.size() - end.size()) == end;
  };

  // open file
  auto json = load_json(filename);

  // check vrsion
  if (!json.contains("asset") || !json.at("asset").contains("version"))
    return load_json_scene_version40(filename, json, scene, noparallel);

  // parse json value
  auto get_opt = [](const json_value& json, const string& key, auto& value) {
    value = json.value(key, value);
  };
  auto get_req = [](const json_value& json, const string& key, auto& value) {
    json.at(key).get(value);
  };
  auto get_orf = [](const json_value& json, const string& key, int& value,
                     const unordered_map<string, int>& map) {
    auto values = json.value(key, string{});
    value       = values.empty() ? -1 : map.at(values);
  };
  auto get_rrf = [](const json_value& json, const string& key, int& value,
                     const unordered_map<string, int>& map) {
    auto values = json.at(key).get<string>();
    value       = map.at(values);
  };

  // references
  auto shape_map    = unordered_map<string, int>{};
  auto texture_map  = unordered_map<string, int>{};
  auto material_map = unordered_map<string, int>{};

  // filenames
  auto shape_datafiles   = vector<string>{};
  auto texture_datafiles = vector<string>{};
  auto subdiv_datafiles  = vector<string>{};

  // // load json instance
  // struct ply_instance {
  //   vector<frame3f> frames = {};
  // };
  // using ply_instance_handle = int;
  // auto ply_instances        = vector<ply_instance>{};
  // auto ply_instances_names  = vector<string>{};
  // auto ply_instance_map     = unordered_map<string, ply_instance_handle>{
  //     {"", invalidid}};
  // auto instance_ply = unordered_map<int, ply_instance_handle>{};
  // auto get_ist      = [&scene, &ply_instances, &ply_instances_names,
  //                    &ply_instance_map, &instance_ply](const json_value&
  //                    json, const string& key, const instance_data& instance)
  //                    {
  //   auto name = json.value(key, string{});
  //   if (name.empty()) return;
  //   auto instance_id = (int)(&instance - scene.instances.data());
  //   auto it          = ply_instance_map.find(name);
  //   if (it != ply_instance_map.end()) {
  //     instance_ply[instance_id] = it->second;
  //   } else {
  //     ply_instances_names.emplace_back(name);
  //     ply_instances.emplace_back(ply_instance());
  //     auto ply_instance_id      = (int)ply_instances.size() - 1;
  //     ply_instance_map[name]    = ply_instance_id;
  //     instance_ply[instance_id] = ply_instance_id;
  //   }
  // };
  // auto get_ply_instance_name = [&ply_instances, &ply_instances_names](
  //                                  const scene_data&   scene,
  //                                  const ply_instance& instance) -> string {
  //   return ply_instances_names[&instance - ply_instances.data()];
  // };

  // prepare data
  auto dirname = path_dirname(filename);

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
        if (element.is_string()) {
          auto filepath = element.get<string>();
          element       = load_json(path_join(dirname, "cameras", filepath));
        }
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
      texture_datafiles.reserve(group.size());
      for (auto& [key, element] : group.items()) {
        [[maybe_unused]] auto& texture = scene.textures.emplace_back();
        scene.texture_names.push_back(key);
        auto& datafile   = texture_datafiles.emplace_back();
        texture_map[key] = (int)scene.textures.size() - 1;
        if (element.is_string()) {
          auto filename = element.get<string>();
          if (!ends_with(filename, ".json")) {
            element             = json_object();
            element["datafile"] = filename;
          } else {
            element = load_json(path_join(dirname, "textures", filename));
          }
        }
        get_req(element, "datafile", datafile);
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
        if (element.is_string()) {
          auto filename = element.get<string>();
          element       = load_json(path_join(dirname, "materials", filename));
        }
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
        get_orf(element, "emission_tex", material.emission_tex, texture_map);
        get_orf(element, "color_tex", material.color_tex, texture_map);
        get_orf(element, "roughness_tex", material.roughness_tex, texture_map);
        get_orf(
            element, "scattering_tex", material.scattering_tex, texture_map);
        get_orf(element, "normal_tex", material.normal_tex, texture_map);
      }
    }
    if (json.contains("shapes")) {
      auto& group = json.at("shapes");
      scene.shapes.reserve(group.size());
      scene.shape_names.reserve(group.size());
      shape_datafiles.reserve(group.size());
      for (auto& [key, element] : group.items()) {
        [[maybe_unused]] auto& shape = scene.shapes.emplace_back();
        scene.shape_names.push_back(key);
        auto& datafile = shape_datafiles.emplace_back();
        shape_map[key] = (int)scene.shapes.size() - 1;
        if (element.is_string()) {
          auto filename = element.get<string>();
          if (!ends_with(filename, ".json")) {
            element             = json_object();
            element["datafile"] = filename;
          } else {
            element = load_json(path_join(dirname, "shapes", filename));
          }
        }
        get_req(element, "datafile", datafile);
      }
    }
    if (json.contains("subdivs")) {
      auto& group = json.at("subdivs");
      scene.subdivs.reserve(group.size());
      scene.subdiv_names.reserve(group.size());
      subdiv_datafiles.reserve(group.size());
      for (auto& [key, element] : group.items()) {
        auto& subdiv = scene.subdivs.emplace_back();
        scene.subdiv_names.emplace_back(key);
        auto& datafile = subdiv_datafiles.emplace_back();
        if (element.is_string()) {
          auto filename = element.get<string>();
          element       = load_json(path_join(dirname, "subdivs", filename));
        }
        get_req(element, "datafile", datafile);
        get_rrf(element, "shape", subdiv.shape, shape_map);
        get_opt(element, "subdivisions", subdiv.subdivisions);
        get_opt(element, "catmullclark", subdiv.catmullclark);
        get_opt(element, "smooth", subdiv.smooth);
        get_opt(element, "displacement", subdiv.displacement);
        get_orf(
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
        if (element.is_string()) {
          auto filename = element.get<string>();
          element       = load_json(path_join(dirname, "instances", filename));
        }
        get_opt(element, "frame", instance.frame);
        get_rrf(element, "shape", instance.shape, shape_map);
        get_rrf(element, "material", instance.material, material_map);
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
        if (element.is_string()) {
          auto filename = element.get<string>();
          element = load_json(path_join(dirname, "environments", filename));
        }
        get_opt(element, "frame", environment.frame);
        get_opt(element, "emission", environment.emission);
        get_orf(element, "emission_tex", environment.emission_tex, texture_map);
        if (element.contains("lookat")) {
          get_opt(element, "lookat", (mat3f&)environment.frame);
          environment.frame = lookat_frame(environment.frame.x,
              environment.frame.y, environment.frame.z, false);
        }
      }
    }
    // if (json.contains("object_instances")) {
    //   for (auto& [key, element] : json.at("object_instances").items()) {
    //     auto& instance = ply_instances.emplace_back();
    //     if (element.is_string()) {
    //       auto filename = element.get<string>();
    //       element       = load_json(path_join(dirname, "objects", filename));
    //     }
    //     get_opt(element, "frame", instance.frame);
    //     get_rrf(element, "shape", instance.shape, shape_map);
    //     get_rrf(element, "material", instance.material, material_map);
    //     if (element.contains("lookat")) {
    //       get_opt(element, "lookat", (mat3f&)instance.frame);
    //       instance.frame = lookat_frame(
    //           instance.frame.x, instance.frame.y, instance.frame.z, false);
    //     }
    //     if (element.contains("ply_instance")) {
    //       get_ist(element, "instance", instance);
    //     }
    //   }
    // }
    // if (json.contains("objects")) {
    //   for (auto& [key, element] : json.at("objects").items()) {
    //     auto& instance = scene.instances.emplace_back();
    //     scene.instance_names.emplace_back(key);
    //     if (element.is_string()) {
    //       auto filename = element.get<string>();
    //       element       = load_json(path_join(dirname, "objects", filename));
    //     }
    //     get_opt(element, "frame", instance.frame);
    //     get_rrf(element, "shape", instance.shape, shape_map);
    //     get_rrf(element, "material", instance.material, material_map);
    //     if (element.contains("lookat")) {
    //       get_opt(element, "lookat", (mat3f&)instance.frame);
    //       instance.frame = lookat_frame(
    //           instance.frame.x, instance.frame.y, instance.frame.z, false);
    //     }
    //     if (element.contains("ply_instance")) {
    //       get_ist(element, "instance", instance);
    //     }
    //   }
    // }
  } catch (const json_error& error) {
    throw io_error::parse_error(
        filename, json_value::get_path(&json, error.where()));
  } catch (const io_error& error) {
    throw io_error::dependent_error(filename, error);
  } catch (...) {
    throw io_error::parse_error(filename);
  }

  // fix paths
  for (auto& datafile : shape_datafiles)
    datafile = path_join(dirname, "shapes", datafile);
  for (auto& datafile : texture_datafiles)
    datafile = path_join(dirname, "textures", datafile);
  for (auto& datafile : subdiv_datafiles)
    datafile = path_join(dirname, "subdivs", datafile);

  // load resources
  try {
    if (noparallel) {
      // load shapes
      for (auto idx : range(scene.shapes.size())) {
        load_shape(shape_datafiles[idx], scene.shapes[idx], true);
      }
      // load subdivs
      for (auto idx : range(scene.subdivs.size())) {
        load_subdiv(subdiv_datafiles[idx], scene.subdivs[idx]);
      }
      // load textures
      for (auto idx : range(scene.textures.size())) {
        load_texture(texture_datafiles[idx], scene.textures[idx]);
      }
      // load instances
      // for (auto& ply_instance : ply_instances) {
      //   auto path = find_path(
      //       get_ply_instance_name(scene, ply_instance), "instances",
      //       {".ply"});
      //   load_instance(path_join(dirname, path), ply_instance.frames);
      // }
    } else {
      // load shapes
      parallel_for(scene.shapes.size(), [&](size_t idx) {
        load_shape(shape_datafiles[idx], scene.shapes[idx], true);
      });
      // load subdivs
      parallel_for(scene.subdivs.size(), [&](size_t idx) {
        load_subdiv(subdiv_datafiles[idx], scene.subdivs[idx]);
      });
      // load textures
      parallel_for(scene.textures.size(), [&](size_t idx) {
        load_texture(texture_datafiles[idx], scene.textures[idx]);
      });
      // // load instances
      // parallel_foreach(ply_instances, [&](auto& ply_instance) {
      //   auto path = find_path(
      //       get_ply_instance_name(scene, ply_instance), "instances",
      //       {".ply"});
      //   load_instance(path_join(dirname, path), ply_instance.frames);
      // });
    }
  } catch (const io_error& exception) {
    throw io_error::dependent_error(filename, exception);
  }

  // apply instances
  // if (!ply_instances.empty()) {
  //   auto instances      = scene.instances;
  //   auto instance_names = scene.instance_names;
  //   scene.instances.clear();
  //   scene.instance_names.clear();
  //   for (auto& instance : instances) {
  //     auto it = instance_ply.find((int)(&instance - instances.data()));
  //     if (it == instance_ply.end()) {
  //       auto& ninstance = scene.instances.emplace_back();
  //       scene.instance_names.emplace_back(
  //           instance_names[&instance - instances.data()]);
  //       ninstance.frame    = instance.frame;
  //       ninstance.shape    = instance.shape;
  //       ninstance.material = instance.material;
  //     } else {
  //       auto& ply_instance = ply_instances[it->second];
  //       auto  instance_id  = 0;
  //       for (auto& frame : ply_instance.frames) {
  //         auto& ninstance = scene.instances.emplace_back();
  //         scene.instance_names.emplace_back(
  //             instance_names[&instance - instances.data()] + "_" +
  //             std::to_string(instance_id++));
  //         ninstance.frame    = frame * instance.frame;
  //         ninstance.shape    = instance.shape;
  //         ninstance.material = instance.material;
  //       }
  //     }
  //   }
  // }

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
    auto& item = json.insert_back(name);
    item       = json_object{};
    return item;
  };
  auto add_value = [](json_value& json, const string& name,
                       const auto& value) -> json_value& {
    auto& item = json.insert_back(name);
    item       = value;
    return item;
  };
  auto set_opt = [](json_value& json, const string& name, const auto& value,
                     const auto& def) {
    if (value == def) return;
    json.insert_back(name) = value;
  };
  auto set_req = [](json_value& json, const string& name, const auto& value) {
    json.insert_back(name) = value;
  };
  auto set_orf = [](json_value& json, const string& name, int value,
                     const vector<string>& names) {
    if (value < 0) return;
    json.insert_back(name) = names.at(value);
  };
  auto set_rrf = [](json_value& json, const string& name, int value,
                     const vector<string>& names) {
    json.insert_back(name) = names.at(value);
  };

  // dirname
  auto dirname = path_dirname(filename);

  // names
  auto camera_names  = make_names(scene.cameras, scene.camera_names, "camera");
  auto texture_names = make_names(
      scene.textures, scene.texture_names, "texture");
  auto material_names = make_names(
      scene.materials, scene.material_names, "material");
  auto shape_names    = make_names(scene.shapes, scene.shape_names, "shape");
  auto subdiv_names   = make_names(scene.subdivs, scene.subdiv_names, "subdiv");
  auto instance_names = make_names(
      scene.instances, scene.instance_names, "instance");

  // filenames
  auto shape_datafiles   = vector<string>(shape_names.size());
  auto texture_datafiles = vector<string>(texture_names.size());
  auto subdiv_datafiles  = vector<string>(subdiv_names.size());
  for (auto idx : range(shape_datafiles.size())) {
    shape_datafiles[idx] = path_join(
        dirname, "shapes", shape_names[idx] + ".ply");
  }
  for (auto idx : range(texture_datafiles.size())) {
    texture_datafiles[idx] = path_join(dirname, "textures",
        texture_names[idx] +
            (scene.textures[idx].pixelsf.empty() ? ".png" : ".hdr"));
  }
  for (auto idx : range(subdiv_datafiles.size())) {
    subdiv_datafiles[idx] = path_join(
        dirname, "subdivs", subdiv_names[idx] + ".obj");
  }

  // save json file
  auto json = json_value{};
  json      = json_object{};

  // asset
  {
    auto& element = add_object(json, "asset");
    set_opt(element, "copyright", scene.copyright, "");
    set_req(element, "generator",
        "Yocto/GL - https://github.com/xelatihy/yocto-gl");
    set_req(element, "version", "4.1");
  }

  if (!scene.cameras.empty()) {
    auto  default_ = sceneio_camera{};
    auto& group    = add_object(json, "cameras");
    group.reserve(scene.cameras.size());
    for (auto&& [idx, camera] : enumerate(scene.cameras)) {
      auto& element = add_object(group, camera_names[idx]);
      set_opt(element, "frame", camera.frame, default_.frame);
      set_opt(
          element, "orthographic", camera.orthographic, default_.orthographic);
      set_opt(element, "lens", camera.lens, default_.lens);
      set_opt(element, "aspect", camera.aspect, default_.aspect);
      set_opt(element, "film", camera.film, default_.film);
      set_opt(element, "focus", camera.focus, default_.focus);
      set_opt(element, "aperture", camera.aperture, default_.aperture);
    }
  }

  if (!scene.textures.empty()) {
    auto& group = add_object(json, "textures");
    group.reserve(scene.textures.size());
    for (auto&& [idx, texture] : enumerate(scene.textures)) {
      add_value(group, texture_names[idx],
          texture_names[idx] + (texture.pixelsf.empty() ? ".png" : ".hdr"));
    }
  }

  if (!scene.materials.empty()) {
    auto  default_ = sceneio_material{};
    auto& group    = add_object(json, "materials");
    group.reserve(scene.materials.size());
    for (auto&& [idx, material] : enumerate(scene.materials)) {
      auto& element = add_object(group, material_names[idx]);
      set_opt(element, "type", material.type, default_.type);
      set_opt(element, "emission", material.emission, default_.emission);
      set_opt(element, "color", material.color, default_.color);
      set_opt(element, "metallic", material.metallic, default_.metallic);
      set_opt(element, "roughness", material.roughness, default_.roughness);
      set_opt(element, "ior", material.ior, default_.ior);
      set_opt(element, "trdepth", material.trdepth, default_.trdepth);
      set_opt(element, "scattering", material.scattering, default_.scattering);
      set_opt(element, "scanisotropy", material.scanisotropy,
          default_.scanisotropy);
      set_opt(element, "opacity", material.opacity, default_.opacity);
      set_orf(element, "emission_tex", material.emission_tex, texture_names);
      set_orf(element, "color_tex", material.color_tex, texture_names);
      set_orf(element, "roughness_tex", material.roughness_tex, texture_names);
      set_orf(
          element, "scattering_tex", material.scattering_tex, texture_names);
      set_orf(element, "normal_tex", material.normal_tex, texture_names);
    }
  }

  if (!scene.shapes.empty()) {
    auto& group = add_object(json, "shapes");
    group.reserve(scene.shapes.size());
    for (auto&& [idx, shape] : enumerate(scene.shapes)) {
      add_value(group, shape_names[idx], shape_names[idx] + ".ply");
    }
  }

  if (!scene.subdivs.empty()) {
    auto  default_ = subdiv_data{};
    auto& group    = add_object(json, "subdivs");
    group.reserve(scene.subdivs.size());
    for (auto&& [idx, subdiv] : enumerate(scene.subdivs)) {
      auto& element = add_object(group, subdiv_names[idx]);
      set_rrf(element, "shape", subdiv.shape, shape_names);
      set_req(element, "datafile", subdiv_names[idx] + ".obj");
      set_opt(
          element, "subdivisions", subdiv.subdivisions, default_.subdivisions);
      set_opt(
          element, "catmullclark", subdiv.catmullclark, default_.subdivisions);
      set_opt(element, "smooth", subdiv.smooth, default_.subdivisions);
      set_opt(
          element, "displacement", subdiv.displacement, default_.subdivisions);
      set_opt(element, "displacement_tex",
          get_texture_name(scene, subdiv.displacement_tex), "");
    }
  }

  if (!scene.instances.empty()) {
    auto  default_ = sceneio_instance{};
    auto& group    = add_object(json, "instances");
    group.reserve(scene.instances.size());
    for (auto& instance : scene.instances) {
      auto& element = add_object(group, get_instance_name(scene, instance));
      set_opt(element, "frame", instance.frame, default_.frame);
      set_rrf(element, "shape", instance.shape, shape_names);
      set_rrf(element, "material", instance.material, material_names);
    }
  }

  if (!scene.environments.empty()) {
    auto  default_ = sceneio_environment{};
    auto& group    = add_object(json, "environments");
    group.reserve(scene.environments.size());
    for (auto& environment : scene.environments) {
      auto& element = add_object(
          group, get_environment_name(scene, environment));
      set_opt(element, "frame", environment.frame, default_.frame);
      set_opt(element, "emission", environment.emission, default_.emission);
      set_orf(element, "emission_tex", environment.emission_tex, texture_names);
    }
  }

  // save json
  save_json(filename, json);

  // dirname
  try {
    if (noparallel) {
      // save shapes
      for (auto idx : range(scene.shapes.size())) {
        save_shape(shape_datafiles[idx], scene.shapes[idx], true);
      }
      // save subdiv
      for (auto idx : range(scene.subdivs.size())) {
        save_subdiv(subdiv_datafiles[idx], scene.subdivs[idx]);
      }
      // save textures
      for (auto idx : range(scene.textures.size())) {
        save_texture(texture_datafiles[idx], scene.textures[idx]);
      }
    } else {
      // save shapes
      parallel_for(scene.shapes.size(), [&](auto idx) {
        save_shape(shape_datafiles[idx], scene.shapes[idx], true);
      });
      // save subdivs
      parallel_for(scene.subdivs.size(), [&](auto idx) {
        save_subdiv(subdiv_datafiles[idx], scene.subdivs[idx]);
      });
      // save textures
      parallel_for(scene.textures.size(), [&](auto idx) {
        save_texture(texture_datafiles[idx], scene.textures[idx]);
      });
    }
  } catch (const io_error& exception) {
    throw io_error::dependent_error(filename, exception);
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
  // dirname
  auto dirname = path_dirname(filename);

  // data files

  // load obj
  auto obj = load_obj(filename, false, true);

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
      material.type      = material_type::metallic;
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

  try {
    if (noparallel) {
      // load textures
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
  } catch (const io_error& exception) {
    throw io_error::dependent_error(filename, exception);
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
        save_texture(path_join(dirname, path), texture);
      });
    }
  } catch (const io_error& exception) {
    throw io_error::dependent_error(filename, exception);
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
  scene.shapes.push_back(load_shape(filename, true));
  scene.instances.push_back({identity3x4f, (int)scene.shapes.size() - 1, -1});

  // fix scene
  add_missing_material(scene);
  add_missing_camera(scene);
  add_missing_radius(scene);
}

static void save_ply_scene(
    const string& filename, const scene_data& scene, bool noparallel) {
  // save shape
  if (scene.shapes.empty()) throw io_error::shape_error(filename);
  save_shape(filename, scene.shapes.front(), false);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// STL CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static void load_stl_scene(
    const string& filename, scene_data& scene, bool noparallel) {
  // load ply mesh and make instance
  scene.shapes.push_back(load_shape(filename, true));
  scene.instances.push_back({identity3x4f, (int)scene.shapes.size() - 1, -1});

  // fix scene
  add_missing_material(scene);
  add_missing_camera(scene);
  add_missing_radius(scene);
}

static void save_stl_scene(
    const string& filename, const scene_data& scene, bool noparallel) {
  // save shape
  if (scene.shapes.empty()) throw io_error::shape_error(filename);
  save_shape(filename, scene.shapes.front(), false);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GLTF CONVESION
// -----------------------------------------------------------------------------
namespace yocto {

// Load a scene
static void load_gltf_scene(
    const string& filename, scene_data& scene, bool noparallel) {
  // load gltf
  auto gltf = load_json(filename);

  // parse buffers
  auto buffers_paths = vector<string>{};
  auto buffers       = vector<vector<byte>>();
  try {
    if (gltf.contains("buffers")) {
      for (auto& gbuffer : gltf.at("buffers")) {
        if (!gbuffer.contains("uri")) throw io_error::parse_error(filename);
        buffers_paths.push_back(gbuffer.value("uri", ""));
        buffers.emplace_back();
      }
    }
  } catch (...) {
    throw io_error::parse_error(filename);
  }

  // dirname
  auto dirname = path_dirname(filename);

  try {
    if (noparallel) {
      // load buffers
      for (auto& buffer : buffers) {
        auto& path = buffers_paths[&buffer - &buffers.front()];
        load_binary(path_join(dirname, path), buffer);
      }
    } else {
      // load buffers
      parallel_foreach(buffers, [&](auto& buffer) {
        auto& path = buffers_paths[&buffer - &buffers.front()];
        load_binary(path_join(dirname, path), buffer);
      });
    }
  } catch (const io_error& exception) {
    throw io_error::dependent_error(filename, exception);
  }

  // convert asset
  if (gltf.contains("asset")) {
    try {
      scene.copyright = gltf.value("copyright", ""s);
    } catch (...) {
      throw io_error::parse_error(filename);
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
          throw io_error::parse_error(filename);
        }
      }
    } catch (...) {
      throw io_error::parse_error(filename);
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
      throw io_error::parse_error(filename);
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
      throw io_error::parse_error(filename);
    }
  }

  // convert meshes
  auto mesh_primitives = vector<vector<sceneio_instance>>{};
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
            if (gaccessor.contains("sparse"))
              throw std::invalid_argument{"sparse accessor"};
            auto& gview =
                gltf.at("bufferViews").at(gaccessor.value("bufferView", -1));
            auto& buffer      = buffers.at(gview.value("buffer", 0));
            auto  components  = type_components.at(gaccessor.value("type", ""));
            auto  dcomponents = components;
            auto  count       = gaccessor.value("count", (size_t)0);
            auto  data        = (float*)nullptr;
            if (gname == "POSITION") {
              if (components != 3)
                throw std::invalid_argument{"invalid accessor"};
              shape.positions.resize(count);
              data = (float*)shape.positions.data();
            } else if (gname == "NORMAL") {
              if (components != 3)
                throw std::invalid_argument{"invalid accessor"};
              shape.normals.resize(count);
              data = (float*)shape.normals.data();
            } else if (gname == "TEXCOORD" || gname == "TEXCOORD_0") {
              if (components != 2)
                throw std::invalid_argument{"invalid accessor"};
              shape.texcoords.resize(count);
              data = (float*)shape.texcoords.data();
            } else if (gname == "COLOR" || gname == "COLOR_0") {
              if (components != 3 && components != 4)
                throw std::invalid_argument{"invalid accessor"};
              shape.colors.resize(count);
              data = (float*)shape.colors.data();
              if (components == 3) {
                dcomponents = 4;
                for (auto& c : shape.colors) c.w = 1;
              }
            } else if (gname == "TANGENT") {
              if (components != 4)
                throw std::invalid_argument{"invalid accessor"};
              shape.tangents.resize(count);
              data = (float*)shape.tangents.data();
            } else if (gname == "RADIUS") {
              if (components != 1)
                throw std::invalid_argument{"invalid accessor"};
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
              throw std::invalid_argument{"invalid accessor"};
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
              // points
              throw io_error{filename, "primitive_error"};
            } else {
              throw io_error{filename, "primitive_error"};
            }
          } else {
            auto& gaccessor =
                gltf.at("accessors").at(gprimitive.value("indices", -1));
            auto& gview =
                gltf.at("bufferViews").at(gaccessor.value("bufferView", -1));
            auto& buffer = buffers.at(gview.value("buffer", 0));
            if (gaccessor.value("type", "") != "SCALAR")
              throw std::invalid_argument{"invalid accessor"};
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
              throw std::invalid_argument{"invalid accessor"};
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
              throw io_error{filename, "primitive_error"};
            } else {
              throw io_error{filename, "primitive_error"};
            }
          }
        }
      }
    } catch (const io_error& error) {
      throw;
    } catch (...) {
      throw io_error::parse_error(filename);
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
      throw io_error::parse_error(filename);
    }
  }

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
  } catch (const io_error& exception) {
    throw io_error::dependent_error(filename, exception);
  }

  // fix scene
  add_missing_material(scene);
  add_missing_camera(scene);
  add_missing_radius(scene);
}

// Load a scene
static void save_gltf_scene(
    const string& filename, const scene_data& scene, bool noparallel) {
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
  save_json(filename, gltf);

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
        auto path = "textures/" + get_texture_name(scene, texture) +
                    (!texture.pixelsf.empty() ? ".hdr" : ".png");
        save_texture(path_join(dirname, path), texture);
      }
    } else {
      // save shapes
      parallel_foreach(scene.shapes, [&](auto& shape) {
        auto path = "shapes/" + get_shape_name(scene, shape) + ".bin";
        save_binshape(path_join(dirname, path), shape);
      });
      // save textures
      parallel_foreach(scene.textures, [&](auto& texture) {
        auto path = "textures/" + get_texture_name(scene, texture) +
                    (!texture.pixelsf.empty() ? ".hdr"s : ".png"s);
        save_texture(path_join(dirname, path), texture);
      });
    }
  } catch (const io_error& exception) {
    throw io_error::dependent_error(filename, exception);
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
      {pbrt_mtype::metal, material_type::metallic},
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
        load_texture(path_join(dirname, path), texture);
      });
    }
  } catch (const io_error& exception) {
    throw io_error::dependent_error(filename, exception);
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
      {material_type::metallic, pbrt_mtype::metal},
      {material_type::refractive, pbrt_mtype::glass},
      {material_type::transparent, pbrt_mtype::thinglass},
      {material_type::subsurface, pbrt_mtype::matte},
      {material_type::volume, pbrt_mtype::matte},
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
        save_shape(path_join(dirname, path), shape, true);
      });
      // save textures
      parallel_foreach(scene.textures, [&](auto& texture) {
        auto path = "textures/" + get_texture_name(scene, texture) +
                    (!texture.pixelsf.empty() ? ".hdr"s : ".png"s);
        save_texture(path_join(dirname, path), texture);
      });
    }
  } catch (const io_error& exception) {
    throw io_error::dependent_error(filename, exception);
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
