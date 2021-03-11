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

#include "ext/json.hpp"
#include "ext/stb_image.h"
#include "ext/stb_image_resize.h"
#include "ext/stb_image_write.h"
#include "ext/tinyexr.h"
#include "yocto_cli.h"
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

// List the contents of a directory
vector<string> list_directory(const string& filename) {
  auto entries = vector<string>{};
  for (auto entry : std::filesystem::directory_iterator(make_path(filename))) {
    entries.push_back(entry.path().generic_u8string());
  }
  return entries;
}

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

// Get the current directory
string path_current() { return std::filesystem::current_path().u8string(); }

}  // namespace yocto

// -----------------------------------------------------------------------------
// FILE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Opens a file with a utf8 file name
static FILE* fopen_utf8(const char* filename, const char* mode) {
#ifdef _WIN32
  auto path8 = std::filesystem::u8path(filename);
  auto wmode = std::wstring(string{mode}.begin(), string{mode}.end());
  return _wfopen(path8.c_str(), wmode.c_str());
#else
  return fopen(filename, mode);
#endif
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
// IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

template <typename T>
static T swap_endian(T value) {
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

// Pfm load
static float* load_pfm(
    const string& filename, int* width, int* height, int* components, int req) {
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
  if (!fs) return nullptr;

  // buffer
  auto buffer = array<char, 4096>{};
  auto toks   = vector<string>();

  // read magic
  if (!fgets(buffer.data(), (int)buffer.size(), fs)) return nullptr;
  toks = split_string(buffer.data());
  if (toks[0] == "Pf") {
    *components = 1;
  } else if (toks[0] == "PF") {
    *components = 3;
  } else {
    return nullptr;
  }

  // read width, height
  if (!fgets(buffer.data(), (int)buffer.size(), fs)) return nullptr;
  toks    = split_string(buffer.data());
  *width  = atoi(toks[0].c_str());
  *height = atoi(toks[1].c_str());

  // read scale
  if (!fgets(buffer.data(), (int)buffer.size(), fs)) return nullptr;
  toks   = split_string(buffer.data());
  auto s = atof(toks[0].c_str());

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
      pixels[i] = swap_endian(pixels[i]);
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
static bool save_pfm(const char* filename, int width, int height,
    int components, const float* pixels) {
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
bool load_image(const string& filename, color_image& image, string& error) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };
  auto preset_error = [filename, &error]() {
    error = filename + ": " + error;
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
    auto pixels = (float*)nullptr;
    if (LoadEXR(&pixels, &image.width, &image.height, filename.c_str(),
            nullptr) != 0)
      return read_error();
    image.linear = true;
    image.pixels = from_linear(pixels, image.width, image.height);
    free(pixels);
    return true;
  } else if (ext == ".pfm" || ext == ".PFM") {
    auto ncomp  = 0;
    auto pixels = load_pfm(
        filename.c_str(), &image.width, &image.height, &ncomp, 4);
    if (!pixels) return read_error();
    image.linear = true;
    image.pixels = from_linear(pixels, image.width, image.height);
    delete[] pixels;
    return true;
  } else if (ext == ".hdr" || ext == ".HDR") {
    auto ncomp  = 0;
    auto pixels = stbi_loadf(
        filename.c_str(), &image.width, &image.height, &ncomp, 4);
    if (!pixels) return read_error();
    image.linear = true;
    image.pixels = from_linear(pixels, image.width, image.height);
    free(pixels);
    return true;
  } else if (ext == ".png" || ext == ".PNG") {
    auto ncomp  = 0;
    auto pixels = stbi_load(
        filename.c_str(), &image.width, &image.height, &ncomp, 4);
    if (!pixels) return read_error();
    image.linear = false;
    image.pixels = from_srgb(pixels, image.width, image.height);
    free(pixels);
    return true;
  } else if (ext == ".jpg" || ext == ".JPG") {
    auto ncomp  = 0;
    auto pixels = stbi_load(
        filename.c_str(), &image.width, &image.height, &ncomp, 4);
    if (!pixels) return read_error();
    image.linear = false;
    image.pixels = from_srgb(pixels, image.width, image.height);
    free(pixels);
    return true;
  } else if (ext == ".tga" || ext == ".TGA") {
    auto ncomp  = 0;
    auto pixels = stbi_load(
        filename.c_str(), &image.width, &image.height, &ncomp, 4);
    if (!pixels) return read_error();
    image.linear = false;
    image.pixels = from_srgb(pixels, image.width, image.height);
    free(pixels);
    return true;
  } else if (ext == ".bmp" || ext == ".BMP") {
    auto ncomp  = 0;
    auto pixels = stbi_load(
        filename.c_str(), &image.width, &image.height, &ncomp, 4);
    if (!pixels) return read_error();
    image.linear = false;
    image.pixels = from_srgb(pixels, image.width, image.height);
    free(pixels);
    return true;
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    // create preset
    if (!make_image_preset(image, path_basename(filename), error))
      return preset_error();
    return true;
  } else {
    return format_error();
  }
}

// Saves an hdr image.
bool save_image(
    const string& filename, const color_image& image, string& error) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto write_error = [filename, &error]() {
    error = filename + ": write error";
    return false;
  };

  // conversion helpers
  auto to_linear = [](const color_image& image) {
    if (image.linear) return image.pixels;
    auto pixelsf = vector<vec4f>(image.pixels.size());
    srgb_to_rgb(pixelsf, image.pixels);
    return pixelsf;
  };
  auto to_srgb = [](const color_image& image) {
    auto pixelsb = vector<vec4b>(image.pixels.size());
    if (image.linear) {
      rgb_to_srgb(pixelsb, image.pixels);
    } else {
      float_to_byte(pixelsb, image.pixels);
    }
    return pixelsb;
  };

  auto ext = path_extension(filename);
  if (ext == ".hdr" || ext == ".HDR") {
    if (!stbi_write_hdr(filename.c_str(), (int)image.width, (int)image.height,
            4, (const float*)to_linear(image).data()))
      return write_error();
    return true;
  } else if (ext == ".pfm" || ext == ".PFM") {
    if (!save_pfm(filename.c_str(), image.width, image.height, 4,
            (const float*)to_linear(image).data()))
      return write_error();
    return true;
  } else if (ext == ".exr" || ext == ".EXR") {
    if (SaveEXR((const float*)to_linear(image).data(), (int)image.width,
            (int)image.height, 4, 1, filename.c_str(), nullptr) < 0)
      return write_error();
    return true;
  } else if (ext == ".png" || ext == ".PNG") {
    if (!stbi_write_png(filename.c_str(), (int)image.width, (int)image.height,
            4, (const byte*)to_srgb(image).data(), (int)image.width * 4))
      return write_error();
    return true;
  } else if (ext == ".jpg" || ext == ".JPG") {
    if (!stbi_write_jpg(filename.c_str(), (int)image.width, (int)image.height,
            4, (const byte*)to_srgb(image).data(), 75))
      return write_error();
    return true;
  } else if (ext == ".tga" || ext == ".TGA") {
    if (!stbi_write_tga(filename.c_str(), (int)image.width, (int)image.height,
            4, (const byte*)to_srgb(image).data()))
      return write_error();
    return true;
  } else if (ext == ".bmp" || ext == ".BMP") {
    if (!stbi_write_bmp(filename.c_str(), (int)image.width, (int)image.height,
            4, (const byte*)to_srgb(image).data()))
      return write_error();
    return true;
  } else {
    return format_error();
  }
}

bool make_image_preset(color_image& image, const string& type_, string& error) {
  auto type = path_basename(type_);

  auto width = 1024, height = 1024;
  if (type.find("sky") != type.npos) width = 2048;
  if (type.find("images2") != type.npos) width = 2048;
  if (type == "grid") {
    image = make_grid(width, height);
  } else if (type == "checker") {
    image = make_checker(width, height);
  } else if (type == "bumps") {
    image = make_bumps(width, height);
  } else if (type == "uvramp") {
    image = make_uvramp(width, height);
  } else if (type == "gammaramp") {
    image = make_gammaramp(width, height);
  } else if (type == "blackbodyramp") {
    image = make_blackbodyramp(width, height);
  } else if (type == "uvgrid") {
    image = make_uvgrid(width, height);
  } else if (type == "colormapramp") {
    image = make_colormapramp(width, height);
  } else if (type == "sky") {
    image = make_sunsky(
        width, height, pif / 4, 3.0, false, 1.0, 1.0, vec3f{0.7, 0.7, 0.7});
  } else if (type == "sunsky") {
    image = make_sunsky(
        width, height, pif / 4, 3.0, true, 1.0, 1.0, vec3f{0.7, 0.7, 0.7});
  } else if (type == "noise") {
    image = make_noisemap(width, height, 1);
  } else if (type == "fbm") {
    image = make_fbmmap(width, height, 1);
  } else if (type == "ridge") {
    image = make_ridgemap(width, height, 1);
  } else if (type == "turbulence") {
    image = make_turbulencemap(width, height, 1);
  } else if (type == "bump-normal") {
    image = make_bumps(width, height);
    // TODO(fabio): fix color space
    // img   = srgb_to_rgb(bump_to_normal(img, 0.05f));
  } else if (type == "images1") {
    auto sub_types = vector<string>{"grid", "uvgrid", "checker", "gammaramp",
        "bumps", "bump-normal", "noise", "fbm", "blackbodyramp"};
    auto sub_imgs  = vector<color_image>(sub_types.size());
    for (auto i = 0; i < sub_imgs.size(); i++) {
      if (!make_image_preset(sub_imgs[i], sub_types[i], error)) return false;
    }
    auto montage_size = zero2i;
    for (auto& sub_img : sub_imgs) {
      montage_size.x += sub_img.width;
      montage_size.y = max(montage_size.y, sub_img.height);
    }
    image    = make_image(montage_size.x, montage_size.y, sub_imgs[0].linear);
    auto pos = 0;
    for (auto& sub_img : sub_imgs) {
      set_region(image, sub_img, pos, 0);
      pos += sub_img.width;
    }
  } else if (type == "images2") {
    auto sub_types = vector<string>{"sky", "sunsky"};
    auto sub_imgs  = vector<color_image>(sub_types.size());
    for (auto i = 0; i < sub_imgs.size(); i++) {
      if (!make_image_preset(sub_imgs[i], sub_types[i], error)) return false;
    }
    auto montage_size = zero2i;
    for (auto& sub_img : sub_imgs) {
      montage_size.x += sub_img.width;
      montage_size.y = max(montage_size.y, sub_img.height);
    }
    image    = make_image(montage_size.x, montage_size.y, sub_imgs[0].linear);
    auto pos = 0;
    for (auto& sub_img : sub_imgs) {
      set_region(image, sub_img, pos, 0);
      pos += sub_img.width;
    }
  } else if (type == "test-floor") {
    image = make_grid(width, height);
    image = add_border(image, 0.0025);
  } else if (type == "test-grid") {
    image = make_grid(width, height);
  } else if (type == "test-checker") {
    image = make_checker(width, height);
  } else if (type == "test-bumps") {
    image = make_bumps(width, height);
  } else if (type == "test-uvramp") {
    image = make_uvramp(width, height);
  } else if (type == "test-gammaramp") {
    image = make_gammaramp(width, height);
  } else if (type == "test-blackbodyramp") {
    image = make_blackbodyramp(width, height);
  } else if (type == "test-colormapramp") {
    image = make_colormapramp(width, height);
    // TODO(fabio): fix color space
    // img   = srgb_to_rgb(img);
  } else if (type == "test-uvgrid") {
    image = make_uvgrid(width, height);
  } else if (type == "test-sky") {
    image = make_sunsky(
        width, height, pif / 4, 3.0, false, 1.0, 1.0, vec3f{0.7, 0.7, 0.7});
  } else if (type == "test-sunsky") {
    image = make_sunsky(
        width, height, pif / 4, 3.0, true, 1.0, 1.0, vec3f{0.7, 0.7, 0.7});
  } else if (type == "test-noise") {
    image = make_noisemap(width, height);
  } else if (type == "test-fbm") {
    image = make_noisemap(width, height);
  } else if (type == "test-bumps-normal") {
    image = make_bumps(width, height);
    image = bump_to_normal(image, 0.05);
  } else if (type == "test-bumps-displacement") {
    image = make_bumps(width, height);
    // TODO(fabio): fix color space
    // img   = srgb_to_rgb(img);
  } else if (type == "test-fbm-displacement") {
    image = make_fbmmap(width, height);
    // TODO(fabio): fix color space
    // img   = srgb_to_rgb(img);
  } else if (type == "test-checker-opacity") {
    image = make_checker(width, height, 1, {1, 1, 1, 1}, {0, 0, 0, 0});
  } else if (type == "test-grid-opacity") {
    image = make_grid(width, height, 1, {1, 1, 1, 1}, {0, 0, 0, 0});
  } else {
    error = "unknown preset";
    image = {};
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
    const string& filename, scene_texture& texture, string& error) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };
  auto preset_error = [filename, &error]() {
    error = filename + ": " + error;
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
  } else if (ext == ".pfm" || ext == ".PFM") {
    auto ncomp  = 0;
    auto pixels = load_pfm(
        filename.c_str(), &texture.width, &texture.height, &ncomp, 4);
    if (!pixels) return read_error();
    texture.linear  = true;
    texture.pixelsf = vector<vec4f>{
        (vec4f*)pixels, (vec4f*)pixels + texture.width * texture.height};
    delete[] pixels;
    return true;
  } else if (ext == ".hdr" || ext == ".HDR") {
    auto ncomp  = 0;
    auto pixels = stbi_loadf(
        filename.c_str(), &texture.width, &texture.height, &ncomp, 4);
    if (!pixels) return read_error();
    texture.linear  = true;
    texture.pixelsf = vector<vec4f>{
        (vec4f*)pixels, (vec4f*)pixels + texture.width * texture.height};
    free(pixels);
    return true;
  } else if (ext == ".png" || ext == ".PNG") {
    auto ncomp  = 0;
    auto pixels = stbi_load(
        filename.c_str(), &texture.width, &texture.height, &ncomp, 4);
    if (!pixels) return read_error();
    texture.linear  = false;
    texture.pixelsb = vector<vec4b>{
        (vec4b*)pixels, (vec4b*)pixels + texture.width * texture.height};
    free(pixels);
    return true;
  } else if (ext == ".jpg" || ext == ".JPG") {
    auto ncomp  = 0;
    auto pixels = stbi_load(
        filename.c_str(), &texture.width, &texture.height, &ncomp, 4);
    if (!pixels) return read_error();
    texture.linear  = false;
    texture.pixelsb = vector<vec4b>{
        (vec4b*)pixels, (vec4b*)pixels + texture.width * texture.height};
    free(pixels);
    return true;
  } else if (ext == ".tga" || ext == ".TGA") {
    auto ncomp  = 0;
    auto pixels = stbi_load(
        filename.c_str(), &texture.width, &texture.height, &ncomp, 4);
    if (!pixels) return read_error();
    texture.linear  = false;
    texture.pixelsb = vector<vec4b>{
        (vec4b*)pixels, (vec4b*)pixels + texture.width * texture.height};
    free(pixels);
    return true;
  } else if (ext == ".bmp" || ext == ".BMP") {
    auto ncomp  = 0;
    auto pixels = stbi_load(
        filename.c_str(), &texture.width, &texture.height, &ncomp, 4);
    if (!pixels) return read_error();
    texture.linear  = false;
    texture.pixelsb = vector<vec4b>{
        (vec4b*)pixels, (vec4b*)pixels + texture.width * texture.height};
    free(pixels);
    return true;
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    // create preset
    if (!make_texture_preset(texture, path_basename(filename), error))
      return preset_error();
    return true;
  } else {
    return format_error();
  }
}

// Saves an hdr image.
bool save_texture(
    const string& filename, const scene_texture& texture, string& error) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto write_error = [filename, &error]() {
    error = filename + ": write error";
    return false;
  };
  auto hdr_error = [filename, &error]() {
    error = filename + ": cannot save hdr texture to ldr file";
    return false;
  };
  auto ldr_error = [filename, &error]() {
    error = filename + ": cannot save ldr texture to hdr file";
    return false;
  };

  // check for correct handling
  if (!texture.pixelsf.empty() && is_ldr_filename(filename)) return hdr_error();
  if (!texture.pixelsb.empty() && is_hdr_filename(filename)) return ldr_error();

  auto ext = path_extension(filename);
  if (ext == ".hdr" || ext == ".HDR") {
    if (!stbi_write_hdr(filename.c_str(), (int)texture.width,
            (int)texture.height, 4, (const float*)texture.pixelsf.data()))
      return write_error();
    return true;
  } else if (ext == ".pfm" || ext == ".PFM") {
    if (!save_pfm(filename.c_str(), texture.width, texture.height, 4,
            (const float*)texture.pixelsf.data()))
      return write_error();
    return true;
  } else if (ext == ".exr" || ext == ".EXR") {
    if (SaveEXR((const float*)texture.pixelsf.data(), (int)texture.width,
            (int)texture.height, 4, 1, filename.c_str(), nullptr) < 0)
      return write_error();
    return true;
  } else if (ext == ".png" || ext == ".PNG") {
    if (!stbi_write_png(filename.c_str(), (int)texture.width,
            (int)texture.height, 4, (const byte*)texture.pixelsb.data(),
            (int)texture.width * 4))
      return write_error();
    return true;
  } else if (ext == ".jpg" || ext == ".JPG") {
    if (!stbi_write_jpg(filename.c_str(), (int)texture.width,
            (int)texture.height, 4, (const byte*)texture.pixelsb.data(), 75))
      return write_error();
    return true;
  } else if (ext == ".tga" || ext == ".TGA") {
    if (!stbi_write_tga(filename.c_str(), (int)texture.width,
            (int)texture.height, 4, (const byte*)texture.pixelsb.data()))
      return write_error();
    return true;
  } else if (ext == ".bmp" || ext == ".BMP") {
    if (!stbi_write_bmp(filename.c_str(), (int)texture.width,
            (int)texture.height, 4, (const byte*)texture.pixelsb.data()))
      return write_error();
    return true;
  } else {
    return format_error();
  }
}

bool make_texture_preset(
    scene_texture& texture, const string& type, string& error) {
  auto image = color_image{};
  if (!make_image_preset(image, type, error)) return false;
  texture = image_to_texture(image);
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load ply mesh
bool load_shape(const string& filename, scene_shape& shape, string& error,
    bool flip_texcoord) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };
  auto preset_error = [filename, &error]() {
    error = filename + ": " + error;
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
    // create preset
    if (!make_shape_preset(shape, path_basename(filename), error))
      return preset_error();
    return true;
  } else {
    return format_error();
  }
}

// Save ply mesh
bool save_shape(const string& filename, const scene_shape& shape, string& error,
    bool flip_texcoord, bool ascii) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };
  auto line_error = [filename, &error]() {
    error = filename + ": unsupported lines";
    return false;
  };
  auto point_error = [filename, &error]() {
    error = filename + ": unsupported points";
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
    return save_ply(filename, ply, error);
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
    return save_obj(filename, obj, error);
  } else if (ext == ".stl" || ext == ".STL") {
    auto stl = stl_model{};
    if (!shape.lines.empty()) return line_error();
    if (!shape.points.empty()) return point_error();
    if (!shape.triangles.empty()) {
      add_triangles(stl, shape.triangles, shape.positions, {});
    } else if (!shape.quads.empty()) {
      add_triangles(stl, quads_to_triangles(shape.quads), shape.positions, {});
    } else {
      return shape_error();
    }
    return save_stl(filename, stl, error);
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
    return save_text(filename, str, error);
  } else {
    return format_error();
  }
}

// Load ply mesh
bool load_fvshape(const string& filename, scene_fvshape& shape, string& error,
    bool flip_texcoord) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };
  auto preset_error = [filename, &error]() {
    error = filename + ": " + error;
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
    // create preset
    if (!make_fvshape_preset(shape, path_basename(filename), error))
      return preset_error();
    return true;
  } else {
    return format_error();
  }
}

// Save ply mesh
bool save_fvshape(const string& filename, const scene_fvshape& shape,
    string& error, bool flip_texcoord, bool ascii) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto shape_error = [filename, &error]() {
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
    return save_ply(filename, ply, error);
  } else if (ext == ".obj" || ext == ".OBJ") {
    auto obj = obj_shape{};
    add_positions(obj, shape.positions);
    add_normals(obj, shape.positions);
    add_texcoords(obj, shape.texcoords, flip_texcoord);
    add_fvquads(obj, shape.quadspos, shape.quadsnorm, shape.quadstexcoord, 0);
    return save_obj(filename, obj, error);
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
    return save_stl(filename, stl, error);
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
    return save_text(filename, str, error);
  } else {
    return format_error();
  }
}

// Shape presets used ofr testing.
bool make_shape_preset(scene_shape& shape, const string& type, string& error) {
  auto set_quads = [&](scene_shape&& shape_) {
    shape.quads     = shape_.quads;
    shape.positions = shape_.positions;
    shape.normals   = shape_.normals;
    shape.texcoords = shape_.texcoords;
  };
  auto set_triangles = [&](scene_shape&& shape_) {
    shape.triangles = shape_.triangles;
    shape.positions = shape_.positions;
    shape.normals   = shape_.normals;
    shape.texcoords = shape_.texcoords;
  };
  auto set_lines = [&](scene_shape&& shape_) {
    shape.lines     = shape_.lines;
    shape.positions = shape_.positions;
    shape.normals   = shape_.normals;
    shape.texcoords = shape_.texcoords;
    shape.radius    = shape_.radius;
  };
  auto set_points = [&](scene_shape&& shape_) {
    shape.points    = shape_.points;
    shape.positions = shape_.positions;
    shape.normals   = shape_.normals;
    shape.texcoords = shape_.texcoords;
    shape.radius    = shape_.radius;
  };
  auto set_fvquads = [&](scene_fvshape&& shape_) {
    shape.quads     = shape_.quadspos;
    shape.positions = shape_.positions;
    shape.normals   = shape_.normals;
    shape.texcoords = shape_.texcoords;
  };

  if (type == "default-quad") {
    set_quads(make_rect());
  } else if (type == "default-quady") {
    set_quads(make_recty());
  } else if (type == "default-cube") {
    set_quads(make_box());
  } else if (type == "default-cube-rounded") {
    set_quads(make_rounded_box());
  } else if (type == "default-sphere") {
    set_quads(make_sphere());
  } else if (type == "default-disk") {
    set_quads(make_disk());
  } else if (type == "default-disk-bulged") {
    set_quads(make_bulged_disk());
  } else if (type == "default-quad-bulged") {
    set_quads(make_bulged_rect());
  } else if (type == "default-uvsphere") {
    set_quads(make_uvsphere());
  } else if (type == "default-uvsphere-flipcap") {
    set_quads(make_capped_uvsphere());
  } else if (type == "default-uvdisk") {
    set_quads(make_uvdisk());
  } else if (type == "default-uvcylinder") {
    set_quads(make_uvcylinder());
  } else if (type == "default-uvcylinder-rounded") {
    set_quads(make_rounded_uvcylinder({32, 32, 32}));
  } else if (type == "default-geosphere") {
    set_triangles(make_geosphere());
  } else if (type == "default-floor") {
    set_quads(make_floor());
  } else if (type == "default-floor-bent") {
    set_quads(make_bent_floor());
  } else if (type == "default-matball") {
    set_quads(make_sphere());
  } else if (type == "default-hairball") {
    auto base = make_sphere(pow2(5), 0.8);
    set_lines(make_hair(base, {4, 65536}, {0.2, 0.2}, {0.002, 0.001}));
  } else if (type == "default-hairball-interior") {
    set_quads(make_sphere(pow2(5), 0.8));
  } else if (type == "default-suzanne") {
    set_quads(make_monkey());
  } else if (type == "default-cube-facevarying") {
    set_fvquads(make_fvbox());
  } else if (type == "default-sphere-facevarying") {
    set_fvquads(make_fvsphere());
  } else if (type == "default-quady-displaced") {
    set_quads(make_recty({256, 256}));
  } else if (type == "default-sphere-displaced") {
    set_quads(make_sphere(128));
  } else if (type == "test-cube") {
    set_quads(make_rounded_box(
        {32, 32, 32}, {0.075f, 0.075f, 0.075f}, {1, 1, 1}, 0.3 * 0.075f));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-uvsphere") {
    set_quads(make_uvsphere({32, 32}, 0.075));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-uvsphere-flipcap") {
    set_quads(make_capped_uvsphere({32, 32}, 0.075, {1, 1}, 0.3 * 0.075));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-sphere") {
    set_quads(make_sphere(32, 0.075f, 1));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-sphere-displaced") {
    set_quads(make_sphere(128, 0.075f, 1));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-disk") {
    set_quads(make_disk(32, 0.075f, 1));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-uvcylinder") {
    set_quads(make_rounded_uvcylinder(
        {32, 32, 32}, {0.075, 0.075}, {1, 1, 1}, 0.3 * 0.075));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-floor") {
    set_quads(make_floor({1, 1}, {2, 2}, {20, 20}));
  } else if (type == "test-quad") {
    set_quads(make_rect({1, 1}, {0.075, 0.075}, {1, 1}));
  } else if (type == "test-quady") {
    set_quads(make_recty({1, 1}, {0.075, 0.075}, {1, 1}));
  } else if (type == "test-quad-displaced") {
    set_quads(make_rect({256, 256}, {0.075, 0.075}, {1, 1}));
  } else if (type == "test-quady-displaced") {
    set_quads(make_recty({256, 256}, {0.075, 0.075}, {1, 1}));
  } else if (type == "test-matball") {
    set_quads(make_sphere(32, 0.075));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-hairball1") {
    auto base = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : base.positions) p += {0, 0.075, 0};
    set_lines(make_hair(base, {4, 65536}, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f}, {0.03, 100}));
  } else if (type == "test-hairball2") {
    auto base = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : base.positions) p += {0, 0.075, 0};
    set_lines(make_hair(base, {4, 65536}, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f}));
  } else if (type == "test-hairball3") {
    auto base = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : base.positions) p += {0, 0.075, 0};
    set_lines(make_hair(base, {4, 65536}, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f}, {0, 0}, {0.5, 128}));
  } else if (type == "test-hairball-interior") {
    set_quads(make_sphere(32, 0.075f * 0.8f, 1));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-suzanne-subdiv") {
    set_quads(make_monkey(0.075f * 0.8f));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-cube-subdiv") {
    // set_quads(make_cube( 0.075f);
    set_fvquads(make_fvcube(0.075f));
    // make_fvbox(quadspos, quadsnorm, quadstexcoord, positions, normals,
    //      texcoords, {1, 1, 1}, {0.075f, 0.075f, 0.075f});
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-arealight1") {
    set_quads(make_rect({1, 1}, {0.2, 0.2}));
  } else if (type == "test-arealight2") {
    set_quads(make_rect({1, 1}, {0.2, 0.2}));
  } else if (type == "test-largearealight1") {
    set_quads(make_rect({1, 1}, {0.4, 0.4}));
  } else if (type == "test-largearealight2") {
    set_quads(make_rect({1, 1}, {0.4, 0.4}));
  } else if (type == "test-pointlight1") {
    set_points(make_point(0));
  } else if (type == "test-pointlight2") {
    set_points(make_point(0));
  } else if (type == "test-point") {
    set_points(make_points(1));
  } else if (type == "test-points") {
    set_points(make_points(4096));
  } else if (type == "test-points-random") {
    set_points(make_random_points(4096, {0.2, 0.2, 0.2}));
  } else if (type == "test-particles") {
    set_points(make_points(4096));
  } else if (type == "test-cloth") {
    set_quads(make_rect({64, 64}, {0.2, 0.2}));
  } else if (type == "test-clothy") {
    set_quads(make_recty({64, 64}, {0.2, 0.2}));
  } else {
    error = "unknown preset";
    return false;
  }
  return true;
}

// Shape presets used for testing.
bool make_fvshape_preset(
    scene_fvshape& shape, const string& type, string& error) {
  auto set_quads = [&](scene_shape&& shape_) {
    shape.quadspos  = shape_.quads;
    shape.positions = shape_.positions;
    if (!shape_.normals.empty()) shape.quadsnorm = shape_.quads;
    shape.normals = shape_.normals;
    if (!shape_.texcoords.empty()) shape.quadstexcoord = shape_.quads;
    shape.texcoords = shape_.texcoords;
  };
  auto set_triangles = [&](scene_shape&& shape) {
    throw std::invalid_argument{"bad shape type"};
  };
  auto set_lines = [&](scene_shape&& shape) {
    throw std::invalid_argument{"bad shape type"};
  };
  auto set_points = [&](scene_shape&& shape) {
    throw std::invalid_argument{"bad shape type"};
  };
  auto set_fvquads = [&](scene_fvshape&& shape_) {
    shape.quadspos      = shape_.quadspos;
    shape.quadsnorm     = shape_.quadsnorm;
    shape.quadstexcoord = shape_.quadstexcoord;
    shape.positions     = shape_.positions;
    shape.normals       = shape_.normals;
    shape.texcoords     = shape_.texcoords;
  };

  if (type == "default-quad") {
    set_quads(make_rect());
  } else if (type == "default-quady") {
    set_quads(make_recty());
  } else if (type == "default-cube") {
    set_quads(make_box());
  } else if (type == "default-cube-rounded") {
    set_quads(make_rounded_box());
  } else if (type == "default-sphere") {
    set_quads(make_sphere());
  } else if (type == "default-disk") {
    set_quads(make_disk());
  } else if (type == "default-disk-bulged") {
    set_quads(make_bulged_disk());
  } else if (type == "default-quad-bulged") {
    set_quads(make_bulged_rect());
  } else if (type == "default-uvsphere") {
    set_quads(make_uvsphere());
  } else if (type == "default-uvsphere-flipcap") {
    set_quads(make_capped_uvsphere());
  } else if (type == "default-uvdisk") {
    set_quads(make_uvdisk());
  } else if (type == "default-uvcylinder") {
    set_quads(make_uvcylinder());
  } else if (type == "default-uvcylinder-rounded") {
    set_quads(make_rounded_uvcylinder({32, 32, 32}));
  } else if (type == "default-geosphere") {
    set_triangles(make_geosphere());
  } else if (type == "default-floor") {
    set_quads(make_floor());
  } else if (type == "default-floor-bent") {
    set_quads(make_bent_floor());
  } else if (type == "default-matball") {
    set_quads(make_sphere());
  } else if (type == "default-hairball") {
    auto base = make_sphere(pow2(5), 0.8);
    set_lines(make_hair(base, {4, 65536}, {0.2, 0.2}, {0.002, 0.001}));
  } else if (type == "default-hairball-interior") {
    set_quads(make_sphere(pow2(5), 0.8));
  } else if (type == "default-suzanne") {
    set_quads(make_monkey());
  } else if (type == "default-cube-facevarying") {
    set_fvquads(make_fvbox());
  } else if (type == "default-sphere-facevarying") {
    set_fvquads(make_fvsphere());
  } else if (type == "default-quady-displaced") {
    set_quads(make_recty({256, 256}));
  } else if (type == "default-sphere-displaced") {
    set_quads(make_sphere(128));
  } else if (type == "test-cube") {
    set_quads(make_rounded_box(
        {32, 32, 32}, {0.075f, 0.075f, 0.075f}, {1, 1, 1}, 0.3 * 0.075f));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-uvsphere") {
    set_quads(make_uvsphere({32, 32}, 0.075));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-uvsphere-flipcap") {
    set_quads(make_capped_uvsphere({32, 32}, 0.075, {1, 1}, 0.3 * 0.075));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-sphere") {
    set_quads(make_sphere(32, 0.075f, 1));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-sphere-displaced") {
    set_quads(make_sphere(128, 0.075f, 1));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-disk") {
    set_quads(make_disk(32, 0.075f, 1));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-uvcylinder") {
    set_quads(make_rounded_uvcylinder(
        {32, 32, 32}, {0.075, 0.075}, {1, 1, 1}, 0.3 * 0.075));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-floor") {
    set_quads(make_floor({1, 1}, {2, 2}, {20, 20}));
  } else if (type == "test-quad") {
    set_quads(make_rect({1, 1}, {0.075, 0.075}, {1, 1}));
  } else if (type == "test-quady") {
    set_quads(make_recty({1, 1}, {0.075, 0.075}, {1, 1}));
  } else if (type == "test-quad-displaced") {
    set_quads(make_rect({256, 256}, {0.075, 0.075}, {1, 1}));
  } else if (type == "test-quady-displaced") {
    set_quads(make_recty({256, 256}, {0.075, 0.075}, {1, 1}));
  } else if (type == "test-matball") {
    set_quads(make_sphere(32, 0.075));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-hairball1") {
    auto base = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : base.positions) p += {0, 0.075, 0};
    set_lines(make_hair(base, {4, 65536}, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f}, {0.03, 100}));
  } else if (type == "test-hairball2") {
    auto base = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : base.positions) p += {0, 0.075, 0};
    set_lines(make_hair(base, {4, 65536}, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f}));
  } else if (type == "test-hairball3") {
    auto base = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : base.positions) p += {0, 0.075, 0};
    set_lines(make_hair(base, {4, 65536}, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f}, {0, 0}, {0.5, 128}));
  } else if (type == "test-hairball-interior") {
    set_quads(make_sphere(32, 0.075f * 0.8f, 1));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-suzanne-subdiv") {
    set_quads(make_monkey(0.075f * 0.8f));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-cube-subdiv") {
    // set_quads(make_cube( 0.075f);
    set_fvquads(make_fvcube(0.075f));
    // make_fvbox(quadspos, quadsnorm, quadstexcoord, positions, normals,
    //      texcoords, {1, 1, 1}, {0.075f, 0.075f, 0.075f});
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-arealight1") {
    set_quads(make_rect({1, 1}, {0.2, 0.2}));
  } else if (type == "test-arealight2") {
    set_quads(make_rect({1, 1}, {0.2, 0.2}));
  } else if (type == "test-largearealight1") {
    set_quads(make_rect({1, 1}, {0.4, 0.4}));
  } else if (type == "test-largearealight2") {
    set_quads(make_rect({1, 1}, {0.4, 0.4}));
  } else if (type == "test-pointlight1") {
    set_points(make_point(0));
  } else if (type == "test-pointlight2") {
    set_points(make_point(0));
  } else if (type == "test-point") {
    set_points(make_points(1));
  } else if (type == "test-points") {
    set_points(make_points(4096));
  } else if (type == "test-points-random") {
    set_points(make_random_points(4096, {0.2, 0.2, 0.2}));
  } else if (type == "test-particles") {
    set_points(make_points(4096));
  } else if (type == "test-cloth") {
    set_quads(make_rect({64, 64}, {0.2, 0.2}));
  } else if (type == "test-clothy") {
    set_quads(make_recty({64, 64}, {0.2, 0.2}));
  } else {
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
    const scene_model& scene, int idx) {
  if (scene.camera_names.empty())
    return get_element_name("camera", idx, scene.cameras.size());
  return scene.camera_names[idx];
}
[[maybe_unused]] static string get_environment_name(
    const scene_model& scene, int idx) {
  if (scene.environment_names.empty())
    return get_element_name("environment", idx, scene.environments.size());
  return scene.environment_names[idx];
}
[[maybe_unused]] static string get_shape_name(
    const scene_model& scene, int idx) {
  if (scene.shape_names.empty())
    return get_element_name("shape", idx, scene.shapes.size());
  return scene.shape_names[idx];
}
[[maybe_unused]] static string get_texture_name(
    const scene_model& scene, int idx) {
  if (scene.texture_names.empty())
    return get_element_name("texture", idx, scene.textures.size());
  return scene.texture_names[idx];
}
[[maybe_unused]] static string get_instance_name(
    const scene_model& scene, int idx) {
  if (scene.instance_names.empty())
    return get_element_name("instance", idx, scene.instances.size());
  return scene.instance_names[idx];
}
[[maybe_unused]] static string get_material_name(
    const scene_model& scene, int idx) {
  if (scene.material_names.empty())
    return get_element_name("material", idx, scene.materials.size());
  return scene.material_names[idx];
}
[[maybe_unused]] static string get_subdiv_name(
    const scene_model& scene, int idx) {
  if (scene.subdiv_names.empty())
    return get_element_name("subdiv", idx, scene.subdivs.size());
  return scene.subdiv_names[idx];
}

[[maybe_unused]] static string get_camera_name(
    const scene_model& scene, const scene_camera& camera) {
  return get_camera_name(scene, (int)(&camera - scene.cameras.data()));
}
[[maybe_unused]] static string get_environment_name(
    const scene_model& scene, const scene_environment& environment) {
  return get_environment_name(
      scene, (int)(&environment - scene.environments.data()));
}
[[maybe_unused]] static string get_shape_name(
    const scene_model& scene, const scene_shape& shape) {
  return get_shape_name(scene, (int)(&shape - scene.shapes.data()));
}
[[maybe_unused]] static string get_texture_name(
    const scene_model& scene, const scene_texture& texture) {
  return get_texture_name(scene, (int)(&texture - scene.textures.data()));
}
[[maybe_unused]] static string get_instance_name(
    const scene_model& scene, const scene_instance& instance) {
  return get_instance_name(scene, (int)(&instance - scene.instances.data()));
}
[[maybe_unused]] static string get_material_name(
    const scene_model& scene, const scene_material& material) {
  return get_material_name(scene, (int)(&material - scene.materials.data()));
}
[[maybe_unused]] static string get_subdiv_name(
    const scene_model& scene, const scene_subdiv& subdiv) {
  return get_subdiv_name(scene, (int)(&subdiv - scene.subdivs.data()));
}

// Add missing cameras.
void add_missing_camera(scene_model& scene) {
  if (!scene.cameras.empty()) return;
  scene.camera_names.emplace_back("camera");
  auto& camera        = scene.cameras.emplace_back();
  camera.orthographic = false;
  camera.film         = 0.036;
  camera.aspect       = (float)16 / (float)9;
  camera.aperture     = 0;
  camera.lens         = 0.050;
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
static void add_missing_radius(scene_model& scene, float radius = 0.001f) {
  for (auto& shape : scene.shapes) {
    if (shape.points.empty() && shape.lines.empty()) continue;
    if (!shape.radius.empty()) continue;
    shape.radius.assign(shape.positions.size(), radius);
  }
}

// Add missing cameras.
void add_missing_material(scene_model& scene) {
  auto default_material = invalidid;
  for (auto& instance : scene.instances) {
    if (instance.material >= 0) continue;
    if (default_material == invalidid) {
      auto& material   = scene.materials.emplace_back();
      material.color   = {0.8, 0.8, 0.8};
      default_material = (int)scene.materials.size() - 1;
    }
    instance.material = default_material;
  }
}

// Reduce memory usage
static void trim_memory(scene_model& scene) {
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
void make_test(scene_model& scene, const test_params& params) {
  return;
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
bool make_scene_preset(scene_model& scene, const string& type, string& error) {
  if (type == "cornellbox") {
    make_cornellbox(scene);
    return true;
  } else if (type == "features1") {
    make_test(
        scene, {test_cameras_type::standard, test_environments_type::sky,
                   test_arealights_type::standard, test_floor_type::standard,
                   test_shapes_type::features1, test_materials_type::features1,
                   test_instance_name_type::material});
    return true;
  } else if (type == "features2") {
    make_test(
        scene, {test_cameras_type::standard, test_environments_type::sky,
                   test_arealights_type::standard, test_floor_type::standard,
                   test_shapes_type::features2, test_materials_type::features2,
                   test_instance_name_type::shape});
    return true;
  } else if (type == "materials1") {
    make_test(
        scene, {test_cameras_type::wide, test_environments_type::sky,
                   test_arealights_type::large, test_floor_type::standard,
                   test_shapes_type::rows, test_materials_type::materials1,
                   test_instance_name_type::material});
    return true;
  } else if (type == "materials2") {
    make_test(
        scene, {test_cameras_type::wide, test_environments_type::sky,
                   test_arealights_type::large, test_floor_type::standard,
                   test_shapes_type::rows, test_materials_type::materials2,
                   test_instance_name_type::material});
    return true;
  } else if (type == "materials3") {
    make_test(
        scene, {test_cameras_type::wide, test_environments_type::sky,
                   test_arealights_type::large, test_floor_type::standard,
                   test_shapes_type::rows, test_materials_type::materials3,
                   test_instance_name_type::material});
    return true;
  } else if (type == "materials4") {
    make_test(
        scene, {test_cameras_type::wide, test_environments_type::sky,
                   test_arealights_type::large, test_floor_type::standard,
                   test_shapes_type::rows, test_materials_type::materials4,
                   test_instance_name_type::material});
    return true;
  } else if (type == "materials5") {
    make_test(
        scene, {test_cameras_type::wide, test_environments_type::sky,
                   test_arealights_type::large, test_floor_type::standard,
                   test_shapes_type::rows, test_materials_type::materials5,
                   test_instance_name_type::material});
    return true;
  } else if (type == "shapes1") {
    make_test(scene, {test_cameras_type::standard, test_environments_type::sky,
                         test_arealights_type::large, test_floor_type::standard,
                         test_shapes_type::shapes1, test_materials_type::uvgrid,
                         test_instance_name_type::shape});
    return true;
  } else if (type == "shapes2") {
    make_test(scene, {test_cameras_type::standard, test_environments_type::sky,
                         test_arealights_type::large, test_floor_type::standard,
                         test_shapes_type::shapes2, test_materials_type::uvgrid,
                         test_instance_name_type::shape});
    return true;
  } else if (type == "shapes3") {
    make_test(scene, {test_cameras_type::standard, test_environments_type::sky,
                         test_arealights_type::large, test_floor_type::standard,
                         test_shapes_type::shapes3, test_materials_type::hair,
                         test_instance_name_type::shape});
    return true;
  } else if (type == "environments1") {
    make_test(scene,
        {test_cameras_type::standard, test_environments_type::sky,
            test_arealights_type::none, test_floor_type::standard,
            test_shapes_type::bunny_sphere, test_materials_type::plastic_metal,
            test_instance_name_type::material});
    return true;
  } else if (type == "environments2") {
    make_test(scene,
        {test_cameras_type::standard, test_environments_type::sunsky,
            test_arealights_type::none, test_floor_type::standard,
            test_shapes_type::bunny_sphere, test_materials_type::plastic_metal,
            test_instance_name_type::material});
    return true;
  } else if (type == "arealights1") {
    make_test(scene,
        {test_cameras_type::standard, test_environments_type::none,
            test_arealights_type::standard, test_floor_type::standard,
            test_shapes_type::bunny_sphere, test_materials_type::plastic_metal,
            test_instance_name_type::material});
    return true;
  } else {
    error = "unknown preset";
    return false;
  }
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GENERIC SCENE LOADING
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a scene in the builtin JSON format.
static bool load_json_scene(
    const string& filename, scene_model& scene, string& error, bool noparallel);
static bool save_json_scene(const string& filename, const scene_model& scene,
    string& error, bool noparallel);

// Load/save a scene from/to OBJ.
static bool load_obj_scene(
    const string& filename, scene_model& scene, string& error, bool noparallel);
static bool save_obj_scene(const string& filename, const scene_model& scene,
    string& error, bool noparallel);

// Load/save a scene from/to PLY. Loads/saves only one mesh with no other data.
static bool load_ply_scene(
    const string& filename, scene_model& scene, string& error, bool noparallel);
static bool save_ply_scene(const string& filename, const scene_model& scene,
    string& error, bool noparallel);

// Load/save a scene from/to STL. Loads/saves only one mesh with no other data.
static bool load_stl_scene(
    const string& filename, scene_model& scene, string& error, bool noparallel);
static bool save_stl_scene(const string& filename, const scene_model& scene,
    string& error, bool noparallel);

// Load/save a scene from/to glTF.
static bool load_gltf_scene(
    const string& filename, scene_model& scene, string& error, bool noparallel);
static bool save_gltf_scene(const string& filename, const scene_model& scene,
    string& error, bool noparallel);

// Load/save a scene from/to pbrt-> This is not robust at all and only
// works on scene that have been previously adapted since the two renderers
// are too different to match.
static bool load_pbrt_scene(
    const string& filename, scene_model& scene, string& error, bool noparallel);
static bool save_pbrt_scene(const string& filename, const scene_model& scene,
    string& error, bool noparallel);

// Load/save a scene preset.
static bool load_preset_scene(
    const string& filename, scene_model& scene, string& error, bool noparallel);

// Load a scene
bool load_scene(const string& filename, scene_model& scene, string& error,
    bool noparallel) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };

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
    return load_preset_scene(filename, scene, error, noparallel);
  } else {
    return format_error();
  }
}

// Save a scene
bool save_scene(const string& filename, const scene_model& scene, string& error,
    bool noparallel) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };

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
    return format_error();
  }
}

// Load/save a scene preset.
static bool load_preset_scene(const string& filename, scene_model& scene,
    string& error, bool noparallel) {
  auto preset_error = [filename, &error]() {
    error = filename + ": " + error;
    return false;
  };

  // make preset
  if (!make_scene_preset(scene, path_basename(filename), error))
    return preset_error();

  // done
  return true;
}

// Make missing scene directories
bool make_scene_directories(
    const string& filename, const scene_model& scene, string& error) {
  // make a directory if needed
  if (!make_directory(path_dirname(filename), error)) return false;
  if (!scene.shapes.empty()) {
    if (!make_directory(path_join(path_dirname(filename), "shapes"), error))
      return false;
  }
  if (!scene.textures.empty()) {
    if (!make_directory(path_join(path_dirname(filename), "textures"), error))
      return false;
  }
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// INDIVIDUAL ELEMENTS
// -----------------------------------------------------------------------------
namespace yocto {

// load instances
static bool load_instance(
    const string& filename, vector<frame3f>& frames, string& error) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = ply_model{};
    if (!load_ply(filename, ply, error)) return false;
    get_values(ply, "instance",
        {"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz", "ox", "oy",
            "oz"},
        frames);
    return true;
  } else {
    return format_error();
  }
}

// save instances
bool save_instance(const string& filename, const vector<frame3f>& frames,
    string& error, bool ascii = false) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = ply_model{};
    add_values(ply, "instance",
        {"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz", "ox", "oy",
            "oz"},
        frames);
    return save_ply(filename, ply, error);
  } else {
    return format_error();
  }
}

// load subdiv
bool load_subdiv(const string& filename, scene_subdiv& subdiv, string& error) {
  auto lsubdiv = scene_fvshape{};
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
    const string& filename, const scene_subdiv& subdiv, string& error) {
  auto ssubdiv          = scene_fvshape{};
  ssubdiv.quadspos      = subdiv.quadspos;
  ssubdiv.quadsnorm     = subdiv.quadsnorm;
  ssubdiv.quadstexcoord = subdiv.quadstexcoord;
  ssubdiv.positions     = subdiv.positions;
  ssubdiv.normals       = subdiv.normals;
  ssubdiv.texcoords     = subdiv.texcoords;
  return save_fvshape(filename, ssubdiv, error, true);
}

// save binary shape
static bool save_binshape(
    const string& filename, const scene_shape& shape, string& error) {
  auto open_error = [filename, &error]() {
    error = filename + ": file not found";
    return false;
  };
  auto write_error = [filename, &error]() {
    error = filename + ": write error";
    return false;
  };

  auto write_values = [](FILE* fs, const auto& values) -> bool {
    return fwrite(values.data(), sizeof(values.front()), values.size(), fs) ==
           values.size();
  };

  auto fs       = fopen_utf8(filename.c_str(), "wb");
  auto fs_guard = unique_ptr<FILE, int (*)(FILE*)>(fs, &fclose);
  if (!fs) return open_error();

  if (!write_values(fs, shape.positions)) return write_error();
  if (!write_values(fs, shape.normals)) return write_error();
  if (!write_values(fs, shape.texcoords)) return write_error();
  if (!write_values(fs, shape.colors)) return write_error();
  if (!write_values(fs, shape.radius)) return write_error();
  if (!write_values(fs, shape.points)) return write_error();
  if (!write_values(fs, shape.lines)) return write_error();
  if (!write_values(fs, shape.triangles)) return write_error();
  if (!write_values(fs, quads_to_triangles(shape.quads))) return write_error();

  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON SUPPORT
// -----------------------------------------------------------------------------
namespace yocto {

using njson = nlohmann::json;
using std::array;

// load/save json
[[maybe_unused]] static bool load_json(
    const string& filename, njson& json, string& error) {
  // error helpers
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error in json";
    return false;
  };
  auto text = ""s;
  if (!load_text(filename, text, error)) return false;
  try {
    json = njson::parse(text);
    return true;
  } catch (std::exception&) {
    return parse_error();
  }
}

[[maybe_unused]] static bool save_json(
    const string& filename, const njson& json, string& error) {
  return save_text(filename, json.dump(2), error);
}

// support for json conversions
inline void to_json(njson& j, const vec3f& value) {
  nlohmann::to_json(j, (const array<float, 3>&)value);
}
inline void to_json(njson& j, const vec4f& value) {
  nlohmann::to_json(j, (const array<float, 4>&)value);
}
inline void to_json(njson& j, const frame3f& value) {
  nlohmann::to_json(j, (const array<float, 12>&)value);
}
inline void to_json(njson& j, const mat4f& value) {
  nlohmann::to_json(j, (const array<float, 16>&)value);
}

inline void from_json(const njson& j, vec3f& value) {
  nlohmann::from_json(j, (array<float, 3>&)value);
}
inline void from_json(const njson& j, vec4f& value) {
  nlohmann::from_json(j, (array<float, 4>&)value);
}
inline void from_json(const njson& j, mat3f& value) {
  nlohmann::from_json(j, (array<float, 9>&)value);
}
inline void from_json(const njson& j, frame3f& value) {
  nlohmann::from_json(j, (array<float, 12>&)value);
}
inline void from_json(const njson& j, mat4f& value) {
  nlohmann::from_json(j, (array<float, 16>&)value);
}

inline void to_json(njson& j, scene_material_type value) {
  j = scene_material_names.at((int)value);
}
inline void from_json(const njson& j, scene_material_type& value) {
  auto values = j.get<string>();
  auto pos    = std::find(
      scene_material_names.begin(), scene_material_names.end(), values);
  if (pos == scene_material_names.end())
    throw std::invalid_argument{"unknown value"};
  value = (scene_material_type)(pos - scene_material_names.begin());
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load a scene in the builtin JSON format.
static bool load_json_scene(const string& filename, scene_model& scene,
    string& error, bool noparallel) {
  auto json_error = [filename]() {
    // error does not need setting
    return false;
  };
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
    error = filename + ": unknow key at " + path;
    return false;
  };
  auto material_error = [filename, &error](const string& name) {
    error = filename + ": missing material " + string{name};
    return false;
  };
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // open file
  auto js = njson{};
  if (!load_json(filename, js, error)) return json_error();

  // parse json value
  auto get_value = [](const njson& js, auto& value) -> bool {
    try {
      value = js;
      return true;
    } catch (...) {
      return false;
    }
  };

  // parse json reference
  auto shape_map = unordered_map<string, int>{};
  auto get_shape = [&scene, &shape_map, &get_value](
                       const njson& js, int& value) -> bool {
    auto name = ""s;
    if (!get_value(js, name)) return false;
    auto it = shape_map.find(name);
    if (it != shape_map.end()) {
      value = it->second;
      return true;
    }
    scene.shape_names.emplace_back(name);
    scene.shapes.emplace_back();
    auto shape_id   = (int)scene.shapes.size() - 1;
    shape_map[name] = shape_id;
    value           = shape_id;
    return true;
  };

  // parse json reference
  auto material_map = unordered_map<string, int>{};
  auto material_set = vector<bool>{};
  auto get_material = [&scene, &material_map, &material_set, &get_value](
                          const njson& js, int& value) -> bool {
    auto name = ""s;
    if (!get_value(js, name)) return false;
    auto it = material_map.find(name);
    if (it != material_map.end()) {
      value = it->second;
      return true;
    }
    scene.material_names.emplace_back(name);
    scene.materials.emplace_back();
    auto material_id   = (int)scene.materials.size() - 1;
    material_map[name] = material_id;
    value              = material_id;
    material_set.push_back(false);
    return true;
  };

  // parse json reference
  auto texture_map = unordered_map<string, int>{};
  auto get_texture = [&scene, &texture_map, &get_value](
                         const njson& js, int& value) -> bool {
    auto name = ""s;
    if (!get_value(js, name)) return false;
    auto it = texture_map.find(name);
    if (it != texture_map.end()) {
      value = it->second;
      return true;
    }
    scene.texture_names.emplace_back(name);
    scene.textures.emplace_back();
    auto texture_id   = (int)scene.textures.size() - 1;
    texture_map[name] = texture_id;
    value             = texture_id;
    return true;
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
  auto instance_ply      = unordered_map<int, ply_instance_handle>{};
  auto get_ply_instances = [&scene, &ply_instances, &ply_instances_names,
                               &ply_instance_map, &instance_ply,
                               &get_value](const njson& js,
                               const scene_instance&    instance) -> bool {
    auto name = ""s;
    if (!get_value(js, name)) return false;
    if (name.empty()) return true;
    auto instance_id = (int)(&instance - scene.instances.data());
    auto it          = ply_instance_map.find(name);
    if (it != ply_instance_map.end()) {
      instance_ply[instance_id] = it->second;
      return true;
    }
    ply_instances_names.emplace_back(name);
    ply_instances.emplace_back(ply_instance());
    auto ply_instance_id      = (int)ply_instances.size() - 1;
    ply_instance_map[name]    = ply_instance_id;
    instance_ply[instance_id] = ply_instance_id;
    return true;
  };
  auto get_ply_instance_name = [&ply_instances, &ply_instances_names](
                                   const scene_model&  scene,
                                   const ply_instance& instance) -> string {
    return ply_instances_names[&instance - ply_instances.data()];
  };

  // helper for iteration
  auto iterate_object = [](const njson& js) { return js.items(); };
  auto check_object   = [](const njson& js) { return js.is_object(); };

  // loop over external dictionaries
  for (auto& [gname, group] : iterate_object(js)) {
    if (gname == "asset") {
      if (!check_object(group)) return parse_error(gname);
      for (auto& [key, value] : iterate_object(group)) {
        if (key == "copyright") {
          if (!get_value(value, scene.copyright))
            return parse_error(gname, key);
        } else if (key == "generator") {
          // skip
        } else {
          return key_error(gname, key);
        }
      }
    } else if (gname == "cameras") {
      if (!check_object(group)) return parse_error(gname);
      for (auto& [name, element] : iterate_object(group)) {
        if (!check_object(element)) return parse_error(gname, name);
        scene.camera_names.emplace_back(name);
        auto& camera = scene.cameras.emplace_back();
        for (auto& [key, value] : iterate_object(element)) {
          if (key == "frame") {
            if (!get_value(value, camera.frame))
              return parse_error(gname, name, key);
          } else if (key == "orthographic") {
            if (!get_value(value, camera.orthographic))
              return parse_error(gname, name, key);
          } else if (key == "ortho") {
            // backward compatibility
            if (!get_value(value, camera.orthographic))
              return parse_error(gname, name, key);
          } else if (key == "lens") {
            if (!get_value(value, camera.lens))
              return parse_error(gname, name, key);
          } else if (key == "aspect") {
            if (!get_value(value, camera.aspect))
              return parse_error(gname, name, key);
          } else if (key == "film") {
            if (!get_value(value, camera.film))
              return parse_error(gname, name, key);
          } else if (key == "focus") {
            if (!get_value(value, camera.focus))
              return parse_error(gname, name, key);
          } else if (key == "aperture") {
            if (!get_value(value, camera.aperture))
              return parse_error(gname, name, key);
          } else if (key == "lookat") {
            if (!get_value(value, (mat3f&)camera.frame))
              return parse_error(gname, name, key);
            camera.focus = length(camera.frame.x - camera.frame.y);
            camera.frame = lookat_frame(
                camera.frame.x, camera.frame.y, camera.frame.z);
          } else {
            return key_error(gname, name, key);
          }
        }
      }
    } else if (gname == "environments") {
      if (!check_object(group)) return parse_error(gname);
      for (auto& [name, element] : iterate_object(group)) {
        if (!check_object(element)) return parse_error(gname, name);
        scene.environment_names.emplace_back(name);
        auto& environment = scene.environments.emplace_back();
        for (auto& [key, value] : iterate_object(element)) {
          if (key == "frame") {
            if (!get_value(value, environment.frame))
              return parse_error(gname, name, key);
          } else if (key == "emission") {
            if (!get_value(value, environment.emission))
              return parse_error(gname, name, key);
          } else if (key == "emission_tex") {
            if (!get_texture(value, environment.emission_tex))
              return parse_error(gname, name, key);
          } else if (key == "lookat") {
            if (!get_value(value, (mat3f&)environment.frame))
              return parse_error(gname, name, key);
            environment.frame = lookat_frame(environment.frame.x,
                environment.frame.y, environment.frame.z, true);
          } else {
            return key_error(gname, name, key);
          }
        }
      }
    } else if (gname == "materials") {
      if (!check_object(group)) return parse_error(gname);
      for (auto& [name, element] : iterate_object(group)) {
        if (!check_object(element)) return parse_error(gname, name);
        if (material_map.find(name) == material_map.end()) {
          scene.material_names.emplace_back(name);
          scene.materials.emplace_back();
          material_map[name] = (int)scene.materials.size() - 1;
          material_set.push_back(false);
        }
        auto& material = scene.materials.at(material_map.at(name));
        material_set[&material - &scene.materials.front()] = true;
        for (auto& [key, value] : iterate_object(element)) {
          if (key == "type") {
            if (!get_value(value, material.type))
              return parse_error(gname, name, key);
          } else if (key == "emission") {
            if (!get_value(value, material.emission))
              return parse_error(gname, name, key);
          } else if (key == "color") {
            if (!get_value(value, material.color))
              return parse_error(gname, name, key);
          } else if (key == "metallic") {
            if (!get_value(value, material.metallic))
              return parse_error(gname, name, key);
          } else if (key == "roughness") {
            if (!get_value(value, material.roughness))
              return parse_error(gname, name, key);
          } else if (key == "ior") {
            if (!get_value(value, material.ior))
              return parse_error(gname, name, key);
          } else if (key == "trdepth") {
            if (!get_value(value, material.trdepth))
              return parse_error(gname, name, key);
          } else if (key == "scattering") {
            if (!get_value(value, material.scattering))
              return parse_error(gname, name, key);
          } else if (key == "scanisotropy") {
            if (!get_value(value, material.scanisotropy))
              return parse_error(gname, name, key);
          } else if (key == "opacity") {
            if (!get_value(value, material.opacity))
              return parse_error(gname, name, key);
          } else if (key == "emission_tex") {
            if (!get_texture(value, material.emission_tex))
              return parse_error(gname, name, key);
          } else if (key == "color_tex") {
            if (!get_texture(value, material.color_tex))
              return parse_error(gname, name, key);
          } else if (key == "roughness_tex") {
            if (!get_texture(value, material.roughness_tex))
              return parse_error(gname, name, key);
          } else if (key == "scattering_tex") {
            if (!get_texture(value, material.scattering_tex))
              return parse_error(gname, name, key);
          } else if (key == "normal_tex") {
            if (!get_texture(value, material.normal_tex))
              return parse_error(gname, name, key);
          } else {
            return key_error(gname, name, key);
          }
        }
      }
    } else if (gname == "instances" || gname == "objects") {
      if (!check_object(group)) return parse_error(gname);
      for (auto& [name, element] : iterate_object(group)) {
        if (!check_object(element)) return parse_error(gname, name);
        scene.instance_names.emplace_back(name);
        auto& instance = scene.instances.emplace_back();
        for (auto [key, value] : iterate_object(element)) {
          if (key == "frame") {
            if (!get_value(value, instance.frame))
              return parse_error(gname, name, key);
          } else if (key == "lookat") {
            if (!get_value(value, (mat3f&)instance.frame))
              return parse_error(gname, name, key);
            instance.frame = lookat_frame(
                instance.frame.x, instance.frame.y, instance.frame.z, true);
          } else if (key == "material") {
            if (!get_material(value, instance.material))
              return parse_error(gname, name, key);
          } else if (key == "shape") {
            if (!get_shape(value, instance.shape))
              return parse_error(gname, name, key);
          } else if (key == "instance") {
            if (!get_ply_instances(value, instance))
              return parse_error(gname, name, key);
          } else {
            return key_error(gname, name, key);
          }
        }
      }
    } else if (gname == "subdivs") {
      if (!check_object(group)) return parse_error(gname);
      for (auto& [name, element] : iterate_object(group)) {
        if (!check_object(element)) return parse_error(gname, name);
        scene.subdiv_names.emplace_back(name);
        auto& subdiv = scene.subdivs.emplace_back();
        for (auto& [key, value] : iterate_object(element)) {
          if (key == "shape") {
            if (!get_shape(value, subdiv.shape))
              return parse_error(gname, name, key);
          } else if (key == "subdivisions") {
            if (!get_value(value, subdiv.subdivisions))
              return parse_error(gname, name, key);
          } else if (key == "catmullcark") {
            if (!get_value(value, subdiv.catmullclark))
              return parse_error(gname, name, key);
          } else if (key == "smooth") {
            if (!get_value(value, subdiv.smooth))
              return parse_error(gname, name, key);
          } else if (key == "displacement") {
            if (!get_value(value, subdiv.displacement))
              return parse_error(gname, name, key);
          } else if (key == "displacement_tex") {
            if (!get_texture(value, subdiv.displacement_tex))
              return parse_error(gname, name, key);
          } else {
            return key_error(gname, name, key);
          }
        }
      }
    } else {
      return key_error(gname);
    }
  }

  // check materials
  for (auto& material : scene.materials) {
    if (!material_set[&material - &scene.materials.front()])
      return material_error(get_material_name(scene, material));
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
  if (noparallel) {
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
    // helpers
    auto mutex = std::mutex{};
    // load shapes
    parallel_foreach(scene.shapes, [&](auto& shape) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
      }
      auto path = find_path(
          get_shape_name(scene, shape), "shapes", {".ply", ".obj"});
      auto err = string{};
      if (!load_shape(path_join(dirname, path), shape, err, true)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
    // load subdivs
    parallel_foreach(scene.subdivs, [&](auto& subdiv) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
      }
      auto path = find_path(
          get_subdiv_name(scene, subdiv), "subdivs", {".ply", ".obj"});
      auto err = string{};
      if (!load_subdiv(path_join(dirname, path), subdiv, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
    // load textures
    parallel_foreach(scene.textures, [&](auto& texture) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
      }
      auto path = find_path(get_texture_name(scene, texture), "textures",
          {".hdr", ".exr", ".png", ".jpg"});
      auto err  = string{};
      if (!load_texture(path_join(dirname, path), texture, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
    // load instances
    parallel_foreach(ply_instances, [&](auto& ply_instance) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
      }
      auto path = find_path(
          get_ply_instance_name(scene, ply_instance), "instances", {".ply"});
      auto err = string{};
      if (!load_instance(path_join(dirname, path), ply_instance.frames, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
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

// Save a scene in the builtin JSON format.
static bool save_json_scene(const string& filename, const scene_model& scene,
    string& error, bool noparallel) {
  auto conversion_error = [filename, &error](const string& message) {
    // should never happen
    throw std::runtime_error{"programmer error"};
    error = filename + ": conversion error (" + message + ")";
    return false;
  };
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // helpers to handel old code paths
  auto insert_object = [](njson& js, const string& name) -> njson& {
    js[name] = njson::object();
    return js[name];
  };
  auto insert_value = [](njson& js, const string& name, const auto& value) {
    js[name] = value;
  };

  // save json file
  auto js = njson::object();

  // asset
  {
    auto& element = insert_object(js, "asset");
    if (!scene.copyright.empty()) {
      insert_value(element, "copyright", scene.copyright);
    }
    insert_value(element, "generator",
        "Yocto/GL - https://github.com/xelatihy/yocto-gl");
  }

  auto def_cam = sceneio_camera{};
  if (!scene.cameras.empty()) {
    auto& group = insert_object(js, "cameras");
    for (auto& camera : scene.cameras) {
      auto& element = insert_object(group, get_camera_name(scene, camera));
      if (camera.frame != def_cam.frame) {
        insert_value(element, "frame", camera.frame);
      }
      if (camera.orthographic != def_cam.orthographic) {
        insert_value(element, "orthographic", camera.orthographic);
      }
      if (camera.lens != def_cam.lens) {
        insert_value(element, "lens", camera.lens);
      }
      if (camera.aspect != def_cam.aspect) {
        insert_value(element, "aspect", camera.aspect);
      }
      if (camera.film != def_cam.film) {
        insert_value(element, "film", camera.film);
      }
      if (camera.focus != def_cam.focus) {
        insert_value(element, "focus", camera.focus);
      }
      if (camera.aperture != def_cam.aperture) {
        insert_value(element, "aperture", camera.aperture);
      }
    }
  }

  auto def_env = sceneio_environment{};
  if (!scene.environments.empty()) {
    auto& group = insert_object(js, "environments");
    for (auto& environment : scene.environments) {
      auto& element = insert_object(
          group, get_environment_name(scene, environment));
      if (environment.frame != def_env.frame) {
        insert_value(element, "frame", environment.frame);
      }
      if (environment.emission != def_env.emission) {
        insert_value(element, "emission", environment.emission);
      }
      if (environment.emission_tex != invalidid) {
        insert_value(element, "emission_tex",
            get_texture_name(scene, environment.emission_tex));
      }
    }
  }

  auto def_material = sceneio_material{};
  if (!scene.materials.empty()) {
    auto& group = insert_object(js, "materials");
    for (auto& material : scene.materials) {
      auto& element = insert_object(group, get_material_name(scene, material));
      if (material.type != def_material.type) {
        insert_value(element, "type", material.type);
      }
      if (material.emission != def_material.emission) {
        insert_value(element, "emission", material.emission);
      }
      if (material.color != def_material.color) {
        insert_value(element, "color", material.color);
      }
      if (material.metallic != def_material.metallic) {
        insert_value(element, "metallic", material.metallic);
      }
      if (material.roughness != def_material.roughness) {
        insert_value(element, "roughness", material.roughness);
      }
      if (material.ior != def_material.ior) {
        insert_value(element, "ior", material.ior);
      }
      if (material.trdepth != def_material.trdepth) {
        insert_value(element, "trdepth", material.trdepth);
      }
      if (material.scattering != def_material.scattering) {
        insert_value(element, "scattering", material.scattering);
      }
      if (material.scanisotropy != def_material.scanisotropy) {
        insert_value(element, "scanisotropy", material.scanisotropy);
      }
      if (material.opacity != def_material.opacity) {
        insert_value(element, "opacity", material.opacity);
      }
      if (material.emission_tex != invalidid) {
        insert_value(element, "emission_tex",
            get_texture_name(scene, material.emission_tex));
      }
      if (material.color_tex != invalidid) {
        insert_value(
            element, "color_tex", get_texture_name(scene, material.color_tex));
      }
      if (material.roughness_tex != invalidid) {
        insert_value(element, "roughness_tex",
            get_texture_name(scene, material.roughness_tex));
      }
      if (material.scattering_tex != invalidid) {
        insert_value(element, "scattering_tex",
            get_texture_name(scene, material.scattering_tex));
      }
      if (material.normal_tex != invalidid) {
        insert_value(element, "normal_tex",
            get_texture_name(scene, material.normal_tex));
      }
    }
  }

  auto def_instance = sceneio_instance{};
  if (!scene.instances.empty()) {
    auto& group = insert_object(js, "instances");
    for (auto& instance : scene.instances) {
      auto& element = insert_object(group, get_instance_name(scene, instance));
      if (instance.frame != def_instance.frame) {
        insert_value(element, "frame", instance.frame);
      }
      if (instance.shape != invalidid) {
        insert_value(element, "shape", get_shape_name(scene, instance.shape));
      }
      if (instance.material != invalidid) {
        insert_value(
            element, "material", get_material_name(scene, instance.material));
      }
    }
  }

  auto def_subdiv = scene_subdiv{};
  if (!scene.subdivs.empty()) {
    auto& group = insert_object(js, "subdivs");
    for (auto& subdiv : scene.subdivs) {
      auto& element = insert_object(group, get_subdiv_name(scene, subdiv));
      if (subdiv.shape != invalidid) {
        insert_value(element, "shape", get_shape_name(scene, subdiv.shape));
      }
      if (subdiv.subdivisions != def_subdiv.subdivisions) {
        insert_value(element, "subdivisions", subdiv.subdivisions);
      }
      if (subdiv.catmullclark != def_subdiv.catmullclark) {
        insert_value(element, "catmullclark", subdiv.catmullclark);
      }
      if (subdiv.smooth != def_subdiv.smooth) {
        insert_value(element, "smooth", subdiv.smooth);
      }
      if (subdiv.displacement != def_subdiv.displacement) {
        insert_value(element, "displacement", subdiv.displacement);
      }
      if (subdiv.displacement_tex != invalidid) {
        insert_value(element, "displacement_tex",
            get_texture_name(scene, subdiv.displacement_tex));
      }
    }
  }

  // save json
  if (!save_json(filename, js, error)) return false;

  // dirname
  auto dirname = path_dirname(filename);

  if (noparallel) {
    // save shapes
    for (auto& shape : scene.shapes) {
      auto path = "shapes/" + get_shape_name(scene, shape) + ".ply";
      if (!save_shape(path_join(dirname, path), shape, error, true))
        return dependent_error();
    }

    // save subdiv
    for (auto& subdiv : scene.subdivs) {
      auto path = "subdivs/" + get_subdiv_name(scene, subdiv) + ".obj";
      if (!save_subdiv(path_join(dirname, path), subdiv, error))
        return dependent_error();
    }

    // save textures
    for (auto& texture : scene.textures) {
      auto path = "textures/" + get_texture_name(scene, texture) +
                  (!texture.pixelsf.empty() ? ".hdr"s : ".png"s);
      if (!save_texture(path_join(dirname, path), texture, error))
        return dependent_error();
    }
  } else {
    // mutex
    auto mutex = std::mutex{};
    // save shapes
    parallel_foreach(scene.shapes, [&](auto& shape) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
      }
      auto path = "shapes/" + get_shape_name(scene, shape) + ".ply";
      auto err  = string{};
      if (!save_shape(path_join(dirname, path), shape, err, true)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
    // save subdivs
    parallel_foreach(scene.subdivs, [&](auto& subdiv) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
      }
      auto path = "subdivs/" + get_subdiv_name(scene, subdiv) + ".obj";
      auto err  = string{};
      if (!save_subdiv(path_join(dirname, path), subdiv, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
    // save textures
    parallel_foreach(scene.textures, [&](auto& texture) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
      }
      auto path = "textures/" + get_texture_name(scene, texture) +
                  (!texture.pixelsf.empty() ? ".hdr"s : ".png"s);
      auto err = string{};
      if (!save_texture(path_join(dirname, path), texture, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
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
static bool load_obj_scene(const string& filename, scene_model& scene,
    string& error, bool noparallel) {
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };
  auto material_error = [filename, &error](const string& name) {
    error = filename + ": missing material " + name;
    return false;
  };
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // load obj
  auto obj = obj_model{};
  if (!load_obj(filename, obj, error, false, true)) return false;

  // convert cameras
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
  for (auto& omaterial : obj.materials) {
    auto& material        = scene.materials.emplace_back();
    material.type         = scene_material_type::gltfpbr;
    material.emission     = omaterial.emission;
    material.emission_tex = omaterial.emission_tex;
    if (max(omaterial.transmission) > 0.1) {
      material.type      = scene_material_type::transparent;
      material.color     = omaterial.transmission;
      material.color_tex = omaterial.transmission_tex;
    } else if (max(omaterial.specular) > 0.2) {
      material.type      = scene_material_type::metallic;
      material.color     = omaterial.specular;
      material.color_tex = omaterial.specular_tex;
    } else if (max(omaterial.specular) > 0) {
      material.type      = scene_material_type::glossy;
      material.color     = omaterial.diffuse;
      material.color_tex = omaterial.diffuse_tex;
    } else {
      material.type      = scene_material_type::matte;
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
  for (auto& oenvironment : obj.environments) {
    auto& environment        = scene.environments.emplace_back();
    environment.frame        = oenvironment.frame;
    environment.emission     = oenvironment.emission;
    environment.emission_tex = oenvironment.emission_tex;
  }

  // dirname
  auto dirname = path_dirname(filename);

  if (noparallel) {
    // load textures
    for (auto& texture : scene.textures) {
      auto& path = texture_paths[&texture - &scene.textures.front()];
      if (!load_texture(path_join(dirname, path), texture, error))
        return dependent_error();
    }
  } else {
    // mutex
    auto mutex = std::mutex{};
    // load textures
    parallel_foreach(scene.textures, [&](auto& texture) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
      }
      auto& path = texture_paths[&texture - &scene.textures.front()];
      auto  err  = string{};
      if (!load_texture(path_join(dirname, path), texture, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
  }

  // fix scene
  add_missing_camera(scene);
  add_missing_radius(scene);

  // done
  return true;
}

static bool save_obj_scene(const string& filename, const scene_model& scene,
    string& error, bool noparallel) {
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

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
  auto dirname = path_dirname(filename);

  if (noparallel) {
    // save textures
    for (auto& texture : scene.textures) {
      auto path = "textures/" + get_texture_name(scene, texture) +
                  (!texture.pixelsf.empty() ? ".hdr"s : ".png"s);
      if (!save_texture(path_join(dirname, path), texture, error))
        return dependent_error();
    }
  } else {
    // mutex
    auto mutex = std::mutex{};
    // save textures
    parallel_foreach(scene.textures, [&](auto& texture) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
      }
      auto path = "textures/" + get_texture_name(scene, texture) +
                  (!texture.pixelsf.empty() ? ".hdr"s : ".png"s);
      auto err = string{};
      if (!save_texture(path_join(dirname, path), texture, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
  }

  // done
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PLY CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static bool load_ply_scene(const string& filename, scene_model& scene,
    string& error, bool noparallel) {
  // load ply mesh
  auto& shape = scene.shapes.emplace_back();
  if (!load_shape(filename, shape, error, true)) return false;

  // create instance
  auto& instance = scene.instances.emplace_back();
  instance.shape = (int)scene.shapes.size() - 1;

  // fix scene
  add_missing_material(scene);
  add_missing_camera(scene);
  add_missing_radius(scene);

  // done
  return true;
}

static bool save_ply_scene(const string& filename, const scene_model& scene,
    string& error, bool noparallel) {
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };

  if (scene.shapes.empty()) return shape_error();

  // save shape
  auto& shape = scene.shapes.front();
  if (!save_shape(filename, shape, error, false)) return false;

  // done
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// STL CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static bool load_stl_scene(const string& filename, scene_model& scene,
    string& error, bool noparallel) {
  // load stl mesh
  auto& shape = scene.shapes.emplace_back();
  if (!load_shape(filename, shape, error, true)) return false;

  // create instance
  auto& instance = scene.instances.emplace_back();
  instance.shape = (int)scene.shapes.size() - 1;

  // fix scene
  add_missing_material(scene);
  add_missing_camera(scene);
  add_missing_radius(scene);

  // done
  return true;
}

static bool save_stl_scene(const string& filename, const scene_model& scene,
    string& error, bool noparallel) {
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };

  if (scene.shapes.empty()) return shape_error();

  // save shape
  auto& shape = scene.shapes.front();
  if (!save_shape(filename, shape, error, false)) return false;

  // done
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GLTF CONVESION
// -----------------------------------------------------------------------------
namespace yocto {

// Load a scene
static bool load_gltf_scene(const string& filename, scene_model& scene,
    string& error, bool noparallel) {
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };
  auto primitive_error = [filename, &error]() {
    error = filename + ": primitive error";
    return false;
  };
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error";
    return false;
  };
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // load gltf
  auto gltf = njson{};
  if (!load_json(filename, gltf, error)) return false;

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
  auto dirname = path_dirname(filename);

  if (noparallel) {
    // load buffers
    for (auto& buffer : buffers) {
      auto& path = buffers_paths[&buffer - &buffers.front()];
      if (!load_binary(path_join(dirname, path), buffer, error))
        return dependent_error();
    }
  } else {
    // mutex
    auto mutex = std::mutex{};
    // load buffers
    parallel_foreach(buffers, [&](auto& buffer) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
      }
      auto& path = buffers_paths[&buffer - &buffers.front()];
      auto  err  = string{};
      if (!load_binary(path_join(dirname, path), buffer, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
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
  auto cameras = vector<scene_camera>{};
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
          camera.film   = 0.036;
        } else if (type == "perspective") {
          auto& gpersp  = gcamera.at("perspective");
          camera.aspect = gpersp.value("aspectRatio", 0.0f);
          auto yfov     = gpersp.value("yfov", radians(45.0f));
          if (camera.aspect == 0) camera.aspect = 16.0f / 9.0f;
          camera.film = 0.036;
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
  auto get_texture = [&gltf](const njson& js, const string& name) -> int {
    if (!js.contains(name)) return invalidid;
    auto& ginfo    = js.at(name);
    auto& gtexture = gltf.at("textures").at(ginfo.value("index", -1));
    return gtexture.value("source", -1);
  };

  // convert textures
  auto texture_paths = vector<string>{};
  if (gltf.contains("images")) {
    try {
      for (auto& gimage : gltf.at("images")) {
        scene.textures.emplace_back();
        texture_paths.push_back(gimage.value("uri", ""));
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
        material.type     = scene_material_type::gltfpbr;
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
            material.type      = scene_material_type::transparent;
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
              return primitive_error();
            } else {
              return primitive_error();
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
              return primitive_error();
            } else {
              return primitive_error();
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
    // mutex
    auto mutex = std::mutex{};
    // load textures
    parallel_foreach(scene.textures, [&](auto& texture) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
      }
      auto& path = texture_paths[&texture - &scene.textures.front()];
      auto  err  = string{};
      if (!load_texture(path_join(dirname, path), texture, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
  }

  // fix scene
  add_missing_material(scene);
  add_missing_camera(scene);
  add_missing_radius(scene);

  // load done
  return true;
}

// Load a scene
static bool save_gltf_scene(const string& filename, const scene_model& scene,
    string& error, bool noparallel) {
  auto write_error = [filename, &error]() {
    error = filename + ": write error";
    return false;
  };
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };
  auto fvshape_error = [filename, &error]() {
    error = filename + ": face-varying not supported";
    return false;
  };

  // convert scene to json
  auto gltf = njson::object();

  // asset
  {
    auto& gasset        = gltf["asset"];
    gasset              = njson::object();
    gasset["version"]   = "2.0";
    gasset["generator"] = "Yocto/GL - https://github.com/xelatihy/yocto-gl";
    if (!scene.copyright.empty()) gasset["copyright"] = scene.copyright;
  }

  // cameras
  if (!scene.cameras.empty()) {
    auto& gcameras = gltf["cameras"];
    gcameras       = njson::array();
    for (auto& camera : scene.cameras) {
      auto& gcamera               = gcameras.emplace_back();
      gcamera                     = njson::object();
      gcamera["name"]             = get_camera_name(scene, camera);
      gcamera["type"]             = "perspective";
      auto& gperspective          = gcamera["perspective"];
      gperspective                = njson::object();
      gperspective["aspectRatio"] = camera.aspect;
      gperspective["yfov"]        = 0.660593;  // TODO(fabio): yfov
      gperspective["znear"]       = 0.001;     // TODO(fabio): configurable?
    }
  }

  // textures
  if (!scene.textures.empty()) {
    auto& gtextures  = gltf["textures"];
    gtextures        = njson::array();
    auto& gsamplers  = gltf["samplers"];
    gsamplers        = njson::array();
    auto& gimages    = gltf["images"];
    gimages          = njson::array();
    auto& gsampler   = gsamplers.emplace_back();
    gsampler         = njson::object();
    gsampler["name"] = "sampler";
    for (auto& texture : scene.textures) {
      auto  name          = get_texture_name(scene, texture);
      auto& gimage        = gimages.emplace_back();
      gimage              = njson::object();
      gimage["name"]      = name;
      gimage["uri"]       = "textures/" + name + ".png";
      auto& gtexture      = gtextures.emplace_back();
      gtexture            = njson::object();
      gtexture["name"]    = name;
      gtexture["sampler"] = 0;
      gtexture["source"]  = (int)gimages.size() - 1;
    }
  }

  // materials
  if (!scene.materials.empty()) {
    auto& gmaterials = gltf["materials"];
    gmaterials       = njson::array();
    for (auto& material : scene.materials) {
      auto& gmaterial             = gmaterials.emplace_back();
      gmaterial                   = njson::object();
      gmaterial["name"]           = get_material_name(scene, material);
      gmaterial["emissiveFactor"] = material.emission;
      auto& gpbr                  = gmaterial["pbrMetallicRoughness"];
      gpbr                        = njson::object();
      gpbr["baseColorFactor"]     = vec4f{material.color.x, material.color.y,
          material.color.z, material.opacity};
      gpbr["metallicFactor"]      = material.metallic;
      gpbr["roughnessFactor"]     = material.roughness;
      if (material.emission_tex != invalidid) {
        gmaterial["emissiveTexture"]          = njson::object();
        gmaterial["emissiveTexture"]["index"] = material.emission_tex;
      }
      if (material.normal_tex != invalidid) {
        gmaterial["normalTexture"]          = njson::object();
        gmaterial["normalTexture"]["index"] = material.normal_tex;
      }
      if (material.color_tex != invalidid) {
        gpbr["baseColorTexture"]          = njson::object();
        gpbr["baseColorTexture"]["index"] = material.color_tex;
      }
      if (material.roughness_tex != invalidid) {
        gpbr["metallicRoughnessTexture"]          = njson::object();
        gpbr["metallicRoughnessTexture"]["index"] = material.roughness_tex;
      }
    }
  }

  // add an accessor
  auto set_view = [](njson& gview, njson& gbuffer, const auto& data,
                      size_t buffer_id) {
    gview                 = njson::object();
    gview["buffer"]       = buffer_id;
    gview["byteLength"]   = data.size() * sizeof(data.front());
    gview["byteOffset"]   = gbuffer["byteLength"];
    gbuffer["byteLength"] = gbuffer.value("byteLength", (size_t)0) +
                            data.size() * sizeof(data.front());
  };
  auto set_vaccessor = [](njson& gaccessor, const auto& data, size_t view_id,
                           bool minmax = false) {
    static auto types = unordered_map<size_t, string>{
        {1, "SCALAR"}, {2, "VEC2"}, {3, "VEC3"}, {4, "VEC4"}};
    gaccessor                  = njson::object();
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
  auto set_iaccessor = [](njson& gaccessor, const auto& data, size_t view_id,
                           bool minmax = false) {
    gaccessor                  = njson::object();
    gaccessor["bufferView"]    = view_id;
    gaccessor["componentType"] = 5125;
    gaccessor["count"] = data.size() * sizeof(data.front()) / sizeof(int);
    gaccessor["type"]  = "SCALAR";
  };

  // meshes
  auto shape_primitives = vector<njson>();
  shape_primitives.reserve(scene.shapes.size());
  if (!scene.shapes.empty()) {
    auto& gaccessors = gltf["accessors"];
    gaccessors       = njson::array();
    auto& gviews     = gltf["bufferViews"];
    gviews           = njson::array();
    auto& gbuffers   = gltf["buffers"];
    gbuffers         = njson::array();
    for (auto& shape : scene.shapes) {
      auto& gbuffer         = gbuffers.emplace_back();
      gbuffer["uri"]        = "shapes/" + get_shape_name(scene, shape) + ".bin";
      gbuffer["byteLength"] = (size_t)0;
      auto& gprimitive      = shape_primitives.emplace_back();
      gprimitive            = njson::object();
      auto& gattributes     = gprimitive["attributes"];
      gattributes           = njson::object();
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
    gmeshes       = njson::array();
    for (auto& instance : scene.instances) {
      auto key = mesh_key{instance.shape, instance.material};
      if (mesh_map.find(key) != mesh_map.end()) continue;
      auto& gmesh   = gmeshes.emplace_back();
      gmesh         = njson::object();
      gmesh["name"] = get_shape_name(scene, instance.shape) + "_" +
                      get_material_name(scene, instance.material);
      gmesh["primitives"] = njson::array();
      gmesh["primitives"].push_back(shape_primitives.at(instance.shape));
      gmesh["primitives"].back()["material"] = instance.material;
      mesh_map[key]                          = gmeshes.size() - 1;
    }
  } else if (!scene.shapes.empty()) {
    auto& gmeshes = gltf["meshes"];
    gmeshes       = njson::array();
    auto shape_id = 0;
    for (auto& primitives : shape_primitives) {
      auto& gmesh         = gmeshes.emplace_back();
      gmesh               = njson::object();
      gmesh["name"]       = get_shape_name(scene, shape_id++);
      gmesh["primitives"] = njson::array();
      gmesh["primitives"].push_back(primitives);
    }
  }

  // nodes
  if (!scene.cameras.empty() || !scene.instances.empty()) {
    auto& gnodes   = gltf["nodes"];
    gnodes         = njson::array();
    auto camera_id = 0;
    for (auto& camera : scene.cameras) {
      auto& gnode     = gnodes.emplace_back();
      gnode           = njson::object();
      gnode["name"]   = get_camera_name(scene, camera);
      gnode["matrix"] = frame_to_mat(camera.frame);
      gnode["camera"] = camera_id++;
    }
    for (auto& instance : scene.instances) {
      auto& gnode     = gnodes.emplace_back();
      gnode           = njson::object();
      gnode["name"]   = get_instance_name(scene, instance);
      gnode["matrix"] = frame_to_mat(instance.frame);
      gnode["mesh"]   = mesh_map.at({instance.shape, instance.material});
    }
    // root children
    auto& groot     = gnodes.emplace_back();
    groot           = njson::object();
    groot["name"]   = "root";
    auto& gchildren = groot["children"];
    gchildren       = njson::array();
    for (auto idx = (size_t)0; idx < gnodes.size() - 1; idx++)
      gchildren.push_back(idx);
    // scene
    auto& gscenes     = gltf["scenes"];
    gscenes           = njson::array();
    auto& gscene      = gscenes.emplace_back();
    gscene            = njson::object();
    auto& gscenenodes = gscene["nodes"];
    gscenenodes       = njson::array();
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
  auto dirname = path_dirname(filename);

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
    // mutex
    auto mutex = std::mutex{};
    // save shapes
    parallel_foreach(scene.shapes, [&](auto& shape) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
      }
      auto path = "shapes/" + get_shape_name(scene, shape) + ".bin";
      auto err  = string{};
      if (!save_binshape(path_join(dirname, path), shape, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
    // save textures
    parallel_foreach(scene.textures, [&](auto& texture) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
      }
      auto path = "textures/" + get_texture_name(scene, texture) +
                  (!texture.pixelsf.empty() ? ".hdr"s : ".png"s);
      auto err = string{};
      if (!save_texture(path_join(dirname, path), texture, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
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
static bool load_pbrt_scene(const string& filename, scene_model& scene,
    string& error, bool noparallel) {
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // load pbrt
  auto pbrt = pbrt_model{};
  if (!load_pbrt(filename, pbrt, error)) return false;

  // convert cameras
  for (auto& pcamera : pbrt.cameras) {
    auto& camera  = scene.cameras.emplace_back();
    camera.frame  = pcamera.frame;
    camera.aspect = pcamera.aspect;
    camera.film   = 0.036;
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
  auto material_type_map = unordered_map<pbrt_mtype, scene_material_type>{
      {pbrt_mtype::matte, scene_material_type::matte},
      {pbrt_mtype::plastic, scene_material_type::glossy},
      {pbrt_mtype::metal, scene_material_type::metallic},
      {pbrt_mtype::glass, scene_material_type::refractive},
      {pbrt_mtype::thinglass, scene_material_type::transparent},
      {pbrt_mtype::subsurface, scene_material_type::matte},
  };

  // convert material
  for (auto& pmaterial : pbrt.materials) {
    auto& material = scene.materials.emplace_back();
    material.type  = material_type_map.at(pmaterial.type);
    if (pmaterial.emission != zero3f) {
      material.type = scene_material_type::matte;
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
    // mutex
    auto mutex = std::mutex{};
    // load shapes
    parallel_foreach(scene.shapes, [&](auto& shape) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
      }
      auto& path = shapes_paths[&shape - &scene.shapes.front()];
      if (path.empty()) return;
      auto err = string{};
      if (!load_shape(path_join(dirname, path), shape, err, true)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
    // load textures
    parallel_foreach(scene.textures, [&](auto& texture) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
      }
      auto& path = texture_paths[&texture - &scene.textures.front()];
      auto  err  = string{};
      if (!load_texture(path_join(dirname, path), texture, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
  }

  // fix scene
  add_missing_camera(scene);
  add_missing_radius(scene);

  // done
  return true;
}

// Save a pbrt scene
static bool save_pbrt_scene(const string& filename, const scene_model& scene,
    string& error, bool noparallel) {
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

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
  auto material_type_map = unordered_map<scene_material_type, pbrt_mtype>{
      {scene_material_type::matte, pbrt_mtype::matte},
      {scene_material_type::glossy, pbrt_mtype::plastic},
      {scene_material_type::metallic, pbrt_mtype::metal},
      {scene_material_type::refractive, pbrt_mtype::glass},
      {scene_material_type::transparent, pbrt_mtype::thinglass},
      {scene_material_type::subsurface, pbrt_mtype::matte},
      {scene_material_type::volume, pbrt_mtype::matte},
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
  auto dirname = path_dirname(filename);

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
    // mutex
    auto mutex = std::mutex{};
    // save shapes
    parallel_foreach(scene.shapes, [&](auto& shape) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
      }
      auto path = "shapes/" + get_shape_name(scene, shape) + ".ply";
      auto err  = string{};
      if (!save_shape(path_join(dirname, path), shape, err, true)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
    // save textures
    parallel_foreach(scene.textures, [&](auto& texture) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
      }
      auto path = "textures/" + get_texture_name(scene, texture) +
                  (!texture.pixelsf.empty() ? ".hdr"s : ".png"s);
      auto err = string{};
      if (!save_texture(path_join(dirname, path), texture, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
  }

  // done
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
