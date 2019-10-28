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

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_imageio.h"
#include "yocto_commonio.h"
#include "yocto_image.h"

#if !defined(_WIN32) && !defined(_WIN64)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wunused-variable"
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#endif

// #ifndef _clang_analyzer__

#define STB_IMAGE_IMPLEMENTATION
#include "ext/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "ext/stb_image_write.h"

#define TINYEXR_IMPLEMENTATION
#include "ext/tinyexr.h"

// #endif

#if !defined(_WIN32) && !defined(_WIN64)
#pragma GCC diagnostic pop
#endif

#include <memory>

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGEIO
// -----------------------------------------------------------------------------
namespace yocto {

// Split a string
static inline vector<string> split_string(const string& str) {
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
static inline float* load_pfm(
    const char* filename, int* w, int* h, int* nc, int req) {
  auto fs = fopen(filename, "rb");
  if (!fs) return nullptr;
  auto fs_guard = std::unique_ptr<FILE, void (*)(FILE*)>{
      fs, [](FILE* f) { fclose(f); }};

  // buffer
  char buffer[4096];
  auto toks = vector<string>();

  // read magic
  if (!fgets(buffer, sizeof(buffer), fs)) return nullptr;
  toks = split_string(buffer);
  if (toks[0] == "Pf")
    *nc = 1;
  else if (toks[0] == "PF")
    *nc = 3;
  else
    return nullptr;

  // read w, h
  if (!fgets(buffer, sizeof(buffer), fs)) return nullptr;
  toks = split_string(buffer);
  *w   = atoi(toks[0].c_str());
  *h   = atoi(toks[1].c_str());

  // read scale
  if (!fgets(buffer, sizeof(buffer), fs)) return nullptr;
  toks   = split_string(buffer);
  auto s = atof(toks[0].c_str());

  // read the data (flip y)
  auto npixels = (size_t)(*w) * (size_t)(*h);
  auto nvalues = npixels * (size_t)(*nc);
  auto nrow    = (size_t)(*w) * (size_t)(*nc);
  auto pixels  = std::unique_ptr<float[]>(new float[nvalues]);
  for (auto j = *h - 1; j >= 0; j--) {
    if (fread(pixels.get() + j * nrow, sizeof(float), nrow, fs) != nrow)
      return nullptr;
  }

  // endian conversion
  if (s > 0) {
    for (auto i = 0; i < nvalues; ++i) {
      auto dta = (uint8_t*)(pixels.get() + i);
      std::swap(dta[0], dta[3]);
      std::swap(dta[1], dta[2]);
    }
  }

  // scale
  auto scl = (s > 0) ? s : -s;
  if (scl != 1) {
    for (auto i = 0; i < nvalues; i++) pixels[i] *= scl;
  }

  // proper number of channels
  if (!req || *nc == req) return pixels.release();

  // pack into channels
  if (req < 0 || req > 4) {
    return nullptr;
  }
  auto cpixels = std::unique_ptr<float[]>(new float[req * npixels]);
  for (auto i = 0ull; i < npixels; i++) {
    auto vp = pixels.get() + i * (*nc);
    auto cp = cpixels.get() + i * req;
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
  return cpixels.release();
}

// save pfm
static inline bool save_pfm(
    const char* filename, int w, int h, int nc, const float* pixels) {
  auto fs = fopen(filename, "wb");
  if (!fs) return false;
  auto fs_guard = std::unique_ptr<FILE, void (*)(FILE*)>{
      fs, [](FILE* f) { fclose(f); }};

  if (fprintf(fs, "%s\n", (nc == 1) ? "Pf" : "PF") < 0) return false;
  if (fprintf(fs, "%d %d\n", w, h) < 0) return false;
  if (fprintf(fs, "-1\n") < 0) return false;
  if (nc == 1 || nc == 3) {
    if (fwrite(pixels, sizeof(float), w * h * nc, fs) != w * h * nc)
      return false;
  } else {
    for (auto i = 0; i < w * h; i++) {
      auto vz = 0.0f;
      auto v  = pixels + i * nc;
      if (fwrite(v + 0, sizeof(float), 1, fs) != 1) return false;
      if (fwrite(v + 1, sizeof(float), 1, fs) != 1) return false;
      if (nc == 2) {
        if (fwrite(&vz, sizeof(float), 1, fs) != 1) return false;
      } else {
        if (fwrite(v + 2, sizeof(float), 1, fs) != 1) return false;
      }
    }
  }

  return true;
}

// load exr image weith tiny exr
static inline const char* get_tinyexr_error(int error) {
  switch (error) {
    case TINYEXR_ERROR_INVALID_MAGIC_NUMBER: return "INVALID_MAGIC_NUMBER";
    case TINYEXR_ERROR_INVALID_EXR_VERSION: return "INVALID_EXR_VERSION";
    case TINYEXR_ERROR_INVALID_ARGUMENT: return "INVALID_ARGUMENT";
    case TINYEXR_ERROR_INVALID_DATA: return "INVALID_DATA";
    case TINYEXR_ERROR_INVALID_FILE: return "INVALID_FILE";
    // case TINYEXR_ERROR_INVALID_PARAMETER: return "INVALID_PARAMETER";
    case TINYEXR_ERROR_CANT_OPEN_FILE: return "CANT_OPEN_FILE";
    case TINYEXR_ERROR_UNSUPPORTED_FORMAT: return "UNSUPPORTED_FORMAT";
    case TINYEXR_ERROR_INVALID_HEADER: return "INVALID_HEADER";
    default: throw std::runtime_error("unknown tinyexr error");
  }
}

// Check if an image is HDR based on filename.
bool is_hdr_filename(const string& filename) {
  auto ext = get_extension(filename);
  return ext == ".hdr" || ext == ".exr" || ext == ".pfm";
}

// Loads an hdr image.
image<vec4f> load_image(const string& filename) {
  auto img = image<vec4f>{};
  load_image(filename, img);
  return img;
}

// Loads an hdr image.
void load_image(const string& filename, image<vec4f>& img) {
  auto ext = get_extension(filename);
  if (ext == ".ypreset") {
    return make_image_preset(img, get_basename(filename));
  }
  if (ext == ".exr" || ext == ".EXR") {
    auto width = 0, height = 0;
    auto pixels = (float*)nullptr;
    if (auto error = LoadEXR(
            &pixels, &width, &height, filename.c_str(), nullptr);
        error < 0)
      throw std::runtime_error("error loading image " + filename + "("s +
                               get_tinyexr_error(error) + ")"s);
    if (!pixels) throw std::runtime_error("error loading image " + filename);
    img = image{{width, height}, (const vec4f*)pixels};
    free(pixels);
  } else if (ext == ".pfm" || ext == ".PFM") {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = load_pfm(filename.c_str(), &width, &height, &ncomp, 4);
    if (!pixels) throw std::runtime_error("error loading image " + filename);
    img = image{{width, height}, (const vec4f*)pixels};
    delete[] pixels;
  } else if (ext == ".hdr" || ext == ".HDR") {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_loadf(filename.c_str(), &width, &height, &ncomp, 4);
    if (!pixels) throw std::runtime_error("error loading image " + filename);
    img = image{{width, height}, (const vec4f*)pixels};
    free(pixels);
  } else if (!is_hdr_filename(filename)) {
    img = srgb_to_rgb(load_imageb(filename));
  } else {
    throw std::runtime_error("unsupported image format " + ext);
  }
}

// Saves an hdr image.
void save_image(const string& filename, const image<vec4f>& img) {
  auto ext = get_extension(filename);
  if (ext == ".hdr" || ext == ".HDR") {
    if (!stbi_write_hdr(filename.c_str(), img.size().x, img.size().y, 4,
            (float*)img.data()))
      throw std::runtime_error("error saving image " + filename);
  } else if (ext == ".pfm" || ext == ".PFM") {
    if (!save_pfm(filename.c_str(), img.size().x, img.size().y, 4,
            (float*)img.data()))
      throw std::runtime_error("error saving image " + filename);
  } else if (ext == ".exr" || ext == ".EXR") {
    if (SaveEXR((float*)img.data(), img.size().x, img.size().y, 4,
            filename.c_str()) < 0)
      throw std::runtime_error("error saving image " + filename);
  } else if (!is_hdr_filename(filename)) {
    save_imageb(filename, rgb_to_srgbb(img));
  } else {
    throw std::runtime_error("unsupported image format " + ext);
  }
}

// Loads an ldr image.
image<vec4b> load_imageb(const string& filename) {
  auto img = image<vec4b>{};
  load_imageb(filename, img);
  return img;
}

// Loads an ldr image.
void load_imageb(const string& filename, image<vec4b>& img) {
  auto ext = get_extension(filename);
  if (ext == ".ypreset") {
    return make_image_preset(img, get_basename(filename));
  }
  if (ext == ".png" || ext == ".PNG" || ext == ".jpg" || ext == ".JPG" ||
      ext == ".tga" || ext == ".TGA" || ext == ".bmp" || ext == ".BMP") {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_load(filename.c_str(), &width, &height, &ncomp, 4);
    if (!pixels) throw std::runtime_error("error loading image " + filename);
    img = image{{width, height}, (const vec4b*)pixels};
    free(pixels);
  } else if (is_hdr_filename(filename)) {
    img = rgb_to_srgbb(load_image(filename));
  } else {
    throw std::runtime_error("unsupported image format " + ext);
  }
}

// Saves an ldr image.
void save_imageb(const string& filename, const image<vec4b>& img) {
  auto ext = get_extension(filename);
  if (ext == ".png" || ext == ".PNG") {
    if (!stbi_write_png(filename.c_str(), img.size().x, img.size().y, 4,
            img.data(), img.size().x * 4))
      throw std::runtime_error("error saving image " + filename);
  } else if (ext == ".jpg" || ext == ".JPG") {
    if (!stbi_write_jpg(
            filename.c_str(), img.size().x, img.size().y, 4, img.data(), 75))
      throw std::runtime_error("error saving image " + filename);
  } else if (ext == ".tga" || ext == ".TGA") {
    if (!stbi_write_tga(
            filename.c_str(), img.size().x, img.size().y, 4, img.data()))
      throw std::runtime_error("error saving image " + filename);
  } else if (ext == ".bmp" || ext == ".BMP") {
    if (!stbi_write_bmp(
            filename.c_str(), img.size().x, img.size().y, 4, img.data()))
      throw std::runtime_error("error saving image " + filename);
  } else if (is_hdr_filename(filename)) {
    save_image(filename, srgb_to_rgb(img));
  } else {
    throw std::runtime_error("unsupported image format " + ext);
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR VOLUME IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

namespace impl {

// Volume load
static inline float* load_yvol(
    const char* filename, int* w, int* h, int* d, int* nc, int req) {
  auto fs = fopen(filename, "rb");
  if (!fs) return nullptr;
  auto fs_guard = std::unique_ptr<FILE, void (*)(FILE*)>{
      fs, [](FILE* f) { fclose(f); }};

  // buffer
  char buffer[4096];
  auto toks = vector<string>();

  // read magic
  if (!fgets(buffer, sizeof(buffer), fs)) return nullptr;
  toks = split_string(buffer);
  if (toks[0] != "YVOL") return nullptr;

  // read w, h
  if (!fgets(buffer, sizeof(buffer), fs)) return nullptr;
  toks = split_string(buffer);
  *w   = atoi(toks[0].c_str());
  *h   = atoi(toks[1].c_str());
  *d   = atoi(toks[2].c_str());
  *nc  = atoi(toks[3].c_str());

  // read data
  auto nvoxels = (size_t)(*w) * (size_t)(*h) * (size_t)(*d);
  auto nvalues = nvoxels * (size_t)(*nc);
  auto voxels  = std::unique_ptr<float[]>(new float[nvalues]);
  if (fread(voxels.get(), sizeof(float), nvalues, fs) != nvalues)
    return nullptr;

  // proper number of channels
  if (!req || *nc == req) return voxels.release();

  // pack into channels
  if (req < 0 || req > 4) {
    return nullptr;
  }
  auto cvoxels = std::unique_ptr<float[]>(new float[req * nvoxels]);
  for (auto i = 0; i < nvoxels; i++) {
    auto vp = voxels.get() + i * (*nc);
    auto cp = cvoxels.get() + i * req;
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
    } else if (*nc == 2) {
      switch (req) {
        case 1: cp[0] = vp[0]; break;
        case 2:
          cp[0] = vp[0];
          cp[1] = vp[1];
          break;
        case 3:
          cp[0] = vp[0];
          cp[1] = vp[1];
          break;
        case 4:
          cp[0] = vp[0];
          cp[1] = vp[1];
          break;
      }
    } else if (*nc == 3) {
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
    } else if (*nc == 4) {
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
          cp[3] = vp[3];
          break;
      }
    }
  }
  return cvoxels.release();
}

// save pfm
static inline bool save_yvol(
    const char* filename, int w, int h, int d, int nc, const float* voxels) {
  auto fs = fopen(filename, "wb");
  if (!fs) return false;
  auto fs_guard = std::unique_ptr<FILE, void (*)(FILE*)>{
      fs, [](FILE* f) { fclose(f); }};

  if (fprintf(fs, "YVOL\n") < 0) return false;
  if (fprintf(fs, "%d %d %d %d\n", w, h, d, nc) < 0) return false;
  auto nvalues = (size_t)w * (size_t)h * (size_t)d * (size_t)nc;
  if (fwrite(voxels, sizeof(float), nvalues, fs) != nvalues) return false;

  return true;
}

// Loads volume data from binary format.
void load_volume(const string& filename, volume<float>& vol) {
  auto width = 0, height = 0, depth = 0, ncomp = 0;
  auto voxels = load_yvol(filename.c_str(), &width, &height, &depth, &ncomp, 1);
  if (!voxels) {
    throw std::runtime_error("error loading volume " + filename);
  }
  vol = volume{{width, height, depth}, (const float*)voxels};
  delete[] voxels;
}

// Saves volume data in binary format.
void save_volume(const string& filename, const volume<float>& vol) {
  if (!save_yvol(filename.c_str(), vol.size().x, vol.size().y, vol.size().z, 1,
          vol.data())) {
    throw std::runtime_error("error saving volume " + filename);
  }
}

}  // namespace impl

// Loads volume data from binary format.
void load_volume(const string& filename, volume<float>& vol) {
  impl::load_volume(filename, vol);
}

// Saves volume data in binary format.
void save_volume(const string& filename, const volume<float>& vol) {
  impl::save_volume(filename, vol);
}

}  // namespace yocto
