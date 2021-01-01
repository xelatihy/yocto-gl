//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//

#include <yocto/yocto_commonio.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>
#if YOCTO_OPENGL == 1
#include <yocto_gui/yocto_imageviewer.h>
#endif
using namespace yocto;

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

bool make_image_preset(const string& type_, image<vec4f>& img, string& error) {
  auto set_region = [](image<vec4f>& img, const image<vec4f>& region,
                        const vec2i& offset) {
    for (auto j = 0; j < region.height(); j++) {
      for (auto i = 0; i < region.width(); i++) {
        if (!img.contains({i, j})) continue;
        img[vec2i{i, j} + offset] = region[{i, j}];
      }
    }
  };

  auto type = path_basename(type_);

  auto size = vec2i{1024, 1024};
  if (type.find("sky") != type.npos) size = {2048, 1024};
  if (type.find("images2") != type.npos) size = {2048, 1024};
  if (type == "grid") {
    img = make_grid(size);
  } else if (type == "checker") {
    img = make_checker(size);
  } else if (type == "bumps") {
    img = make_bumps(size);
  } else if (type == "uvramp") {
    img = make_uvramp(size);
  } else if (type == "gammaramp") {
    img = make_gammaramp(size);
  } else if (type == "blackbodyramp") {
    img = make_blackbodyramp(size);
  } else if (type == "uvgrid") {
    img = make_uvgrid(size);
  } else if (type == "colormap") {
    img = make_colormapramp(size);
    img = srgb_to_rgb(img);
  } else if (type == "sky") {
    img = make_sunsky(
        size, pif / 4, 3.0, false, 1.0, 1.0, vec3f{0.7, 0.7, 0.7});
  } else if (type == "sunsky") {
    img = make_sunsky(size, pif / 4, 3.0, true, 1.0, 1.0, vec3f{0.7, 0.7, 0.7});
  } else if (type == "noise") {
    img = make_noisemap(size, 1);
  } else if (type == "fbm") {
    img = make_fbmmap(size, 1);
  } else if (type == "ridge") {
    img = make_ridgemap(size, 1);
  } else if (type == "turbulence") {
    img = make_turbulencemap(size, 1);
  } else if (type == "bump-normal") {
    img = make_bumps(size);
    img = srgb_to_rgb(bump_to_normal(img, 0.05f));
  } else if (type == "images1") {
    auto sub_types = vector<string>{"grid", "uvgrid", "checker", "gammaramp",
        "bumps", "bump-normal", "noise", "fbm", "blackbodyramp"};
    auto sub_imgs  = vector<image<vec4f>>(sub_types.size());
    for (auto i = 0; i < sub_imgs.size(); i++) {
      if (!make_image_preset(sub_types[i], sub_imgs[i], error)) return false;
    }
    auto montage_size = zero2i;
    for (auto& sub_img : sub_imgs) {
      montage_size.x += sub_img.width();
      montage_size.y = max(montage_size.y, sub_img.height());
    }
    img      = image<vec4f>(montage_size);
    auto pos = 0;
    for (auto& sub_img : sub_imgs) {
      set_region(img, sub_img, {pos, 0});
      pos += sub_img.width();
    }
  } else if (type == "images2") {
    auto sub_types = vector<string>{"sky", "sunsky"};
    auto sub_imgs  = vector<image<vec4f>>(sub_types.size());
    for (auto i = 0; i < sub_imgs.size(); i++) {
      if (!make_image_preset(sub_types[i], sub_imgs[i], error)) return false;
    }
    auto montage_size = zero2i;
    for (auto& sub_img : sub_imgs) {
      montage_size.x += sub_img.width();
      montage_size.y = max(montage_size.y, sub_img.height());
    }
    img      = image<vec4f>(montage_size);
    auto pos = 0;
    for (auto& sub_img : sub_imgs) {
      set_region(img, sub_img, {pos, 0});
      pos += sub_img.width();
    }
  } else if (type == "test-floor") {
    img = make_grid(size);
    img = add_border(img, 0.0025);
  } else if (type == "test-grid") {
    img = make_grid(size);
  } else if (type == "test-checker") {
    img = make_checker(size);
  } else if (type == "test-bumps") {
    img = make_bumps(size);
  } else if (type == "test-uvramp") {
    img = make_uvramp(size);
  } else if (type == "test-gammaramp") {
    img = make_gammaramp(size);
  } else if (type == "test-blackbodyramp") {
    img = make_blackbodyramp(size);
  } else if (type == "test-colormapramp") {
    img = make_colormapramp(size);
    img = srgb_to_rgb(img);
  } else if (type == "test-uvgrid") {
    img = make_uvgrid(size);
  } else if (type == "test-sky") {
    img = make_sunsky(
        size, pif / 4, 3.0, false, 1.0, 1.0, vec3f{0.7, 0.7, 0.7});
  } else if (type == "test-sunsky") {
    img = make_sunsky(size, pif / 4, 3.0, true, 1.0, 1.0, vec3f{0.7, 0.7, 0.7});
  } else if (type == "test-noise") {
    img = make_noisemap(size);
  } else if (type == "test-fbm") {
    img = make_noisemap(size);
  } else if (type == "test-bumps-normal") {
    img = make_bumps(size);
    img = bump_to_normal(img, 0.05);
  } else if (type == "test-bumps-displacement") {
    img = make_bumps(size);
    img = srgb_to_rgb(img);
  } else if (type == "test-fbm-displacement") {
    img = make_fbmmap(size);
    img = srgb_to_rgb(img);
  } else if (type == "test-checker-opacity") {
    img = make_checker(size, 1, {1, 1, 1, 1}, {0, 0, 0, 0});
  } else if (type == "test-grid-opacity") {
    img = make_grid(size, 1, {1, 1, 1, 1}, {0, 0, 0, 0});
  } else {
    error = "unknown preset";
    img   = {};
    return false;
  }
  return true;
}

bool make_image_preset(const string& type, image<vec4b>& img, string& error) {
  auto imgf = image<vec4f>{};
  if (make_image_preset(type, imgf, error)) return false;
  img = rgb_to_srgbb(imgf);
  return true;
}

bool is_preset_filename(const string& filename) {
  return path_extension(filename) == ".ypreset";
}

bool is_preset_hdr(const string& type_) {
  auto type = path_basename(type_);
  return type.find("sky") != string::npos && type.find("sun") != string::npos;
}

}  // namespace yocto

// convert params
struct convert_params {
  string image    = "image.png";
  string output   = "out.png";
  bool   logo     = false;
  float  exposure = 0;
  bool   filmic   = false;
  int    width    = 0;
  int    height   = 0;
};

// Json IO
void serialize_value(json_mode mode, json_value& json, convert_params& value,
    const string& description) {
  serialize_object(mode, json, value, description);
  serialize_property(mode, json, value.image, "image", "Input image.", true);
  serialize_property(mode, json, value.output, "output", "Output image.");
  serialize_property(
      mode, json, value.exposure, "exposure", "Tonemap exposure.");
  serialize_property(mode, json, value.filmic, "filmic", "Tonemap filmic.");
  serialize_property(mode, json, value.width, "width", "Resize width.");
  serialize_property(mode, json, value.height, "height", "Resize height.");
  serialize_property(mode, json, value.logo, "logo", "Add logo.");
  serialize_clipositionals(mode, json, {"image"});
  serialize_clialternates(mode, json,
      {{"output", "o"}, {"exposure", "e"}, {"filmic", "f"}, {"width", "w"},
          {"height", "h"}, {"logo", "L"}});
}

// convert images
int run_convert(const convert_params& params) {
  // load
  auto hdr     = image<vec4f>{};
  auto ldr     = image<vec4b>{};
  auto ioerror = string{};
  if (is_preset_filename(params.image)) {
    if (is_preset_hdr(params.image)) {
      if (!make_image_preset(path_basename(params.image), hdr, ioerror))
        return print_fatal(ioerror);
    } else {
      if (!make_image_preset(path_basename(params.image), ldr, ioerror))
        return print_fatal(ioerror);
    }
  } else if (is_hdr_filename(params.image)) {
    if (!load_image(params.image, hdr, ioerror)) return print_fatal(ioerror);
  } else {
    if (!load_image(params.image, ldr, ioerror)) return print_fatal(ioerror);
  }

  // resize if needed
  if (params.width != 0 || params.height != 0) {
    if (!hdr.empty()) {
      hdr = resize_image(hdr, params.width, params.height);
    } else {
      ldr = resize_image(ldr, params.width, params.height);
    }
  }

  // tonemap if needed
  if (!hdr.empty() && !is_hdr_filename(params.output)) {
    ldr = tonemap_imageb(hdr, params.exposure, params.filmic);
    hdr = {};
  }
  if (!ldr.empty() && is_hdr_filename(params.output)) {
    hdr = srgb_to_rgb(ldr);
    ldr = {};
  }

  // apply logo
  if (params.logo) {
    if (!hdr.empty()) {
      hdr = add_logo(hdr);
    } else {
      ldr = add_logo(ldr);
    }
  }

  // save
  if (!hdr.empty()) {
    if (!save_image(params.output, hdr, ioerror)) return print_fatal(ioerror);
  } else {
    if (!save_image(params.output, ldr, ioerror)) return print_fatal(ioerror);
  }

  // done
  return 0;
}

// view params
struct view_params {
  vector<string> images = {"image.png"};
  string         output = "out.png";
  bool           logo   = false;
};

// Json IO
void serialize_value(json_mode mode, json_value& json, view_params& value,
    const string& description) {
  serialize_object(mode, json, value, description);
  serialize_property(mode, json, value.images, "images", "Input images.", true);
  serialize_property(mode, json, value.output, "output", "Output image.");
  serialize_clipositionals(mode, json, {"images"});
  serialize_clialternates(mode, json, {{"output", "o"}});
}

#ifndef YOCTO_OPENGL

// view images
int run_view(const view_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

// view images
int run_view(const view_params& params) {
  // open viewer
  auto viewer_guard = make_imageview("yimage");
  auto viewer       = viewer_guard.get();

  // set image
  for (auto& filename : params.images) {
    // load
    auto hdr     = image<vec4f>{};
    auto ldr     = image<vec4b>{};
    auto ioerror = string{};
    if (is_preset_filename(filename)) {
      if (is_preset_hdr(filename)) {
        if (!make_image_preset(path_basename(filename), hdr, ioerror))
          return print_fatal(ioerror);
      } else {
        if (!make_image_preset(path_basename(filename), ldr, ioerror))
          return print_fatal(ioerror);
      }
    } else if (is_hdr_filename(filename)) {
      if (!load_image(filename, hdr, ioerror)) return print_fatal(ioerror);
    } else {
      if (!load_image(filename, ldr, ioerror)) return print_fatal(ioerror);
    }

    // push image to the viewer
    if (!hdr.empty()) {
      set_image(viewer, filename, hdr);
    } else {
      set_image(viewer, filename, ldr);
    }
  }

  // run view
  run_view(viewer);

  // done
  return 0;
}

#endif

// grade params
struct grade_params : colorgrade_params {
  string image  = "image.png";
  string output = "out.png";
  bool   logo   = false;
};

// Json IO
void serialize_value(json_mode mode, json_value& json, colorgrade_params& value,
    const string& description) {
  serialize_object(mode, json, value, description);
  serialize_property(mode, json, value.exposure, "exposure", "Hdr exposure");
  serialize_property(
      mode, json, (array<float, 3>&)value.tint, "tint", "Hdr tint");
  serialize_property(
      mode, json, value.lincontrast, "lincontrast", "Hdr lin contrast");
  serialize_property(
      mode, json, value.logcontrast, "logcontrast", "Hdr log contrast");
  serialize_property(
      mode, json, value.linsaturation, "linsaturation", "Hdr saturation");
  serialize_property(mode, json, value.filmic, "filmic", "Hdr filmic curve");
  serialize_property(mode, json, value.srgb, "srgb", "sRGB coversion");
  serialize_property(mode, json, value.contrast, "contrast", "Ldr contrast");
  serialize_property(
      mode, json, value.saturation, "saturation", "Ldr saturation");
  serialize_property(mode, json, value.shadows, "shadows", "Ldr shadows");
  serialize_property(mode, json, value.midtones, "midtones", "Ldr midtones");
  serialize_property(
      mode, json, value.highlights, "highlights", "Ldr highlights");
  serialize_property(mode, json, (array<float, 3>&)value.shadows_color,
      "shadows_color", "Ldr shadows color");
  serialize_property(mode, json, (array<float, 3>&)value.midtones_color,
      "midtones_color", "Ldr sidtones color");
  serialize_property(mode, json, (array<float, 3>&)value.highlights_color,
      "highlights_color", "Ldr sighlights color");
}

// Json IO
void serialize_value(json_mode mode, json_value& json, grade_params& value,
    const string& description) {
  serialize_object(mode, json, value, description);
  serialize_property(mode, json, value.image, "image", "Input image.", true);
  serialize_property(mode, json, value.output, "output", "Output image.");
  serialize_value(mode, json, (colorgrade_params&)value, description);
  serialize_clipositionals(mode, json, {"image"});
  serialize_clialternates(mode, json, {{"output", "o"}});
}

#ifndef YOCTO_OPENGL

// grade images
int run_grade(const grade_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

// grade images
int run_grade(const grade_params& params) {
  // open viewer
  auto viewer_guard = make_imageview("yimage");
  auto viewer       = viewer_guard.get();

  // load image
  auto img     = image<vec4f>{};
  auto ioerror = string{};
  if (is_preset_filename(params.image)) {
    if (!make_image_preset(path_basename(params.image), img, ioerror))
      return print_fatal(ioerror);
  } else {
    if (!load_image(params.image, img, ioerror)) return print_fatal(ioerror);
  }

  // grade image
  auto graded = image<vec4b>{img.imsize()};
  colorgrade_image_mt(graded, img, true, params);

  // set view
  set_image(viewer, params.image, graded);
  set_params(
      viewer, params.image, to_json(params), to_schema(params, "Color grade"));

  // set callback
  set_callback(viewer, [&params, &graded, &img, viewer](const string& name,
                           const json_value& uiparams, const gui_input&) {
    if (uiparams.is_null()) return;
    serialize_value(json_mode::from_json, (json_value&)uiparams,
        (colorgrade_params&)params, "");
    colorgrade_image_mt(graded, img, true, params);
    set_image(viewer, name, graded);
  });

  // run view
  run_view(viewer);

  // done
  return 0;
}

#endif

// resize params
struct diff_params {
  string image1    = "image1.png";
  string image2    = "image2.png";
  string output    = "";
  bool   logo      = false;
  bool   signal    = false;
  float  threshold = 0;
};

// Json IO
void serialize_value(json_mode mode, json_value& json, diff_params& value,
    const string& description) {
  serialize_object(mode, json, value, description);
  serialize_property(
      mode, json, value.image1, "image1", "Input image 1.", true);
  serialize_property(
      mode, json, value.image2, "image2", "Input image 2.", true);
  serialize_property(mode, json, value.output, "output", "Output image.");
  serialize_property(mode, json, value.signal, "signal", "Error on diff.");
  serialize_property(
      mode, json, value.threshold, "threshold", "Diff threshold.");
  serialize_property(mode, json, value.logo, "logo", "Add logo.");
  serialize_clipositionals(mode, json, {"image1", "image2"});
  serialize_clialternates(mode, json, {{"output", "o"}});
}

// resize images
int run_diff(const diff_params& params) {
  // load
  auto ioerror = string{};
  auto img1 = image<vec4f>{}, img2 = image<vec4f>{};
  if (path_extension(params.image1) == ".ypreset") {
    if (!make_image_preset(path_basename(params.image1), img1, ioerror))
      return false;
  } else {
    if (!load_image(params.image1, img1, ioerror)) return false;
  }
  if (path_extension(params.image2) == ".ypreset") {
    if (!make_image_preset(path_basename(params.image2), img2, ioerror))
      return false;
  } else {
    if (!load_image(params.image2, img2, ioerror)) return false;
  }

  // check sizes
  if (img1.imsize() != img2.imsize()) {
    ioerror = "image sizes are different";
    return print_fatal(ioerror);
  }

  // compute diff
  auto diff = image_difference(img1, img2, true);

  // save
  if (params.output != "") {
    if (!save_image(
            params.output, params.logo ? add_logo(diff) : diff, ioerror))
      return print_fatal(ioerror);
  }

  // check diff
  if (params.signal) {
    for (auto& c : diff) {
      if (max(xyz(c)) > params.threshold) {
        ioerror = "image content differs";
        return print_fatal(ioerror);
      }
    }
  }

  // done
  return 0;
}

// setalpha params
struct setalpha_params {
  string image      = "image.png";
  string alpha      = "alpha.png";
  string output     = "out.png";
  bool   logo       = false;
  bool   from_color = false;
  bool   to_color   = false;
};

// Json IO
void serialize_value(json_mode mode, json_value& json, setalpha_params& value,
    const string& description) {
  serialize_object(mode, json, value, description);
  serialize_property(mode, json, value.image, "image", "Input image.", true);
  serialize_property(mode, json, value.alpha, "alpha", "Alpha image.", true);
  serialize_property(mode, json, value.output, "output", "Output image.");
  serialize_property(
      mode, json, value.from_color, "from-color", "Alpha from color.");
  serialize_property(
      mode, json, value.to_color, "to-color", "Color from alpha.");
  serialize_property(mode, json, value.logo, "logo", "Add logo.");
  serialize_clipositionals(mode, json, {"image", "alpha"});
  serialize_clialternates(mode, json, {{"output", "o"}});
}

// setalpha images
int run_setalpha(const setalpha_params& params) {
  // load
  auto ioerror = string{};
  auto img = image<vec4f>{}, alpha = image<vec4f>{};
  if (path_extension(params.image) == ".ypreset") {
    if (!make_image_preset(path_basename(params.image), img, ioerror))
      return print_fatal(ioerror);
  } else {
    if (!load_image(params.image, img, ioerror)) return print_fatal(ioerror);
  }
  if (path_extension(params.alpha) == ".ypreset") {
    if (!make_image_preset(path_basename(params.alpha), alpha, ioerror))
      return print_fatal(ioerror);
  } else {
    if (!load_image(params.alpha, alpha, ioerror)) return print_fatal(ioerror);
  }

  // set alpha
  if (params.from_color) {
    if (img.imsize() != alpha.imsize()) print_fatal("bad image size");
    for (auto j = 0; j < img.height(); j++)
      for (auto i = 0; i < img.width(); i++)
        img[{i, j}].w = mean(xyz(alpha[{i, j}]));
  } else {
    if (img.imsize() != alpha.imsize()) print_fatal("bad image size");
    for (auto j = 0; j < img.height(); j++)
      for (auto i = 0; i < img.width(); i++) img[{i, j}].w = alpha[{i, j}].w;
  }

  // set color from alpha
  if (params.to_color) {
    for (auto& c : img) c = vec4f{c.w, c.w, c.w, c.w};
  }

  // edit alpha
  auto out = image<vec4f>{img.imsize()};
  for (auto idx = (size_t)0; idx < out.count(); idx++) {
    auto a = params.from_color ? mean(xyz(img[idx])) : img[idx].w;
    if (params.to_color) {
      out[idx] = {a, a, a, a};
    } else {
      out[idx].w = a;
    }
  }

  // save
  if (!save_image(params.output, params.logo ? add_logo(out) : out, ioerror))
    return print_fatal(ioerror);

  // done
  return 0;
}

struct app_params {
  string          command  = "convert";
  convert_params  convert  = {};
  view_params     view     = {};
  grade_params    grade    = {};
  diff_params     diff     = {};
  setalpha_params setalpha = {};
};

// Json IO
void serialize_value(json_mode mode, json_value& json, app_params& value,
    const string& description) {
  serialize_object(mode, json, value, description);
  serialize_command(mode, json, value.command, "command", "Command.");
  serialize_property(mode, json, value.convert, "convert", "Convert images.");
  serialize_property(mode, json, value.view, "view", "View images.");
  serialize_property(mode, json, value.grade, "grade", "Grade images.");
  serialize_property(mode, json, value.diff, "diff", "Diff two images.");
  serialize_property(
      mode, json, value.setalpha, "setalpha", "Set alpha in images.");
}

int main(int argc, const char* argv[]) {
  // command line parameters
  auto params = app_params{};
  parse_cli(params, "Process and view images", argc, argv);

  // dispatch commands
  if (params.command == "convert") {
    return run_convert(params.convert);
  } else if (params.command == "view") {
    return run_view(params.view);
  } else if (params.command == "grade") {
    return run_grade(params.grade);
  } else if (params.command == "diff") {
    return run_diff(params.diff);
  } else if (params.command == "setalpha") {
    return run_setalpha(params.setalpha);
  } else {
    return print_fatal("unknown command " + params.command);
  }
}
