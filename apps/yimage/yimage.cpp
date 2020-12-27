//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
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

bool make_image_preset(const string& type, image<vec4f>& img, string& error) {
  auto set_region = [](image<vec4f>& img, const image<vec4f>& region,
                        const vec2i& offset) {
    for (auto j = 0; j < region.height(); j++) {
      for (auto i = 0; i < region.width(); i++) {
        if (!img.contains({i, j})) continue;
        img[vec2i{i, j} + offset] = region[{i, j}];
      }
    }
  };

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

}  // namespace yocto

// convert params
struct convert_params {
  string image  = "";
  string output = "";
  bool   logo   = false;
};

// convert images
bool convert_images(const convert_params& params, string& error) {
  // load
  auto img = image<vec4f>{};
  if (path_extension(params.image) == ".ypreset") {
    if (!make_image_preset(path_basename(params.image), img, error))
      return false;
  } else {
    if (!load_image(params.image, img, error)) return false;
  }

  // save
  if (!save_image(params.output, params.logo ? add_logo(img) : img, error))
    return false;

  // done
  return true;
}

// tonemap params
struct tonemap_params {
  string image    = "";
  string output   = "";
  bool   logo     = false;
  float  exposure = 0;
  bool   filmic   = false;
  bool   srgb     = true;
};

// tonemap images
bool tonemap_images(const tonemap_params& params, string& error) {
  // load
  auto img = image<vec4f>{};
  if (path_extension(params.image) == ".ypreset") {
    if (!make_image_preset(path_basename(params.image), img, error))
      return false;
  } else {
    if (!load_image(params.image, img, error)) return false;
  }

  // tonemap
  auto out = tonemap_imageb(img, params.exposure, params.filmic, params.srgb);

  // save
  if (!save_image(params.output, params.logo ? add_logo(out) : out, error))
    return false;

  // done
  return true;
}

// resize params
struct resize_params {
  string image  = "";
  string output = "";
  bool   logo   = false;
  int    width  = 0;
  int    height = 0;
};

// resize images
bool resize_images(const resize_params& params, string& error) {
  // check for hdr values
  if (is_hdr_filename(params.image) != is_hdr_filename(params.output)) {
    error = "image and output should be both hdr or ldr";
    return false;
  }

  // check for hdr values
  if (params.width == 0 && params.height == 0) {
    error = "set width or height";
    return false;
  }

  // deterime whether is hdr
  if (is_hdr_filename(params.image)) {
    // load
    auto img = image<vec4f>{};
    if (path_extension(params.image) == ".ypreset") {
      if (!make_image_preset(path_basename(params.image), img, error))
        return false;
    } else {
      if (!load_image(params.image, img, error)) return false;
    }

    // resize
    auto out = resize_image(img, {params.width, params.height});

    // save
    if (!save_image(params.output, params.logo ? add_logo(out) : out, error))
      return false;

    // done
    return true;
  } else {
    // load
    auto img = image<vec4b>{};
    if (path_extension(params.image) == ".ypreset") {
      if (!make_image_preset(path_basename(params.image), img, error))
        return false;
    } else {
      if (!load_image(params.image, img, error)) return false;
    }

    // resize
    auto out = resize_image(img, {params.width, params.height});

    // save
    if (!save_image(params.output, params.logo ? add_logo(out) : out, error))
      return false;

    // done
    return true;
  }
}

// resize params
struct diff_params {
  string image1    = "";
  string image2    = "";
  string output    = "";
  bool   logo      = false;
  int    width     = 0;
  int    height    = 0;
  bool   signal    = false;
  float  threshold = 0;
};

// resize images
bool diff_images(const diff_params& params, string& error) {
  // load
  auto img1 = image<vec4f>{}, img2 = image<vec4f>{};
  if (path_extension(params.image1) == ".ypreset") {
    if (!make_image_preset(path_basename(params.image1), img1, error))
      return false;
  } else {
    if (!load_image(params.image1, img1, error)) return false;
  }
  if (path_extension(params.image2) == ".ypreset") {
    if (!make_image_preset(path_basename(params.image2), img2, error))
      return false;
  } else {
    if (!load_image(params.image2, img2, error)) return false;
  }

  // check sizes
  if (img1.imsize() != img2.imsize()) {
    error = "image sizes are different";
    return false;
  }

  // compute diff
  auto diff = image_difference(img1, img2, true);

  // save
  if (!save_image(params.output, params.logo ? add_logo(diff) : diff, error))
    return false;

  // check diff
  if (params.signal) {
    for (auto& c : diff) {
      if (max(xyz(c)) > params.threshold) {
        error = "image content differs";
        return false;
      }
    }
  }

  // done
  return true;
}

// setalpha params
struct setalpha_params {
  string image      = "";
  string alpha      = "";
  string output     = "";
  bool   logo       = false;
  bool   from_color = false;
  bool   to_color   = false;
};

// setalpha images
bool setalpha_images(const setalpha_params& params, string& error) {
  // load
  auto img = image<vec4f>{}, alpha = image<vec4f>{};
  if (path_extension(params.image) == ".ypreset") {
    if (!make_image_preset(path_basename(params.image), img, error))
      return false;
  } else {
    if (!load_image(params.image, img, error)) return false;
  }
  if (path_extension(params.alpha) == ".ypreset") {
    if (!make_image_preset(path_basename(params.alpha), alpha, error))
      return false;
  } else {
    if (!load_image(params.alpha, alpha, error)) return false;
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
  if (!save_image(params.output, params.logo ? add_logo(out) : out, error))
    return false;

  // done
  return true;
}

// bilateral params
struct bilateral_params {
  string image         = "";
  string output        = "";
  bool   logo          = false;
  float  spatial_sigma = 0.0f;
  float  range_sigma   = 0.0f;
};

// bilateral images
bool bilateral_images(const bilateral_params& params, string& error) {
  // load
  auto img = image<vec4f>{};
  if (path_extension(params.image) == ".ypreset") {
    if (!make_image_preset(path_basename(params.image), img, error))
      return false;
  } else {
    if (!load_image(params.image, img, error)) return false;
  }

  // edit alpha
  auto out = filter_bilateral(img, params.spatial_sigma, params.range_sigma);

  // save
  if (!save_image(params.output, params.logo ? add_logo(out) : out, error))
    return false;

  // done
  return true;
}

int main(int argc, const char* argv[]) {
  // command line parameters
  auto convert   = convert_params{};
  auto tonemap   = tonemap_params{};
  auto resize    = resize_params{};
  auto diff      = diff_params{};
  auto setalpha  = setalpha_params{};
  auto bilateral = bilateral_params{};

  // parse command line
  auto  cli         = make_cli("yimage", "Transform images");
  auto& cli_convert = add_command(cli, "convert", "Convert images");
  add_option(cli_convert, "--logo/--no-logo", convert.logo, "Add logo");
  add_option(cli_convert, "--output,-o", convert.output, "Output image");
  add_option(cli_convert, "image", convert.image, "Input image", true);
  auto& cli_tonemap = add_command(cli, "tonemap", "Tonemap images");
  add_option(cli_tonemap, "--exposure,-e", tonemap.exposure, "Exposure");
  add_option(
      cli_tonemap, "--filmic/--no-filmic", tonemap.filmic, "Filmic curve");
  add_option(cli_tonemap, "--srgb/--no-srgb", tonemap.srgb, "Srgb curve");
  add_option(cli_tonemap, "--logo/--no-logo", tonemap.logo, "Add logo");
  add_option(cli_tonemap, "--output,-o", tonemap.output, "Output image");
  add_option(cli_tonemap, "image", tonemap.image, "Input image", true);
  auto& cli_resize = add_command(cli, "resize", "Resize images");
  add_option(cli_resize, "--width", resize.width,
      "resize size (0 to maintain aspect)");
  add_option(cli_resize, "--height", resize.height,
      "resize size (0 to maintain aspect)");
  add_option(cli_resize, "--logo/--no-logo", resize.logo, "Add logo");
  add_option(cli_resize, "--output,-o", resize.output, "Output image");
  add_option(cli_resize, "image", resize.image, "Input image", true);
  auto& cli_diff = add_command(cli, "diff", "Diff two images");
  add_option(cli_diff, "--signal", diff.signal, "Signal a diff as error");
  add_option(cli_diff, "--threshold,", diff.threshold, "Diff threshold");
  add_option(cli_diff, "--logo/--no-logo", diff.logo, "Add logo");
  add_option(cli_diff, "--output,-o", diff.output, "Output image");
  add_option(cli_diff, "image1", diff.image1, "Input image", true);
  add_option(cli_diff, "image2", diff.image2, "Input image", true);
  auto& cli_setalpha = add_command(cli, "setalpha", "Set alpha in images");
  add_option(cli_setalpha, "--alpha", setalpha.alpha, "Alpha filename", true);
  add_option(cli_setalpha, "--from-color/--no-from-color", setalpha.from_color,
      "Get alpha from color");
  add_option(cli_setalpha, "--to-color/--no-to-color", setalpha.to_color,
      "Set color as alpha");
  add_option(cli_setalpha, "--logo/--no-logo", tonemap.logo, "Add logo");
  add_option(cli_setalpha, "--output,-o", tonemap.output, "Output image");
  add_option(cli_setalpha, "image", tonemap.image, "Input image", true);
  auto& cli_bilateral = add_command(
      cli, "bilateral", "Apply bilateral filtering to images");
  add_option(cli_bilateral, "--spatial-sigma", bilateral.spatial_sigma,
      "blur spatial sigma");
  add_option(cli_bilateral, "--range-sigma", bilateral.range_sigma,
      "bilateral blur range sigma");
  add_option(cli_bilateral, "--logo/--no-logo", resize.logo, "Add logo");
  add_option(cli_bilateral, "--output,-o", resize.output, "Output image");
  add_option(cli_bilateral, "image", resize.image, "Input image", true);

  // parse cli
  parse_cli(cli, argc, argv);

  // dispatch commands
  auto command = get_command(cli);
  auto error   = string{};
  if (command == "convert") {
    if (!convert_images(convert, error)) print_fatal(error);
  } else if (command == "tonemap") {
    if (!tonemap_images(tonemap, error)) print_fatal(error);
  } else if (command == "resize") {
    if (!resize_images(resize, error)) print_fatal(error);
  } else if (command == "diff") {
    if (!diff_images(diff, error)) print_fatal(error);
  } else if (command == "setalpha") {
    if (!setalpha_images(setalpha, error)) print_fatal(error);
  } else if (command == "bilateral") {
    if (!bilateral_images(bilateral, error)) print_fatal(error);
  } else {
    print_fatal("unknown command " + command);
  }

  // done
  return 0;
}
