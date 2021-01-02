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

#include <yocto/yocto_color.h>
#include <yocto/yocto_commonio.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>
#if YOCTO_OPENGL == 1
#include <yocto_gui/yocto_glview.h>
#endif
using namespace yocto;

namespace yocto {

bool is_preset_filename(const string& filename) {
  return path_extension(filename) == ".ypreset";
}

bool is_preset_hdr(const string& type_) {
  auto type = path_basename(type_);
  return type.find("sky") != string::npos && type.find("sun") != string::npos;
}

bool make_image_preset(const string& type_, image_data& image, string& error) {
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
  } else if (type == "colormap") {
    image = make_colormapramp(width, height);
    // TODO(fabio): fix color space
    // image   = srgb_to_rgb(image);
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
    auto sub_imgs  = vector<image_data>(sub_types.size());
    for (auto i = 0; i < sub_imgs.size(); i++) {
      if (!make_image_preset(sub_types[i], sub_imgs[i], error)) return false;
    }
    auto montage_size = zero2i;
    for (auto& sub_img : sub_imgs) {
      montage_size.x += sub_img.width;
      montage_size.y = max(montage_size.y, sub_img.height);
    }
    image    = make_image(montage_size.x, montage_size.y, is_hdr(sub_imgs[0]));
    auto pos = 0;
    for (auto& sub_img : sub_imgs) {
      set_region(image, sub_img, pos, 0);
      pos += sub_img.width;
    }
  } else if (type == "images2") {
    auto sub_types = vector<string>{"sky", "sunsky"};
    auto sub_imgs  = vector<image_data>(sub_types.size());
    for (auto i = 0; i < sub_imgs.size(); i++) {
      if (!make_image_preset(sub_types[i], sub_imgs[i], error)) return false;
    }
    auto montage_size = zero2i;
    for (auto& sub_img : sub_imgs) {
      montage_size.x += sub_img.width;
      montage_size.y = max(montage_size.y, sub_img.height);
    }
    image    = make_image(montage_size.x, montage_size.y, is_hdr(sub_imgs[0]));
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
  auto image   = image_data{};
  auto ioerror = string{};
  if (is_preset_filename(params.image)) {
    if (!make_image_preset(path_basename(params.image), image, ioerror))
      return print_fatal(ioerror);
  } else {
    if (!load_image(params.image, image, ioerror)) return print_fatal(ioerror);
  }

  // resize if needed
  if (params.width != 0 || params.height != 0) {
    image = resize_image(image, params.width, params.height);
  }

  // tonemap if needed
  if (is_hdr(image) && is_ldr_filename(params.output)) {
    image = tonemap_image(image, params.exposure, params.filmic);
  }

  // apply logo
  if (params.logo) {
    image = add_logo(image);
  }

  // save
  if (!save_image(params.output, image, ioerror)) return print_fatal(ioerror);

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
  auto viewer_guard = make_imageviewer("yimage");
  auto viewer       = viewer_guard.get();

  // set image
  for (auto& filename : params.images) {
    // load
    auto image   = image_data{};
    auto ioerror = string{};
    if (is_preset_filename(filename)) {
      if (!make_image_preset(path_basename(filename), image, ioerror))
        return print_fatal(ioerror);
    } else {
      if (!load_image(filename, image, ioerror)) return print_fatal(ioerror);
    }

    // push image to the viewer
    set_image(viewer, filename, image);
  }

  // run view
  run_viewer(viewer);

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
  auto viewer_guard = make_imageviewer("yimage");
  auto viewer       = viewer_guard.get();

  // load image
  auto image   = image_data{};
  auto ioerror = string{};
  if (is_preset_filename(params.image)) {
    if (!make_image_preset(path_basename(params.image), image, ioerror))
      return print_fatal(ioerror);
  } else {
    if (!load_image(params.image, image, ioerror)) return print_fatal(ioerror);
  }

  // grade image
  auto graded = make_image(image.width, image.height, false);
  colorgrade_image_mt(graded, image, params);

  // set view
  set_image(viewer, params.image, graded);
  set_widgets(
      viewer, params.image, to_json(params), to_schema(params, "Color grade"));

  // set callback
  set_callback(viewer, [&params, &graded, &image, viewer](const string& name,
                           const json_value& uiparams, const gui_input&) {
    if (uiparams.is_null()) return;
    serialize_value(json_mode::from_json, (json_value&)uiparams,
        (colorgrade_params&)params, "");
    colorgrade_image_mt(graded, image, params);
    set_image(viewer, name, graded);
  });

  // run view
  run_viewer(viewer);

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
  auto image1 = image_data{}, image2 = image_data{};
  if (path_extension(params.image1) == ".ypreset") {
    if (!make_image_preset(path_basename(params.image1), image1, ioerror))
      return false;
  } else {
    if (!load_image(params.image1, image1, ioerror)) return false;
  }
  if (path_extension(params.image2) == ".ypreset") {
    if (!make_image_preset(path_basename(params.image2), image2, ioerror))
      return false;
  } else {
    if (!load_image(params.image2, image2, ioerror)) return false;
  }

  // check sizes
  if (image1.width != image2.width || image1.height != image2.height) {
    ioerror = "image sizes are different";
    return print_fatal(ioerror);
  }

  // check types
  if (is_hdr(image1) != is_hdr(image2) || is_ldr(image1) != is_ldr(image2)) {
    ioerror = "image types are different";
    return print_fatal(ioerror);
  }

  // compute diff
  auto diff = image_difference(image1, image2, true);

  // save
  if (params.output != "") {
    if (!save_image(
            params.output, params.logo ? add_logo(diff) : diff, ioerror))
      return print_fatal(ioerror);
  }

  // check diff
  if (params.signal) {
    for (auto& c : diff.hdr) {
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
  auto image = image_data{}, alpha = image_data{};
  if (path_extension(params.image) == ".ypreset") {
    if (!make_image_preset(path_basename(params.image), image, ioerror))
      return print_fatal(ioerror);
  } else {
    if (!load_image(params.image, image, ioerror)) return print_fatal(ioerror);
  }
  if (path_extension(params.alpha) == ".ypreset") {
    if (!make_image_preset(path_basename(params.alpha), alpha, ioerror))
      return print_fatal(ioerror);
  } else {
    if (!load_image(params.alpha, alpha, ioerror)) return print_fatal(ioerror);
  }

  // check sizes
  if (image.width != alpha.width || image.height != alpha.height) {
    ioerror = "image sizes are different";
    return print_fatal(ioerror);
  }

  // check types
  if (is_hdr(image) != is_hdr(alpha) || is_ldr(image) != is_ldr(alpha)) {
    ioerror = "image types are different";
    return print_fatal(ioerror);
  }

  // edit alpha
  auto out = make_image(image.width, image.height, is_hdr(image));
  for (auto idx = 0; idx < image.width * image.height; idx++) {
    if (is_hdr(image)) {
      auto a = params.from_color ? mean(xyz(image.hdr[idx])) : image.hdr[idx].w;
      if (params.to_color) {
        out.hdr[idx] = {a, a, a, a};
      } else {
        out.hdr[idx].w = a;
      }
    } else {
      auto a = params.from_color
                   ? float_to_byte(mean(xyz(byte_to_float(image.ldr[idx]))))
                   : image.ldr[idx].w;
      if (params.to_color) {
        out.ldr[idx] = {a, a, a, a};
      } else {
        out.ldr[idx].w = a;
      }
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
