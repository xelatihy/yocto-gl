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

// Cli
void add_command(cli_command& cli, const string& name, convert_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_positional(cmd, "image", value.image, "Input image.");
  add_optional(cmd, "output", value.output, "Output image.", {}, "o");
  add_optional(
      cmd, "exposure", value.exposure, "Tonemap exposure.", {-100, +100}, "e");
  add_optional(cmd, "filmic", value.filmic, "Tonemap filmic.", {}, "f");
  add_optional(cmd, "width", value.width, "Resize width.", {1, int_max}, "w");
  add_optional(
      cmd, "height", value.height, "Resize height.", {1, int_max}, "h");
  add_optional(cmd, "logo", value.logo, "Add logo.");
}

// convert images
int run_convert(const convert_params& params) {
  // load
  auto image   = image_data{};
  auto ioerror = string{};
  if (!load_image(params.image, image, ioerror)) return print_fatal(ioerror);

  // resize if needed
  if (params.width != 0 || params.height != 0) {
    image = resize_image(image, params.width, params.height);
  }

  // tonemap if needed
  if (image.linear && is_ldr_filename(params.output)) {
    image = tonemap_image(image, params.exposure, params.filmic, true);
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

// Cli
void add_command(cli_command& cli, const string& name, view_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_positional(cmd, "images", value.images, "Input images.");
  add_optional(cmd, "output", value.output, "Output image.", {}, "o");
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
  auto viewer = make_imageviewer("yimage");

  // set image
  for (auto& filename : params.images) {
    // load
    auto image   = image_data{};
    auto ioerror = string{};
    if (!load_image(filename, image, ioerror)) return print_fatal(ioerror);

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

// Cli
void add_command(cli_command& cli, const string& name, grade_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_positional(cmd, "image", value.image, "Input image.");
  add_optional(cmd, "output", value.output, "Output image.", {}, "o");
}

#ifndef YOCTO_OPENGL

// grade images
int run_grade(const grade_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

// Parameter conversions
void from_params(const gui_params& uiparams, colorgrade_params& params) {
  params.exposure         = uiparams.at("exposure");
  params.tint             = uiparams.at("tint");
  params.lincontrast      = uiparams.at("lincontrast");
  params.logcontrast      = uiparams.at("logcontrast");
  params.linsaturation    = uiparams.at("linsaturation");
  params.filmic           = uiparams.at("filmic");
  params.srgb             = uiparams.at("srgb");
  params.contrast         = uiparams.at("contrast");
  params.saturation       = uiparams.at("saturation");
  params.shadows          = uiparams.at("shadows");
  params.midtones         = uiparams.at("midtones");
  params.highlights       = uiparams.at("highlights");
  params.shadows_color    = uiparams.at("shadows_color");
  params.midtones_color   = uiparams.at("midtones_color");
  params.highlights_color = uiparams.at("highlights_color");
}
void to_params(gui_params& uiparams, const colorgrade_params& params) {
  uiparams["exposure"]         = {params.exposure, {-10, 10}};
  uiparams["tint"]             = {params.tint, true};
  uiparams["lincontrast"]      = {params.lincontrast, {0, 1}};
  uiparams["logcontrast"]      = {params.logcontrast, {0, 1}};
  uiparams["linsaturation"]    = {params.linsaturation, {0, 1}};
  uiparams["filmic"]           = {params.filmic};
  uiparams["srgb"]             = {params.srgb};
  uiparams["contrast"]         = {params.contrast, {0, 1}};
  uiparams["saturation"]       = {params.saturation, {0, 1}};
  uiparams["shadows"]          = {params.shadows, {0, 1}};
  uiparams["midtones"]         = {params.midtones, {0, 1}};
  uiparams["highlights"]       = {params.highlights, {0, 1}};
  uiparams["shadows_color"]    = {params.shadows_color, true};
  uiparams["midtones_color"]   = {params.midtones_color, true};
  uiparams["highlights_color"] = {params.highlights_color, true};
}

// grade images
int run_grade(const grade_params& params) {
  // open viewer
  auto viewer = make_imageviewer("yimage");

  // load image
  auto image   = image_data{};
  auto ioerror = string{};
  if (!load_image(params.image, image, ioerror)) return print_fatal(ioerror);

  // grade image
  auto graded = make_image(image.width, image.height, false, true);
  colorgrade_image(graded, image, params);

  // set view
  set_image(viewer, params.image, graded);
  auto uiparams = gui_params{};
  to_params(uiparams, params);
  set_params(viewer, params.image, "Color grade", uiparams);

  // set callback
  set_params_callback(
      viewer, [&](const string& name, const gui_params& uiparams) {
        if (uiparams.empty()) return;
        auto gparams = params;
        from_params(uiparams, gparams);
        colorgrade_image_mt(graded, image, gparams);
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

// Cli
void add_command(cli_command& cli, const string& name, diff_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_positional(cmd, "image1", value.image1, "Input image 1.");
  add_positional(cmd, "image2", value.image2, "Input image 2.");
  add_optional(cmd, "output", value.output, "Output image.", {}, "o");
  add_optional(cmd, "signal", value.signal, "Error on diff.");
  add_optional(cmd, "threshold", value.threshold, "Diff threshold.");
  add_optional(cmd, "logo", value.logo, "Add logo.");
}

// resize images
int run_diff(const diff_params& params) {
  // load
  auto ioerror = string{};
  auto image1 = image_data{}, image2 = image_data{};
  if (!load_image(params.image1, image1, ioerror)) return false;
  if (!load_image(params.image2, image2, ioerror)) return false;

  // check sizes
  if (image1.width != image2.width || image1.height != image2.height) {
    ioerror = "image sizes are different";
    return print_fatal(ioerror);
  }

  // check types
  if (!image1.pixelsf.empty() != !image2.pixelsf.empty() ||
      !image1.pixelsf.empty() != !image2.pixelsf.empty()) {
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
    for (auto& c : diff.pixelsf) {
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

// Cli
void add_command(cli_command& cli, const string& name, setalpha_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_positional(cmd, "image", value.image, "Input image.");
  add_positional(cmd, "alpha", value.alpha, "Alpha image.");
  add_optional(cmd, "output", value.output, "Output image.", {}, "o");
  add_optional(cmd, "from-color", value.from_color, "Alpha from color.");
  add_optional(cmd, "to-color", value.to_color, "Color from alpha.");
  add_optional(cmd, "logo", value.logo, "Add logo.");
}

// setalpha images
int run_setalpha(const setalpha_params& params) {
  // load
  auto ioerror = string{};
  auto image = image_data{}, alpha = image_data{};
  if (!load_image(params.image, image, ioerror)) return print_fatal(ioerror);
  if (!load_image(params.alpha, alpha, ioerror)) return print_fatal(ioerror);

  // check sizes
  if (image.width != alpha.width || image.height != alpha.height) {
    ioerror = "image sizes are different";
    return print_fatal(ioerror);
  }

  // check types
  if (!image.pixelsf.empty() != !alpha.pixelsf.empty() ||
      !image.pixelsf.empty() != !alpha.pixelsf.empty()) {
    ioerror = "image types are different";
    return print_fatal(ioerror);
  }

  // check types
  if (image.linear != alpha.linear || !image.linear != !alpha.linear) {
    ioerror = "image types are different";
    return print_fatal(ioerror);
  }

  // edit alpha
  auto out = make_image(
      image.width, image.height, image.linear, !image.pixelsf.empty());
  for (auto j = 0; j < image.height; j++) {
    for (auto i = 0; i < image.width; i++) {
      auto calpha = get_pixel(alpha, i, j);
      auto alpha_ = params.from_color ? mean(xyz(calpha)) : calpha.w;
      if (params.to_color) {
        set_pixel(out, i, j, {alpha_, alpha_, alpha_, alpha_});
      } else {
        auto color = get_pixel(image, i, j);
        color.w    = alpha_;
        set_pixel(out, i, j, color);
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

// Cli
void add_commands(cli_command& cli, const string& name, app_params& value,
    const string& usage) {
  cli = make_cli(name, usage);
  add_command_name(cli, "command", value.command, "Command.");
  add_command(cli, "convert", value.convert, "Convert images.");
  add_command(cli, "view", value.view, "View images.");
  add_command(cli, "grade", value.grade, "Grade images.");
  add_command(cli, "diff", value.diff, "Diff two images.");
  add_command(cli, "setalpha", value.setalpha, "Set alpha in images.");
}

// Parse cli
void parse_cli(app_params& params, int argc, const char** argv) {
  auto cli = cli_command{};
  add_commands(cli, "yimage", params, "Process and view images.");
  parse_cli(cli, argc, argv);
}

int main(int argc, const char* argv[]) {
  // command line parameters
  auto params = app_params{};
  parse_cli(params, argc, argv);

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
