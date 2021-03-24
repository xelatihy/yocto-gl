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

#include <yocto/yocto_cli.h>
#include <yocto/yocto_color.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_scene.h>
#include <yocto/yocto_sceneio.h>
#if YOCTO_OPENGL == 1
#include <yocto_gui/yocto_glview.h>
#endif
using namespace yocto;

// convert params
struct convert_params {
  string image    = "image.png";
  string output   = "out.png";
  float  exposure = 0;
  bool   filmic   = false;
  int    width    = 0;
  int    height   = 0;
};

// Cli
void add_command(cli_command& cli, const string& name, convert_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_argument(cmd, "image", value.image, "Input image.");
  add_option(cmd, "output", value.output, "Output image.");
  add_option(
      cmd, "exposure", value.exposure, "Tonemap exposure.", {-100, +100});
  add_option(cmd, "filmic", value.filmic, "Tonemap filmic.");
  add_option(cmd, "width", value.width, "Resize width.", {1, int_max});
  add_option(cmd, "height", value.height, "Resize height.", {1, int_max});
}

// convert images
int run_convert(const convert_params& params) {
  // load
  auto image   = color_image{};
  auto ioerror = string{};
  if (!load_image(params.image, image, ioerror)) return print_fatal(ioerror);

  // resize if needed
  if (params.width != 0 || params.height != 0) {
    image = resize_image(image, params.width, params.height);
  }

  // tonemap if needed
  if (image.linear && is_ldr_filename(params.output)) {
    image = tonemap_image(image, params.exposure, params.filmic);
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
};

// Cli
void add_command(cli_command& cli, const string& name, view_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_argument(cmd, "images", value.images, "Input images.");
  add_option(cmd, "output", value.output, "Output image.");
}

#ifndef YOCTO_OPENGL

// view images
int run_view(const view_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

// view images
int run_view(const view_params& params) {
  // load
  auto images = vector<color_image>(params.images.size());
  for (auto idx = 0; idx < (int)params.images.size(); idx++) {
    auto ioerror = string{};
    if (!load_image(params.images[idx], images[idx], ioerror))
      return print_fatal(ioerror);
  }

  // run viewer
  view_images("yimage", params.images, images);

  // done
  return 0;
}

#endif

// grade params
struct grade_params : colorgrade_params {
  string image  = "image.png";
  string output = "out.png";
};

// Cli
void add_command(cli_command& cli, const string& name, grade_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_argument(cmd, "image", value.image, "Input image.");
  add_option(cmd, "output", value.output, "Output image.");
}

#ifndef YOCTO_OPENGL

// grade images
int run_grade(const grade_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

// grade images
int run_grade(const grade_params& params) {
  // load image
  auto image   = color_image{};
  auto ioerror = string{};
  if (!load_image(params.image, image, ioerror)) return print_fatal(ioerror);

  // run viewer
  colorgrade_image("yimage", params.image, image);

  // done
  return 0;
}

#endif

// resize params
struct diff_params {
  string image1    = "image1.png";
  string image2    = "image2.png";
  string output    = "";
  bool   signal    = false;
  float  threshold = 0;
};

// Cli
void add_command(cli_command& cli, const string& name, diff_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_argument(cmd, "image1", value.image1, "Input image 1.");
  add_argument(cmd, "image2", value.image2, "Input image 2.");
  add_option(cmd, "output", value.output, "Output image.");
  add_option(cmd, "signal", value.signal, "Error on diff.");
  add_option(cmd, "threshold", value.threshold, "Diff threshold.");
}

// resize images
int run_diff(const diff_params& params) {
  // load
  auto ioerror = string{};
  auto image1 = color_image{}, image2 = color_image{};
  if (!load_image(params.image1, image1, ioerror)) return false;
  if (!load_image(params.image2, image2, ioerror)) return false;

  // check sizes
  if (image1.width != image2.width || image1.height != image2.height) {
    ioerror = "image sizes are different";
    return print_fatal(ioerror);
  }

  // check types
  if (image1.linear != image2.linear) {
    ioerror = "image types are different";
    return print_fatal(ioerror);
  }

  // compute diff
  auto diff = image_difference(image1, image2, true);

  // save
  if (params.output != "") {
    if (!save_image(params.output, diff, ioerror)) return print_fatal(ioerror);
  }

  // check diff
  if (params.signal) {
    for (auto& c : diff.pixels) {
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
  bool   from_color = false;
  bool   from_black = false;
  bool   to_color   = false;
};

// Cli
void add_command(cli_command& cli, const string& name, setalpha_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_argument(cmd, "image", value.image, "Input image.");
  add_argument(cmd, "alpha", value.alpha, "Alpha image.");
  add_option(cmd, "output", value.output, "Output image.");
  add_option(cmd, "from-color", value.from_color, "Alpha from color.");
  add_option(cmd, "from-black", value.from_black, "Alpha from black.");
  add_option(cmd, "to-color", value.to_color, "Color from alpha.");
}

// setalpha images
int run_setalpha(const setalpha_params& params) {
  // load
  auto ioerror = string{};
  auto image = color_image{}, alpha = color_image{};
  if (!load_image(params.image, image, ioerror)) return print_fatal(ioerror);
  if (!load_image(params.alpha, alpha, ioerror)) return print_fatal(ioerror);

  // check sizes
  if (image.width != alpha.width || image.height != alpha.height) {
    ioerror = "image sizes are different";
    return print_fatal(ioerror);
  }

  // check types
  if (image.linear != alpha.linear) {
    ioerror = "image types are different";
    return print_fatal(ioerror);
  }

  // check types
  if (image.linear != alpha.linear || !image.linear != !alpha.linear) {
    ioerror = "image types are different";
    return print_fatal(ioerror);
  }

  // edit alpha
  auto out = make_image(image.width, image.height, image.linear);
  for (auto idx = (size_t)0; idx < image.pixels.size(); idx++) {
    auto calpha = alpha.pixels[idx];
    auto alpha_ = params.from_color   ? mean(xyz(calpha))
                  : params.from_black ? (mean(xyz(calpha)) > 0.01 ? 1.0f : 0.0f)
                                      : calpha.w;
    if (params.to_color) {
      out.pixels[idx] = {alpha_, alpha_, alpha_, alpha_};
    } else {
      auto color      = image.pixels[idx];
      color.w         = alpha_;
      out.pixels[idx] = color;
    }
  }

  // save
  if (!save_image(params.output, out, ioerror)) return print_fatal(ioerror);

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
