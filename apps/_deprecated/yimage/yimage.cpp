//
// LICENSE:
//
// Copyright (c) 2016 -- 2022 Fabio Pellacini
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
#include <yocto/yocto_gui.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_scene.h>
#include <yocto/yocto_sceneio.h>

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
void add_options(cli_command& cli, convert_params& params) {
  add_option(cli, "image", params.image, "Input image.");
  add_option(cli, "output", params.output, "Output image.");
  add_option(cli, "exposure", params.exposure, "Tonemap exposure.");
  add_option(cli, "filmic", params.filmic, "Tonemap filmic.");
  add_option(cli, "width", params.width, "Resize width.");
  add_option(cli, "height", params.height, "Resize height.");
}

// convert images
void run_convert(const convert_params& params) {
  // load
  auto image = load_image(params.image);

  // resize if needed
  if (params.width != 0 || params.height != 0) {
    image = resize_image(image, params.width, params.height);
  }

  // tonemap if needed
  if (image.linear && is_ldr_filename(params.output)) {
    image = tonemap_image(image, params.exposure, params.filmic);
  }

  // save
  save_image(params.output, image);
}

// view params
struct view_params {
  vector<string> images = {"image.png"};
  string         output = "out.png";
};

// Cli
void add_options(cli_command& cli, view_params& params) {
  add_option(cli, "images", params.images, "Input images.");
  add_option(cli, "output", params.output, "Output image.");
}

// view images
void run_view(const view_params& params) {
  // load
  auto images = vector<image_data>{};
  for (auto& image : params.images) images.push_back(load_image(image));

  // run viewer
  show_image_gui("yimage", params.images, images);
}

// grade params
struct grade_params : colorgrade_params {
  string image  = "image.png";
  string output = "out.png";
};

// Cli
void add_options(cli_command& cli, grade_params& params) {
  add_option(cli, "image", params.image, "Input image.");
  add_option(cli, "output", params.output, "Output image.");
}

// grade images
void run_grade(const grade_params& params) {
  // load image
  auto image = load_image(params.image);

  // run viewer
  show_colorgrade_gui("yimage", params.image, image);
}

// resize params
struct diff_params {
  string image1    = "image1.png";
  string image2    = "image2.png";
  string output    = "";
  bool   signal    = false;
  float  threshold = 0;
};

// Cli
void add_options(cli_command& cli, diff_params& params) {
  add_option(cli, "image1", params.image1, "Input image 1.");
  add_option(cli, "image2", params.image2, "Input image 2.");
  add_option(cli, "output", params.output, "Output image.");
  add_option(cli, "signal", params.signal, "Error on diff.");
  add_option(cli, "threshold", params.threshold, "Diff threshold.");
}

// resize images
void run_diff(const diff_params& params) {
  // load
  auto image1 = load_image(params.image1);
  auto image2 = load_image(params.image2);

  // check sizes
  if (image1.width != image2.width || image1.height != image2.height)
    throw io_error("different image sizes");

  // check types
  if (image1.linear != image2.linear) throw io_error("different image types");

  // compute diff
  auto diff = image_difference(image1, image2, true);

  // save
  if (params.output != "") save_image(params.output, diff);

  // check diff
  if (params.signal) {
    for (auto& c : diff.pixels) {
      if (max(xyz(c)) > params.threshold) {
        throw io_error("image content differ");
      }
    }
  }
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
void add_options(cli_command& cli, setalpha_params& params) {
  add_option(cli, "image", params.image, "Input image");
  add_option(cli, "alpha", params.alpha, "Alpha image");
  add_option(cli, "output", params.output, "Output image");
  add_option(cli, "from-color", params.from_color, "Alpha from color");
  add_option(cli, "from-black", params.from_black, "Alpha from black");
  add_option(cli, "to-color", params.to_color, "color from alpha");
}

// setalpha images
void run_setalpha(const setalpha_params& params) {
  // load
  auto image = load_image(params.image);
  auto alpha = load_image(params.alpha);

  // check sizes
  if (image.width != alpha.width || image.height != alpha.height)
    throw io_error("different image sizes");

  // check types
  if (image.linear != alpha.linear) throw io_error("different image types");

  // edit alpha
  auto out = make_image(image.width, image.height, image.linear);
  for (auto idx : range(image.pixels.size())) {
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
  save_image(params.output, out);
}

struct app_params {
  string          command  = "convert";
  convert_params  convert  = {};
  view_params     view     = {};
  grade_params    grade    = {};
  diff_params     diff     = {};
  setalpha_params setalpha = {};
};

// Main
int main(int argc, const char* argv[]) {
  try {
    // command line parameters
    auto params = app_params{};
    auto cli    = make_cli("yimage", "process and view images");
    add_command_var(cli, params.command);
    add_command(cli, "convert", params.convert, "convert images");
    add_command(cli, "view", params.view, "view images");
    add_command(cli, "grade", params.grade, "grade images");
    add_command(cli, "diff", params.diff, "diff two images");
    add_command(cli, "setalpha", params.setalpha, "set images alpha");
    parse_cli(cli, argc, argv);

    // dispatch commands
    if (params.command == "convert") {
      run_convert(params.convert);
    } else if (params.command == "view") {
      run_view(params.view);
    } else if (params.command == "grade") {
      run_grade(params.grade);
    } else if (params.command == "diff") {
      run_diff(params.diff);
    } else if (params.command == "setalpha") {
      run_setalpha(params.setalpha);
    } else {
      throw io_error("unknown command");
    }
  } catch (const std::exception& error) {
    print_error(error.what());
    return 1;
  }

  // done
  return 0;
}
