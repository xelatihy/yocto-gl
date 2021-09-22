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
#include <yocto/yocto_gui.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_scene.h>
#include <yocto/yocto_sceneio.h>

#include "ext/CLI11.hpp"

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
void add_options(CLI::App& cli, convert_params& params) {
  cli.add_option("image", params.image, "Input image.");
  cli.add_option("--output", params.output, "Output image.");
  cli.add_option("--exposure", params.exposure, "Tonemap exposure.")
      ->check(CLI::Range(-100, +100));
  cli.add_flag("--filmic", params.filmic, "Tonemap filmic.");
  cli.add_option("--width", params.width, "Resize width.")
      ->check(CLI::Range(1, int_max));
  cli.add_option("--height", params.height, "Resize height.")
      ->check(CLI::Range(1, int_max));
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
void add_options(CLI::App& cli, view_params& params) {
  cli.add_option("images", params.images, "Input images.");
  cli.add_option("--output", params.output, "Output image.");
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
void add_options(CLI::App& cli, grade_params& params) {
  cli.add_option("image", params.image, "Input image.");
  cli.add_option("--output", params.output, "Output image.");
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
void add_options(CLI::App& cli, diff_params& params) {
  cli.add_option("image1", params.image1, "Input image 1.");
  cli.add_option("image2", params.image2, "Input image 2.");
  cli.add_option("--output", params.output, "Output image.");
  cli.add_flag("--signal", params.signal, "Error on diff.");
  cli.add_option("--threshold", params.threshold, "Diff threshold.");
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
void add_options(CLI::App& cli, setalpha_params& params) {
  cli.add_option("image", params.image, "Input image.");
  cli.add_option("alpha", params.alpha, "Alpha image.");
  cli.add_option("--output", params.output, "Output image.");
  cli.add_flag("--from-color", params.from_color, "Alpha from color.");
  cli.add_flag("--from-black", params.from_black, "Alpha from black.");
  cli.add_flag("--to-color", params.to_color, "Color from alpha.");
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
  // command line parameters
  auto params = app_params{};
  auto cli    = CLI::App("Process and view images");
  add_options(
      *cli.add_subcommand("convert", "Convert images."), params.convert);
  add_options(*cli.add_subcommand("view", "View images."), params.view);
  add_options(*cli.add_subcommand("grade", "Grade images."), params.grade);
  add_options(*cli.add_subcommand("diff", "Diff two images."), params.diff);
  add_options(
      *cli.add_subcommand("setalpha", "Set alpha in images."), params.setalpha);
  cli.require_subcommand(1);
  try {
    cli.parse(argc, argv);
    params.command = cli.get_subcommands().front()->get_name();
  } catch (const CLI::ParseError& e) {
    return cli.exit(e);
  }

  // dispatch commands
  try {
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
  } catch (const io_error& error) {
    std::cerr << "error: " << error.what() << "\n";
    return 1;
  }

  // done
  return 0;
}
