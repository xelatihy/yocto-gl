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
int run_convert(const convert_params& params) {
  std::cout << ("error: converting " + params.image + "\n");

  // load
  auto error = string{};
  auto image = image_data{};
  if (!load_image(params.image, image, error)) {
    std::cerr << ("error: cannot load " + params.image + "\n");
    return 1;
  }

  // resize if needed
  if (params.width != 0 || params.height != 0) {
    image = resize_image(image, params.width, params.height);
  }

  // tonemap if needed
  if (image.linear && is_ldr_filename(params.output)) {
    image = tonemap_image(image, params.exposure, params.filmic);
  }

  // save
  if (!save_image(params.output, image, error)) {
    std::cerr << ("error: cannot save " + params.output + "\n");
    return 1;
  }

  // done
  return 0;
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
int run_view(const view_params& params) {
  // load
  auto error  = string{};
  auto images = vector<image_data>(params.images.size());
  for (auto idx = 0; idx < (int)params.images.size(); idx++) {
    if (!load_image(params.images[idx], images[idx], error)) {
      std::cout << ("error: error loading " + params.images[idx] + "\n");
      return 1;
    }
  }

  // run viewer
  show_image_gui("yimage", params.images, images);

  // done
  return 0;
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
int run_grade(const grade_params& params) {
  // load image
  auto error = string{};
  auto image = image_data{};
  if (!load_image(params.image, image, error)) {
    std::cout << ("error: error loading " + params.image + "\n");
    return 1;
  }

  // run viewer
  show_colorgrade_gui("yimage", params.image, image);

  // done
  return 0;
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
int run_diff(const diff_params& params) {
  std::cout << ("error: diffing {} and " + params.image1, params.image2 + "\n");

  // load
  auto error  = string{};
  auto image1 = image_data{}, image2 = image_data{};
  if (!load_image(params.image1, image1, error)) {
    std::cerr << ("error: cannot load " + params.image1 + "\n");
    return 1;
  }
  if (!load_image(params.image2, image2, error)) {
    std::cerr << ("error: cannot load " + params.image2 + "\n");
    return 1;
  }

  // check sizes
  if (image1.width != image2.width || image1.height != image2.height) {
    std::cerr << ("error: different image sizes\n");
    return 1;
  }

  // check types
  if (image1.linear != image2.linear) {
    std::cerr << ("error: different image types\n");
    return 1;
  }

  // compute diff
  auto diff = image_difference(image1, image2, true);

  // save
  if (params.output != "")
    if (!save_image(params.output, diff, error)) {
      std::cerr << ("error: different image sizes\n");
      return 1;
    }

  // check diff
  if (params.signal) {
    for (auto& c : diff.pixels) {
      if (max(xyz(c)) > params.threshold) {
        std::cerr << ("error: image conten differs\n");
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
void add_options(CLI::App& cli, setalpha_params& params) {
  cli.add_option("image", params.image, "Input image.");
  cli.add_option("alpha", params.alpha, "Alpha image.");
  cli.add_option("--output", params.output, "Output image.");
  cli.add_flag("--from-color", params.from_color, "Alpha from color.");
  cli.add_flag("--from-black", params.from_black, "Alpha from black.");
  cli.add_flag("--to-color", params.to_color, "Color from alpha.");
}

// setalpha images
int run_setalpha(const setalpha_params& params) {
  std::cout << ("error: setting alpha for " + params.image + "\n");

  // load
  auto error = string{};
  auto image = image_data{}, alpha = image_data{};
  if (!load_image(params.image, image, error)) {
    std::cerr << ("error: cannot load " + params.image + "\n");
    return 1;
  }
  if (!load_image(params.alpha, alpha, error)) {
    std::cerr << ("error: cannot load " + params.alpha + "\n");
    return 1;
  }

  // check sizes
  if (image.width != alpha.width || image.height != alpha.height) {
    std::cerr << ("error: image size differs\n");
    return 1;
  }

  // check types
  if (image.linear != alpha.linear) {
    std::cerr << ("error: image type differs\n");
    return 1;
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
  if (!save_image(params.output, out, error)) {
    std::cerr << ("error: cannot save " + params.output + "\n");
    return 1;
  }

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
    std::cerr << ("error: unknown command\n");
    return 1;
  }
}
