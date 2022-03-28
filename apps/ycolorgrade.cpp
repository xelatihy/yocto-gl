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
using namespace std::string_literals;

// main function
void run(const vector<string>& args) {
  // parameters
  auto imagename   = "image.png"s;
  auto outname     = "out.png"s;
  auto paramsname  = ""s;
  auto interactive = true;
  auto dumpname    = ""s;
  auto params      = colorgrade_params{};

  // parse command line
  auto cli = make_cli("ycolorgrade", "adjust image colors");
  add_option(cli, "image", imagename, "Input image.");
  add_option(cli, "output", outname, "Output image.");
  add_option(cli, "params", paramsname, "params filename");
  add_option(cli, "interactive", interactive, "Run interactively.");
  add_option(cli, "dumpparams", dumpname, "dump params filename");
  parse_cli(cli, args);

  // load config
  if (!paramsname.empty()) {
    update_colorgrade_params(paramsname, params);
  }

  // dump config
  if (!dumpname.empty()) {
    save_colorgrade_params(dumpname, params);
  }

  // load image
  auto image = load_image(imagename);

  // switch between interactive and offline
  if (!interactive) {
    // apply color grade
    image = colorgrade_image(image, params);

    // save image
    save_image(outname, image);
  } else {
#ifdef YOCTO_OPENGL

    // color grading parameters
    auto params = colorgrade_params{};

    // display image
    auto display = make_image(image.width, image.height, false);
    colorgrade_image_mt(display, image, params);

    // opengl image
    auto glimage  = glimage_state{};
    auto glparams = glimage_params{};

    // callbacks
    auto callbacks = gui_callbacks{};
    callbacks.init = [&](const gui_input& input) {
      init_image(glimage);
      set_image(glimage, display);
    };
    callbacks.clear = [&](const gui_input& input) { clear_image(glimage); };
    callbacks.draw  = [&](const gui_input& input) {
      update_image_params(input, image, glparams);
      draw_image(glimage, glparams);
    };
    callbacks.widgets = [&](const gui_input& input) {
      if (draw_gui_header("colorgrade")) {
        auto edited = 0;
        edited += draw_gui_slider("exposure", params.exposure, -5, 5);
        edited += draw_gui_coloredit("tint", params.tint);
        edited += draw_gui_slider("lincontrast", params.lincontrast, 0, 1);
        edited += draw_gui_slider("logcontrast", params.logcontrast, 0, 1);
        edited += draw_gui_slider("linsaturation", params.linsaturation, 0, 1);
        edited += draw_gui_checkbox("filmic", params.filmic);
        continue_gui_line();
        edited += draw_gui_checkbox("srgb", params.srgb);
        edited += draw_gui_slider("contrast", params.contrast, 0, 1);
        edited += draw_gui_slider("saturation", params.saturation, 0, 1);
        edited += draw_gui_slider("shadows", params.shadows, 0, 1);
        edited += draw_gui_slider("midtones", params.midtones, 0, 1);
        edited += draw_gui_slider("highlights", params.highlights, 0, 1);
        edited += draw_gui_coloredit("shadows color", params.shadows_color);
        edited += draw_gui_coloredit("midtones color", params.midtones_color);
        edited += draw_gui_coloredit(
            "highlights color", params.highlights_color);
        end_gui_header();
        if (edited) {
          colorgrade_image_mt(display, image, params);
          set_image(glimage, display);
        }
      }
      draw_image_widgets(input, image, display, glparams);
    };
    callbacks.uiupdate = [&glparams](const gui_input& input) {
      uiupdate_image_params(input, glparams);
    };

    // run ui
    show_gui_window({1280 + 320, 720}, "ycolorgrade - " + imagename, callbacks);
#else
    throw io_error{"Interactive requires OpenGL"};
#endif
  }
}

// Main
int main(int argc, const char* argv[]) {
  try {
    run({argv, argv + argc});
    return 0;
  } catch (const std::exception& error) {
    print_error(error.what());
    return 1;
  }
}
