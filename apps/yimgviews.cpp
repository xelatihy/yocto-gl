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

#include "../yocto/yocto_cli.h"
#include "../yocto/yocto_image.h"
#include "../yocto_gui/yocto_gui.h"
using namespace ym;

#include <future>
using namespace std::string_literals;

struct app_state {
  // original data
  std::string filename = "image.png";
  std::string outname  = "out.png";

  // image data
  yimg::image<vec4f> source = {};

  // diplay data
  yimg::image<vec4f>      display    = {};
  float                   exposure   = 0;
  bool                    filmic     = false;
  yimg::colorgrade_params params     = {};
  bool                    colorgrade = false;

  // viewing properties
  ygui::image*       glimage  = new ygui::image{};
  ygui::image_params glparams = {};

  ~app_state() {
    if (glimage) delete glimage;
  }
};

void update_display(app_state* app) {
  if (app->display.size() != app->source.size()) app->display = app->source;
  if (app->colorgrade) {
    colorgrade_image_mt(app->display, app->source, true, app->params);
  } else {
    tonemap_image_mt(app->display, app->source, app->exposure, app->filmic);
  }
}

int main(int argc, const char* argv[]) {
  // prepare application
  auto app_guard = std::make_unique<app_state>();
  auto app       = app_guard.get();
  auto filenames = std::vector<std::string>{};

  // command line options
  auto cli = ycli::make_cli("yimgviews", "view images");
  add_option(cli, "--output,-o", app->outname, "image output");
  add_option(cli, "image", app->filename, "image filename", true);
  parse_cli(cli, argc, argv);

  // load image
  auto ioerror = ""s;
  if (!load_image(app->filename, app->source, ioerror)) {
    ycli::print_fatal(ioerror);
    return 1;
  }

  // update display
  update_display(app);

  // callbacks
  auto callbacks    = ygui::ui_callbacks{};
  callbacks.draw_cb = [app](ygui::window* win, const ygui::input& input) {
    app->glparams.window      = input.window_size;
    app->glparams.framebuffer = input.framebuffer_viewport;
    if (!is_initialized(app->glimage)) {
      init_glimage(app->glimage);
      set_glimage(app->glimage, app->display, false, false);
    }
    update_imview(app->glparams.center, app->glparams.scale,
        app->display.size(), app->glparams.window, app->glparams.fit);
    draw_glimage(app->glimage, app->glparams);
  };
  callbacks.widgets_cb = [app](ygui::window* win, const ygui::input& input) {
    auto edited = 0;
    if (begin_header(win, "tonemap")) {
      edited += draw_slider(win, "exposure", app->exposure, -5, 5);
      edited += draw_checkbox(win, "filmic", app->filmic);
      end_header(win);
    }
    if (begin_header(win, "colorgrade")) {
      auto& params = app->params;
      edited += draw_checkbox(win, "apply colorgrade", app->colorgrade);
      edited += draw_slider(win, "exposure", params.exposure, -5, 5);
      edited += draw_coloredit(win, "tint", params.tint);
      edited += draw_slider(win, "lincontrast", params.lincontrast, 0, 1);
      edited += draw_slider(win, "logcontrast", params.logcontrast, 0, 1);
      edited += draw_slider(win, "linsaturation", params.linsaturation, 0, 1);
      edited += draw_checkbox(win, "filmic", params.filmic);
      continue_line(win);
      edited += draw_checkbox(win, "srgb", params.srgb);
      edited += draw_slider(win, "contrast", params.contrast, 0, 1);
      edited += draw_slider(win, "saturation", params.saturation, 0, 1);
      edited += draw_slider(win, "shadows", params.shadows, 0, 1);
      edited += draw_slider(win, "midtones", params.midtones, 0, 1);
      edited += draw_slider(win, "highlights", params.highlights, 0, 1);
      edited += draw_coloredit(win, "shadows color", params.shadows_color);
      edited += draw_coloredit(win, "midtones color", params.midtones_color);
      edited += draw_coloredit(
          win, "highlights color", params.highlights_color);
      end_header(win);
    }
    if (begin_header(win, "inspect")) {
      draw_slider(win, "zoom", app->glparams.scale, 0.1, 10);
      draw_checkbox(win, "fit", app->glparams.fit);
      auto ij = get_image_coords(input.mouse_pos, app->glparams.center,
          app->glparams.scale, app->source.size());
      draw_dragger(win, "mouse", ij);
      auto img_pixel = zero4f, display_pixel = zero4f;
      if (ij.x >= 0 && ij.x < app->source.size().x && ij.y >= 0 &&
          ij.y < app->source.size().y) {
        img_pixel     = app->source[{ij.x, ij.y}];
        display_pixel = app->display[{ij.x, ij.y}];
      }
      draw_coloredit(win, "image", img_pixel);
      draw_dragger(win, "display", display_pixel);
      end_header(win);
    }
    if (edited) {
      update_display(app);
      if (!is_initialized(app->glimage)) init_glimage(app->glimage);
      set_glimage(app->glimage, app->display, false, false);
    }
  };
  callbacks.uiupdate_cb = [app](ygui::window* win, const ygui::input& input) {
    // handle mouse
    if (input.mouse_left && !input.widgets_active) {
      app->glparams.center += input.mouse_pos - input.mouse_last;
    }
    if (input.mouse_right && !input.widgets_active) {
      app->glparams.scale *= powf(
          2, (input.mouse_pos.x - input.mouse_last.x) * 0.001f);
    }
  };

  // run ui
  run_ui({1280, 720}, "yimgviews", callbacks);

  // done
  return 0;
}
