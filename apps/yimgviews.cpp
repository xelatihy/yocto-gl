//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#include "../yocto/yocto_common.h"
#include "../yocto/yocto_commonio.h"
#include "../yocto/yocto_image.h"
#include "yocto_opengl.h"
using namespace yocto;

struct app_state {
  // original data
  string filename = "image.png";
  string outname  = "out.png";

  // image data
  image<vec4f> source = {};

  // diplay data
  image<vec4f>      display          = {};
  tonemap_params    tonemap_prms     = {};
  colorgrade_params colorgrade_prms  = {};
  bool              apply_colorgrade = false;

  // viewing properties
  opengl_image        gl_image  = {};
  draw_glimage_params draw_prms = {};
};

void update_display(app_state& app) {
  if (app.display.size() != app.source.size()) app.display = app.source;
  auto regions = make_image_regions(app.source.size(), 128);
  parallel_foreach(regions, [&app](const image_region& region) {
    tonemap_region(app.display, app.source, region, app.tonemap_prms);
    if (app.apply_colorgrade) {
      colorgrade_region(app.display, app.display, region, app.colorgrade_prms);
    }
  });
}

void draw(const opengl_window& win) {
  auto& app = *(app_state*)get_gluser_pointer(win);
  clear_glframebuffer(vec4f{0.15f, 0.15f, 0.15f, 1.0f});
  if (!app.gl_image) update_glimage(app.gl_image, app.display, false, false);
  app.draw_prms.window      = get_glwindow_size(win);
  app.draw_prms.framebuffer = get_glframebuffer_viewport(win);
  update_imview(app.draw_prms.center, app.draw_prms.scale, app.display.size(),
      app.draw_prms.window, app.draw_prms.fit);
  draw_glimage(app.gl_image, app.draw_prms);
  swap_glbuffers(win);
}

void run_ui(app_state& app) {
  // window
  auto win = opengl_window();
  init_glwindow(win, {1280, 720}, "yimview", &app, draw);

  // window values
  auto mouse_pos = zero2f, last_pos = zero2f;
  while (!should_glwindow_close(win)) {
    last_pos         = mouse_pos;
    mouse_pos        = get_glmouse_pos(win);
    auto mouse_left  = get_glmouse_left(win);
    auto mouse_right = get_glmouse_right(win);

    // handle mouse
    if (mouse_left) {
      app.draw_prms.center += mouse_pos - last_pos;
    }
    if (mouse_right) {
      app.draw_prms.scale *= powf(2, (mouse_pos.x - last_pos.x) * 0.001f);
    }

    // draw
    draw(win);

    // event hadling
    process_glevents(win);
  }

  // cleanup
  delete_glwindow(win);
}

int main(int argc, const char* argv[]) {
  // prepare application
  auto app       = app_state();
  auto filenames = vector<string>{};

  // command line options
  auto cli = make_cli("yimgview", "view images");
  add_cli_option(cli, "--output,-o", app.outname, "image output");
  add_cli_option(cli, "image", app.filename, "image filename", true);
  if (!parse_cli(cli, argc, argv)) exit(1);

  // load image
  load_image(app.filename, app.source);
  update_display(app);

  // run ui
  run_ui(app);

  // done
  return 0;
}
