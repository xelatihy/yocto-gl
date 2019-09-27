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

#include "../yocto/yocto_image.h"
#include "../yocto/yocto_utils.h"
#include "yocto_opengl.h"
using namespace yocto;

#include <atomic>
#include <future>
#include <thread>

struct app_state {
  // original data
  string filename = "image.png";
  string outname  = "out.png";

  // image data
  image<vec4f> img = {};

  // diplay data
  image<vec4f>   display = {};
  tonemap_params    tonemap_prms    = {};
  colorgrade_params colorgrade_prms = {};
  bool              colorgrade      = false;

  // viewing properties
  vec2f image_center = zero2f;
  float image_scale  = 1;
  bool  zoom_to_fit  = false;
  opengl_texture gl_txt  = {};
};

void update_display(app_state& app) {
  if (app.display.size() != app.img.size()) app.display = app.img;
  auto regions = make_regions(app.img.size(), 128);
  parallel_foreach(regions, [&app](const image_region& region) {
    tonemap(app.display, app.img, region, app.tonemap_prms);
    if (app.colorgrade) {
      colorgrade(app.display, app.display, region, app.colorgrade_prms);
    }
  });
}

void draw(const opengl_window& win) {
  auto& app      = *(app_state*)get_gluser_pointer(win);
  auto  win_size = get_glwindow_size(win);
  auto  fb_view  = get_glframebuffer_viewport(win);
  set_glviewport(fb_view);
  clear_glframebuffer(vec4f{0.15f, 0.15f, 0.15f, 1.0f});
  if (!app.gl_txt) {
    init_gltexture(app.gl_txt, app.display, false, false, false);
  }
  update_imview(app.image_center, app.image_scale, app.display.size(),
      win_size, app.zoom_to_fit);
  draw_glimage_background(
      app.gl_txt, win_size.x, win_size.y, app.image_center, app.image_scale);
  set_glblending(true);
  draw_glimage(
      app.gl_txt, win_size.x, win_size.y, app.image_center, app.image_scale);
  set_glblending(false);
  swap_glbuffers(win);
}

void run_ui(app_state& app) {
  // window
  auto win = opengl_window();
  init_glwindow(win, {1280, 720}, "yimview", &app, draw);

  // window values
  auto mouse_pos = zero2f, last_pos = zero2f;
  while (!should_glwindow_close(win)) {
    last_pos            = mouse_pos;
    mouse_pos           = get_glmouse_pos(win);
    auto mouse_left     = get_glmouse_left(win);
    auto mouse_right    = get_glmouse_right(win);

    // handle mouse
    if (mouse_left) {
      app.image_center += mouse_pos - last_pos;
    }
    if (mouse_right) {
      app.image_scale *= powf(2, (mouse_pos.x - last_pos.x) * 0.001f);
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
  load_image(app.filename, app.img);
  update_display(app);

  // run ui
  run_ui(app);

  // done
  return 0;
}
