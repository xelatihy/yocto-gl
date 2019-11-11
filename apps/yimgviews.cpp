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

#include "../yocto/yocto_commonio.h"
#include "../yocto/yocto_image.h"
#include "yocto_opengl.h"
using namespace yocto;

#include <future>

struct app_state {
  // original data
  string filename = "image.png";
  string outname  = "out.png";

  // image data
  image<vec4f> source = {};

  // diplay data
  image<vec4f>      display    = {};
  float             exposure   = 0;
  bool              filmic     = false;
  colorgrade_params params     = {};
  bool              colorgrade = false;

  // viewing properties
  opengl_image        glimage  = {};
  draw_glimage_params glparams = {};
};

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename Func>
inline void parallel_for(const vec2i& size, Func&& func) {
  auto             futures  = vector<std::future<void>>{};
  auto             nthreads = std::thread::hardware_concurrency();
  std::atomic<int> next_idx(0);
  for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
    futures.emplace_back(
        std::async(std::launch::async, [&func, &next_idx, size]() {
          while (true) {
            auto j = next_idx.fetch_add(1);
            if (j >= size.y) break;
            for (auto i = 0; i < size.x; i++) func({i, j});
          }
        }));
  }
  for (auto& f : futures) f.get();
}

void update_display(app_state& app) {
  if (app.display.size() != app.source.size()) app.display = app.source;
  parallel_for(app.source.size(), [&app](const vec2i& ij) {
    if (app.colorgrade) {
      app.display[ij] = colorgrade(app.source[ij], true, app.params);
    } else {
      app.display[ij] = tonemap(app.source[ij], app.exposure, app.filmic);
    }
  });
}

void draw(const opengl_window& win) {
  auto& app = *(app_state*)get_gluser_pointer(win);
  clear_glframebuffer(vec4f{0.15f, 0.15f, 0.15f, 1.0f});
  if (!app.glimage) update_glimage(app.glimage, app.display, false, false);
  app.glparams.window      = get_glwindow_size(win);
  app.glparams.framebuffer = get_glframebuffer_viewport(win);
  update_imview(app.glparams.center, app.glparams.scale, app.display.size(),
      app.glparams.window, app.glparams.fit);
  draw_glimage(app.glimage, app.glparams);
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
      app.glparams.center += mouse_pos - last_pos;
    }
    if (mouse_right) {
      app.glparams.scale *= powf(2, (mouse_pos.x - last_pos.x) * 0.001f);
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
  if (!load_image(app.filename, app.source))
    print_fatal("cannot load " + app.filename);
  update_display(app);

  // run ui
  run_ui(app);

  // done
  return 0;
}
