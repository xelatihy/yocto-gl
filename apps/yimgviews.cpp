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
#include "yocto_opengl.h"
using namespace yocto;

#include <future>
using namespace std;

#include "ext/CLI11.hpp"

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
  shared_ptr<opengl_image> glimage  = nullptr;
  draw_glimage_params      glparams = {};
};

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the pixel index as a vec2i.
template <typename Func>
inline void parallel_for(const vec2i& size, Func&& func) {
  auto        futures  = vector<future<void>>{};
  auto        nthreads = thread::hardware_concurrency();
  atomic<int> next_idx(0);
  for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
    futures.emplace_back(async(launch::async, [&func, &next_idx, size]() {
      while (true) {
        auto j = next_idx.fetch_add(1);
        if (j >= size.y) break;
        for (auto i = 0; i < size.x; i++) func({i, j});
      }
    }));
  }
  for (auto& f : futures) f.get();
}

void update_display(shared_ptr<app_state> app) {
  if (app->display.size() != app->source.size()) app->display = app->source;
  parallel_for(app->source.size(), [app](const vec2i& ij) {
    if (app->colorgrade) {
      app->display[ij] = colorgrade(app->source[ij], true, app->params);
    } else {
      app->display[ij] = tonemap(app->source[ij], app->exposure, app->filmic);
    }
  });
}

int run_app(int argc, const char* argv[]) {
  // prepare application
  auto app       = make_shared<app_state>();
  auto filenames = vector<string>{};

  // command line options
  auto cli = CLI::App{"view images"};
  cli.add_option("--output,-o", app->outname, "image output");
  cli.add_option("image", app->filename, "image filename")->required();
  try {
    cli.parse(argc, argv);
  } catch (CLI::ParseError& e) {
    return cli.exit(e);
  }

  // load image
  load_image(app->filename, app->source);

  // update display
  update_display(app);

  // create window
  auto win = make_glwindow({1280, 720}, "yimgviews", false);

  // set callbacks
  set_draw_glcallback(
      win, [app](shared_ptr<opengl_window> win, const opengl_input& input) {
        app->glparams.window      = input.window_size;
        app->glparams.framebuffer = input.framebuffer_viewport;
        if (!app->glimage) {
          app->glimage = make_glimage();
          set_glimage(app->glimage, app->display, false, false);
        }
        update_imview(app->glparams.center, app->glparams.scale,
            app->display.size(), app->glparams.window, app->glparams.fit);
        draw_glimage(app->glimage, app->glparams);
      });
  set_uiupdate_glcallback(
      win, [app](shared_ptr<opengl_window> win, const opengl_input& input) {
        // handle mouse
        if (input.mouse_left) {
          app->glparams.center += input.mouse_pos - input.mouse_last;
        }
        if (input.mouse_right) {
          app->glparams.scale *= powf(
              2, (input.mouse_pos.x - input.mouse_last.x) * 0.001f);
        }
      });

  // run ui
  run_ui(win);

  // cleanup
  clear_glwindow(win);

  // done
  return 0;
}

int main(int argc, const char* argv[]) {
  try {
    return run_app(argc, argv);
  } catch (std::exception& e) {
    fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
