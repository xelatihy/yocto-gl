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

#include <deque>
#include <future>
using namespace std;

#include "ext/CLI11.hpp"
#include "ext/filesystem.hpp"
namespace fs = ghc::filesystem;

struct image_stats {
  vec4f         min       = zero4f;
  vec4f         max       = zero4f;
  vec4f         average   = zero4f;
  vector<vec3f> histogram = {};
};

struct app_state {
  // original data
  string name     = "";
  string filename = "";
  string outname  = "";

  // image data
  image<vec4f> source  = {};
  image<vec4f> display = {};

  // image stats
  image_stats source_stats  = {};
  image_stats display_stats = {};

  // tonemapping values
  float             exposure   = 0;
  bool              filmic     = false;
  colorgrade_params params     = {};
  bool              colorgrade = false;

  // viewing properties
  shared_ptr<opengl_image> glimage   = nullptr;
  draw_glimage_params      glparams  = {};
  bool                     glupdated = true;

  // loading status
  atomic<bool> ok     = false;
  future<void> loader = {};
  string       status = "";
  string       error  = "";
};

// app states
struct app_states {
  // data
  vector<shared_ptr<app_state>> states   = {};
  shared_ptr<app_state>         selected = nullptr;
  deque<shared_ptr<app_state>>  loading  = {};

  // default options
  float             exposure = 0;
  bool              filmic   = false;
  colorgrade_params params   = {};
};

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
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

// compute min/max
void compute_stats(
    image_stats& stats, const image<vec4f>& img, bool linear_hdr) {
  auto max_histo = linear_hdr ? 8 : 1;
  stats.min      = vec4f{flt_max};
  stats.max      = vec4f{flt_min};
  stats.average  = zero4f;
  stats.histogram.assign(256, zero3f);
  for (auto& p : img) {
    stats.min = min(stats.min, p);
    stats.max = max(stats.max, p);
    stats.average += p;
    stats.histogram[(int)(clamp(p.x / max_histo, 0.f, 1.f) * 255)].x += 1;
    stats.histogram[(int)(clamp(p.y / max_histo, 0.f, 1.f) * 255)].y += 1;
    stats.histogram[(int)(clamp(p.z / max_histo, 0.f, 1.f) * 255)].z += 1;
  }
  auto num_pixels = (size_t)img.size().x * (size_t)img.size().y;
  for (auto& v : stats.histogram) v /= num_pixels;
  stats.average /= num_pixels;
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
  compute_stats(app->display_stats, app->display, false);
  app->glupdated = true;
}

// add a new image
void load_image_async(shared_ptr<app_states> apps, const string& filename) {
  auto app      = make_shared<app_state>();
  app->filename = filename;
  app->outname  = fs::path(filename).replace_extension(".display.png").string();
  app->name     = fs::path(filename).filename();
  app->exposure = apps->exposure;
  app->filmic   = apps->filmic;
  app->params   = apps->params;
  app->status   = "loading";
  app->loader   = async(launch::async, [app]() {
    load_image(app->filename, app->source);
    compute_stats(
        app->source_stats, app->source, is_hdr_filename(app->filename));
    if (app->colorgrade) {
      app->display = colorgrade_image(app->display, true, app->params);
    } else {
      app->display = tonemap_image(app->source, app->exposure, app->filmic);
    }
    compute_stats(app->display_stats, app->display, false);
  });
  apps->states.push_back(app);
  apps->loading.push_back(app);
  if (!apps->selected) apps->selected = apps->states.front();
}

void draw_glwidgets(shared_ptr<opengl_window> win, shared_ptr<app_states> apps,
    const opengl_input& input) {
  static string load_path = "", save_path = "", error_message = "";
  if (draw_glfiledialog_button(win, "load", true, "load image", load_path,
          false, "./", "", "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    load_image_async(apps, load_path);
    load_path = "";
  }
  continue_glline(win);
  if (draw_glfiledialog_button(win, "save",
          apps->selected && apps->selected->ok, "save image", save_path, true,
          fs::path(save_path).parent_path(), fs::path(save_path).filename(),
          "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    auto app     = apps->selected;
    app->outname = save_path;
    try {
      save_image(app->outname, app->display);
    } catch (exception& e) {
      push_glmessage(win, "cannot save " + app->outname);
      log_glinfo(win, "cannot save " + app->outname);
      log_glinfo(win, e.what());
    }
    save_path = "";
  }
  continue_glline(win);
  if (draw_glbutton(win, "close", (bool)apps->selected)) {
    if (apps->selected->loader.valid()) return;
    apps->states.erase(
        std::find(apps->states.begin(), apps->states.end(), apps->selected));
    apps->selected = apps->states.empty() ? nullptr : apps->states.front();
  }
  continue_glline(win);
  if (draw_glbutton(win, "quit")) {
    set_close(win, true);
  }
  draw_glcombobox(win, "image", apps->selected, apps->states, false);
  if (!apps->selected) return;
  auto app = apps->selected;
  if (app->status != "") draw_gllabel(win, "status", app->status);
  if (app->error != "") draw_gllabel(win, "error", app->error);
  if (!app->ok) return;
  if (begin_glheader(win, "tonemap")) {
    auto edited = 0;
    edited += draw_glslider(win, "exposure", app->exposure, -5, 5);
    edited += draw_glcheckbox(win, "filmic", app->filmic);
    if (edited) update_display(app);
    end_glheader(win);
  }
  if (begin_glheader(win, "colorgrade")) {
    auto& params = app->params;
    auto  edited = 0;
    edited += draw_glcheckbox(win, "apply colorgrade", app->colorgrade);
    edited += draw_glslider(win, "exposure", params.exposure, -5, 5);
    edited += draw_glcoloredit(win, "tint", params.tint);
    edited += draw_glslider(win, "lincontrast", params.lincontrast, 0, 1);
    edited += draw_glslider(win, "logcontrast", params.logcontrast, 0, 1);
    edited += draw_glslider(win, "linsaturation", params.linsaturation, 0, 1);
    edited += draw_glcheckbox(win, "filmic", params.filmic);
    continue_glline(win);
    edited += draw_glcheckbox(win, "srgb", params.srgb);
    continue_glline(win);
    if (draw_glbutton(win, "auto wb")) {
      auto wb     = 1 / xyz(app->source_stats.average);
      params.tint = wb / max(wb);
      edited += 1;
    }
    edited += draw_glslider(win, "contrast", params.contrast, 0, 1);
    edited += draw_glslider(win, "saturation", params.saturation, 0, 1);
    edited += draw_glslider(win, "shadows", params.shadows, 0, 1);
    edited += draw_glslider(win, "midtones", params.midtones, 0, 1);
    edited += draw_glslider(win, "highlights", params.highlights, 0, 1);
    edited += draw_glcoloredit(win, "shadows color", params.shadows_color);
    edited += draw_glcoloredit(win, "midtones color", params.midtones_color);
    edited += draw_glcoloredit(
        win, "highlights color", params.highlights_color);
    if (edited) update_display(app);
    end_glheader(win);
  }
  if (begin_glheader(win, "inspect")) {
    draw_gllabel(win, "image", fs::path(app->filename).filename());
    draw_gllabel(win, "filename", app->filename);
    draw_gllabel(win, "outname", app->outname);
    draw_gllabel(win, "image",
        to_string(app->source.size().x) + " x " +
            to_string(app->source.size().y));
    draw_glslider(win, "zoom", app->glparams.scale, 0.1, 10);
    draw_glcheckbox(win, "fit", app->glparams.fit);
    auto ij = get_image_coords(input.mouse_pos, app->glparams.center,
        app->glparams.scale, app->source.size());
    draw_gldragger(win, "mouse", ij);
    auto img_pixel = zero4f, display_pixel = zero4f;
    if (ij.x >= 0 && ij.x < app->source.size().x && ij.y >= 0 &&
        ij.y < app->source.size().y) {
      img_pixel     = app->source[{ij.x, ij.y}];
      display_pixel = app->display[{ij.x, ij.y}];
    }
    draw_glcoloredit(win, "image", img_pixel);
    draw_gldragger(win, "display", display_pixel);
    draw_gldragger(win, "image min", app->source_stats.min);
    draw_gldragger(win, "image max", app->source_stats.max);
    draw_gldragger(win, "image avg", app->source_stats.average);
    draw_glhistogram(win, "image histo", app->source_stats.histogram);
    draw_gldragger(win, "display min", app->display_stats.min);
    draw_gldragger(win, "display max", app->display_stats.max);
    draw_gldragger(win, "display avg", app->display_stats.average);
    draw_glhistogram(win, "display histo", app->display_stats.histogram);
    end_glheader(win);
  }
}

void draw(shared_ptr<opengl_window> win, shared_ptr<app_states> apps,
    const opengl_input& input) {
  if (!apps->selected || !apps->selected->ok) return;
  auto app                  = apps->selected;
  app->glparams.window      = input.window_size;
  app->glparams.framebuffer = input.framebuffer_viewport;
  if (!app->glimage) app->glimage = make_glimage();
  if (app->glupdated) {
    set_glimage(app->glimage, app->display, false, false);
    app->glupdated = false;
  }
  update_imview(app->glparams.center, app->glparams.scale, app->display.size(),
      app->glparams.window, app->glparams.fit);
  draw_glimage(app->glimage, app->glparams);
}

void update(shared_ptr<opengl_window> win, shared_ptr<app_states> apps) {
  auto is_ready = [](const future<void>& result) -> bool {
    return result.valid() &&
           result.wait_for(chrono::microseconds(0)) == future_status::ready;
  };

  while (!apps->loading.empty()) {
    auto app = apps->loading.front();
    if (!is_ready(app->loader)) break;
    apps->loading.pop_front();
    try {
      app->loader.get();
      update_display(app);
      app->ok     = true;
      app->status = "ok";
    } catch (std::exception& e) {
      app->status = "";
      app->error  = e.what();
    }
  }
}

int run_app(int argc, const char* argv[]) {
  // prepare application
  auto apps      = make_shared<app_states>();
  auto filenames = vector<string>{};

  // command line options
  auto cli = CLI::App{"view images"};
  cli.add_option("images", filenames, "image filenames")->required();
  try {
    cli.parse(argc, argv);
  } catch (CLI::ParseError& e) {
    return cli.exit(e);
  }

  // loading images
  for (auto filename : filenames) load_image_async(apps, filename);

  // window
  auto win = make_glwindow({1280 + 320, 720}, "yimview", true);

  // callbacks
  set_update_glcallback(
      win, [apps](shared_ptr<opengl_window> win, const opengl_input& input) {
        update(win, apps);
      });
  set_draw_glcallback(
      win, [apps](shared_ptr<opengl_window> win, const opengl_input& input) {
        draw(win, apps, input);
      });
  set_widgets_glcallback(
      win, [apps](shared_ptr<opengl_window> win, const opengl_input& input) {
        draw_glwidgets(win, apps, input);
      });
  set_uiupdate_glcallback(
      win, [apps](shared_ptr<opengl_window> win, const opengl_input& input) {
        if (!apps->selected) return;
        auto app = apps->selected;
        // handle mouse
        if (input.mouse_left && !input.widgets_active) {
          app->glparams.center += input.mouse_pos - input.mouse_last;
        }
        if (input.mouse_right && !input.widgets_active) {
          app->glparams.scale *= powf(
              2, (input.mouse_pos.x - input.mouse_last.x) * 0.001f);
        }
      });
  set_drop_glcallback(
      win, [apps](shared_ptr<opengl_window> win, const vector<string>& paths,
               const opengl_input& input) {
        for (auto path : paths) load_image_async(apps, path);
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
    std::cerr << e.what() << "\n";
    return 1;
  }
}
