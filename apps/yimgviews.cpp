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

#include "ext/filesystem.hpp"
namespace fs = ghc::filesystem;

struct image_stats {
  vec4f         min       = zero4f;
  vec4f         max       = zero4f;
  vec4f         average   = zero4f;
  vector<vec3f> histogram = {};
};

enum struct app_task_type { none, load, save, display, close };

struct app_state {
  // original data
  string name     = "";
  string filename = "";
  string outname  = "";

  // image data
  image<vec4f> img = {};

  // diplay image
  image<vec4f>   display = {};
  opengl_texture gl_txt  = {};

  // image stats
  image_stats image_stats, display_stats;

  // tonemapping values
  tonemap_params    tonemap_prms    = {};
  colorgrade_params colorgrade_prms = {};

  // computation futures
  bool                      load_done = false, display_done = false;
  std::deque<app_task_type> task_queue;

  // viewing properties
  vec2f image_center = zero2f;
  float image_scale  = 1;
  bool  zoom_to_fit  = false;
};

// compute min/max
void compute_image_stats(
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

void update_app_display(const string& filename, const image<vec4f>& img,
    image<vec4f>& display, image_stats& stats,
    const tonemap_params&    tonemap_prms,
    const colorgrade_params& colorgrade_prms, std::atomic<bool>* cancel,
    std::deque<image_region>& queue, std::mutex& queuem) {
  auto                regions  = make_regions(img.size(), 128);
  auto                futures  = vector<std::future<void>>{};
  auto                nthreads = std::thread::hardware_concurrency();
  std::atomic<size_t> next_idx(0);
  for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
    futures.emplace_back(std::async(std::launch::async,
        [&next_idx, cancel, &regions, &queue, &queuem, &display, &img,
            tonemap_prms, colorgrade_prms]() {
          while (true) {
            if (cancel && *cancel) break;
            auto idx = next_idx.fetch_add(1);
            if (idx >= regions.size()) break;
            auto region = regions[idx];
            tonemap(display, img, region, tonemap_prms);
            if (colorgrade_prms != colorgrade_params{}) {
              colorgrade(display, display, region, colorgrade_prms);
            }
            {
              std::lock_guard guard{queuem};
              queue.push_back(region);
            }
          }
        }));
  }
  for (auto& f : futures) f.get();
  compute_image_stats(stats, display, false);
}

// add a new image
void add_new_image(app_state& app, const string& filename) {
  app.filename        = filename;
  app.outname         = fs::path(filename).replace_extension(".display.png");
  app.name            = fs::path(filename).filename();
  app.tonemap_prms    = app.tonemap_prms;
  app.colorgrade_prms = app.colorgrade_prms;
  app.load_done       = false;
  app.display_done    = false;
  app.task_queue.emplace_back(app_task_type::load);
}

void draw_glwidgets(const opengl_window& win) {
  auto&         app = *(app_state*)get_gluser_pointer(win);
  if (!begin_glwidgets_window(win, "yimview")) return;
  if (draw_glbutton(win, "save", app.display_done)) {
    // TODO: save
    // app.images[app.selected].outname;
  }
  continue_glline(win);
  if (draw_glbutton(win, "quit")) set_glwindow_close(win, true);
  if (begin_glheader(win, "tonemap")) {
    auto options = app.tonemap_prms;
    draw_glslider(win, "exposure", options.exposure, -5, 5);
    draw_glcoloredit(win, "tint", options.tint);
    draw_glslider(win, "contrast", options.contrast, 0, 1);
    draw_glslider(win, "logcontrast", options.logcontrast, 0, 1);
    draw_glslider(win, "saturation", options.saturation, 0, 1);
    draw_glcheckbox(win, "filmic", options.filmic);
    continue_glline(win);
    draw_glcheckbox(win, "srgb", options.srgb);
    continue_glline(win);
    if (draw_glbutton(win, "auto wb")) {
      auto wb      = 1 / xyz(app.image_stats.average);
      options.tint = wb / max(wb);
    }
    if (options != app.tonemap_prms) {
      app.tonemap_prms = options;
      if (app.load_done) app.task_queue.emplace_back(app_task_type::display);
    }
    end_glheader(win);
  }
  if (begin_glheader(win, "colorgrade")) {
    auto options = app.colorgrade_prms;
    draw_glslider(win, "contrast", options.contrast, 0, 1);
    draw_glslider(win, "ldr shadows", options.shadows, 0, 1);
    draw_glslider(win, "ldr midtones", options.midtones, 0, 1);
    draw_glslider(win, "highlights", options.highlights, 0, 1);
    draw_glcoloredit(win, "shadows color", options.shadows_color);
    draw_glcoloredit(win, "midtones color", options.midtones_color);
    draw_glcoloredit(win, "highlights color", options.highlights_color);
    if (options != app.colorgrade_prms) {
      app.colorgrade_prms = options;
      if (app.load_done) app.task_queue.emplace_back(app_task_type::display);
    }
    end_glheader(win);
  }
  if (begin_glheader(win, "inspect")) {
    draw_gllabel(win, "image", fs::path(app.filename).filename());
    draw_gllabel(win, "filename", app.filename);
    draw_gllabel(win, "outname", app.outname);
    draw_gllabel(win, "image", "%d x %d", app.img.size().x, app.img.size().y);
    draw_glslider(win, "zoom", app.image_scale, 0.1, 10);
    draw_glcheckbox(win, "zoom to fit", app.zoom_to_fit);
    auto mouse_pos = get_glmouse_pos(win);
    auto ij        = get_image_coords(
        mouse_pos, app.image_center, app.image_scale, app.img.size());
    draw_gldragger(win, "mouse", ij);
    auto img_pixel = zero4f, display_pixel = zero4f;
    if (ij.x >= 0 && ij.x < app.img.size().x && ij.y >= 0 &&
        ij.y < app.img.size().y) {
      img_pixel     = app.img[{ij.x, ij.y}];
      display_pixel = app.display[{ij.x, ij.y}];
    }
    draw_glcoloredit(win, "image", img_pixel);
    draw_gldragger(win, "display", display_pixel);
    auto img_stats = (app.load_done) ? app.image_stats : image_stats{};
    draw_gldragger(win, "image min", img_stats.min);
    draw_gldragger(win, "image max", img_stats.max);
    draw_gldragger(win, "image avg", img_stats.average);
    draw_glhistogram(win, "image histo", img_stats.histogram);
    auto display_stats = (app.load_done) ? app.display_stats : image_stats{};
    draw_gldragger(win, "display min", display_stats.min);
    draw_gldragger(win, "display max", display_stats.max);
    draw_gldragger(win, "display avg", display_stats.average);
    draw_glhistogram(win, "display histo", display_stats.histogram);
    end_glheader(win);
  }
  if (begin_glheader(win, "log")) {
    draw_gllog(win);
    end_glheader(win);
  }
}

void draw(const opengl_window& win) {
  auto& app      = *(app_state*)get_gluser_pointer(win);
  auto  win_size = get_glwindow_size(win);
  auto  fb_view  = get_glframebuffer_viewport(win);
  set_glviewport(fb_view);
  clear_glframebuffer(vec4f{0.15f, 0.15f, 0.15f, 1.0f});
  if (app.load_done && app.gl_txt) {
    update_imview(app.image_center, app.image_scale, app.display.size(),
        win_size, app.zoom_to_fit);
    draw_glimage_background(app.gl_txt, win_size.x, win_size.y,
        app.image_center, app.image_scale);
    set_glblending(true);
    draw_glimage(app.gl_txt, win_size.x, win_size.y, app.image_center,
        app.image_scale);
    set_glblending(false);
  }
  begin_glwidgets(win);
  draw_glwidgets(win);
  end_glwidgets(win);
  swap_glbuffers(win);
}

void update(const opengl_window& win, app_state& app) {
  // handle tasks
    if (app.task_queue.empty()) return;
    auto cmd = app.task_queue.front();
    if(cmd == app_task_type::close) return;
    app.task_queue.pop_back();
    switch (cmd) {
      case app_task_type::load: {
        log_glinfo(win, "start loading " + app.filename);
        app.load_done = false;
        app.img       = {};
        try {
          load_image(app.filename, app.img);
        } catch (...) {
          log_glinfo(win, "cannot load " + app.filename);
        }
        compute_image_stats(
            app.image_stats, app.img, is_hdr_filename(app.filename));
        app.load_done = true;
        app.name      = fs::path(app.filename).filename().string() + " [" +
                   std::to_string(app.img.size().x) + "x" +
                   std::to_string(app.img.size().y) + "]";
        app.display = app.img;
        log_glinfo(win, "done loading " + app.filename);
        init_gltexture(app.gl_txt, app.display, false, false, false);
        app.task_queue.emplace_back(app_task_type::display);
      } break;
      case app_task_type::save: {
        log_glinfo(win, "start saving " + app.outname);
        try {
          if (!is_hdr_filename(app.outname)) {
            auto ldr = image<vec4b>{};
            float_to_byte(ldr, app.display);
            save_imageb(app.outname, ldr);
          } else {
            auto aux = image<vec4f>{};
            srgb_to_rgb(aux, app.display);
            save_image(app.outname, aux);
          }
        } catch (...) {
          log_glinfo(win, "cannot save " + app.outname);
        }
        log_glinfo(win, "done saving " + app.outname);
      } break;
      case app_task_type::display: {
        auto                regions  = make_regions(app.img.size(), 128);
        auto                futures  = vector<std::future<void>>{};
        auto                nthreads = std::thread::hardware_concurrency();
        std::atomic<size_t> next_idx(0);
        for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
          futures.emplace_back(std::async(std::launch::async,
              [&next_idx, &regions, &app]() {
                while (true) {
                  auto idx = next_idx.fetch_add(1);
                  if (idx >= regions.size()) break;
                  auto region = regions[idx];
                  tonemap(app.display, app.img, region, app.tonemap_prms);
                  if (app.colorgrade_prms != colorgrade_params{}) {
                    colorgrade(app.display, app.display, region, app.colorgrade_prms);
                  }
                }
              }));
        }
        for (auto& f : futures) f.get();
        compute_image_stats(app.display_stats, app.display, false);
        update_gltexture(app.gl_txt, app.display, false);
        app.display_done = true;
      } break;
      default: break;
    }
}

void run_ui(app_state& app) {
  // window
  auto win = opengl_window();
  init_glwindow(win, {1280 + 320, 720}, "yimview", &app, draw);

  // init widgets
  init_glwidgets(win);

  // window values
  auto mouse_pos = zero2f, last_pos = zero2f;
  while (!should_glwindow_close(win)) {
    last_pos            = mouse_pos;
    mouse_pos           = get_glmouse_pos(win);
    auto mouse_left     = get_glmouse_left(win);
    auto mouse_right    = get_glmouse_right(win);
    auto widgets_active = get_glwidgets_active(win);

    // handle mouse
    if (mouse_left && !widgets_active) {
      app.image_center += mouse_pos - last_pos;
    }
    if (mouse_right && !widgets_active) {
      app.image_scale *= powf(2, (mouse_pos.x - last_pos.x) * 0.001f);
    }

    // update
    update(win, app);

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
  add_cli_option(cli, "images", filenames, "image filenames", true);
  if (!parse_cli(cli, argc, argv)) exit(1);

  // loading images
  for (auto filename : filenames) add_new_image(app, filename);

  // run ui
  run_ui(app);

  // done
  return 0;
}
