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

struct image_stats {
  vec4f         min       = zero4f;
  vec4f         max       = zero4f;
  vec4f         average   = zero4f;
  vector<vec3f> histogram = {};
};

struct app_state {
  // original data
  string name     = "image";
  string filename = "image.png";
  string outname  = "out.png";

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
  bool              colorgrade      = false;

  // viewing properties
  vec2f image_center = zero2f;
  float image_scale  = 1;
  bool  zoom_to_fit  = false;

  // computation
  deque<pair<string, int>> updates     = {};
  bool                     load_done   = false;
  string                   error       = "";
  std::future<void>        load_worker = {};
};

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

void update_display(app_state& app) {
  if (app.display.size() != app.img.size()) app.display = app.img;
  auto regions = make_regions(app.img.size(), 128);
  parallel_foreach(regions, [&app](const image_region& region) {
    tonemap(app.display, app.img, region, app.tonemap_prms);
    if (app.colorgrade) {
      colorgrade(app.display, app.display, region, app.colorgrade_prms);
    }
  });
  compute_stats(app.display_stats, app.display, false);
}

void load_image(app_state& app) {
  app.load_done = false;
  load_image(app.filename, app.img);
  compute_stats(app.image_stats, app.img, is_hdr_filename(app.filename));
  update_display(app);
  app.load_done = true;
}

void draw_glwidgets(const opengl_window& win) {
  auto& app = *(app_state*)get_gluser_pointer(win);
  if (!begin_glwidgets_window(win, "yimview")) return;
  if(!app.error.empty()) {
    if(!is_glmodal_open(win, "error")) {
      open_glmodal(win, "error");
    } else if (!draw_glmessage(win, "error", app.error)) {
      app.error = "";
    }
  }
  if (begin_glheader(win, "yimview")) {
    draw_gllabel(win, "image",
        get_filename(app.filename) + " @ " + to_string(app.img.size().x) +
            " x " + to_string(app.img.size().y));
    draw_gllabel(win, "filename", app.filename);
    draw_gllabel(win, "outname", app.outname);
    if (draw_glbutton(win, "save")) {
      // TODO: save
      // app.images[app.selected].outname;
    }
    continue_glline(win);
    if (draw_glbutton(win, "quit")) set_glwindow_close(win, true);
    end_glheader(win);
  }
  if (app.load_done && begin_glheader(win, "tonemap")) {
    auto& tonemap = app.tonemap_prms;
    auto  edited  = 0;
    edited += (int)draw_glslider(win, "exposure", tonemap.exposure, -5, 5);
    edited += (int)draw_glcoloredit(win, "tint", tonemap.tint);
    edited += (int)draw_glslider(win, "contrast", tonemap.contrast, 0, 1);
    edited += (int)draw_glslider(win, "logcontrast", tonemap.logcontrast, 0, 1);
    edited += (int)draw_glslider(win, "saturation", tonemap.saturation, 0, 1);
    edited += (int)draw_glcheckbox(win, "filmic", tonemap.filmic);
    continue_glline(win);
    edited += (int)draw_glcheckbox(win, "srgb", tonemap.srgb);
    continue_glline(win);
    if (draw_glbutton(win, "auto wb")) {
      auto wb      = 1 / xyz(app.image_stats.average);
      tonemap.tint = wb / max(wb);
      edited += 1;
    }
    if (edited) app.updates.push_back({"tonemap", -1});
    end_glheader(win);
  }
  if (app.load_done && begin_glheader(win, "colorgrade")) {
    auto& colorgrade = app.colorgrade_prms;
    auto  edited     = 0;
    edited += (int)draw_glcheckbox(win, "enable colorgrade", app.colorgrade);
    edited += (int)draw_glslider(win, "contrast", colorgrade.contrast, 0, 1);
    edited += (int)draw_glslider(win, "ldr shadows", colorgrade.shadows, 0, 1);
    edited += (int)draw_glslider(
        win, "ldr midtones", colorgrade.midtones, 0, 1);
    edited += (int)draw_glslider(
        win, "highlights", colorgrade.highlights, 0, 1);
    edited += (int)draw_glcoloredit(
        win, "shadows color", colorgrade.shadows_color);
    edited += (int)draw_glcoloredit(
        win, "midtones color", colorgrade.midtones_color);
    edited += (int)draw_glcoloredit(
        win, "highlights color", colorgrade.highlights_color);
    if (edited) app.updates.push_back({"colorgrade", -1});
    end_glheader(win);
  }
  if (app.load_done && begin_glheader(win, "inspect")) {
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
    draw_gldragger(win, "image min", app.image_stats.min);
    draw_gldragger(win, "image max", app.image_stats.max);
    draw_gldragger(win, "image avg", app.image_stats.average);
    draw_glhistogram(win, "image histo", app.image_stats.histogram);
    draw_gldragger(win, "display min", app.display_stats.min);
    draw_gldragger(win, "display max", app.display_stats.max);
    draw_gldragger(win, "display avg", app.display_stats.average);
    draw_glhistogram(win, "display histo", app.display_stats.histogram);
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
  if (app.gl_txt) {
    update_imview(app.image_center, app.image_scale, app.display.size(),
        win_size, app.zoom_to_fit);
    draw_glimage_background(
        app.gl_txt, win_size.x, win_size.y, app.image_center, app.image_scale);
    set_glblending(true);
    draw_glimage(
        app.gl_txt, win_size.x, win_size.y, app.image_center, app.image_scale);
    set_glblending(false);
  }
  begin_glwidgets(win);
  draw_glwidgets(win);
  end_glwidgets(win);
  swap_glbuffers(win);
}

void update(const opengl_window& win, app_state& app) {
  if (is_valid(app.load_worker)) {
    if(!is_ready(app.load_worker)) return;
    try {
      app.load_worker.get();
      app.load_done = true;
    } catch(const std::exception& e) {
      app.error = "cannot load image "s + e.what();
      log_glinfo(win, "cannot load image " + app.filename);
      log_glinfo(win, e.what());
    }
  }
  if (!app.load_done) return;
  if (app.gl_txt.size != app.img.size()) {
    init_gltexture(app.gl_txt, app.display, false, false, false);
  }
  if (app.updates.empty()) return;
  update_display(app);
  update_gltexture(app.gl_txt, app.display, false);
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
  add_cli_option(cli, "--output,-o", app.outname, "image output");
  add_cli_option(cli, "image", app.filename, "image filename", true);
  if (!parse_cli(cli, argc, argv)) exit(1);

  // load image
  app.load_worker = run_async([&app]() { load_image(app); });

  // run ui
  run_ui(app);

  // done
  return 0;
}
