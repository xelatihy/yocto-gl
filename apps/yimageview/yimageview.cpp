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

#include <yocto/yocto_commonio.h>
#include <yocto/yocto_image.h>
#include <yocto_gui/yocto_gui.h>
using namespace yocto;

#include <atomic>
#include <deque>
#include <future>

#include "ext/filesystem.hpp"
namespace sfs = ghc::filesystem;

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
  ogl_image*       glimage   = new ogl_image{};
  ogl_image_params glparams  = {};
  bool             glupdated = true;

  // loading status
  std::atomic<bool> ok           = false;
  std::future<void> loader       = {};
  string            status       = "";
  string            error        = "";
  string            loader_error = "";

  // cleanup
  ~app_state() {
    if (glimage) delete glimage;
  }
};

// app states
struct app_states {
  // data
  vector<app_state*>     states   = {};
  app_state*             selected = nullptr;
  std::deque<app_state*> loading  = {};

  // default options
  float             exposure = 0;
  bool              filmic   = false;
  colorgrade_params params   = {};

  // cleanup
  ~app_states() {
    for (auto state : states) delete state;
  }
};

// compute min/max
void compute_stats(
    image_stats& stats, const image<vec4f>& img, bool linear_hdr) {
  auto max_histo = linear_hdr ? 8 : 1;
  stats.min      = vec4f{flt_max, flt_max, flt_max, flt_max};
  stats.max      = vec4f{flt_min, flt_min, flt_min, flt_min};
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

void update_display(app_state* app) {
  if (app->display.size() != app->source.size()) app->display = app->source;
  if (app->colorgrade) {
    colorgrade_image_mt(app->display, app->source, true, app->params);
  } else {
    tonemap_image_mt(app->display, app->source, app->exposure, app->filmic);
  }
  compute_stats(app->display_stats, app->display, false);
  app->glupdated = true;
}

// add a new image
void load_image_async(app_states* apps, const string& filename) {
  auto app      = apps->states.emplace_back(new app_state{});
  app->filename = filename;
  app->outname = sfs::path(filename).replace_extension(".display.png").string();
  app->name    = sfs::path(filename).filename();
  app->exposure = apps->exposure;
  app->filmic   = apps->filmic;
  app->params   = apps->params;
  app->status   = "loading";
  app->loader   = std::async(std::launch::async, [app]() {
    if (!load_image(app->filename, app->source, app->loader_error)) return;
    compute_stats(
        app->source_stats, app->source, is_hdr_filename(app->filename));
    if (app->colorgrade) {
      app->display = colorgrade_image(app->display, true, app->params);
    } else {
      app->display = tonemap_image(app->source, app->exposure, app->filmic);
    }
    compute_stats(app->display_stats, app->display, false);
  });
  apps->loading.push_back(app);
  if (!apps->selected) apps->selected = apps->states.front();
}

void draw_widgets(gui_window* win, app_states* apps, const gui_input& input) {
  static string load_path = "", save_path = "", error_message = "";
  if (draw_filedialog_button(win, "load", true, "load image", load_path, false,
          "./", "", "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    load_image_async(apps, load_path);
    load_path = "";
  }
  continue_line(win);
  if (draw_filedialog_button(win, "save", apps->selected && apps->selected->ok,
          "save image", save_path, true, sfs::path(save_path).parent_path(),
          sfs::path(save_path).filename(),
          "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    auto app     = apps->selected;
    app->outname = save_path;
    save_image(app->outname, app->display, app->error);
    save_path = "";
  }
  continue_line(win);
  if (draw_button(win, "close", (bool)apps->selected)) {
    if (apps->selected->loader.valid()) return;
    delete apps->selected;
    apps->states.erase(
        std::find(apps->states.begin(), apps->states.end(), apps->selected));
    apps->selected = apps->states.empty() ? nullptr : apps->states.front();
  }
  continue_line(win);
  if (draw_button(win, "quit")) {
    set_close(win, true);
  }
  draw_combobox(win, "image", apps->selected, apps->states, false);
  if (!apps->selected) return;
  auto app = apps->selected;
  if (app->status != "") draw_label(win, "status", app->status);
  if (app->error != "") draw_label(win, "error", app->error);
  if (!app->ok) return;
  if (begin_header(win, "tonemap")) {
    auto edited = 0;
    edited += draw_slider(win, "exposure", app->exposure, -5, 5);
    edited += draw_checkbox(win, "filmic", app->filmic);
    if (edited) update_display(app);
    end_header(win);
  }
  if (begin_header(win, "colorgrade")) {
    auto& params = app->params;
    auto  edited = 0;
    edited += draw_checkbox(win, "apply colorgrade", app->colorgrade);
    edited += draw_slider(win, "exposure", params.exposure, -5, 5);
    edited += draw_coloredit(win, "tint", params.tint);
    edited += draw_slider(win, "lincontrast", params.lincontrast, 0, 1);
    edited += draw_slider(win, "logcontrast", params.logcontrast, 0, 1);
    edited += draw_slider(win, "linsaturation", params.linsaturation, 0, 1);
    edited += draw_checkbox(win, "filmic", params.filmic);
    continue_line(win);
    edited += draw_checkbox(win, "srgb", params.srgb);
    continue_line(win);
    if (draw_button(win, "auto wb")) {
      auto wb     = 1 / xyz(app->source_stats.average);
      params.tint = wb / max(wb);
      edited += 1;
    }
    edited += draw_slider(win, "contrast", params.contrast, 0, 1);
    edited += draw_slider(win, "saturation", params.saturation, 0, 1);
    edited += draw_slider(win, "shadows", params.shadows, 0, 1);
    edited += draw_slider(win, "midtones", params.midtones, 0, 1);
    edited += draw_slider(win, "highlights", params.highlights, 0, 1);
    edited += draw_coloredit(win, "shadows color", params.shadows_color);
    edited += draw_coloredit(win, "midtones color", params.midtones_color);
    edited += draw_coloredit(win, "highlights color", params.highlights_color);
    if (edited) update_display(app);
    end_header(win);
  }
  if (begin_header(win, "inspect")) {
    draw_label(win, "image", sfs::path(app->filename).filename());
    draw_label(win, "filename", app->filename);
    draw_label(win, "outname", app->outname);
    draw_label(win, "image",
        std::to_string(app->source.size().x) + " x " +
            std::to_string(app->source.size().y));
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
    draw_dragger(win, "image min", app->source_stats.min);
    draw_dragger(win, "image max", app->source_stats.max);
    draw_dragger(win, "image avg", app->source_stats.average);
    draw_histogram(win, "image histo", app->source_stats.histogram);
    draw_dragger(win, "display min", app->display_stats.min);
    draw_dragger(win, "display max", app->display_stats.max);
    draw_dragger(win, "display avg", app->display_stats.average);
    draw_histogram(win, "display histo", app->display_stats.histogram);
    end_header(win);
  }
}

void draw(gui_window* win, app_states* apps, const gui_input& input) {
  if (!apps->selected || !apps->selected->ok) return;
  auto app                  = apps->selected;
  app->glparams.window      = input.window_size;
  app->glparams.framebuffer = input.framebuffer_viewport;
  if (!is_initialized(app->glimage)) init_image(app->glimage);
  if (app->glupdated) {
    set_image(app->glimage, app->display, false, false);
    app->glupdated = false;
  }
  update_imview(app->glparams.center, app->glparams.scale, app->display.size(),
      app->glparams.window, app->glparams.fit);
  draw_image(app->glimage, app->glparams);
}

void update(gui_window* win, app_states* apps) {
  auto is_ready = [](const std::future<void>& result) -> bool {
    return result.valid() && result.wait_for(std::chrono::microseconds(0)) ==
                                 std::future_status::ready;
  };

  while (!apps->loading.empty()) {
    auto app = apps->loading.front();
    if (!is_ready(app->loader)) break;
    apps->loading.pop_front();
    app->loader.get();
    if (app->loader_error.empty()) {
      update_display(app);
      app->ok     = true;
      app->status = "ok";
    } else {
      app->status = "";
      app->error  = app->loader_error;
    }
  }
}

int main(int argc, const char* argv[]) {
  // prepare application
  auto apps_guard = std::make_unique<app_states>();
  auto apps       = apps_guard.get();
  auto filenames  = vector<string>{};

  // command line options
  auto cli = make_cli("yimgview", "view images");
  add_option(cli, "images", filenames, "image filenames", true);
  parse_cli(cli, argc, argv);

  // loading images
  for (auto filename : filenames) load_image_async(apps, filename);

  // callbacks
  auto callbacks     = gui_callbacks{};
  callbacks.clear_cb = [apps](gui_window* win, const gui_input& input) {
    for (auto app : apps->states) clear_image(app->glimage);
  };
  callbacks.update_cb = [apps](gui_window* win, const gui_input& input) {
    update(win, apps);
  };
  callbacks.draw_cb = [apps](gui_window* win, const gui_input& input) {
    draw(win, apps, input);
  };
  callbacks.widgets_cb = [apps](gui_window* win, const gui_input& input) {
    draw_widgets(win, apps, input);
  };
  callbacks.uiupdate_cb = [apps](gui_window* win, const gui_input& input) {
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
  };
  callbacks.drop_cb = [apps](gui_window* win, const vector<string>& paths,
                          const gui_input& input) {
    for (auto path : paths) load_image_async(apps, path);
  };

  // run ui
  run_ui({1280 + 320, 720}, "yimview", callbacks);

  // done
  return 0;
}
