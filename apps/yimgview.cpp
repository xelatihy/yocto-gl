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
#include <list>

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
  opengl_image        glimage   = {};
  draw_glimage_params glparams  = {};
  bool                glupdated = true;

  // error
  string error = "";
};

struct app_states {
  // data
  std::list<app_state>         states   = {};
  int                          selected = -1;
  std::list<app_state>         loading  = {};
  std::list<std::future<bool>> loaders  = {};

  // get image
  app_state& get_selected() {
    auto it = states.begin();
    std::advance(it, selected);
    return *it;
  }
  const app_state& get_selected() const {
    auto it = states.begin();
    std::advance(it, selected);
    return *it;
  }

  // default options
  float             exposure = 0;
  bool              filmic   = false;
  colorgrade_params params   = {};
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
  if (app.display.size() != app.source.size()) app.display = app.source;
  parallel_for(app.source.size(), [&app](const vec2i& ij) {
    if (app.colorgrade) {
      app.display[ij] = colorgrade(app.source[ij], true, app.params);
    } else {
      app.display[ij] = tonemap(app.source[ij], app.exposure, app.filmic);
    }
  });
  compute_stats(app.display_stats, app.display, false);
  app.glupdated = true;
}

// add a new image
void load_image_async(app_states& apps, const string& filename) {
  auto& app     = apps.loading.emplace_back();
  app.filename  = filename;
  app.outname   = replace_extension(filename, ".display.png");
  app.name      = get_filename(filename);
  app.exposure  = apps.exposure;
  app.filmic    = apps.filmic;
  app.params    = apps.params;
  apps.selected = (int)apps.states.size() - 1;
  apps.loaders.push_back(std::async(std::launch::async, [&app]() -> bool {
    if (!load_image(app.filename, app.source)) {
      app.error = "cannot load " + app.filename;
      return false;
    }
    compute_stats(app.source_stats, app.source, is_hdr_filename(app.filename));
    if (app.colorgrade) {
      app.display = colorgrade_image(app.display, true, app.params);
    } else {
      app.display = tonemap_image(app.source, app.exposure, app.filmic);
    }
    compute_stats(app.display_stats, app.display, false);
    return true;
  }));
}

void draw_glwidgets(const opengl_window& win) {
  static string load_path = "", save_path = "", error_message = "";
  auto&         apps     = *(app_states*)get_gluser_pointer(win);
  auto          image_ok = !apps.states.empty() && apps.selected >= 0;
  if (!begin_glwidgets_window(win, "yimview")) return;
  draw_glmessages(win);
  if (draw_glfiledialog_button(win, "load", true, "load image", load_path,
          false, "./", "", "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    load_image_async(apps, load_path);
    load_path = "";
  }
  continue_glline(win);
  if (draw_glfiledialog_button(win, "save", image_ok, "save image", save_path,
          true, get_dirname(save_path), get_filename(save_path),
          "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    auto& app   = apps.get_selected();
    app.outname = save_path;
    try {
      save_image(app.outname, app.display);
    } catch (std::exception& e) {
      push_glmessage("cannot save " + app.outname);
      log_glinfo(win, "cannot save " + app.outname);
      log_glinfo(win, e.what());
    }
    save_path = "";
  }
  continue_glline(win);
  if (draw_glbutton(win, "close", image_ok)) {
    auto it = apps.states.begin();
    std::advance(it, apps.selected);
    apps.states.erase(it);
    apps.selected = apps.states.empty() ? -1 : 0;
  }
  continue_glline(win);
  if (draw_glbutton(win, "quit")) {
    set_glwindow_close(win, true);
  }
  draw_glcombobox(
      win, "image", apps.selected, (int)apps.states.size(),
      [&apps](int idx) {
        auto it = apps.states.begin();
        std::advance(it, idx);
        return it->name.c_str();
      },
      false);
  if (image_ok && begin_glheader(win, "tonemap")) {
    auto& app    = apps.get_selected();
    auto  edited = 0;
    edited += draw_glslider(win, "exposure", app.exposure, -5, 5);
    edited += draw_glcheckbox(win, "filmic", app.filmic);
    if (edited) update_display(app);
    end_glheader(win);
  }
  if (image_ok && begin_glheader(win, "colorgrade")) {
    auto& app    = apps.get_selected();
    auto& params = app.params;
    auto  edited = 0;
    edited += draw_glcheckbox(win, "apply colorgrade", app.colorgrade);
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
      auto wb     = 1 / xyz(app.source_stats.average);
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
  if (image_ok && begin_glheader(win, "inspect")) {
    auto& app = apps.get_selected();
    draw_gllabel(win, "image", get_filename(app.filename));
    draw_gllabel(win, "filename", app.filename);
    draw_gllabel(win, "outname", app.outname);
    draw_gllabel(win, "image",
        std::to_string(app.source.size().x) + " x " +
            std::to_string(app.source.size().y));
    draw_glslider(win, "zoom", app.glparams.scale, 0.1, 10);
    draw_glcheckbox(win, "fit", app.glparams.fit);
    auto mouse_pos = get_glmouse_pos(win);
    auto ij        = get_image_coords(
        mouse_pos, app.glparams.center, app.glparams.scale, app.source.size());
    draw_gldragger(win, "mouse", ij);
    auto img_pixel = zero4f, display_pixel = zero4f;
    if (ij.x >= 0 && ij.x < app.source.size().x && ij.y >= 0 &&
        ij.y < app.source.size().y) {
      img_pixel     = app.source[{ij.x, ij.y}];
      display_pixel = app.display[{ij.x, ij.y}];
    }
    draw_glcoloredit(win, "image", img_pixel);
    draw_gldragger(win, "display", display_pixel);
    draw_gldragger(win, "image min", app.source_stats.min);
    draw_gldragger(win, "image max", app.source_stats.max);
    draw_gldragger(win, "image avg", app.source_stats.average);
    draw_glhistogram(win, "image histo", app.source_stats.histogram);
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
  auto& apps = *(app_states*)get_gluser_pointer(win);
  if (!apps.states.empty() && apps.selected >= 0) {
    auto& app                = apps.get_selected();
    app.glparams.window      = get_glwindow_size(win);
    app.glparams.framebuffer = get_glframebuffer_viewport(win);
    if (!app.glimage || app.glupdated) {
      update_glimage(app.glimage, app.display, false, false);
      app.glupdated = false;
    }
    update_imview(app.glparams.center, app.glparams.scale, app.display.size(),
        app.glparams.window, app.glparams.fit);
    draw_glimage(app.glimage, app.glparams);
  }
  begin_glwidgets(win);
  draw_glwidgets(win);
  end_glwidgets(win);
  swap_glbuffers(win);
}

void update(const opengl_window& win, app_states& app) {
  auto is_ready = [](const std::future<bool>& result) -> bool {
    return result.valid() && result.wait_for(std::chrono::microseconds(0)) ==
                                 std::future_status::ready;
  };

  while (!app.loaders.empty() && is_ready(app.loaders.front())) {
    if (!app.loaders.front().get()) {
      push_glmessage(win, "cannot load image " + app.loading.front().filename);
      log_glinfo(win, "cannot load image " + app.loading.front().filename);
      log_glinfo(win, app.loading.front().error);
      break;
    }
    app.states.splice(app.states.end(), app.loading, app.loading.begin());
    app.loaders.pop_front();
    update_display(app.states.back());
    if (app.selected < 0) app.selected = (int)app.states.size() - 1;
  }
}

void run_ui(app_states& apps) {
  // window
  auto win = opengl_window();
  init_glwindow(win, {1280 + 320, 720}, "yimview", &apps, draw);
  set_drop_glcallback(
      win, [](const opengl_window& win, const vector<string>& paths) {
        auto& app = *(app_states*)get_gluser_pointer(win);
        for (auto path : paths) load_image_async(app, path);
      });

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
      auto& img = apps.get_selected();
      img.glparams.center += mouse_pos - last_pos;
    }
    if (mouse_right && !widgets_active) {
      auto& img = apps.get_selected();
      img.glparams.scale *= powf(2, (mouse_pos.x - last_pos.x) * 0.001f);
    }

    // update
    update(win, apps);

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
  auto apps      = app_states();
  auto filenames = vector<string>{};

  // command line options
  auto cli = make_cli("yimgview", "view images");
  add_cli_option(cli, "images", filenames, "image filenames", true);
  if (!parse_cli(cli, argc, argv)) exit(1);

  // loading images
  for (auto filename : filenames) load_image_async(apps, filename);

  // run ui
  run_ui(apps);

  // done
  return 0;
}
