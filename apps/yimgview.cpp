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

#include <list>
using std::list;

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
  tonemap_params    tonemap_prms     = {};
  colorgrade_params colorgrade_prms  = {};
  bool              apply_colorgrade = false;

  // viewing properties
  opengl_image        gl_image  = {};
  draw_glimage_params draw_prms = {};

  // rendering properties
  int                  render_region  = 0;
  vector<image_region> render_regions = {};
  bool                 render_stats   = false;
};

struct app_states {
  // data
  std::list<app_state>         states       = {};
  int                          selected     = -1;
  std::list<app_state>         loading      = {};
  std::list<std::future<void>> loaders = {};

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
  tonemap_params    tonemap_prms    = {};
  colorgrade_params colorgrade_prms = {};
};

// reset display
void reset_display(app_state& app) {
  if (app.display.size() != app.source.size()) app.display = app.source;
  app.render_region  = 0;
  app.render_regions = make_image_regions(app.source.size(), 256);
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

// add a new image
void load_image_async(app_states& apps, const string& filename) {
  auto& app           = apps.loading.emplace_back();
  app.filename        = filename;
  app.outname         = replace_extension(filename, ".display.png");
  app.name            = get_filename(filename);
  app.tonemap_prms    = app.tonemap_prms;
  app.colorgrade_prms = app.colorgrade_prms;
  apps.selected       = (int)apps.states.size() - 1;
  apps.loaders.push_back(run_async([&app]() {
    load_image(app.filename, app.source);
    compute_stats(app.source_stats, app.source, is_hdr_filename(app.filename));
    app.display = tonemap_image(app.source, app.tonemap_prms);
    if (app.apply_colorgrade)
      app.display = colorgrade_image(app.display, app.colorgrade_prms);
    compute_stats(app.display_stats, app.display, false);
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
    auto& params = app.tonemap_prms;
    auto  edited = 0;
    edited += draw_glslider(win, "exposure", params.exposure, -5, 5);
    edited += draw_glcoloredit(win, "tint", params.tint);
    edited += draw_glslider(win, "contrast", params.contrast, 0, 1);
    edited += draw_glslider(win, "logcontrast", params.logcontrast, 0, 1);
    edited += draw_glslider(win, "saturation", params.saturation, 0, 1);
    edited += draw_glcheckbox(win, "filmic", params.filmic);
    continue_glline(win);
    edited += draw_glcheckbox(win, "srgb", params.srgb);
    continue_glline(win);
    if (draw_glbutton(win, "auto wb")) {
      auto wb     = 1 / xyz(app.source_stats.average);
      params.tint = wb / max(wb);
      edited += 1;
    }
    if (edited) reset_display(app);
    end_glheader(win);
  }
  if (image_ok && begin_glheader(win, "colorgrade")) {
    auto& app    = apps.get_selected();
    auto& params = app.colorgrade_prms;
    auto  edited = 0;
    edited += draw_glcheckbox(win, "apply colorgrade", app.apply_colorgrade);
    edited += draw_glslider(win, "contrast", params.contrast, 0, 1);
    edited += draw_glslider(win, "ldr shadows", params.shadows, 0, 1);
    edited += draw_glslider(win, "ldr midtones", params.midtones, 0, 1);
    edited += draw_glslider(win, "highlights", params.highlights, 0, 1);
    edited += draw_glcoloredit(win, "shadows color", params.shadows_color);
    edited += draw_glcoloredit(win, "midtones color", params.midtones_color);
    edited += draw_glcoloredit(
        win, "highlights color", params.highlights_color);
    if (edited) reset_display(app);
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
    draw_glslider(win, "zoom", app.draw_prms.scale, 0.1, 10);
    draw_glcheckbox(win, "fit", app.draw_prms.fit);
    auto mouse_pos = get_glmouse_pos(win);
    auto ij        = get_image_coords(mouse_pos, app.draw_prms.center,
        app.draw_prms.scale, app.source.size());
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
    auto& app                 = apps.get_selected();
    app.draw_prms.window      = get_glwindow_size(win);
    app.draw_prms.framebuffer = get_glframebuffer_viewport(win);
    if (!app.gl_image) update_glimage(app.gl_image, app.display, false, false);
    update_imview(app.draw_prms.center, app.draw_prms.scale, app.display.size(),
        app.draw_prms.window, app.draw_prms.fit);
    draw_glimage(app.gl_image, app.draw_prms);
  }
  begin_glwidgets(win);
  draw_glwidgets(win);
  end_glwidgets(win);
  swap_glbuffers(win);
}

void update(const opengl_window& win, app_states& app) {
  while (!app.loaders.empty() && is_ready(app.loaders.front())) {
    try {
      app.loaders.front().get();
    } catch (const std::exception& e) {
      push_glmessage(win, "cannot load image " + app.loading.front().filename);
      log_glinfo(win, "cannot load image " + app.loading.front().filename);
      log_glinfo(win, e.what());
      break;
    }
    app.states.splice(app.states.end(), app.loading, app.loading.begin());
    app.loaders.pop_front();
    reset_display(app.states.back());
    if (app.selected < 0) app.selected = (int)app.states.size() - 1;
  }
  for (auto& app : app.states) {
    if (app.render_region < app.render_regions.size()) {
      auto num_regions = min(12, app.render_regions.size() - app.render_region);
      parallel_for(num_regions, [&app](int idx) {
        auto& region = app.render_regions[app.render_region + idx];
        tonemap_region(app.display, app.source,
            app.render_regions[app.render_region + idx], app.tonemap_prms);
        if (app.apply_colorgrade) {
          colorgrade_region(
              app.display, app.display, region, app.colorgrade_prms);
        }
      });
      if (!app.gl_image) {
        update_glimage(app.gl_image, app.display, false, false);
      } else {
        for (auto idx = 0; idx < num_regions; idx++)
          update_glimage_region(app.gl_image, app.display,
              app.render_regions[app.render_region + idx]);
      }
      app.render_region += num_regions;
    } else if (app.render_stats) {
      compute_stats(app.display_stats, app.display, false);
    }
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
      img.draw_prms.center += mouse_pos - last_pos;
    }
    if (mouse_right && !widgets_active) {
      auto& img = apps.get_selected();
      img.draw_prms.scale *= powf(2, (mouse_pos.x - last_pos.x) * 0.001f);
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
