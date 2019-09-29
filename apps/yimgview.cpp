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
  vec2f          image_center = zero2f;
  float          image_scale  = 1;
  bool           zoom_to_fit  = false;
  opengl_texture gl_txt       = {};

  // rendering properties
  int                  render_region  = 0;
  vector<image_region> render_regions = {};
  bool                 render_stats   = false;
};

struct app_states {
  // data
  std::list<app_state>    states       = {};
  int                     selected     = -1;
  std::list<app_state>    loading      = {};
  std::list<future<void>> load_workers = {};

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
void reset_display(app_state& state) {
  if (state.display.size() != state.source.size()) state.display = state.source;
  state.render_region  = 0;
  state.render_regions = make_regions(state.source.size(), 256);
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
void load_image_async(app_states& app, const string& filename) {
  auto& state           = app.loading.emplace_back();
  state.filename        = filename;
  state.outname         = replace_extension(filename, ".display.png");
  state.name            = get_filename(filename);
  state.tonemap_prms    = app.tonemap_prms;
  state.colorgrade_prms = app.colorgrade_prms;
  app.selected          = (int)app.states.size() - 1;
  app.load_workers.push_back(run_async([&state]() {
    load_image(state.filename, state.source);
    compute_stats(
        state.source_stats, state.source, is_hdr_filename(state.filename));
    state.display = tonemap(state.source, state.tonemap_prms);
    if (state.apply_colorgrade)
      state.display = colorgrade(state.display, state.colorgrade_prms);
    compute_stats(state.display_stats, state.display, false);
  }));
}

void draw_glwidgets(const opengl_window& win) {
  static string load_path = "", save_path = "", error_message = "";
  auto&         app      = *(app_states*)get_gluser_pointer(win);
  auto          image_ok = !app.states.empty() && app.selected >= 0;
  if (!begin_glwidgets_window(win, "yimview")) return;
  draw_glmessages(win);
  if (draw_glfiledialog(win, "load image", load_path, false, "./", "",
          "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    load_image_async(app, load_path);
  }
  if (draw_glfiledialog(win, "save image", save_path, true,
          get_dirname(save_path), get_filename(save_path),
          "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    // app.images[app.selected].outname = save_path;
    // app.images[app.selected].task_queue.emplace_back(app_task_type::save);
    save_path = "";
    // TODO> implement save
  }
  if (draw_glbutton(win, "load")) {
    open_glmodal(win, "load image");
  }
  continue_glline(win);
  if (draw_glbutton(win, "save", image_ok)) {
    // TODO: implement save
    // save_path = app.images[app.selected].outname;
    // open_glmodal(win, "save image");
  }
  continue_glline(win);
  if (draw_glbutton(win, "close", image_ok)) {
    auto it = app.states.begin();
    std::advance(it, app.selected);
    app.states.erase(it);
    app.selected = app.states.empty() ? -1 : 0;
  }
  continue_glline(win);
  if (draw_glbutton(win, "quit")) {
    set_glwindow_close(win, true);
  }
  if (app.states.empty()) return;
  draw_glcombobox(
      win, "image", app.selected, (int)app.states.size(),
      [&app](int idx) {
        auto it = app.states.begin();
        std::advance(it, idx);
        return it->name.c_str();
      },
      false);
  if (image_ok && begin_glheader(win, "tonemap")) {
    auto& state  = app.get_selected();
    auto& params = state.tonemap_prms;
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
      auto wb     = 1 / xyz(state.source_stats.average);
      params.tint = wb / max(wb);
      edited += 1;
    }
    if (edited) reset_display(state);
    end_glheader(win);
  }
  if (image_ok && begin_glheader(win, "colorgrade")) {
    auto& state  = app.get_selected();
    auto& params = state.colorgrade_prms;
    auto  edited = 0;
    edited += draw_glcheckbox(win, "apply colorgrade", state.apply_colorgrade);
    edited += draw_glslider(win, "contrast", params.contrast, 0, 1);
    edited += draw_glslider(win, "ldr shadows", params.shadows, 0, 1);
    edited += draw_glslider(win, "ldr midtones", params.midtones, 0, 1);
    edited += draw_glslider(win, "highlights", params.highlights, 0, 1);
    edited += draw_glcoloredit(win, "shadows color", params.shadows_color);
    edited += draw_glcoloredit(win, "midtones color", params.midtones_color);
    edited += draw_glcoloredit(
        win, "highlights color", params.highlights_color);
    if (edited) reset_display(state);
    end_glheader(win);
  }
  if (image_ok && begin_glheader(win, "inspect")) {
    auto& state = app.get_selected();
    draw_gllabel(win, "image", get_filename(state.filename));
    draw_gllabel(win, "filename", state.filename);
    draw_gllabel(win, "outname", state.outname);
    draw_gllabel(
        win, "image", "%d x %d", state.source.size().x, state.source.size().y);
    draw_glslider(win, "zoom", state.image_scale, 0.1, 10);
    draw_glcheckbox(win, "zoom to fit", state.zoom_to_fit);
    auto mouse_pos = get_glmouse_pos(win);
    auto ij        = get_image_coords(
        mouse_pos, state.image_center, state.image_scale, state.source.size());
    draw_gldragger(win, "mouse", ij);
    auto img_pixel = zero4f, display_pixel = zero4f;
    if (ij.x >= 0 && ij.x < state.source.size().x && ij.y >= 0 &&
        ij.y < state.source.size().y) {
      img_pixel     = state.source[{ij.x, ij.y}];
      display_pixel = state.display[{ij.x, ij.y}];
    }
    draw_glcoloredit(win, "image", img_pixel);
    draw_gldragger(win, "display", display_pixel);
    draw_gldragger(win, "image min", state.source_stats.min);
    draw_gldragger(win, "image max", state.source_stats.max);
    draw_gldragger(win, "image avg", state.source_stats.average);
    draw_glhistogram(win, "image histo", state.source_stats.histogram);
    draw_gldragger(win, "display min", state.display_stats.min);
    draw_gldragger(win, "display max", state.display_stats.max);
    draw_gldragger(win, "display avg", state.display_stats.average);
    draw_glhistogram(win, "display histo", state.display_stats.histogram);
    end_glheader(win);
  }
  if (begin_glheader(win, "log")) {
    draw_gllog(win);
    end_glheader(win);
  }
}

void draw(const opengl_window& win) {
  auto& app      = *(app_states*)get_gluser_pointer(win);
  auto  win_size = get_glwindow_size(win);
  auto  fb_view  = get_glframebuffer_viewport(win);
  set_glviewport(fb_view);
  clear_glframebuffer(vec4f{0.15f, 0.15f, 0.15f, 1.0f});
  if (!app.states.empty() && app.selected >= 0) {
    auto& state = app.get_selected();
    if (!state.gl_txt || state.gl_txt.size != state.display.size())
      init_gltexture(state.gl_txt, state.display, false, false, false);
    update_imview(state.image_center, state.image_scale, state.display.size(),
        win_size, state.zoom_to_fit);
    draw_glimage_background(state.gl_txt, win_size.x, win_size.y,
        state.image_center, state.image_scale);
    set_glblending(true);
    draw_glimage(state.gl_txt, win_size.x, win_size.y, state.image_center,
        state.image_scale);
    set_glblending(false);
  }
  begin_glwidgets(win);
  draw_glwidgets(win);
  end_glwidgets(win);
  swap_glbuffers(win);
}

void update(const opengl_window& win, app_states& app) {
  while (app.load_workers.empty() && is_ready(app.load_workers.front())) {
    try {
      app.load_workers.front().get();
    } catch (const std::exception& e) {
      push_glmessage(win, "cannot load image " + app.loading.front().filename);
      log_glinfo(win, "cannot load image " + app.loading.front().filename);
      log_glinfo(win, e.what());
      break;
    }
    app.states.splice(app.states.end(), app.loading, app.loading.begin());
    reset_display(app.states.back());
    if (app.selected < 0) app.selected = (int)app.states.size() - 1;
  }
  for (auto& state : app.states) {
    if (state.render_region < state.render_regions.size()) {
      auto num_regions = min(
          12, state.render_regions.size() - state.render_region);
      parallel_for(num_regions, [&state](int idx) {
        auto& region = state.render_regions[state.render_region + idx];
        tonemap(state.display, state.source,
            state.render_regions[state.render_region + idx],
            state.tonemap_prms);
        if (state.apply_colorgrade) {
          colorgrade(
              state.display, state.display, region, state.colorgrade_prms);
        }
      });
      if (!state.gl_txt || state.gl_txt.size != state.display.size()) {
        init_gltexture(state.gl_txt, state.display, false, false, false);
      } else {
        for (auto idx = 0; idx < num_regions; idx++)
          update_gltexture_region(state.gl_txt, state.display,
              state.render_regions[state.render_region + idx], false);
      }
      state.render_region += num_regions;
    } else if (state.render_stats) {
      compute_stats(state.display_stats, state.display, false);
    }
  }
}

void run_ui(app_states& states) {
  // window
  auto win = opengl_window();
  init_glwindow(win, {1280 + 320, 720}, "yimview", &states, draw);
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
      auto& img = states.get_selected();
      img.image_center += mouse_pos - last_pos;
    }
    if (mouse_right && !widgets_active) {
      auto& img = states.get_selected();
      img.image_scale *= powf(2, (mouse_pos.x - last_pos.x) * 0.001f);
    }

    // update
    update(win, states);

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
  auto states       = app_states();
  auto filenames = vector<string>{};

  // command line options
  auto cli = make_cli("yimgview", "view images");
  add_cli_option(cli, "images", filenames, "image filenames", true);
  if (!parse_cli(cli, argc, argv)) exit(1);

  // loading images
  for (auto filename : filenames) load_image_async(states, filename);

  // run ui
  run_ui(states);

  // done
  return 0;
}
