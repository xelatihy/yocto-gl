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
#include <yocto/yocto_parallel.h>
#include <yocto_gui/yocto_imgui.h>
#include <yocto_gui/yocto_opengl.h>
#include <yocto_gui/yocto_window.h>
using namespace yocto;

#include <deque>

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

  // imgui
  gui_widgets widgets = {};

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
  auto num_pixels = (size_t)img.width() * (size_t)img.height();
  for (auto& v : stats.histogram) v /= num_pixels;
  stats.average /= num_pixels;
}

void update_display(app_state* app) {
  if (app->display.imsize() != app->source.imsize()) app->display = app->source;
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
  app->outname  = replace_extension(filename, ".display.png");
  app->name     = path_filename(filename);
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

void draw_widgets(app_states* apps, const gui_input& input) {
  auto widgets = &apps->widgets;
  begin_imgui(widgets, "yimageview");

  static string load_path = "", save_path = "", error_message = "";
  if (draw_filedialog_button(widgets, "load", true, "load image", load_path,
          false, "./", "", "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    load_image_async(apps, load_path);
    load_path = "";
  }
  continue_line(widgets);
  if (draw_filedialog_button(widgets, "save",
          apps->selected && apps->selected->ok, "save image", save_path, true,
          path_dirname(save_path), path_filename(save_path),
          "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    auto app     = apps->selected;
    app->outname = save_path;
    save_image(app->outname, app->display, app->error);
    save_path = "";
  }
  continue_line(widgets);
  if (draw_button(widgets, "close", (bool)apps->selected)) {
    if (apps->selected->loader.valid()) {
      end_imgui(widgets);
      return;
    }
    delete apps->selected;
    apps->states.erase(
        std::find(apps->states.begin(), apps->states.end(), apps->selected));
    apps->selected = apps->states.empty() ? nullptr : apps->states.front();
  }
  continue_line(widgets);
  if (draw_button(widgets, "quit")) {
    set_close(widgets->window, true);
  }
  draw_combobox(widgets, "image", apps->selected, apps->states, false);
  if (!apps->selected) {
    end_imgui(widgets);
    return;
  }
  auto app = apps->selected;
  if (app->status != "") draw_label(widgets, "status", app->status);
  if (app->error != "") draw_label(widgets, "error", app->error);
  if (!app->ok) {
    end_imgui(widgets);
    return;
  }
  if (begin_header(widgets, "tonemap")) {
    auto edited = 0;
    edited += draw_slider(widgets, "exposure", app->exposure, -5, 5);
    edited += draw_checkbox(widgets, "filmic", app->filmic);
    if (edited) update_display(app);
    end_header(widgets);
  }
  if (begin_header(widgets, "colorgrade")) {
    auto& params = app->params;
    auto  edited = 0;
    edited += draw_checkbox(widgets, "apply colorgrade", app->colorgrade);
    edited += draw_slider(widgets, "exposure", params.exposure, -5, 5);
    edited += draw_coloredit(widgets, "tint", params.tint);
    edited += draw_slider(widgets, "lincontrast", params.lincontrast, 0, 1);
    edited += draw_slider(widgets, "logcontrast", params.logcontrast, 0, 1);
    edited += draw_slider(widgets, "linsaturation", params.linsaturation, 0, 1);
    edited += draw_checkbox(widgets, "filmic", params.filmic);
    continue_line(widgets);
    edited += draw_checkbox(widgets, "srgb", params.srgb);
    continue_line(widgets);
    if (draw_button(widgets, "auto wb")) {
      auto wb     = 1 / xyz(app->source_stats.average);
      params.tint = wb / max(wb);
      edited += 1;
    }
    edited += draw_slider(widgets, "contrast", params.contrast, 0, 1);
    edited += draw_slider(widgets, "saturation", params.saturation, 0, 1);
    edited += draw_slider(widgets, "shadows", params.shadows, 0, 1);
    edited += draw_slider(widgets, "midtones", params.midtones, 0, 1);
    edited += draw_slider(widgets, "highlights", params.highlights, 0, 1);
    edited += draw_coloredit(widgets, "shadows color", params.shadows_color);
    edited += draw_coloredit(widgets, "midtones color", params.midtones_color);
    edited += draw_coloredit(
        widgets, "highlights color", params.highlights_color);
    if (edited) update_display(app);
    end_header(widgets);
  }
  if (begin_header(widgets, "inspect")) {
    draw_label(widgets, "image", app->name);
    draw_label(widgets, "filename", app->filename);
    draw_label(widgets, "outname", app->outname);
    draw_label(widgets, "image",
        std::to_string(app->source.width()) + " x " +
            std::to_string(app->source.height()));
    draw_slider(widgets, "zoom", app->glparams.scale, 0.1, 10);
    draw_checkbox(widgets, "fit", app->glparams.fit);
    auto ij = image_coords(input.mouse_pos, app->glparams.center,
        app->glparams.scale, app->source.imsize());
    draw_dragger(widgets, "mouse", ij);
    auto img_pixel = zero4f, display_pixel = zero4f;
    if (ij.x >= 0 && ij.x < app->source.width() && ij.y >= 0 &&
        ij.y < app->source.height()) {
      img_pixel     = app->source[{ij.x, ij.y}];
      display_pixel = app->display[{ij.x, ij.y}];
    }
    draw_coloredit(widgets, "image", img_pixel);
    draw_dragger(widgets, "display", display_pixel);
    draw_dragger(widgets, "image min", app->source_stats.min);
    draw_dragger(widgets, "image max", app->source_stats.max);
    draw_dragger(widgets, "image avg", app->source_stats.average);
    draw_histogram(widgets, "image histo", app->source_stats.histogram);
    draw_dragger(widgets, "display min", app->display_stats.min);
    draw_dragger(widgets, "display max", app->display_stats.max);
    draw_dragger(widgets, "display avg", app->display_stats.average);
    draw_histogram(widgets, "display histo", app->display_stats.histogram);
    end_header(widgets);
  }
  end_imgui(widgets);
}

void draw_image(app_states* apps, const gui_input& input) {
  if (!apps->selected || !apps->selected->ok) return;
  auto app                  = apps->selected;
  app->glparams.window      = input.window_size;
  app->glparams.framebuffer = input.framebuffer_viewport;
  if (!is_initialized(app->glimage)) init_image(app->glimage);
  if (app->glupdated) {
    set_image(app->glimage, app->display, false, false);
    app->glupdated = false;
  }
  std::tie(app->glparams.center, app->glparams.scale) = camera_imview(
      app->glparams.center, app->glparams.scale, app->display.imsize(),
      app->glparams.window, app->glparams.fit);
  draw_image(app->glimage, app->glparams);
}

void update(app_states* apps) {
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

void update_view(app_states* apps, const gui_input& input) {
  if (is_active(&apps->widgets)) return;

  if (!apps->selected) return;
  auto app = apps->selected;
  // handle mouse
  if (input.mouse_left) {
    app->glparams.center += input.mouse_pos - input.mouse_last;
  }
  if (input.mouse_right) {
    app->glparams.scale *= powf(
        2, (input.mouse_pos.x - input.mouse_last.x) * 0.001f);
  }
}

void drop(app_states* apps, const gui_input& input) {
  if (input.dropped.size()) {
    for (auto& path : input.dropped) load_image_async(apps, path);
  }
}

void update_app(const gui_input& input, void* data) {
  auto apps = (app_states*)data;
  update_view(apps, input);
  update(apps);
  draw_image(apps, input);
  draw_widgets(apps, input);
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

  auto window = new gui_window{};
  init_window(window, {1280 + 320, 720}, "yimageviews", true);
  window->user_data = apps;
  apps->widgets     = create_imgui(window);

  // run ui
  run_ui(window, update_app);

  // clear
  for (auto app : apps->states) clear_image(app->glimage);

  // done
  return 0;
}
