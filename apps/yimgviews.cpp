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
#include <any>
#include <typeindex>
using std::any;
using std::any_cast;
using std::type_index;

struct image_stats {
  vec4f         min       = zero4f;
  vec4f         max       = zero4f;
  vec4f         average   = zero4f;
  vector<vec3f> histogram = {};
};

struct app_image {
  // original data
  string name     = "";
  string filename = "";
  string outname  = "";

  // image data
  image<vec4f> source = {};
  image<vec4f>   display = {};

  // image stats
  image_stats source_stats = {};
  image_stats display_stats = {};

  // tonemapping values
  tonemap_params    tonemap_prms    = {};
  colorgrade_params colorgrade_prms = {};
  bool apply_colorgrade = false;

  // computation futures
  bool load_done = false;
  bool display_done = false;
  future<void> load_worker = {};
  future<void> display_worker = {};
  string error = "";

  // viewing properties
  vec2f image_center = zero2f;
  float image_scale  = 1;
  bool  zoom_to_fit  = false;
  opengl_texture gl_txt  = {};
};

struct app_state {
  // data
  std::deque<app_image> images;
  int                   selected = -1;
  std::deque<string>    errors;

  // default options
  tonemap_params    tonemap_prms    = {};
  colorgrade_params colorgrade_prms = {};
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

void update_display(app_image& image) {
  if (image.display.size() != image.source.size()) image.display = image.source;
  auto regions = make_regions(image.source.size(), 128);
  parallel_foreach(regions, [&image](const image_region& region) {
    tonemap(image.display, image.source, region, image.tonemap_prms);
    if (image.apply_colorgrade) {
      colorgrade(image.display, image.display, region, image.colorgrade_prms);
    }
  });
  compute_stats(image.display_stats, image.display, false);
}

void update_texture(app_image& image) {
  if(!image.load_done) return;
  if(!image.gl_txt) {
    init_gltexture(image.gl_txt, image.display, false, false, false);
  }
  update_gltexture(image.gl_txt, image.display, false);
}

// add a new image
void load_image_async(app_state& app, const string& filename) {
  auto& image           = app.images.emplace_back();
  image.filename        = filename;
  image.outname         = replace_extension(filename, ".display.png");
  image.name            = get_filename(filename);
  image.tonemap_prms    = app.tonemap_prms;
  image.colorgrade_prms = app.colorgrade_prms;
  image.load_done       = false;
  image.display_done    = false;
  app.selected = (int)app.images.size() - 1;
  image.load_worker = run_async([&image]() { 
    load_image(image.filename, image.source);
    compute_stats(image.source_stats, image.source, is_hdr_filename(image.filename));
    update_display(image);
  });
}

void close_image(app_state& app) {
  app.images.erase(app.images.begin() + app.selected);
  app.selected = app.images.empty() ? -1 : 0;
}

void draw_glwidgets(const opengl_window& win) {
  static string load_path = "", save_path = "", error_message = "";
  auto&         app = *(app_state*)get_gluser_pointer(win);
  auto image_ok = !app.images.empty() && app.selected >= 0 && app.images[app.selected].load_done;
  if (!begin_glwidgets_window(win, "yimview")) return;
  draw_glmessages(win);
  if (draw_glfiledialog(win, "load image", load_path, false, "./", "",
          "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    load_image_async(app, load_path);
  }
  if (draw_glfiledialog(win, "save image", save_path, true,
          get_dirname(save_path), get_filename(save_path),
          "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    app.images[app.selected].outname = save_path;
    // app.images[app.selected].task_queue.emplace_back(app_task_type::save);
    save_path = "";
    // TODO> implement save
  }
  if (draw_glbutton(win, "load")) {
    open_glmodal(win, "load image");
  }
  continue_glline(win);
  if (draw_glbutton(win, "save",
          app.selected >= 0 && app.images[app.selected].display_done)) {
    save_path = app.images[app.selected].outname;
    open_glmodal(win, "save image");
  }
  continue_glline(win);
  if (draw_glbutton(win, "close", app.selected >= 0)) {
    close_image(app);
  }
  continue_glline(win);
  if (draw_glbutton(win, "quit")) {
    set_glwindow_close(win, true);
  }
  if (app.images.empty()) return;
  draw_glcombobox(
      win, "image", app.selected, (int)app.images.size(),
      [&app](int idx) { return app.images[idx].name.c_str(); }, false);
  if (image_ok && begin_glheader(win, "tonemap")) {
    auto& image = app.images[app.selected];
    auto params = image.tonemap_prms;
    auto edited = 0;
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
      auto wb      = 1 / xyz(image.source_stats.average);
      params.tint = wb / max(wb);
      edited += 1;
    }
    if (edited) {
      image.tonemap_prms = params;
      update_display(image);
      update_texture(image);
    }
    end_glheader(win);
  }
  if (image_ok && begin_glheader(win, "colorgrade")) {
    auto& image = app.images[app.selected];
    auto apply_colorgrade = image.apply_colorgrade;
    auto params = image.colorgrade_prms;
    auto edited = 0;
    edited += draw_glcheckbox(win, "apply colorgrade", apply_colorgrade);
    edited += draw_glslider(win, "contrast", params.contrast, 0, 1);
    edited += draw_glslider(win, "ldr shadows", params.shadows, 0, 1);
    edited += draw_glslider(win, "ldr midtones", params.midtones, 0, 1);
    edited += draw_glslider(win, "highlights", params.highlights, 0, 1);
    edited += draw_glcoloredit(win, "shadows color", params.shadows_color);
    edited += draw_glcoloredit(win, "midtones color", params.midtones_color);
    edited += draw_glcoloredit(win, "highlights color", params.highlights_color);
    if (edited) {
      image.apply_colorgrade = apply_colorgrade;
      image.colorgrade_prms = params;
      update_display(image);
      update_texture(image);
    }
    end_glheader(win);
  }
  if (image_ok && begin_glheader(win, "inspect")) {
    auto& image = app.images[app.selected];
    draw_gllabel(win, "image", get_filename(image.filename));
    draw_gllabel(win, "filename", image.filename);
    draw_gllabel(win, "outname", image.outname);
    draw_gllabel(win, "image", "%d x %d", image.source.size().x, image.source.size().y);
    draw_glslider(win, "zoom", image.image_scale, 0.1, 10);
    draw_glcheckbox(win, "zoom to fit", image.zoom_to_fit);
    auto mouse_pos = get_glmouse_pos(win);
    auto ij        = get_image_coords(
        mouse_pos, image.image_center, image.image_scale, image.source.size());
    draw_gldragger(win, "mouse", ij);
    auto img_pixel = zero4f, display_pixel = zero4f;
    if (ij.x >= 0 && ij.x < image.source.size().x && ij.y >= 0 &&
        ij.y < image.source.size().y) {
      img_pixel     = image.source[{ij.x, ij.y}];
      display_pixel = image.display[{ij.x, ij.y}];
    }
    draw_glcoloredit(win, "image", img_pixel);
    draw_gldragger(win, "display", display_pixel);
    auto img_stats = (image.load_done) ? image.source_stats : image_stats{};
    draw_gldragger(win, "image min", img_stats.min);
    draw_gldragger(win, "image max", img_stats.max);
    draw_gldragger(win, "image avg", img_stats.average);
    draw_glhistogram(win, "image histo", img_stats.histogram);
    auto display_stats = (image.load_done) ? image.display_stats : image_stats{};
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
  auto image_ok = !app.images.empty() && app.selected >= 0 && app.images[app.selected].load_done;
  if(image_ok) {
    auto& image = app.images.at(app.selected);
    if(!image.gl_txt) update_texture(image);
    update_imview(image.image_center, image.image_scale, image.display.size(),
        win_size, image.zoom_to_fit);
    draw_glimage_background(image.gl_txt, win_size.x, win_size.y,
        image.image_center, image.image_scale);
    set_glblending(true);
    draw_glimage(image.gl_txt, win_size.x, win_size.y, image.image_center,
        image.image_scale);
    set_glblending(false);
  }
  begin_glwidgets(win);
  draw_glwidgets(win);
  end_glwidgets(win);
  swap_glbuffers(win);
}

void update(const opengl_window& win, app_state& app) {
  for(auto& image : app.images) {
    if (is_valid(image.load_worker)) {
      if(!is_ready(image.load_worker)) return;
      try {
        image.load_worker.get();
        image.load_done = true;
      } catch(const std::exception& e) {
        push_glmessage(win, "cannot load image " + image.filename);
        log_glinfo(win, "cannot load image " + image.filename);
        log_glinfo(win, e.what());
        break;
      }
    }
  }
}

void drop_callback(const opengl_window& win, const vector<string>& paths) {
  auto& app = *(app_state*)get_gluser_pointer(win);
  for (auto path : paths) load_image_async(app, path);
}

void run_ui(app_state& app) {
  // window
  auto win = opengl_window();
  init_glwindow(win, {1280 + 320, 720}, "yimview", &app, draw);
  set_drop_glcallback(win, drop_callback);

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
      auto& img = app.images.at(app.selected);
      img.image_center += mouse_pos - last_pos;
    }
    if (mouse_right && !widgets_active) {
      auto& img = app.images.at(app.selected);
      img.image_scale *= powf(2, (mouse_pos.x - last_pos.x) * 0.001f);
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
  for (auto filename : filenames) load_image_async(app, filename);
  app.selected = 0;

  // run ui
  run_ui(app);

  // done
  return 0;
}
