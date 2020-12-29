//
// Simpler image viewer.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//

#include "yocto_imageviewer.h"

#include <yocto/yocto_commonio.h>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::make_unique;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE VIEWER API
// -----------------------------------------------------------------------------
namespace yocto {

// make an image viewer
unique_ptr<imageview_state> make_imageview(const string& title) {
  auto state = make_unique<imageview_state>();
  // state->name = title;
  return state;
}

// Open and image viewer
unique_ptr<imageview_state> open_imageview(const string& title) {
  auto state = make_unique<imageview_state>();
  // state->name   = title;
  state->runner = std::async(std::launch::async, run_view, state.get());
  return state;
}

// Wait for the viewer to close
void wait_view(imageview_state* state) { state->runner.wait(); }

// Close viewer
void close_view(imageview_state* state) {
  state->queue.push(imageview_command{imageview_command_type::quit});
}

// Set image
void set_image(
    imageview_state* state, const string& name, const image<vec4f>& img) {
  state->queue.push(
      imageview_command{imageview_command_type::set, name, img, {}});
}
void set_image(
    imageview_state* state, const string& name, const image<vec4b>& img) {
  state->queue.push(
      imageview_command{imageview_command_type::set, name, {}, img});
}
void close_image(imageview_state* state, const string& name) {
  state->queue.push(imageview_command{imageview_command_type::close, name});
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE VIEWER INTERNALS
// -----------------------------------------------------------------------------
namespace yocto {

static void update_display(imageview_image* img) {
  if (img->display.imsize() != img->source.imsize()) img->display = img->source;
  tonemap_image_mt(img->display, img->source, img->exposure, img->filmic);
  // if (app->colorgrade) {
  //   colorgrade_image_mt(app->display, app->source, true, app->params);
  // } else {
  //   tonemap_image_mt(app->display, app->source, app->exposure, app->filmic);
  // }
}

void draw_widgets(
    gui_window* win, imageview_state* state, const gui_input& input) {
  static string load_path = "", save_path = "", error_message = "";
  if (draw_filedialog_button(win, "load", true, "load image", load_path, false,
          "./", "", "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    // load_image_async(state, load_path);
    load_path = "";
  }
  continue_line(win);
  if (draw_filedialog_button(win, "save", state->selected, "save image",
          save_path, true, path_dirname(save_path), path_filename(save_path),
          "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    // state->selected->outname = save_path;
    // save_image(img->outname, img->display, img->error);
    save_path = "";
  }
  continue_line(win);
  if (draw_button(win, "close", (bool)state->selected)) {
    close_image(state, state->selected->name);
  }
  continue_line(win);
  if (draw_button(win, "quit")) {
    close_view(state);
  }
  // draw_combobox(win, "image", state->selected, state->images, false);
  if (!state->selected) return;
  // if (begin_header(win, "tonemap")) {
  //   auto edited = 0;
  //   edited += draw_slider(win, "exposure", app->exposure, -5, 5);
  //   edited += draw_checkbox(win, "filmic", app->filmic);
  //   if (edited) update_display(app);
  //   end_header(win);
  // }
  // if (begin_header(win, "colorgrade")) {
  //   auto& params = app->params;
  //   auto  edited = 0;
  //   edited += draw_checkbox(win, "apply colorgrade", app->colorgrade);
  //   edited += draw_slider(win, "exposure", params.exposure, -5, 5);
  //   edited += draw_coloredit(win, "tint", params.tint);
  //   edited += draw_slider(win, "lincontrast", params.lincontrast, 0, 1);
  //   edited += draw_slider(win, "logcontrast", params.logcontrast, 0, 1);
  //   edited += draw_slider(win, "linsaturation", params.linsaturation, 0, 1);
  //   edited += draw_checkbox(win, "filmic", params.filmic);
  //   continue_line(win);
  //   edited += draw_checkbox(win, "srgb", params.srgb);
  //   continue_line(win);
  //   if (draw_button(win, "auto wb")) {
  //     auto wb     = 1 / xyz(app->source_stats.average);
  //     params.tint = wb / max(wb);
  //     edited += 1;
  //   }
  //   edited += draw_slider(win, "contrast", params.contrast, 0, 1);
  //   edited += draw_slider(win, "saturation", params.saturation, 0, 1);
  //   edited += draw_slider(win, "shadows", params.shadows, 0, 1);
  //   edited += draw_slider(win, "midtones", params.midtones, 0, 1);
  //   edited += draw_slider(win, "highlights", params.highlights, 0, 1);
  //   edited += draw_coloredit(win, "shadows color", params.shadows_color);
  //   edited += draw_coloredit(win, "midtones color", params.midtones_color);
  //   edited += draw_coloredit(win, "highlights color",
  //   params.highlights_color); if (edited) update_display(app);
  //   end_header(win);
  // }
  // if (begin_header(win, "inspect")) {
  //   draw_label(win, "image", app->name);
  //   draw_label(win, "filename", app->filename);
  //   draw_label(win, "outname", app->outname);
  //   draw_label(win, "image",
  //       std::to_string(app->source.width()) + " x " +
  //           std::to_string(app->source.height()));
  //   draw_slider(win, "zoom", app->glparams.scale, 0.1, 10);
  //   draw_checkbox(win, "fit", app->glparams.fit);
  //   auto ij = image_coords(input.mouse_pos, app->glparams.center,
  //       app->glparams.scale, app->source.imsize());
  //   draw_dragger(win, "mouse", ij);
  //   auto img_pixel = zero4f, display_pixel = zero4f;
  //   if (ij.x >= 0 && ij.x < app->source.width() && ij.y >= 0 &&
  //       ij.y < app->source.height()) {
  //     img_pixel     = app->source[{ij.x, ij.y}];
  //     display_pixel = app->display[{ij.x, ij.y}];
  //   }
  //   draw_coloredit(win, "image", img_pixel);
  //   draw_dragger(win, "display", display_pixel);
  //   draw_dragger(win, "image min", app->source_stats.min);
  //   draw_dragger(win, "image max", app->source_stats.max);
  //   draw_dragger(win, "image avg", app->source_stats.average);
  //   draw_histogram(win, "image histo", app->source_stats.histogram);
  //   draw_dragger(win, "display min", app->display_stats.min);
  //   draw_dragger(win, "display max", app->display_stats.max);
  //   draw_dragger(win, "display avg", app->display_stats.average);
  //   draw_histogram(win, "display histo", app->display_stats.histogram);
  //   end_header(win);
  // }
}

void draw(gui_window* win, imageview_state* state, const gui_input& input) {
  if (!state->selected) {
    clear_ogl_framebuffer(ogl_image_params{}.background);
    return;
  }
  auto img                  = state->selected;
  img->glparams.window      = input.window_size;
  img->glparams.framebuffer = input.framebuffer_viewport;
  if (!is_initialized(img->glimage)) init_image(img->glimage);
  std::tie(img->glparams.center, img->glparams.scale) = camera_imview(
      img->glparams.center, img->glparams.scale, img->display.imsize(),
      img->glparams.window, img->glparams.fit);
  draw_image(img->glimage, img->glparams);
}

void update(gui_window* win, imageview_state* state, const gui_input& input) {
  auto command     = imageview_command{};
  auto has_command = state->queue.try_pop(command);
  if (has_command) {
    switch (command.type) {
      case imageview_command_type::set: {
        auto img = (imageview_image*)nullptr;
        for (auto& image : state->images) {
          if (image->name == command.name) {
            img = image.get();
            break;
          }
        }
        if (img == nullptr) {
          state->images.emplace_back(make_unique<imageview_image>());
          img  = state->images.back().get();
          img->name = command.name;
          if (state->selected == nullptr) state->selected = img;
        }
        img->source = command.hdr;
        update_display(img);
        if (!is_initialized(img->glimage)) init_image(img->glimage);
        set_image(img->glimage, img->display, false, false);
      } break;
      case imageview_command_type::close: {
        auto position = -1;
        for (auto pos = (size_t)0; pos < state->images.size(); pos++) {
          if (state->images[pos]->name == command.name) position = (int)pos;
        }
        if (position >= 0) {
          auto fix_selected = state->selected == state->images[position].get();
          state->images.erase(state->images.begin() + position);
          if (fix_selected)
            state->selected =
                state->images.empty() ? nullptr : state->images.front().get();
        }
      } break;
      case imageview_command_type::quit: {
        set_close(win, true);
      } break;
    }
  }
}

// Run application
void run_view(imageview_state* state) {
  // callbacks
  auto callbacks     = gui_callbacks{};
  callbacks.clear_cb = [state](gui_window* win, const gui_input& input) {
    for (auto& image : state->images) clear_image(image->glimage);
  };
  callbacks.update_cb = [state](gui_window* win, const gui_input& input) {
    update(win, state, input);
  };
  callbacks.draw_cb = [state](gui_window* win, const gui_input& input) {
    draw(win, state, input);
  };
  callbacks.widgets_cb = [state](gui_window* win, const gui_input& input) {
    draw_widgets(win, state, input);
  };
  callbacks.uiupdate_cb = [state](gui_window* win, const gui_input& input) {
    if (!state->selected) return;
    // auto app = state->selected;
    // // handle mouse
    // if (input.mouse_left && !input.widgets_active) {
    //   app->glparams.center += input.mouse_pos - input.mouse_last;
    // }
    // if (input.mouse_right && !input.widgets_active) {
    //   app->glparams.scale *= powf(
    //       2, (input.mouse_pos.x - input.mouse_last.x) * 0.001f);
    // }
  };

  // run ui
  run_ui({1280 + 320, 720}, "yimview", callbacks);
}

}  // namespace yocto
