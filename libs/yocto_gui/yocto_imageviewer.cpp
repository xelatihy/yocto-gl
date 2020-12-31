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

// grab input
static imageview_image* get_image(imageview_state* state, const string& name);
static imageview_input* get_input(imageview_state* state, const string& name);

// make an image viewer
unique_ptr<imageview_state> make_imageview(const string& title) {
  auto state = make_unique<imageview_state>();
  // state->name = title;
  return state;
}

// Set image
void set_image(imageview_state* state, const string& name,
    const image<vec4f>& img, float exposure, bool filmic) {
  auto lock  = std::lock_guard{state->input_mutex};
  auto input = get_input(state, name);
  if (!input) {
    input =
        state->inputs.emplace_back(std::make_unique<imageview_input>()).get();
  }
  input->name     = name;
  input->hdr      = img;
  input->ldr      = {};
  input->exposure = exposure;
  input->filmic   = filmic;
  input->ichanged = true;
}
void set_image(
    imageview_state* state, const string& name, const image<vec4b>& img) {
  auto lock  = std::lock_guard{state->input_mutex};
  auto input = get_input(state, name);
  if (!input) {
    input =
        state->inputs.emplace_back(std::make_unique<imageview_input>()).get();
  }
  input->name     = name;
  input->hdr      = {};
  input->ldr      = img;
  input->exposure = 0;
  input->filmic   = false;
  input->ichanged = true;
}
// Close image
void close_image(imageview_state* state, const string& name) {
  auto lock  = std::lock_guard{state->input_mutex};
  auto input = get_input(state, name);
  if (!input) return;
  input->close = true;
}

// Set params
void set_param(imageview_state* state, const string& name, const string& pname,
    const json_value& param, const json_value& schema) {
  auto lock  = std::lock_guard{state->input_mutex};
  auto input = get_input(state, name);
  if (!input) return;
  // TODO(fabio): implement this
}
void set_params(imageview_state* state, const string& name,
    const json_value& params, const json_value& schema) {
  auto lock  = std::lock_guard{state->input_mutex};
  auto input = get_input(state, name);
  if (!input) return;
  input->params = params;
  input->schema = schema;
}

// Callback
void set_callback(imageview_state* state, const imageview_callback& callback) {
  state->callback = callback;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE VIEWER INTERNALS
// -----------------------------------------------------------------------------
namespace yocto {

// grab input
static imageview_image* get_image(imageview_state* state, const string& name) {
  for (auto& img : state->images)
    if (img->name == name) return img.get();
  return nullptr;
}
static imageview_input* get_input(imageview_state* state, const string& name) {
  for (auto& img : state->inputs)
    if (img->name == name) return img.get();
  return nullptr;
}

static void update_display(imageview_image* img) {
  if (!img->hdr.empty()) {
    if (img->display.imsize() != img->hdr.imsize())
      img->display.resize(img->hdr.imsize());
    tonemap_image_mt(img->display, img->hdr, img->exposure, img->filmic);
  } else if (!img->ldr.empty()) {
    img->display = img->ldr;
  } else {
    // TODO(fabio): decide about empty images
  }
  if (!is_initialized(img->glimage)) init_image(img->glimage);
  set_image(img->glimage, img->display, false, false);
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
    set_close(win, true);
  }
  draw_combobox(win, "image", state->selected, state->images, false);
  if (!state->selected) return;
  if (begin_header(win, "inspect")) {
    auto img = state->selected;
    draw_label(win, "name", img->name);
    auto size = img->display.imsize();
    draw_dragger(win, "size", size);
    draw_slider(win, "zoom", img->glparams.scale, 0.1, 10);
    draw_checkbox(win, "fit", img->glparams.fit);
    auto ij = image_coords(input.mouse_pos, img->glparams.center,
        img->glparams.scale, img->display.imsize());
    draw_dragger(win, "mouse", ij);
    auto hdr_pixel     = zero4f;
    auto ldr_pixel     = zero4b;
    auto display_pixel = zero4b;
    if (ij.x >= 0 && ij.x < img->display.width() && ij.y >= 0 &&
        ij.y < img->display.height()) {
      hdr_pixel     = !img->hdr.empty() ? img->hdr[{ij.x, ij.y}] : zero4f;
      ldr_pixel     = !img->ldr.empty() ? img->ldr[{ij.x, ij.y}] : zero4b;
      display_pixel = img->display[{ij.x, ij.y}];
    }
    if (!img->hdr.empty()) {
      draw_coloredit(win, "source", hdr_pixel);
    } else {
      draw_coloredit(win, "source", ldr_pixel);
    }
    draw_coloredit(win, "display", display_pixel);
    end_header(win);
  }
  if (!state->selected->hdr.empty()) {
    if (begin_header(win, "tonemap")) {
      auto img    = state->selected;
      auto edited = 0;
      edited += draw_slider(win, "exposure", img->exposure, -5, 5);
      edited += draw_checkbox(win, "filmic", img->filmic);
      if (edited) update_display(img);
      end_header(win);
    }
  }
  if (!state->selected->params.empty()) {
    if (draw_params(win, "params", state->selected->params,
            state->selected->schema, false)) {
      if (state->callback)
        state->callback(state->selected->name, state->selected->params, {});
    }
  }
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
  // process inputs
  auto lock = std::lock_guard{state->input_mutex};

  // close images
  for (auto idx = (size_t)0; idx < state->inputs.size(); idx++) {
    if (!state->inputs[idx]->close) continue;
    if (state->selected == state->images[idx].get()) state->selected = nullptr;
    state->inputs.erase(state->inputs.begin() + idx);
    state->images.erase(state->images.begin() + idx);
    idx--;
  }

  // add images
  for (auto idx = (size_t)0; idx < state->inputs.size(); idx++) {
    if (idx >= state->images.size()) {
      state->images.emplace_back(std::make_unique<imageview_image>());
      state->images[idx]->name = state->inputs[idx]->name;
    }
  }

  // update images
  for (auto idx = (size_t)0; idx < state->inputs.size(); idx++) {
    if (state->inputs[idx]->ichanged) {
      state->images[idx]->hdr      = state->inputs[idx]->hdr;
      state->images[idx]->ldr      = state->inputs[idx]->ldr;
      state->inputs[idx]->ichanged = false;
      update_display(state->images[idx].get());
    }
    if (state->inputs[idx]->pchanged) {
      state->images[idx]->params   = state->inputs[idx]->params;
      state->images[idx]->schema   = state->inputs[idx]->schema;
      state->inputs[idx]->pchanged = false;
    }
  }

  // selected
  if (state->selected == nullptr && !state->images.empty())
    state->selected = state->images[0].get();
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
    if (input.widgets_active) return;
    auto img = state->selected;
    // handle mouse
    if (input.modifier_alt) {
      if (input.mouse_left) {
        img->glparams.center += input.mouse_pos - input.mouse_last;
      }
      if (input.mouse_right) {
        img->glparams.scale *= powf(
            2, (input.mouse_pos.x - input.mouse_last.x) * 0.001f);
      }
    } else {
      if (state->callback) state->callback(state->selected->name, {}, input);
    }
  };

  // run ui
  run_ui({1280 + 320, 720}, "yimview", callbacks);
}

}  // namespace yocto
