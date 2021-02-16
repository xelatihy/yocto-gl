//
// Simpler image viewer.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
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

#include "yocto_glview.h"

#include <yocto/yocto_cli.h>
#include <yocto/yocto_geometry.h>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::make_unique;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE AND TRACE VIEW
// -----------------------------------------------------------------------------
namespace yocto {

// Open a window and show an image
void view_image(const string& title, const string& name,
    const image_data& image, const progress_callback& progress_cb) {
  // display image
  if (progress_cb) progress_cb("tonemap image", 0, 1);
  auto  display  = make_image(image.width, image.height, false, true);
  float exposure = 0;
  bool  filmic   = false;
  tonemap_image_mt(display, image, exposure, filmic);
  if (progress_cb) progress_cb("tonemap image", 1, 1);

  // opengl image
  auto glimage  = ogl_image{};
  auto glparams = ogl_image_params{};

  // top level combo
  auto names    = vector<string>{name};
  auto selected = 0;

  // callbacks
  auto callbacks    = gui_callbacks{};
  callbacks.init_cb = [&](gui_window* win, const gui_input& input) {
    init_image(glimage);
    set_image(glimage, display, false, false);
  };
  callbacks.clear_cb = [&](gui_window* win, const gui_input& input) {
    clear_image(glimage);
  };
  callbacks.draw_cb = [&](gui_window* win, const gui_input& input) {
    glparams.window                           = input.window_size;
    glparams.framebuffer                      = input.framebuffer_viewport;
    std::tie(glparams.center, glparams.scale) = camera_imview(glparams.center,
        glparams.scale, {display.width, display.height}, glparams.window,
        glparams.fit);
    draw_image(glimage, glparams);
  };
  callbacks.widgets_cb = [&](gui_window* win, const gui_input& input) {
    auto edited = 0;
    draw_combobox(win, "name", selected, names);
    if (begin_header(win, "tonemap")) {
      edited += draw_slider(win, "exposure", exposure, -5, 5);
      edited += draw_checkbox(win, "filmic", filmic);
      end_header(win);
      if (edited) {
        if (progress_cb) progress_cb("tonemap image", 0, 1);
        tonemap_image_mt(display, image, exposure, filmic);
        set_image(glimage, display, false, false);
        if (progress_cb) progress_cb("tonemap image", 0, 1);
      }
    }
    if (begin_header(win, "inspect")) {
      draw_slider(win, "zoom", glparams.scale, 0.1, 10);
      draw_checkbox(win, "fit", glparams.fit);
      auto [i, j] = image_coords(input.mouse_pos, glparams.center,
          glparams.scale, {image.width, image.height});
      auto ij     = vec2i{i, j};
      draw_dragger(win, "mouse", ij);
      auto hdr_pixel     = zero4f;
      auto ldr_pixel     = zero4b;
      auto display_pixel = zero4b;
      if (i >= 0 && i < image.width && j >= 0 && j < image.height) {
        display_pixel = image.pixelsb[j * image.width + i];
        if (!image.pixelsf.empty())
          hdr_pixel = image.pixelsf[j * image.width + i];
        if (!image.pixelsb.empty())
          ldr_pixel = image.pixelsb[j * image.width + i];
      }
      if (!image.pixelsf.empty()) {
        draw_coloredit(win, "image", hdr_pixel);
      } else {
        draw_coloredit(win, "image", ldr_pixel);
      }
      draw_coloredit(win, "display", display_pixel);
      end_header(win);
    }
  };
  callbacks.uiupdate_cb = [&glparams](gui_window* win, const gui_input& input) {
    // handle mouse
    if (input.mouse_left && !input.widgets_active) {
      glparams.center += input.mouse_pos - input.mouse_last;
    }
    if (input.mouse_right && !input.widgets_active) {
      glparams.scale *= powf(
          2, (input.mouse_pos.x - input.mouse_last.x) * 0.001f);
    }
  };

  // run ui
  run_ui({1280 + 320, 720}, title, callbacks);
}

// Open a window and show an image
void view_images(const string& title, const vector<string>& names,
    const vector<image_data>& images, const progress_callback& progress_cb) {
  // display image
  if (progress_cb) progress_cb("tonemap image", 0, (int)images.size());
  auto displays  = vector<image_data>(images.size());
  auto exposures = vector<float>(images.size(), 0);
  auto filmics   = vector<bool>(images.size(), false);
  for (auto idx = 0; idx < (int)images.size(); idx++) {
    if (progress_cb) progress_cb("tonemap image", idx, (int)images.size());
    displays[idx] = make_image(
        images[idx].width, images[idx].height, false, true);
    tonemap_image_mt(displays[idx], images[idx], exposures[idx], filmics[idx]);
  }
  if (progress_cb)
    progress_cb("tonemap image", (int)images.size(), (int)images.size());

  // opengl image
  auto glimages  = vector<ogl_image>(images.size());
  auto glparamss = vector<ogl_image_params>(images.size());

  // selection
  auto selected = 0;

  // callbacks
  auto callbacks    = gui_callbacks{};
  callbacks.init_cb = [&glimages, &displays, &images](
                          gui_window* win, const gui_input& input) {
    for (auto idx = 0; idx < (int)images.size(); idx++) {
      init_image(glimages[idx]);
      set_image(glimages[idx], displays[idx], false, false);
    }
  };
  callbacks.clear_cb = [&glimages, &images](
                           gui_window* win, const gui_input& input) {
    for (auto idx = 0; idx < (int)images.size(); idx++) {
      clear_image(glimages[idx]);
    }
  };
  callbacks.draw_cb = [&glimages, &glparamss, &displays, &selected](
                          gui_window* win, const gui_input& input) {
    auto& glparams                            = glparamss[selected];
    auto& display                             = displays[selected];
    auto& glimage                             = glimages[selected];
    glparams.window                           = input.window_size;
    glparams.framebuffer                      = input.framebuffer_viewport;
    std::tie(glparams.center, glparams.scale) = camera_imview(glparams.center,
        glparams.scale, {display.width, display.height}, glparams.window,
        glparams.fit);
    draw_image(glimage, glparams);
  };
  callbacks.widgets_cb = [&](gui_window* win, const gui_input& input) {
    auto edited = 0;
    if (draw_combobox(win, "name", selected, names)) {
    }
    auto& display  = displays[selected];
    auto& glimage  = glimages[selected];
    auto& image    = images[selected];
    auto& exposure = exposures[selected];
    auto  filmic   = (bool)filmics[selected];  // vector of bool ...
    if (begin_header(win, "tonemap")) {
      edited += draw_slider(win, "exposure", exposure, -5, 5);
      edited += draw_checkbox(win, "filmic", filmic);
      end_header(win);
      if (edited) {
        if (progress_cb) progress_cb("tonemap image", 0, 1);
        filmics[selected] = filmic;
        tonemap_image_mt(display, image, exposure, filmic);
        set_image(glimage, display, false, false);
        if (progress_cb) progress_cb("tonemap image", 1, 1);
      }
    }
    auto& glparams = glparamss[selected];
    if (begin_header(win, "inspect")) {
      draw_slider(win, "zoom", glparams.scale, 0.1, 10);
      draw_checkbox(win, "fit", glparams.fit);
      auto [i, j] = image_coords(input.mouse_pos, glparams.center,
          glparams.scale, {image.width, image.height});
      auto ij     = vec2i{i, j};
      draw_dragger(win, "mouse", ij);
      auto hdr_pixel     = zero4f;
      auto ldr_pixel     = zero4b;
      auto display_pixel = zero4b;
      if (i >= 0 && i < image.width && j >= 0 && j < image.height) {
        display_pixel = image.pixelsb[j * image.width + i];
        if (!image.pixelsf.empty())
          hdr_pixel = image.pixelsf[j * image.width + i];
        if (!image.pixelsb.empty())
          ldr_pixel = image.pixelsb[j * image.width + i];
      }
      if (!image.pixelsf.empty()) {
        draw_coloredit(win, "image", hdr_pixel);
      } else {
        draw_coloredit(win, "image", ldr_pixel);
      }
      draw_coloredit(win, "display", display_pixel);
      end_header(win);
    }
  };
  callbacks.uiupdate_cb = [&](gui_window* win, const gui_input& input) {
    // handle mouse
    auto& glparams = glparamss[selected];
    if (input.mouse_left && !input.widgets_active) {
      glparams.center += input.mouse_pos - input.mouse_last;
    }
    if (input.mouse_right && !input.widgets_active) {
      glparams.scale *= powf(
          2, (input.mouse_pos.x - input.mouse_last.x) * 0.001f);
    }
  };

  // run ui
  run_ui({1280 + 320, 720}, title, callbacks);
}

// Open a window and show an image
void colorgrade_image(const string& title, const string& name,
    const image_data& image, const progress_callback& progress_cb) {
  // color grading parameters
  auto params = colorgrade_params{};

  // display image
  if (progress_cb) progress_cb("colorgrade image", 0, 1);
  auto display = make_image(image.width, image.height, false, true);
  colorgrade_image_mt(display, image, params);
  if (progress_cb) progress_cb("colorgrade image", 1, 1);

  // opengl image
  auto glimage  = ogl_image{};
  auto glparams = ogl_image_params{};

  // top level combo
  auto names    = vector<string>{name};
  auto selected = 0;

  // callbacks
  auto callbacks    = gui_callbacks{};
  callbacks.init_cb = [&](gui_window* win, const gui_input& input) {
    init_image(glimage);
    set_image(glimage, display, false, false);
  };
  callbacks.clear_cb = [&](gui_window* win, const gui_input& input) {
    clear_image(glimage);
  };
  callbacks.draw_cb = [&](gui_window* win, const gui_input& input) {
    glparams.window                           = input.window_size;
    glparams.framebuffer                      = input.framebuffer_viewport;
    std::tie(glparams.center, glparams.scale) = camera_imview(glparams.center,
        glparams.scale, {display.width, display.height}, glparams.window,
        glparams.fit);
    draw_image(glimage, glparams);
  };
  callbacks.widgets_cb = [&](gui_window* win, const gui_input& input) {
    draw_combobox(win, "name", selected, names);
    if (begin_header(win, "colorgrade")) {
      auto edited = 0;
      edited += draw_slider(win, "exposure", params.exposure, -5, 5);
      edited += draw_coloredit(win, "tint", params.tint);
      edited += draw_slider(win, "lincontrast", params.lincontrast, 0, 1);
      edited += draw_slider(win, "logcontrast", params.logcontrast, 0, 1);
      edited += draw_slider(win, "linsaturation", params.linsaturation, 0, 1);
      edited += draw_checkbox(win, "filmic", params.filmic);
      continue_line(win);
      edited += draw_checkbox(win, "srgb", params.srgb);
      edited += draw_slider(win, "contrast", params.contrast, 0, 1);
      edited += draw_slider(win, "saturation", params.saturation, 0, 1);
      edited += draw_slider(win, "shadows", params.shadows, 0, 1);
      edited += draw_slider(win, "midtones", params.midtones, 0, 1);
      edited += draw_slider(win, "highlights", params.highlights, 0, 1);
      edited += draw_coloredit(win, "shadows color", params.shadows_color);
      edited += draw_coloredit(win, "midtones color", params.midtones_color);
      edited += draw_coloredit(
          win, "highlights color", params.highlights_color);
      end_header(win);
      if (edited) {
        if (progress_cb) progress_cb("colorgrade image", 0, 1);
        colorgrade_image_mt(display, image, params);
        set_image(glimage, display, false, false);
        if (progress_cb) progress_cb("colorgrade image", 0, 1);
      }
    }
    if (begin_header(win, "inspect")) {
      draw_slider(win, "zoom", glparams.scale, 0.1, 10);
      draw_checkbox(win, "fit", glparams.fit);
      auto [i, j] = image_coords(input.mouse_pos, glparams.center,
          glparams.scale, {image.width, image.height});
      auto ij     = vec2i{i, j};
      draw_dragger(win, "mouse", ij);
      auto hdr_pixel     = zero4f;
      auto ldr_pixel     = zero4b;
      auto display_pixel = zero4b;
      if (i >= 0 && i < image.width && j >= 0 && j < image.height) {
        display_pixel = image.pixelsb[j * image.width + i];
        if (!image.pixelsf.empty())
          hdr_pixel = image.pixelsf[j * image.width + i];
        if (!image.pixelsb.empty())
          ldr_pixel = image.pixelsb[j * image.width + i];
      }
      if (!image.pixelsf.empty()) {
        draw_coloredit(win, "image", hdr_pixel);
      } else {
        draw_coloredit(win, "image", ldr_pixel);
      }
      draw_coloredit(win, "display", display_pixel);
      end_header(win);
    }
  };
  callbacks.uiupdate_cb = [&glparams](gui_window* win, const gui_input& input) {
    // handle mouse
    if (input.mouse_left && !input.widgets_active) {
      glparams.center += input.mouse_pos - input.mouse_last;
    }
    if (input.mouse_right && !input.widgets_active) {
      glparams.scale *= powf(
          2, (input.mouse_pos.x - input.mouse_last.x) * 0.001f);
    }
  };

  // run ui
  run_ui({1280 + 320, 720}, title, callbacks);
}

// Open a window and show a shape via path tracing
void view_shape(const string& title, const string& name,
    const shape_data& shape, bool addsky,
    const progress_callback& progress_cb) {
  // initialize path tracer scene
  auto scene = scene_scene{};
  if (progress_cb) progress_cb("create scene", 0, 1);
  scene.shape_names.emplace_back("shape");
  scene.shapes.emplace_back(shape);
  scene.material_names.emplace_back("material");
  auto& material = scene.materials.emplace_back();
  material.color = {0.8, 0.8, 0.8};
  scene.instance_names.emplace_back("instance");
  auto& instance    = scene.instances.emplace_back();
  instance.shape    = 0;
  instance.material = 0;
  add_camera(scene);
  if (addsky) add_sky(scene);
  if (progress_cb) progress_cb("create scene", 0, 1);

  // run view
  view_scene(title, name, scene, progress_cb);
}

// Open a window and show an scene via path tracing
void view_scene(const string& title, const string& name, scene_scene& scene,
    const progress_callback& progress_cb) {
  return view_scene(title, name, scene, "", progress_cb);
}

// Open a window and show an scene via path tracing
void view_scene(const string& title, const string& name, scene_scene& scene,
    const string& camname, const progress_callback& progress_cb) {
  // rendering params
  auto params     = trace_params{};
  auto has_lights = std::any_of(scene.instances.begin(), scene.instances.end(),
                        [&scene](sceneio_instance& instance) {
                          auto& material = scene.materials[instance.material];
                          return material.emission != zero3f;
                        }) ||
                    std::any_of(scene.environments.begin(),
                        scene.environments.end(),
                        [](const sceneio_environment& environment) {
                          return environment.emission != zero3f;
                        });
  if (!has_lights) params.sampler = trace_sampler_type::eyelight;

  // get camera
  params.camera = find_camera(scene, camname);

  // run viewer
  view_scene(title, name, scene, params, progress_cb);
}

#if 1

// Open a window and show an scene via path tracing
void view_scene(const string& title, const string& name, scene_scene& scene,
    const trace_params& params_, const progress_callback& progress_cb) {
  // copy params and camera
  auto params = params_;

  // build bvh
  auto bvh = make_bvh(scene, params, progress_cb);

  // init renderer
  auto lights = make_lights(scene, params, progress_cb);

  // fix renderer type if no lights
  if (lights.lights.empty() && is_sampler_lit(params)) {
    print_info("no lights presents --- switching to eyelight");
    params.sampler = trace_sampler_type::eyelight;
  }

  // init images
  auto trace_size = [](const scene_scene&   scene,
                        const trace_params& params) -> vec2i {
    auto& camera = scene.cameras[params.camera];
    if (camera.aspect >= 1) {
      return {params.resolution, (int)round(params.resolution / camera.aspect)};
    } else {
      return {(int)round(params.resolution * camera.aspect), params.resolution};
    }
  };
  auto [width, height] = trace_size(scene, params);
  auto image           = make_image(width, height, true, false);
  auto display         = make_image(width, height, false, true);
  tonemap_image_mt(display, image, params.exposure);

  // opengl image
  auto glimage  = ogl_image{};
  auto glparams = ogl_image_params{};

  // top level combo
  auto names    = vector<string>{name};
  auto selected = 0;

  // camera names
  auto camera_names = scene.camera_names;
  if (camera_names.empty()) {
    for (auto idx = 0; idx < (int)scene.cameras.size(); idx++) {
      camera_names.push_back("camera" + std::to_string(idx + 1));
    }
  }

  // init state
  auto worker = trace_worker{};
  auto state  = trace_state{};

  // renderer update
  auto render_update  = std::atomic<bool>{};
  auto render_current = std::atomic<int>{};
  auto render_mutex   = std::mutex{};
  auto reset_display  = [&]() {
    // stop render
    trace_stop(worker);

    // start render
    trace_start(
        worker, state, scene, bvh, lights, params,
        [&](const string& message, int sample, int samples) {
          if (progress_cb) progress_cb("render sample", sample, params.samples);
          render_current = sample;
        },
        [&](const image_data& render, int sample, int samples) {
          // if (current > 0) return;
          auto lock = std::lock_guard{render_mutex};
          image     = render;
          tonemap_image_mt(display, image, params.exposure);
          render_update = true;
        });
  };

  // start rendeting
  reset_display();

  // callbacks
  auto callbacks    = gui_callbacks{};
  callbacks.init_cb = [&](gui_window* win, const gui_input& input) {
    auto lock = std::lock_guard{render_mutex};
    init_image(glimage);
    set_image(glimage, display, false, false);
  };
  callbacks.clear_cb = [&](gui_window* win, const gui_input& input) {
    clear_image(glimage);
  };
  callbacks.draw_cb = [&](gui_window* win, const gui_input& input) {
    // update image
    if (render_update) {
      auto lock = std::lock_guard{render_mutex};
      set_image(glimage, display, false, false);
      render_update = false;
    }
    // draw image
    glparams.window                           = input.window_size;
    glparams.framebuffer                      = input.framebuffer_viewport;
    std::tie(glparams.center, glparams.scale) = camera_imview(glparams.center,
        glparams.scale, {display.width, display.height}, glparams.window,
        glparams.fit);
    draw_image(glimage, glparams);
  };
  callbacks.widgets_cb = [&](gui_window* win, const gui_input& input) {
    auto edited = 0;
    draw_combobox(win, "name", selected, names);
    // draw_progressbar(win, "render", app->current, app->total);
    if (begin_header(win, "render")) {
      auto edited  = 0;
      auto tparams = params;
      edited += draw_combobox(win, "camera", tparams.camera, camera_names);
      edited += draw_slider(win, "resolution", tparams.resolution, 180, 4096);
      edited += draw_slider(win, "samples", tparams.samples, 16, 4096);
      edited += draw_combobox(
          win, "tracer", (int&)tparams.sampler, trace_sampler_names);
      edited += draw_combobox(
          win, "false color", (int&)tparams.falsecolor, trace_falsecolor_names);
      edited += draw_slider(win, "bounces", tparams.bounces, 1, 128);
      edited += draw_checkbox(win, "envhidden", tparams.envhidden);
      continue_line(win);
      edited += draw_checkbox(win, "filter", tparams.tentfilter);
      edited += draw_slider(win, "pratio", tparams.pratio, 1, 64);
      // edited += draw_slider(win, "exposure", tparams.exposure, -5, 5);
      end_header(win);
      if (edited) {
        trace_stop(worker);
        // check image changes
        auto [width, height]   = trace_size(scene, params);
        auto [twidth, theight] = trace_size(scene, tparams);
        if (width != twidth || height != theight) {
          image   = make_image(twidth, theight, true, false);
          display = make_image(twidth, theight, false, true);
          set_image(glimage, display, false, false);
        }
        params = tparams;
        reset_display();
      }
    }
    if (begin_header(win, "tonemap")) {
      edited += draw_slider(win, "exposure", params.exposure, -5, 5);
      end_header(win);
      if (edited) {
        if (progress_cb) progress_cb("tonemap image", 0, 1);
        tonemap_image_mt(display, image, params.exposure);
        set_image(glimage, display, false, false);
        if (progress_cb) progress_cb("tonemap image", 0, 1);
      }
    }
    if (begin_header(win, "inspect")) {
      draw_slider(win, "zoom", glparams.scale, 0.1, 10);
      draw_checkbox(win, "fit", glparams.fit);
      auto [i, j] = image_coords(input.mouse_pos, glparams.center,
          glparams.scale, {image.width, image.height});
      auto ij     = vec2i{i, j};
      draw_dragger(win, "mouse", ij);
      auto hdr_pixel     = zero4f;
      auto ldr_pixel     = zero4b;
      auto display_pixel = zero4b;
      if (i >= 0 && i < image.width && j >= 0 && j < image.height) {
        display_pixel = image.pixelsb[j * image.width + i];
        if (!image.pixelsf.empty())
          hdr_pixel = image.pixelsf[j * image.width + i];
        if (!image.pixelsb.empty())
          ldr_pixel = image.pixelsb[j * image.width + i];
      }
      if (!image.pixelsf.empty()) {
        draw_coloredit(win, "image", hdr_pixel);
      } else {
        draw_coloredit(win, "image", ldr_pixel);
      }
      draw_coloredit(win, "display", display_pixel);
      end_header(win);
    }
  };
  callbacks.uiupdate_cb = [&](gui_window* win, const gui_input& input) {
    if ((input.mouse_left || input.mouse_right) && !input.modifier_alt &&
        !input.widgets_active) {
      auto dolly  = 0.0f;
      auto pan    = zero2f;
      auto rotate = zero2f;
      if (input.mouse_left && !input.modifier_shift)
        rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
      if (input.mouse_right)
        dolly = (input.mouse_pos.x - input.mouse_last.x) / 100.0f;
      auto camera = scene.cameras[params.camera];
      if (input.mouse_left && input.modifier_shift)
        pan = (input.mouse_pos - input.mouse_last) * camera.focus / 200.0f;
      pan.x                                = -pan.x;
      std::tie(camera.frame, camera.focus) = camera_turntable(
          camera.frame, camera.focus, rotate, dolly, pan);
      if (camera.frame != scene.cameras[params.camera].frame ||
          camera.focus != scene.cameras[params.camera].focus) {
        trace_stop(worker);
        scene.cameras[params.camera] = camera;
        reset_display();
      }
    }
  };

  // run ui
  run_ui({1280 + 320, 720}, title, callbacks);

  // done
  trace_stop(worker);
}

#else

// Parameter conversions
void from_params(const gui_params& uiparams, trace_params& params) {
  params.camera     = uiparams.at("camera");
  params.resolution = uiparams.at("resolution");
  params.sampler    = uiparams.at("sampler");
  params.falsecolor = uiparams.at("falsecolor");
  params.samples    = uiparams.at("samples");
  params.bounces    = uiparams.at("bounces");
  params.clamp      = uiparams.at("clamp");
  params.nocaustics = uiparams.at("nocaustics");
  params.envhidden  = uiparams.at("envhidden");
  params.tentfilter = uiparams.at("tentfilter");
}
void to_params(gui_params& uiparams, const trace_params& params,
    const vector<string>& camera_names) {
  uiparams["camera"]     = {params.camera, camera_names};
  uiparams["resolution"] = {params.resolution, {128, 4096}};
  uiparams["sampler"]    = {params.sampler, trace_sampler_names};
  uiparams["falsecolor"] = {params.falsecolor, trace_falsecolor_names};
  uiparams["samples"]    = {params.samples, {1, 4096}};
  uiparams["bounces"]    = {params.bounces, {1, 64}};
  uiparams["clamp"]      = {params.clamp, {10, 1000}};
  uiparams["nocaustics"] = params.nocaustics;
  uiparams["envhidden"]  = params.envhidden;
  uiparams["tentfilter"] = params.tentfilter;
}

// Open a window and show an scene via path tracing
void view_scene(const string& title, const string& name, scene_scene& scene,
    const trace_params& params_, const progress_callback& progress_cb) {
  // open viewer
  auto viewer = make_imageviewer(title);

  // copy params and camera
  auto params = params_;

  // build bvh
  auto bvh = make_bvh(scene, params, progress_cb);

  // init renderer
  auto lights = make_lights(scene, params, progress_cb);

  // fix renderer type if no lights
  if (lights.lights.empty() && is_sampler_lit(params)) {
    // TODO(fabio): fix this
  }

  // init state
  auto worker = trace_worker{};
  auto state  = trace_state{};

  // render start
  trace_start(
      worker, state, scene, bvh, lights, params,
      [&viewer, &progress_cb, name](
          const string& message, int sample, int nsamples) {
        set_param(viewer, name, "sample", {sample, {1, 4096}, true});
        if (progress_cb) progress_cb(message, sample, nsamples);
      },
      [&viewer, name](const image_data& render, int current, int total) {
        set_image(viewer, name, render);
      });

  // create camera names if missing
  if (scene.camera_names.empty()) {
    for (auto& camera : scene.cameras) {
      scene.camera_names.push_back(
          "camera" + std::to_string(&camera - scene.cameras.data()));
    }
  }

  // show rendering params
  auto uiparams = gui_params();
  to_params(uiparams, params, scene.camera_names);
  set_params(viewer, name, "render", uiparams);

  // set callback
  set_params_callback(
      viewer, [&](const string& name, const gui_params& uiparams) {
        trace_stop(worker);
        from_params(uiparams, params);
        trace_start(
            worker, state, scene, bvh, lights, params,
            [&viewer, &progress_cb, name](
                const string& message, int sample, int nsamples) {
              set_param(viewer, name, "sample", {sample, {1, 4096}, true});
              if (progress_cb) progress_cb(message, sample, nsamples);
            },
            [&viewer, name](const image_data& render, int current, int total) {
              set_image(viewer, name, render);
            });
      });
  set_input_callback(viewer, [&](const string& name, const gui_input& input) {
    if ((input.mouse_left || input.mouse_right) &&
        input.mouse_pos != input.mouse_last) {
      trace_stop(worker);
      auto dolly  = 0.0f;
      auto pan    = zero2f;
      auto rotate = zero2f;
      if (input.mouse_left && !input.modifier_shift)
        rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
      if (input.mouse_right)
        dolly = (input.mouse_pos.x - input.mouse_last.x) / 100.0f;
      auto& camera = scene.cameras[params.camera];
      if (input.mouse_left && input.modifier_shift)
        pan = (input.mouse_pos - input.mouse_last) * camera.focus / 200.0f;
      pan.x                                = -pan.x;
      std::tie(camera.frame, camera.focus) = camera_turntable(
          camera.frame, camera.focus, rotate, dolly, pan);
      trace_start(
          worker, state, scene, bvh, lights, params,
          [&viewer, &progress_cb, name](
              const string& message, int sample, int nsamples) {
            set_param(viewer, name, "sample", {sample, {1, 4096}, true});
            if (progress_cb) progress_cb(message, sample, nsamples);
          },
          [&viewer, name](const image_data& render, int current, int total) {
            set_image(viewer, name, render);
          });
    }
  });

  // run view
  run_viewer(viewer);

  // stop
  trace_stop(worker);
}

#endif

static void init_glscene(shade_scene& glscene, const sceneio_scene& ioscene,
    progress_callback progress_cb) {
  // handle progress
  auto progress = vec2i{
      0, (int)ioscene.cameras.size() + (int)ioscene.materials.size() +
             (int)ioscene.textures.size() + (int)ioscene.shapes.size() +
             (int)ioscene.instances.size()};

  // init scene
  init_scene(glscene);

  // camera
  for (auto& iocamera : ioscene.cameras) {
    if (progress_cb) progress_cb("convert camera", progress.x++, progress.y);
    auto& camera = glscene.cameras.at(add_camera(glscene));
    set_frame(camera, iocamera.frame);
    set_lens(camera, iocamera.lens, iocamera.aspect, iocamera.film);
    set_nearfar(camera, 0.001, 10000);
  }

  // textures
  for (auto& iotexture : ioscene.textures) {
    if (progress_cb) progress_cb("convert texture", progress.x++, progress.y);
    auto  handle    = add_texture(glscene);
    auto& gltexture = glscene.textures[handle];
    if (!iotexture.pixelsf.empty()) {
      set_texture(
          gltexture, iotexture.width, iotexture.height, iotexture.pixelsf);
    } else if (!iotexture.pixelsb.empty()) {
      set_texture(
          gltexture, iotexture.width, iotexture.height, iotexture.pixelsb);
    }
  }

  // material
  for (auto& iomaterial : ioscene.materials) {
    if (progress_cb) progress_cb("convert material", progress.x++, progress.y);
    auto  handle     = add_material(glscene);
    auto& glmaterial = glscene.materials[handle];
    set_emission(glmaterial, iomaterial.emission, iomaterial.emission_tex);
    set_opacity(glmaterial, iomaterial.opacity, invalid_handle);
    set_normalmap(glmaterial, iomaterial.normal_tex);
    switch (iomaterial.type) {
      case material_type::matte: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 0, invalid_handle);
        set_metallic(glmaterial, 0, invalid_handle);
        set_roughness(glmaterial, 0, invalid_handle);
      } break;
      case material_type::plastic: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 1, invalid_handle);
        set_metallic(glmaterial, 0, invalid_handle);
        set_roughness(glmaterial, iomaterial.roughness, invalid_handle);
      } break;
      case material_type::metal: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 0, invalid_handle);
        set_metallic(glmaterial, 1, invalid_handle);
        set_roughness(
            glmaterial, iomaterial.roughness, iomaterial.roughness_tex);
      } break;
      case material_type::metallic: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 1, invalid_handle);
        set_metallic(glmaterial, iomaterial.metallic, invalid_handle);
        set_roughness(glmaterial, iomaterial.roughness, invalid_handle);
      } break;
      default: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 0, invalid_handle);
        set_metallic(glmaterial, 0, invalid_handle);
        set_roughness(
            glmaterial, iomaterial.roughness, iomaterial.roughness_tex);
      } break;
    }
  }

  // shapes
  for (auto& ioshape : ioscene.shapes) {
    if (progress_cb) progress_cb("convert shape", progress.x++, progress.y);
    add_shape(glscene, ioshape.points, ioshape.lines, ioshape.triangles,
        ioshape.quads, ioshape.positions, ioshape.normals, ioshape.texcoords,
        ioshape.colors);
  }

  // shapes
  for (auto& ioinstance : ioscene.instances) {
    if (progress_cb) progress_cb("convert instance", progress.x++, progress.y);
    auto  handle     = add_instance(glscene);
    auto& glinstance = glscene.instances[handle];
    set_frame(glinstance, ioinstance.frame);
    set_shape(glinstance, ioinstance.shape);
    set_material(glinstance, ioinstance.material);
  }

  // environments
  for (auto& ioenvironment : ioscene.environments) {
    auto  handle        = add_environment(glscene);
    auto& glenvironment = glscene.environments[handle];
    set_frame(glenvironment, ioenvironment.frame);
    set_emission(
        glenvironment, ioenvironment.emission, ioenvironment.emission_tex);
  }

  // init environments
  init_environments(glscene);

  // done
  if (progress_cb) progress_cb("convert done", progress.x++, progress.y);
}

using glview_scene_callback = std::function<void(gui_window* win,
    const gui_input& input, scene_scene& scene, shade_scene& glscene)>;

void glview_scene(scene_scene& scene, const string& name, const string& camname,
    const progress_callback&     progress_cb,
    const glview_scene_callback& widgets_callback,
    const glview_scene_callback& uiupdate_callback) {
  // glscene
  auto glscene = shade_scene{};

  // draw params
  auto params = shade_params{};

  // callbacks
  auto callbacks    = gui_callbacks{};
  callbacks.init_cb = [&glscene, &scene, &progress_cb](
                          gui_window* win, const gui_input& input) {
    init_glscene(glscene, scene, progress_cb);
  };
  callbacks.clear_cb = [&glscene](gui_window* win, const gui_input& input) {
    clear_scene(glscene);
  };
  callbacks.draw_cb = [&glscene, &params](
                          gui_window* win, const gui_input& input) {
    draw_scene(
        glscene, glscene.cameras.at(0), input.framebuffer_viewport, params);
  };
  callbacks.widgets_cb = [&glscene, &scene, &params, &widgets_callback](
                             gui_window* win, const gui_input& input) {
    draw_checkbox(win, "wireframe", params.wireframe);
    continue_line(win);
    draw_checkbox(win, "faceted", params.faceted);
    continue_line(win);
    draw_checkbox(win, "double sided", params.double_sided);
    draw_combobox(win, "lighting", (int&)params.lighting, shade_lighting_names);
    draw_slider(win, "exposure", params.exposure, -10, 10);
    draw_slider(win, "gamma", params.gamma, 0.1f, 4);
    draw_slider(win, "near", params.near, 0.01f, 1.0f);
    draw_slider(win, "far", params.far, 1000.0f, 10000.0f);
    if (widgets_callback) widgets_callback(win, input, scene, glscene);
  };
  callbacks.update_cb = [](gui_window* win, const gui_input& input) {
    // update(win, apps);
  };
  callbacks.uiupdate_cb = [&glscene, &scene, &uiupdate_callback](
                              gui_window* win, const gui_input& input) {
    // handle mouse and keyboard for navigation
    if (uiupdate_callback) uiupdate_callback(win, input, scene, glscene);
    if ((input.mouse_left || input.mouse_right) && !input.modifier_alt &&
        !input.modifier_ctrl && !input.widgets_active) {
      auto dolly  = 0.0f;
      auto pan    = zero2f;
      auto rotate = zero2f;
      if (input.mouse_left && !input.modifier_shift)
        rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
      if (input.mouse_right)
        dolly = (input.mouse_pos.x - input.mouse_last.x) / 100.0f;
      if (input.mouse_left && input.modifier_shift)
        pan = (input.mouse_pos - input.mouse_last) / 100.0f;
      auto& camera                         = scene.cameras.at(0);
      std::tie(camera.frame, camera.focus) = camera_turntable(
          camera.frame, camera.focus, rotate, dolly, pan);
      set_frame(glscene.cameras.at(0), camera.frame);
    }
  };

  // run ui
  run_ui({1280 + 320, 720}, "yshade", callbacks);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE VIEWER API
// -----------------------------------------------------------------------------
namespace yocto {

// grab input
// static imageview_image* get_image(imageview_state* viewer, const string&
// name);
static ogl_imageinput* get_input(ogl_imageviewer& viewer, const string& name);

// make an image viewer
ogl_imageviewer make_imageviewer(const string& title) {
  return ogl_imageviewer{};
}

// Set image
void set_image(
    ogl_imageviewer& viewer, const string& name, const image_data& image) {
  auto lock  = std::lock_guard{viewer.input_mutex};
  auto input = get_input(viewer, name);
  if (!input) {
    input =
        viewer.inputs.emplace_back(std::make_unique<ogl_imageinput>()).get();
  }
  input->name     = name;
  input->image    = image;
  input->exposure = 0;
  input->filmic   = false;
  input->ichanged = true;
}
// Close image
void close_image(ogl_imageviewer& viewer, const string& name) {
  auto lock  = std::lock_guard{viewer.input_mutex};
  auto input = get_input(viewer, name);
  if (!input) return;
  input->close = true;
}

// Set params
void set_param(ogl_imageviewer& viewer, const string& name, const string& pname,
    const gui_param& param) {
  auto lock  = std::lock_guard{viewer.input_mutex};
  auto input = get_input(viewer, name);
  if (!input) return;
  input->params[pname] = param;
  input->pchanged      = true;
}
void set_params(ogl_imageviewer& viewer, const string& name,
    const string& pname, const gui_params& params) {
  auto lock  = std::lock_guard{viewer.input_mutex};
  auto input = get_input(viewer, name);
  if (!input) return;
  input->pname    = pname;
  input->params   = params;
  input->pchanged = true;
}

// Callback
void set_params_callback(
    ogl_imageviewer& viewer, const ogl_imageviewer_pcallback& callback) {
  viewer.pcallback = callback;
}
void set_input_callback(
    ogl_imageviewer& viewer, const ogl_imageviewer_icallback& callback) {
  viewer.icallback = callback;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE VIEWER INTERNALS
// -----------------------------------------------------------------------------
namespace yocto {

// grab input
// static imageview_image* get_image(imageview_state* viewer, const string&
// name)
// {
//   for (auto& img : viewer.images)
//     if (view->name == name) return img.get();
//   return nullptr;
// }
static ogl_imageinput* get_input(ogl_imageviewer& viewer, const string& name) {
  for (auto& input : viewer.inputs)
    if (input->name == name) return input.get();
  return nullptr;
}

static void update_display(ogl_imageview* view) {
  if (view->display.width != view->image.width ||
      view->display.height != view->image.height) {
    view->display = make_image(
        view->image.width, view->image.height, false, true);
  }
  if (view->image.linear) {
    tonemap_image_mt(view->display, view->image, view->exposure, view->filmic);
  } else if (!view->image.pixelsf.empty() || !view->image.pixelsf.empty()) {
    convert_image(view->display, view->image);
  } else if (!view->image.pixelsb.empty() || !view->image.pixelsb.empty()) {
    convert_image(view->display, view->image);
  } else {
    // TODO(fabio): decide about empty images
  }
  if (!is_initialized(view->glimage)) init_image(view->glimage);
  set_image(view->glimage, view->display, false, false);
}

void draw_widgets(
    gui_window* win, ogl_imageviewer& viewer, const gui_input& input) {
  static string load_path = "", save_path = "", error_message = "";
  if (draw_filedialog_button(win, "load", true, "load image", load_path, false,
          "./", "", "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    // load_image_async(viewer, load_path);
    load_path = "";
  }
  continue_line(win);
  if (draw_filedialog_button(win, "save", viewer.selected, "save image",
          save_path, true, path_dirname(save_path) + "/",
          path_filename(save_path), "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    auto error = string{};
    save_image(save_path, viewer.selected->display, error);
    save_path = "";
  }
  continue_line(win);
  if (draw_button(win, "close", (bool)viewer.selected)) {
    close_image(viewer, viewer.selected->name);
  }
  continue_line(win);
  if (draw_button(win, "quit")) {
    set_close(win, true);
  }
  draw_combobox(win, "image", viewer.selected, viewer.views, false);
  if (!viewer.selected) return;
  if (begin_header(win, "inspect")) {
    auto view = viewer.selected;
    draw_label(win, "name", view->name);
    auto size = vec2i{view->display.width, view->display.height};
    draw_dragger(win, "size", size);
    draw_slider(win, "zoom", view->glparams.scale, 0.1, 10);
    draw_checkbox(win, "fit", view->glparams.fit);
    auto [i, j] = image_coords(input.mouse_pos, view->glparams.center,
        view->glparams.scale, vec2i{view->display.width, view->display.height});
    auto ij     = vec2i{i, j};
    draw_dragger(win, "mouse", ij);
    auto hdr_pixel     = zero4f;
    auto ldr_pixel     = zero4b;
    auto display_pixel = zero4b;
    auto width = view->display.width, height = view->display.height;
    if (i >= 0 && i < width && j >= 0 && j < height) {
      hdr_pixel     = !view->image.pixelsf.empty()
                          ? view->image.pixelsf[j * width + i]
                          : zero4f;
      ldr_pixel     = !view->image.pixelsb.empty()
                          ? view->image.pixelsb[j * width + i]
                          : zero4b;
      display_pixel = view->display.pixelsb[j * width + i];
    }
    if (!view->image.pixelsf.empty()) {
      draw_coloredit(win, "source", hdr_pixel);
    } else {
      draw_coloredit(win, "source", ldr_pixel);
    }
    draw_coloredit(win, "display", display_pixel);
    end_header(win);
  }
  if (!viewer.selected->image.pixelsf.empty()) {
    if (begin_header(win, "tonemap")) {
      auto view   = viewer.selected;
      auto edited = 0;
      edited += draw_slider(win, "exposure", view->exposure, -5, 5);
      edited += draw_checkbox(win, "filmic", view->filmic);
      if (edited) update_display(view);
      end_header(win);
    }
  }
  if (!viewer.selected->params.empty()) {
    if (draw_params(win, viewer.selected->pname, viewer.selected->params)) {
      for (auto pos = (size_t)0; pos < viewer.views.size(); pos++) {
        if (viewer.views[pos].get() == viewer.selected) {
          auto lock                  = std::lock_guard{viewer.input_mutex};
          viewer.inputs[pos]->params = viewer.selected->params;
        }
      }
      if (viewer.pcallback)
        viewer.pcallback(viewer.selected->name, viewer.selected->params);
    }
  }
}

void draw(gui_window* win, ogl_imageviewer& viewer, const gui_input& input) {
  if (!viewer.selected) {
    clear_ogl_framebuffer(ogl_image_params{}.background);
    return;
  }
  auto view                  = viewer.selected;
  view->glparams.window      = input.window_size;
  view->glparams.framebuffer = input.framebuffer_viewport;
  if (!is_initialized(view->glimage)) init_image(view->glimage);
  std::tie(view->glparams.center, view->glparams.scale) = camera_imview(
      view->glparams.center, view->glparams.scale,
      {view->display.width, view->display.height}, view->glparams.window,
      view->glparams.fit);
  draw_image(view->glimage, view->glparams);
}

void update(gui_window* win, ogl_imageviewer& viewer, const gui_input& input) {
  // process inputs
  auto lock = std::lock_guard{viewer.input_mutex};

  // close images
  for (auto idx = (size_t)0; idx < viewer.inputs.size(); idx++) {
    if (!viewer.inputs[idx]->close) continue;
    if (viewer.selected == viewer.views[idx].get()) viewer.selected = nullptr;
    viewer.inputs.erase(viewer.inputs.begin() + idx);
    viewer.views.erase(viewer.views.begin() + idx);
    idx--;
  }

  // add images
  for (auto idx = (size_t)0; idx < viewer.inputs.size(); idx++) {
    if (idx >= viewer.views.size()) {
      viewer.views.emplace_back(std::make_unique<ogl_imageview>());
      viewer.views[idx]->name = viewer.inputs[idx]->name;
    }
  }

  // update images
  for (auto idx = (size_t)0; idx < viewer.inputs.size(); idx++) {
    if (viewer.inputs[idx]->ichanged) {
      viewer.views[idx]->image     = viewer.inputs[idx]->image;
      viewer.inputs[idx]->ichanged = false;
      update_display(viewer.views[idx].get());
    }
    if (viewer.inputs[idx]->pchanged) {
      viewer.views[idx]->params    = viewer.inputs[idx]->params;
      viewer.views[idx]->pname     = viewer.inputs[idx]->pname;
      viewer.inputs[idx]->pchanged = false;
    }
  }

  // selected
  if (viewer.selected == nullptr && !viewer.views.empty())
    viewer.selected = viewer.views[0].get();
}

// Run application
void run_viewer(ogl_imageviewer& viewer) {
  // callbacks
  auto callbacks     = gui_callbacks{};
  callbacks.clear_cb = [&viewer](gui_window* win, const gui_input& input) {
    for (auto& image : viewer.views) clear_image(image->glimage);
  };
  callbacks.update_cb = [&viewer](gui_window* win, const gui_input& input) {
    update(win, viewer, input);
  };
  callbacks.draw_cb = [&viewer](gui_window* win, const gui_input& input) {
    draw(win, viewer, input);
  };
  callbacks.widgets_cb = [&viewer](gui_window* win, const gui_input& input) {
    draw_widgets(win, viewer, input);
  };
  callbacks.uiupdate_cb = [&viewer](gui_window* win, const gui_input& input) {
    if (!viewer.selected) return;
    if (input.widgets_active) return;
    auto view = viewer.selected;
    // handle mouse
    if (input.modifier_alt) {
      if (input.mouse_left) {
        view->glparams.center += input.mouse_pos - input.mouse_last;
      }
      if (input.mouse_right) {
        view->glparams.scale *= powf(
            2, (input.mouse_pos.x - input.mouse_last.x) * 0.001f);
      }
    } else {
      if (viewer.icallback) viewer.icallback(viewer.selected->name, input);
    }
  };

  // run ui
  run_ui({1280 + 320, 720}, "yimview", callbacks);
}

}  // namespace yocto
