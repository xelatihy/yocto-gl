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

#include <yocto/yocto_commonio.h>
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
void view_image(
    const string& title, const string& name, const image<vec4f>& img) {
  // open viewer
  auto viewer = make_imageviewer(title);

  // set view
  set_image(viewer, name, img);

  // run view
  run_viewer(viewer);
}
void view_image(
    const string& title, const string& name, const image<vec4b>& img) {
  // open viewer
  auto viewer = make_imageviewer(title);

  // set view
  set_image(viewer, name, img);

  // run view
  run_viewer(viewer);
}
void view_image(
    const string& title, const string& name, const image_data& img) {
  // open viewer
  auto viewer = make_imageviewer(title);

  // set view
  set_image(viewer, name, img);

  // run view
  run_viewer(viewer);
}

// Open a window and show a shape via path tracing
void view_shape(const string& title, const string& name,
    const shape_data& shape, bool addsky,
    const progress_callback& progress_cb) {
  // initialize path tracer scene
  auto scene = scene_scene{};
  print_progress("create scene", 0, 1);
  scene.shape_names.emplace_back("shape");
  auto& ioshape     = scene.shapes.emplace_back();
  ioshape.points    = shape.points;
  ioshape.lines     = shape.lines;
  ioshape.triangles = shape.triangles;
  ioshape.quads     = shape.quads;
  ioshape.positions = shape.positions;
  ioshape.normals   = shape.normals;
  ioshape.texcoords = shape.texcoords;
  ioshape.colors    = shape.colors;
  ioshape.radius    = shape.radius;
  scene.instance_names.emplace_back("instance");
  auto& instance = scene.instances.emplace_back();
  instance.shape = 0;
  add_cameras(scene);
  add_materials(scene);
  if (addsky) add_sky(scene);
  print_progress("create scene", 0, 1);

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
  auto bvh = make_bvh(scene, params, print_progress);

  // init renderer
  auto lights = make_lights(scene, params, print_progress);

  // fix renderer type if no lights
  if (lights.lights.empty() && is_sampler_lit(params)) {
    print_info("no lights presents, image will be black");
  }

  // init state
  auto worker = trace_worker{};
  auto state  = trace_state{};

  // render start
  trace_start(
      worker, state, scene, bvh, lights, params,
      [&viewer, name](const string& message, int sample, int nsamples) {
        set_param(viewer, name, "sample", {sample, {1, 4096}, true});
        print_progress(message, sample, nsamples);
      },
      [&viewer, name](const image<vec4f>& render, int current, int total) {
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
            [&viewer, name](const string& message, int sample, int nsamples) {
              set_param(viewer, name, "sample", {sample, {1, 4096}, true});
              print_progress(message, sample, nsamples);
            },
            [&viewer, name](const image<vec4f>& render, int current,
                int total) { set_image(viewer, name, render); });
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
          [&viewer, name](const string& message, int sample, int nsamples) {
            set_param(viewer, name, "sample", {sample, {1, 4096}, true});
            print_progress(message, sample, nsamples);
          },
          [&viewer, name](const image<vec4f>& render, int current, int total) {
            set_image(viewer, name, render);
          });
    }
  });

  // run view
  run_viewer(viewer);

  // stop
  trace_stop(worker);
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
void set_image(ogl_imageviewer& viewer, const string& name,
    const image<vec4f>& img, float exposure, bool filmic) {
  auto lock  = std::lock_guard{viewer.input_mutex};
  auto input = get_input(viewer, name);
  if (!input) {
    input =
        viewer.inputs.emplace_back(std::make_unique<ogl_imageinput>()).get();
  }
  input->name     = name;
  input->image    = make_image(img.width(), img.height(), true, img.data());
  input->exposure = exposure;
  input->filmic   = filmic;
  input->ichanged = true;
}
void set_image(
    ogl_imageviewer& viewer, const string& name, const image<vec4b>& img) {
  auto lock  = std::lock_guard{viewer.input_mutex};
  auto input = get_input(viewer, name);
  if (!input) {
    input =
        viewer.inputs.emplace_back(std::make_unique<ogl_imageinput>()).get();
  }
  input->name     = name;
  input->image    = make_image(img.width(), img.height(), false, img.data());
  input->exposure = 0;
  input->filmic   = false;
  input->ichanged = true;
}
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
    auto error = ""s;
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
