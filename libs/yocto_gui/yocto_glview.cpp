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
  auto viewer_guard = make_imageviewer(title);
  auto viewer       = viewer_guard.get();

  // set view
  set_image(viewer, name, img);

  // run view
  run_viewer(viewer);
}
void view_image(
    const string& title, const string& name, const image<vec4b>& img) {
  // open viewer
  auto viewer_guard = make_imageviewer(title);
  auto viewer       = viewer_guard.get();

  // set view
  set_image(viewer, name, img);

  // run view
  run_viewer(viewer);
}

// Open a window and show a shape via path tracing
void view_shape(const string& title, const string& name,
    const generic_shape& shape, bool addsky,
    const progress_callback& progress_cb) {
  // initialize path tracer scene
  auto scene_guard = std::make_unique<sceneio_scene>();
  auto scene       = scene_guard.get();
  print_progress("create scene", 0, 1);
  auto ioshape       = add_shape(scene, "shape");
  ioshape->points    = shape.points;
  ioshape->lines     = shape.lines;
  ioshape->triangles = shape.triangles;
  ioshape->quads     = shape.quads;
  ioshape->positions = shape.positions;
  ioshape->normals   = shape.normals;
  ioshape->texcoords = shape.texcoords;
  ioshape->colors    = shape.colors;
  ioshape->radius    = shape.radius;
  auto instance      = add_instance(scene, "instance");
  instance->shape    = ioshape;
  add_cameras(scene);
  add_materials(scene);
  if (addsky) add_sky(scene);
  print_progress("create scene", 0, 1);

  // run view
  view_scene(title, name, scene, scene->cameras[0]->name, progress_cb);
}

// Construct a scene from io
static void init_scene(trace_scene* scene, const sceneio_scene* ioscene,
    trace_camera*& camera, const sceneio_camera* iocamera,
    progress_callback progress_cb = {}) {
  // handle progress
  auto progress = vec2i{
      0, (int)ioscene->cameras.size() + (int)ioscene->environments.size() +
             (int)ioscene->materials.size() + (int)ioscene->textures.size() +
             (int)ioscene->shapes.size() + (int)ioscene->instances.size()};

  auto camera_map     = unordered_map<const sceneio_camera*, trace_camera*>{};
  camera_map[nullptr] = nullptr;
  for (auto iocamera : ioscene->cameras) {
    if (progress_cb)
      progress_cb("converting cameras", progress.x++, progress.y);
    auto camera          = add_camera(scene);
    camera->frame        = iocamera->frame;
    camera->lens         = iocamera->lens;
    camera->aspect       = iocamera->aspect;
    camera->film         = iocamera->film;
    camera->orthographic = iocamera->orthographic;
    camera->aperture     = iocamera->aperture;
    camera->focus        = iocamera->focus;
    camera_map[iocamera] = camera;
  }

  auto texture_map = unordered_map<const sceneio_texture*, trace_texture*>{};
  texture_map[nullptr] = nullptr;
  for (auto iotexture : ioscene->textures) {
    if (progress_cb)
      progress_cb("converting textures", progress.x++, progress.y);
    auto texture           = add_texture(scene);
    texture->hdr           = iotexture->hdr;
    texture->ldr           = iotexture->ldr;
    texture_map[iotexture] = texture;
  }

  auto material_map = unordered_map<const sceneio_material*, trace_material*>{};
  material_map[nullptr] = nullptr;
  for (auto iomaterial : ioscene->materials) {
    if (progress_cb)
      progress_cb("converting materials", progress.x++, progress.y);
    auto material              = add_material(scene);
    material->emission         = iomaterial->emission;
    material->color            = iomaterial->color;
    material->specular         = iomaterial->specular;
    material->roughness        = iomaterial->roughness;
    material->metallic         = iomaterial->metallic;
    material->ior              = iomaterial->ior;
    material->spectint         = iomaterial->spectint;
    material->coat             = iomaterial->coat;
    material->transmission     = iomaterial->transmission;
    material->translucency     = iomaterial->translucency;
    material->scattering       = iomaterial->scattering;
    material->scanisotropy     = iomaterial->scanisotropy;
    material->trdepth          = iomaterial->trdepth;
    material->opacity          = iomaterial->opacity;
    material->thin             = iomaterial->thin;
    material->emission_tex     = texture_map.at(iomaterial->emission_tex);
    material->color_tex        = texture_map.at(iomaterial->color_tex);
    material->specular_tex     = texture_map.at(iomaterial->specular_tex);
    material->metallic_tex     = texture_map.at(iomaterial->metallic_tex);
    material->roughness_tex    = texture_map.at(iomaterial->roughness_tex);
    material->transmission_tex = texture_map.at(iomaterial->transmission_tex);
    material->translucency_tex = texture_map.at(iomaterial->translucency_tex);
    material->spectint_tex     = texture_map.at(iomaterial->spectint_tex);
    material->scattering_tex   = texture_map.at(iomaterial->scattering_tex);
    material->coat_tex         = texture_map.at(iomaterial->coat_tex);
    material->opacity_tex      = texture_map.at(iomaterial->opacity_tex);
    material->normal_tex       = texture_map.at(iomaterial->normal_tex);
    material_map[iomaterial]   = material;
  }

  auto shape_map     = unordered_map<const sceneio_shape*, trace_shape*>{};
  shape_map[nullptr] = nullptr;
  for (auto ioshape : ioscene->shapes) {
    if (progress_cb) progress_cb("converting shapes", progress.x++, progress.y);
    auto shape              = add_shape(scene);
    shape->points           = ioshape->points;
    shape->lines            = ioshape->lines;
    shape->triangles        = ioshape->triangles;
    shape->quads            = ioshape->quads;
    shape->quadspos         = ioshape->quadspos;
    shape->quadsnorm        = ioshape->quadsnorm;
    shape->quadstexcoord    = ioshape->quadstexcoord;
    shape->positions        = ioshape->positions;
    shape->normals          = ioshape->normals;
    shape->texcoords        = ioshape->texcoords;
    shape->colors           = ioshape->colors;
    shape->radius           = ioshape->radius;
    shape->tangents         = ioshape->tangents;
    shape->subdivisions     = ioshape->subdivisions;
    shape->catmullclark     = ioshape->catmullclark;
    shape->smooth           = ioshape->smooth;
    shape->displacement     = ioshape->displacement;
    shape->displacement_tex = texture_map.at(ioshape->displacement_tex);
    shape_map[ioshape]      = shape;
  }

  for (auto ioinstance : ioscene->instances) {
    if (progress_cb)
      progress_cb("converting instances", progress.x++, progress.y);
    auto instance      = add_instance(scene);
    instance->frame    = ioinstance->frame;
    instance->shape    = shape_map.at(ioinstance->shape);
    instance->material = material_map.at(ioinstance->material);
  }

  for (auto ioenvironment : ioscene->environments) {
    if (progress_cb)
      progress_cb("converting environments", progress.x++, progress.y);
    auto environment          = add_environment(scene);
    environment->frame        = ioenvironment->frame;
    environment->emission     = ioenvironment->emission;
    environment->emission_tex = texture_map.at(ioenvironment->emission_tex);
  }

  // done
  if (progress_cb) progress_cb("converting done", progress.x++, progress.y);

  // get camera
  camera = camera_map.at(iocamera);
}

// Open a window and show an scene via path tracing
void view_scene(const string& title, const string& name,
    const sceneio_scene* scene, const string& camera_,
    const progress_callback& progress_cb) {
  // get camera
  auto camera = get_camera(scene, camera_);

  // scene conversion
  auto ptscene_guard = std::make_unique<trace_scene>();
  auto ptscene       = ptscene_guard.get();
  auto ptcamera      = (trace_camera*)nullptr;
  init_scene(ptscene, scene, ptcamera, camera, progress_cb);

  // rendering params
  auto params = trace_params{};
  auto has_lights =
      std::any_of(scene->instances.begin(), scene->instances.end(),
          [](sceneio_instance* instance) {
            return instance->material->emission != zero3f;
          }) ||
      std::any_of(scene->environments.begin(), scene->environments.end(),
          [](sceneio_environment* environment) {
            return environment->emission != zero3f;
          });
  if (!has_lights) params.sampler = trace_sampler_type::eyelight;

  // run viewer
  view_scene(title, name, ptscene, ptcamera, params, progress_cb);
}

// Open a window and show an scene via path tracing
void view_scene(const string& title, const string& name,
    const trace_scene* scene, const trace_camera* camera_,
    const trace_params& params_, const progress_callback& progress_cb) {
  // open viewer
  auto viewer_guard = make_imageviewer(title);
  auto viewer       = viewer_guard.get();

  // copy params and camera
  auto params       = params_;
  auto camera_guard = make_unique<trace_camera>();
  auto camera       = camera_guard.get();
  *camera           = *camera_;

  // build bvh
  auto bvh_guard = std::make_unique<trace_bvh>();
  auto bvh       = bvh_guard.get();
  init_bvh(bvh, scene, params, print_progress);

  // init renderer
  auto lights_guard = std::make_unique<trace_lights>();
  auto lights       = lights_guard.get();
  init_lights(lights, scene, params, print_progress);

  // fix renderer type if no lights
  if (lights->lights.empty() && is_sampler_lit(params)) {
    print_info("no lights presents, image will be black");
  }

  // init state
  auto state_guard = std::make_unique<trace_state>();
  auto state       = state_guard.get();

  // render start
  trace_start(
      state, scene, camera, bvh, lights, params,
      [viewer, name](const string& message, int sample, int nsamples) {
        set_widget(viewer, name, "sample", to_json(sample),
            to_schema(sample, "Current sample"));
        print_progress(message, sample, nsamples);
      },
      [viewer, name](const image<vec4f>& render, int current, int total) {
        set_image(viewer, name, render);
      });

  // show rendering params
  set_widgets(
      viewer, name, to_json(params), to_schema(params, "Render params"));

  // set callback
  set_callback(viewer,
      [state, scene, camera, bvh, lights, viewer, &params](const string& name,
          const json_value& uiparams, const gui_input& input) {
        if (name != name) return;
        if (!uiparams.is_null()) {
          trace_stop(state);
          params = from_json<trace_params>(uiparams);
          // show rendering params
          set_widgets(viewer, name, to_json(params),
              to_schema(params, "Render params"));
          trace_start(
              state, scene, camera, bvh, lights, params,
              [viewer, name](const string& message, int sample, int nsamples) {
                set_widget(viewer, name, "sample", to_json(sample),
                    to_schema(sample, "Current sample"));
                print_progress(message, sample, nsamples);
              },
              [viewer, name](const image<vec4f>& render, int current,
                  int total) { set_image(viewer, name, render); });
        } else if ((input.mouse_left || input.mouse_right) &&
                   input.mouse_pos != input.mouse_last) {
          trace_stop(state);
          auto dolly  = 0.0f;
          auto pan    = zero2f;
          auto rotate = zero2f;
          if (input.mouse_left && !input.modifier_shift)
            rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
          if (input.mouse_right)
            dolly = (input.mouse_pos.x - input.mouse_last.x) / 100.0f;
          if (input.mouse_left && input.modifier_shift)
            pan = (input.mouse_pos - input.mouse_last) * camera->focus / 200.0f;
          pan.x                                  = -pan.x;
          std::tie(camera->frame, camera->focus) = camera_turntable(
              camera->frame, camera->focus, rotate, dolly, pan);
          trace_start(
              state, scene, camera, bvh, lights, params,
              [viewer, name](const string& message, int sample, int nsamples) {
                set_widget(viewer, name, "sample", to_json(sample),
                    to_schema(sample, "Current sample"));
                print_progress(message, sample, nsamples);
              },
              [viewer, name](const image<vec4f>& render, int current,
                  int total) { set_image(viewer, name, render); });
        }
      });

  // run view
  run_viewer(viewer);

  // stop
  trace_stop(state);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE VIEWER API
// -----------------------------------------------------------------------------
namespace yocto {

// grab input
// static imageview_image* get_image(imageview_state* viewer, const string&
// name);
static ogl_imageinput* get_input(ogl_imageviewer* viewer, const string& name);

// make an image viewer
unique_ptr<ogl_imageviewer> make_imageviewer(const string& title) {
  auto viewer = make_unique<ogl_imageviewer>();
  // viewer->name = title;
  return viewer;
}

// Set image
void set_image(ogl_imageviewer* viewer, const string& name,
    const image<vec4f>& img, float exposure, bool filmic) {
  auto lock  = std::lock_guard{viewer->input_mutex};
  auto input = get_input(viewer, name);
  if (!input) {
    input =
        viewer->inputs.emplace_back(std::make_unique<ogl_imageinput>()).get();
  }
  input->name     = name;
  input->hdr      = img;
  input->ldr      = {};
  input->exposure = exposure;
  input->filmic   = filmic;
  input->ichanged = true;
}
void set_image(
    ogl_imageviewer* viewer, const string& name, const image<vec4b>& img) {
  auto lock  = std::lock_guard{viewer->input_mutex};
  auto input = get_input(viewer, name);
  if (!input) {
    input =
        viewer->inputs.emplace_back(std::make_unique<ogl_imageinput>()).get();
  }
  input->name     = name;
  input->hdr      = {};
  input->ldr      = img;
  input->exposure = 0;
  input->filmic   = false;
  input->ichanged = true;
}
// Close image
void close_image(ogl_imageviewer* viewer, const string& name) {
  auto lock  = std::lock_guard{viewer->input_mutex};
  auto input = get_input(viewer, name);
  if (!input) return;
  input->close = true;
}

// Set params
void set_widget(ogl_imageviewer* viewer, const string& name,
    const string& pname, const json_value& param, const json_value& schema) {
  auto lock  = std::lock_guard{viewer->input_mutex};
  auto input = get_input(viewer, name);
  if (!input) return;
  input->widgets[pname]              = param;
  input->schema["properties"][pname] = schema;
  input->wchanged                    = true;
}
void set_widgets(ogl_imageviewer* viewer, const string& name,
    const json_value& params, const json_value& schema) {
  auto lock  = std::lock_guard{viewer->input_mutex};
  auto input = get_input(viewer, name);
  if (!input) return;
  input->widgets  = params;
  input->schema   = schema;
  input->wchanged = true;
}

// Callback
void set_callback(
    ogl_imageviewer* viewer, const ogl_imageviewer_callback& callback) {
  viewer->callback = callback;
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
//   for (auto& img : viewer->images)
//     if (view->name == name) return img.get();
//   return nullptr;
// }
static ogl_imageinput* get_input(ogl_imageviewer* viewer, const string& name) {
  for (auto& input : viewer->inputs)
    if (input->name == name) return input.get();
  return nullptr;
}

static void update_display(ogl_imageview* view) {
  if (!view->hdr.empty()) {
    if (view->display.imsize() != view->hdr.imsize())
      view->display.resize(view->hdr.imsize());
    tonemap_image_mt(view->display, view->hdr, view->exposure, view->filmic);
  } else if (!view->ldr.empty()) {
    view->display = view->ldr;
  } else {
    // TODO(fabio): decide about empty images
  }
  if (!is_initialized(view->glimage)) init_image(view->glimage);
  set_image(view->glimage, view->display, false, false);
}

void draw_widgets(
    gui_window* win, ogl_imageviewer* viewer, const gui_input& input) {
  static string load_path = "", save_path = "", error_message = "";
  if (draw_filedialog_button(win, "load", true, "load image", load_path, false,
          "./", "", "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    // load_image_async(viewer, load_path);
    load_path = "";
  }
  continue_line(win);
  if (draw_filedialog_button(win, "save", viewer->selected, "save image",
          save_path, true, path_dirname(save_path), path_filename(save_path),
          "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    // viewer->selected->outname = save_path;
    // save_image(view->outname, view->display, view->error);
    save_path = "";
  }
  continue_line(win);
  if (draw_button(win, "close", (bool)viewer->selected)) {
    close_image(viewer, viewer->selected->name);
  }
  continue_line(win);
  if (draw_button(win, "quit")) {
    set_close(win, true);
  }
  draw_combobox(win, "image", viewer->selected, viewer->views, false);
  if (!viewer->selected) return;
  if (begin_header(win, "inspect")) {
    auto view = viewer->selected;
    draw_label(win, "name", view->name);
    auto size = view->display.imsize();
    draw_dragger(win, "size", size);
    draw_slider(win, "zoom", view->glparams.scale, 0.1, 10);
    draw_checkbox(win, "fit", view->glparams.fit);
    auto ij = image_coords(input.mouse_pos, view->glparams.center,
        view->glparams.scale, view->display.imsize());
    draw_dragger(win, "mouse", ij);
    auto hdr_pixel     = zero4f;
    auto ldr_pixel     = zero4b;
    auto display_pixel = zero4b;
    if (ij.x >= 0 && ij.x < view->display.width() && ij.y >= 0 &&
        ij.y < view->display.height()) {
      hdr_pixel     = !view->hdr.empty() ? view->hdr[{ij.x, ij.y}] : zero4f;
      ldr_pixel     = !view->ldr.empty() ? view->ldr[{ij.x, ij.y}] : zero4b;
      display_pixel = view->display[{ij.x, ij.y}];
    }
    if (!view->hdr.empty()) {
      draw_coloredit(win, "source", hdr_pixel);
    } else {
      draw_coloredit(win, "source", ldr_pixel);
    }
    draw_coloredit(win, "display", display_pixel);
    end_header(win);
  }
  if (!viewer->selected->hdr.empty()) {
    if (begin_header(win, "tonemap")) {
      auto view   = viewer->selected;
      auto edited = 0;
      edited += draw_slider(win, "exposure", view->exposure, -5, 5);
      edited += draw_checkbox(win, "filmic", view->filmic);
      if (edited) update_display(view);
      end_header(win);
    }
  }
  if (!viewer->selected->widgets.empty()) {
    if (draw_params(win, "params", viewer->selected->widgets,
            viewer->selected->schema, false)) {
      if (viewer->callback)
        viewer->callback(viewer->selected->name, viewer->selected->widgets, {});
    }
  }
}

void draw(gui_window* win, ogl_imageviewer* viewer, const gui_input& input) {
  if (!viewer->selected) {
    clear_ogl_framebuffer(ogl_image_params{}.background);
    return;
  }
  auto view                  = viewer->selected;
  view->glparams.window      = input.window_size;
  view->glparams.framebuffer = input.framebuffer_viewport;
  if (!is_initialized(view->glimage)) init_image(view->glimage);
  std::tie(view->glparams.center, view->glparams.scale) = camera_imview(
      view->glparams.center, view->glparams.scale, view->display.imsize(),
      view->glparams.window, view->glparams.fit);
  draw_image(view->glimage, view->glparams);
}

void update(gui_window* win, ogl_imageviewer* viewer, const gui_input& input) {
  // process inputs
  auto lock = std::lock_guard{viewer->input_mutex};

  // close images
  for (auto idx = (size_t)0; idx < viewer->inputs.size(); idx++) {
    if (!viewer->inputs[idx]->close) continue;
    if (viewer->selected == viewer->views[idx].get())
      viewer->selected = nullptr;
    viewer->inputs.erase(viewer->inputs.begin() + idx);
    viewer->views.erase(viewer->views.begin() + idx);
    idx--;
  }

  // add images
  for (auto idx = (size_t)0; idx < viewer->inputs.size(); idx++) {
    if (idx >= viewer->views.size()) {
      viewer->views.emplace_back(std::make_unique<ogl_imageview>());
      viewer->views[idx]->name = viewer->inputs[idx]->name;
    }
  }

  // update images
  for (auto idx = (size_t)0; idx < viewer->inputs.size(); idx++) {
    if (viewer->inputs[idx]->ichanged) {
      viewer->views[idx]->hdr       = viewer->inputs[idx]->hdr;
      viewer->views[idx]->ldr       = viewer->inputs[idx]->ldr;
      viewer->inputs[idx]->ichanged = false;
      update_display(viewer->views[idx].get());
    }
    if (viewer->inputs[idx]->wchanged) {
      viewer->views[idx]->widgets   = viewer->inputs[idx]->widgets;
      viewer->views[idx]->schema    = viewer->inputs[idx]->schema;
      viewer->inputs[idx]->wchanged = false;
    }
  }

  // selected
  if (viewer->selected == nullptr && !viewer->views.empty())
    viewer->selected = viewer->views[0].get();
}

// Run application
void run_viewer(ogl_imageviewer* viewer) {
  // callbacks
  auto callbacks     = gui_callbacks{};
  callbacks.clear_cb = [viewer](gui_window* win, const gui_input& input) {
    for (auto& image : viewer->views) clear_image(image->glimage);
  };
  callbacks.update_cb = [viewer](gui_window* win, const gui_input& input) {
    update(win, viewer, input);
  };
  callbacks.draw_cb = [viewer](gui_window* win, const gui_input& input) {
    draw(win, viewer, input);
  };
  callbacks.widgets_cb = [viewer](gui_window* win, const gui_input& input) {
    draw_widgets(win, viewer, input);
  };
  callbacks.uiupdate_cb = [viewer](gui_window* win, const gui_input& input) {
    if (!viewer->selected) return;
    if (input.widgets_active) return;
    auto view = viewer->selected;
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
      if (viewer->callback) viewer->callback(viewer->selected->name, {}, input);
    }
  };

  // run ui
  run_ui({1280 + 320, 720}, "yimview", callbacks);
}

}  // namespace yocto
