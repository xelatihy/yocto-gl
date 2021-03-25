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
#include <yocto/yocto_parallel.h>

#include <cassert>
#include <stdexcept>

#include "ext/glad/glad.h"

#ifdef _WIN32
#undef near
#undef far
#endif

// -----------------------------------------------------------------------------
// VIEW HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

static void update_image_params(gui_window* win, const gui_input& input,
    const color_image& image, glimage_params& glparams) {
  glparams.window                           = input.window_size;
  glparams.framebuffer                      = input.framebuffer_viewport;
  std::tie(glparams.center, glparams.scale) = camera_imview(glparams.center,
      glparams.scale, {image.width, image.height}, glparams.window,
      glparams.fit);
}

static bool uiupdate_image_params(
    gui_window* win, const gui_input& input, glimage_params& glparams) {
  // handle mouse
  if (input.mouse_left && !input.widgets_active) {
    glparams.center += input.mouse_pos - input.mouse_last;
    return true;
  }
  if (input.mouse_right && !input.widgets_active) {
    glparams.scale *= pow(
        2.0f, (input.mouse_pos.x - input.mouse_last.x) * 0.001f);
    return true;
  }
  return false;
}

static bool uiupdate_camera_params(
    gui_window* win, const gui_input& input, scene_camera& camera) {
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
      pan = (input.mouse_pos - input.mouse_last) * camera.focus / 200.0f;
    pan.x               = -pan.x;
    auto [frame, focus] = camera_turntable(
        camera.frame, camera.focus, rotate, dolly, pan);
    if (camera.frame != frame || camera.focus != focus) {
      camera.frame = frame;
      camera.focus = focus;
      return true;
    }
  }
  return false;
}

static bool draw_tonemap_params(
    gui_window* win, const gui_input& input, float& exposure, bool& filmic) {
  auto edited = 0;
  if (begin_header(win, "tonemap")) {
    edited += draw_slider(win, "exposure", exposure, -5, 5);
    edited += draw_checkbox(win, "filmic", filmic);
    end_header(win);
  }
  return (bool)edited;
}

static bool draw_image_inspector(gui_window* win, const gui_input& input,
    const color_image& image, const color_image& display,
    glimage_params& glparams) {
  if (begin_header(win, "inspect")) {
    draw_slider(win, "zoom", glparams.scale, 0.1, 10);
    draw_checkbox(win, "fit", glparams.fit);
    auto [i, j] = image_coords(input.mouse_pos, glparams.center, glparams.scale,
        {image.width, image.height});
    auto ij     = vec2i{i, j};
    draw_dragger(win, "mouse", ij);
    auto image_pixel   = zero4f;
    auto display_pixel = zero4f;
    if (i >= 0 && i < image.width && j >= 0 && j < image.height) {
      image_pixel   = image.pixels[j * image.width + i];
      display_pixel = image.pixels[j * image.width + i];
    }
    draw_coloredit(win, "image", image_pixel);
    draw_coloredit(win, "display", display_pixel);
    end_header(win);
  }
  return false;
}

struct scene_selection {
  int camera      = 0;
  int instance    = 0;
  int environment = 0;
  int shape       = 0;
  int texture     = 0;
  int material    = 0;
  int subdiv      = 0;
};

static bool draw_scene_editor(gui_window* win, scene_model& scene,
    scene_selection& selection, const function<void()>& before_edit) {
  auto edited = 0;
  if (begin_header(win, "cameras")) {
    draw_combobox(win, "camera", selection.camera, scene.camera_names);
    auto camera = scene.cameras.at(selection.camera);
    edited += draw_checkbox(win, "ortho", camera.orthographic);
    edited += draw_slider(win, "lens", camera.lens, 0.001, 1);
    edited += draw_slider(win, "aspect", camera.aspect, 0.1, 5);
    edited += draw_slider(win, "film", camera.film, 0.1, 0.5);
    edited += draw_slider(win, "focus", camera.focus, 0.001, 100);
    edited += draw_slider(win, "aperture", camera.aperture, 0, 1);
    //   frame3f frame        = identity3x4f;
    if (edited) {
      if (before_edit) before_edit();
      scene.cameras.at(selection.camera) = camera;
    }
    end_header(win);
  }
  if (begin_header(win, "environments")) {
    draw_combobox(
        win, "environment", selection.environment, scene.environment_names);
    auto environment = scene.environments.at(selection.environment);
    edited += draw_hdrcoloredit(win, "emission", environment.emission);
    edited += draw_combobox(win, "emission_tex", environment.emission_tex,
        scene.texture_names, true);
    //   frame3f frame        = identity3x4f;
    if (edited) {
      if (before_edit) before_edit();
      scene.environments.at(selection.environment) = environment;
    }
    end_header(win);
  }
  if (begin_header(win, "instances")) {
    draw_combobox(win, "instance", selection.instance, scene.instance_names);
    auto instance = scene.instances.at(selection.instance);
    edited += draw_combobox(win, "shape", instance.shape, scene.shape_names);
    edited += draw_combobox(
        win, "material", instance.material, scene.material_names);
    //   frame3f frame        = identity3x4f;
    if (edited) {
      if (before_edit) before_edit();
      scene.instances.at(selection.instance) = instance;
    }
    end_header(win);
  }
  if (begin_header(win, "materials")) {
    draw_combobox(win, "material", selection.material, scene.material_names);
    auto material = scene.materials.at(selection.material);
    edited += draw_hdrcoloredit(win, "emission", material.emission);
    edited += draw_combobox(
        win, "emission_tex", material.emission_tex, scene.texture_names, true);
    edited += draw_hdrcoloredit(win, "color", material.color);
    edited += draw_combobox(
        win, "color_tex", material.color_tex, scene.texture_names, true);
    edited += draw_slider(win, "roughness", material.roughness, 0, 1);
    edited += draw_combobox(win, "roughness_tex", material.roughness_tex,
        scene.texture_names, true);
    edited += draw_slider(win, "metallic", material.metallic, 0, 1);
    edited += draw_slider(win, "ior", material.ior, 0.1, 5);
    if (edited) {
      if (before_edit) before_edit();
      scene.materials.at(selection.material) = material;
    }
    end_header(win);
  }
  if (begin_header(win, "shapes")) {
    draw_combobox(win, "shape", selection.shape, scene.shape_names);
    auto& shape = scene.shapes.at(selection.shape);
    draw_label(win, "points", (int)shape.points.size());
    draw_label(win, "lines", (int)shape.lines.size());
    draw_label(win, "triangles", (int)shape.triangles.size());
    draw_label(win, "quads", (int)shape.quads.size());
    draw_label(win, "positions", (int)shape.positions.size());
    draw_label(win, "normals", (int)shape.normals.size());
    draw_label(win, "texcoords", (int)shape.texcoords.size());
    draw_label(win, "colors", (int)shape.colors.size());
    draw_label(win, "radius", (int)shape.radius.size());
    draw_label(win, "tangents", (int)shape.tangents.size());
    end_header(win);
  }
  if (begin_header(win, "textures")) {
    draw_combobox(win, "texture", selection.texture, scene.texture_names);
    auto& texture = scene.textures.at(selection.texture);
    draw_label(win, "width", texture.width);
    draw_label(win, "height", texture.height);
    draw_label(win, "linear", texture.linear);
    draw_label(win, "byte", !texture.pixelsb.empty());
    end_header(win);
  }
  if (begin_header(win, "subdivs")) {
    draw_combobox(win, "subdiv", selection.subdiv, scene.subdiv_names);
    auto& subdiv = scene.subdivs.at(selection.subdiv);
    draw_label(win, "quadspos", (int)subdiv.quadspos.size());
    draw_label(win, "quadsnorm", (int)subdiv.quadsnorm.size());
    draw_label(win, "quadstexcoord", (int)subdiv.quadstexcoord.size());
    draw_label(win, "positions", (int)subdiv.positions.size());
    draw_label(win, "normals", (int)subdiv.normals.size());
    draw_label(win, "texcoords", (int)subdiv.texcoords.size());
    end_header(win);
  }
  return (bool)edited;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE AND TRACE VIEW
// -----------------------------------------------------------------------------
namespace yocto {

// Open a window and show an image
void view_image(
    const string& title, const string& name, const color_image& image) {
  // display image
  auto  display  = make_image(image.width, image.height, false);
  float exposure = 0;
  bool  filmic   = false;
  tonemap_image_mt(display, image, exposure, filmic);

  // opengl image
  auto glimage  = glimage_state{};
  auto glparams = glimage_params{};

  // top level combo
  auto names    = vector<string>{name};
  auto selected = 0;

  // callbacks
  auto callbacks    = gui_callbacks{};
  callbacks.init_cb = [&](gui_window* win, const gui_input& input) {
    init_image(glimage);
    set_image(glimage, display);
  };
  callbacks.clear_cb = [&](gui_window* win, const gui_input& input) {
    clear_image(glimage);
  };
  callbacks.draw_cb = [&](gui_window* win, const gui_input& input) {
    update_image_params(win, input, image, glparams);
    draw_image(glimage, glparams);
  };
  callbacks.widgets_cb = [&](gui_window* win, const gui_input& input) {
    draw_combobox(win, "name", selected, names);
    if (draw_tonemap_params(win, input, exposure, filmic)) {
      tonemap_image_mt(display, image, exposure, filmic);
      set_image(glimage, display);
    }
    draw_image_inspector(win, input, image, display, glparams);
  };
  callbacks.uiupdate_cb = [&](gui_window* win, const gui_input& input) {
    uiupdate_image_params(win, input, glparams);
  };

  // run ui
  run_ui({1280 + 320, 720}, title, callbacks);
}

// Open a window and show an image
void view_images(const string& title, const vector<string>& names,
    const vector<color_image>& images) {
  // display image
  auto displays  = vector<color_image>(images.size());
  auto exposures = vector<float>(images.size(), 0);
  auto filmics   = vector<bool>(images.size(), false);
  for (auto idx = 0; idx < (int)images.size(); idx++) {
    displays[idx] = make_image(images[idx].width, images[idx].height, false);
    tonemap_image_mt(displays[idx], images[idx], exposures[idx], filmics[idx]);
  }

  // opengl image
  auto glimages  = vector<glimage_state>(images.size());
  auto glparamss = vector<glimage_params>(images.size());

  // selection
  auto selected = 0;

  // callbacks
  auto callbacks    = gui_callbacks{};
  callbacks.init_cb = [&](gui_window* win, const gui_input& input) {
    for (auto idx = 0; idx < (int)images.size(); idx++) {
      init_image(glimages[idx]);
      set_image(glimages[idx], displays[idx]);
    }
  };
  callbacks.clear_cb = [&](gui_window* win, const gui_input& input) {
    for (auto idx = 0; idx < (int)images.size(); idx++) {
      clear_image(glimages[idx]);
    }
  };
  callbacks.draw_cb = [&](gui_window* win, const gui_input& input) {
    update_image_params(win, input, displays[selected], glparamss[selected]);
    draw_image(glimages[selected], glparamss[selected]);
  };
  callbacks.widgets_cb = [&](gui_window* win, const gui_input& input) {
    draw_combobox(win, "name", selected, names);
    auto filmic = (bool)filmics[selected];  // vector of bool ...
    if (draw_tonemap_params(win, input, exposures[selected], filmic)) {
      filmics[selected] = filmic;
      tonemap_image_mt(displays[selected], images[selected],
          exposures[selected], filmics[selected]);
      set_image(glimages[selected], displays[selected]);
    }
    draw_image_inspector(
        win, input, images[selected], displays[selected], glparamss[selected]);
  };
  callbacks.uiupdate_cb = [&](gui_window* win, const gui_input& input) {
    uiupdate_image_params(win, input, glparamss[selected]);
  };

  // run ui
  run_ui({1280 + 320, 720}, title, callbacks);
}

// Open a window and show an image
void colorgrade_image(
    const string& title, const string& name, const color_image& image) {
  // color grading parameters
  auto params = colorgrade_params{};

  // display image
  auto display = make_image(image.width, image.height, false);
  colorgrade_image_mt(display, image, params);

  // opengl image
  auto glimage  = glimage_state{};
  auto glparams = glimage_params{};

  // top level combo
  auto names    = vector<string>{name};
  auto selected = 0;

  // callbacks
  auto callbacks    = gui_callbacks{};
  callbacks.init_cb = [&](gui_window* win, const gui_input& input) {
    init_image(glimage);
    set_image(glimage, display);
  };
  callbacks.clear_cb = [&](gui_window* win, const gui_input& input) {
    clear_image(glimage);
  };
  callbacks.draw_cb = [&](gui_window* win, const gui_input& input) {
    update_image_params(win, input, image, glparams);
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
        colorgrade_image_mt(display, image, params);
        set_image(glimage, display);
      }
    }
    draw_image_inspector(win, input, image, display, glparams);
  };
  callbacks.uiupdate_cb = [&glparams](gui_window* win, const gui_input& input) {
    uiupdate_image_params(win, input, glparams);
  };

  // run ui
  run_ui({1280 + 320, 720}, title, callbacks);
}

// Open a window and show an scene via path tracing
void view_scene(const string& title, const string& name, scene_model& scene,
    const trace_params& params_, bool print, bool edit) {
  // copy params and camera
  auto params = params_;

  // build bvh
  if (print) print_progress_begin("build bvh");
  auto bvh = make_bvh(scene, params);
  if (print) print_progress_end();

  // init renderer
  if (print) print_progress_begin("init lights");
  auto lights = make_lights(scene, params);
  if (print) print_progress_end();

  // fix renderer type if no lights
  if (lights.lights.empty() && is_sampler_lit(params)) {
    if (print) print_info("no lights presents --- switching to eyelight");
    params.sampler = trace_sampler_type::eyelight;
  }

  // init state
  if (print) print_progress_begin("init state");
  auto state   = make_state(scene, params);
  auto image   = make_image(state.width, state.height, true);
  auto display = make_image(state.width, state.height, false);
  auto render  = make_image(state.width, state.height, true);
  if (print) print_progress_end();

  // opengl image
  auto glimage  = glimage_state{};
  auto glparams = glimage_params{};

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

  // renderer update
  auto render_update  = std::atomic<bool>{};
  auto render_current = std::atomic<int>{};
  auto render_mutex   = std::mutex{};
  auto render_worker  = future<void>{};
  auto render_stop    = atomic<bool>{};
  auto reset_display  = [&]() {
    // stop render
    render_stop = true;
    if (render_worker.valid()) render_worker.get();

    state   = make_state(scene, params);
    image   = make_image(state.width, state.height, true);
    display = make_image(state.width, state.height, false);
    render  = make_image(state.width, state.height, true);

    render_worker = {};
    render_stop   = false;

    // preview
    auto pparams = params;
    pparams.resolution /= params.pratio;
    pparams.samples = 1;
    auto pstate     = make_state(scene, pparams);
    trace_samples(pstate, scene, bvh, lights, pparams);
    auto preview = get_render(pstate);
    for (auto idx = 0; idx < state.width * state.height; idx++) {
      auto i = idx % render.width, j = idx / render.width;
      auto pi            = clamp(i / params.pratio, 0, preview.width - 1),
           pj            = clamp(j / params.pratio, 0, preview.height - 1);
      render.pixels[idx] = preview.pixels[pj * preview.width + pi];
    }
    // if (current > 0) return;
    {
      auto lock      = std::lock_guard{render_mutex};
      render_current = 0;
      image          = render;
      tonemap_image_mt(display, image, params.exposure, params.filmic);
      render_update = true;
    }

    // start renderer
    render_worker = std::async(std::launch::async, [&]() {
      for (auto sample = 0; sample < params.samples; sample++) {
        if (render_stop) return;
        parallel_for(state.width, state.height, [&](int i, int j) {
          if (render_stop) return;
          trace_sample(state, scene, bvh, lights, i, j, params);
        });
        state.samples += 1;
        if (!render_stop) {
          auto lock      = std::lock_guard{render_mutex};
          render_current = 0;
          if (!params.denoise || render_stop) {
            get_render(render, state);
          } else {
            get_denoised(render, state);
          }
          image = render;
          tonemap_image_mt(display, image, params.exposure, params.filmic);
          render_update = true;
        }
      }
    });
  };

  // stop render
  auto stop_render = [&]() {
    render_stop = true;
    if (render_worker.valid()) render_worker.get();
  };

  // start rendering
  reset_display();

  // prepare selection
  auto selection = scene_selection{};

  // callbacks
  auto callbacks    = gui_callbacks{};
  callbacks.init_cb = [&](gui_window* win, const gui_input& input) {
    auto lock = std::lock_guard{render_mutex};
    init_image(glimage);
    set_image(glimage, display);
  };
  callbacks.clear_cb = [&](gui_window* win, const gui_input& input) {
    clear_image(glimage);
  };
  callbacks.draw_cb = [&](gui_window* win, const gui_input& input) {
    // update image
    if (render_update) {
      auto lock = std::lock_guard{render_mutex};
      set_image(glimage, display);
      render_update = false;
    }
    update_image_params(win, input, image, glparams);
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
      edited += draw_slider(win, "clamp", tparams.clamp, 10, 1000);
      edited += draw_checkbox(win, "envhidden", tparams.envhidden);
      continue_line(win);
      edited += draw_checkbox(win, "filter", tparams.tentfilter);
      edited += draw_slider(win, "pratio", tparams.pratio, 1, 64);
      // edited += draw_slider(win, "exposure", tparams.exposure, -5, 5);
      end_header(win);
      if (edited) {
        stop_render();
        params = tparams;
        reset_display();
      }
    }
    if (begin_header(win, "tonemap")) {
      edited += draw_slider(win, "exposure", params.exposure, -5, 5);
      edited += draw_checkbox(win, "filmic", params.filmic);
      edited += draw_checkbox(win, "denoise", params.denoise);
      end_header(win);
      if (edited) {
        tonemap_image_mt(display, image, params.exposure, params.filmic);
        set_image(glimage, display);
      }
    }
    draw_image_inspector(win, input, image, display, glparams);
    if (edit) {
      if (draw_scene_editor(win, scene, selection, [&]() { stop_render(); })) {
        reset_display();
      }
    }
  };
  callbacks.uiupdate_cb = [&](gui_window* win, const gui_input& input) {
    auto camera = scene.cameras[params.camera];
    if (uiupdate_camera_params(win, input, camera)) {
      stop_render();
      scene.cameras[params.camera] = camera;
      reset_display();
    }
  };

  // run ui
  run_ui({1280 + 320, 720}, title, callbacks);

  // done
  stop_render();
}

void glview_scene(const string& title, const string& name, scene_model& scene,
    const glscene_params& params_, const glview_callback& widgets_callback,
    const glview_callback& uiupdate_callback,
    const glview_callback& update_callback) {
  // glscene
  auto glscene = glscene_state{};

  // draw params
  auto params = params_;

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

  // gpu updates
  auto updated_shapes   = vector<int>{};
  auto updated_textures = vector<int>{};

  // callbacks
  auto callbacks    = gui_callbacks{};
  callbacks.init_cb = [&](gui_window* win, const gui_input& input) {
    init_glscene(glscene, scene);
  };
  callbacks.clear_cb = [&](gui_window* win, const gui_input& input) {
    clear_scene(glscene);
  };
  callbacks.draw_cb = [&](gui_window* win, const gui_input& input) {
    draw_scene(glscene, scene, input.framebuffer_viewport, params);
  };
  callbacks.widgets_cb = [&](gui_window* win, const gui_input& input) {
    draw_combobox(win, "name", selected, names);
    if (begin_header(win, "shade")) {
      draw_combobox(win, "camera", params.camera, camera_names);
      draw_checkbox(win, "wireframe", params.wireframe);
      continue_line(win);
      draw_checkbox(win, "faceted", params.faceted);
      continue_line(win);
      draw_checkbox(win, "double sided", params.double_sided);
      draw_combobox(
          win, "lighting", (int&)params.lighting, glscene_lighting_names);
      draw_slider(win, "exposure", params.exposure, -10, 10);
      draw_slider(win, "gamma", params.gamma, 0.1f, 4);
      draw_slider(win, "near", params.near, 0.01f, 1.0f);
      draw_slider(win, "far", params.far, 1000.0f, 10000.0f);
      end_header(win);
    }
    // draw_scene_editor(win, scene, selection, {});
    if (widgets_callback) {
      widgets_callback(win, input, updated_shapes, updated_textures);
      if (!updated_shapes.empty() || !updated_textures.empty()) {
        update_glscene(glscene, scene, updated_shapes, updated_textures);
        updated_shapes.clear();
        updated_textures.clear();
      }
    }
  };
  callbacks.update_cb = [&](gui_window* win, const gui_input& input) {
    if (update_callback) {
      update_callback(win, input, updated_shapes, updated_textures);
      if (!updated_shapes.empty() || !updated_textures.empty()) {
        update_glscene(glscene, scene, updated_shapes, updated_textures);
        updated_shapes.clear();
        updated_textures.clear();
      }
    }
  };
  callbacks.uiupdate_cb = [&](gui_window* win, const gui_input& input) {
    // handle mouse and keyboard for navigation
    if (uiupdate_callback) {
      uiupdate_callback(win, input, updated_shapes, updated_textures);
      if (!updated_shapes.empty() || !updated_textures.empty()) {
        update_glscene(glscene, scene, updated_shapes, updated_textures);
        updated_shapes.clear();
        updated_textures.clear();
      }
    }
    auto camera = scene.cameras.at(params.camera);
    if (uiupdate_camera_params(win, input, camera)) {
      scene.cameras.at(params.camera) = camera;
    }
  };

  // run ui
  run_ui({1280 + 320, 720}, "yshade", callbacks);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// OPENGL HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

// assert on error
[[maybe_unused]] static GLenum _assert_ogl_error() {
  auto error_code = glGetError();
  if (error_code != GL_NO_ERROR) {
    auto error = ""s;
    switch (error_code) {
      case GL_INVALID_ENUM: error = "INVALID_ENUM"; break;
      case GL_INVALID_VALUE: error = "INVALID_VALUE"; break;
      case GL_INVALID_OPERATION: error = "INVALID_OPERATION"; break;
      // case GL_STACK_OVERFLOW: error = "STACK_OVERFLOW"; break;
      // case GL_STACK_UNDERFLOW: error = "STACK_UNDERFLOW"; break;
      case GL_OUT_OF_MEMORY: error = "OUT_OF_MEMORY"; break;
      case GL_INVALID_FRAMEBUFFER_OPERATION:
        error = "INVALID_FRAMEBUFFER_OPERATION";
        break;
    }
    printf("\n    OPENGL ERROR: %s\n\n", error.c_str());
  }
  return error_code;
}
static void assert_glerror() { assert(_assert_ogl_error() == GL_NO_ERROR); }

// initialize program
void set_program(uint& program_id, uint& vertex_id, uint& fragment_id,
    const string& vertex, const string& fragment) {
  // error
  auto program_error = [&](const char* message, const char* log) {
    if (program_id) glDeleteProgram(program_id);
    if (vertex_id) glDeleteShader(program_id);
    if (fragment_id) glDeleteShader(program_id);
    program_id  = 0;
    vertex_id   = 0;
    fragment_id = 0;
    printf("%s\n", message);
    printf("%s\n", log);
  };

  const char* ccvertex   = vertex.data();
  const char* ccfragment = fragment.data();
  auto        errflags   = 0;
  auto        errbuf     = array<char, 10000>{};

  assert_glerror();

  // create vertex
  vertex_id = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vertex_id, 1, &ccvertex, NULL);
  glCompileShader(vertex_id);
  glGetShaderiv(vertex_id, GL_COMPILE_STATUS, &errflags);
  if (errflags == 0) {
    glGetShaderInfoLog(vertex_id, 10000, 0, errbuf.data());
    return program_error("vertex shader not compiled", errbuf.data());
  }
  assert_glerror();

  // create fragment
  fragment_id = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fragment_id, 1, &ccfragment, NULL);
  glCompileShader(fragment_id);
  glGetShaderiv(fragment_id, GL_COMPILE_STATUS, &errflags);
  if (errflags == 0) {
    glGetShaderInfoLog(fragment_id, 10000, 0, errbuf.data());
    return program_error("fragment shader not compiled", errbuf.data());
  }
  assert_glerror();

  // create program
  program_id = glCreateProgram();
  glAttachShader(program_id, vertex_id);
  glAttachShader(program_id, fragment_id);
  glLinkProgram(program_id);
  glGetProgramiv(program_id, GL_LINK_STATUS, &errflags);
  if (errflags == 0) {
    glGetProgramInfoLog(program_id, 10000, 0, errbuf.data());
    return program_error("program not linked", errbuf.data());
  }
// TODO(fabio): Apparently validation must be done just before drawing.
//    https://community.khronos.org/t/samplers-of-different-types-use-the-same-textur/66329
// If done here, validation fails when using cubemaps and textures in the
// same shader. We should create a function validate_program() anc call it
// separately.
#if 0
  glValidateProgram(program_id);
  glGetProgramiv(program_id, GL_VALIDATE_STATUS, &errflags);
  if (!errflags) {
    glGetProgramInfoLog(program_id, 10000, 0, errbuf.data());
    return program_error("program not validated", errbuf.data());
  }
  assert_glerror();
#endif
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// HIGH-LEVEL OPENGL IMAGE DRAWING
// -----------------------------------------------------------------------------
namespace yocto {

static auto glimage_vertex =
    R"(
#version 330
in vec2 positions;
out vec2 frag_texcoord;
uniform vec2 window_size, image_size;
uniform vec2 image_center;
uniform float image_scale;
void main() {
    vec2 pos = (positions * 0.5) * image_size * image_scale + image_center;
    gl_Position = vec4(2 * pos.x / window_size.x - 1, 1 - 2 * pos.y / window_size.y, 0, 1);
    frag_texcoord = positions * 0.5 + 0.5;
}
)";
#if 0
static auto glimage_vertex = R"(
#version 330
in vec2 positions;
out vec2 frag_texcoord;
uniform vec2 window_size, image_size, border_size;
uniform vec2 image_center;
uniform float image_scale;
void main() {
    vec2 pos = (positions * 0.5) * (image_size + border_size*2) * image_scale + image_center;
    gl_Position = vec4(2 * pos.x / window_size.x - 1, 1 - 2 * pos.y / window_size.y, 0.1, 1);
    frag_texcoord = positions * 0.5 + 0.5;
}
)";
#endif
static auto glimage_fragment =
    R"(
#version 330
in vec2 frag_texcoord;
out vec4 frag_color;
uniform sampler2D txt;
void main() {
    frag_color = texture(txt, frag_texcoord);
}
)";
#if 0
static auto glimage_fragment = R"(
#version 330
in vec2 frag_texcoord;
out vec4 frag_color;
uniform vec2 image_size, border_size;
uniform float image_scale;
void main() {
    ivec2 imcoord = ivec2(frag_texcoord * (image_size + border_size*2) - border_size);
    ivec2 tilecoord = ivec2(frag_texcoord * (image_size + border_size*2) * image_scale - border_size);
    ivec2 tile = tilecoord / 16;
    if(imcoord.x <= 0 || imcoord.y <= 0 || 
        imcoord.x >= image_size.x || imcoord.y >= image_size.y) frag_color = vec4(0,0,0,1);
    else if((tile.x + tile.y) % 2 == 0) frag_color = vec4(0.1,0.1,0.1,1);
    else frag_color = vec4(0.3,0.3,0.3,1);
}
)";
#endif

// init image program
bool init_image(glimage_state& glimage) {
  // program
  set_program(glimage.program, glimage.vertex, glimage.fragment, glimage_vertex,
      glimage_fragment);

  // vertex arrays
  glGenVertexArrays(1, &glimage.vertexarray);
  glBindVertexArray(glimage.vertexarray);

  // buffers
  auto positions = vector<vec3f>{
      {-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0}};
  glGenBuffers(1, &glimage.positions);
  glBindBuffer(GL_ARRAY_BUFFER, glimage.positions);
  glBufferData(GL_ARRAY_BUFFER, sizeof(vec3f) * positions.size(),
      positions.data(), GL_STATIC_DRAW);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, false, 0, nullptr);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  auto triangles = vector<vec3i>{{0, 1, 3}, {3, 2, 1}};
  glGenBuffers(1, &glimage.triangles);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, glimage.triangles);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(vec3i) * triangles.size(),
      triangles.data(), GL_STATIC_DRAW);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

  // done
  // glBindVertexArray(0);
  return true;
}

// clear an opengl image
void clear_image(glimage_state& glimage) {
  if (glimage.texture) glDeleteTextures(1, &glimage.texture);
  if (glimage.program) glDeleteProgram(glimage.program);
  if (glimage.vertex) glDeleteProgram(glimage.vertex);
  if (glimage.fragment) glDeleteProgram(glimage.fragment);
  if (glimage.vertexarray) glDeleteVertexArrays(1, &glimage.vertexarray);
  if (glimage.positions) glDeleteBuffers(1, &glimage.positions);
  if (glimage.triangles) glDeleteBuffers(1, &glimage.triangles);
  glimage = {};
}

void set_image(glimage_state& glimage, const color_image& img) {
  if (!glimage.texture || glimage.width != img.width ||
      glimage.height != img.height) {
    if (!glimage.texture) glGenTextures(1, &glimage.texture);
    glBindTexture(GL_TEXTURE_2D, glimage.texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img.width, img.height, 0, GL_RGBA,
        GL_FLOAT, img.pixels.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  } else {
    glBindTexture(GL_TEXTURE_2D, glimage.texture);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img.width, img.height, GL_RGBA,
        GL_FLOAT, img.pixels.data());
  }
  glimage.width  = img.width;
  glimage.height = img.height;
}

// draw image
void draw_image(glimage_state& glimage, const glimage_params& params) {
  // check errors
  assert_glerror();

  // viewport and framebuffer
  glViewport(params.framebuffer.x, params.framebuffer.y, params.framebuffer.z,
      params.framebuffer.w);
  glClearColor(params.background.x, params.background.y, params.background.z,
      params.background.w);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);

  // bind program and params
  glUseProgram(glimage.program);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, glimage.texture);
  glUniform1i(glGetUniformLocation(glimage.program, "txt"), 0);
  glUniform2f(glGetUniformLocation(glimage.program, "window_size"),
      (float)params.window.x, (float)params.window.y);
  glUniform2f(glGetUniformLocation(glimage.program, "image_size"),
      (float)glimage.width, (float)glimage.height);
  glUniform2f(glGetUniformLocation(glimage.program, "image_center"),
      params.center.x, params.center.y);
  glUniform1f(
      glGetUniformLocation(glimage.program, "image_scale"), params.scale);
  assert_glerror();

  // draw
  glBindVertexArray(glimage.vertexarray);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, glimage.triangles);
  glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, nullptr);
  glBindVertexArray(0);
  assert_glerror();

  // unbind program
  glUseProgram(0);
  assert_glerror();
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE DRAWING
// -----------------------------------------------------------------------------
namespace yocto {

#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverlength-strings"
#endif

static const char* glscene_vertex =
    R"(
#version 330

layout(location = 0) in vec3 positions;           // vertex position (in mesh coordinate frame)
layout(location = 1) in vec3 normals;             // vertex normal (in mesh coordinate frame)
layout(location = 2) in vec2 texcoords;           // vertex texcoords
layout(location = 3) in vec4 colors;              // vertex color
layout(location = 4) in vec4 tangents;            // vertex tangent space

uniform mat4 frame;             // shape transform
uniform mat4 frameit;           // shape transform

uniform mat4 view;              // inverse of the camera frame (as a matrix)
uniform mat4 projection;        // camera projection

out vec3 position;              // [to fragment shader] vertex position (in world coordinate)
out vec3 normal;                // [to fragment shader] vertex normal (in world coordinate)
out vec2 texcoord;              // [to fragment shader] vertex texture coordinates
out vec4 scolor;                // [to fragment shader] vertex color
out vec4 tangsp;                // [to fragment shader] vertex tangent space

// main function
void main() {
  // copy values
  position = positions;
  normal = normals;
  tangsp = tangents;
  texcoord = texcoords;
  scolor = colors;

  // world projection
  position = (frame * vec4(position,1)).xyz;
  normal = (frameit * vec4(normal,0)).xyz;
  tangsp.xyz = (frame * vec4(tangsp.xyz,0)).xyz;

  // clip
  gl_Position = projection * view * vec4(position,1);
}
)";

static const char* glscene_fragment =
    R"(
#version 330

in vec3 position;  // [from vertex shader] position in world space
in vec3 normal;    // [from vertex shader] normal in world space
in vec2 texcoord;  // [from vertex shader] texcoord
in vec4 scolor;    // [from vertex shader] color
in vec4 tangsp;    // [from vertex shader] tangent space

uniform int element;
uniform bool unlit;
uniform bool faceted;
uniform vec4 highlight;
uniform bool double_sided;

uniform vec3 emission;            // material ke
uniform vec3 color;               // material kd
uniform float specular;           // material ks
uniform float metallic;           // material km
uniform float roughness;          // material rs
uniform float opacity;            // material op

uniform bool emission_tex_on;     // material ke texture on
uniform sampler2D emission_tex;   // material ke texture
uniform bool color_tex_on;        // material kd texture on
uniform sampler2D color_tex;      // material kd texture
uniform bool roughness_tex_on;    // material rs texture on
uniform sampler2D roughness_tex;  // material rs texture
uniform bool normalmap_tex_on;    // material normal texture on
uniform sampler2D normalmap_tex;  // material normal texture

uniform int  lighting;            // eyelight shading
uniform vec3 ambient;             // ambient light
uniform int  lights_num;          // number of lights
uniform vec3 lights_direction[16];// light positions
uniform vec3 lights_emission[16]; // light intensities

uniform mat4 frame;              // shape transform
uniform mat4 frameit;            // shape transform

uniform vec3 eye;              // camera position
uniform mat4 view;             // inverse of the camera frame (as a matrix)
uniform mat4 projection;       // camera projection

uniform float exposure; 
uniform float gamma;

out vec4 frag_color;      

float pif = 3.14159265;

struct shade_brdf {
  vec3  emission;
  vec3  diffuse;
  vec3  specular;
  float roughness;
  float opacity;
};

vec3 eval_brdf_color(vec3 value, sampler2D tex, bool tex_on) {
  vec3 result = value;
  if (tex_on) result *= texture(tex, texcoord).rgb;
  return result;
}
float eval_brdf_value(float value, sampler2D tex, bool tex_on) {
  float result = value;
  if (tex_on) result *= texture(tex, texcoord).r;
  return result;
}

shade_brdf eval_brdf() {
  vec4 emission_t = vec4(emission, 1);
  if (emission_tex_on) emission_t *= texture(emission_tex, texcoord);
  vec4 base_t = scolor * vec4(color, opacity);
  if (color_tex_on) base_t *= pow(texture(color_tex, texcoord), vec4(2.2,2.2,2.2,1));
  float metallic_t = metallic;
  float roughness_t = roughness;
  roughness_t = roughness_t * roughness_t;
  if (roughness_t < 0.03 * 0.03) roughness_t = 0.03 * 0.03;
  float specular_t = specular;

  // color?
  shade_brdf brdf;
  brdf.emission  = emission_t.xyz;
  brdf.diffuse   = base_t.xyz * (1 - metallic_t);
  brdf.specular  = specular_t * (base_t.xyz * metallic_t + vec3(0.04) * (1 - metallic_t));
  brdf.roughness = roughness_t;
  brdf.opacity   = base_t.w;
  return brdf;
}

vec3 eval_brdfcos(shade_brdf brdf, vec3 n, vec3 incoming, vec3 outgoing) {
  vec3 halfway = normalize(incoming+outgoing);
  float ndi = dot(incoming,n), ndo = dot(outgoing,n), ndh = dot(halfway,n);
  if(ndi<=0 || ndo <=0) return vec3(0);
  vec3 diff = ndi * brdf.diffuse / pif;
  if(ndh<=0) return diff;
  float cos2 = ndh * ndh;
  float tan2 = (1 - cos2) / cos2;
  float alpha2 = brdf.roughness * brdf.roughness;
  float d = alpha2 / (pif * cos2 * cos2 * (alpha2 + tan2) * (alpha2 + tan2));
  float lambda_o = (-1 + sqrt(1 + (1 - ndo * ndo) / (ndo * ndo))) / 2;
  float lambda_i = (-1 + sqrt(1 + (1 - ndi * ndi) / (ndi * ndi))) / 2;
  float g = 1 / (1 + lambda_o + lambda_i);
  vec3 spec = ndi * brdf.specular * d * g / (4*ndi*ndo);
  return diff+spec;
}

vec3 apply_normal_map(vec2 texcoord, vec3 normal, vec4 tangsp) {
  if(!normalmap_tex_on) return normal;
  vec3 tangu = normalize((frame * vec4(normalize(tangsp.xyz),0)).xyz);
  vec3 tangv = normalize(cross(normal, tangu));
  if(tangsp.w < 0) tangv = -tangv;
  vec3 texture = 2 * texture(normalmap_tex,texcoord).xyz - 1;
  // texture.y = -texture.y;
  return normalize( tangu * texture.x + tangv * texture.y + normal * texture.z );
}

vec3 triangle_normal(vec3 position) {
  vec3 fdx = dFdx(position); 
  vec3 fdy = dFdy(position); 
  return normalize((frame * vec4(normalize(cross(fdx, fdy)), 0)).xyz);
}

#define element_points 1
#define element_lines 2
#define element_triangles 3

vec3 eval_normal(vec3 outgoing) {
  vec3 norm;
  if (element == element_triangles) {
    if (faceted) {
      norm = triangle_normal(position);
    } else {
      norm = normalize(normal);
    }
    // apply normal map
    norm = apply_normal_map(texcoord, norm, tangsp);
    if (double_sided) norm = faceforward(norm, -outgoing, norm);
  }

  if (element == element_lines) {
    vec3 tangent = normalize(normal);
    norm         = normalize(outgoing - tangent * dot(outgoing, tangent));
  }

  return norm;
}

#define lighting_eyelight 0
#define lighting_camlight 1

// main
void main() {
  // view vector
  vec3 outgoing = normalize(eye - position);
  vec3 n = eval_normal(outgoing);

  // get material color from textures
  shade_brdf brdf = eval_brdf();
  if(brdf.opacity < 0.005) discard;

  if(unlit) {
    frag_color = vec4(brdf.emission + brdf.diffuse, brdf.opacity);
    return; 
  }

  // emission
  vec3 radiance = brdf.emission;

  // check early exit
  if(brdf.diffuse != vec3(0,0,0) || brdf.specular != vec3(0,0,0)) {
    // eyelight shading
    if(lighting == lighting_eyelight) {
      vec3 incoming = outgoing;
      radiance += pif * eval_brdfcos(brdf, n, incoming, outgoing);
    }
    if(lighting == lighting_camlight) {
      // accumulate ambient
      radiance += ambient * brdf.diffuse;
      // foreach light
      for(int lid = 0; lid < lights_num; lid ++) {
        radiance += lights_emission[lid] * 
          eval_brdfcos(brdf, n, lights_direction[lid], outgoing);
      }
    }
  }

  // final color correction
  radiance = pow(radiance * pow(2,exposure), vec3(1/gamma));

  // highlighting
  if(highlight.w > 0) {
    if(mod(int(gl_FragCoord.x)/4 + int(gl_FragCoord.y)/4, 2)  == 0)
        radiance = highlight.xyz * highlight.w + radiance * (1-highlight.w);
  }

  // output final color by setting gl_FragColor
  frag_color = vec4(radiance, brdf.opacity);
}
)";

#ifndef _WIN32
#pragma GCC diagnostic pop
#endif

// Create texture
void set_texture(glscene_texture& gltexture, const scene_texture& texture) {
  if (!gltexture.texture || gltexture.width != texture.width ||
      gltexture.height != texture.height) {
    if (!gltexture.texture) glGenTextures(1, &gltexture.texture);
    glBindTexture(GL_TEXTURE_2D, gltexture.texture);
    if (!texture.pixelsb.empty()) {
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texture.width, texture.height, 0,
          GL_RGBA, GL_UNSIGNED_BYTE, texture.pixelsb.data());
    } else {
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texture.width, texture.height, 0,
          GL_RGBA, GL_FLOAT, texture.pixelsf.data());
    }
    glGenerateMipmap(GL_TEXTURE_2D);
    glTexParameteri(
        GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  } else {
    glBindTexture(GL_TEXTURE_2D, gltexture.texture);
    if (!texture.pixelsb.empty()) {
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, texture.width, texture.height,
          GL_RGBA, GL_UNSIGNED_BYTE, texture.pixelsb.data());
    } else {
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, texture.width, texture.height,
          GL_RGBA, GL_FLOAT, texture.pixelsf.data());
    }
    glGenerateMipmap(GL_TEXTURE_2D);
  }
}

// Clean texture
void clear_texture(glscene_texture& gltexture) {
  if (gltexture.texture) {
    glDeleteTextures(1, &gltexture.texture);
    gltexture.texture = 0;
  }
}

// Create shape
void set_shape(glscene_shape& glshape, const scene_shape& shape) {
  auto set_vertex = [](uint& buffer, int& num, const auto& data,
                        const auto& def, int location) {
    if (data.empty()) {
      if (buffer) glDeleteBuffers(1, &buffer);
      buffer = 0;
      num    = 0;
      glDisableVertexAttribArray(location);
      if constexpr (sizeof(def) == sizeof(float))
        glVertexAttrib1f(location, (float)def);
      if constexpr (sizeof(def) == sizeof(vec2f))
        glVertexAttrib2fv(location, (float*)&def.x);
      if constexpr (sizeof(def) == sizeof(vec3f))
        glVertexAttrib3fv(location, (float*)&def.x);
      if constexpr (sizeof(def) == sizeof(vec4f))
        glVertexAttrib4fv(location, (float*)&def.x);
    } else {
      if (!buffer || (int)data.size() != num) {
        if (buffer) glDeleteBuffers(1, &buffer);
        glGenBuffers(1, &buffer);
        glBindBuffer(GL_ARRAY_BUFFER, buffer);
        glBufferData(GL_ARRAY_BUFFER, data.size() * sizeof(data.front()),
            data.data(), GL_STATIC_DRAW);
        num = (int)data.size();
      } else {
        // we have enough space
        glBindBuffer(GL_ARRAY_BUFFER, buffer);
        glBufferSubData(GL_ARRAY_BUFFER, 0, data.size() * sizeof(data.front()),
            data.data());
      }
      glBindBuffer(GL_ARRAY_BUFFER, buffer);
      glEnableVertexAttribArray(location);
      glVertexAttribPointer(location, sizeof(data.front()) / sizeof(float),
          GL_FLOAT, false, 0, nullptr);
    }
  };

  auto set_indices = [](uint& buffer, int& num, const auto& data) {
    if (data.empty()) {
      if (buffer) glDeleteBuffers(1, &buffer);
      buffer = 0;
      num    = 0;
    } else {
      if (!buffer || (int)data.size() != num) {
        if (buffer) glDeleteBuffers(1, &buffer);
        glGenBuffers(1, &buffer);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,
            data.size() * sizeof(data.front()), data.data(), GL_STATIC_DRAW);
        num = (int)data.size();
      } else {
        // we have enough space
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer);
        glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0,
            data.size() * sizeof(data.front()), data.data());
      }
    }
  };

  if (!glshape.vertexarray) glGenVertexArrays(1, &glshape.vertexarray);
  glBindVertexArray(glshape.vertexarray);
  set_indices(glshape.points, glshape.num_points, shape.points);
  set_indices(glshape.lines, glshape.num_lines, shape.lines);
  set_indices(glshape.triangles, glshape.num_triangles, shape.triangles);
  set_indices(
      glshape.quads, glshape.num_quads, quads_to_triangles(shape.quads));
  set_vertex(glshape.positions, glshape.num_positions, shape.positions,
      vec3f{0, 0, 0}, 0);
  set_vertex(
      glshape.normals, glshape.num_normals, shape.normals, vec3f{0, 0, 1}, 1);
  set_vertex(glshape.texcoords, glshape.num_texcoords, shape.texcoords,
      vec2f{0, 0}, 2);
  set_vertex(
      glshape.colors, glshape.num_colors, shape.colors, vec4f{1, 1, 1, 1}, 3);
  set_vertex(glshape.tangents, glshape.num_tangents, shape.tangents,
      vec4f{0, 0, 1, 1}, 4);
  glBindVertexArray(0);
}

// Clean shape
void clear_shape(glscene_shape& glshape) {
  if (glshape.vertexarray) glDeleteVertexArrays(1, &glshape.vertexarray);
  if (glshape.positions) glDeleteBuffers(1, &glshape.positions);
  if (glshape.normals) glDeleteBuffers(1, &glshape.normals);
  if (glshape.texcoords) glDeleteBuffers(1, &glshape.texcoords);
  if (glshape.colors) glDeleteBuffers(1, &glshape.colors);
  if (glshape.tangents) glDeleteBuffers(1, &glshape.tangents);
  if (glshape.points) glDeleteBuffers(1, &glshape.points);
  if (glshape.lines) glDeleteBuffers(1, &glshape.lines);
  if (glshape.triangles) glDeleteBuffers(1, &glshape.triangles);
  if (glshape.quads) glDeleteBuffers(1, &glshape.quads);
  glshape = {};
  assert_glerror();
}

// init scene
void init_glscene(glscene_state& glscene, const scene_model& ioscene) {
  // program
  set_program(glscene.program, glscene.vertex, glscene.fragment, glscene_vertex,
      glscene_fragment);

  // textures
  for (auto& iotexture : ioscene.textures) {
    auto& gltexture = glscene.textures.emplace_back();
    set_texture(gltexture, iotexture);
  }

  // shapes
  for (auto& ioshape : ioscene.shapes) {
    auto& glshape = glscene.shapes.emplace_back();
    set_shape(glshape, ioshape);
  }
}

// update scene
void update_glscene(glscene_state& glscene, const scene_model& scene,
    const vector<int>& updated_shapes, const vector<int>& updated_textures) {
  for (auto shape_id : updated_shapes) {
    set_shape(glscene.shapes[shape_id], scene.shapes[shape_id]);
  }
  for (auto texture_id : updated_textures) {
    set_texture(glscene.textures[texture_id], scene.textures[texture_id]);
  }
}

// Clear an OpenGL scene
void clear_scene(glscene_state& glscene) {
  for (auto& texture : glscene.textures) clear_texture(texture);
  for (auto& shape : glscene.shapes) clear_shape(shape);
  if (glscene.program) glDeleteProgram(glscene.program);
  if (glscene.vertex) glDeleteProgram(glscene.vertex);
  if (glscene.fragment) glDeleteProgram(glscene.fragment);
}

[[maybe_unused]] static void draw_shape(glscene_shape& shape) {
  if (shape.vertexarray == 0) return;
  glBindVertexArray(shape.vertexarray);

  if (shape.points) {
    glPointSize(shape.point_size);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, shape.points);
    glDrawElements(
        GL_POINTS, (GLsizei)shape.num_points * 1, GL_UNSIGNED_INT, nullptr);
    glPointSize(shape.point_size);
  }
  if (shape.lines) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, shape.lines);
    glDrawElements(
        GL_LINES, (GLsizei)shape.num_lines * 2, GL_UNSIGNED_INT, nullptr);
  }
  if (shape.triangles) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, shape.triangles);
    glDrawElements(GL_TRIANGLES, (GLsizei)shape.num_triangles * 3,
        GL_UNSIGNED_INT, nullptr);
  }
  if (shape.quads) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, shape.quads);
    glDrawElements(
        GL_TRIANGLES, (GLsizei)shape.num_quads * 3, GL_UNSIGNED_INT, nullptr);
  }

  glBindVertexArray(0);
  assert_glerror();
}

void draw_scene(glscene_state& glscene, const scene_model& scene,
    const vec4i& viewport, const glscene_params& params) {
  // check errors
  assert_glerror();

  // viewport and framebuffer
  glViewport(viewport.x, viewport.y, viewport.z, viewport.w);
  glClearColor(params.background.x, params.background.y, params.background.z,
      params.background.w);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);

  // set program
  auto& program = glscene.program;
  glUseProgram(program);

  // camera
  auto& camera        = scene.cameras.at(params.camera);
  auto  camera_aspect = (float)viewport.z / (float)viewport.w;
  auto  camera_yfov =
      camera_aspect >= 0
           ? (2 * atan(camera.film / (camera_aspect * 2 * camera.lens)))
           : (2 * atan(camera.film / (2 * camera.lens)));
  auto view_matrix       = frame_to_mat(inverse(camera.frame));
  auto projection_matrix = perspective_mat(
      camera_yfov, camera_aspect, params.near, params.far);
  glUniform3f(glGetUniformLocation(program, "eye"), camera.frame.o.x,
      camera.frame.o.y, camera.frame.o.z);
  glUniformMatrix4fv(
      glGetUniformLocation(program, "view"), 1, false, &view_matrix.x.x);
  glUniformMatrix4fv(glGetUniformLocation(program, "projection"), 1, false,
      &projection_matrix.x.x);

  // params
  glUniform1f(glGetUniformLocation(program, "exposure"), params.exposure);
  glUniform1f(glGetUniformLocation(program, "gamma"), params.gamma);
  glUniform1i(glGetUniformLocation(program, "double_sided"),
      params.double_sided ? 1 : 0);

  static auto lights_direction = vector<vec3f>{normalize(vec3f{1, 1, 1}),
      normalize(vec3f{-1, 1, 1}), normalize(vec3f{-1, -1, 1}),
      normalize(vec3f{0.1, 0.5, -1})};
  static auto lights_emission  = vector<vec3f>{vec3f{pif / 2, pif / 2, pif / 2},
      vec3f{pif / 2, pif / 2, pif / 2}, vec3f{pif / 4, pif / 4, pif / 4},
      vec3f{pif / 4, pif / 4, pif / 4}};
  if (params.lighting == glscene_lighting_type::camlight) {
    glUniform1i(glGetUniformLocation(program, "lighting"), 1);
    glUniform3f(glGetUniformLocation(program, "ambient"), 0, 0, 0);
    glUniform1i(glGetUniformLocation(program, "lights_num"),
        (int)lights_direction.size());
    for (auto lid = 0; lid < lights_direction.size(); lid++) {
      auto is        = std::to_string(lid);
      auto direction = transform_direction(camera.frame, lights_direction[lid]);
      glUniform3f(glGetUniformLocation(
                      program, ("lights_direction[" + is + "]").c_str()),
          direction.x, direction.y, direction.z);
      glUniform3f(glGetUniformLocation(
                      program, ("lights_emission[" + is + "]").c_str()),
          lights_emission[lid].x, lights_emission[lid].y,
          lights_emission[lid].z);
    }
  } else if (params.lighting == glscene_lighting_type::eyelight) {
    glUniform1i(glGetUniformLocation(program, "lighting"), 0);
    glUniform1i(glGetUniformLocation(program, "lights_num"), 0);
  } else {
    throw std::invalid_argument{"unknown lighting type"};
  }

  // helper
  auto set_texture = [&glscene](uint program, const char* name,
                         const char* name_on, int texture_idx, int unit) {
    if (texture_idx >= 0) {
      auto& gltexture = glscene.textures.at(texture_idx);
      glActiveTexture(GL_TEXTURE0 + unit);
      glBindTexture(GL_TEXTURE_2D, gltexture.texture);
      glUniform1i(glGetUniformLocation(program, name), unit);
      glUniform1i(glGetUniformLocation(program, name_on), 1);
    } else {
      glActiveTexture(GL_TEXTURE0 + unit);
      glBindTexture(GL_TEXTURE_2D, 0);
      glUniform1i(glGetUniformLocation(program, name), unit);
      glUniform1i(glGetUniformLocation(program, name_on), 0);
    }
  };

  // draw instances
  if (params.wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  for (auto& instance : scene.instances) {
    auto& glshape  = glscene.shapes.at(instance.shape);
    auto& material = scene.materials.at(instance.material);

    auto shape_xform     = frame_to_mat(instance.frame);
    auto shape_inv_xform = transpose(
        frame_to_mat(inverse(instance.frame, params.non_rigid_frames)));
    glUniformMatrix4fv(
        glGetUniformLocation(program, "frame"), 1, false, &shape_xform.x.x);
    glUniformMatrix4fv(glGetUniformLocation(program, "frameit"), 1, false,
        &shape_inv_xform.x.x);
    glUniform1i(glGetUniformLocation(program, "faceted"),
        (params.faceted || glshape.normals == 0) ? 1 : 0);

    glUniform1i(glGetUniformLocation(program, "unlit"), 0);
    glUniform3f(glGetUniformLocation(program, "emission"), material.emission.x,
        material.emission.y, material.emission.z);
    glUniform3f(glGetUniformLocation(program, "color"), material.color.x,
        material.color.y, material.color.z);
    glUniform1f(glGetUniformLocation(program, "specular"), 1);
    glUniform1f(glGetUniformLocation(program, "metallic"), material.metallic);
    glUniform1f(glGetUniformLocation(program, "roughness"), material.roughness);
    glUniform1f(glGetUniformLocation(program, "opacity"), material.opacity);
    if (material.type == scene_material_type::matte ||
        material.type == scene_material_type::transparent ||
        material.type == scene_material_type::refractive ||
        material.type == scene_material_type::subsurface ||
        material.type == scene_material_type::volume) {
      glUniform1f(glGetUniformLocation(program, "specular"), 0);
    }
    if (material.type == scene_material_type::metallic) {
      glUniform1f(glGetUniformLocation(program, "metallic"), 1);
    }
    glUniform1f(glGetUniformLocation(program, "double_sided"),
        params.double_sided ? 1 : 0);
    set_texture(
        program, "emission_tex", "emission_tex_on", material.emission_tex, 0);
    set_texture(program, "color_tex", "color_tex_on", material.color_tex, 1);
    set_texture(program, "roughness_tex", "roughness_tex_on",
        material.roughness_tex, 3);
    set_texture(
        program, "normalmap_tex", "normalmap_tex_on", material.normal_tex, 5);
    assert_glerror();

    if (glshape.points)
      glUniform1i(glGetUniformLocation(program, "element"), 1);
    if (glshape.lines) glUniform1i(glGetUniformLocation(program, "element"), 2);
    if (glshape.triangles)
      glUniform1i(glGetUniformLocation(program, "element"), 3);
    if (glshape.quads) glUniform1i(glGetUniformLocation(program, "element"), 3);
    assert_glerror();

    glBindVertexArray(glshape.vertexarray);
    if (glshape.points) {
      glPointSize(glshape.point_size);
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, glshape.points);
      glDrawElements(
          GL_POINTS, (GLsizei)glshape.num_points * 1, GL_UNSIGNED_INT, nullptr);
      glPointSize(glshape.point_size);
    }
    if (glshape.lines) {
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, glshape.lines);
      glDrawElements(
          GL_LINES, (GLsizei)glshape.num_lines * 2, GL_UNSIGNED_INT, nullptr);
    }
    if (glshape.triangles) {
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, glshape.triangles);
      glDrawElements(GL_TRIANGLES, (GLsizei)glshape.num_triangles * 3,
          GL_UNSIGNED_INT, nullptr);
    }
    if (glshape.quads) {
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, glshape.quads);
      glDrawElements(GL_TRIANGLES, (GLsizei)glshape.num_quads * 3,
          GL_UNSIGNED_INT, nullptr);
    }

    glBindVertexArray(0);
    assert_glerror();
  }
  if (params.wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  // done
  glUseProgram(0);
}

}  // namespace yocto
