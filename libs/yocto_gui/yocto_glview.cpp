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

void update_image_params(gui_window* win, const gui_input& input,
    const color_image& image, ogl_image_params& glparams) {
  glparams.window                           = input.window_size;
  glparams.framebuffer                      = input.framebuffer_viewport;
  std::tie(glparams.center, glparams.scale) = camera_imview(glparams.center,
      glparams.scale, {image.width, image.height}, glparams.window,
      glparams.fit);
}

void update_image_params(gui_window* win, const gui_input& input,
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
    glparams.scale *= powf(
        2, (input.mouse_pos.x - input.mouse_last.x) * 0.001f);
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

static void init_glscene(shade_scene& glscene, const scene_model& ioscene) {
  // init scene
  init_scene(glscene);

  // camera
  for (auto& iocamera : ioscene.cameras) {
    auto& camera = glscene.cameras.at(add_camera(glscene));
    set_frame(camera, iocamera.frame);
    set_lens(camera, iocamera.lens, iocamera.aspect, iocamera.film);
    set_nearfar(camera, 0.001, 10000);
  }

  // textures
  for (auto& iotexture : ioscene.textures) {
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
    auto  handle     = add_material(glscene);
    auto& glmaterial = glscene.materials[handle];
    set_emission(glmaterial, iomaterial.emission, iomaterial.emission_tex);
    set_opacity(glmaterial, iomaterial.opacity, invalidid);
    set_normalmap(glmaterial, iomaterial.normal_tex);
    switch (iomaterial.type) {
      case scene_material_type::matte: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 0, invalidid);
        set_metallic(glmaterial, 0, invalidid);
        set_roughness(glmaterial, 0, invalidid);
      } break;
      case scene_material_type::glossy: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 1, invalidid);
        set_metallic(glmaterial, 0, invalidid);
        set_roughness(glmaterial, iomaterial.roughness, invalidid);
      } break;
      case scene_material_type::metallic: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 0, invalidid);
        set_metallic(glmaterial, 1, invalidid);
        set_roughness(
            glmaterial, iomaterial.roughness, iomaterial.roughness_tex);
      } break;
      case scene_material_type::gltfpbr: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 1, invalidid);
        set_metallic(glmaterial, iomaterial.metallic, invalidid);
        set_roughness(glmaterial, iomaterial.roughness, invalidid);
      } break;
      default: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 0, invalidid);
        set_metallic(glmaterial, 0, invalidid);
        set_roughness(
            glmaterial, iomaterial.roughness, iomaterial.roughness_tex);
      } break;
    }
  }

  // shapes
  for (auto& ioshape : ioscene.shapes) {
    add_shape(glscene, ioshape.points, ioshape.lines, ioshape.triangles,
        ioshape.quads, ioshape.positions, ioshape.normals, ioshape.texcoords,
        ioshape.colors);
  }

  // shapes
  for (auto& ioinstance : ioscene.instances) {
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
}

static void update_glscene(shade_scene& glscene, const scene_model& scene,
    const vector<int>& updated_shapes, const vector<int>& updated_textures) {
  for (auto shape_id : updated_shapes) {
    auto& shape   = scene.shapes.at(shape_id);
    auto& glshape = glscene.shapes.at(shape_id);
    if (!shape.points.empty()) set_points(glshape, shape.points);
    if (!shape.lines.empty()) set_lines(glshape, shape.lines);
    if (!shape.triangles.empty()) set_triangles(glshape, shape.triangles);
    if (!shape.quads.empty()) set_quads(glshape, shape.quads);
    if (!shape.positions.empty()) set_positions(glshape, shape.positions);
    if (!shape.normals.empty()) set_normals(glshape, shape.normals);
    if (!shape.texcoords.empty()) set_texcoords(glshape, shape.texcoords);
    if (!shape.colors.empty()) set_colors(glshape, shape.colors);
    if (!shape.tangents.empty()) set_tangents(glshape, shape.tangents);
  }
  for (auto texture_id : updated_textures) {
    auto& texture   = scene.textures.at(texture_id);
    auto& gltexture = glscene.textures.at(texture_id);
    if (!texture.pixelsf.empty()) {
      set_texture(gltexture, texture.width, texture.height, texture.pixelsf);
    } else if (!texture.pixelsb.empty()) {
      set_texture(gltexture, texture.width, texture.height, texture.pixelsb);
    }
  }
}

void glview_scene(const string& title, const string& name, scene_model& scene,
    const shade_params& params_, const glview_callback& widgets_callback,
    const glview_callback& uiupdate_callback,
    const glview_callback& update_callback) {
  // glscene
  auto glscene = shade_scene{};

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
    draw_scene(glscene, input.framebuffer_viewport, params);
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
          win, "lighting", (int&)params.lighting, shade_lighting_names);
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
    auto camera = scene.cameras.at(0);
    if (uiupdate_camera_params(win, input, camera)) {
      scene.cameras.at(0) = camera;
      set_frame(glscene.cameras.at(0), camera.frame);
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
static void assert_ogl_error_() { assert(_assert_ogl_error() == GL_NO_ERROR); }

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

  assert_ogl_error_();

  // create vertex
  vertex_id = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vertex_id, 1, &ccvertex, NULL);
  glCompileShader(vertex_id);
  glGetShaderiv(vertex_id, GL_COMPILE_STATUS, &errflags);
  if (errflags == 0) {
    glGetShaderInfoLog(vertex_id, 10000, 0, errbuf.data());
    return program_error("vertex shader not compiled", errbuf.data());
  }
  assert_ogl_error_();

  // create fragment
  fragment_id = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fragment_id, 1, &ccfragment, NULL);
  glCompileShader(fragment_id);
  glGetShaderiv(fragment_id, GL_COMPILE_STATUS, &errflags);
  if (errflags == 0) {
    glGetShaderInfoLog(fragment_id, 10000, 0, errbuf.data());
    return program_error("fragment shader not compiled", errbuf.data());
  }
  assert_ogl_error_();

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
  assert_ogl_error_();
#endif
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// HIGH-LEVEL OPENGL IMAGE DRAWING
// -----------------------------------------------------------------------------
namespace yocto {

static auto ogl_image_vertex =
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
static auto ogl_image_vertex = R"(
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
static auto ogl_image_fragment =
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
static auto ogl_image_fragment = R"(
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
  set_program(glimage.program, glimage.vertex, glimage.fragment,
      ogl_image_vertex, ogl_image_fragment);

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
  assert_ogl_error_();

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
  assert_ogl_error_();

  // draw
  glBindVertexArray(glimage.vertexarray);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, glimage.triangles);
  glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, nullptr);
  glBindVertexArray(0);
  assert_ogl_error_();

  // unbind program
  glUseProgram(0);
  assert_ogl_error_();
}

}  // namespace yocto
