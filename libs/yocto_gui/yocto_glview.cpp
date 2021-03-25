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

static void init_glscene(glscene_state& glscene, const scene_model& ioscene) {
  // init scene
  init_scene(glscene);

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

static void update_glscene(glscene_state& glscene, const scene_model& scene,
    const vector<int>& updated_shapes, const vector<int>& updated_textures) {
  for (auto shape_id : updated_shapes) {
    set_shape(glscene.shapes[shape_id], scene.shapes[shape_id]);
  }
  for (auto texture_id : updated_textures) {
    set_texture(glscene.textures[texture_id], scene.textures[texture_id]);
  }
}

void glview_scene(const string& title, const string& name, scene_model& scene,
    const shade_params& params_, const glview_callback& widgets_callback,
    const glview_callback& uiupdate_callback,
    const glview_callback& update_callback) {
  // glscene
  auto glscene = glscene_state{};

  // draw params
  auto params = (glscene_params&)params_;

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

// -----------------------------------------------------------------------------
// SCENE DRAWING
// -----------------------------------------------------------------------------
namespace yocto {

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
  glshape.vertexarray = 0;
  if (glshape.positions) glDeleteBuffers(1, &glshape.positions);
  if (glshape.normals) glDeleteBuffers(1, &glshape.normals);
  if (glshape.texcoords) glDeleteBuffers(1, &glshape.texcoords);
  if (glshape.colors) glDeleteBuffers(1, &glshape.colors);
  if (glshape.tangents) glDeleteBuffers(1, &glshape.tangents);
  if (glshape.points) glDeleteBuffers(1, &glshape.points);
  if (glshape.lines) glDeleteBuffers(1, &glshape.lines);
  if (glshape.triangles) glDeleteBuffers(1, &glshape.triangles);
  if (glshape.quads) glDeleteBuffers(1, &glshape.quads);
  assert_ogl_error_();
}

glscene_state::~glscene_state() { clear_scene(*this); }

static const char* shade_instance_vertex();
static const char* shade_instanced_vertex();
static const char* shade_instance_fragment();

static const char* shade_environment_fragment();

static const char* precompute_brdflut_vertex();
static const char* precompute_brdflut_fragment();

static const char* precompute_cubemap_vertex();
static const char* precompute_environment_fragment();
static const char* precompute_irradiance_fragment();
static const char* precompute_reflections_fragment();

static void init_environment(
    glscene_state& scene, glscene_environment& environment);
static void init_envlight(
    glscene_state& scene, glscene_environment& environment);

// Initialize an OpenGL scene
void init_scene(glscene_state& scene) {
  if (is_initialized(scene.instance_program)) return;
  set_program(scene.instance_program, shade_instance_vertex(),
      shade_instance_fragment(), true);
  // set_program(scene.envlight_program, shade_instance_vertex(),
  //     shade_instance_fragment(), true);
  set_program(scene.environment_program, precompute_cubemap_vertex(),
      shade_environment_fragment(), true);
}

// Initialize data for environment lighting
void init_environments(glscene_state& scene, bool precompute_envlight) {
  for (auto& environment : scene.environments) {
    init_environment(scene, environment);
    if (precompute_envlight) init_envlight(scene, environment);
  }
}

// Check if we have an envlight
bool has_envlight(const glscene_state& scene) {
  return !scene.environments.empty() && !scene.envlight_shapes.empty();
}

// Clear an OpenGL scene
void clear_scene(glscene_state& scene) {
  for (auto& texture : scene.textures) clear_texture(texture);
  for (auto& shape : scene.shapes) clear_shape(shape);
  for (auto& shape : scene.envlight_shapes) clear_shape(shape);
  for (auto& cubemap : scene.envlight_cubemaps) clear_cubemap(cubemap);
  for (auto& cubemap : scene.envlight_diffuses) clear_cubemap(cubemap);
  for (auto& cubemap : scene.envlight_speculars) clear_cubemap(cubemap);
  for (auto& texture : scene.envlight_brdfluts) clear_texture(texture);
  clear_program(scene.environment_program);
  clear_program(scene.instance_program);
}

// environment properties
glenvironment_handle add_environment(glscene_state& scene) {
  scene.environments.emplace_back();
  return (int)scene.environments.size() - 1;
}
void set_frame(glscene_environment& environment, const frame3f& frame) {
  environment.frame = frame;
}
void set_emission(glscene_environment& environment, const vec3f& emission,
    gltexture_handle emission_tex) {
  environment.emission     = emission;
  environment.emission_tex = emission_tex;
}

glenvironment_handle add_environment(glscene_state& scene, const frame3f& frame,
    const vec3f& emission, gltexture_handle emission_tex) {
  auto  handle      = add_environment(scene);
  auto& environment = scene.environments[handle];
  set_frame(environment, frame);
  set_emission(environment, emission, emission_tex);
  return handle;
}

struct shade_view {
  frame3f camera_frame      = {};
  mat4f   view_matrix       = {};
  mat4f   projection_matrix = {};
};

void set_view_uniforms(ogl_program& program, const shade_view& view) {
  set_uniform(program, "eye", view.camera_frame.o);
  set_uniform(program, "view", view.view_matrix);
  set_uniform(program, "projection", view.projection_matrix);
}

void set_params_uniforms(ogl_program& program, const glscene_params& params) {
  set_uniform(program, "exposure", params.exposure);
  set_uniform(program, "gamma", params.gamma);
  set_uniform(program, "double_sided", params.double_sided);
}

// Draw a shape
void set_instance_uniforms(const glscene_state& scene, ogl_program& program,
    const frame3f& frame, const glscene_shape& shape,
    const scene_material& material, const glscene_params& params) {
  auto shape_xform     = frame_to_mat(frame);
  auto shape_inv_xform = transpose(
      frame_to_mat(inverse(frame, params.non_rigid_frames)));
  set_uniform(program, "frame", shape_xform);
  set_uniform(program, "frameit", shape_inv_xform);
  set_uniform(program, "offset", 0.0f);
  set_uniform(program, "faceted", params.faceted || shape.normals == 0);

  auto set_texture = [&scene](ogl_program& program, const char* name,
                         const char* name_on, int texture_idx, int unit) {
    if (texture_idx >= 0) {
      auto& texture = scene.textures.at(texture_idx);
      glActiveTexture(GL_TEXTURE0 + unit);
      glBindTexture(GL_TEXTURE_2D, texture.texture);
      glUniform1i(glGetUniformLocation(program.program_id, name), unit);
      glUniform1i(glGetUniformLocation(program.program_id, name_on), 1);
    } else {
      glActiveTexture(GL_TEXTURE0 + unit);
      glBindTexture(GL_TEXTURE_2D, 0);
      glUniform1i(glGetUniformLocation(program.program_id, name), unit);
      glUniform1i(glGetUniformLocation(program.program_id, name_on), 0);
    }
  };

  set_uniform(program, "unlit", false);
  set_uniform(program, "emission", material.emission);
  set_uniform(program, "diffuse", material.color);
  set_uniform(program, "specular",
      vec3f{material.metallic, material.metallic, material.metallic});
  set_uniform(program, "roughness", material.roughness);
  set_uniform(program, "opacity", material.opacity);
  set_uniform(program, "double_sided", params.double_sided);
  set_texture(
      program, "emission_tex", "emission_tex_on", material.emission_tex, 0);
  set_texture(program, "diffuse_tex", "diffuse_tex_on", material.color_tex, 1);
  set_texture(program, "specular_tex", "specular_tex_on", -1, 2);
  set_texture(
      program, "roughness_tex", "roughness_tex_on", material.roughness_tex, 3);
  set_texture(program, "opacity_tex", "opacity_tex_on", -1, 4);
  set_texture(
      program, "normalmap_tex", "normalmap_tex_on", material.normal_tex, 5);

  assert_ogl_error();

  if (shape.points) set_uniform(program, "element", 1);
  if (shape.lines) set_uniform(program, "element", 2);
  if (shape.triangles) set_uniform(program, "element", 3);
  if (shape.quads) set_uniform(program, "element", 3);
  assert_ogl_error_();
}

static void draw_shape(glscene_shape& shape) {
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
  assert_ogl_error_();
}

void draw_environments(glscene_state& scene, const shade_view& view,
    const glscene_params& params) {
  if (params.hide_environment) return;
  auto& program = scene.environment_program;
  if (!is_initialized(program)) return;
  bind_program(program);
  set_view_uniforms(program, view);
  set_params_uniforms(program, params);
  for (auto& environment : scene.environments) {
    if (environment.envlight_cubemap == glinvalid_handle) continue;
    set_uniform(program, "emission", environment.emission);
    set_uniform(program, "emission_tex",
        scene.envlight_cubemaps[environment.envlight_cubemap], 0);
    draw_shape(scene.envlight_shapes[environment.envlight_shape]);
  }
  unbind_program();
}

void set_lighting_uniforms(ogl_program& program, const glscene_state& scene,
    const shade_view& view, const glscene_params& params) {
  struct gui_light {
    vec3f position = {0, 0, 0};
    vec3f emission = {0, 0, 0};
    bool  camera   = false;
  };

  static auto camera_light0 = gui_light{
      normalize(vec3f{1, 1, 1}), vec3f{pif / 2, pif / 2, pif / 2}, true};
  static auto camera_light1 = gui_light{
      normalize(vec3f{-1, 1, 1}), vec3f{pif / 2, pif / 2, pif / 2}, true};
  static auto camera_light2 = gui_light{
      normalize(vec3f{-1, -1, 1}), vec3f{pif / 4, pif / 4, pif / 4}, true};
  static auto camera_light3 = gui_light{
      normalize(vec3f{0.1, 0.5, -1}), vec3f{pif / 4, pif / 4, pif / 4}, true};
  static auto camera_lights = vector<gui_light*>{
      &camera_light0, &camera_light1, &camera_light2, &camera_light3};

  auto lighting = params.lighting;
  if (lighting == glscene_lighting_type::envlight && !has_envlight(scene))
    lighting = glscene_lighting_type::camlight;
  if (lighting == glscene_lighting_type::envlight && has_envlight(scene)) {
    if (!has_envlight(scene)) return;
    auto environment = scene.environments.front();
    set_uniform(program, "lighting", 2);
    set_uniform(program, "lights_num", 0);
    set_uniform(program, "envlight_scale", environment.emission);
    set_uniform(program, "envlight_irradiance",
        scene.envlight_diffuses[environment.envlight_diffuse_], 6);
    set_uniform(program, "envlight_reflection",
        scene.envlight_speculars[environment.envlight_specular_], 7);
    set_uniform(program, "envlight_brdflut",
        scene.envlight_brdfluts[environment.envlight_brdflut_], 8);
  } else if (lighting == glscene_lighting_type::camlight) {
    auto& lights = camera_lights;
    set_uniform(program, "lighting", 1);
    set_uniform(program, "ambient", vec3f{0, 0, 0});
    set_uniform(program, "lights_num", (int)lights.size());
    auto lid = 0;
    for (auto light : lights) {
      auto is = std::to_string(lid);
      if (light->camera) {
        auto position = transform_direction(view.camera_frame, light->position);
        set_uniform(program, ("lights_position[" + is + "]").c_str(), position);
      } else {
        set_uniform(
            program, ("lights_position[" + is + "]").c_str(), light->position);
      }
      set_uniform(
          program, ("lights_emission[" + is + "]").c_str(), light->emission);
      set_uniform(program, ("lights_type[" + is + "]").c_str(), 1);
      lid++;
    }
    set_uniform(program, "envlight_irradiance", ogl_cubemap{}, 6);
    set_uniform(program, "envlight_reflection", ogl_cubemap{}, 7);
    set_uniform(program, "envlight_brdflut", ogl_texture{}, 8);
  } else if (lighting == glscene_lighting_type::eyelight) {
    set_uniform(program, "lighting", 0);
    set_uniform(program, "lights_num", 0);
    set_uniform(program, "envlight_irradiance", ogl_cubemap{}, 6);
    set_uniform(program, "envlight_reflection", ogl_cubemap{}, 7);
    set_uniform(program, "envlight_brdflut", ogl_texture{}, 8);
  } else {
    throw std::invalid_argument{"unknown lighting type"};
  }
  assert_ogl_error();
}

void draw_instances(glscene_state& glscene, const scene_model& scene,
    const shade_view& view, const glscene_params& params) {
  // set program
  auto& program = glscene.instance_program;
  bind_program(program);

  // set scene uniforms
  set_view_uniforms(program, view);
  set_params_uniforms(program, params);

  // set lighting uniforms
  set_lighting_uniforms(program, glscene, view, params);

  set_ogl_wireframe(params.wireframe);
  for (auto& instance : scene.instances) {
    // if (instance.hidden) continue;
    set_instance_uniforms(glscene, program, instance.frame,
        glscene.shapes.at(instance.shape),
        scene.materials.at(instance.material), params);
    draw_shape(glscene.shapes.at(instance.shape));
  }
  unbind_program();
}

static shade_view make_scene_view(const scene_camera& camera,
    const vec4i& viewport, const glscene_params& params) {
  auto camera_aspect = (float)viewport.z / (float)viewport.w;
  auto camera_yfov =
      camera_aspect >= 0
          ? (2 * atan(camera.film / (camera_aspect * 2 * camera.lens)))
          : (2 * atan(camera.film / (2 * camera.lens)));
  auto view_matrix       = frame_to_mat(inverse(camera.frame));
  auto projection_matrix = perspective_mat(
      camera_yfov, camera_aspect, params.near, params.far);

  auto view              = shade_view{};
  view.camera_frame      = camera.frame;
  view.view_matrix       = view_matrix;
  view.projection_matrix = projection_matrix;
  return view;
}

void draw_scene(glscene_state& glscene, const scene_model& scene,
    const vec4i& viewport, const glscene_params& params) {
  clear_ogl_framebuffer(params.background);
  set_ogl_viewport(viewport);

  auto& camera = scene.cameras.at(params.camera);
  auto  view   = make_scene_view(camera, viewport, params);
  draw_instances(glscene, scene, view, params);
  draw_environments(glscene, view, params);
}

// image based lighting

// Using 6 render passes, precompute a cubemap given a sampler for the
// environment. The input sampler can be either a cubemap or a latlong texture.
template <typename Sampler>
static void precompute_cubemap(ogl_cubemap& cubemap, const Sampler& environment,
    ogl_program& program, int size, int num_mipmap_levels = 1) {
  // init cubemap with no data
  set_cubemap(cubemap, size, 3,
      array<float*, 6>{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
      true, true, true);
  auto cube = ogl_shape{};
  set_cube_shape(cube);

  auto framebuffer = ogl_framebuffer{};
  set_framebuffer(framebuffer, {size, size});

  auto cameras = array<frame3f, 6>{
      lookat_frame({0, 0, 0}, {1, 0, 0}, {0, 1, 0}),
      lookat_frame({0, 0, 0}, {-1, 0, 0}, {0, 1, 0}),
      lookat_frame({0, 0, 0}, {0, -1, 0}, {0, 0, -1}),
      lookat_frame({0, 0, 0}, {0, 1, 0}, {0, 0, 1}),
      lookat_frame({0, 0, 0}, {0, 0, -1}, {0, 1, 0}),
      lookat_frame({0, 0, 0}, {0, 0, 1}, {0, 1, 0}),
  };

  bind_framebuffer(framebuffer);
  bind_program(program);
  for (int mipmap_level = 0; mipmap_level < num_mipmap_levels; mipmap_level++) {
    // resize render buffer and viewport
    set_framebuffer(framebuffer, {size, size});
    set_ogl_viewport(vec2i{size, size});

    for (auto i = 0; i < 6; ++i) {
      // perspective_mat(fov, aspect, near, far)
      auto camera_proj = perspective_mat(radians(90.0f), 1, 1, 100);
      auto camera_view = frame_to_mat(inverse(cameras[i]));

      set_framebuffer_texture(framebuffer, cubemap, i, mipmap_level);
      clear_ogl_framebuffer({0, 0, 0, 0}, true);

      set_uniform(program, "view", camera_view);
      set_uniform(program, "projection", camera_proj);
      set_uniform(program, "eye", vec3f{0, 0, 0});
      set_uniform(program, "mipmap_level", mipmap_level);
      set_uniform(program, "environment", environment, 0);

      draw_shape(cube);
    }
    size /= 2;
  }
  unbind_program();
  unbind_framebuffer();
}

// Using 6 render passes, precompute a cubemap given a sampler for the
// environment. The input sampler can be either a cubemap or a latlong texture.
static void precompute_cubemap(ogl_cubemap& cubemap,
    const glscene_texture& environment, ogl_program& program, int size,
    int num_mipmap_levels = 1) {
  // init cubemap with no data
  set_cubemap(cubemap, size, 3,
      array<float*, 6>{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
      true, true, true);
  auto cube = ogl_shape{};
  set_cube_shape(cube);

  auto framebuffer = ogl_framebuffer{};
  set_framebuffer(framebuffer, {size, size});

  auto cameras = array<frame3f, 6>{
      lookat_frame({0, 0, 0}, {1, 0, 0}, {0, 1, 0}),
      lookat_frame({0, 0, 0}, {-1, 0, 0}, {0, 1, 0}),
      lookat_frame({0, 0, 0}, {0, -1, 0}, {0, 0, -1}),
      lookat_frame({0, 0, 0}, {0, 1, 0}, {0, 0, 1}),
      lookat_frame({0, 0, 0}, {0, 0, -1}, {0, 1, 0}),
      lookat_frame({0, 0, 0}, {0, 0, 1}, {0, 1, 0}),
  };

  bind_framebuffer(framebuffer);
  bind_program(program);
  for (int mipmap_level = 0; mipmap_level < num_mipmap_levels; mipmap_level++) {
    // resize render buffer and viewport
    set_framebuffer(framebuffer, {size, size});
    set_ogl_viewport(vec2i{size, size});

    for (auto i = 0; i < 6; ++i) {
      // perspective_mat(fov, aspect, near, far)
      auto camera_proj = perspective_mat(radians(90.0f), 1, 1, 100);
      auto camera_view = frame_to_mat(inverse(cameras[i]));

      set_framebuffer_texture(framebuffer, cubemap, i, mipmap_level);
      clear_ogl_framebuffer({0, 0, 0, 0}, true);

      set_uniform(program, "view", camera_view);
      set_uniform(program, "projection", camera_proj);
      set_uniform(program, "eye", vec3f{0, 0, 0});
      set_uniform(program, "mipmap_level", mipmap_level);

      glActiveTexture(GL_TEXTURE0 + 0);
      glBindTexture(GL_TEXTURE_2D, environment.texture);
      glUniform1i(glGetUniformLocation(program.program_id, "environment"), 0);

      draw_shape(cube);
    }
    size /= 2;
  }
  unbind_program();
  unbind_framebuffer();
}

static void precompute_brdflut(ogl_texture& texture) {
  auto size        = vec2i{512, 512};
  auto screen_quad = ogl_shape{};
  set_quad_shape(screen_quad);

  auto program = ogl_program{};
  set_program(program, precompute_brdflut_vertex(),
      precompute_brdflut_fragment(), true);

  set_texture(
      texture, size.x, size.y, 3, (float*)nullptr, true, true, false, false);

  auto framebuffer = ogl_framebuffer{};
  set_framebuffer(framebuffer, size);
  set_framebuffer_texture(framebuffer, texture, 0);

  bind_framebuffer(framebuffer);
  bind_program(program);

  set_ogl_viewport(size);
  clear_ogl_framebuffer({0, 0, 0, 0}, true);

  draw_shape(screen_quad);

  unbind_program();
  unbind_framebuffer();
}

static void init_environment(
    glscene_state& scene, glscene_environment& environment) {
  // init drawing data
  if (environment.envlight_cubemap == glinvalid_handle) {
    scene.envlight_cubemaps.emplace_back();
    environment.envlight_cubemap = (int)scene.envlight_cubemaps.size() - 1;
  }
  if (environment.envlight_shape == glinvalid_handle) {
    scene.envlight_shapes.emplace_back();
    environment.envlight_shape = (int)scene.envlight_shapes.size() - 1;
  }

  // init program and shape for drawing the environment
  set_cube_shape(scene.envlight_shapes[environment.envlight_shape]);

  // precompute cubemap from environment texture
  auto size    = scene.textures[environment.emission_tex].height;
  auto program = ogl_program{};
  set_program(program, precompute_cubemap_vertex(),
      precompute_environment_fragment(), true);
  precompute_cubemap(scene.envlight_cubemaps[environment.envlight_cubemap],
      scene.textures[environment.emission_tex], program, size, 1);
}

void init_envlight(glscene_state& scene, glscene_environment& environment) {
  // init drawing data
  if (environment.envlight_diffuse_ == glinvalid_handle) {
    scene.envlight_diffuses.emplace_back();
    environment.envlight_diffuse_ = (int)scene.envlight_diffuses.size() - 1;
  }
  if (environment.envlight_cubemap == glinvalid_handle) {
    scene.envlight_cubemaps.emplace_back();
    environment.envlight_cubemap = (int)scene.envlight_cubemaps.size() - 1;
  }
  if (environment.envlight_specular_ == glinvalid_handle) {
    scene.envlight_speculars.emplace_back();
    environment.envlight_specular_ = (int)scene.envlight_speculars.size() - 1;
  }
  if (environment.envlight_brdflut_ == glinvalid_handle) {
    scene.envlight_brdfluts.emplace_back();
    environment.envlight_brdflut_ = (int)scene.envlight_brdfluts.size() - 1;
  }
  // precompute irradiance map
  auto diffuse_program = ogl_program{};
  set_program(diffuse_program, precompute_cubemap_vertex(),
      precompute_irradiance_fragment(), true);
  precompute_cubemap(scene.envlight_diffuses[environment.envlight_diffuse_],
      scene.envlight_cubemaps[environment.envlight_cubemap], diffuse_program,
      64);

  // precompute specular map
  auto specular_program = ogl_program{};
  set_program(specular_program, precompute_cubemap_vertex(),
      precompute_reflections_fragment(), true);
  precompute_cubemap(scene.envlight_speculars[environment.envlight_specular_],
      scene.envlight_cubemaps[environment.envlight_cubemap], specular_program,
      256, 6);

  // precompute lookup texture for specular brdf
  precompute_brdflut(scene.envlight_brdfluts[environment.envlight_brdflut_]);
}

#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverlength-strings"
#endif

static const char* shade_instance_vertex() {
  static const char* code =
      R"(
#version 330

layout(location = 0) in vec3 positions;           // vertex position (in mesh coordinate frame)
layout(location = 1) in vec3 normals;             // vertex normal (in mesh coordinate frame)
layout(location = 2) in vec2 texcoords;           // vertex texcoords
layout(location = 3) in vec4 colors;              // vertex color
layout(location = 4) in vec4 tangents;            // vertex tangent space

uniform mat4 frame;             // shape transform
uniform mat4 frameit;           // shape transform
uniform float offset;           // shape normal offset

uniform mat4 view;              // inverse of the camera frame (as a matrix)
uniform mat4 projection;        // camera projection

out vec3 position;              // [to fragment shader] vertex position (in world coordinate)
out vec3 normal;                // [to fragment shader] vertex normal (in world coordinate)
out vec2 texcoord;              // [to fragment shader] vertex texture coordinates
out vec4 color;                 // [to fragment shader] vertex color
out vec4 tangsp;                // [to fragment shader] vertex tangent space

// main function
void main() {
  // copy values
  position = positions;
  normal = normals;
  tangsp = tangents;
  texcoord = texcoords;
  color = colors;

  // normal offset
  if(offset != 0) {
    position += offset * normal;
  }

  // world projection
  position = (frame * vec4(position,1)).xyz;
  normal = (frameit * vec4(normal,0)).xyz;
  tangsp.xyz = (frame * vec4(tangsp.xyz,0)).xyz;

  // clip
  gl_Position = projection * view * vec4(position,1);
}
)";
  return code;
}

[[maybe_unused]] static const char* shade_instanced_vertex() {
  static const char* code = R"(
#version 330

layout(location = 0) in vec3 positions;
layout(location = 1) in vec3 normals;
layout(location = 2) in vec2 texcoords;
layout(location = 3) in vec4 colors;
layout(location = 4) in vec4 tangents;
layout(location = 5) in vec3 instance_from;
layout(location = 6) in vec3 instance_to;

uniform mat4  frame;
uniform mat4  frameit;
uniform float offset = 0;

uniform mat4 view;
uniform mat4 projection;

out vec3 position;
out vec3 normal;
out vec2 texcoord;
out vec4 color;
out vec4 tangsp;

// main function
void main() {
  // copy values
  position = positions;
  normal   = normals;
  tangsp   = tangents;
  texcoord = texcoords;
  color    = colors;

  // normal offset
  if (offset != 0) {
    position += offset * normal;
  }

  // world projection
  position   = (frame * vec4(position, 1)).xyz;
  normal     = (frameit * vec4(normal, 0)).xyz;
  tangsp.xyz = (frame * vec4(tangsp.xyz, 0)).xyz;

  if (instance_from != instance_to) {
    vec3 dir = instance_to - instance_from;

    vec3 up = abs(dir.z) < 0.999 ? vec3(0.0, 0.0, 1.0) : vec3(1.0, 0.0, 0.0);
    vec3 tangent   = normalize(cross(up, dir));
    vec3 bitangent = normalize(cross(dir, tangent));

    mat3 mat;
    mat[2]    = dir;
    mat[0]    = tangent;
    mat[1]    = bitangent;
    position  = mat * position;
    normal    = mat * normal;
    tangent   = mat * tangent;
    bitangent = mat * bitangent;
  }
  position += instance_from;

  // clip
  gl_Position = projection * view * vec4(position, 1);
}
)";
  return code;
}

static const char* shade_instance_fragment() {
  static const char* code =
      R"(
#version 330

in vec3 position;  // [from vertex shader] position in world space
in vec3 normal;    // [from vertex shader] normal in world space
in vec2 texcoord;  // [from vertex shader] texcoord
in vec4 color;     // [from vertex shader] color
in vec4 tangsp;    // [from vertex shader] tangent space

uniform int element;
uniform bool unlit;
uniform bool faceted;
uniform vec4 highlight;
uniform bool double_sided;

uniform vec3 emission;            // material ke
uniform vec3 diffuse;             // material kd
uniform vec3 specular;            // material ks
uniform float roughness;          // material rs
uniform float opacity;            // material op

uniform bool emission_tex_on;     // material ke texture on
uniform sampler2D emission_tex;   // material ke texture
uniform bool diffuse_tex_on;      // material kd texture on
uniform sampler2D diffuse_tex;    // material kd texture
uniform bool specular_tex_on;     // material ks texture on
uniform sampler2D specular_tex;   // material ks texture
uniform bool roughness_tex_on;    // material rs texture on
uniform sampler2D roughness_tex;  // material rs texture
uniform bool opacity_tex_on;      // material op texture on
uniform sampler2D opacity_tex;    // material op texture
uniform bool normalmap_tex_on;    // material normal texture on
uniform sampler2D normalmap_tex;  // material normal texture

uniform int  lighting;            // eyelight shading
uniform vec3 ambient;             // ambient light
uniform int  lights_num;                // number of lights
uniform int  lights_type[16];     // light type (0 -> point, 1 -> directional)
uniform vec3 lights_position[16];            // light positions
uniform vec3 lights_emission[16]; // light intensities

// precomputed textures for image based lighting
uniform vec3        envlight_scale;
uniform samplerCube envlight_irradiance;
uniform samplerCube envlight_reflection;
uniform sampler2D   envlight_brdflut;

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
  // color?
  shade_brdf brdf;
  brdf.emission  = eval_brdf_color(emission, emission_tex, emission_tex_on);
  brdf.diffuse   = eval_brdf_color(diffuse, diffuse_tex, diffuse_tex_on);
  brdf.specular  = eval_brdf_color(specular, specular_tex, specular_tex_on);
  brdf.roughness = eval_brdf_value(roughness, roughness_tex, roughness_tex_on);
  brdf.opacity   = eval_brdf_value(opacity, opacity_tex, opacity_tex_on);
  vec3 base = brdf.diffuse;
  float metallic = brdf.specular.x;
  brdf.diffuse = base * (1 - metallic);
  brdf.specular = base * metallic + vec3(0.04) * (1 - metallic);
  brdf.roughness = brdf.roughness * brdf.roughness;
  return brdf;
}

void eval_light(int lid, vec3 position, out vec3 radiance, out vec3 incoming) {
  radiance = vec3(0,0,0);
  incoming = vec3(0,0,0);
  if(lights_type[lid] == 0) {
    // compute point light color at position
    radiance = lights_emission[lid] / pow(length(lights_position[lid]-position),2);
    // compute light direction at position
    incoming = normalize(lights_position[lid]-position);
  }
  else if(lights_type[lid] == 1) {
    // compute light color
    radiance = lights_emission[lid];
    // compute light direction
    incoming = normalize(lights_position[lid]);
  }
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
  vec3 texture = 2 * pow(texture(normalmap_tex,texcoord).xyz, vec3(1/2.2)) - 1;
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
  }

  if (element == element_lines) {
    vec3 tangent = normalize(normal);
    norm         = normalize(outgoing - tangent * dot(outgoing, tangent));
  }

  // apply normal map
  norm = apply_normal_map(texcoord, norm, tangsp);

  // use faceforward to ensure the normals points toward us
  if (double_sided) norm = faceforward(norm, -outgoing, norm);
  return norm;
}
    
vec3 sample_prefiltered_refleciton(vec3 incoming, float roughness) {
  int   MAX_REFLECTION_LOD = 5;
  float lod                = sqrt(roughness) * MAX_REFLECTION_LOD;
  return textureLod(envlight_reflection, incoming, lod).rgb;
}

#define lighting_eyelight 0
#define lighting_camlight 1
#define lighting_envlight 2

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
        vec3 cl = vec3(0,0,0); vec3 incoming = vec3(0,0,0);
        eval_light(lid, position, cl, incoming);
        radiance += cl * eval_brdfcos(brdf, n, incoming, outgoing);
      }
    }
    if (lighting == lighting_envlight) {
      // diffuse
      radiance += brdf.diffuse * envlight_scale * textureLod(envlight_irradiance, n, 0).rgb;
      // specular
      vec3 incoming   = normalize(reflect(-outgoing, n));
      vec3 reflection = envlight_scale * sample_prefiltered_refleciton(incoming, brdf.roughness);
      vec2 env_brdf   = texture(envlight_brdflut, vec2(max(dot(n, outgoing), 0.0), roughness)).rg;
      radiance += reflection * (brdf.specular * env_brdf.x + env_brdf.y);
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
  return code;
}

#if 0
static const char* shade_envlight_fragment() {
  static const char* code = R"(
#version 330

in vec3 position;  // [from vertex shader] position in world space
in vec3 normal;    // [from vertex shader] normal in world space
in vec2 texcoord;  // [from vertex shader] texcoord
in vec4 color;     // [from vertex shader] color
in vec4 tangsp;    // [from vertex shader] tangent space

uniform int  element;
uniform bool unlit;
uniform bool faceted;
uniform vec4 highlight;
uniform bool double_sided;

uniform vec3  emission;   // material ke
uniform vec3  diffuse;    // material kd
uniform vec3  specular;   // material ks
uniform float roughness;  // material rs
uniform float opacity;    // material op

uniform bool      emission_tex_on;   // material ke texture on
uniform sampler2D emission_tex;      // material ke texture
uniform bool      diffuse_tex_on;    // material kd texture on
uniform sampler2D diffuse_tex;       // material kd texture
uniform bool      specular_tex_on;   // material ks texture on
uniform sampler2D specular_tex;      // material ks texture
uniform bool      roughness_tex_on;  // material rs texture on
uniform sampler2D roughness_tex;     // material rs texture
uniform bool      opacity_tex_on;    // material op texture on
uniform sampler2D opacity_tex;       // material op texture
uniform bool      normalmap_tex_on;  // material normal texture on
uniform sampler2D normalmap_tex;     // material normal texture

// precomputed textures for image based lighting
uniform samplerCube envlight_irradiance;
uniform samplerCube envlight_reflection;
uniform sampler2D   envlight_brdflut;

uniform mat4 frame;    // shape transform
uniform mat4 frameit;  // shape transform

uniform vec3 eye;         // camera position
uniform mat4 view;        // inverse of the camera frame (as a matrix)
uniform mat4 projection;  // camera projection

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
  shade_brdf brdf;
  brdf.emission  = eval_brdf_color(emission, emission_tex, emission_tex_on);
  brdf.diffuse   = eval_brdf_color(diffuse, diffuse_tex, diffuse_tex_on);
  brdf.specular  = eval_brdf_color(specular, specular_tex, specular_tex_on);
  brdf.roughness = eval_brdf_value(roughness, roughness_tex, roughness_tex_on);
  brdf.opacity   = eval_brdf_value(opacity, opacity_tex, opacity_tex_on);
  vec3 base = brdf.diffuse;
  float metallic = brdf.specular.x;
  brdf.diffuse = base * (1 - metallic);
  brdf.specular = base * metallic + vec3(0.04) * (1 - metallic);
  brdf.roughness = brdf.roughness * brdf.roughness;
  return brdf;
}

vec3 apply_normal_map(vec2 texcoord, vec3 normal, vec4 tangsp) {
  if (!normalmap_tex_on) return normal;
  vec3 tangu = normalize((frame * vec4(normalize(tangsp.xyz), 0)).xyz);
  vec3 tangv = normalize(cross(normal, tangu));
  if (tangsp.w < 0) tangv = -tangv;
  vec3 texture = 2 * pow(texture(normalmap_tex, texcoord).xyz, vec3(1 / 2.2)) -
                 1;
  // texture.y = -texture.y;
  return normalize(tangu * texture.x + tangv * texture.y + normal * texture.z);
}

vec3 triangle_normal(vec3 position) {
  vec3 fdx = dFdx(position);
  vec3 fdy = dFdy(position);
  return normalize((frame * vec4(normalize(cross(fdx, fdy)), 0)).xyz);
}

#define element_points 1
#define element_lines 2
#define element_triangles 3
#define etype_quads 3

vec3 compute_normal(vec3 V) {
  vec3 N;
  if (element == element_triangles) {
    if (faceted) {
      N = triangle_normal(position);
    } else {
      N = normalize(normal);
    }
  }

  if (element == element_lines) {
    // normal of lines is coplanar with view vector and direction tangent to the
    // line
    vec3 tangent = normalize(normal);
    N            = normalize(V - tangent * dot(V, tangent));
  }

  // apply normal map
  N = apply_normal_map(texcoord, N, tangsp);

  // use faceforward to ensure the normals points toward us
  if (double_sided) N = faceforward(N, -V, N);
  return N;
}

vec3 sample_prefiltered_refleciton(vec3 L, float roughness) {
  int   MAX_REFLECTION_LOD = 5;
  float lod                = sqrt(roughness) * MAX_REFLECTION_LOD;
  return textureLod(envlight_reflection, L, lod).rgb;
}

// main
void main() {
  vec3 V = normalize(eye - position);
  vec3 N = compute_normal(V);

  shade_brdf brdf = eval_brdf();
  if (brdf.opacity < 0.005) discard;

  if(unlit) {
    frag_color = vec4(brdf.emission + brdf.diffuse, brdf.opacity);
    return; 
  }

  // emission
  vec3 radiance = brdf.emission;

  // diffuse
  radiance += brdf.diffuse * textureLod(envlight_irradiance, N, 0).rgb;

  // specular
  vec3 L          = normalize(reflect(-V, N));
  vec3 reflection = sample_prefiltered_refleciton(L, brdf.roughness);
  vec2 env_brdf   = texture(envlight_brdflut, vec2(max(dot(N, V), 0.0), roughness)).rg;
  radiance += reflection * (brdf.specular * env_brdf.x + env_brdf.y);

  // final color correction
  radiance = pow(radiance * pow(2, exposure), vec3(1 / gamma));

  // output final color by setting gl_FragColor
  frag_color = vec4(radiance, brdf.opacity);
}
)";
  return code;
}

#endif

static const char* precompute_brdflut_vertex() {
  static const char* code = R"(
#version 330

layout(location = 0) in vec3 positions;  // vertex position

out vec3 position;  // vertex position (in world coordinate)

// main function
void main() {
  position = positions;

  gl_Position = vec4(position, 1);
}
)";
  return code;
}

static const char* precompute_brdflut_fragment() {
  static const char* code = R"(
#version 330

out vec3 frag_color;

in vec3 position;  //  position in world space

const float pif = 3.14159265359;

float radical_inverse(uint bits) {
  bits = (bits << 16u) | (bits >> 16u);
  bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
  bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
  bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
  bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
  return float(bits) * 2.3283064365386963e-10;  // / 0x100000000
}

vec2 hammersley(uint i, uint N) {
  return vec2(float(i) / float(N), radical_inverse(i));
}

float geometry_schlick_ggx(float NdotV, float roughness) {
  float a = roughness;
  float k = (a * a) / 2.0;

  float nom   = NdotV;
  float denom = NdotV * (1.0 - k) + k;

  return nom / denom;
}

float geometry_smith(vec3 N, vec3 V, vec3 L, float roughness) {
  float NdotV = max(dot(N, V), 0.0);
  float NdotL = max(dot(N, L), 0.0);
  float ggx2  = geometry_schlick_ggx(NdotV, roughness);
  float ggx1  = geometry_schlick_ggx(NdotL, roughness);

  return ggx1 * ggx2;
}

vec3 importance_sample_ggx(vec2 Xi, vec3 N, float roughness) {
  float a = roughness * roughness;

  float phi      = 2.0 * pif * Xi.x;
  float cosTheta = sqrt((1.0 - Xi.y) / (1.0 + (a * a - 1.0) * Xi.y));
  float sinTheta = sqrt(1.0 - cosTheta * cosTheta);

  // from spherical coordinates to cartesian coordinates
  vec3 H;
  H.x = cos(phi) * sinTheta;
  H.y = sin(phi) * sinTheta;
  H.z = cosTheta;

  // from tangent-space vector to world-space sample vector
  vec3 up        = abs(N.z) < 0.999 ? vec3(0.0, 0.0, 1.0) : vec3(1.0, 0.0, 0.0);
  vec3 tangent   = normalize(cross(up, N));
  vec3 bitangent = cross(N, tangent);

  vec3 sampleVec = tangent * H.x + bitangent * H.y + N * H.z;
  return normalize(sampleVec);
}

vec2 integrate_brdf(float NdotV, float roughness) {
  vec3 V;
  V.x = sqrt(1.0 - NdotV * NdotV);
  V.y = 0.0;
  V.z = NdotV;

  float A = 0.0;
  float B = 0.0;

  vec3 N = vec3(0.0, 0.0, 1.0);

  const uint SAMPLE_COUNT = 1024u;
  for (uint i = 0u; i < SAMPLE_COUNT; ++i) {
    vec2 Xi = hammersley(i, SAMPLE_COUNT);
    vec3 H  = importance_sample_ggx(Xi, N, roughness);
    vec3 L  = normalize(2.0 * dot(V, H) * H - V);

    float NdotL = max(L.z, 0.0);
    float NdotH = max(H.z, 0.0);
    float VdotH = max(dot(V, H), 0.0);

    if (NdotL > 0.0) {
      float G     = geometry_smith(N, V, L, roughness);
      float G_Vis = (G * VdotH) / (NdotH * NdotV);
      float Fc    = pow(1.0 - VdotH, 5.0);

      A += (1.0 - Fc) * G_Vis;
      B += Fc * G_Vis;
    }
  }
  A /= float(SAMPLE_COUNT);
  B /= float(SAMPLE_COUNT);
  return vec2(A, B);
}

void main() {
  vec2 uv              = position.xy * 0.5 + 0.5;
  vec2 integrated_brdf = integrate_brdf(uv.x, uv.y);
  frag_color           = vec3(integrated_brdf, 0.0);
}
)";
  return code;
}

static const char* precompute_cubemap_vertex() {
  static const char* code = R"(
#version 330

layout(location = 0) in vec3 positions;  // vertex position

uniform mat4 view;        // inverse of the camera frame (as a matrix)
uniform mat4 projection;  // camera projection

out vec3 position;  // vertex position (in world coordinate)
out vec2 texcoord;  // vertex texture coordinates

// main function
void main() {
  // copy values
  position = positions;

  // clip
  vec3 view_no_transform = (view * vec4(position * 100.0, 0)).xyz;
  gl_Position            = projection * vec4(view_no_transform, 1);
}
)";
  return code;
}

static const char* shade_environment_fragment() {
  static const char* code = R"(
#version 330

out vec3 frag_color;

in vec3 position;  //  position in world space

uniform vec3 eye;         // camera position
uniform mat4 view;        // inverse of the camera frame (as a matrix)
uniform mat4 projection;  // camera projection

uniform float exposure;
uniform float gamma;

uniform vec3 emission;
uniform samplerCube emission_tex;

void main() {
  vec3 direction = normalize(position);
  vec3 radiance = emission * texture(emission_tex, direction).rgb;

  // final color correction
  radiance   = pow(radiance * pow(2, exposure), vec3(1 / gamma));
  frag_color = radiance;
}
)";
  return code;
}

static const char* precompute_environment_fragment() {
  static const char* code = R"(
#version 330

out vec3 frag_color;

in vec3 position;  //  position in world space

uniform vec3 eye;         // camera position
uniform mat4 view;        // inverse of the camera frame (as a matrix)
uniform mat4 projection;  // camera projection

uniform sampler2D environment;

vec2 sample_spherical_map(vec3 v) {
  vec2 uv = vec2(atan(v.z, v.x), asin(v.y));
  uv *= vec2(0.1591, 0.3183);  // inv atan
  uv += 0.5;
  uv.x = 1 - uv.x;
  return uv;
}

void main() {
  vec3 normal = normalize(position);
  vec2 uv     = sample_spherical_map(normal);
  vec3 color  = texture(environment, uv).rgb;
  frag_color = color;
}
)";
  return code;
}

static const char* precompute_irradiance_fragment() {
  static const char* code = R"(
#version 330

out vec3 frag_color;

in vec3 position;  //  position in world space

uniform vec3 eye;         // camera position
uniform mat4 view;        // inverse of the camera frame (as a matrix)
uniform mat4 projection;  // camera projection

uniform samplerCube environment;

const float pif = 3.14159265359;

vec3 direction(float phi, float theta) {
  return vec3(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
}

void main() {
  vec3 normal = normalize(position);
  // TODO: Why do we need these flips?
  normal.z *= -1;
  normal.y *= -1;

  vec3 up    = vec3(0.0, 1.0, 0.0);
  vec3 right = normalize(cross(up, normal));
  up         = normalize(cross(normal, right));
  mat3 rot   = mat3(right, up, normal);

  vec3 irradiance = vec3(0.0);

  int phi_samples   = 256;
  int theta_samples = 128;
  for (int x = 0; x < phi_samples; x++) {
    for (int y = 0; y < theta_samples; y++) {
      float phi    = (2.0 * pif * x) / phi_samples;
      float theta  = (0.5 * pif * y) / theta_samples;
      vec3  sample = rot * direction(phi, theta);
      // TODO: Artifacts on Mac if we don't force the LOD.
      vec3 environment = textureLod(environment, sample, 0).rgb;
      irradiance += environment * cos(theta) * sin(theta);
    }
  }
  irradiance *= pif / (phi_samples * theta_samples);
  frag_color = irradiance;
}
)";
  return code;
}

static const char* precompute_reflections_fragment() {
  static const char* code = R"(
#version 330

out vec3 frag_color;

in vec3 position;  //  position in world space

uniform vec3 eye;         // camera position
uniform mat4 view;        // inverse of the camera frame (as a matrix)
uniform mat4 projection;  // camera projection
uniform int  mipmap_level;
uniform int  num_samples = 1024;

uniform samplerCube environment;

const float pif = 3.14159265359;

float radical_inverse(uint bits) {
  bits = (bits << 16u) | (bits >> 16u);
  bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
  bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
  bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
  bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
  return float(bits) * 2.3283064365386963e-10;  // / 0x100000000
}

vec2 hammersley(uint i, int N) {
  return vec2(float(i) / float(N), radical_inverse(i));
}

vec3 sample_ggx(vec2 rn, vec3 N, float roughness) {
  float a = roughness * roughness;

  float phi       = 2.0 * pif * rn.x;
  float cos_theta = sqrt((1.0 - rn.y) / (1.0 + (a * a - 1.0) * rn.y));
  float sin_theta = sqrt(1.0 - cos_theta * cos_theta);

  // from spherical coordinates to cartesian coordinates
  vec3 H;
  H.x = cos(phi) * sin_theta;
  H.y = sin(phi) * sin_theta;
  H.z = cos_theta;

  // from tangent-space vector to world-space sample vector
  vec3 up        = abs(N.z) < 0.999 ? vec3(0.0, 0.0, 1.0) : vec3(1.0, 0.0, 0.0);
  vec3 tangent   = normalize(cross(up, N));
  vec3 bitangent = normalize(cross(N, tangent));

  vec3 result = tangent * H.x + bitangent * H.y + N * H.z;
  return normalize(result);
}

void main() {
  vec3 N = normalize(position);
  N.z *= -1;
  N.y *= -1;
  vec3 R = N;
  vec3 V = N;

  float roughness = float(mipmap_level) / 5.0;
  roughness *= roughness;

  float total_weight = 0.0;
  vec3  result       = vec3(0.0);
  for (uint i = 0u; i < uint(num_samples); i++) {
    vec2  rn    = hammersley(i, num_samples);
    vec3  H     = sample_ggx(rn, N, roughness);
    vec3  L     = normalize(reflect(-V, H));
    float NdotL = dot(N, L);
    if (NdotL > 0.0) {
      result += textureLod(environment, L, 0).rgb * NdotL;
      total_weight += NdotL;
    }
  }
  result = result / total_weight;

  frag_color = vec3(result);
}
)";
  return code;
}

#ifndef _WIN32
#pragma GCC diagnostic pop
#endif

}  // namespace yocto
