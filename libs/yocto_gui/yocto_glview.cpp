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

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include <GLFW/glfw3.h>

#include "ext/imgui/imgui.h"
#include "ext/imgui/imgui_impl_glfw.h"
#include "ext/imgui/imgui_impl_opengl3.h"
#include "ext/imgui/imgui_internal.h"

#ifdef _WIN32
#undef near
#undef far
#endif

// -----------------------------------------------------------------------------
// VIEW HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

static void update_image_params(const glinput_state& input,
    const color_image& image, glimage_params& glparams) {
  glparams.window                           = input.window_size;
  glparams.framebuffer                      = input.framebuffer_viewport;
  std::tie(glparams.center, glparams.scale) = camera_imview(glparams.center,
      glparams.scale, {image.width, image.height}, glparams.window,
      glparams.fit);
}

static bool uiupdate_image_params(
    const glinput_state& input, glimage_params& glparams) {
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
    const glinput_state& input, scene_camera& camera) {
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
    const glinput_state& input, float& exposure, bool& filmic) {
  auto edited = 0;
  if (begin_glheader("tonemap")) {
    edited += draw_glslider("exposure", exposure, -5, 5);
    edited += draw_glcheckbox("filmic", filmic);
    end_glheader();
  }
  return (bool)edited;
}

static bool draw_image_inspector(const glinput_state& input,
    const color_image& image, const color_image& display,
    glimage_params& glparams) {
  if (begin_glheader("inspect")) {
    draw_glslider("zoom", glparams.scale, 0.1, 10);
    draw_glcheckbox("fit", glparams.fit);
    auto [i, j] = image_coords(input.mouse_pos, glparams.center, glparams.scale,
        {image.width, image.height});
    auto ij     = vec2i{i, j};
    draw_gldragger("mouse", ij);
    auto image_pixel   = zero4f;
    auto display_pixel = zero4f;
    if (i >= 0 && i < image.width && j >= 0 && j < image.height) {
      image_pixel   = image.pixels[j * image.width + i];
      display_pixel = image.pixels[j * image.width + i];
    }
    draw_glcoloredit("image", image_pixel);
    draw_glcoloredit("display", display_pixel);
    end_glheader();
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

static bool draw_scene_editor(scene_model& scene, scene_selection& selection,
    const function<void()>& before_edit) {
  auto edited = 0;
  if (begin_glheader("cameras")) {
    draw_glcombobox("camera", selection.camera, scene.camera_names);
    auto camera = scene.cameras.at(selection.camera);
    edited += draw_glcheckbox("ortho", camera.orthographic);
    edited += draw_glslider("lens", camera.lens, 0.001, 1);
    edited += draw_glslider("aspect", camera.aspect, 0.1, 5);
    edited += draw_glslider("film", camera.film, 0.1, 0.5);
    edited += draw_glslider("focus", camera.focus, 0.001, 100);
    edited += draw_glslider("aperture", camera.aperture, 0, 1);
    //   frame3f frame        = identity3x4f;
    if (edited) {
      if (before_edit) before_edit();
      scene.cameras.at(selection.camera) = camera;
    }
    end_glheader();
  }
  if (begin_glheader("environments")) {
    draw_glcombobox(
        "environment", selection.environment, scene.environment_names);
    auto environment = scene.environments.at(selection.environment);
    edited += draw_glcoloredithdr("emission", environment.emission);
    edited += draw_glcombobox(
        "emission_tex", environment.emission_tex, scene.texture_names, true);
    //   frame3f frame        = identity3x4f;
    if (edited) {
      if (before_edit) before_edit();
      scene.environments.at(selection.environment) = environment;
    }
    end_glheader();
  }
  if (begin_glheader("instances")) {
    draw_glcombobox("instance", selection.instance, scene.instance_names);
    auto instance = scene.instances.at(selection.instance);
    edited += draw_glcombobox("shape", instance.shape, scene.shape_names);
    edited += draw_glcombobox(
        "material", instance.material, scene.material_names);
    //   frame3f frame        = identity3x4f;
    if (edited) {
      if (before_edit) before_edit();
      scene.instances.at(selection.instance) = instance;
    }
    end_glheader();
  }
  if (begin_glheader("materials")) {
    draw_glcombobox("material", selection.material, scene.material_names);
    auto material = scene.materials.at(selection.material);
    edited += draw_glcoloredithdr("emission", material.emission);
    edited += draw_glcombobox(
        "emission_tex", material.emission_tex, scene.texture_names, true);
    edited += draw_glcoloredithdr("color", material.color);
    edited += draw_glcombobox(
        "color_tex", material.color_tex, scene.texture_names, true);
    edited += draw_glslider("roughness", material.roughness, 0, 1);
    edited += draw_glcombobox(
        "roughness_tex", material.roughness_tex, scene.texture_names, true);
    edited += draw_glslider("metallic", material.metallic, 0, 1);
    edited += draw_glslider("ior", material.ior, 0.1, 5);
    if (edited) {
      if (before_edit) before_edit();
      scene.materials.at(selection.material) = material;
    }
    end_glheader();
  }
  if (begin_glheader("shapes")) {
    draw_glcombobox("shape", selection.shape, scene.shape_names);
    auto& shape = scene.shapes.at(selection.shape);
    draw_gllabel("points", (int)shape.points.size());
    draw_gllabel("lines", (int)shape.lines.size());
    draw_gllabel("triangles", (int)shape.triangles.size());
    draw_gllabel("quads", (int)shape.quads.size());
    draw_gllabel("positions", (int)shape.positions.size());
    draw_gllabel("normals", (int)shape.normals.size());
    draw_gllabel("texcoords", (int)shape.texcoords.size());
    draw_gllabel("colors", (int)shape.colors.size());
    draw_gllabel("radius", (int)shape.radius.size());
    draw_gllabel("tangents", (int)shape.tangents.size());
    end_glheader();
  }
  if (begin_glheader("textures")) {
    draw_glcombobox("texture", selection.texture, scene.texture_names);
    auto& texture = scene.textures.at(selection.texture);
    draw_gllabel("width", texture.width);
    draw_gllabel("height", texture.height);
    draw_gllabel("linear", texture.linear);
    draw_gllabel("byte", !texture.pixelsb.empty());
    end_glheader();
  }
  if (begin_glheader("subdivs")) {
    draw_glcombobox("subdiv", selection.subdiv, scene.subdiv_names);
    auto& subdiv = scene.subdivs.at(selection.subdiv);
    draw_gllabel("quadspos", (int)subdiv.quadspos.size());
    draw_gllabel("quadsnorm", (int)subdiv.quadsnorm.size());
    draw_gllabel("quadstexcoord", (int)subdiv.quadstexcoord.size());
    draw_gllabel("positions", (int)subdiv.positions.size());
    draw_gllabel("normals", (int)subdiv.normals.size());
    draw_gllabel("texcoords", (int)subdiv.texcoords.size());
    end_glheader();
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
  auto callbacks    = glwindow_callbacks{};
  callbacks.init_cb = [&](const glinput_state& input) {
    init_image(glimage);
    set_image(glimage, display);
  };
  callbacks.clear_cb = [&](const glinput_state& input) {
    clear_image(glimage);
  };
  callbacks.draw_cb = [&](const glinput_state& input) {
    update_image_params(input, image, glparams);
    draw_image(glimage, glparams);
  };
  callbacks.widgets_cb = [&](const glinput_state& input) {
    draw_glcombobox("name", selected, names);
    if (draw_tonemap_params(input, exposure, filmic)) {
      tonemap_image_mt(display, image, exposure, filmic);
      set_image(glimage, display);
    }
    draw_image_inspector(input, image, display, glparams);
  };
  callbacks.uiupdate_cb = [&](const glinput_state& input) {
    uiupdate_image_params(input, glparams);
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
  auto callbacks    = glwindow_callbacks{};
  callbacks.init_cb = [&](const glinput_state& input) {
    for (auto idx = 0; idx < (int)images.size(); idx++) {
      init_image(glimages[idx]);
      set_image(glimages[idx], displays[idx]);
    }
  };
  callbacks.clear_cb = [&](const glinput_state& input) {
    for (auto idx = 0; idx < (int)images.size(); idx++) {
      clear_image(glimages[idx]);
    }
  };
  callbacks.draw_cb = [&](const glinput_state& input) {
    update_image_params(input, displays[selected], glparamss[selected]);
    draw_image(glimages[selected], glparamss[selected]);
  };
  callbacks.widgets_cb = [&](const glinput_state& input) {
    draw_glcombobox("name", selected, names);
    auto filmic = (bool)filmics[selected];  // vector of bool ...
    if (draw_tonemap_params(input, exposures[selected], filmic)) {
      filmics[selected] = filmic;
      tonemap_image_mt(displays[selected], images[selected],
          exposures[selected], filmics[selected]);
      set_image(glimages[selected], displays[selected]);
    }
    draw_image_inspector(
        input, images[selected], displays[selected], glparamss[selected]);
  };
  callbacks.uiupdate_cb = [&](const glinput_state& input) {
    uiupdate_image_params(input, glparamss[selected]);
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
  auto callbacks    = glwindow_callbacks{};
  callbacks.init_cb = [&](const glinput_state& input) {
    init_image(glimage);
    set_image(glimage, display);
  };
  callbacks.clear_cb = [&](const glinput_state& input) {
    clear_image(glimage);
  };
  callbacks.draw_cb = [&](const glinput_state& input) {
    update_image_params(input, image, glparams);
    draw_image(glimage, glparams);
  };
  callbacks.widgets_cb = [&](const glinput_state& input) {
    draw_glcombobox("name", selected, names);
    if (begin_glheader("colorgrade")) {
      auto edited = 0;
      edited += draw_glslider("exposure", params.exposure, -5, 5);
      edited += draw_glcoloredit("tint", params.tint);
      edited += draw_glslider("lincontrast", params.lincontrast, 0, 1);
      edited += draw_glslider("logcontrast", params.logcontrast, 0, 1);
      edited += draw_glslider("linsaturation", params.linsaturation, 0, 1);
      edited += draw_glcheckbox("filmic", params.filmic);
      continue_glline();
      edited += draw_glcheckbox("srgb", params.srgb);
      edited += draw_glslider("contrast", params.contrast, 0, 1);
      edited += draw_glslider("saturation", params.saturation, 0, 1);
      edited += draw_glslider("shadows", params.shadows, 0, 1);
      edited += draw_glslider("midtones", params.midtones, 0, 1);
      edited += draw_glslider("highlights", params.highlights, 0, 1);
      edited += draw_glcoloredit("shadows color", params.shadows_color);
      edited += draw_glcoloredit("midtones color", params.midtones_color);
      edited += draw_glcoloredit("highlights color", params.highlights_color);
      end_glheader();
      if (edited) {
        colorgrade_image_mt(display, image, params);
        set_image(glimage, display);
      }
    }
    draw_image_inspector(input, image, display, glparams);
  };
  callbacks.uiupdate_cb = [&glparams](const glinput_state& input) {
    uiupdate_image_params(input, glparams);
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
          render_current = state.samples;
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
  auto callbacks    = glwindow_callbacks{};
  callbacks.init_cb = [&](const glinput_state& input) {
    auto lock = std::lock_guard{render_mutex};
    init_image(glimage);
    set_image(glimage, display);
  };
  callbacks.clear_cb = [&](const glinput_state& input) {
    clear_image(glimage);
  };
  callbacks.draw_cb = [&](const glinput_state& input) {
    // update image
    if (render_update) {
      auto lock = std::lock_guard{render_mutex};
      set_image(glimage, display);
      render_update = false;
    }
    update_image_params(input, image, glparams);
    draw_image(glimage, glparams);
  };
  callbacks.widgets_cb = [&](const glinput_state& input) {
    auto edited = 0;
    draw_glcombobox("name", selected, names);
    auto current = (int)render_current;
    draw_glprogressbar("sample", current, params.samples);
    if (begin_glheader("render")) {
      auto edited  = 0;
      auto tparams = params;
      edited += draw_glcombobox("camera", tparams.camera, camera_names);
      edited += draw_glslider("resolution", tparams.resolution, 180, 4096);
      edited += draw_glslider("samples", tparams.samples, 16, 4096);
      edited += draw_glcombobox(
          "tracer", (int&)tparams.sampler, trace_sampler_names);
      edited += draw_glcombobox(
          "false color", (int&)tparams.falsecolor, trace_falsecolor_names);
      edited += draw_glslider("bounces", tparams.bounces, 1, 128);
      edited += draw_glslider("clamp", tparams.clamp, 10, 1000);
      edited += draw_glcheckbox("envhidden", tparams.envhidden);
      continue_glline();
      edited += draw_glcheckbox("filter", tparams.tentfilter);
      edited += draw_glslider("pratio", tparams.pratio, 1, 64);
      // edited += draw_glslider("exposure", tparams.exposure, -5, 5);
      end_glheader();
      if (edited) {
        stop_render();
        params = tparams;
        reset_display();
      }
    }
    if (begin_glheader("tonemap")) {
      edited += draw_glslider("exposure", params.exposure, -5, 5);
      edited += draw_glcheckbox("filmic", params.filmic);
      edited += draw_glcheckbox("denoise", params.denoise);
      end_glheader();
      if (edited) {
        tonemap_image_mt(display, image, params.exposure, params.filmic);
        set_image(glimage, display);
      }
    }
    draw_image_inspector(input, image, display, glparams);
    if (edit) {
      if (draw_scene_editor(scene, selection, [&]() { stop_render(); })) {
        reset_display();
      }
    }
  };
  callbacks.uiupdate_cb = [&](const glinput_state& input) {
    auto camera = scene.cameras[params.camera];
    if (uiupdate_camera_params(input, camera)) {
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
  auto callbacks    = glwindow_callbacks{};
  callbacks.init_cb = [&](const glinput_state& input) {
    init_glscene(glscene, scene);
  };
  callbacks.clear_cb = [&](const glinput_state& input) {
    clear_scene(glscene);
  };
  callbacks.draw_cb = [&](const glinput_state& input) {
    draw_scene(glscene, scene, input.framebuffer_viewport, params);
  };
  callbacks.widgets_cb = [&](const glinput_state& input) {
    draw_glcombobox("name", selected, names);
    if (begin_glheader("shade")) {
      draw_glcombobox("camera", params.camera, camera_names);
      draw_glcheckbox("wireframe", params.wireframe);
      continue_glline();
      draw_glcheckbox("faceted", params.faceted);
      continue_glline();
      draw_glcheckbox("double sided", params.double_sided);
      draw_glcombobox(
          "lighting", (int&)params.lighting, glscene_lighting_names);
      draw_glslider("exposure", params.exposure, -10, 10);
      draw_glslider("gamma", params.gamma, 0.1f, 4);
      draw_glslider("near", params.near, 0.01f, 1.0f);
      draw_glslider("far", params.far, 1000.0f, 10000.0f);
      end_glheader();
    }
    // draw_scene_editor(scene, selection, {});
    if (widgets_callback) {
      widgets_callback(input, updated_shapes, updated_textures);
      if (!updated_shapes.empty() || !updated_textures.empty()) {
        update_glscene(glscene, scene, updated_shapes, updated_textures);
        updated_shapes.clear();
        updated_textures.clear();
      }
    }
  };
  callbacks.update_cb = [&](const glinput_state& input) {
    if (update_callback) {
      update_callback(input, updated_shapes, updated_textures);
      if (!updated_shapes.empty() || !updated_textures.empty()) {
        update_glscene(glscene, scene, updated_shapes, updated_textures);
        updated_shapes.clear();
        updated_textures.clear();
      }
    }
  };
  callbacks.uiupdate_cb = [&](const glinput_state& input) {
    // handle mouse and keyboard for navigation
    if (uiupdate_callback) {
      uiupdate_callback(input, updated_shapes, updated_textures);
      if (!updated_shapes.empty() || !updated_textures.empty()) {
        update_glscene(glscene, scene, updated_shapes, updated_textures);
        updated_shapes.clear();
        updated_textures.clear();
      }
    }
    auto camera = scene.cameras.at(params.camera);
    if (uiupdate_camera_params(input, camera)) {
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

// -----------------------------------------------------------------------------
// WINDOW
// -----------------------------------------------------------------------------
namespace yocto {

// OpenGL window wrapper
struct glwindow_state {
  string              title         = "";
  init_glcallback     init_cb       = {};
  clear_glcallback    clear_cb      = {};
  draw_glcallback     draw_cb       = {};
  widgets_glcallback  widgets_cb    = {};
  update_glcallback   update_cb     = {};
  uiupdate_glcallback uiupdate_cb   = {};
  int                 widgets_width = 0;
  bool                widgets_left  = true;
  glinput_state       input         = {};
  vec2i               window_size   = {0, 0};
  vec4f               background    = {0.15f, 0.15f, 0.15f, 1.0f};
};

static void draw_window(glwindow_state& state) {
  glClearColor(state.background.x, state.background.y, state.background.z,
      state.background.w);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  if (state.draw_cb) state.draw_cb(state.input);
  if (state.widgets_cb) {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    auto window = state.window_size;
    if (state.widgets_left) {
      ImGui::SetNextWindowPos({0, 0});
      ImGui::SetNextWindowSize({(float)state.widgets_width, (float)window.y});
    } else {
      ImGui::SetNextWindowPos({(float)(window.x - state.widgets_width), 0});
      ImGui::SetNextWindowSize({(float)state.widgets_width, (float)window.y});
    }
    ImGui::SetNextWindowCollapsed(false);
    ImGui::SetNextWindowBgAlpha(1);
    if (ImGui::Begin(state.title.c_str(), nullptr,
            ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse |
                ImGuiWindowFlags_NoSavedSettings)) {
      state.widgets_cb(state.input);
    }
    ImGui::End();
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
  }
}

// run the user interface with the give callbacks
void run_ui(const vec2i& size, const string& title,
    const glwindow_callbacks& callbacks, int widgets_width, bool widgets_left) {
  // init glfw
  if (!glfwInit())
    throw std::runtime_error("cannot initialize windowing system");
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#if __APPLE__
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

  // create state
  auto state        = glwindow_state{};
  state.title       = title;
  state.init_cb     = callbacks.init_cb;
  state.clear_cb    = callbacks.clear_cb;
  state.draw_cb     = callbacks.draw_cb;
  state.widgets_cb  = callbacks.widgets_cb;
  state.update_cb   = callbacks.update_cb;
  state.uiupdate_cb = callbacks.uiupdate_cb;

  // create window
  auto window = glfwCreateWindow(
      size.x, size.y, title.c_str(), nullptr, nullptr);
  if (window == nullptr)
    throw std::runtime_error{"cannot initialize windowing system"};
  glfwMakeContextCurrent(window);
  glfwSwapInterval(1);  // Enable vsync

  // set user data
  glfwSetWindowUserPointer(window, &state);

  // set callbacks
  glfwSetWindowRefreshCallback(window, [](GLFWwindow* window) {
    auto& state = *(glwindow_state*)glfwGetWindowUserPointer(window);
    glfwGetWindowSize(window, &state.window_size.x, &state.window_size.y);
    draw_window(state);
    glfwSwapBuffers(window);
  });
  glfwSetWindowSizeCallback(
      window, [](GLFWwindow* window, int width, int height) {
        auto& state = *(glwindow_state*)glfwGetWindowUserPointer(window);
        glfwGetWindowSize(
            window, &state.input.window_size.x, &state.input.window_size.y);
        if (state.widgets_width)
          state.input.window_size.x -= state.widgets_width;
        glfwGetFramebufferSize(window, &state.input.framebuffer_viewport.z,
            &state.input.framebuffer_viewport.w);
        state.input.framebuffer_viewport.x = 0;
        state.input.framebuffer_viewport.y = 0;
        if (state.widgets_width) {
          auto win_size = zero2i;
          glfwGetWindowSize(window, &win_size.x, &win_size.y);
          auto offset = (int)(state.widgets_width *
                              (float)state.input.framebuffer_viewport.z /
                              win_size.x);
          state.input.framebuffer_viewport.z -= offset;
          if (state.widgets_left) state.input.framebuffer_viewport.x += offset;
        }
      });

  // init gl extensions
  if (!gladLoadGL())
    throw std::runtime_error{"cannot initialize OpenGL extensions"};

  // widgets
  if (callbacks.widgets_cb) {
    ImGui::CreateContext();
    ImGui::GetIO().IniFilename       = nullptr;
    ImGui::GetStyle().WindowRounding = 0;
    ImGui_ImplGlfw_InitForOpenGL(window, true);
#ifndef __APPLE__
    ImGui_ImplOpenGL3_Init();
#else
    ImGui_ImplOpenGL3_Init("#version 330");
#endif
    ImGui::StyleColorsDark();
    state.widgets_width = widgets_width;
    state.widgets_left  = widgets_left;
  }

  // init
  if (state.init_cb) state.init_cb(state.input);

  // run ui
  while (!glfwWindowShouldClose(window)) {
    // update input
    state.input.mouse_last = state.input.mouse_pos;
    auto mouse_posx = 0.0, mouse_posy = 0.0;
    glfwGetCursorPos(window, &mouse_posx, &mouse_posy);
    state.input.mouse_pos = vec2f{(float)mouse_posx, (float)mouse_posy};
    if (state.widgets_width && state.widgets_left)
      state.input.mouse_pos.x -= state.widgets_width;
    state.input.mouse_left = glfwGetMouseButton(
                                 window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
    state.input.mouse_right =
        glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
    state.input.modifier_alt =
        glfwGetKey(window, GLFW_KEY_LEFT_ALT) == GLFW_PRESS ||
        glfwGetKey(window, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS;
    state.input.modifier_shift =
        glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
        glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;
    state.input.modifier_ctrl =
        glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS ||
        glfwGetKey(window, GLFW_KEY_RIGHT_CONTROL) == GLFW_PRESS;
    glfwGetWindowSize(
        window, &state.input.window_size.x, &state.input.window_size.y);
    if (state.widgets_width) state.input.window_size.x -= state.widgets_width;
    glfwGetFramebufferSize(window, &state.input.framebuffer_viewport.z,
        &state.input.framebuffer_viewport.w);
    state.input.framebuffer_viewport.x = 0;
    state.input.framebuffer_viewport.y = 0;
    if (state.widgets_width) {
      auto win_size = zero2i;
      glfwGetWindowSize(window, &win_size.x, &win_size.y);
      auto offset = (int)(state.widgets_width *
                          (float)state.input.framebuffer_viewport.z /
                          win_size.x);
      state.input.framebuffer_viewport.z -= offset;
      if (state.widgets_left) state.input.framebuffer_viewport.x += offset;
    }
    if (state.widgets_width) {
      auto io                    = &ImGui::GetIO();
      state.input.widgets_active = io->WantTextInput || io->WantCaptureMouse ||
                                   io->WantCaptureKeyboard;
    }

    // time
    state.input.clock_last = state.input.clock_now;
    state.input.clock_now =
        std::chrono::high_resolution_clock::now().time_since_epoch().count();
    state.input.time_now = (double)state.input.clock_now / 1000000000.0;
    state.input.time_delta =
        (double)(state.input.clock_now - state.input.clock_last) / 1000000000.0;

    // update ui
    if (state.uiupdate_cb && !state.input.widgets_active)
      state.uiupdate_cb(state.input);

    // update
    if (state.update_cb) state.update_cb(state.input);

    // draw
    glfwGetWindowSize(window, &state.window_size.x, &state.window_size.y);
    draw_window(state);
    glfwSwapBuffers(window);

    // event hadling
    glfwPollEvents();
  }

  // clear
  if (state.clear_cb) state.clear_cb(state.input);

  // clear
  glfwDestroyWindow(window);
  glfwTerminate();
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// OPENGL WIDGETS
// -----------------------------------------------------------------------------
namespace yocto {

bool begin_glheader(const char* lbl) {
  if (!ImGui::CollapsingHeader(lbl)) return false;
  ImGui::PushID(lbl);
  return true;
}
void end_glheader() { ImGui::PopID(); }

bool draw_glbutton(const char* lbl, bool enabled) {
  if (enabled) {
    return ImGui::Button(lbl);
  } else {
    ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
    ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
    auto ok = ImGui::Button(lbl);
    ImGui::PopItemFlag();
    ImGui::PopStyleVar();
    return ok;
  }
}

void draw_gllabel(const char* lbl, const string& label) {
  ImGui::LabelText(lbl, "%s", label.c_str());
}
void draw_gllabel(const char* lbl, int value) {
  ImGui::LabelText(lbl, "%s", std::to_string(value).c_str());
}
void draw_gllabel(const char* lbl, bool value) {
  ImGui::LabelText(lbl, "%s", value ? "true" : "false");
}

void draw_glseparator() { ImGui::Separator(); }

void continue_glline() { ImGui::SameLine(); }

bool draw_gltextinput(const char* lbl, string& value) {
  auto buffer = array<char, 4096>{};
  auto num    = 0;
  for (auto c : value) buffer[num++] = c;
  buffer[num] = 0;
  auto edited = ImGui::InputText(lbl, buffer.data(), buffer.size());
  if (edited) value = buffer.data();
  return edited;
}

bool draw_glslider(const char* lbl, float& value, float min, float max) {
  return ImGui::SliderFloat(lbl, &value, min, max);
}
bool draw_glslider(const char* lbl, vec2f& value, float min, float max) {
  return ImGui::SliderFloat2(lbl, &value.x, min, max);
}
bool draw_glslider(const char* lbl, vec3f& value, float min, float max) {
  return ImGui::SliderFloat3(lbl, &value.x, min, max);
}
bool draw_glslider(const char* lbl, vec4f& value, float min, float max) {
  return ImGui::SliderFloat4(lbl, &value.x, min, max);
}

bool draw_glslider(const char* lbl, int& value, int min, int max) {
  return ImGui::SliderInt(lbl, &value, min, max);
}
bool draw_glslider(const char* lbl, vec2i& value, int min, int max) {
  return ImGui::SliderInt2(lbl, &value.x, min, max);
}
bool draw_glslider(const char* lbl, vec3i& value, int min, int max) {
  return ImGui::SliderInt3(lbl, &value.x, min, max);
}
bool draw_glslider(const char* lbl, vec4i& value, int min, int max) {
  return ImGui::SliderInt4(lbl, &value.x, min, max);
}

bool draw_gldragger(
    const char* lbl, float& value, float speed, float min, float max) {
  return ImGui::DragFloat(lbl, &value, speed, min, max);
}
bool draw_gldragger(
    const char* lbl, vec2f& value, float speed, float min, float max) {
  return ImGui::DragFloat2(lbl, &value.x, speed, min, max);
}
bool draw_gldragger(
    const char* lbl, vec3f& value, float speed, float min, float max) {
  return ImGui::DragFloat3(lbl, &value.x, speed, min, max);
}
bool draw_gldragger(
    const char* lbl, vec4f& value, float speed, float min, float max) {
  return ImGui::DragFloat4(lbl, &value.x, speed, min, max);
}

bool draw_gldragger(
    const char* lbl, int& value, float speed, int min, int max) {
  return ImGui::DragInt(lbl, &value, speed, min, max);
}
bool draw_gldragger(
    const char* lbl, vec2i& value, float speed, int min, int max) {
  return ImGui::DragInt2(lbl, &value.x, speed, min, max);
}
bool draw_gldragger(
    const char* lbl, vec3i& value, float speed, int min, int max) {
  return ImGui::DragInt3(lbl, &value.x, speed, min, max);
}
bool draw_gldragger(
    const char* lbl, vec4i& value, float speed, int min, int max) {
  return ImGui::DragInt4(lbl, &value.x, speed, min, max);
}

bool draw_gldragger(const char* lbl, array<float, 2>& value, float speed,
    float min, float max) {
  return ImGui::DragFloat2(lbl, value.data(), speed, min, max);
}
bool draw_gldragger(const char* lbl, array<float, 3>& value, float speed,
    float min, float max) {
  return ImGui::DragFloat3(lbl, value.data(), speed, min, max);
}
bool draw_gldragger(const char* lbl, array<float, 4>& value, float speed,
    float min, float max) {
  return ImGui::DragFloat4(lbl, value.data(), speed, min, max);
}

bool draw_gldragger(
    const char* lbl, array<int, 2>& value, float speed, int min, int max) {
  return ImGui::DragInt2(lbl, value.data(), speed, min, max);
}
bool draw_gldragger(
    const char* lbl, array<int, 3>& value, float speed, int min, int max) {
  return ImGui::DragInt3(lbl, value.data(), speed, min, max);
}
bool draw_gldragger(
    const char* lbl, array<int, 4>& value, float speed, int min, int max) {
  return ImGui::DragInt4(lbl, value.data(), speed, min, max);
}

bool draw_glcheckbox(const char* lbl, bool& value) {
  return ImGui::Checkbox(lbl, &value);
}
bool draw_glcheckbox(const char* lbl, bool& value, bool invert) {
  if (!invert) {
    return draw_glcheckbox(lbl, value);
  } else {
    auto inverted = !value;
    auto edited   = ImGui::Checkbox(lbl, &inverted);
    if (edited) value = !inverted;
    return edited;
  }
}

bool draw_glcoloredit(const char* lbl, vec3f& value) {
  auto flags = ImGuiColorEditFlags_Float;
  return ImGui::ColorEdit3(lbl, &value.x, flags);
}

bool draw_glcoloredit(const char* lbl, vec4f& value) {
  auto flags = ImGuiColorEditFlags_Float;
  return ImGui::ColorEdit4(lbl, &value.x, flags);
}

bool draw_glcoloredithdr(const char* lbl, vec3f& value) {
  auto color    = value;
  auto exposure = 0.0f;
  auto scale    = max(color);
  if (scale > 1) {
    color /= scale;
    exposure = log2(scale);
  }
  auto edit_exposure = draw_glslider(
      (lbl + " [exp]"s).c_str(), exposure, 0, 10);
  auto edit_color = draw_glcoloredit((lbl + " [col]"s).c_str(), color);
  if (edit_exposure || edit_color) {
    value = color * exp2(exposure);
    return true;
  } else {
    return false;
  }
}
bool draw_glcoloredithdr(const char* lbl, vec4f& value) {
  auto color    = value;
  auto exposure = 0.0f;
  auto scale    = max(xyz(color));
  if (scale > 1) {
    color.x /= scale;
    color.y /= scale;
    color.z /= scale;
    exposure = log2(scale);
  }
  auto edit_exposure = draw_glslider(
      (lbl + " [exp]"s).c_str(), exposure, 0, 10);
  auto edit_color = draw_glcoloredit((lbl + " [col]"s).c_str(), color);
  if (edit_exposure || edit_color) {
    value.x = color.x * exp2(exposure);
    value.y = color.y * exp2(exposure);
    value.z = color.z * exp2(exposure);
    value.w = color.w;
    return true;
  } else {
    return false;
  }
}

bool draw_glcoloredit(const char* lbl, vec4b& value) {
  auto valuef = byte_to_float(value);
  if (ImGui::ColorEdit4(lbl, &valuef.x)) {
    value = float_to_byte(valuef);
    return true;
  } else {
    return false;
  }
}

bool draw_glcombobox(const char* lbl, int& value, const vector<string>& labels,
    bool include_null) {
  if (!ImGui::BeginCombo(lbl, value >= 0 ? labels.at(value).c_str() : "<none>"))
    return false;
  auto old_val = value;
  if (include_null) {
    ImGui::PushID(100000);
    if (ImGui::Selectable("<none>", value < 0)) value = -1;
    if (value < 0) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  for (auto i = 0; i < labels.size(); i++) {
    ImGui::PushID(i);
    if (ImGui::Selectable(labels[i].c_str(), value == i)) value = i;
    if (value == i) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  ImGui::EndCombo();
  return value != old_val;
}

bool draw_glcombobox(const char* lbl, string& value,
    const vector<string>& labels, bool include_null) {
  if (!ImGui::BeginCombo(lbl, value.c_str())) return false;
  auto old_val = value;
  if (include_null) {
    ImGui::PushID(100000);
    if (ImGui::Selectable("<none>", value.empty())) value = "";
    if (value.empty()) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  for (auto i = 0; i < labels.size(); i++) {
    ImGui::PushID(i);
    if (ImGui::Selectable(labels[i].c_str(), value == labels[i]))
      value = labels[i];
    if (value == labels[i]) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  ImGui::EndCombo();
  return value != old_val;
}

bool draw_glcombobox(const char* lbl, int& idx, int num,
    const function<string(int)>& labels, bool include_null) {
  if (num <= 0) idx = -1;
  if (!ImGui::BeginCombo(lbl, idx >= 0 ? labels(idx).c_str() : "<none>"))
    return false;
  auto old_idx = idx;
  if (include_null) {
    ImGui::PushID(100000);
    if (ImGui::Selectable("<none>", idx < 0)) idx = -1;
    if (idx < 0) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  for (auto i = 0; i < num; i++) {
    ImGui::PushID(i);
    if (ImGui::Selectable(labels(i).c_str(), idx == i)) idx = i;
    if (idx == i) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  ImGui::EndCombo();
  return idx != old_idx;
}

void draw_glprogressbar(const char* lbl, float fraction) {
  ImGui::PushStyleColor(ImGuiCol_PlotHistogram, ImVec4(0.5, 0.5, 1, 0.25));
  ImGui::ProgressBar(fraction, ImVec2(0.0f, 0.0f));
  ImGui::SameLine(0.0f, ImGui::GetStyle().ItemInnerSpacing.x);
  ImGui::Text(lbl, ImVec2(0.0f, 0.0f));
  ImGui::PopStyleColor(1);
}

void draw_glprogressbar(const char* lbl, int current, int total) {
  auto overlay = std::to_string(current) + "/" + std::to_string(total);
  ImGui::PushStyleColor(ImGuiCol_PlotHistogram, ImVec4(0.5, 0.5, 1, 0.25));
  ImGui::ProgressBar(
      (float)current / (float)total, ImVec2(0.0f, 0.0f), overlay.c_str());
  ImGui::SameLine(0.0f, ImGui::GetStyle().ItemInnerSpacing.x);
  ImGui::Text(lbl, ImVec2(0.0f, 0.0f));
  ImGui::PopStyleColor(1);
}

void draw_histogram(const char* lbl, const float* values, int count) {
  ImGui::PlotHistogram(lbl, values, count);
}
void draw_histogram(const char* lbl, const vector<float>& values) {
  ImGui::PlotHistogram(lbl, values.data(), (int)values.size(), 0, nullptr,
      flt_max, flt_max, {0, 0}, 4);
}
void draw_histogram(const char* lbl, const vector<vec2f>& values) {
  ImGui::PlotHistogram((lbl + " x"s).c_str(), (const float*)values.data() + 0,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec2f));
  ImGui::PlotHistogram((lbl + " y"s).c_str(), (const float*)values.data() + 1,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec2f));
}
void draw_histogram(const char* lbl, const vector<vec3f>& values) {
  ImGui::PlotHistogram((lbl + " x"s).c_str(), (const float*)values.data() + 0,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec3f));
  ImGui::PlotHistogram((lbl + " y"s).c_str(), (const float*)values.data() + 1,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec3f));
  ImGui::PlotHistogram((lbl + " z"s).c_str(), (const float*)values.data() + 2,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec3f));
}
void draw_histogram(const char* lbl, const vector<vec4f>& values) {
  ImGui::PlotHistogram((lbl + " x"s).c_str(), (const float*)values.data() + 0,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec4f));
  ImGui::PlotHistogram((lbl + " y"s).c_str(), (const float*)values.data() + 1,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec4f));
  ImGui::PlotHistogram((lbl + " z"s).c_str(), (const float*)values.data() + 2,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec4f));
  ImGui::PlotHistogram((lbl + " w"s).c_str(), (const float*)values.data() + 3,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec4f));
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// OPENGL WIDGETS
// -----------------------------------------------------------------------------
namespace yocto {

enum struct glwidgets_param_type {
  // clang-format off
  value1f, value2f, value3f, value4f, 
  value1i, value2i, value3i, value4i, 
  value1s, value1b
  // clang-format on
};

struct glwidgets_param {
  // constructors
  glwidgets_param()
      : type{glwidgets_param_type::value1f}
      , valuef{0, 0, 0, 0}
      , minmaxf{0, 0}
      , readonly{true} {}
  glwidgets_param(
      float value, const vec2f& minmax = {0, 0}, bool readonly = false)
      : type{glwidgets_param_type::value1f}
      , valuef{value, 0, 0, 0}
      , minmaxf{minmax}
      , readonly{readonly} {}
  glwidgets_param(
      vec2f value, const vec2f& minmax = {0, 0}, bool readonly = false)
      : type{glwidgets_param_type::value2f}
      , valuef{value.x, value.y, 0, 0}
      , minmaxf{minmax}
      , readonly{readonly} {}
  glwidgets_param(
      vec3f value, const vec2f& minmax = {0, 0}, bool readonly = false)
      : type{glwidgets_param_type::value3f}
      , valuef{value.x, value.y, value.z}
      , minmaxf{minmax}
      , readonly{readonly} {}
  glwidgets_param(
      vec4f value, const vec2f& minmax = {0, 0}, bool readonly = false)
      : type{glwidgets_param_type::value4f}
      , valuef{value.x, value.y, value.z, value.w}
      , minmaxf{minmax}
      , readonly{readonly} {}
  glwidgets_param(vec3f value, bool color, bool readonly = false)
      : type{glwidgets_param_type::value3f}
      , valuef{value.x, value.y, value.z, 1}
      , color{color}
      , readonly{readonly} {}
  glwidgets_param(vec4f value, bool color, bool readonly = false)
      : type{glwidgets_param_type::value4f}
      , valuef{value.x, value.y, value.z, value.w}
      , color{color}
      , readonly{readonly} {}
  glwidgets_param(
      int value, const vec2i& minmax = {0, 0}, bool readonly = false)
      : type{glwidgets_param_type::value1i}
      , valuei{value, 0, 0, 0}
      , minmaxi{minmax}
      , readonly{readonly} {}
  glwidgets_param(
      vec2i value, const vec2i& minmax = {0, 0}, bool readonly = false)
      : type{glwidgets_param_type::value2i}
      , valuei{value.x, value.y, 0, 0}
      , minmaxi{minmax}
      , readonly{readonly} {}
  glwidgets_param(
      vec3i value, const vec2i& minmax = {0, 0}, bool readonly = false)
      : type{glwidgets_param_type::value3i}
      , valuei{value.x, value.y, value.z, 0}
      , minmaxi{minmax}
      , readonly{readonly} {}
  glwidgets_param(
      vec4i value, const vec2i& minmax = {0, 0}, bool readonly = false)
      : type{glwidgets_param_type::value4i}
      , valuei{value.x, value.y, value.z, value.w}
      , minmaxi{minmax}
      , readonly{readonly} {}
  glwidgets_param(bool value, bool readonly = false)
      : type{glwidgets_param_type::value1b}
      , valueb{value}
      , readonly{readonly} {}
  glwidgets_param(const string& value, bool readonly = false)
      : type{glwidgets_param_type::value1s}
      , values{value}
      , readonly{readonly} {}
  glwidgets_param(
      const string& value, const vector<string>& labels, bool readonly = false)
      : type{glwidgets_param_type::value1s}
      , values{value}
      , labels{labels}
      , readonly{readonly} {}
  glwidgets_param(
      int value, const vector<string>& labels, bool readonly = false)
      : type{glwidgets_param_type::value1i}
      , valuei{value, 0, 0, 0}
      , labels{labels}
      , readonly{readonly} {}
  template <typename T, typename = std::enable_if_t<std::is_enum_v<T>>>
  glwidgets_param(T value, const vector<string>& labels, bool readonly = false)
      : type{glwidgets_param_type::value1i}
      , valuei{(int)value, 0, 0, 0}
      , labels{labels}
      , readonly{readonly} {}

  // conversions
  operator float() const {
    check_type(glwidgets_param_type::value1f);
    return valuef.x;
  }
  operator vec2f() const {
    check_type(glwidgets_param_type::value2f);
    return {valuef.x, valuef.y};
  }
  operator vec3f() const {
    check_type(glwidgets_param_type::value3f);
    return {valuef.x, valuef.y, valuef.z};
  }
  operator vec4f() const {
    check_type(glwidgets_param_type::value4f);
    return {valuef.x, valuef.y, valuef.z, valuef.w};
  }
  operator int() const {
    check_type(glwidgets_param_type::value1i);
    return valuei.x;
  }
  operator vec2i() const {
    check_type(glwidgets_param_type::value2i);
    return {valuei.x, valuei.y};
  }
  operator vec3i() const {
    check_type(glwidgets_param_type::value3i);
    return {valuei.x, valuei.y, valuei.z};
  }
  operator vec4i() const {
    check_type(glwidgets_param_type::value4i);
    return {valuei.x, valuei.y, valuei.z, valuei.w};
  }
  operator bool() const {
    check_type(glwidgets_param_type::value1b);
    return valueb;
  }
  operator string() const {
    check_type(glwidgets_param_type::value1s);
    return values;
  }
  template <typename T, typename = std::enable_if_t<std::is_enum_v<T>>>
  operator T() const {
    check_type(glwidgets_param_type::value1i);
    return (T)valuei.x;
  }

  // type checking
  void check_type(glwidgets_param_type type) const {
    if (type != this->type) throw std::invalid_argument{"bad gui type"};
  }

  // value
  glwidgets_param_type type   = glwidgets_param_type::value1f;
  vec4f                valuef = {0, 0, 0, 0};
  vec4i                valuei = {0, 0, 0, 0};
  bool                 valueb = false;
  string               values = "";

  // display properties
  vec2f          minmaxf  = {0, 0};
  vec2i          minmaxi  = {0, 0};
  bool           color    = false;
  vector<string> labels   = {};
  bool           readonly = false;
};

struct glwidgets_params {
  using container      = vector<pair<string, glwidgets_param>>;
  using iterator       = container::iterator;
  using const_iterator = container::const_iterator;

  glwidgets_params() {}

  bool   empty() const { return items.empty(); }
  size_t size() const { return items.size(); }

  glwidgets_param& operator[](const string& key) {
    auto item = find(key);
    if (item == end()) return items.emplace_back(key, glwidgets_param{}).second;
    return item->second;
  }
  const glwidgets_param& operator[](const string& key) const { return at(key); }

  glwidgets_param& at(const string& key) {
    auto item = find(key);
    if (item == end()) throw std::out_of_range{"key not found " + key};
    return item->second;
  }
  const glwidgets_param& at(const string& key) const {
    auto item = find(key);
    if (item == end()) throw std::out_of_range{"key not found " + key};
    return item->second;
  }

  iterator find(const string& key) {
    for (auto iterator = items.begin(); iterator != items.end(); ++iterator) {
      if (iterator->first == key) return iterator;
    }
    return items.end();
  }
  const_iterator find(const string& key) const {
    for (auto iterator = items.begin(); iterator != items.end(); ++iterator) {
      if (iterator->first == key) return iterator;
    }
    return items.end();
  }

  iterator       begin() { return items.begin(); }
  iterator       end() { return items.end(); }
  const_iterator begin() const { return items.begin(); }
  const_iterator end() const { return items.end(); }

 private:
  vector<pair<string, glwidgets_param>> items;
};

// draw param
bool draw_glparam(const string& name, glwidgets_param& param) {
  auto copy = param;
  switch (param.type) {
    case glwidgets_param_type::value1f:
      if (param.minmaxf.x == param.minmaxf.y) {
        return draw_gldragger(name.c_str(), param.readonly
                                                ? (float&)copy.valuef
                                                : (float&)param.valuef) &&
               !param.readonly;
      } else {
        return draw_glslider(name.c_str(),
                   param.readonly ? (float&)copy.valuef : (float&)param.valuef,
                   param.minmaxf.x, param.minmaxf.y) &&
               !param.readonly;
      }
      break;
    case glwidgets_param_type::value2f:
      if (param.minmaxf.x == param.minmaxf.y) {
        return draw_gldragger(name.c_str(), param.readonly
                                                ? (vec2f&)copy.valuef
                                                : (vec2f&)param.valuef) &&
               !param.readonly;
      } else {
        return draw_glslider(name.c_str(),
                   param.readonly ? (vec2f&)copy.valuef : (vec2f&)param.valuef,
                   param.minmaxf.x, param.minmaxf.y) &&
               !param.readonly;
      }
      break;
    case glwidgets_param_type::value3f:
      if (param.color) {
        return draw_glcoloredit(name.c_str(), param.readonly
                                                  ? (vec3f&)copy.valuef
                                                  : (vec3f&)param.valuef) &&
               !param.readonly;
      } else if (param.minmaxf.x == param.minmaxf.y) {
        return draw_gldragger(name.c_str(), param.readonly
                                                ? (vec3f&)copy.valuef
                                                : (vec3f&)param.valuef) &&
               !param.readonly;
      } else {
        return draw_glslider(name.c_str(),
                   param.readonly ? copy.valuef : param.valuef, param.minmaxf.x,
                   param.minmaxf.y) &&
               !param.readonly;
      }
      break;
    case glwidgets_param_type::value4f:
      if (param.color) {
        return draw_glcoloredit(name.c_str(), param.readonly
                                                  ? (vec4f&)copy.valuef
                                                  : (vec4f&)param.valuef) &&
               !param.readonly;
      } else if (param.minmaxf.x == param.minmaxf.y) {
        return draw_gldragger(name.c_str(), param.readonly
                                                ? (vec4f&)copy.valuef
                                                : (vec4f&)param.valuef) &&
               !param.readonly;
      } else {
        return draw_glslider(name.c_str(),
                   param.readonly ? (vec4f&)copy.valuef : (vec4f&)param.valuef,
                   param.minmaxf.x, param.minmaxf.y) &&
               !param.readonly;
      }
      break;
    case glwidgets_param_type::value1i:
      if (!param.labels.empty()) {
        return draw_glcombobox(name.c_str(),
                   param.readonly ? (int&)copy.valuei : (int&)param.valuei,
                   param.labels) &&
               !param.readonly;
      } else if (param.minmaxi.x == param.minmaxi.y) {
        return draw_gldragger(name.c_str(),
                   param.readonly ? (int&)copy.valuei : (int&)param.valuei) &&
               !param.readonly;
      } else {
        return draw_glslider(name.c_str(),
                   param.readonly ? (int&)copy.valuei : (int&)param.valuei,
                   param.minmaxi.x, param.minmaxi.y) &&
               !param.readonly;
      }
      break;
    case glwidgets_param_type::value2i:
      if (param.minmaxi.x == param.minmaxi.y) {
        return draw_gldragger(name.c_str(), param.readonly
                                                ? (vec2i&)copy.valuei
                                                : (vec2i&)param.valuei) &&
               !param.readonly;
      } else {
        return draw_glslider(name.c_str(),
                   param.readonly ? (vec2i&)copy.valuei : (vec2i&)param.valuei,
                   param.minmaxi.x, param.minmaxi.y) &&
               !param.readonly;
      }
      break;
    case glwidgets_param_type::value3i:
      if (param.minmaxi.x == param.minmaxi.y) {
        return draw_gldragger(name.c_str(), param.readonly
                                                ? (vec3i&)copy.valuei
                                                : (vec3i&)param.valuei) &&
               !param.readonly;
      } else {
        return draw_glslider(name.c_str(),
                   param.readonly ? (vec3i&)copy.valuei : (vec3i&)param.valuei,
                   param.minmaxi.x, param.minmaxi.y) &&
               !param.readonly;
      }
      break;
    case glwidgets_param_type::value4i:
      if (param.minmaxi.x == param.minmaxi.y) {
        return draw_gldragger(name.c_str(), param.readonly
                                                ? (vec4i&)copy.valuei
                                                : (vec4i&)param.valuei) &&
               !param.readonly;
      } else {
        return draw_glslider(name.c_str(),
                   param.readonly ? (vec4i&)copy.valuei : (vec4i&)param.valuei,
                   param.minmaxi.x, param.minmaxi.y) &&
               !param.readonly;
      }
      break;
    case glwidgets_param_type::value1s:
      if (!param.labels.empty()) {
        return draw_glcombobox(name.c_str(),
                   param.readonly ? copy.values : param.values, param.labels) &&
               !param.readonly;
      } else {
        return draw_gltextinput(
                   name.c_str(), param.readonly ? copy.values : param.values) &&
               !param.readonly;
      }
      break;
    case glwidgets_param_type::value1b:
      if (!param.labels.empty()) {
        // maybe we should implement something different here
        return draw_glcheckbox(
                   name.c_str(), param.readonly ? copy.valueb : param.valueb) &&
               !param.readonly;
      } else {
        return draw_glcheckbox(
                   name.c_str(), param.readonly ? copy.valueb : param.valueb) &&
               !param.readonly;
      }
      break;
    default: return false;
  }
}

// draw params
bool draw_glparams(const string& name, glwidgets_params& params) {
  auto edited = false;
  if (begin_glheader(name.c_str())) {
    for (auto& [name, param] : params) {
      auto pedited = draw_glparam(name, param);
      edited       = edited || pedited;
    }
    end_glheader();
  }
  return edited;
}

}  // namespace yocto
