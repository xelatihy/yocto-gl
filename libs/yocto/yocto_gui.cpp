//
// Simpler image viewer.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2022 Fabio Pellacini
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

#include "yocto_gui.h"

#ifdef YOCTO_OPENGL

#include <glad/glad.h>

#include <cassert>
#include <cstdlib>
#include <future>
#include <stdexcept>

#include "yocto_cutrace.h"
#include "yocto_geometry.h"

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include <GLFW/glfw3.h>
#include <imgui/backends/imgui_impl_glfw.h>
#include <imgui/backends/imgui_impl_opengl3.h>
#include <imgui/imgui.h>
#include <imgui_internal.h>

#ifdef _WIN32
#undef near
#undef far
#endif

// -----------------------------------------------------------------------------
// PARALLEL HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the two integer indices.
template <typename T, typename Func>
inline void parallel_for(T num1, T num2, Func&& func) {
  auto              futures  = vector<std::future<void>>{};
  auto              nthreads = std::thread::hardware_concurrency();
  std::atomic<T>    next_idx(0);
  std::atomic<bool> has_error(false);
  for (auto thread_id = 0; thread_id < (int)nthreads; thread_id++) {
    futures.emplace_back(std::async(
        std::launch::async, [&func, &next_idx, &has_error, num1, num2]() {
          try {
            while (true) {
              auto j = next_idx.fetch_add(1);
              if (j >= num2) break;
              if (has_error) break;
              for (auto i = (T)0; i < num1; i++) func(i, j);
            }
          } catch (...) {
            has_error = true;
            throw;
          }
        }));
  }
  for (auto& f : futures) f.get();
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE DRAWING
// -----------------------------------------------------------------------------
namespace yocto {

// Opengl texture
struct glscene_texture {
  // texture properties
  int width  = 0;
  int height = 0;

  // opengl state
  uint texture = 0;
};

// Create texture
static void set_texture(
    glscene_texture& gltexture, const texture_data& texture);

// Clean texture
static void clear_texture(glscene_texture& gltexture);

// Opengl shape
struct glscene_shape {
  // Shape properties
  int num_positions = 0;
  int num_normals   = 0;
  int num_texcoords = 0;
  int num_colors    = 0;
  int num_tangents  = 0;
  int num_points    = 0;
  int num_lines     = 0;
  int num_triangles = 0;
  int num_quads     = 0;

  // OpenGl state
  uint  vertexarray = 0;
  uint  positions   = 0;
  uint  normals     = 0;
  uint  texcoords   = 0;
  uint  colors      = 0;
  uint  tangents    = 0;
  uint  points      = 0;
  uint  lines       = 0;
  uint  triangles   = 0;
  uint  quads       = 0;
  float point_size  = 1;
};

// Create shape
static void set_shape(glscene_shape& glshape, const shape_data& shape);

// Clean shape
static void clear_shape(glscene_shape& glshape);

// Opengl scene
struct glscene_state {
  // scene objects
  vector<glscene_shape>   shapes   = {};
  vector<glscene_texture> textures = {};

  // programs
  uint program  = 0;
  uint vertex   = 0;
  uint fragment = 0;
};

// init scene
static void init_glscene(glscene_state& glscene, const scene_data& scene);

// update scene
static void update_glscene(glscene_state& glscene, const scene_data& scene,
    const vector<int>& updated_shapes, const vector<int>& updated_textures);

// Clear an OpenGL scene
static void clear_scene(glscene_state& scene);

// draw scene
static void draw_scene(glscene_state& glscene, const scene_data& scene,
    const vec4i& viewport, const shade_params& params);

}  // namespace yocto

// -----------------------------------------------------------------------------
// VIEW HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

void update_image_params(
    const gui_input& input, const image_data& image, glimage_params& glparams) {
  glparams.window                           = input.window;
  glparams.framebuffer                      = input.framebuffer;
  std::tie(glparams.center, glparams.scale) = camera_imview(glparams.center,
      glparams.scale, {image.width, image.height}, glparams.window,
      glparams.fit);
}

bool uiupdate_image_params(const gui_input& input, glimage_params& glparams) {
  // handle mouse
  if (input.mouse.x && input.modifiers.x && !input.onwidgets) {
    if (input.modifiers.z) {
      glparams.scale *= pow(2.0f, (input.cursor.y - input.last.y) * 0.001f);
      return true;
    } else {
      glparams.center += input.cursor - input.last;
      return true;
    }
  }
  return false;
}

bool uiupdate_camera_params(const gui_input& input, camera_data& camera) {
  if (input.mouse.x && input.modifiers.x && !input.onwidgets) {
    auto dolly  = 0.0f;
    auto pan    = vec2f{0, 0};
    auto rotate = vec2f{0, 0};
    if (input.modifiers.y) {
      pan   = (input.cursor - input.last) * camera.focus / 200.0f;
      pan.x = -pan.x;
    } else if (input.modifiers.z) {
      dolly = (input.cursor.y - input.last.y) / 100.0f;
    } else {
      rotate = (input.cursor - input.last) / 100.0f;
    }
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

bool draw_tonemap_widgets(
    const gui_input& input, float& exposure, bool& filmic) {
  auto edited = 0;
  if (draw_gui_header("tonemap")) {
    edited += draw_gui_slider("exposure", exposure, -5, 5);
    edited += draw_gui_checkbox("filmic", filmic);
    end_gui_header();
  }
  return (bool)edited;
}

bool draw_image_widgets(const gui_input& input, const image_data& image,
    const image_data& display, glimage_params& glparams) {
  if (draw_gui_header("inspect")) {
    draw_gui_slider("zoom", glparams.scale, 0.1f, 10);
    draw_gui_checkbox("fit", glparams.fit);
    draw_gui_coloredit("background", glparams.background);
    auto [i, j] = image_coords(input.cursor, glparams.center, glparams.scale,
        vec2i{image.width, image.height});
    auto ij     = vec2i{i, j};
    draw_gui_dragger("mouse", ij);
    auto image_pixel   = vec4f{0, 0, 0, 0};
    auto display_pixel = vec4f{0, 0, 0, 0};
    if (i >= 0 && i < image.width && j >= 0 && j < image.height) {
      image_pixel   = image.pixels[j * image.width + i];
      display_pixel = image.pixels[j * image.width + i];
    }
    draw_gui_coloredit("image", image_pixel);
    draw_gui_coloredit("display", display_pixel);
    end_gui_header();
  }
  return false;
}

bool draw_image_widgets(
    const gui_input& input, const image_data& image, glimage_params& glparams) {
  if (draw_gui_header("inspect")) {
    draw_gui_slider("zoom", glparams.scale, 0.1f, 10);
    draw_gui_checkbox("fit", glparams.fit);
    draw_gui_coloredit("background", glparams.background);
    auto [i, j] = image_coords(input.cursor, glparams.center, glparams.scale,
        vec2i{image.width, image.height});
    auto ij     = vec2i{i, j};
    draw_gui_dragger("mouse", ij);
    auto image_pixel   = vec4f{0, 0, 0, 0};
    auto display_pixel = vec4f{0, 0, 0, 0};
    if (i >= 0 && i < image.width && j >= 0 && j < image.height) {
      image_pixel   = image.pixels[j * image.width + i];
      display_pixel = tonemap(
          image_pixel, glparams.exposure, glparams.filmic, glparams.srgb);
    }
    draw_gui_coloredit("image", image_pixel);
    draw_gui_coloredit("display", display_pixel);
    end_gui_header();
  }
  return false;
}

// draw trace params
bool draw_trace_widgets(const gui_input& input, int sample,
    trace_params& params, const vector<string>& camera_names) {
  auto edited = 0;
  draw_gui_progressbar("sample", sample, params.samples);
  if (draw_gui_header("render")) {
    edited += draw_gui_combobox("camera", params.camera, camera_names);
    edited += draw_gui_slider("resolution", params.resolution, 180, 4096);
    edited += draw_gui_slider("samples", params.samples, 16, 4096);
    edited += draw_gui_combobox(
        "tracer", (int&)params.sampler, trace_sampler_names);
    edited += draw_gui_combobox(
        "false color", (int&)params.falsecolor, trace_falsecolor_names);
    edited += draw_gui_slider("bounces", params.bounces, 1, 128);
    edited += draw_gui_slider("batch", params.batch, 1, 16);
    edited += draw_gui_slider("clamp", params.clamp, 10, 1000);
    edited += draw_gui_checkbox("envhidden", params.envhidden);
    continue_gui_line();
    edited += draw_gui_checkbox("filter", params.tentfilter);
    edited += draw_gui_slider("pratio", params.pratio, 1, 64);
    edited += draw_gui_checkbox("denoise", params.denoise);
    end_gui_header();
  }
  return (bool)edited;
}

bool draw_scene_widgets(scene_data& scene, scene_selection& selection,
    const function<void()>& before_edit) {
  auto edited = 0;
  if (draw_gui_header("cameras")) {
    draw_gui_combobox("camera", selection.camera, scene.camera_names);
    auto camera = scene.cameras.at(selection.camera);
    edited += draw_gui_checkbox("ortho", camera.orthographic);
    edited += draw_gui_slider("lens", camera.lens, 0.001f, 1);
    edited += draw_gui_slider("aspect", camera.aspect, 0.1f, 5);
    edited += draw_gui_slider("film", camera.film, 0.1f, 0.5f);
    edited += draw_gui_slider("focus", camera.focus, 0.001f, 100);
    edited += draw_gui_slider("aperture", camera.aperture, 0, 1);
    if (edited) {
      if (before_edit) before_edit();
      scene.cameras.at(selection.camera) = camera;
    }
    end_gui_header();
  }
  if (draw_gui_header("environments")) {
    draw_gui_combobox(
        "environment", selection.environment, scene.environment_names);
    auto environment = scene.environments.at(selection.environment);
    edited += draw_gui_coloredithdr("emission", environment.emission);
    edited += draw_gui_combobox(
        "emission_tex", environment.emission_tex, scene.texture_names, true);
    if (edited) {
      if (before_edit) before_edit();
      scene.environments.at(selection.environment) = environment;
    }
    end_gui_header();
  }
  if (draw_gui_header("instances")) {
    draw_gui_combobox("instance", selection.instance, scene.instance_names);
    auto instance = scene.instances.at(selection.instance);
    edited += draw_gui_combobox("shape", instance.shape, scene.shape_names);
    edited += draw_gui_combobox(
        "material", instance.material, scene.material_names);
    if (edited) {
      if (before_edit) before_edit();
      scene.instances.at(selection.instance) = instance;
    }
    end_gui_header();
  }
  if (draw_gui_header("materials")) {
    draw_gui_combobox("material", selection.material, scene.material_names);
    auto material = scene.materials.at(selection.material);
    edited += draw_gui_coloredithdr("emission", material.emission);
    edited += draw_gui_combobox(
        "emission_tex", material.emission_tex, scene.texture_names, true);
    edited += draw_gui_coloredithdr("color", material.color);
    edited += draw_gui_combobox(
        "color_tex", material.color_tex, scene.texture_names, true);
    edited += draw_gui_slider("roughness", material.roughness, 0, 1);
    edited += draw_gui_combobox(
        "roughness_tex", material.roughness_tex, scene.texture_names, true);
    edited += draw_gui_slider("metallic", material.metallic, 0, 1);
    edited += draw_gui_slider("ior", material.ior, 0.1f, 5);
    if (edited) {
      if (before_edit) before_edit();
      scene.materials.at(selection.material) = material;
    }
    end_gui_header();
  }
  if (draw_gui_header("shapes")) {
    draw_gui_combobox("shape", selection.shape, scene.shape_names);
    auto& shape = scene.shapes.at(selection.shape);
    draw_gui_label("points", (int)shape.points.size());
    draw_gui_label("lines", (int)shape.lines.size());
    draw_gui_label("triangles", (int)shape.triangles.size());
    draw_gui_label("quads", (int)shape.quads.size());
    draw_gui_label("positions", (int)shape.positions.size());
    draw_gui_label("normals", (int)shape.normals.size());
    draw_gui_label("texcoords", (int)shape.texcoords.size());
    draw_gui_label("colors", (int)shape.colors.size());
    draw_gui_label("radius", (int)shape.radius.size());
    draw_gui_label("tangents", (int)shape.tangents.size());
    end_gui_header();
  }
  if (draw_gui_header("textures")) {
    draw_gui_combobox("texture", selection.texture, scene.texture_names);
    auto& texture = scene.textures.at(selection.texture);
    draw_gui_label("width", texture.width);
    draw_gui_label("height", texture.height);
    draw_gui_label("linear", texture.linear);
    draw_gui_label("byte", !texture.pixelsb.empty());
    end_gui_header();
  }
  if (draw_gui_header("subdivs")) {
    draw_gui_combobox("subdiv", selection.subdiv, scene.subdiv_names);
    auto& subdiv = scene.subdivs.at(selection.subdiv);
    draw_gui_label("quadspos", (int)subdiv.quadspos.size());
    draw_gui_label("quadsnorm", (int)subdiv.quadsnorm.size());
    draw_gui_label("quadstexcoord", (int)subdiv.quadstexcoord.size());
    draw_gui_label("positions", (int)subdiv.positions.size());
    draw_gui_label("normals", (int)subdiv.normals.size());
    draw_gui_label("texcoords", (int)subdiv.texcoords.size());
    end_gui_header();
  }
  return (bool)edited;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE AND TRACE VIEW
// -----------------------------------------------------------------------------
namespace yocto {

// Open a window and show an image
void show_image_gui(
    const string& title, const string& name, const image_data& image) {
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
  auto callbacks = gui_callbacks{};
  callbacks.init = [&](const gui_input& input) {
    init_image(glimage);
    set_image(glimage, display);
  };
  callbacks.clear = [&](const gui_input& input) { clear_image(glimage); };
  callbacks.draw  = [&](const gui_input& input) {
    update_image_params(input, image, glparams);
    draw_image(glimage, glparams);
  };
  callbacks.widgets = [&](const gui_input& input) {
    draw_gui_combobox("name", selected, names);
    if (draw_tonemap_widgets(input, exposure, filmic)) {
      tonemap_image_mt(display, image, exposure, filmic);
      set_image(glimage, display);
    }
    draw_image_widgets(input, image, display, glparams);
  };
  callbacks.uiupdate = [&](const gui_input& input) {
    uiupdate_image_params(input, glparams);
  };

  // run ui
  show_gui_window({1280 + 320, 720}, title, callbacks);
}

// Open a window and show an image
void show_image_gui(const string& title, const vector<string>& names,
    const vector<image_data>& images) {
  // display image
  auto displays  = vector<image_data>(images.size());
  auto exposures = vector<float>(images.size(), 0);
  auto filmics   = vector<bool>(images.size(), false);
  for (auto idx : range((int)images.size())) {
    displays[idx] = make_image(images[idx].width, images[idx].height, false);
    tonemap_image_mt(displays[idx], images[idx], exposures[idx], filmics[idx]);
  }

  // opengl image
  auto glimages  = vector<glimage_state>(images.size());
  auto glparamss = vector<glimage_params>(images.size());

  // selection
  auto selected = 0;

  // callbacks
  auto callbacks = gui_callbacks{};
  callbacks.init = [&](const gui_input& input) {
    for (auto idx : range(images.size())) {
      init_image(glimages[idx]);
      set_image(glimages[idx], displays[idx]);
    }
  };
  callbacks.clear = [&](const gui_input& input) {
    for (auto idx : range(images.size())) {
      clear_image(glimages[idx]);
    }
  };
  callbacks.draw = [&](const gui_input& input) {
    update_image_params(input, displays[selected], glparamss[selected]);
    draw_image(glimages[selected], glparamss[selected]);
  };
  callbacks.widgets = [&](const gui_input& input) {
    draw_gui_combobox("name", selected, names);
    auto filmic = (bool)filmics[selected];  // vector of bool ...
    if (draw_tonemap_widgets(input, exposures[selected], filmic)) {
      filmics[selected] = filmic;
      tonemap_image_mt(displays[selected], images[selected],
          exposures[selected], filmics[selected]);
      set_image(glimages[selected], displays[selected]);
    }
    draw_image_widgets(
        input, images[selected], displays[selected], glparamss[selected]);
  };
  callbacks.uiupdate = [&](const gui_input& input) {
    uiupdate_image_params(input, glparamss[selected]);
  };

  // run ui
  show_gui_window({1280 + 320, 720}, title, callbacks);
}

// Open a window and show an image
void show_colorgrade_gui(
    const string& title, const string& name, const image_data& image) {
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
  auto callbacks = gui_callbacks{};
  callbacks.init = [&](const gui_input& input) {
    init_image(glimage);
    set_image(glimage, display);
  };
  callbacks.clear = [&](const gui_input& input) { clear_image(glimage); };
  callbacks.draw  = [&](const gui_input& input) {
    update_image_params(input, image, glparams);
    draw_image(glimage, glparams);
  };
  callbacks.widgets = [&](const gui_input& input) {
    draw_gui_combobox("name", selected, names);
    if (draw_gui_header("colorgrade")) {
      auto edited = 0;
      edited += draw_gui_slider("exposure", params.exposure, -5, 5);
      edited += draw_gui_coloredit("tint", params.tint);
      edited += draw_gui_slider("lincontrast", params.lincontrast, 0, 1);
      edited += draw_gui_slider("logcontrast", params.logcontrast, 0, 1);
      edited += draw_gui_slider("linsaturation", params.linsaturation, 0, 1);
      edited += draw_gui_checkbox("filmic", params.filmic);
      continue_gui_line();
      edited += draw_gui_checkbox("srgb", params.srgb);
      edited += draw_gui_slider("contrast", params.contrast, 0, 1);
      edited += draw_gui_slider("saturation", params.saturation, 0, 1);
      edited += draw_gui_slider("shadows", params.shadows, 0, 1);
      edited += draw_gui_slider("midtones", params.midtones, 0, 1);
      edited += draw_gui_slider("highlights", params.highlights, 0, 1);
      edited += draw_gui_coloredit("shadows color", params.shadows_color);
      edited += draw_gui_coloredit("midtones color", params.midtones_color);
      edited += draw_gui_coloredit("highlights color", params.highlights_color);
      end_gui_header();
      if (edited) {
        colorgrade_image_mt(display, image, params);
        set_image(glimage, display);
      }
    }
    draw_image_widgets(input, image, display, glparams);
  };
  callbacks.uiupdate = [&glparams](const gui_input& input) {
    uiupdate_image_params(input, glparams);
  };

  // run ui
  show_gui_window({1280 + 320, 720}, title, callbacks);
}

// Open a window and show an scene via path tracing
void show_trace_gui(const string& title, const string& name, scene_data& scene,
    const trace_params& params_, bool print, bool edit) {
  // copy params and camera
  auto params = params_;

  // rendering context
  auto context = make_trace_context(params);

  // build bvh
  auto bvh = make_trace_bvh(scene, params);

  // init renderer
  auto lights = make_trace_lights(scene, params);

  // fix renderer type if no lights
  if (lights.lights.empty() && is_sampler_lit(params)) {
    params.sampler = trace_sampler_type::eyelight;
  }

  // init state
  auto state = make_trace_state(scene, params);
  auto image = make_image(state.width, state.height, true);

  // opengl image
  auto glimage     = glimage_state{};
  auto glparams    = glimage_params{};
  glparams.tonemap = true;

  // top level combo
  auto names    = vector<string>{name};
  auto selected = 0;

  // camera names
  auto camera_names = scene.camera_names;
  if (camera_names.empty()) {
    for (auto idx : range(scene.cameras.size())) {
      camera_names.push_back("camera" + std::to_string(idx + 1));
    }
  }

  // render previews
  auto render_preview = [&]() -> bool {
    // preview
    auto pparams = params;
    pparams.resolution /= params.pratio;
    pparams.samples = 1;
    auto pstate     = make_trace_state(scene, pparams);
    trace_samples(pstate, scene, bvh, lights, pparams);
    auto preview = get_image(pstate);
    for (auto idx = 0; idx < state.width * state.height; idx++) {
      auto i = idx % image.width, j = idx / image.width;
      auto pi           = clamp(i / params.pratio, 0, preview.width - 1),
           pj           = clamp(j / params.pratio, 0, preview.height - 1);
      image.pixels[idx] = preview.pixels[pj * preview.width + pi];
    }
    return true;
  };

  // reset rendering
  auto render_reset = [&]() {
    // make sure we can start
    trace_cancel(context);
    state = make_trace_state(scene, params);
    if (image.width != state.width || image.height != state.height)
      image = make_image(state.width, state.height, true);
  };

  // start rendering batch
  auto render_start = [&]() {
    trace_cancel(context);
    trace_start(context, state, scene, bvh, lights, params);
  };

  // cancel
  auto render_cancel = [&]() { trace_cancel(context); };

  // check if batch is done and update image
  auto render_done = [&]() {
    if (context.done) {
      get_image(image, state);
      return true;
    } else {
      return false;
    }
  };

  // prepare selection
  auto selection = scene_selection{};

  // callbacks
  auto callbacks = gui_callbacks{};
  callbacks.init = [&](const gui_input& input) {
    init_image(glimage);
    render_reset();
    render_preview();
    set_image(glimage, image);
    render_start();
  };
  callbacks.clear = [&](const gui_input& input) { clear_image(glimage); };
  callbacks.draw  = [&](const gui_input& input) {
    // update image
    if (render_done()) {
      set_image(glimage, image);
      render_start();
    }
    update_image_params(input, image, glparams);
    draw_image(glimage, glparams);
  };
  callbacks.widgets = [&](const gui_input& input) {
    auto edited = 0;
    draw_gui_combobox("name", selected, names);
    auto current = (int)state.samples;  // TODO maybe an atomic?
    draw_gui_progressbar("sample", current, params.samples);
    if (draw_gui_header("render")) {
      auto edited  = 0;
      auto tparams = params;
      edited += draw_gui_combobox("camera", tparams.camera, camera_names);
      edited += draw_gui_slider("resolution", tparams.resolution, 180, 4096);
      edited += draw_gui_slider("samples", tparams.samples, 16, 4096);
      edited += draw_gui_combobox(
          "tracer", (int&)tparams.sampler, trace_sampler_names);
      edited += draw_gui_combobox(
          "false color", (int&)tparams.falsecolor, trace_falsecolor_names);
      edited += draw_gui_slider("bounces", tparams.bounces, 1, 128);
      edited += draw_gui_slider("batch", tparams.batch, 1, 16);
      edited += draw_gui_slider("clamp", tparams.clamp, 10, 1000);
      edited += draw_gui_checkbox("envhidden", tparams.envhidden);
      continue_gui_line();
      edited += draw_gui_checkbox("filter", tparams.tentfilter);
      edited += draw_gui_slider("pratio", tparams.pratio, 1, 64);
      edited += draw_gui_checkbox("denoise", tparams.denoise);
      end_gui_header();
      if (edited) {
        render_cancel();
        params = tparams;
        render_reset();
        if (render_preview()) set_image(glimage, image);
        render_start();
      }
    }
    if (draw_gui_header("tonemap")) {
      edited += draw_gui_slider("exposure", glparams.exposure, -5, 5);
      edited += draw_gui_checkbox("filmic", glparams.filmic);
      end_gui_header();
      if (edited) set_image(glimage, image);
    }
    draw_image_widgets(input, image, glparams);
    if (edit) {
      if (draw_scene_widgets(scene, selection, [&]() { render_cancel(); })) {
        render_reset();
        if (render_preview()) set_image(glimage, image);
        render_start();
      }
    }
  };
  callbacks.uiupdate = [&](const gui_input& input) {
    auto camera = scene.cameras[params.camera];
    if (uiupdate_camera_params(input, camera)) {
      render_cancel();
      scene.cameras[params.camera] = camera;
      render_reset();
      if (render_preview()) set_image(glimage, image);
      render_start();
    }
  };

  // run ui
  show_gui_window({1280 + 320, 720}, title, callbacks);

  // done
  render_cancel();
}

#if YOCTO_CUDA

// Open a window and show an scene via path tracing
void show_cutrace_gui(const string& title, const string& name,
    scene_data& scene, const trace_params& params_, bool print, bool edit) {
  // copy params and camera
  auto params = params_;

  // initialize context
  auto context = make_cutrace_context(params);

  // upload scene to the gpu
  auto cuscene = make_cutrace_scene(context, scene, params);

  // build bvh
  auto bvh = make_cutrace_bvh(context, cuscene, params);

  // init lights
  auto lights = make_cutrace_lights(context, scene, params);

  // fix renderer type if no lights
  // if (lights.lights.empty() && is_sampler_lit(params)) {
  //   params.sampler = trace_sampler_type::eyelight;
  // }

  // state
  auto state = make_cutrace_state(context, scene, params);

  // preview state
  auto pparams = params;
  pparams.resolution /= params.pratio;
  pparams.samples = 1;
  auto pstate     = make_cutrace_state(context, scene, pparams);

  // init state
  auto image = make_image(state.width, state.height, true);

  // opengl image
  auto glimage     = glimage_state{};
  auto glparams    = glimage_params{};
  glparams.tonemap = true;

  // top level combo
  auto names    = vector<string>{name};
  auto selected = 0;

  // camera names
  auto camera_names = scene.camera_names;
  if (camera_names.empty()) {
    for (auto idx : range(scene.cameras.size())) {
      camera_names.push_back("camera" + std::to_string(idx + 1));
    }
  }

  // render preview
  auto render_preview = [&]() -> bool {
    auto pparams = params;
    pparams.resolution /= params.pratio;
    pparams.samples = 1;
    reset_cutrace_state(context, pstate, scene, pparams);
    trace_start(context, pstate, cuscene, bvh, lights, scene, pparams);
    trace_samples(context, pstate, cuscene, bvh, lights, scene, pparams);
    auto preview = get_rendered_image(pstate);
    for (auto idx = 0; idx < state.width * state.height; idx++) {
      auto i = idx % image.width, j = idx / image.width;
      auto pi           = clamp(i / params.pratio, 0, preview.width - 1),
           pj           = clamp(j / params.pratio, 0, preview.height - 1);
      image.pixels[idx] = preview.pixels[pj * preview.width + pi];
    }
    return true;
  };

  // reset renderer
  auto render_reset = [&]() {
    reset_cutrace_state(context, state, scene, params);
    if (image.width != state.width || image.height != state.height)
      image = make_image(state.width, state.height, true);
  };

  // render samples synchronously
  auto render_samples = [&]() {
    if (state.samples >= params.samples) return false;
    if (state.samples == 0) {
      trace_start(context, state, cuscene, bvh, lights, scene, params);
    }
    trace_samples(context, state, cuscene, bvh, lights, scene, params);
    get_image(image, state);
    return true;
  };

  // prepare selection
  auto selection = scene_selection{};

  // callbacks
  auto callbacks = gui_callbacks{};
  callbacks.init = [&](const gui_input& input) {
    init_image(glimage);
    render_reset();
    if (render_preview()) set_image(glimage, image);
  };
  callbacks.clear = [&](const gui_input& input) { clear_image(glimage); };
  callbacks.draw  = [&](const gui_input& input) {
    update_image_params(input, image, glparams);
    draw_image(glimage, glparams);
    if (render_samples()) set_image(glimage, image);
  };
  callbacks.widgets = [&](const gui_input& input) {
    auto edited = 0;
    draw_gui_combobox("name", selected, names);
    auto current = state.samples;
    draw_gui_progressbar("sample", current, params.samples);
    if (draw_gui_header("render")) {
      auto edited  = 0;
      auto tparams = params;
      edited += draw_gui_combobox("camera", tparams.camera, camera_names);
      edited += draw_gui_slider("resolution", tparams.resolution, 180, 4096);
      edited += draw_gui_slider("samples", tparams.samples, 16, 4096);
      edited += draw_gui_combobox(
          "tracer", (int&)tparams.sampler, trace_sampler_names);
      edited += draw_gui_combobox(
          "false color", (int&)tparams.falsecolor, trace_falsecolor_names);
      edited += draw_gui_slider("bounces", tparams.bounces, 1, 128);
      edited += draw_gui_slider("batch", tparams.batch, 1, 16);
      edited += draw_gui_slider("clamp", tparams.clamp, 10, 1000);
      edited += draw_gui_checkbox("envhidden", tparams.envhidden);
      continue_gui_line();
      edited += draw_gui_checkbox("filter", tparams.tentfilter);
      edited += draw_gui_slider("pratio", tparams.pratio, 1, 64);
      edited += draw_gui_checkbox("denoise", tparams.denoise);
      end_gui_header();
      if (edited) {
        params = tparams;
        render_reset();
        if (render_preview()) set_image(glimage, image);
      }
    }
    if (draw_gui_header("tonemap")) {
      edited += draw_gui_slider("exposure", glparams.exposure, -5, 5);
      edited += draw_gui_checkbox("filmic", glparams.filmic);
      end_gui_header();
    }
    draw_image_widgets(input, image, glparams);
    if (edit) {
      if (draw_scene_widgets(scene, selection, [&]() {})) {
        render_reset();
        if (render_preview()) set_image(glimage, image);
      }
    }
  };
  callbacks.uiupdate = [&](const gui_input& input) {
    auto camera = scene.cameras[params.camera];
    if (uiupdate_camera_params(input, camera)) {
      scene.cameras[params.camera] = camera;
      update_cutrace_cameras(context, cuscene, scene, params);
      render_reset();
      if (render_preview()) set_image(glimage, image);
    }
  };

  // run ui
  show_gui_window({1280 + 320, 720}, title, callbacks);
}

#else

static void exit_nocuda() { throw std::runtime_error{"Cuda not linked"}; }

// Open a window and show an scene via path tracing
void show_cutrace_gui(const string& title, const string& name,
    scene_data& scene, const trace_params& params, bool print, bool edit) {
  exit_nocuda();
}

#endif

void show_shade_gui(const string& title, const string& name, scene_data& scene,
    const shade_params& params_, const glview_callback& widgets_callback,
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
    for (auto idx : range((int)scene.cameras.size())) {
      camera_names.push_back("camera" + std::to_string(idx + 1));
    }
  }

  // gpu updates
  auto updated_shapes   = vector<int>{};
  auto updated_textures = vector<int>{};

  // callbacks
  auto callbacks = gui_callbacks{};
  callbacks.init = [&](const gui_input& input) {
    init_glscene(glscene, scene);
  };
  callbacks.clear = [&](const gui_input& input) { clear_scene(glscene); };
  callbacks.draw  = [&](const gui_input& input) {
    draw_scene(glscene, scene, input.framebuffer, params);
  };
  callbacks.widgets = [&](const gui_input& input) {
    draw_gui_combobox("name", selected, names);
    if (draw_gui_header("shade")) {
      draw_gui_combobox("camera", params.camera, camera_names);
      draw_gui_checkbox("wireframe", params.wireframe);
      continue_gui_line();
      draw_gui_checkbox("faceted", params.faceted);
      continue_gui_line();
      draw_gui_checkbox("double sided", params.double_sided);
      draw_gui_combobox(
          "lighting", (int&)params.lighting, shade_lighting_names);
      draw_gui_slider("exposure", params.exposure, -10, 10);
      draw_gui_slider("gamma", params.gamma, 0.1f, 4);
      draw_gui_slider("near", params.near, 0.01f, 1.0f);
      draw_gui_slider("far", params.far, 1000.0f, 10000.0f);
      draw_gui_coloredit("background", params.background);
      end_gui_header();
    }
    // draw_scene_widgets(scene, selection, {});
    if (widgets_callback) {
      widgets_callback(input, updated_shapes, updated_textures);
      if (!updated_shapes.empty() || !updated_textures.empty()) {
        update_glscene(glscene, scene, updated_shapes, updated_textures);
        updated_shapes.clear();
        updated_textures.clear();
      }
    }
  };
  callbacks.update = [&](const gui_input& input) {
    if (update_callback) {
      update_callback(input, updated_shapes, updated_textures);
      if (!updated_shapes.empty() || !updated_textures.empty()) {
        update_glscene(glscene, scene, updated_shapes, updated_textures);
        updated_shapes.clear();
        updated_textures.clear();
      }
    }
  };
  callbacks.uiupdate = [&](const gui_input& input) {
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
  show_gui_window({1280 + 320, 720}, "yshade", callbacks);
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
    auto error = string{};
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
static void set_program(uint& program_id, uint& vertex_id, uint& fragment_id,
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
uniform vec2 window, image_size;
uniform vec2 image_center;
uniform float image_scale;
void main() {
    vec2 pos = (positions * 0.5) * image_size * image_scale + image_center;
    gl_Position = vec4(2 * pos.x / window.x - 1, 1 - 2 * pos.y / window.y, 0, 1);
    frag_texcoord = positions * 0.5 + 0.5;
}
)";
#if 0
static auto glimage_vertex = R"(
#version 330
in vec2 positions;
out vec2 frag_texcoord;
uniform vec2 window, image_size, border_size;
uniform vec2 image_center;
uniform float image_scale;
void main() {
    vec2 pos = (positions * 0.5) * (image_size + border_size*2) * image_scale + image_center;
    gl_Position = vec4(2 * pos.x / window.x - 1, 1 - 2 * pos.y / window.y, 0.1, 1);
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
uniform vec4 background;
uniform bool tonemap;
uniform float exposure;
uniform bool srgb;
uniform bool filmic;

float rgb_to_srgb(float rgb) {
  return (rgb <= 0.0031308f) ? 12.92f * rgb
                             : (1 + 0.055f) * pow(rgb, 1 / 2.4f) - 0.055f;
}

vec3 tonemap_filmic(vec3 hdr_) {
  // https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
  vec3 hdr = hdr_ * 0.6f;  // brings it back to ACES range
  vec3 ldr = (hdr * hdr * 2.51f + hdr * 0.03f) /
              (hdr * hdr * 2.43f + hdr * 0.59f + 0.14f);
  return max(vec3(0,0,0), ldr);
}

void main() {
  if(tonemap) {
    vec4 color =  texture(txt, frag_texcoord);
    color.xyz = color.xyz * pow(2, exposure);
    if (filmic) color.xyz = tonemap_filmic(color.xyz);
    if (srgb) color.xyz = vec3(rgb_to_srgb(color.x),rgb_to_srgb(color.y),rgb_to_srgb(color.z));
    frag_color = color;
  } else { 
    frag_color = texture(txt, frag_texcoord);
  }
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

void set_image(glimage_state& glimage, const image_data& image) {
  if (!glimage.texture || glimage.width != image.width ||
      glimage.height != image.height) {
    if (!glimage.texture) glGenTextures(1, &glimage.texture);
    glBindTexture(GL_TEXTURE_2D, glimage.texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, image.width, image.height, 0,
        GL_RGBA, GL_FLOAT, image.pixels.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  } else {
    glBindTexture(GL_TEXTURE_2D, glimage.texture);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, image.width, image.height, GL_RGBA,
        GL_FLOAT, image.pixels.data());
  }
  glimage.width  = image.width;
  glimage.height = image.height;
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

  // blend
  glEnable(GL_BLEND);
  glBlendEquationSeparate(GL_FUNC_ADD, GL_FUNC_ADD);
  glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ZERO);

  // bind program and params
  glUseProgram(glimage.program);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, glimage.texture);
  glUniform1i(glGetUniformLocation(glimage.program, "txt"), 0);
  glUniform2f(glGetUniformLocation(glimage.program, "window"),
      (float)params.window.x, (float)params.window.y);
  glUniform2f(glGetUniformLocation(glimage.program, "image_size"),
      (float)glimage.width, (float)glimage.height);
  glUniform2f(glGetUniformLocation(glimage.program, "image_center"),
      params.center.x, params.center.y);
  glUniform1f(
      glGetUniformLocation(glimage.program, "image_scale"), params.scale);
  glUniform4f(glGetUniformLocation(glimage.program, "background"),
      params.background.x, params.background.y, params.background.z,
      params.background.w);
  glUniform1i(
      glGetUniformLocation(glimage.program, "tonemap"), params.tonemap ? 1 : 0);
  glUniform1f(
      glGetUniformLocation(glimage.program, "exposure"), params.exposure);
  glUniform1i(
      glGetUniformLocation(glimage.program, "srgb"), params.srgb ? 1 : 0);
  glUniform1i(
      glGetUniformLocation(glimage.program, "filmic"), params.filmic ? 1 : 0);
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

  // blend
  glDisable(GL_BLEND);
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
static void set_texture(
    glscene_texture& gltexture, const texture_data& texture) {
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
static void clear_texture(glscene_texture& gltexture) {
  if (gltexture.texture) {
    glDeleteTextures(1, &gltexture.texture);
    gltexture.texture = 0;
  }
}

// Create shape
static void set_shape(glscene_shape& glshape, const shape_data& shape) {
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
static void clear_shape(glscene_shape& glshape) {
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
static void init_glscene(glscene_state& glscene, const scene_data& ioscene) {
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
static void update_glscene(glscene_state& glscene, const scene_data& scene,
    const vector<int>& updated_shapes, const vector<int>& updated_textures) {
  for (auto shape_id : updated_shapes) {
    set_shape(glscene.shapes[shape_id], scene.shapes[shape_id]);
  }
  for (auto texture_id : updated_textures) {
    set_texture(glscene.textures[texture_id], scene.textures[texture_id]);
  }
}

// Clear an OpenGL scene
static void clear_scene(glscene_state& glscene) {
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

static void draw_scene(glscene_state& glscene, const scene_data& scene,
    const vec4i& viewport, const shade_params& params) {
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
      normalize(vec3f{0.1f, 0.5f, -1})};
  static auto lights_emission  = vector<vec3f>{vec3f{pif / 2, pif / 2, pif / 2},
      vec3f{pif / 2, pif / 2, pif / 2}, vec3f{pif / 4, pif / 4, pif / 4},
      vec3f{pif / 4, pif / 4, pif / 4}};
  if (params.lighting == shade_lighting::camlight) {
    glUniform1i(glGetUniformLocation(program, "lighting"), 1);
    glUniform3f(glGetUniformLocation(program, "ambient"), 0, 0, 0);
    glUniform1i(glGetUniformLocation(program, "lights_num"),
        (int)lights_direction.size());
    for (auto lid : range((int)lights_direction.size())) {
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
  } else if (params.lighting == shade_lighting::eyelight) {
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
    if (material.type == material_type::matte ||
        material.type == material_type::transparent ||
        material.type == material_type::refractive ||
        material.type == material_type::subsurface ||
        material.type == material_type::volumetric) {
      glUniform1f(glGetUniformLocation(program, "specular"), 0);
    }
    if (material.type == material_type::reflective) {
      glUniform1f(glGetUniformLocation(program, "metallic"), 1);
    }
    glUniform1i(glGetUniformLocation(program, "double_sided"),
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
  string       title         = "";
  gui_callback init          = {};
  gui_callback clear         = {};
  gui_callback draw          = {};
  gui_callback widgets       = {};
  gui_callback update        = {};
  gui_callback uiupdate      = {};
  int          widgets_width = 0;
  bool         widgets_left  = true;
  gui_input    input         = {};
  vec2i        window        = {0, 0};
  vec4f        background    = {0.15f, 0.15f, 0.15f, 1.0f};
};

static void draw_window(glwindow_state& state) {
  glClearColor(state.background.x, state.background.y, state.background.z,
      state.background.w);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  if (state.draw) state.draw(state.input);
  if (state.widgets) {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    auto window = state.window;
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
      state.widgets(state.input);
    }
    ImGui::End();
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
  }
}

// run the user interface with the give callbacks
void show_gui_window(const vec2i& size, const string& title,
    const gui_callbacks& callbacks, int widgets_width, bool widgets_left) {
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
  auto state     = glwindow_state{};
  state.title    = title;
  state.init     = callbacks.init;
  state.clear    = callbacks.clear;
  state.draw     = callbacks.draw;
  state.widgets  = callbacks.widgets;
  state.update   = callbacks.update;
  state.uiupdate = callbacks.uiupdate;

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
    glfwGetWindowSize(window, &state.window.x, &state.window.y);
    draw_window(state);
    glfwSwapBuffers(window);
  });
  glfwSetWindowSizeCallback(
      window, [](GLFWwindow* window, int width, int height) {
        auto& state = *(glwindow_state*)glfwGetWindowUserPointer(window);
        glfwGetWindowSize(window, &state.input.window.x, &state.input.window.y);
        if (state.widgets_width) state.input.window.x -= state.widgets_width;
        glfwGetFramebufferSize(
            window, &state.input.framebuffer.z, &state.input.framebuffer.w);
        state.input.framebuffer.x = 0;
        state.input.framebuffer.y = 0;
        if (state.widgets_width) {
          auto win_size = vec2i{0, 0};
          glfwGetWindowSize(window, &win_size.x, &win_size.y);
          auto offset = (int)(state.widgets_width *
                              (float)state.input.framebuffer.z / win_size.x);
          state.input.framebuffer.z -= offset;
          if (state.widgets_left) state.input.framebuffer.x += offset;
        }
      });

  // init gl extensions
  if (!gladLoadGL())
    throw std::runtime_error{"cannot initialize OpenGL extensions"};

  // widgets
  if (callbacks.widgets) {
    ImGui::CreateContext();
    ImGui::GetIO().IniFilename       = nullptr;
    ImGui::GetStyle().WindowRounding = 0;
    ImGui_ImplGlfw_InitForOpenGL(window, true);
#ifndef __APPLE__
    ImGui_ImplOpenGL3_Init();
#else
    ImGui_ImplOpenGL3_Init("#version 330");
#endif
    // ImGui::StyleColorsDark();
    ImGui::StyleColorsLight();
    state.widgets_width = widgets_width;
    state.widgets_left  = widgets_left;
  }

  // init
  if (state.init) state.init(state.input);

  // run ui
  while (!glfwWindowShouldClose(window)) {
    // update input
    state.input.last = state.input.cursor;
    auto mouse_posx = 0.0, mouse_posy = 0.0;
    glfwGetCursorPos(window, &mouse_posx, &mouse_posy);
    state.input.cursor = vec2f{(float)mouse_posx, (float)mouse_posy};
    if (state.widgets_width && state.widgets_left)
      state.input.cursor.x -= state.widgets_width;
    state.input.mouse = {
        glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS ? 1
                                                                         : 0,
        glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS ? 1
                                                                          : 0,
        glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS ? 1
                                                                           : 0,
    };
    state.input.modifiers = {
        (glfwGetKey(window, GLFW_KEY_LEFT_ALT) == GLFW_PRESS ||
            glfwGetKey(window, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS)
            ? 1
            : 0,
        (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
            glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS)
            ? 1
            : 0,
        (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS ||
            glfwGetKey(window, GLFW_KEY_RIGHT_CONTROL) == GLFW_PRESS)
            ? 1
            : 0};
    glfwGetWindowSize(window, &state.input.window.x, &state.input.window.y);
    if (state.widgets_width) state.input.window.x -= state.widgets_width;
    glfwGetFramebufferSize(
        window, &state.input.framebuffer.z, &state.input.framebuffer.w);
    state.input.framebuffer.x = 0;
    state.input.framebuffer.y = 0;
    if (state.widgets_width) {
      auto win_size = vec2i{0, 0};
      glfwGetWindowSize(window, &win_size.x, &win_size.y);
      auto offset = (int)(state.widgets_width *
                          (float)state.input.framebuffer.z / win_size.x);
      state.input.framebuffer.z -= offset;
      if (state.widgets_left) state.input.framebuffer.x += offset;
    }
    if (state.widgets_width) {
      auto io               = &ImGui::GetIO();
      state.input.onwidgets = io->WantTextInput || io->WantCaptureMouse ||
                              io->WantCaptureKeyboard;
    }

    // update ui
    if (state.uiupdate && !state.input.onwidgets) state.uiupdate(state.input);

    // update
    if (state.update) state.update(state.input);

    // draw
    glfwGetWindowSize(window, &state.window.x, &state.window.y);
    draw_window(state);
    glfwSwapBuffers(window);

    // event hadling
    glfwPollEvents();
  }

  // clear
  if (state.clear) state.clear(state.input);

  // clear
  glfwDestroyWindow(window);
  glfwTerminate();
}  // namespace yocto

}  // namespace yocto

// -----------------------------------------------------------------------------
// OPENGL WIDGETS
// -----------------------------------------------------------------------------
namespace yocto {

bool draw_gui_header(const char* lbl) {
  if (!ImGui::CollapsingHeader(lbl)) return false;
  ImGui::PushID(lbl);
  return true;
}
void end_gui_header() { ImGui::PopID(); }

bool draw_gui_button(const char* lbl, bool enabled) {
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

void draw_gui_label(const char* lbl, const string& label) {
  ImGui::LabelText(lbl, "%s", label.c_str());
}
void draw_gui_label(const char* lbl, int value) {
  ImGui::LabelText(lbl, "%s", std::to_string(value).c_str());
}
void draw_gui_label(const char* lbl, bool value) {
  ImGui::LabelText(lbl, "%s", value ? "true" : "false");
}

void draw_gui_separator() { ImGui::Separator(); }

void continue_gui_line() { ImGui::SameLine(); }

bool draw_gui_textinput(const char* lbl, string& value) {
  auto buffer = array<char, 4096>{};
  auto num    = 0;
  for (auto c : value) buffer[num++] = c;
  buffer[num] = 0;
  auto edited = ImGui::InputText(lbl, buffer.data(), buffer.size());
  if (edited) value = buffer.data();
  return edited;
}

bool draw_gui_slider(const char* lbl, float& value, float min, float max) {
  return ImGui::SliderFloat(lbl, &value, min, max);
}
bool draw_gui_slider(const char* lbl, vec2f& value, float min, float max) {
  return ImGui::SliderFloat2(lbl, &value.x, min, max);
}
bool draw_gui_slider(const char* lbl, vec3f& value, float min, float max) {
  return ImGui::SliderFloat3(lbl, &value.x, min, max);
}
bool draw_gui_slider(const char* lbl, vec4f& value, float min, float max) {
  return ImGui::SliderFloat4(lbl, &value.x, min, max);
}
bool draw_gui_slider(const char* lbl, int& value, int min, int max) {
  return ImGui::SliderInt(lbl, &value, min, max);
}
bool draw_gui_slider(const char* lbl, vec2i& value, int min, int max) {
  return ImGui::SliderInt2(lbl, &value.x, min, max);
}
bool draw_gui_slider(const char* lbl, vec3i& value, int min, int max) {
  return ImGui::SliderInt3(lbl, &value.x, min, max);
}
bool draw_gui_slider(const char* lbl, vec4i& value, int min, int max) {
  return ImGui::SliderInt4(lbl, &value.x, min, max);
}

bool draw_gui_dragger(
    const char* lbl, float& value, float speed, float min, float max) {
  return ImGui::DragFloat(lbl, &value, speed, min, max);
}
bool draw_gui_dragger(
    const char* lbl, vec2f& value, float speed, float min, float max) {
  return ImGui::DragFloat2(lbl, &value.x, speed, min, max);
}
bool draw_gui_dragger(
    const char* lbl, vec3f& value, float speed, float min, float max) {
  return ImGui::DragFloat3(lbl, &value.x, speed, min, max);
}
bool draw_gui_dragger(
    const char* lbl, vec4f& value, float speed, float min, float max) {
  return ImGui::DragFloat4(lbl, &value.x, speed, min, max);
}
bool draw_gui_dragger(
    const char* lbl, int& value, float speed, int min, int max) {
  return ImGui::DragInt(lbl, &value, speed, min, max);
}
bool draw_gui_dragger(
    const char* lbl, vec2i& value, float speed, int min, int max) {
  return ImGui::DragInt2(lbl, &value.x, speed, min, max);
}
bool draw_gui_dragger(
    const char* lbl, vec3i& value, float speed, int min, int max) {
  return ImGui::DragInt3(lbl, &value.x, speed, min, max);
}
bool draw_gui_dragger(
    const char* lbl, vec4i& value, float speed, int min, int max) {
  return ImGui::DragInt4(lbl, &value.x, speed, min, max);
}

bool draw_gui_dragger(const char* lbl, array<float, 2>& value, float speed,
    float min, float max) {
  return ImGui::DragFloat2(lbl, value.data(), speed, min, max);
}
bool draw_gui_dragger(const char* lbl, array<float, 3>& value, float speed,
    float min, float max) {
  return ImGui::DragFloat3(lbl, value.data(), speed, min, max);
}
bool draw_gui_dragger(const char* lbl, array<float, 4>& value, float speed,
    float min, float max) {
  return ImGui::DragFloat4(lbl, value.data(), speed, min, max);
}

bool draw_gui_dragger(
    const char* lbl, array<int, 2>& value, float speed, int min, int max) {
  return ImGui::DragInt2(lbl, value.data(), speed, min, max);
}
bool draw_gui_dragger(
    const char* lbl, array<int, 3>& value, float speed, int min, int max) {
  return ImGui::DragInt3(lbl, value.data(), speed, min, max);
}
bool draw_gui_dragger(
    const char* lbl, array<int, 4>& value, float speed, int min, int max) {
  return ImGui::DragInt4(lbl, value.data(), speed, min, max);
}

bool draw_gui_checkbox(const char* lbl, bool& value) {
  return ImGui::Checkbox(lbl, &value);
}
bool draw_gui_checkbox(const char* lbl, bool& value, bool invert) {
  if (!invert) {
    return draw_gui_checkbox(lbl, value);
  } else {
    auto inverted = !value;
    auto edited   = ImGui::Checkbox(lbl, &inverted);
    if (edited) value = !inverted;
    return edited;
  }
}

bool draw_gui_coloredit(const char* lbl, vec3f& value) {
  auto flags = ImGuiColorEditFlags_Float;
  return ImGui::ColorEdit3(lbl, &value.x, flags);
}
bool draw_gui_coloredit(const char* lbl, vec4f& value) {
  auto flags = ImGuiColorEditFlags_Float;
  return ImGui::ColorEdit4(lbl, &value.x, flags);
}

bool draw_gui_coloredithdr(const char* lbl, vec3f& value) {
  auto color    = value;
  auto exposure = 0.0f;
  auto scale    = max(color);
  if (scale > 1) {
    color /= scale;
    exposure = log2(scale);
  }
  auto edit_exposure = draw_gui_slider(
      (string{lbl} + " [exp]").c_str(), exposure, 0, 10);
  auto edit_color = draw_gui_coloredit((string{lbl} + " [col]").c_str(), color);
  if (edit_exposure || edit_color) {
    value = color * exp2(exposure);
    return true;
  } else {
    return false;
  }
}
bool draw_gui_coloredithdr(const char* lbl, vec4f& value) {
  auto color    = value;
  auto exposure = 0.0f;
  auto scale    = max(xyz(color));
  if (scale > 1) {
    color.x /= scale;
    color.y /= scale;
    color.z /= scale;
    exposure = log2(scale);
  }
  auto edit_exposure = draw_gui_slider(
      (string{lbl} + " [exp]").c_str(), exposure, 0, 10);
  auto edit_color = draw_gui_coloredit((string{lbl} + " [col]").c_str(), color);
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

bool draw_gui_coloredit(const char* lbl, vec4b& value) {
  auto valuef = byte_to_float(value);
  if (ImGui::ColorEdit4(lbl, &valuef.x)) {
    value = float_to_byte(valuef);
    return true;
  } else {
    return false;
  }
}

bool draw_gui_combobox(const char* lbl, int& value,
    const vector<string>& labels, bool include_null) {
  if (!ImGui::BeginCombo(lbl, value >= 0 ? labels.at(value).c_str() : "<none>"))
    return false;
  auto old_val = value;
  if (include_null) {
    ImGui::PushID(100000);
    if (ImGui::Selectable("<none>", value < 0)) value = -1;
    if (value < 0) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  for (auto i : range((int)labels.size())) {
    ImGui::PushID(i);
    if (ImGui::Selectable(labels[i].c_str(), value == i)) value = i;
    if (value == i) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  ImGui::EndCombo();
  return value != old_val;
}

bool draw_gui_combobox(const char* lbl, string& value,
    const vector<string>& labels, bool include_null) {
  if (!ImGui::BeginCombo(lbl, value.c_str())) return false;
  auto old_val = value;
  if (include_null) {
    ImGui::PushID(100000);
    if (ImGui::Selectable("<none>", value.empty())) value = "";
    if (value.empty()) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  for (auto i : range(labels.size())) {
    ImGui::PushID((int)i);
    if (ImGui::Selectable(labels[i].c_str(), value == labels[i]))
      value = labels[i];
    if (value == labels[i]) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  ImGui::EndCombo();
  return value != old_val;
}

bool draw_gui_combobox(const char* lbl, int& idx, int num,
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
  for (auto i : range(num)) {
    ImGui::PushID(i);
    if (ImGui::Selectable(labels(i).c_str(), idx == i)) idx = i;
    if (idx == i) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  ImGui::EndCombo();
  return idx != old_idx;
}

void draw_gui_progressbar(const char* lbl, float fraction) {
  ImGui::PushStyleColor(ImGuiCol_PlotHistogram, ImVec4(0.5, 0.5, 1, 0.25));
  ImGui::ProgressBar(fraction, ImVec2(0.0f, 0.0f));
  ImGui::SameLine(0.0f, ImGui::GetStyle().ItemInnerSpacing.x);
  ImGui::Text(lbl, ImVec2(0.0f, 0.0f));
  ImGui::PopStyleColor(1);
}

void draw_gui_progressbar(const char* lbl, int current, int total) {
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
  ImGui::PlotHistogram((string{lbl} + " x").c_str(),
      (const float*)values.data() + 0, (int)values.size(), 0, nullptr, flt_max,
      flt_max, {0, 0}, sizeof(vec2f));
  ImGui::PlotHistogram((string{lbl} + " y").c_str(),
      (const float*)values.data() + 1, (int)values.size(), 0, nullptr, flt_max,
      flt_max, {0, 0}, sizeof(vec2f));
}
void draw_histogram(const char* lbl, const vector<vec3f>& values) {
  ImGui::PlotHistogram((string{lbl} + " x").c_str(),
      (const float*)values.data() + 0, (int)values.size(), 0, nullptr, flt_max,
      flt_max, {0, 0}, sizeof(vec3f));
  ImGui::PlotHistogram((string{lbl} + " y").c_str(),
      (const float*)values.data() + 1, (int)values.size(), 0, nullptr, flt_max,
      flt_max, {0, 0}, sizeof(vec3f));
  ImGui::PlotHistogram((string{lbl} + " z").c_str(),
      (const float*)values.data() + 2, (int)values.size(), 0, nullptr, flt_max,
      flt_max, {0, 0}, sizeof(vec3f));
}
void draw_histogram(const char* lbl, const vector<vec4f>& values) {
  ImGui::PlotHistogram((string{lbl} + " x").c_str(),
      (const float*)values.data() + 0, (int)values.size(), 0, nullptr, flt_max,
      flt_max, {0, 0}, sizeof(vec4f));
  ImGui::PlotHistogram((string{lbl} + " y").c_str(),
      (const float*)values.data() + 1, (int)values.size(), 0, nullptr, flt_max,
      flt_max, {0, 0}, sizeof(vec4f));
  ImGui::PlotHistogram((string{lbl} + " z").c_str(),
      (const float*)values.data() + 2, (int)values.size(), 0, nullptr, flt_max,
      flt_max, {0, 0}, sizeof(vec4f));
  ImGui::PlotHistogram((string{lbl} + " w").c_str(),
      (const float*)values.data() + 3, (int)values.size(), 0, nullptr, flt_max,
      flt_max, {0, 0}, sizeof(vec4f));
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
bool draw_gui_param(const string& name, glwidgets_param& param) {
  auto copy = param;
  switch (param.type) {
    case glwidgets_param_type::value1f:
      if (param.minmaxf.x == param.minmaxf.y) {
        return draw_gui_dragger(name.c_str(), param.readonly
                                                  ? (float&)copy.valuef
                                                  : (float&)param.valuef) &&
               !param.readonly;
      } else {
        return draw_gui_slider(name.c_str(),
                   param.readonly ? (float&)copy.valuef : (float&)param.valuef,
                   param.minmaxf.x, param.minmaxf.y) &&
               !param.readonly;
      }
      break;
    case glwidgets_param_type::value2f:
      if (param.minmaxf.x == param.minmaxf.y) {
        return draw_gui_dragger(name.c_str(), param.readonly
                                                  ? (vec2f&)copy.valuef
                                                  : (vec2f&)param.valuef) &&
               !param.readonly;
      } else {
        return draw_gui_slider(name.c_str(),
                   param.readonly ? (vec2f&)copy.valuef : (vec2f&)param.valuef,
                   param.minmaxf.x, param.minmaxf.y) &&
               !param.readonly;
      }
      break;
    case glwidgets_param_type::value3f:
      if (param.color) {
        return draw_gui_coloredit(name.c_str(), param.readonly
                                                    ? (vec3f&)copy.valuef
                                                    : (vec3f&)param.valuef) &&
               !param.readonly;
      } else if (param.minmaxf.x == param.minmaxf.y) {
        return draw_gui_dragger(name.c_str(), param.readonly
                                                  ? (vec3f&)copy.valuef
                                                  : (vec3f&)param.valuef) &&
               !param.readonly;
      } else {
        return draw_gui_slider(name.c_str(),
                   param.readonly ? copy.valuef : param.valuef, param.minmaxf.x,
                   param.minmaxf.y) &&
               !param.readonly;
      }
      break;
    case glwidgets_param_type::value4f:
      if (param.color) {
        return draw_gui_coloredit(name.c_str(), param.readonly
                                                    ? (vec4f&)copy.valuef
                                                    : (vec4f&)param.valuef) &&
               !param.readonly;
      } else if (param.minmaxf.x == param.minmaxf.y) {
        return draw_gui_dragger(name.c_str(), param.readonly
                                                  ? (vec4f&)copy.valuef
                                                  : (vec4f&)param.valuef) &&
               !param.readonly;
      } else {
        return draw_gui_slider(name.c_str(),
                   param.readonly ? (vec4f&)copy.valuef : (vec4f&)param.valuef,
                   param.minmaxf.x, param.minmaxf.y) &&
               !param.readonly;
      }
      break;
    case glwidgets_param_type::value1i:
      if (!param.labels.empty()) {
        return draw_gui_combobox(name.c_str(),
                   param.readonly ? (int&)copy.valuei : (int&)param.valuei,
                   param.labels) &&
               !param.readonly;
      } else if (param.minmaxi.x == param.minmaxi.y) {
        return draw_gui_dragger(name.c_str(),
                   param.readonly ? (int&)copy.valuei : (int&)param.valuei) &&
               !param.readonly;
      } else {
        return draw_gui_slider(name.c_str(),
                   param.readonly ? (int&)copy.valuei : (int&)param.valuei,
                   param.minmaxi.x, param.minmaxi.y) &&
               !param.readonly;
      }
      break;
    case glwidgets_param_type::value2i:
      if (param.minmaxi.x == param.minmaxi.y) {
        return draw_gui_dragger(name.c_str(), param.readonly
                                                  ? (vec2i&)copy.valuei
                                                  : (vec2i&)param.valuei) &&
               !param.readonly;
      } else {
        return draw_gui_slider(name.c_str(),
                   param.readonly ? (vec2i&)copy.valuei : (vec2i&)param.valuei,
                   param.minmaxi.x, param.minmaxi.y) &&
               !param.readonly;
      }
      break;
    case glwidgets_param_type::value3i:
      if (param.minmaxi.x == param.minmaxi.y) {
        return draw_gui_dragger(name.c_str(), param.readonly
                                                  ? (vec3i&)copy.valuei
                                                  : (vec3i&)param.valuei) &&
               !param.readonly;
      } else {
        return draw_gui_slider(name.c_str(),
                   param.readonly ? (vec3i&)copy.valuei : (vec3i&)param.valuei,
                   param.minmaxi.x, param.minmaxi.y) &&
               !param.readonly;
      }
      break;
    case glwidgets_param_type::value4i:
      if (param.minmaxi.x == param.minmaxi.y) {
        return draw_gui_dragger(name.c_str(), param.readonly
                                                  ? (vec4i&)copy.valuei
                                                  : (vec4i&)param.valuei) &&
               !param.readonly;
      } else {
        return draw_gui_slider(name.c_str(),
                   param.readonly ? (vec4i&)copy.valuei : (vec4i&)param.valuei,
                   param.minmaxi.x, param.minmaxi.y) &&
               !param.readonly;
      }
      break;
    case glwidgets_param_type::value1s:
      if (!param.labels.empty()) {
        return draw_gui_combobox(name.c_str(),
                   param.readonly ? copy.values : param.values, param.labels) &&
               !param.readonly;
      } else {
        return draw_gui_textinput(
                   name.c_str(), param.readonly ? copy.values : param.values) &&
               !param.readonly;
      }
      break;
    case glwidgets_param_type::value1b:
      if (!param.labels.empty()) {
        // maybe we should implement something different here
        return draw_gui_checkbox(
                   name.c_str(), param.readonly ? copy.valueb : param.valueb) &&
               !param.readonly;
      } else {
        return draw_gui_checkbox(
                   name.c_str(), param.readonly ? copy.valueb : param.valueb) &&
               !param.readonly;
      }
      break;
    default: return false;
  }
}

// draw params
bool draw_gui_params(const string& name, glwidgets_params& params) {
  auto edited = false;
  if (draw_gui_header(name.c_str())) {
    for (auto& [name, param] : params) {
      auto pedited = draw_gui_param(name, param);
      edited       = edited || pedited;
    }
    end_gui_header();
  }
  return edited;
}

}  // namespace yocto

#endif