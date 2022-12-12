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
  vec2i extents = {0, 0};

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

void update_image_params(const gui_input& input, const array2d<vec4f>& image,
    glimage_params& glparams) {
  glparams.window                           = input.window;
  glparams.framebuffer                      = input.framebuffer;
  std::tie(glparams.center, glparams.scale) = camera_imview(glparams.center,
      glparams.scale, (vec2i)image.extents(), glparams.window, glparams.fit);
}

bool uiupdate_image_params(const gui_input& input, glimage_params& glparams) {
  // handle mouse
  if (input.mouse.left && input.modifiers.alt && !input.onwidgets) {
    if (input.modifiers.control) {
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
  if (input.mouse.left && input.modifiers.alt && !input.onwidgets) {
    auto dolly  = 0.0f;
    auto pan    = vec2f{0, 0};
    auto rotate = vec2f{0, 0};
    if (input.modifiers.shift) {
      pan = (input.cursor - input.last) * camera.focus / 200.0f;
      pan = lerp(vec2f{1, 0}, vec2f{0, 1}, pan);  // flip x
    } else if (input.modifiers.control) {
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

bool draw_image_widgets(const gui_input& input, const array2d<vec4f>& image,
    const array2d<vec4f>& display, glimage_params& glparams) {
  if (draw_gui_header("inspect")) {
    draw_gui_slider("zoom", glparams.scale, 0.1f, 10);
    draw_gui_checkbox("fit", glparams.fit);
    draw_gui_coloredit("background", glparams.background);
    auto ij = image_coords(
        input.cursor, glparams.center, glparams.scale, (vec2i)image.extents());
    draw_gui_dragger("mouse", ij);
    auto image_pixel = image[ij], display_pixel = image[ij];
    draw_gui_coloredit("image", image_pixel);
    draw_gui_coloredit("display", display_pixel);
    end_gui_header();
  }
  return false;
}

bool draw_image_widgets(const gui_input& input, const array2d<vec4f>& image,
    glimage_params& glparams) {
  if (draw_gui_header("inspect")) {
    draw_gui_slider("zoom", glparams.scale, 0.1f, 10);
    draw_gui_checkbox("fit", glparams.fit);
    draw_gui_coloredit("background", glparams.background);
    auto ij = image_coords(
        input.cursor, glparams.center, glparams.scale, (vec2i)image.extents());
    draw_gui_dragger("mouse", ij);
    auto image_pixel   = image[ij];
    auto display_pixel = tonemap(
        image_pixel, glparams.exposure, glparams.filmic, glparams.srgb);
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
    draw_gui_label("width",
        (int)max(texture.pixelsf.extent(0), texture.pixelsb.extent(0)));
    draw_gui_label("height",
        (int)max(texture.pixelsf.extent(1), texture.pixelsb.extent(1)));
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
void show_colorgrade_gui(const string& title, const string& name,
    const array2d<vec4f>& image, bool linear) {
  // color grading parameters
  auto params = colorgrade_params{};

  // display image
  auto display = image;
  colorgrade_image_mt(display, image, linear, params);

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
        colorgrade_image_mt(display, image, linear, params);
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
    const trace_params& params_, float scale, bool print, bool edit) {
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
  auto image = array2d<vec4f>{state.image.extents()};

  // opengl image
  auto glimage     = glimage_state{};
  auto glparams    = glimage_params{};
  glparams.tonemap = true;
  if (scale > 0) {
    glparams.scale = scale;
    glparams.fit   = false;
  }

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
    for (auto ij : range(state.image.extents())) {
      auto pij  = min(ij / params.pratio, preview.extents() - 1);
      image[ij] = preview[pij];
    }
    return true;
  };

  // reset rendering
  auto render_reset = [&]() {
    // make sure we can start
    trace_cancel(context);
    state = make_trace_state(scene, params);
    if (image.extents() != state.image.extents())
      image = array2d<vec4f>{state.image.extents()};
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
  auto image = array2d<vec4f>{state.image.extents()};

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
    for (auto ij : range(image.extents())) {
      auto pij  = min(ij / params.pratio, preview.extents() - 1);
      image[ij] = preview[pij];
    }
    return true;
  };

  // reset renderer
  auto render_reset = [&]() {
    reset_cutrace_state(context, state, scene, params);
    if (image.extents() != state.image.extents())
      image = array2d<vec4f>{state.image.extents()};
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
inline void clear_program(uint& program, uint& vertex, uint& fragment) {
  if (program) glDeleteProgram(program);
  if (vertex) glDeleteShader(vertex);
  if (fragment) glDeleteShader(fragment);
  program  = 0;
  vertex   = 0;
  fragment = 0;
}
inline void bind_program(uint program) { glUseProgram(program); }

// Viewport and framebuffer
inline void bind_viewport(const vec4i& framebuffer) {
  auto [x, y, w, h] = framebuffer;
  glViewport(x, y, w, h);
}
inline void clear_framebuffer(const vec4f& background) {
  auto [r, g, b, a] = background;
  glClearColor(r, g, b, a);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);
}
inline void bind_blending(bool enable) {
  if (enable) {
    glEnable(GL_BLEND);
    glBlendEquationSeparate(GL_FUNC_ADD, GL_FUNC_ADD);
    glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ZERO);
  } else {
    glDisable(GL_BLEND);
  }
}
inline void bind_wireframe(bool enabled) {
  if (enabled) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  } else {
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  }
}

// Vertex arrays
inline void set_vertexarrays(uint& vertexarray) {
  if (!vertexarray) glGenVertexArrays(1, &vertexarray);
  glBindVertexArray(vertexarray);
}
inline void clear_vertexarrays(uint& vertexarray) {
  glDeleteVertexArrays(1, &vertexarray);
  vertexarray = 0;
}
inline void bind_vertexarrays(uint vertexarray) {
  glBindVertexArray(vertexarray);
}

// Vertex buffers
template <typename T>
inline void set_vertex(uint& buffer, int& num, const vector<T>& values,
    const T& def, int location) {
  if (values.empty()) {
    if (buffer) glDeleteBuffers(1, &buffer);
    buffer = 0;
    num    = 0;
    glDisableVertexAttribArray(location);
    if constexpr (sizeof(def) == sizeof(float))
      glVertexAttrib1f(location, (float)def);
    if constexpr (sizeof(def) == sizeof(vec2f))
      glVertexAttrib2fv(location, data(def));
    if constexpr (sizeof(def) == sizeof(vec3f))
      glVertexAttrib3fv(location, data(def));
    if constexpr (sizeof(def) == sizeof(vec4f))
      glVertexAttrib4fv(location, data(def));
  } else {
    if (!buffer || (int)values.size() != num) {
      if (buffer) glDeleteBuffers(1, &buffer);
      glGenBuffers(1, &buffer);
      glBindBuffer(GL_ARRAY_BUFFER, buffer);
      glBufferData(GL_ARRAY_BUFFER, values.size() * sizeof(values.front()),
          values.data(), GL_STATIC_DRAW);
      num = (int)values.size();
    } else {
      // we have enough space
      glBindBuffer(GL_ARRAY_BUFFER, buffer);
      glBufferSubData(GL_ARRAY_BUFFER, 0,
          values.size() * sizeof(values.front()), values.data());
    }
    glBindBuffer(GL_ARRAY_BUFFER, buffer);
    glEnableVertexAttribArray(location);
    glVertexAttribPointer(location, sizeof(values.front()) / sizeof(float),
        GL_FLOAT, false, 0, nullptr);
  }
}
template <typename T>
inline void set_vertex(uint& buffer, const vector<T>& data, int location) {
  auto num = 0;
  auto def = T{};
  set_vertex(buffer, num, data, def, location);
}
inline void clear_vertex(uint& buffer) {
  if (buffer) glDeleteBuffers(1, &buffer);
  buffer = 0;
}

// Primitive indices
template <typename T>
inline void set_indices(uint& buffer, int& num, const vector<T>& data) {
  if (data.empty()) {
    if (buffer) glDeleteBuffers(1, &buffer);
    buffer = 0;
    num    = 0;
  } else {
    if (!buffer || (int)data.size() != num) {
      if (buffer) glDeleteBuffers(1, &buffer);
      glGenBuffers(1, &buffer);
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer);
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, data.size() * sizeof(data.front()),
          data.data(), GL_STATIC_DRAW);
      num = (int)data.size();
    } else {
      // we have enough space
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer);
      glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0,
          data.size() * sizeof(data.front()), data.data());
    }
  }
}
template <typename T>
inline void set_indices(uint& buffer, const vector<T>& data) {
  auto num = 0;
  set_indices(buffer, num, data);
}
inline void clear_indices(uint& buffer) {
  if (buffer) glDeleteBuffers(1, &buffer);
  buffer = 0;
}
inline void draw_points(uint buffer, int num) {
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer);
  glDrawElements(GL_POINTS, (GLsizei)num, GL_UNSIGNED_INT, nullptr);
}
inline void draw_points(uint buffer, int num, float size) {
  glPointSize(size);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer);
  glDrawElements(GL_POINTS, (GLsizei)num, GL_UNSIGNED_INT, nullptr);
  glPointSize(1);
}
inline void draw_lines(uint buffer, int num) {
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer);
  glDrawElements(GL_LINES, (GLsizei)num * 2, GL_UNSIGNED_INT, nullptr);
}
inline void draw_triangles(uint buffer, int num) {
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer);
  glDrawElements(GL_TRIANGLES, num * 3, GL_UNSIGNED_INT, nullptr);
}

// Buffers
inline void clear_buffer(uint& buffer) {
  if (buffer) glDeleteBuffers(1, &buffer);
  buffer = 0;
}

// Textures
inline void set_texture(
    uint& texture, vec2i& extents, const array2d<vec4f>& image) {
  auto [width, height] = (vec2i)image.extents();
  if (!texture || extents != (vec2i)image.extents()) {
    if (!texture) glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0, GL_RGBA,
        GL_FLOAT, image.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  } else {
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexSubImage2D(
        GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGBA, GL_FLOAT, image.data());
  }
  extents = (vec2i)image.extents();
}
static void set_texture(uint& texture, vec2i& extents,
    const array2d<vec4f>& imagef, const array2d<vec4b>& imageb, bool mipmap) {
  auto imextents       = (vec2i)max(imagef.extents(), imageb.extents());
  auto [width, height] = imextents;
  if (!texture || extents != imextents) {
    if (!texture) glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    if (!imageb.empty()) {
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA,
          GL_UNSIGNED_BYTE, imageb.data());
    } else {
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA,
          GL_FLOAT, imagef.data());
    }
    if (mipmap) {
      glGenerateMipmap(GL_TEXTURE_2D);
      glTexParameteri(
          GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    } else {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    }
  } else {
    glBindTexture(GL_TEXTURE_2D, texture);
    if (!imageb.empty()) {
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGBA,
          GL_UNSIGNED_BYTE, imageb.data());
    } else {
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGBA, GL_FLOAT,
          imagef.data());
    }
    if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
  }
  extents = imextents;
}
inline void clear_texture(uint& texture) {
  if (texture) glDeleteTextures(1, &texture);
  texture = 0;
}

// Parameter setting
inline void bind_uniform(int loc, float v) { glUniform1f(loc, v); }
inline void bind_uniform(int loc, vec2f v) { glUniform2fv(loc, 1, &v[0]); }
inline void bind_uniform(int loc, vec3f v) { glUniform3fv(loc, 1, &v[0]); }
inline void bind_uniform(int loc, vec4f v) { glUniform4fv(loc, 1, &v[0]); }
inline void bind_uniform(int loc, int v) { glUniform1i(loc, v); }
inline void bind_uniform(int loc, vec2i v) { glUniform2iv(loc, 1, &v[0]); }
inline void bind_uniform(int loc, vec3i v) { glUniform3iv(loc, 1, &v[0]); }
inline void bind_uniform(int loc, vec4i v) { glUniform4iv(loc, 1, &v[0]); }
inline void bind_uniform(int loc, bool v) { glUniform1i(loc, v ? 1 : 0); }
inline void bind_uniform(int loc, const mat4f& v) {
  glUniformMatrix4fv(loc, 1, false, &v[0][0]);
}
template <typename T>
inline void bind_uniform(uint program, const char* name, const T& v) {
  return bind_uniform(glGetUniformLocation(program, name), v);
}
inline void bind_texture(int loc, uint texture, int unit) {
  glActiveTexture(GL_TEXTURE0 + unit);
  glBindTexture(GL_TEXTURE_2D, texture);
  bind_uniform(loc, unit);
}
inline void bind_texture(
    uint program, const char* name, uint texture, int unit) {
  return bind_texture(glGetUniformLocation(program, name), texture, unit);
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
  set_vertexarrays(glimage.vertexarray);

  // buffers
  auto positions = vector<vec3f>{
      {-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0}};
  set_vertex(glimage.positions, positions, 0);
  auto triangles = vector<vec3i>{{0, 1, 3}, {3, 2, 1}};
  set_indices(glimage.triangles, triangles);

  return true;
}

// clear an opengl image
void clear_image(glimage_state& glimage) {
  clear_texture(glimage.texture);
  clear_program(glimage.program, glimage.vertex, glimage.fragment);
  clear_vertexarrays(glimage.vertexarray);
  clear_buffer(glimage.positions);
  clear_buffer(glimage.triangles);
  glimage = {};
}

void set_image(glimage_state& glimage, const array2d<vec4f>& image) {
  set_texture(glimage.texture, glimage.extents, image);
}

// draw image
void draw_image(glimage_state& glimage, const glimage_params& params) {
  // check errors
  assert_glerror();

  // viewport and framebuffer
  bind_viewport(params.framebuffer);
  clear_framebuffer(params.background);

  // blend
  bind_blending(true);

  // bind program and params
  bind_program(glimage.program);
  bind_texture(glimage.program, "txt", glimage.texture, 0);
  bind_uniform(glimage.program, "window", (vec2f)params.window);
  bind_uniform(glimage.program, "image_size", (vec2f)glimage.extents);
  bind_uniform(glimage.program, "image_center", params.center);
  bind_uniform(glimage.program, "image_scale", params.scale);
  bind_uniform(glimage.program, "background", params.background);
  bind_uniform(glimage.program, "tonemap", params.tonemap);
  bind_uniform(glimage.program, "exposure", params.exposure);
  bind_uniform(glimage.program, "srgb", params.srgb);
  bind_uniform(glimage.program, "filmic", params.filmic);
  assert_glerror();

  // draw
  bind_vertexarrays(glimage.vertexarray);
  draw_triangles(glimage.triangles, 2);
  bind_vertexarrays(0);
  assert_glerror();

  // unbind program
  bind_program(0);
  assert_glerror();

  // blend
  bind_blending(false);
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
  set_texture(gltexture.texture, gltexture.extents, texture.pixelsf,
      texture.pixelsb, true);
}

// Clean texture
static void clear_texture(glscene_texture& gltexture) {
  clear_texture(gltexture.texture);
}

// Create shape
static void set_shape(glscene_shape& glshape, const shape_data& shape) {
  set_vertexarrays(glshape.vertexarray);
  bind_vertexarrays(glshape.vertexarray);
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
  bind_vertexarrays(0);
}

// Clean shape
static void clear_shape(glscene_shape& glshape) {
  clear_vertexarrays(glshape.vertexarray);
  clear_vertex(glshape.positions);
  clear_vertex(glshape.normals);
  clear_vertex(glshape.texcoords);
  clear_vertex(glshape.colors);
  clear_vertex(glshape.tangents);
  clear_indices(glshape.points);
  clear_indices(glshape.lines);
  clear_indices(glshape.triangles);
  clear_indices(glshape.quads);
  glshape = {};
  assert_glerror();
}

// init scene
static void init_glscene(glscene_state& glscene, const scene_data& ioscene) {
  set_program(glscene.program, glscene.vertex, glscene.fragment, glscene_vertex,
      glscene_fragment);
  for (auto& iotexture : ioscene.textures) {
    set_texture(glscene.textures.emplace_back(), iotexture);
  }
  for (auto& ioshape : ioscene.shapes) {
    set_shape(glscene.shapes.emplace_back(), ioshape);
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
  clear_program(glscene.program, glscene.vertex, glscene.fragment);
}

[[maybe_unused]] static void draw_shape(glscene_shape& shape) {
  if (shape.vertexarray == 0) return;
  bind_vertexarrays(shape.vertexarray);
  if (shape.points)
    draw_points(shape.points, shape.num_points, shape.point_size);
  if (shape.lines) draw_lines(shape.lines, shape.num_lines);
  if (shape.triangles) draw_triangles(shape.triangles, shape.num_triangles);
  if (shape.quads) draw_triangles(shape.quads, shape.num_quads);
  bind_vertexarrays(0);
  assert_glerror();
}

static void draw_scene(glscene_state& glscene, const scene_data& scene,
    const vec4i& viewport, const shade_params& params) {
  // check errors
  assert_glerror();

  // viewport and framebuffer
  bind_viewport(viewport);
  clear_framebuffer(params.background);

  // set program
  auto& program = glscene.program;
  bind_program(program);

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
  bind_uniform(program, "eye", camera.frame.o);
  bind_uniform(program, "view", view_matrix);
  bind_uniform(program, "projection", projection_matrix);

  // params
  bind_uniform(program, "exposure", params.exposure);
  bind_uniform(program, "gamma", params.gamma);
  bind_uniform(program, "double_sided", params.double_sided);

  static auto lights_direction = vector<vec3f>{normalize(vec3f{1, 1, 1}),
      normalize(vec3f{-1, 1, 1}), normalize(vec3f{-1, -1, 1}),
      normalize(vec3f{0.1f, 0.5f, -1})};
  static auto lights_emission  = vector<vec3f>{vec3f{pif / 2, pif / 2, pif / 2},
       vec3f{pif / 2, pif / 2, pif / 2}, vec3f{pif / 4, pif / 4, pif / 4},
       vec3f{pif / 4, pif / 4, pif / 4}};
  if (params.lighting == shade_lighting::camlight) {
    bind_uniform(program, "lighting", 1);
    bind_uniform(program, "ambient", vec3f{0, 0, 0});
    bind_uniform(program, "lights_num", (int)lights_direction.size());
    for (auto lid : range((int)lights_direction.size())) {
      auto is        = std::to_string(lid);
      auto direction = transform_direction(camera.frame, lights_direction[lid]);
      bind_uniform(
          program, ("lights_direction[" + is + "]").c_str(), direction);
      bind_uniform(program, ("lights_emission[" + is + "]").c_str(),
          lights_emission[lid]);
    }
  } else if (params.lighting == shade_lighting::eyelight) {
    bind_uniform(program, "lighting", 0);
    bind_uniform(program, "lights_num", 0);
  } else {
    throw std::invalid_argument{"unknown lighting type"};
  }

  // helper
  auto bind_scene_texture =
      [](uint program, const char* name, const char* name_on,
          const glscene_state& glscene, int texture_idx, int unit) {
        if (texture_idx >= 0) {
          auto& gltexture = glscene.textures.at(texture_idx);
          bind_texture(program, name, gltexture.texture, unit);
          bind_uniform(program, name_on, 1);
        } else {
          bind_texture(program, name, 0, unit);
          bind_uniform(program, name_on, 0);
        }
      };

  // draw instances
  bind_wireframe(params.wireframe);
  for (auto& instance : scene.instances) {
    auto& glshape  = glscene.shapes.at(instance.shape);
    auto& material = scene.materials.at(instance.material);

    auto shape_xform     = frame_to_mat(instance.frame);
    auto shape_inv_xform = transpose(
        frame_to_mat(inverse(instance.frame, params.non_rigid_frames)));
    bind_uniform(program, "frame", shape_xform);
    bind_uniform(program, "frameit", shape_inv_xform);
    bind_uniform(program, "faceted", params.faceted || glshape.normals == 0);

    bind_uniform(program, "unlit", 0);
    bind_uniform(program, "emission", material.emission);
    bind_uniform(program, "color", material.color);
    bind_uniform(program, "specular", 1.0f);
    bind_uniform(program, "metallic", material.metallic);
    bind_uniform(program, "roughness", material.roughness);
    bind_uniform(program, "opacity", material.opacity);
    if (material.type == material_type::matte ||
        material.type == material_type::transparent ||
        material.type == material_type::refractive ||
        material.type == material_type::subsurface ||
        material.type == material_type::volumetric) {
      bind_uniform(program, "specular", 0.0f);
    }
    if (material.type == material_type::reflective) {
      bind_uniform(program, "metallic", 1.0f);
    }
    bind_uniform(program, "double_sided", params.double_sided);
    bind_scene_texture(program, "emission_tex", "emission_tex_on", glscene,
        material.emission_tex, 0);
    bind_scene_texture(
        program, "color_tex", "color_tex_on", glscene, material.color_tex, 1);
    bind_scene_texture(program, "roughness_tex", "roughness_tex_on", glscene,
        material.roughness_tex, 3);
    bind_scene_texture(program, "normalmap_tex", "normalmap_tex_on", glscene,
        material.normal_tex, 5);
    assert_glerror();

    if (glshape.points) bind_uniform(program, "element", 1);
    if (glshape.lines) bind_uniform(program, "element", 2);
    if (glshape.triangles) bind_uniform(program, "element", 3);
    if (glshape.quads) bind_uniform(program, "element", 3);
    assert_glerror();

    bind_vertexarrays(glshape.vertexarray);
    if (glshape.points) {
      draw_points(glshape.points, glshape.num_points, glshape.point_size);
    }
    if (glshape.lines) {
      draw_lines(glshape.lines, glshape.num_lines);
    }
    if (glshape.triangles) {
      draw_triangles(glshape.triangles, glshape.num_triangles);
    }
    if (glshape.quads) {
      draw_triangles(glshape.quads, glshape.num_quads);
    }

    bind_vertexarrays(0);
    assert_glerror();
  }
  bind_wireframe(false);

  // done
  bind_program(0);
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
  auto [r, g, b, a] = state.background;
  glClearColor(r, g, b, a);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  if (state.draw) state.draw(state.input);
  if (state.widgets) {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    auto window          = state.window;
    auto [width, height] = window;
    if (state.widgets_left) {
      ImGui::SetNextWindowPos({0, 0});
      ImGui::SetNextWindowSize({(float)state.widgets_width, (float)height});
    } else {
      ImGui::SetNextWindowPos({(float)(width - state.widgets_width), 0});
      ImGui::SetNextWindowSize({(float)state.widgets_width, (float)height});
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

// helpers
vec2i get_window_size(GLFWwindow* window) {
  auto width = 0, height = 0;
  glfwGetWindowSize(window, &width, &height);
  return {width, height};
}
vec2i get_framebuffer_size(GLFWwindow* window) {
  auto width = 0, height = 0;
  glfwGetFramebufferSize(window, &width, &height);
  return {width, height};
}
vec4i get_framebuffer_viewport(GLFWwindow* window) {
  auto width = 0, height = 0;
  glfwGetFramebufferSize(window, &width, &height);
  return {0, 0, width, height};
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
  auto [width, height] = size;
  auto window          = glfwCreateWindow(
      width, height, title.c_str(), nullptr, nullptr);
  if (window == nullptr)
    throw std::runtime_error{"cannot initialize windowing system"};
  glfwMakeContextCurrent(window);
  glfwSwapInterval(1);  // Enable vsync

  // set user data
  glfwSetWindowUserPointer(window, &state);

  // set callbacks
  glfwSetWindowRefreshCallback(window, [](GLFWwindow* window) {
    auto& state  = *(glwindow_state*)glfwGetWindowUserPointer(window);
    state.window = get_window_size(window);
    draw_window(state);
    glfwSwapBuffers(window);
  });
  glfwSetWindowSizeCallback(
      window, [](GLFWwindow* window, int width, int height) {
        auto& state        = *(glwindow_state*)glfwGetWindowUserPointer(window);
        state.input.window = get_window_size(window);
        if (state.widgets_width)
          state.input.window -= vec2i{state.widgets_width, 0};
        state.input.framebuffer = get_framebuffer_viewport(window);
        if (state.widgets_width) {
          auto [win_width, _]          = get_window_size(window);
          auto [framebuffer_width, __] = get_framebuffer_size(window);
          auto offset =
              (int)(state.widgets_width * (float)framebuffer_width / win_width);
          state.input.framebuffer -= vec4i{0, 0, offset, 0};
          if (state.widgets_left)
            state.input.framebuffer += vec4i{offset, 0, 0, 0};
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
      state.input.cursor -= vec2i{state.widgets_width, 0};
    state.input.mouse = {
        glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS,
        glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS,
        glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS};
    state.input.modifiers = {
        glfwGetKey(window, GLFW_KEY_LEFT_ALT) == GLFW_PRESS ||
            glfwGetKey(window, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS,
        glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
            glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS,
        glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS ||
            glfwGetKey(window, GLFW_KEY_RIGHT_CONTROL) == GLFW_PRESS};
    state.input.window = get_window_size(window);
    if (state.widgets_width)
      state.input.window -= vec2i{state.widgets_width, 0};
    state.input.framebuffer = get_framebuffer_viewport(window);
    if (state.widgets_width) {
      auto [win_width, _]          = get_window_size(window);
      auto [framebuffer_width, __] = get_framebuffer_size(window);
      auto offset =
          (int)(state.widgets_width * (float)framebuffer_width / win_width);
      state.input.framebuffer -= vec4i{0, 0, offset, 0};
      if (state.widgets_left) state.input.framebuffer += vec4i{offset, 0, 0, 0};
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
    state.window = get_window_size(window);
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
  return ImGui::SliderFloat2(lbl, data(value), min, max);
}
bool draw_gui_slider(const char* lbl, vec3f& value, float min, float max) {
  return ImGui::SliderFloat3(lbl, data(value), min, max);
}
bool draw_gui_slider(const char* lbl, vec4f& value, float min, float max) {
  return ImGui::SliderFloat4(lbl, data(value), min, max);
}
bool draw_gui_slider(const char* lbl, int& value, int min, int max) {
  return ImGui::SliderInt(lbl, &value, min, max);
}
bool draw_gui_slider(const char* lbl, vec2i& value, int min, int max) {
  return ImGui::SliderInt2(lbl, data(value), min, max);
}
bool draw_gui_slider(const char* lbl, vec3i& value, int min, int max) {
  return ImGui::SliderInt3(lbl, data(value), min, max);
}
bool draw_gui_slider(const char* lbl, vec4i& value, int min, int max) {
  return ImGui::SliderInt4(lbl, data(value), min, max);
}

bool draw_gui_dragger(
    const char* lbl, float& value, float speed, float min, float max) {
  return ImGui::DragFloat(lbl, &value, speed, min, max);
}
bool draw_gui_dragger(
    const char* lbl, vec2f& value, float speed, float min, float max) {
  return ImGui::DragFloat2(lbl, data(value), speed, min, max);
}
bool draw_gui_dragger(
    const char* lbl, vec3f& value, float speed, float min, float max) {
  return ImGui::DragFloat3(lbl, data(value), speed, min, max);
}
bool draw_gui_dragger(
    const char* lbl, vec4f& value, float speed, float min, float max) {
  return ImGui::DragFloat4(lbl, data(value), speed, min, max);
}
bool draw_gui_dragger(
    const char* lbl, int& value, float speed, int min, int max) {
  return ImGui::DragInt(lbl, &value, speed, min, max);
}
bool draw_gui_dragger(
    const char* lbl, vec2i& value, float speed, int min, int max) {
  return ImGui::DragInt2(lbl, data(value), speed, min, max);
}
bool draw_gui_dragger(
    const char* lbl, vec3i& value, float speed, int min, int max) {
  return ImGui::DragInt3(lbl, data(value), speed, min, max);
}
bool draw_gui_dragger(
    const char* lbl, vec4i& value, float speed, int min, int max) {
  return ImGui::DragInt4(lbl, data(value), speed, min, max);
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
  return ImGui::ColorEdit3(lbl, data(value), flags);
}
bool draw_gui_coloredit(const char* lbl, vec4f& value) {
  auto flags = ImGuiColorEditFlags_Float;
  return ImGui::ColorEdit4(lbl, data(value), flags);
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
    color    = {xyz(color) / scale, alpha(color)};
    exposure = log2(scale);
  }
  auto edit_exposure = draw_gui_slider(
      (string{lbl} + " [exp]").c_str(), exposure, 0, 10);
  auto edit_color = draw_gui_coloredit((string{lbl} + " [col]").c_str(), color);
  if (edit_exposure || edit_color) {
    value = {xyz(color) * exp2(exposure), alpha(color)};
    return true;
  } else {
    return false;
  }
}

bool draw_gui_coloredit(const char* lbl, vec4b& value) {
  auto valuef = byte_to_float(value);
  if (ImGui::ColorEdit4(lbl, data(valuef))) {
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

#endif