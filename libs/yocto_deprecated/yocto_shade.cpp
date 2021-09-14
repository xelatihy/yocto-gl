//
// Utilities for real-time rendering of a scene.
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

#include "yocto_shade.h"

#include <array>
#include <memory>
#include <stdexcept>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::make_unique;
using namespace std::string_literals;

}  // namespace yocto

namespace yocto {

#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverlength-strings"
#endif

#ifndef _WIN32
#pragma GCC diagnostic pop
#endif

bool has_normals(const shade_shape& shape) {
  if (shape.vertex_buffers.size() <= 1) return false;
  return is_initialized(shape.vertex_buffers[1]);
}
const ogl_arraybuffer& get_positions(const shade_shape& shape) {
  if (shape.vertex_buffers.size() <= 0)
    throw std::out_of_range{"no vertex buffer"};
  return shape.vertex_buffers[0];
}
const ogl_arraybuffer& get_normals(const shade_shape& shape) {
  if (shape.vertex_buffers.size() <= 1)
    throw std::out_of_range{"no vertex buffer"};
  return shape.vertex_buffers[1];
}
const ogl_arraybuffer& get_texcoords(const shade_shape& shape) {
  if (shape.vertex_buffers.size() <= 2)
    throw std::out_of_range{"no vertex buffer"};
  return shape.vertex_buffers[2];
}
const ogl_arraybuffer& get_colors(const shade_shape& shape) {
  if (shape.vertex_buffers.size() <= 3)
    throw std::out_of_range{"no vertex buffer"};
  return shape.vertex_buffers[3];
}
const ogl_arraybuffer& get_tangents(const shade_shape& shape) {
  if (shape.vertex_buffers.size() <= 4)
    throw std::out_of_range{"no vertex buffer"};
  return shape.vertex_buffers[4];
}

void set_positions(shade_shape& shape, const vector<vec3f>& positions) {
  if (positions.empty()) {
    set_vertex_buffer(shape, vec3f{0, 0, 0}, 0);
  } else {
    set_vertex_buffer(shape, positions, 0);
  }
}
void set_normals(shade_shape& shape, const vector<vec3f>& normals) {
  if (normals.empty()) {
    set_vertex_buffer(shape, vec3f{0, 0, 1}, 1);
  } else {
    set_vertex_buffer(shape, normals, 1);
  }
}
void set_texcoords(shade_shape& shape, const vector<vec2f>& texcoords) {
  if (texcoords.empty()) {
    set_vertex_buffer(shape, vec2f{0, 0}, 2);
  } else {
    set_vertex_buffer(shape, texcoords, 2);
  }
}
void set_colors(shade_shape& shape, const vector<vec4f>& colors) {
  if (colors.empty()) {
    set_vertex_buffer(shape, vec4f{1, 1, 1, 1}, 3);
  } else {
    set_vertex_buffer(shape, colors, 3);
  }
}
void set_tangents(shade_shape& shape, const vector<vec4f>& tangents) {
  if (tangents.empty()) {
    set_vertex_buffer(shape, vec4f{0, 0, 1, 1}, 4);
  } else {
    set_vertex_buffer(shape, tangents, 4);
  }
}
void set_instances(
    shade_shape& shape, const vector<vec3f>& froms, const vector<vec3f>& tos) {
  if (froms.empty()) {
    set_vertex_buffer(shape, vec3f{0, 0, 0}, 5);
    set_instance_buffer(shape, 5, false);
  } else {
    set_vertex_buffer(shape, froms, 5);
    set_instance_buffer(shape, 5, true);
  }
  if (tos.empty()) {
    set_vertex_buffer(shape, vec3f{0, 0, 0}, 5);
    set_instance_buffer(shape, 6, false);
  } else {
    set_vertex_buffer(shape, tos, 6);
    set_instance_buffer(shape, 6, true);
  }
}

void set_points(shade_shape& shape, const vector<int>& points) {
  set_index_buffer(shape, points);
}
void set_lines(shade_shape& shape, const vector<vec2i>& lines) {
  set_index_buffer(shape, lines);
}
void set_triangles(shade_shape& shape, const vector<vec3i>& triangles) {
  set_index_buffer(shape, triangles);
}
void set_quads(shade_shape& shape, const vector<vec4i>& quads) {
  auto triangles = vector<vec3i>{};
  triangles.reserve(quads.size() * 2);
  for (auto& q : quads) {
    triangles.push_back({q.x, q.y, q.w});
    if (q.z != q.w) triangles.push_back({q.z, q.w, q.y});
  }
  set_index_buffer(shape, triangles);
}

// set point size
void set_point_size(shade_shape& shape, float point_size) {
  set_point_size((ogl_shape&)shape, point_size);
}

shade_scene::~shade_scene() { clear_scene(*this); }

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
    shade_scene& scene, shade_environment& environment);
static void init_envlight(shade_scene& scene, shade_environment& environment);

// Initialize an OpenGL scene
void init_scene(shade_scene& scene, bool instanced_drawing) {
  if (is_initialized(scene.instance_program)) return;
  set_program(scene.instance_program,
      instanced_drawing ? shade_instanced_vertex() : shade_instance_vertex(),
      shade_instance_fragment(), true);
  // set_program(scene.envlight_program, shade_instance_vertex(),
  //     shade_instance_fragment(), true);
  set_program(scene.environment_program, precompute_cubemap_vertex(),
      shade_environment_fragment(), true);
}

// Initialize data for environment lighting
void init_environments(shade_scene& scene, bool precompute_envlight) {
  for (auto& environment : scene.environments) {
    init_environment(scene, environment);
    if (precompute_envlight) init_envlight(scene, environment);
  }
}

// Check if we have an envlight
bool has_envlight(const shade_scene& scene) {
  return !scene.environments.empty() && !scene.envlight_shapes.empty();
}

// Clear an OpenGL scene
void clear_scene(shade_scene& scene) {
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

// add camera
glcamera_handle add_camera(shade_scene& scene) {
  scene.cameras.emplace_back();
  return (int)scene.cameras.size() - 1;
}
void set_frame(shade_camera& camera, const frame3f& frame) {
  camera.frame = frame;
}
void set_lens(shade_camera& camera, float lens, float aspect, float film) {
  camera.lens   = lens;
  camera.aspect = aspect;
  camera.film   = film;
}
void set_nearfar(shade_camera& camera, float near, float far) {
  camera.near = near;
  camera.far  = far;
}

// add texture
gltexture_handle add_texture(shade_scene& scene) {
  scene.textures.emplace_back();
  return (int)scene.textures.size() - 1;
}

// add shape
glshape_handle add_shape(shade_scene& scene) {
  scene.shapes.emplace_back();
  return (int)scene.shapes.size() - 1;
}

// add instance
glinstance_handle add_instance(shade_scene& scene) {
  scene.instances.emplace_back();
  return (int)scene.instances.size() - 1;
}
void set_frame(shade_instance& instance, const frame3f& frame) {
  instance.frame = frame;
}
void set_shape(shade_instance& instance, glshape_handle shape) {
  instance.shape = shape;
}
void set_material(shade_instance& instance, glmaterial_handle material) {
  instance.material = material;
}
void set_hidden(shade_instance& instance, bool hidden) {
  instance.hidden = hidden;
}
void set_highlighted(shade_instance& instance, bool highlighted) {
  instance.highlighted = highlighted;
}

// add material
glmaterial_handle add_material(shade_scene& scene) {
  scene.materials.emplace_back();
  return (int)scene.materials.size() - 1;
}
void set_emission(shade_material& material, const vec3f& emission,
    gltexture_handle emission_tex) {
  material.emission     = emission;
  material.emission_tex = emission_tex;
}
void set_color(
    shade_material& material, const vec3f& color, gltexture_handle color_tex) {
  material.color     = color;
  material.color_tex = color_tex;
}
void set_specular(
    shade_material& material, float specular, gltexture_handle specular_tex) {
  material.specular     = specular;
  material.specular_tex = specular_tex;
}
void set_roughness(
    shade_material& material, float roughness, gltexture_handle roughness_tex) {
  material.roughness     = roughness;
  material.roughness_tex = roughness_tex;
}
void set_opacity(
    shade_material& material, float opacity, gltexture_handle opacity_tex) {
  material.opacity = opacity;
}
void set_metallic(
    shade_material& material, float metallic, gltexture_handle metallic_tex) {
  material.metallic     = metallic;
  material.metallic_tex = metallic_tex;
}
void set_normalmap(shade_material& material, gltexture_handle normal_tex) {
  material.normal_tex = normal_tex;
}
void set_unlit(shade_material& material, bool unlit) { material.unlit = unlit; }

glshape_handle add_shape(shade_scene& scene, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<vec4f>& colors, bool edges) {
  auto  handle = add_shape(scene);
  auto& shape  = scene.shapes[handle];
  if (points.size() != 0) {
    set_points(shape, points);
  } else if (lines.size() != 0) {
    set_lines(shape, lines);
  } else if (triangles.size() != 0) {
    set_triangles(shape, triangles);
  } else if (quads.size() != 0) {
    set_quads(shape, quads);
  }
  set_positions(shape, positions);
  set_normals(shape, normals);
  set_texcoords(shape, texcoords);
  set_colors(shape, colors);
  return handle;
}

// environment properties
glenvironment_handle add_environment(shade_scene& scene) {
  scene.environments.emplace_back();
  return (int)scene.environments.size() - 1;
}
void set_frame(shade_environment& environment, const frame3f& frame) {
  environment.frame = frame;
}
void set_emission(shade_environment& environment, const vec3f& emission,
    gltexture_handle emission_tex) {
  environment.emission     = emission;
  environment.emission_tex = emission_tex;
}

// shortcuts
glcamera_handle add_camera(shade_scene& scene, const frame3f& frame, float lens,
    float aspect, float film, float near, float far) {
  auto  handle = add_camera(scene);
  auto& camera = scene.cameras[handle];
  set_frame(camera, frame);
  set_lens(camera, lens, aspect, film);
  set_nearfar(camera, near, far);
  return handle;
}
glmaterial_handle add_material(shade_scene& scene, const vec3f& emission,
    const vec3f& color, float specular, float metallic, float roughness,
    gltexture_handle emission_tex, gltexture_handle color_tex,
    gltexture_handle specular_tex, gltexture_handle metallic_tex,
    gltexture_handle roughness_tex, gltexture_handle normalmap_tex) {
  auto  handle   = add_material(scene);
  auto& material = scene.materials[handle];
  set_emission(material, emission, emission_tex);
  set_color(material, color, color_tex);
  set_specular(material, specular, specular_tex);
  set_metallic(material, metallic, metallic_tex);
  set_roughness(material, roughness, roughness_tex);
  set_normalmap(material, normalmap_tex);
  return handle;
}

glinstance_handle add_instance(shade_scene& scene, const frame3f& frame,
    glshape_handle shape, glmaterial_handle material, bool hidden,
    bool highlighted) {
  auto  handle   = add_instance(scene);
  auto& instance = scene.instances.back();
  set_frame(instance, frame);
  set_shape(instance, shape);
  set_material(instance, material);
  set_hidden(instance, hidden);
  set_highlighted(instance, highlighted);
  return handle;
}

glenvironment_handle add_environment(shade_scene& scene, const frame3f& frame,
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

void set_params_uniforms(ogl_program& program, const shade_params& params) {
  set_uniform(program, "exposure", params.exposure);
  set_uniform(program, "gamma", params.gamma);
  set_uniform(program, "double_sided", params.double_sided);
}

// Draw a shape
void set_instance_uniforms(const shade_scene& scene, ogl_program& program,
    const frame3f& frame, const shade_shape& shape,
    const shade_material& material, const shade_params& params) {
  auto shape_xform     = frame_to_mat(frame);
  auto shape_inv_xform = transpose(
      frame_to_mat(inverse(frame, params.non_rigid_frames)));
  set_uniform(program, "frame", shape_xform);
  set_uniform(program, "frameit", shape_inv_xform);
  set_uniform(program, "offset", 0.0f);
  set_uniform(program, "faceted", params.faceted || !has_normals(shape));
  //  if (instance.highlighted) {
  //    set_uniform(program, "highlight", vec4f{1, 1, 0, 1});
  //  } else {
  //    set_uniform(program, "highlight", vec4f{0, 0, 0, 0});
  //  }

  auto set_texture = [&scene](ogl_program& program, const char* name,
                         const char* name_on, gltexture_handle texture,
                         int unit) {
    if (texture == glinvalid_handle) {
      auto otexture = ogl_texture{};
      set_uniform(program, name, name_on, otexture, unit);
    } else {
      set_uniform(program, name, name_on, scene.textures[texture], unit);
    }
  };

  set_uniform(program, "unlit", material.unlit);
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
  set_texture(
      program, "specular_tex", "specular_tex_on", material.metallic_tex, 2);
  set_texture(
      program, "roughness_tex", "roughness_tex_on", material.roughness_tex, 3);
  set_texture(
      program, "opacity_tex", "opacity_tex_on", material.opacity_tex, 4);
  set_texture(
      program, "normalmap_tex", "normalmap_tex_on", material.normal_tex, 5);

  assert_ogl_error();

  switch (shape.elements) {
    case ogl_element_type::points: set_uniform(program, "element", 1); break;
    case ogl_element_type::line_strip:
    case ogl_element_type::lines: set_uniform(program, "element", 2); break;
    case ogl_element_type::triangle_strip:
    case ogl_element_type::triangle_fan:
    case ogl_element_type::triangles: set_uniform(program, "element", 3); break;
  }
  assert_ogl_error();
}

static void draw_shape(shade_shape& shape) { draw_shape((ogl_shape&)shape); }

void draw_environments(
    shade_scene& scene, const shade_view& view, const shade_params& params) {
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

void set_lighting_uniforms(ogl_program& program, const shade_scene& scene,
    const shade_view& view, const shade_params& params) {
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
  if (lighting == shade_lighting_type::envlight && !has_envlight(scene))
    lighting = shade_lighting_type::camlight;
  if (lighting == shade_lighting_type::envlight && has_envlight(scene)) {
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
  } else if (lighting == shade_lighting_type::camlight) {
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
  } else if (lighting == shade_lighting_type::eyelight) {
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

void draw_instances(
    shade_scene& scene, const shade_view& view, const shade_params& params) {
  // set program
  auto& program = scene.instance_program;
  bind_program(program);

  // set scene uniforms
  set_view_uniforms(program, view);
  set_params_uniforms(program, params);

  // set lighting uniforms
  set_lighting_uniforms(program, scene, view, params);

  set_ogl_wireframe(params.wireframe);
  for (auto& instance : scene.instances) {
    if (instance.hidden) continue;
    set_instance_uniforms(scene, program, instance.frame,
        scene.shapes.at(instance.shape), scene.materials.at(instance.material),
        params);
    draw_shape(scene.shapes.at(instance.shape));
  }
  unbind_program();
}

static shade_view make_scene_view(const shade_camera& camera,
    const vec4i& viewport, const shade_params& params) {
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

void draw_scene(
    shade_scene& scene, const vec4i& viewport, const shade_params& params) {
  clear_ogl_framebuffer(params.background);
  set_ogl_viewport(viewport);

  auto& camera = scene.cameras.at(params.camera);
  auto  view   = make_scene_view(camera, viewport, params);
  draw_instances(scene, view, params);
  draw_environments(scene, view, params);
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
    shade_scene& scene, shade_environment& environment) {
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

void init_envlight(shade_scene& scene, shade_environment& environment) {
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

static const char* shade_instanced_vertex() {
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

}  // namespace yocto
