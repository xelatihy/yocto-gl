//
// Utilities for real-time rendering of a scene.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
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

#include "yocto_draw.h"

#include <yocto/yocto_commonio.h>

#include <array>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
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

ogl_arraybuffer* get_positions(gui_shape* shape) {
  if (shape->vertex_buffers.size() <= 0) return nullptr;
  return &shape->vertex_buffers[0];
}
ogl_arraybuffer* get_normals(gui_shape* shape) {
  if (shape->vertex_buffers.size() <= 1) return nullptr;
  return &shape->vertex_buffers[1];
}
ogl_arraybuffer* get_texcoords(gui_shape* shape) {
  if (shape->vertex_buffers.size() <= 2) return nullptr;
  return &shape->vertex_buffers[2];
}
ogl_arraybuffer* get_colors(gui_shape* shape) {
  if (shape->vertex_buffers.size() <= 3) return nullptr;
  return &shape->vertex_buffers[3];
}
ogl_arraybuffer* get_tangents(gui_shape* shape) {
  if (shape->vertex_buffers.size() <= 4) return nullptr;
  return &shape->vertex_buffers[4];
}

void set_positions(gui_shape* shape, const vector<vec3f>& positions) {
  if (positions.empty())
    set_vertex_buffer(shape, vec3f{0, 0, 0}, 0);
  else
    set_vertex_buffer(shape, positions, 0);
}
void set_normals(gui_shape* shape, const vector<vec3f>& normals) {
  if (normals.empty())
    set_vertex_buffer(shape, vec3f{0, 0, 1}, 1);
  else
    set_vertex_buffer(shape, normals, 1);
}
void set_texcoords(gui_shape* shape, const vector<vec2f>& texcoords) {
  if (texcoords.empty())
    set_vertex_buffer(shape, vec2f{0, 0}, 2);
  else
    set_vertex_buffer(shape, texcoords, 2);
}
void set_colors(gui_shape* shape, const vector<vec4f>& colors) {
  if (colors.empty())
    set_vertex_buffer(shape, vec4f{1, 1, 1, 1}, 3);
  else
    set_vertex_buffer(shape, colors, 3);
}
void set_tangents(gui_shape* shape, const vector<vec4f>& tangents) {
  if (tangents.empty())
    set_vertex_buffer(shape, vec4f{0, 0, 1, 1}, 4);
  else
    set_vertex_buffer(shape, tangents, 4);
}

void set_points(gui_shape* shape, const vector<int>& points) {
  set_index_buffer(shape, points);
}
void set_lines(gui_shape* shape, const vector<vec2i>& lines) {
  set_index_buffer(shape, lines);
}
void set_triangles(gui_shape* shape, const vector<vec3i>& triangles) {
  set_index_buffer(shape, triangles);
}
void set_quads(gui_shape* shape, const vector<vec4i>& quads) {
  auto triangles = vector<vec3i>{};
  triangles.reserve(quads.size() * 2);
  for (auto& q : quads) {
    triangles.push_back({q.x, q.y, q.w});
    if (q.z != q.w) triangles.push_back({q.z, q.w, q.y});
  }
  set_index_buffer(shape, triangles);
}

gui_scene::~gui_scene() {
  clear_scene(this);
  for (auto camera : cameras) delete camera;
  for (auto shape : shapes) delete shape;
  for (auto material : materials) delete material;
  for (auto texture : textures) delete texture;
  delete environment_shape;
  delete environment_cubemap;
  delete environment_program;
  delete diffuse_cubemap;
  delete specular_cubemap;
  delete brdf_lut;
  delete eyelight_program;
  delete ibl_program;
}

static const char* bake_brdf_vertex_code();
static const char* bake_brdf_fragment_code();

static const char* bake_cubemap_vertex_code();
static const char* bake_environment_fragment_code();
static const char* bake_irradiance_fragment_code();
static const char* bake_reflections_fragment_code();

static void init_environment(gui_scene* scene,
    const gui_texture* environment_tex, const vec3f& environment_emission);

// Initialize an OpenGL scene
void init_scene(gui_scene* scene, const gui_texture* environment_tex,
    const vec3f& environment_emission) {
  if (is_initialized(scene->eyelight_program)) return;
  auto error = ""s, errorlog = ""s;
  auto vert = draw_instances_vertex_code();
  auto frag = draw_instances_eyelight_fragment_code();
  init_program(scene->eyelight_program, vert, frag, error, errorlog);

  if (environment_tex && environment_emission != vec3f{0, 0, 0}) {
    init_environment(scene, environment_tex, environment_emission);
    init_ibl_data(scene);
  }
}

bool is_initialized(gui_scene* scene) {
  return scene && is_initialized(scene->eyelight_program);
}

// Clear an OpenGL scene
void clear_scene(gui_scene* scene) {
  for (auto texture : scene->textures) clear_texture(texture);
  for (auto shape : scene->shapes) clear_shape(shape);
  clear_shape(scene->environment_shape);
  clear_cubemap(scene->environment_cubemap);
  clear_program(scene->environment_program);
  clear_cubemap(scene->diffuse_cubemap);
  clear_cubemap(scene->specular_cubemap);
  clear_texture(scene->brdf_lut);
  clear_program(scene->eyelight_program);
  clear_program(scene->ibl_program);
}

// add camera
gui_camera* add_camera(gui_scene* scene) {
  return scene->cameras.emplace_back(new gui_camera{});
}
void set_frame(gui_camera* camera, const frame3f& frame) {
  camera->frame = frame;
}
void set_lens(gui_camera* camera, float lens, float aspect, float film) {
  camera->lens   = lens;
  camera->aspect = aspect;
  camera->film   = film;
}
void set_nearfar(gui_camera* camera, float near, float far) {
  camera->near = near;
  camera->far  = far;
}

// add texture
gui_texture* add_texture(gui_scene* scene) {
  return scene->textures.emplace_back(new gui_texture{});
}

// add shape
gui_shape* add_shape(gui_scene* scene) {
  auto shape = new gui_shape{};
  shape->vertex_buffers.resize(5);
  set_shape(shape);
  scene->shapes.push_back(shape);
  return shape;
}

// add instance
gui_instance* add_instance(gui_scene* scene) {
  return scene->instances.emplace_back(new gui_instance{});
}
void set_frame(gui_instance* instance, const frame3f& frame) {
  instance->frame = frame;
}
void set_shape(gui_instance* instance, gui_shape* shape) {
  instance->shape = shape;
}
void set_material(gui_instance* instance, gui_material* material) {
  instance->material = material;
}
void set_hidden(gui_instance* instance, bool hidden) {
  instance->hidden = hidden;
}
void set_highlighted(gui_instance* instance, bool highlighted) {
  instance->highlighted = highlighted;
}

// add material
gui_material* add_material(gui_scene* scene) {
  return scene->materials.emplace_back(new gui_material{});
}
void set_emission(
    gui_material* material, const vec3f& emission, gui_texture* emission_tex) {
  material->emission     = emission;
  material->emission_tex = emission_tex;
}
void set_color(
    gui_material* material, const vec3f& color, gui_texture* color_tex) {
  material->color     = color;
  material->color_tex = color_tex;
}
void set_specular(
    gui_material* material, float specular, gui_texture* specular_tex) {
  material->specular     = specular;
  material->specular_tex = specular_tex;
}
void set_roughness(
    gui_material* material, float roughness, gui_texture* roughness_tex) {
  material->roughness     = roughness;
  material->roughness_tex = roughness_tex;
}
void set_opacity(
    gui_material* material, float opacity, gui_texture* opacity_tex) {
  material->opacity = opacity;
}
void set_metallic(
    gui_material* material, float metallic, gui_texture* metallic_tex) {
  material->metallic     = metallic;
  material->metallic_tex = metallic_tex;
}
void set_normalmap(gui_material* material, gui_texture* normal_tex) {
  material->normal_tex = normal_tex;
}

gui_shape* add_shape(gui_scene* scene, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<vec4f>& colors, bool edges) {
  auto shape = add_shape(scene);
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
  return shape;
}

// shortcuts
gui_camera* add_camera(gui_scene* scene, const frame3f& frame, float lens,
    float aspect, float film, float near, float far) {
  auto camera = add_camera(scene);
  set_frame(camera, frame);
  set_lens(camera, lens, aspect, film);
  set_nearfar(camera, near, far);
  return camera;
}
gui_material* add_material(gui_scene* scene, const vec3f& emission,
    const vec3f& color, float specular, float metallic, float roughness,
    gui_texture* emission_tex, gui_texture* color_tex,
    gui_texture* specular_tex, gui_texture* metallic_tex,
    gui_texture* roughness_tex, gui_texture* normalmap_tex) {
  auto material = add_material(scene);
  set_emission(material, emission, emission_tex);
  set_color(material, color, color_tex);
  set_specular(material, specular, specular_tex);
  set_metallic(material, metallic, metallic_tex);
  set_roughness(material, roughness, roughness_tex);
  set_normalmap(material, normalmap_tex);
  return material;
}

gui_instance* add_instance(gui_scene* scene, const frame3f& frame,
    gui_shape* shape, gui_material* material, bool hidden, bool highlighted) {
  auto instance = add_instance(scene);
  set_frame(instance, frame);
  set_shape(instance, shape);
  set_material(instance, material);
  set_hidden(instance, hidden);
  set_highlighted(instance, highlighted);
  return instance;
}

void set_scene_view_uniforms(ogl_program* program, const gui_scene_view& view) {
  set_uniform(program, "eye", view.camera_frame.o);
  set_uniform(program, "view", view.view_matrix);
  set_uniform(program, "projection", view.projection_matrix);
  set_uniform(program, "exposure", view.params.exposure);
  set_uniform(program, "gamma", view.params.gamma);
  set_uniform(program, "faceted", view.params.faceted);
  set_uniform(program, "double_sided", view.params.double_sided);
}

// Draw a shape
void set_instance_uniforms(ogl_program* program, const frame3f& frame,
    const gui_shape* shape, const gui_material* material,
    gui_shading_type shading, bool double_sided, bool non_rigid_frames) {
  auto shape_xform     = frame_to_mat(frame);
  auto shape_inv_xform = transpose(
      frame_to_mat(inverse(frame, non_rigid_frames)));
  set_uniform(program, "frame", shape_xform);
  set_uniform(program, "frameit", shape_inv_xform);
  set_uniform(program, "offset", 0.0f);
  //  if (instance->highlighted) {
  //    set_uniform(program, "highlight", vec4f{1, 1, 0, 1});
  //  } else {
  //    set_uniform(program, "highlight", vec4f{0, 0, 0, 0});
  //  }

  auto mtype = (int)shading;
  set_uniform(program, "mtype", mtype);
  set_uniform(program, "emission", material->emission);
  set_uniform(program, "diffuse", material->color);
  set_uniform(program, "specular",
      vec3f{material->metallic, material->metallic, material->metallic});
  set_uniform(program, "roughness", material->roughness);
  set_uniform(program, "opacity", material->opacity);
  set_uniform(program, "double_sided", double_sided);
  set_uniform(
      program, "emission_tex", "emission_tex_on", material->emission_tex, 0);
  set_uniform(program, "diffuse_tex", "diffuse_tex_on", material->color_tex, 1);
  set_uniform(
      program, "specular_tex", "specular_tex_on", material->metallic_tex, 2);
  set_uniform(
      program, "roughness_tex", "roughness_tex_on", material->roughness_tex, 3);
  set_uniform(
      program, "opacity_tex", "opacity_tex_on", material->opacity_tex, 4);
  set_uniform(
      program, "mat_norm_tex", "mat_norm_tex_on", material->normal_tex, 5);

  assert_ogl_error();

  auto type = shape->elements;
  if (type == ogl_element_type::points) set_uniform(program, "etype", 1);
  if (type == ogl_element_type::lines) set_uniform(program, "etype", 2);
  if (type == ogl_element_type::triangles) set_uniform(program, "etype", 3);
  assert_ogl_error();
}

void draw_environment(gui_scene* scene, const gui_scene_view& view) {
  auto program = scene->environment_program;
  if (program->program_id == 0) return;

  bind_program(program);

  set_scene_view_uniforms(program, view);
  set_uniform(program, "environment", scene->environment_cubemap, 0);

  draw_shape(scene->environment_shape);

  unbind_program();
}

void set_eyelight_uniforms(ogl_program* program, const gui_scene_view& view) {
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

  auto& lights = camera_lights;
  set_uniform(program, "lamb", vec3f{0, 0, 0});
  set_uniform(program, "lnum", (int)lights.size());
  auto lid = 0;
  for (auto light : lights) {
    auto is = std::to_string(lid);
    if (light->camera) {
      auto position = transform_direction(view.camera_frame, light->position);
      set_uniform(program, ("lpos[" + is + "]").c_str(), position);
    } else {
      set_uniform(program, ("lpos[" + is + "]").c_str(), light->position);
    }
    set_uniform(program, ("lke[" + is + "]").c_str(), light->emission);
    set_uniform(program, ("ltype[" + is + "]").c_str(), 1);
    lid++;
  }
  assert_ogl_error();
}

void set_ibl_uniforms(ogl_program* program, const gui_scene* scene) {
  set_uniform(program, "irradiance_cubemap", scene->diffuse_cubemap, 6);
  set_uniform(program, "reflection_cubemap", scene->specular_cubemap, 7);
  set_uniform(program, "brdf_lut", scene->brdf_lut, 8);
}

void draw_instances(gui_scene* scene, const gui_scene_view& view) {
  auto program = scene->eyelight_program;
  if (is_initialized(scene->environment_cubemap) &&
      view.params.lighting == gui_lighting_type::environment) {
    program = scene->ibl_program;
  }

  bind_program(program);

  // set scene uniforms
  set_scene_view_uniforms(program, view);

  // set lighting uniforms
  if (view.params.lighting == gui_lighting_type::eyelight) {
    set_eyelight_uniforms(program, view);
  } else {
    set_ibl_uniforms(program, scene);
  }

  set_ogl_wireframe(view.params.wireframe);
  for (auto instance : scene->instances) {
    if (instance->hidden) continue;
    set_instance_uniforms(program, instance->frame, instance->shape,
        instance->material, instance->shading, view.params.double_sided,
        view.params.non_rigid_frames);
    draw_shape(instance->shape);
  }
  unbind_program();
}

gui_scene_view make_scene_view(
    gui_camera* camera, const vec4i& viewport, const gui_scene_params& params) {
  auto camera_aspect = (float)viewport.z / (float)viewport.w;
  auto camera_yfov =
      camera_aspect >= 0
          ? (2 * atan(camera->film / (camera_aspect * 2 * camera->lens)))
          : (2 * atan(camera->film / (2 * camera->lens)));
  auto view_matrix       = frame_to_mat(inverse(camera->frame));
  auto projection_matrix = perspective_mat(
      camera_yfov, camera_aspect, params.near, params.far);

  auto view              = gui_scene_view{};
  view.camera_frame      = camera->frame;
  view.view_matrix       = view_matrix;
  view.projection_matrix = projection_matrix;
  view.params            = params;
  return view;
}

void draw_scene(gui_scene* scene, gui_camera* camera, const vec4i& viewport,
    const gui_scene_params& params) {
  clear_ogl_framebuffer(params.background);
  set_ogl_viewport(viewport);

  auto view = make_scene_view(camera, viewport, params);
  draw_instances(scene, view);
  draw_environment(scene, view);
}

// image based lighting

// Using 6 render passes, bake a cubemap given a sampler for the environment.
// The input sampler can be either a cubemap or a latlong texture.
template <typename Sampler>
inline void bake_cubemap(ogl_cubemap* cubemap, const Sampler* environment,
    ogl_program* program, int size, int num_mipmap_levels = 1,
    const vec3f& emission = {1, 1, 1}) {
  // init cubemap with no data
  set_cubemap<float>(cubemap, size, 3, true, true, true);
  auto cube = ogl_shape{};
  set_cube_shape(&cube);

  auto framebuffer = ogl_framebuffer{};
  set_framebuffer(&framebuffer, {size, size});

  auto cameras = array<frame3f, 6>{
      lookat_frame({0, 0, 0}, {1, 0, 0}, {0, 1, 0}),
      lookat_frame({0, 0, 0}, {-1, 0, 0}, {0, 1, 0}),
      lookat_frame({0, 0, 0}, {0, -1, 0}, {0, 0, -1}),
      lookat_frame({0, 0, 0}, {0, 1, 0}, {0, 0, 1}),
      lookat_frame({0, 0, 0}, {0, 0, -1}, {0, 1, 0}),
      lookat_frame({0, 0, 0}, {0, 0, 1}, {0, 1, 0}),
  };

  bind_framebuffer(&framebuffer);
  bind_program(program);
  for (int mipmap_level = 0; mipmap_level < num_mipmap_levels; mipmap_level++) {
    // resize render buffer and viewport
    set_framebuffer(&framebuffer, {size, size});
    set_ogl_viewport(vec2i{size, size});

    for (auto i = 0; i < 6; ++i) {
      // perspective_mat(fov, aspect, near, far)
      auto camera_proj = perspective_mat(radians(90), 1, 1, 100);
      auto camera_view = frame_to_mat(inverse(cameras[i]));

      set_framebuffer_texture(&framebuffer, cubemap, i, mipmap_level);
      clear_ogl_framebuffer({0, 0, 0, 0}, true);

      set_uniform(program, "view", camera_view);
      set_uniform(program, "projection", camera_proj);
      set_uniform(program, "eye", vec3f{0, 0, 0});
      set_uniform(program, "mipmap_level", mipmap_level);
      set_uniform(program, "emission", emission);
      set_uniform(program, "environment", environment, 0);

      draw_shape(&cube);
    }
    size /= 2;
  }
  unbind_program();
  unbind_framebuffer();
  clear_framebuffer(&framebuffer);
}

inline void bake_specular_brdf_texture(gui_texture* texture) {
  auto size        = vec2i{512, 512};
  auto framebuffer = ogl_framebuffer{};
  auto screen_quad = ogl_shape{};
  set_quad_shape(&screen_quad);

  auto program = ogl_program{};
  auto error = ""s, errorlog = ""s;
  auto vert = bake_brdf_vertex_code();
  auto frag = bake_brdf_fragment_code();
  init_program(&program, vert, frag, error, errorlog);

  set_texture(texture, size, 3, (float*)nullptr, true, true, false, false);

  set_framebuffer(&framebuffer, size);
  set_framebuffer_texture(&framebuffer, texture, 0);

  bind_framebuffer(&framebuffer);
  bind_program(&program);

  set_ogl_viewport(size);
  clear_ogl_framebuffer({0, 0, 0, 0}, true);

  draw_shape(&screen_quad);

  unbind_program();
  unbind_framebuffer();
  clear_framebuffer(&framebuffer);
  clear_program(&program);
  clear_shape(&screen_quad);
}

static void init_environment(gui_scene* scene,
    const gui_texture* environment_tex, const vec3f& environment_emission) {
  set_cube_shape(scene->environment_shape);

  // init program for drawing the environment
  {
    auto vert = bake_cubemap_vertex_code();
    auto frag = draw_enivronment_fragment_code();
    init_program(scene->environment_program, vert, frag);
  }

  // bake cubemap from environment texture
  {
    auto size    = environment_tex->size.y;
    auto program = new ogl_program{};
    auto vert    = bake_cubemap_vertex_code();
    auto frag    = bake_environment_fragment_code();
    init_program(program, vert, frag);
    bake_cubemap(scene->environment_cubemap, environment_tex, program, size, 1,
        environment_emission);
    clear_program(program);
    delete program;
  }
}

void init_ibl_data(gui_scene* scene) {
  // bake irradiance map
  {
    auto program = new ogl_program{};
    auto vert    = bake_cubemap_vertex_code();
    auto frag    = bake_irradiance_fragment_code();
    init_program(program, vert, frag);
    bake_cubemap(
        scene->diffuse_cubemap, scene->environment_cubemap, program, 64);
    clear_program(program);
    delete program;
  }

  // bake specular map
  {
    auto program = new ogl_program{};
    auto vert    = bake_cubemap_vertex_code();
    auto frag    = bake_reflections_fragment_code();
    init_program(program, vert, frag);
    bake_cubemap(
        scene->specular_cubemap, scene->environment_cubemap, program, 256, 6);
    clear_program(program);
    delete program;
  }

  // bake lookup texture for specular brdf
  bake_specular_brdf_texture(scene->brdf_lut);

  // init shader for IBL shading
  {
    auto vert = draw_instances_vertex_code();
    auto frag = draw_instances_ibl_fragment_code();
    init_program(scene->ibl_program, vert, frag);
  }
}

const char* draw_instances_vertex_code() {
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

const char* draw_instances_eyelight_fragment_code() {
  static const char* code =
      R"(
#version 330

float pif = 3.14159265;

uniform bool eyelight;         // eyelight shading
uniform vec3 lamb;             // ambient light
uniform int  lnum;             // number of lights
uniform int  ltype[16];        // light type (0 -> point, 1 -> directional)
uniform vec3 lpos[16];         // light positions
uniform vec3 lke[16];          // light intensities

void evaluate_light(int lid, vec3 position, out vec3 cl, out vec3 wi) {
  cl = vec3(0,0,0);
  wi = vec3(0,0,0);
  if(ltype[lid] == 0) {
    // compute point light color at position
    cl = lke[lid] / pow(length(lpos[lid]-position),2);
    // compute light direction at position
    wi = normalize(lpos[lid]-position);
  }
  else if(ltype[lid] == 1) {
    // compute light color
    cl = lke[lid];
    // compute light direction
    wi = normalize(lpos[lid]);
  }
}

vec3 brdfcos(int etype, vec3 ke, vec3 kd, vec3 ks, float rs, float op,
    vec3 n, vec3 wi, vec3 wo) {
  if(etype == 0) return vec3(0);
  vec3 wh = normalize(wi+wo);
  float ns = 2/(rs*rs)-2;
  float ndi = dot(wi,n), ndo = dot(wo,n), ndh = dot(wh,n);
  if(etype == 1) {
      return ((1+dot(wo,wi))/2) * kd/pif;
  } else if(etype == 2) {
      float si = sqrt(1-ndi*ndi);
      float so = sqrt(1-ndo*ndo);
      float sh = sqrt(1-ndh*ndh);
      if(si <= 0) return vec3(0);
      vec3 diff = si * kd / pif;
      if(sh<=0) return diff;
      float d = ((2+ns)/(2*pif)) * pow(si,ns);
      vec3 spec = si * ks * d / (4*si*so);
      return diff+spec;
  } else if(etype == 3) {
      if(ndi<=0 || ndo <=0) return vec3(0);
      vec3 diff = ndi * kd / pif;
      if(ndh<=0) return diff;
      float cos2 = ndh * ndh;
      float tan2 = (1 - cos2) / cos2;
      float alpha2 = rs * rs;
      float d = alpha2 / (pif * cos2 * cos2 * (alpha2 + tan2) * (alpha2 + tan2));
      float lambda_o = (-1 + sqrt(1 + (1 - ndo * ndo) / (ndo * ndo))) / 2;
      float lambda_i = (-1 + sqrt(1 + (1 - ndi * ndi) / (ndi * ndi))) / 2;
      float g = 1 / (1 + lambda_o + lambda_i);
      vec3 spec = ndi * ks * d * g / (4*ndi*ndo);
      return diff+spec;
  }
}

uniform int etype;
uniform bool faceted;
uniform vec4 highlight;           // highlighted color

uniform int mtype;                // material type
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

uniform bool mat_norm_tex_on;     // material normal texture on
uniform sampler2D mat_norm_tex;   // material normal texture

uniform bool double_sided;        // double sided rendering

uniform mat4 frame;              // shape transform
uniform mat4 frameit;            // shape transform

bool evaluate_material(vec2 texcoord, vec4 color, out vec3 ke, 
                    out vec3 kd, out vec3 ks, out float rs, out float op) {
  ke = color.xyz * emission;
  kd = color.xyz * diffuse;
  ks = color.xyz * specular;
  rs = roughness;
  op = color.w * opacity;

  vec4 ke_tex = (emission_tex_on) ? texture(emission_tex,texcoord) : vec4(1,1,1,1);
  vec4 kd_tex = (diffuse_tex_on) ? texture(diffuse_tex,texcoord) : vec4(1,1,1,1);
  vec4 ks_tex = (specular_tex_on) ? texture(specular_tex,texcoord) : vec4(1,1,1,1);
  vec4 rs_tex = (roughness_tex_on) ? texture(roughness_tex,texcoord) : vec4(1,1,1,1);
  vec4 op_tex = (opacity_tex_on) ? texture(opacity_tex,texcoord) : vec4(1,1,1,1);

  // get material color from textures and adjust values
  ke *= ke_tex.xyz;
  vec3 kb = kd * kd_tex.xyz;
  float km = ks.x * ks_tex.z;
  kd = kb * (1 - km);
  ks = kb * km + vec3(0.04) * (1 - km);
  rs *= ks_tex.y;
  rs = rs*rs;
  op *= kd_tex.w;

  return true;
}

vec3 apply_normal_map(vec2 texcoord, vec3 normal, vec4 tangsp) {
    if(!mat_norm_tex_on) return normal;
  vec3 tangu = normalize((frame * vec4(normalize(tangsp.xyz),0)).xyz);
  vec3 tangv = normalize(cross(normal, tangu));
  if(tangsp.w < 0) tangv = -tangv;
  vec3 texture = 2 * pow(texture(mat_norm_tex,texcoord).xyz, vec3(1/2.2)) - 1;
  // texture.y = -texture.y;
  return normalize( tangu * texture.x + tangv * texture.y + normal * texture.z );
}

in vec3 position;              // [from vertex shader] position in world space
in vec3 normal;                // [from vertex shader] normal in world space (need normalization)
in vec2 texcoord;              // [from vertex shader] texcoord
in vec4 color;                 // [from vertex shader] color
in vec4 tangsp;                // [from vertex shader] tangent space

uniform vec3 eye;              // camera position
uniform mat4 view;             // inverse of the camera frame (as a matrix)
uniform mat4 projection;       // camera projection

uniform float exposure; 
uniform float gamma;

out vec4 frag_color;      

vec3 triangle_normal(vec3 position) {
  vec3 fdx = dFdx(position); 
  vec3 fdy = dFdy(position); 
  return normalize((frame * vec4(normalize(cross(fdx, fdy)), 0)).xyz);
}

// main
void main() {
  if(mtype == 0) {
    frag_color = vec4(emission, 1);
    return; 
  }

  // view vector
  vec3 wo = normalize(eye - position);

  // prepare normals
  vec3 n;
  if(faceted) {
    n = triangle_normal(position);
  } else {
    n = normalize(normal);
  }

  // apply normal map
  n = apply_normal_map(texcoord, n, tangsp);

  // use faceforward to ensure the normals points toward us
  if(double_sided) n = faceforward(n,-wo,n);

  // get material color from textures
  vec3 brdf_ke, brdf_kd, brdf_ks; float brdf_rs, brdf_op;
  bool has_brdf = evaluate_material(texcoord, color, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op);

  // exit if needed
  if(brdf_op < 0.005) discard;

  // check const color
  if(etype == 0) {
    frag_color = vec4(brdf_ke,brdf_op);
    return;
  }

  // emission
  vec3 c = brdf_ke;

  // check early exit
  if(brdf_kd != vec3(0,0,0) || brdf_ks != vec3(0,0,0)) {
    // eyelight shading
    if(eyelight) {
      vec3 wi = wo;
      c += pif * brdfcos((has_brdf) ? etype : 0, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op, n,wi,wo);
    } else {
      // accumulate ambient
      c += lamb * brdf_kd;
      // foreach light
      for(int lid = 0; lid < lnum; lid ++) {
        vec3 cl = vec3(0,0,0); vec3 wi = vec3(0,0,0);
        evaluate_light(lid, position, cl, wi);
        c += cl * brdfcos((has_brdf) ? etype : 0, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op, n,wi,wo);
      }
    }
  }

  // final color correction
  c = pow(c * pow(2,exposure), vec3(1/gamma));

  // highlighting
  if(highlight.w > 0) {
    if(mod(int(gl_FragCoord.x)/4 + int(gl_FragCoord.y)/4, 2)  == 0)
        c = highlight.xyz * highlight.w + c * (1-highlight.w);
  }

  // output final color by setting gl_FragColor
  frag_color = vec4(c,brdf_op);
}
)";
  return code;
}

const char* draw_instances_ibl_fragment_code() {
  static const char* code = R"(
#version 330

in vec3 position;  // [from vertex shader] position in world space
in vec3 normal;    // [from vertex shader] normal in world space
in vec2 texcoord;  // [from vertex shader] texcoord
in vec4 color;     // [from vertex shader] color
in vec4 tangsp;    // [from vertex shader] tangent space

float pif = 3.14159265;

uniform int  etype;
uniform int  mtype;
uniform bool faceted;
uniform int  shading_type;

uniform vec3  emission;   // material ke
uniform vec3  diffuse;    // material kd
uniform vec3  specular;   // material ks
uniform float roughness;  // material rs
uniform float opacity;    // material op

// baked textures for image based lighting
uniform samplerCube irradiance_cubemap;
uniform samplerCube reflection_cubemap;
uniform sampler2D   brdf_lut;

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

uniform bool      mat_norm_tex_on;  // material normal texture on
uniform sampler2D mat_norm_tex;     // material normal texture

uniform bool double_sided;  // double sided rendering

uniform mat4 frame;    // shape transform
uniform mat4 frameit;  // shape transform

uniform vec3 eye;         // camera position
uniform mat4 view;        // inverse of the camera frame (as a matrix)
uniform mat4 projection;  // camera projection

uniform float exposure;
uniform float gamma;

out vec4 frag_color;

struct brdf_struct {
  vec3  emission;
  vec3  diffuse;
  vec3  specular;
  float roughness;
  float opacity;
} brdf;

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

brdf_struct compute_brdf() {
  brdf_struct brdf;
  brdf.emission  = eval_brdf_color(emission, emission_tex, emission_tex_on);
  brdf.diffuse   = eval_brdf_color(diffuse, diffuse_tex, diffuse_tex_on);
  brdf.specular  = eval_brdf_color(specular, specular_tex, specular_tex_on);
  brdf.roughness = eval_brdf_value(roughness, roughness_tex, roughness_tex_on);
  brdf.opacity   = eval_brdf_value(opacity, opacity_tex, opacity_tex_on);
  return brdf;
}

vec3 apply_normal_map(vec2 texcoord, vec3 normal, vec4 tangsp) {
  if (!mat_norm_tex_on) return normal;
  vec3 tangu = normalize((frame * vec4(normalize(tangsp.xyz), 0)).xyz);
  vec3 tangv = normalize(cross(normal, tangu));
  if (tangsp.w < 0) tangv = -tangv;
  vec3 texture = 2 * pow(texture(mat_norm_tex, texcoord).xyz, vec3(1 / 2.2)) -
                 1;
  // texture.y = -texture.y;
  return normalize(tangu * texture.x + tangv * texture.y + normal * texture.z);
}

vec3 triangle_normal(vec3 position) {
  vec3 fdx = dFdx(position);
  vec3 fdy = dFdy(position);
  return normalize((frame * vec4(normalize(cross(fdx, fdy)), 0)).xyz);
}

#define etype_points 1
#define etype_lines 2
#define etype_triangles 3
#define etype_quads 3

vec3 compute_normal(vec3 V) {
  vec3 N;
  if (etype == etype_triangles) {
    if (faceted) {
      N = triangle_normal(position);
    } else {
      N = normalize(normal);
    }
  }

  if (etype == etype_lines) {
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
  return textureLod(reflection_cubemap, L, lod).rgb;
}

// main
void main() {
  if(mtype == 0) {
    frag_color = vec4(emission, 1);
    return; 
  }

  vec3 V = normalize(eye - position);
  vec3 N = compute_normal(V);

  brdf_struct brdf = compute_brdf();
  if (brdf.opacity < 0.005) discard;

  // emission
  vec3 radiance = brdf.emission;

  // diffuse
  radiance += brdf.diffuse * textureLod(irradiance_cubemap, N, 0).rgb;

  // specular
  vec3 L          = normalize(reflect(-V, N));
  vec3 reflection = sample_prefiltered_refleciton(L, brdf.roughness);
  vec2 env_brdf   = texture(brdf_lut, vec2(max(dot(N, V), 0.0), roughness)).rg;
  radiance += reflection * (brdf.specular * env_brdf.x + env_brdf.y);

  // final color correction
  radiance = pow(radiance * pow(2, exposure), vec3(1 / gamma));

  // output final color by setting gl_FragColor
  frag_color = vec4(radiance, brdf.opacity);
}
)";
  return code;
}

static const char* bake_brdf_vertex_code() {
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

static const char* bake_brdf_fragment_code() {
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

static const char* bake_cubemap_vertex_code() {
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

static const char* bake_environment_fragment_code() {
  static const char* code = R"(
#version 330

out vec3 frag_color;

in vec3 position;  //  position in world space

uniform vec3 eye;         // camera position
uniform mat4 view;        // inverse of the camera frame (as a matrix)
uniform mat4 projection;  // camera projection

uniform sampler2D environment;
uniform vec3 emission = vec3(1);

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

  // TODO(giacomo): We skip gamma correction, assuming the environment is stored
  // in linear space. Is it always true? Probably not.
  frag_color = emission * color;
}
)";
  return code;
}

static const char* bake_irradiance_fragment_code() {
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

static const char* bake_reflections_fragment_code() {
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

const char* draw_enivronment_fragment_code() {
  static const char* code = R"(
#version 330

out vec3 frag_color;

in vec3 position;  //  position in world space

uniform vec3 eye;         // camera position
uniform mat4 view;        // inverse of the camera frame (as a matrix)
uniform mat4 projection;  // camera projection

uniform float exposure;
uniform float gamma;

uniform samplerCube environment;

void main() {
  vec3 v = normalize(position);

  vec3 radiance = texture(environment, v).rgb;

  // final color correction
  radiance   = pow(radiance * pow(2, exposure), vec3(1 / gamma));
  frag_color = radiance;
}
)";
  return code;
}

}  // namespace yocto
